!******************************************************************************!
!                    dk_structure_io_xyz submodule
!------------------------------------------------------------------------------!
!> Specific procedures to read and write XYZ and extended XYZ structure files.
!
!  Part of the dk_io library version 0.1.
!  Written by Paul Fossati, <paul.fossati@cea.fr>
!  Copyright (c) 2024-2025 CEA
!
!------------------------------------------------------------------------------!
! Redistribution and use in source and binary forms, with or without           !
! modification, are permitted provided that the following conditions are met:  !
!                                                                              !
!     * Redistributions of source code must retain the above copyright notice, !
!       this list of conditions and the following disclaimer.                  !
!                                                                              !
!     * Redistributions in binary form must reproduce the above copyright      !
!       notice, this list of conditions and the following disclaimer in the    !
!       documentation and/or other materials provided with the distribution.   !
!                                                                              !
!     * The name of the author may not be used to endorse or promote products  !
!      derived from this software without specific prior written permission    !
!      from the author.                                                        !
!                                                                              !
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"  !
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    !
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE   !
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE     !
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR          !
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF         !
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS     !
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      !
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)      !
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE   !
! POSSIBILITY OF SUCH DAMAGE.                                                  !
!******************************************************************************!
submodule (dk_structure_io) dk_structure_io_xyz
    use dk_math,       only: matrix_inverse
    use dk_string,     only: FString
    use dk_datatable,  only: DataTable
    use dk_datacolumn, only: DataColumn

    use dk_exception, only: FException, FExceptionDescription
    use dk_token, only: Token, get_token, TokenError
    use ext_character

    implicit none(external, type)

contains
!******************************************************************************!
!> Read a structure from a XYZ file.
!******************************************************************************!
    module subroutine read_xyz(input_file, box, positions, tags, masses, velocities, forces, aux, info, stat)
        type(FileObject), intent(inout) :: input_file
        real(e), dimension(3,3), intent(out) :: box
        real(e), dimension(:,:), allocatable, intent(out) :: positions
        character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
        real(e), dimension(:),   allocatable, intent(out), optional :: masses
        real(e), dimension(:,:), allocatable, intent(out), optional :: velocities
        real(e), dimension(:,:), allocatable, intent(out), optional :: forces
        type(DataTable), intent(out), optional :: aux
        type(FileInfo), intent(out), optional :: info
        type(FException), intent(out), optional :: stat
!------
        integer :: a, i, iostat, nFields, atoms_count, start, end
        integer :: x, y, z, vx, vy, vz, fx, fy, fz, mass, read
        character(:), allocatable :: buffer
        type(FException) :: stat_, istat
        real(e), dimension(:), allocatable :: data_array
        real(e), dimension(6) :: parameters
        real(e), dimension(9) :: box_definition
        real(e), dimension(3,3) :: ibox
        type(Token) :: error_token
!------

        nFields = 0

        x = 0
        y = 0
        z = 0
        vx = 0
        vy = 0
        vz = 0
        fx = 0
        fy = 0
        fz = 0
        mass = 0

        body: block

!------ Get the number of atoms
            call input_file%read_line(buffer, stat=stat_)
            read(buffer,*,iostat=iostat) atoms_count
            if (iostat /= 0) then
                call istat%raise("The first line is not a valid number of atoms.")
                exit body
            end if

!------ Process the header
            call input_file%read_line(buffer, stat=stat_)
            read(buffer,*,iostat=iostat) box_definition
            if (iostat == 0) then
                box(:,1) = box_definition(1:3)
                box(:,2) = box_definition(4:6)
                box(:,3) = box_definition(7:9)
            else
                read(buffer,*,iostat=iostat) parameters
                if (iostat == 0) then
                    box = box_with_parameters(parameters)
                else
                    call parse_line(buffer, istat)
                    if (istat /= 0) then
                        exit body
                    end if

                    ! Count the number of auxiliary fields

                end if
            end if
            if (nFields == 0) then
                ! Read the file as standard xyz
                x = 1
                y = 2
                z = 3
                nFields = 3
            end if

!------ Read the atomic positions
            allocate(tags(atoms_count))
            allocate(positions(3,atoms_count))
            if (present(velocities)) allocate(velocities(3,atoms_count))
            if (present(forces)) allocate(forces(3,atoms_count))
            if (present(masses)) allocate(masses(atoms_count))
            allocate(data_array(nFields))
            do i=1, atoms_count
                call input_file%read_line(buffer, stat=stat_)
                if (stat_ /= 0) then
                    call istat%raise("Unexpected end of file.")
                    exit body
                end if

                ! Read atom tag
                start = 1
                do while(start < len(buffer) .and. is_whitespace(buffer(start:start)))
                    start = start + 1
                end do
                end = start + 1
                do while (end < len(buffer) .and. .not. is_whitespace(buffer(end:end)))
                    end = end + 1
                end do

                tags(i) = buffer(start:end-1)
                buffer(start:end-1) = " "

                call read_array(buffer, data_array, read)

                if (read < size(data_array)) then
                    call get_token(buffer, read+1, error_token)
                    error_token%line_number = input_file%current_line
                    if (.not.allocated(error_token%value)) then
                        error_token%start = len_trim(buffer)
                        error_token%line = trim(buffer)
                        error_token%value = buffer(len_trim(buffer):len_trim(buffer))
                        call stat%raise(TokenError("Missing data on line.", &
                            "In XYZ file "//bold(input_file%name)//":", error_token))
                        exit body
                    else
                        call stat%raise(TokenError("Invalid number value.", &
                            "In XYZ file "//bold(input_file%name)//":", error_token))
                        exit body
                    end if
                end if

                positions(1,i) = data_array(x)
                positions(2,i) = data_array(y)
                positions(3,i) = data_array(z)

                if (present(velocities)) then
                    velocities(1,i) = data_array(vx)
                    velocities(2,i) = data_array(vy)
                    velocities(3,i) = data_array(vz)
                end if

                if (present(forces)) then
                    forces(1,i) = data_array(fx)
                    forces(2,i) = data_array(fy)
                    forces(3,i) = data_array(fz)
                end if

                if (present(masses) .and. mass > 0) masses(i) = data_array(mass)
            end do

            ! Get fractional coordinates
            ibox = matrix_inverse(box)
            do i=1, atoms_count
                positions(:,i) = matmul(ibox, positions(:,i))
            end do

            ! Normalise shell names
            do i=1, atoms_count
                a = index(tags(i), "-Shell")
                if (a /= 0) then
                    tags(i)(a:a) = "_"
                end if
            end do

            if (present(info)) then
                info%data_format = trim(STRUCTURE_XYZ)
                info%num_frames = 1
                info%compression_format = input_file%compression_format()
            end if
        end block body

        ! Cleanup
        call stat_%discard
        if (present(stat)) call stat%transfer(istat)
    contains
!******************************************************************************!
!> Parse properties from the property line of an extended XYZ file.
!******************************************************************************!
        subroutine parse_line(line, pstat)
            character(*), intent(in) :: line
            type(FException), intent(out) :: pstat
!------
            integer :: cursor
!------
            cursor = 1

            do while(cursor <= len(line))
                call parse_property(line, cursor, pstat)
                cursor = cursor + 1
                if (pstat /= 0) return
                if (cursor > len(line)) exit

                do while(is_whitespace(line(cursor:cursor)))
                    cursor = cursor + 1
                    if (cursor > len(line)) exit
                end do
            end do
        end subroutine
!******************************************************************************!
!> Parse a property name and value from the property line of an extended XYZ
!> file.
!******************************************************************************!
        subroutine parse_property(line, cursor, pstat)
            character(*), intent(in) :: line
            integer, intent(inout) :: cursor
            type(FException), intent(out) :: pstat
!------
            integer :: start, n, i, j
            character(:), allocatable :: mode, name, value
            character(:), dimension(:), allocatable :: fields
            type(Token) :: error_token
!------

!------ Get the name of the property
            start = cursor
            do while(cursor <= len_trim(line))
                if (line(cursor:cursor) == "=" .or. is_whitespace(line(cursor:cursor))) exit
                cursor = cursor + 1
            end do
            name = line(start:cursor-1)

            ! Skip the whitespace between the property name and the = sign
            do while (cursor <= len_trim(line))
                if (.not.is_whitespace(line(cursor:cursor))) exit
                cursor = cursor + 1
            end do

            ! No = sign was found on the line
            if (cursor > len_trim(line)) then
                call pstat%raise("Expected '=' sign after property name.")
                return
            end if

            ! Got random non-whitespace character instead of an = sign
            if (line(cursor:cursor) /= "=") then
                call pstat%raise("Expected '=' sign after property name.")
                return
            end if
            cursor = cursor + 1

            ! Skip any whitespace between the = sign and the property value
            do while (cursor <= len_trim(line))
                if (.not.is_whitespace(line(cursor:cursor))) exit
                cursor = cursor + 1
            end do

            ! No value was found after the = sign
            if (cursor > len_trim(line)) then
                call pstat%raise("Expected value after '=' sign.")
                return
            end if

!------ Get the value of the property
            if (line(cursor:cursor) == """") then
                mode = "double"
                cursor = cursor + 1
            else
                mode = "bareword"
            end if

            start = cursor
            do while(cursor < len_trim(line))
                cursor = cursor + 1
                if (mode=="double" .and. line(cursor:cursor)=="""") then
                    value = line(start:cursor-1)
                    exit
                end if
                if (mode=="bareword" .and. is_whitespace(line(cursor:cursor))) then
                    value = line(start:cursor-1)
                    exit
                end if
                if (mode=="bareword" .and. cursor == len_trim(line)) then
                    value = line(start:cursor)
                    exit
                end if
            end do

            if (cursor == len_trim(line) .and. mode == "double" .and. line(cursor:cursor) /= """") then
                error_token%start = start
                error_token%line = line
                error_token%value = line(start:cursor)
                error_token%line_number = input_file%current_line
                call stat%raise(TokenError("Unterminated double-quoted string.", &
                    "In XYZ file "//bold(input_file%name)//":", error_token))
                return
            end if

            if (.not.allocated(value)) return

!------ Process the property
            select case(name)
            case ("Lattice")
                read(value,*,iostat=iostat) box_definition
                if (iostat /= 0) then
                    error_token%start = start
                    error_token%line = line
                    error_token%value = value
                    error_token%line_number = input_file%current_line
                    call stat%raise(TokenError("Invalid cell shape.", &
                        "In XYZ file "//bold(input_file%name)//":", error_token))
                    return
                end if
                box(:,1) = box_definition(1:3)
                box(:,2) = box_definition(4:6)
                box(:,3) = box_definition(7:9)
            case ("Origin")
            case ("Properties")
                if (value(:12) /= "species:S:1:") then
                    error_token%start = start
                    error_token%line = line
                    error_token%value = line(start:start+index(line(start:cursor),":")-1)
                    error_token%line_number = input_file%current_line
                    call stat%raise(TokenError("Extended XYZ files are supported only if the first column contains the species names.", &
                        "In XYZ file "//bold(input_file%name)//":", error_token))
                    return
                end if
                call split(value(13:), ":", fields)
                n = 0
                do i=1, size(fields)/3

                    read(fields((i-1)*3+3),*,iostat=iostat) j
                    if (iostat /= 0) then
                        call pstat%raise("Syntax error in comment line")
                        return
                    end if

                    if (fields((i-1)*3+2) /= "R") then
                        call pstat%raise("Only real auxiliary fields are supported in extended XYZ files.")
                        return
                    end if

                    select case(lower_case(fields((i-1)*3+1)))
                    case ("velo", "velocity", "velocities", "vel")
                        vx = n+1
                        vy = n+2
                        vz = n+3
                    case ("force", "forces")
                        fx = n+1
                        fy = n+2
                        fz = n+3
                    case ("pos", "position", "positions")
                        x = n+1
                        y = n+2
                        z = n+3
                    case ("mass", "masses")
                        mass = n+1
                    case default
                        if (j > 1) then
                            call pstat%raise("Only scalar auxiliary fields are supported in extended XYZ files.")
                            return
                        end if
                    end select

                    n = n + j
                end do

                nFields = n
            end select
        end subroutine
    end subroutine
!******************************************************************************!
!> Write a structure to a XYZ file.
!******************************************************************************!
    module subroutine write_xyz(filename, box, positions, tags, masses, velocities, forces, aux, include, scale, stat)
        character(*), intent(in) :: filename
        real(e), dimension(3,3), intent(in) :: box
        real(e), dimension(:,:), intent(in) :: positions
        character(TAG_LENGTH), dimension(:), intent(in) :: tags
        real(e), dimension(:), intent(in), optional :: masses
        real(e), dimension(:,:), intent(in), optional :: velocities
        real(e), dimension(:,:), intent(in), optional :: forces
        type(DataTable), intent(in), optional :: aux
        logical, dimension(:), intent(in), optional :: include
        real(e), intent(in), optional :: scale
        type(FException), intent(out), optional :: stat
!------
        integer :: i, c, cols, atoms_count, n, fields_count
        logical, dimension(:), allocatable :: mask
        character(10000) :: buffer, properties
        real(c_double), dimension(:), allocatable :: tmpdat
        type(FileObject) :: output_file
        real(e), dimension(:), allocatable :: row
        type(FException) :: istat
!------

        body: block

!------ Prepare the file
            call output_file%init(filename, "w", istat)

            ! Return if the file could not be opened
            if (istat /= 0) then
                exit body
            end if

            if (size(tags) /= size(positions, 2)) then
                call istat%raise(FExceptionDescription("The tags array has the wrong size.", &
                    "When writing the extended XYZ file " // bold(filename) // ":"))
                exit body
            end if

            ! Set whether we need to save atomic velocities
            if (present(velocities)) then
                if (size(velocities, 2) /= size(positions, 2) .or. size(velocities, 1) /= 3) then
                    call istat%raise(FExceptionDescription("The velocities array has the wrong size.", &
                        "When writing the XYZ file " // bold(filename) // ":"))
                    exit body
                end if
            end if

            if (present(forces)) then
                if (size(forces, 2) /= size(positions, 2) .or. size(forces, 1) /= 3) then
                    call istat%raise(FExceptionDescription("The forces array has the wrong size.", &
                        "When writing the XYZ file " // bold(filename) // ":"))
                    exit body
                end if
            end if

            if (present(masses)) then
                if (size(masses) /= size(positions, 2)) then
                    call istat%raise(FExceptionDescription("The masses array has the wrong size.", &
                        "When writing the XYZ file " // bold(filename) // ":"))
                    exit body
                end if
            end if

!------ Build the inclusion mask
            if (present(include)) then
                if (size(include) /= size(positions, 2)) then
                    call istat%raise(FExceptionDescription("The atomic mask has the wrong size.", &
                        "When writing the extended XYZ file " // bold(filename) // ":"))
                    exit body
                end if
                allocate(mask, source=include)
            else
                allocate(mask(size(positions, 2)))
                mask = .true.
            end if
            atoms_count = count(mask)

            fields_count = 3  ! Number of coordinates
            if (present(velocities)) fields_count = fields_count + 3  ! Velocities
            if (present(forces)) fields_count = fields_count + 3      ! Forces
            if (present(masses)) fields_count = fields_count + 1      ! Atomic masses
            if (present(aux)) cols = fields_count + aux%num_columns() ! Auxiliary fields
            allocate(tmpdat(fields_count))

!------ Print header
            call output_file%write_line("" // atoms_count)
            properties = 'Properties=species:S:1:pos:R:3'
            if (present(velocities)) properties = trim(properties) // ":velo:R:3"
            if (present(forces)) properties = trim(properties)     // ":force:R:3"
            if (present(masses)) properties = trim(properties)     // ":mass:R:1"
            if (present(aux)) then
                do c=1, aux%num_columns()
                    properties = trim(properties) // ":" // aux%column_name(c) // ":R:1"
                end do
            end if
            write(buffer,'(a,g0,8(1x,g0),a)') 'Lattice="', box(:,1), box(:,2), box(:,3), '" ' // trim(properties)
            call output_file%write_line(trim(buffer))

!------ Special case: empty configuration
            if (atoms_count == 0) then
                call output_file%close
                exit body
            end if

!------ Print atomic data
            do i=1, size(positions, 2)
                if (.not.mask(i)) cycle

                tmpdat(:3) = matmul(box, positions(:,i))
                n = 4
                if (present(velocities)) then
                    tmpdat(n:n+2) = velocities(:,i)
                    n = n + 3
                end if
                if (present(forces)) then
                    tmpdat(n:n+2) = forces(:,i)
                    n = n + 3
                end if
                if (present(masses)) then
                    tmpdat(n:n) = masses(i)
                    n = n + 1
                end if
                if (present(aux)) then
                    call aux%get_row(i, row)
                    if (allocated(row)) tmpdat(n+1:) = row(:)
                end if

                n = writeArray(buffer, tmpdat, int(size(tmpdat), c_size_t), len(buffer))
                call output_file%write_line(tags(i) // " " // trim(buffer))
            end do

            call output_file%close
        end block body

        if (present(stat)) call stat%transfer(istat)
    end subroutine
end submodule
