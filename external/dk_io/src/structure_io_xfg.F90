!******************************************************************************!
!                    dk_structure_io_xfg submodule
!------------------------------------------------------------------------------!
!> Specific procedures to read and write structure files in Atomeye's extended
!> CFG format.
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
submodule (dk_structure_io) dk_structure_io_xfg
    use dk_datatable, only: DataTable
    use dk_datacolumn, only: DataColumn
    use dk_string, only: FString
    
    use dk_exception, only: FException, FExceptionDescription
    use ext_character, only: operator(//), bold, colour, upper_case, is_whitespace

    implicit none(external, type)

contains

    subroutine get_error_context(line, value, line_number, context, start)
        character(*), intent(in) :: line
        character(*), intent(in) :: value
        integer, intent(in) :: line_number
        character(:), allocatable, intent(out) :: context
        integer, intent(in), optional :: start
!------
        character(100) :: line_number_buffer, empty, caret
        integer :: line_number_len, line_len, total_len, i, start_
!------
        ! Get the size of the terminal

        ! If the output is not a terminal, then use a size of 100 characters
        total_len = 80

        ! Get the number of positions for the line number
        line_number_buffer = " "
        empty = " "
        write(line_number_buffer,'(i0)') line_number
        line_number_len = len_trim(line_number_buffer)

        ! Get the length of the line to print
        line_len = total_len - 4 - line_number_len

        ! Get the caret to indicate the error location
        caret = " "
        if (len_trim(value)==1) then
            caret = colour("^", "red")
            start_ = index(line, value)
        else if (len_trim(value) > 1) then
            do i=1, len_trim(value)
                caret(i:i) = "~"
            end do
            caret = colour(trim(caret), "red")
            start_ = index(line, value)
        else if (len_trim(line) > 0) then
            do i=1, len_trim(line)
                caret(i:i) = "~"
            end do
            caret = colour(trim(caret), "red")
            start_ = 1
        end if
        if (present(start)) start_ = start

        ! Get the context
        context = " "//trim(line_number_buffer)//" | "//trim(line)//new_line("")// &
            empty(:1+line_number_len)//" | "//empty(:start_-1)//trim(caret)
    end subroutine
!******************************************************************************!
!> Get a human-readable description of an error that happened while reading an
!> extended CFG file.
!******************************************************************************!
    function XFGFileError(file, message, line, value)
        type(FileObject), intent(in) :: file
        character(*), intent(in) :: message
        character(*), intent(in) :: line
        character(*), intent(in), optional :: value
        type(FExceptionDescription) :: XFGFileError
!------
        character(:), allocatable :: context
!------
        if (present(value)) then
            call get_error_context(line, value, file%current_line, context)
        else
            call get_error_context(line, "", file%current_line, context)
        end if
        XFGFileError = FExceptionDescription(message, &
            "In extended CFG file "//bold(file%full_name())//":"//new_line("")//context)
    end function
!******************************************************************************!
!> Read a structure from an extended CFG file.
!******************************************************************************!
    module subroutine read_xfg(input_file, box, positions, tags, masses, velocities, forces, aux, info, stat)
        type(FileObject), intent(inout) :: input_file
        real(e), dimension(3,3), intent(out) :: box
        real(e), dimension(:,:), allocatable, intent(out) :: positions
        character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
        real(e), dimension(:), allocatable, intent(out), optional :: masses
        real(e), dimension(:,:), allocatable, intent(out), optional :: velocities
        real(e), dimension(:,:), allocatable, intent(out), optional :: forces
        type(DataTable), intent(out), optional :: aux
        type(FileInfo), intent(out), optional :: info
        type(FException), intent(out), optional :: stat
!------
        integer :: a, i, j, iostat, n, nFields, atoms_count, begin, end
        integer :: fx, fy, fz, vx, vy, vz, read
        character(:), allocatable :: buffer, comp, keyword, value
        type(FException) :: stat_, istat
        real(e) :: x
        logical :: start
        logical :: velocityMode, aux_mode
        real(e), dimension(:), allocatable :: dataArray
        character(TAG_LENGTH) :: current_tag
        real(e) :: current_mass, current_aux
        integer, dimension(:), allocatable :: aux_index
        type(FString), dimension(:), allocatable :: headers
!------

        velocityMode = .true.
        nFields = 0
        aux_mode = present(aux)
        fx = 0
        fy = 0
        fz = 0
        vx = 0
        vy = 0
        vz = 0
        x = 1

        body: block

!------ Process the header
            call input_file%read_line(buffer, stat=stat_)
            start = .true.
            current_aux = 0

            headerLoop: do while(stat_ == 0)

                ! Get a copy of the buffer
                comp = buffer

                ! Remove comments
                do i=1, len(comp)
                    if (comp(i:i) == "#") then
                        comp(i:) = ""
                        exit
                    end if
                end do

                ! Ignore empty lines
                if (comp == "") then
                    call input_file%read_line(buffer, stat=stat_)
                    cycle
                end if

                ! Set velocity mode
                if (upper_case(comp) == ".NO_VELOCITY.") then
                    velocityMode = .false.
                    call input_file%read_line(buffer, stat=stat_)
                    cycle
                end if

                ! Read a keyword-value pair
                if (index(comp, "=") > 0) then
                    keyword = trim(comp(:index(comp, "=")-1))
                    value = trim(adjustl(comp(index(comp, "=")+1:)))

                    ! The first non-comment, non-empty line needs to start with "number of particles"
                    if (start) then
                        if (upper_case(keyword) == "NUMBER OF PARTICLES") then

                            read(value,*,iostat=iostat) atoms_count
                            if (iostat /= 0) then
                                call stat%raise(XFGFileError(input_file, &
                                    "Invalid value for the number of particles.", buffer, &
                                    value))
                                exit body
                            end if
                            start = .false.

                            call input_file%read_line(buffer, stat=stat_)
                            cycle
                        else
                            call stat%raise(XFGFileError(input_file, &
                                "The first non-comment line does not have the number of particles.", buffer))
                            exit body
                        end if
                    end if

                    ! Set auxiliary names
                    if (index(keyword, "auxiliary[") /= 0) then
                        if (.not.allocated(headers)) then
                            call stat%raise(XFGFileError(input_file, &
                                "Missing entry_count line.", buffer))
                            exit body
                        end if
                        i = index(keyword, "[")
                        j = index(keyword, "]")

                        read(keyword(i+1:j-1),*,iostat=iostat) a
                        if (iostat /= 0) then
                            call stat%raise(XFGFileError(input_file, &
                                "Invalid auxiliary index.", buffer, keyword(i+1:j-1)))
                            exit body
                        end if
                        if (a > size(headers)+1) then
                            call stat%raise(XFGFileError(input_file, &
                                "Auxiliary index out of bounds.", buffer, keyword(i+1:j-1)))
                            exit body
                        end if

                        headers(a+1) = value

                        call input_file%read_line(buffer, stat=stat_)
                        cycle
                    end if

                    ! Process other keyword-value pairs
                    select case(upper_case(keyword))
                    case ("H0(1,1)")
                        read(value,*,iostat=iostat) box(1,1)
                        if (iostat /= 0) then
                            call stat%raise(XFGFileError(input_file, &
                                "Invalid value for H0(1,1).", buffer, &
                                value))
                            exit body
                        end if

                    case ("H0(1,2)")
                        read(value,*,iostat=iostat) box(1,2)
                        if (iostat /= 0) then
                            call stat%raise(XFGFileError(input_file, &
                                "Invalid value for H0(1,2).", buffer, &
                                value))
                            exit body
                        end if

                    case ("H0(1,3)")
                        read(value,*,iostat=iostat) box(1,3)
                        if (iostat /= 0) then
                            call stat%raise(XFGFileError(input_file, &
                                "Invalid value for H0(1,3).", buffer, &
                                value))
                            exit body
                        end if

                    case ("H0(2,1)")
                        read(value,*,iostat=iostat) box(2,1)
                        if (iostat /= 0) then
                            call stat%raise(XFGFileError(input_file, &
                                "Invalid value for H0(2,1).", buffer, &
                                value))
                            exit body
                        end if

                    case ("H0(2,2)")
                        read(value,*,iostat=iostat) box(2,2)
                        if (iostat /= 0) then
                            call stat%raise(XFGFileError(input_file, &
                                "Invalid value for H0(2,2).", buffer, &
                                value))
                            exit body
                        end if

                    case ("H0(2,3)")
                        read(value,*,iostat=iostat) box(2,3)
                        if (iostat /= 0) then
                            call stat%raise(XFGFileError(input_file, &
                                "Invalid value for H0(2,3).", buffer, &
                                value))
                            exit body
                        end if

                    case ("H0(3,1)")
                        read(value,*,iostat=iostat) box(3,1)
                        if (iostat /= 0) then
                            call stat%raise(XFGFileError(input_file, &
                                "Invalid value for H0(3,1).", buffer, &
                                value))
                            exit body
                        end if

                    case ("H0(3,2)")
                        read(value,*,iostat=iostat) box(3,2)
                        if (iostat /= 0) then
                            call stat%raise(XFGFileError(input_file, &
                                "Invalid value for H0(3,2).", buffer, &
                                value))
                            exit body
                        end if

                    case ("H0(3,3)")
                        read(value,*,iostat=iostat) box(3,3)
                        if (iostat /= 0) then
                            call stat%raise(XFGFileError(input_file, &
                                "Invalid value for H0(3,3).", buffer, &
                                value))
                            exit body
                        end if

                    case ("A")
                        read(value,*,iostat=iostat) x
                        if (iostat /= 0) then
                            call stat%raise(XFGFileError(input_file, &
                                "Invalid value for A.", buffer, &
                                value))
                            exit body
                        end if

                    case ("ENTRY_COUNT")
                        read(value,*,iostat=iostat) nFields
                        if (iostat /= 0 .or. nFields < 0) then
                            call stat%raise(XFGFileError(input_file, &
                                "Invalid value for entry_count.", buffer, &
                                value))
                            exit body
                        end if
                        if (velocitymode) then
                            if (nFields < 6) then
                                call stat%raise(XFGFileError(input_file, &
                                    "At least 6 fields are necessary for atomic positions and velocities.", buffer, &
                                    value))
                                exit body
                            end if
                            allocate(headers(nFields-6))
                        else
                            if (nFields < 3) then
                                call stat%raise(XFGFileError(input_file, &
                                    "At least 3 fields are necessary for atomic positions.", buffer, &
                                    value))
                                exit body
                            end if
                            allocate(headers(nFields-3))
                        end if
                        if (present(aux)) allocate(aux_index(size(headers)))

                    case default
                        call stat%raise(XFGFileError(input_file, &
                            "Invalid keyword.", buffer, &
                            keyword))
                        exit body
                    end select

                    call input_file%read_line(buffer, stat=stat_)
                    cycle
                end if

                ! Check whether this is the definition of the first atom
                read(comp,*,iostat=iostat) current_mass
                if (iostat == 0) then
                    exit
                else
                    call stat%raise(XFGFileError(input_file, &
                        "Invalid mass definition for the first particle.", buffer))
                    exit body
                end if
            end do headerLoop

            ! Multiplication by the scale factor
            box=x*box

            ! Make sure the entry_count line is present
            if (nFields == 0) then
                call stat%raise(XFGFileError(input_file, &
                    "Missing entry_count line.", buffer))
                exit body
            end if

            ! Get the indices of the forces and velocities columns
            if (velocitymode) then
                vx = 4
                vy = 5
                vz = 6
            end if
            if (allocated(headers)) then
                do i=1, size(headers)
                    select case(lower_case(headers(i)%string()))
                    case ("force.x", "forces.x", "fx")
                        fx = i+3
                        if (velocityMode) fx = fx + 3
                    case ("force.y", "forces.y", "fy")
                        fy = i+3
                        if (velocityMode) fy = fy + 3
                    case ("force.z", "forces.z", "fz")
                        fz = i+3
                        if (velocityMode) fz = fz + 3
                    case ("v.x", "vel.x", "velo.x", "velocity.x", "velocities.x", "vx")
                        vx = i+3
                        if (velocityMode) vx = vx + 3
                    case ("v.y", "vel.y", "velo.y", "velocity.y", "velocities.y", "vy")
                        vy = i+3
                        if (velocityMode) vy = vy + 3
                    case ("v.z", "vel.z", "velo.z", "velocity.z", "velocities.z", "vz")
                        vz = i+3
                        if (velocityMode) vz = vz + 3
                    end select
                end do
            end if

            ! Prepare the data arrays
            allocate(positions(3,atoms_count))
            allocate(tags(atoms_count))
            if (atoms_count > 0) then
                positions(:,:) = 0
                tags(:) = " "
            end if
            if (present(masses)) allocate(masses(atoms_count))
            if (present(velocities)) allocate(velocities(3,atoms_count))
            if (present(forces)) allocate(forces(3,atoms_count))
            if (present(aux)) then
                do i=1, size(headers)
                    call aux%add_column(headers(i)%string(), "", length=atoms_count)
                end do
            end if

!------ Process the atomic data

            ! Check whether the configuration is empty
            if (stat_ == "End of file.") then
                exit body
            end if

            ! Prepare the data arrays
            allocate(dataArray(nFields))

            ! Read the species of the first atom
            call input_file%read_line(buffer, stat=stat_)
            begin = 1
            do while(begin < len(buffer) .and. is_whitespace(buffer(begin:begin)))
                begin = begin + 1
            end do
            end = begin + 1
            do while (end < len(buffer) .and. .not. is_whitespace(buffer(end:end)))
                end = end + 1
            end do
            current_tag = buffer(begin:end-1)

            ! Read all the atoms
            i = 0
            call input_file%read_line(buffer, stat=stat_)
            do while(stat_ == 0 .and. i < atoms_count)

                ! Handle errors when getting the current line
                if (stat_ == "End of file.") then
                    write(*,*) " Reached end of file"
                    call istat%discard
                    call istat%raise(FileError(input_file, &
                        message="Missing atoms.", context="Expecting " // n // " atoms (" // i // " present)"))
                end if
                if (stat_ /= 0) exit

                ! Make sure the buffer is NULL_terminated
                if (len_trim(buffer) < len(buffer)) then
                    buffer(len_trim(buffer)+1:len_trim(buffer)+1) = C_NULL_CHAR
                else
                    buffer = trim(buffer) // C_NULL_CHAR
                end if

                ! Try to read an array of numbers
                call read_array(buffer, dataArray, read)

                if (read == 1) then
                    ! There is only a single numerical value, which should be the mass of the next atom
                    current_mass = dataArray(1)

                    ! The next line should be the species of the next atom
                    call input_file%read_line(buffer, stat=stat_)
                    begin = 1
                    do while(begin < len(buffer) .and. is_whitespace(buffer(begin:begin)))
                        begin = begin + 1
                    end do
                    end = begin + 1
                    do while (end < len(buffer) .and. .not. is_whitespace(buffer(end:end)))
                        end = end + 1
                    end do
                    current_tag = buffer(begin:end-1)

                    ! The line after that should be the definition of the next atom
                    call input_file%read_line(buffer, stat=stat_)
                    if (len_trim(buffer) < len(buffer)) then
                        buffer(len_trim(buffer)+1:len_trim(buffer)+1) = C_NULL_CHAR
                    else
                        buffer = trim(buffer) // C_NULL_CHAR
                    end if
                    call read_array(buffer, dataArray, read)
                end if

                ! There is a syntax error if there is not the right amount of data on the current line
                if (read /= size(dataArray)) then
                    write(*,*) ":"//trim(buffer)//":"
                    write(*,*) read, size(dataArray)
                    stop
                end if

                ! Set the next atom
                i = i + 1
                tags(i) = current_tag
                positions(:,i) = dataArray(1:3)
                if (present(masses)) masses(i) = current_mass
                if (present(forces) .and. fx>0) forces(1,i) = dataArray(fx)
                if (present(forces) .and. fy>0) forces(2,i) = dataArray(fy)
                if (present(forces) .and. fz>0) forces(3,i) = dataArray(fz)
                if (present(velocities) .and. vx>0) velocities(1,i) = dataArray(vx)
                if (present(velocities) .and. vy>0) velocities(2,i) = dataArray(vy)
                if (present(velocities) .and. vz>0) velocities(3,i) = dataArray(vz)
                if (present(aux) .and. read > 3) then
                    do a=1, int(read-3)
                        call aux%append(a, dataArray(3+a))
                    end do
                end if

                call input_file%read_line(buffer, stat=stat_)
            end do

            ! Normalise shell names
            do i=1, atoms_count
                a = index(tags(i), "-Shell")
                if (a /= 0) then
                    tags(i)(a:a) = "_"
                end if
            end do

            if (present(info)) then
                info%data_format = trim(STRUCTURE_XFG)
                info%num_frames = 1
                info%compression_format = input_file%compression_format()
            end if
        end block body

        ! Cleanup
        call stat_%discard
    end subroutine
!******************************************************************************!
!> Write a structure to an extended CFG file.
!******************************************************************************!
    module subroutine write_xfg(filename, box, positions, tags, masses, velocities, forces, aux, include, scale, stat)
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
        integer :: i, c, cols, atoms_count, n
        logical :: has_auxiliaries, has_velocities, has_forces
        logical, dimension(:), allocatable :: mask
        character(:), allocatable :: last_tag
        character(1000) :: buffer
        real(c_double), dimension(:), allocatable :: tmpdat
        type(FileObject) :: output_file
        real(e), dimension(:), allocatable :: row
        type(FException) :: istat
        type(FString), dimension(:), allocatable :: headers
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
                    "When writing the extended CFG file " // bold(filename) // ":"))
                exit body
            end if

            ! Set whether we need to save atomic velocities
            has_velocities = present(velocities)
            if (has_velocities) then
                if (size(velocities, 2) /= size(positions, 2) .or. size(velocities, 1) /= 3) then
                    call istat%raise(FExceptionDescription("The velocities array has the wrong size.", &
                        "When writing the extended CFG file " // bold(filename) // ":"))
                    exit body
                end if
            end if

            has_forces = present(forces)
            if (has_forces) then
                if (size(forces, 2) /= size(positions, 2) .or. size(forces, 1) /= 3) then
                    call istat%raise(FExceptionDescription("The forces array has the wrong size.", &
                        "When writing the extended CFG file " // bold(filename) // ":"))
                    exit body
                end if
            end if

            if (present(masses)) then
                if (size(masses) /= size(positions, 2)) then
                    call istat%raise(FExceptionDescription("The masses array has the wrong size.", &
                        "When writing the extended CFG file " // bold(filename) // ":"))
                    exit body
                end if
            end if

!------ Build the inclusion mask
            if (present(include)) then
                if (size(include) /= size(positions, 2)) then
                    call istat%raise(FExceptionDescription("The atomic mask has the wrong size.", &
                        "When writing the extended CFG file " // bold(filename) // ":"))
                    exit body
                end if
                allocate(mask, source=include)
            else
                allocate(mask(size(positions, 2)))
                mask = .true.
            end if
            atoms_count = count(mask)

            has_auxiliaries = present(aux)
            has_auxiliaries = .false.
            cols = 0
            if (has_forces) cols = cols + 3
            if (has_auxiliaries) cols = cols + aux%num_columns()
            allocate(headers(cols))
            n = 0
            if (has_forces) then
                headers(1) = "fx"
                headers(2) = "fy"
                headers(3) = "fz"
                n = 3
            end if
            if (has_auxiliaries) then
                do c=1, aux%num_columns()
                    n = n + 1
                    headers(n) = aux%column_name(c)
                end do
            end if
            cols = cols + 3
            if (has_velocities) cols = cols + 3
            allocate(tmpdat(cols))

!------ Print header
            call output_file%write_line("Number of particles = " // atoms_count)
            call output_file%write_line(get_mark())
            if (present(scale)) call output_file%write_line("A = " // scale)
            call output_file%write_line("H0(1,1) = " // box(1,1))
            call output_file%write_line("H0(1,2) = " // box(1,2))
            call output_file%write_line("H0(1,3) = " // box(1,3))
            call output_file%write_line("H0(2,1) = " // box(2,1))
            call output_file%write_line("H0(2,2) = " // box(2,2))
            call output_file%write_line("H0(2,3) = " // box(2,3))
            call output_file%write_line("H0(3,1) = " // box(3,1))
            call output_file%write_line("H0(3,2) = " // box(3,2))
            call output_file%write_line("H0(3,3) = " // box(3,3))
            if (.not.has_velocities) call output_file%write_line(".NO_VELOCITY.")
            call output_file%write_line("entry_count = " // cols)
            do c=1, size(headers)
                call output_file%write_line("auxiliary[" // (c-1) // "] = " // headers(c)%string())
            end do

!------ Special case: empty configuration
            if (atoms_count == 0) then
                call output_file%close
                exit body
            end if

!------ Print first atom data
            last_tag = "#INVALID"

            do i=1, size(positions, 2)
                if (.not.mask(i)) cycle

                ! Write the atom's element and mass if it is not the same as the previous atom
                if (tags(i) /= last_tag) then
                    if (present(masses)) then
                        call output_file%write_line("" // masses(i))
                    else
                        call output_file%write_line("0.0")
                    end if
                    call output_file%write_line(trim(tags(i)))
                    last_tag = tags(i)
                end if

                ! Get the atom's data
                tmpdat(1:3) = positions(:,i)
                c = 3
                if (has_velocities) then
                    tmpdat(c+1:c+3) = velocities(:,i)
                    c = c + 3
                end if
                if (has_forces) then
                    tmpdat(c+1:c+3) = forces(:,i)
                    c = c + 3
                end if
                if (has_auxiliaries) then
                    call aux%get_row(i, row)
                    tmpdat(c+1:) = row
                end if
                n = writeArray(buffer, tmpdat, int(size(tmpdat), c_size_t), len(buffer))
                call output_file%write_line(trim(buffer))
            end do

            call output_file%close
        end block body

        if (present(stat)) call stat%transfer(istat)
    end subroutine
end submodule
