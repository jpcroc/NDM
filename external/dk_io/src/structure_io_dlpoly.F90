!******************************************************************************!
!                    dk_structure_io_dlpoly submodule
!------------------------------------------------------------------------------!
!> Specific procedures to read and write DL_POLY structure files.
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
submodule (dk_structure_io) dk_structure_io_dlpoly
    use ext_character

    use dk_token,      only: Token, TokenError, get_token
    use dk_math,       only: matrix_inverse
    use dk_exception,  only: FException, FExceptionDescription

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
        start_ = 0
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
        if (start_ > 0) then
            context = " "//trim(line_number_buffer)//" | "//trim(line)//new_line("")// &
                empty(:1+line_number_len)//" | "//empty(:start_-1)//trim(caret)
        else
            context = ""
        end if
    end subroutine
!******************************************************************************!
!> Get a human-readable description of an error that happened while reading a
!> DL_POLY file.
!******************************************************************************!
    function DlpolyFileError(file, message, line, value)
        type(FileObject), intent(in) :: file
        character(*), intent(in) :: message
        character(*), intent(in), optional :: line
        character(*), intent(in), optional :: value
        type(FExceptionDescription) :: DlpolyFileError
!------
        character(:), allocatable :: context
!------
        if (present(line)) then
            if (present(value)) then
                call get_error_context(line, value, file%current_line, context)
            else
                call get_error_context(line, "", file%current_line, context)
            end if

            DlpolyFileError = FExceptionDescription(message, &
                "In DL_POLY file "//bold(file%full_name())//":"//file%current_line//":"//new_line("")//context)
        else
            DlpolyFileError = FExceptionDescription(message, &
                "In DL_POLY file "//bold(file%full_name())//":"//file%current_line//":")
        end if
    end function
!******************************************************************************!
!> Read a structure from a DL_POLY CONFIG or HISTORY file.
!******************************************************************************!
    module subroutine read_dlpoly(input_file, box, positions, tags, velocities, forces, info, stat)
        type(FileObject), intent(inout) :: input_file
        real(e), dimension(3,3), intent(out) :: box
        real(e), dimension(:,:), allocatable, intent(out) :: positions
        character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
        real(e), dimension(:,:), allocatable, intent(out), optional :: velocities
        real(e), dimension(:,:), allocatable, intent(out), optional :: forces
        type(FileInfo), intent(out), optional :: info
        type(FException), intent(out), optional :: stat
!------
        type(FException) :: istat
        character(:), allocatable :: line, config_keys_line, history_keys_line, timestep_line
        integer :: iostat, num_fields, pbc, num_atoms, num_frames, num_records, i, j, read_atoms
        integer :: n, start, end
        real(e), dimension(3,3) :: ibox
        type(Token) :: error_token
!------
        num_fields = -1
        pbc = -1
        num_atoms = -1
        num_frames = -1
        num_records = -1

        body: block

!------ Read file header

            ! The first line is a description
            call input_file%read_line(line, stat=istat)
            if (istat /= 0) exit body

            ! The second line has parameters
            call input_file%read_line(line, stat=istat)
            if (istat /= 0) exit body

            read(line,*,iostat=iostat) num_fields, pbc, num_atoms, num_frames, num_records
            if (iostat == 0) then
                history_keys_line = trim(line)
            else
                read(line,*,iostat=iostat) num_fields, pbc, num_atoms
                if (iostat /= 0) then
                    ! SYNTAX ERROR: the keys line is inconsistent with both the CONFIG and HISTORY formats
                    call istat%raise("Syntax error on record 2")
                    exit body
                end if
                config_keys_line = trim(line)
            end if

!------ Read frame data
            ! The first line is either the timestep record, or the simulation box
            call input_file%read_line(line, stat=istat)
            if (istat /= 0) exit body
            i = index(line, "timestep")
            if (i > 0) then

                ! Check the consistency of the file format
                if (num_frames == -1) then
                    ! SYNTAX ERROR: the keys line is inconsistent with both the CONFIG and HISTORY formats
                    call istat%raise("missing number of frames")
                    exit body
                end if
                if (num_records == -1) then
                    ! SYNTAX ERROR: the keys line is inconsistent with both the CONFIG and HISTORY formats
                    call istat%raise("missing number of records")
                    exit body
                end if

                ! Read the parameters for the current timestep
                timestep_line = trim(line)

                call input_file%read_line(line, stat=istat)
                if (istat /= 0) exit body
            end if

            ! Read the simulation box
            if (pbc > 0) then

                call read_array(line, box(:,1), n)
                if (n < 3) then
                    call get_token(line, n+1, error_token)
                    error_token%line_number = input_file%current_line
                    call istat%raise(TokenError("Invalid element for cell vector a.", &
                        "In DL_POLY file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                    exit body
                end if

                call input_file%read_line(line, stat=istat)
                if (istat /= 0) then
                    ! SYNTAX ERROR: the keys line is incosistent with both the CONFIG and HISTORY formats

                end if

                call read_array(line, box(:,2), n)
                if (n < 3) then
                    call get_token(line, n+1, error_token)
                    error_token%line_number = input_file%current_line
                    call istat%raise(TokenError("Invalid element for cell vector b.", &
                        "In DL_POLY file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                    exit body
                end if

                call input_file%read_line(line, stat=istat)
                if (istat /= 0) then
                    ! SYNTAX ERROR: the keys line is incosistent with both the CONFIG and HISTORY formats

                end if

                call read_array(line, box(:,3), n)
                if (n < 3) then
                    call get_token(line, n+1, error_token)
                    error_token%line_number = input_file%current_line
                    call istat%raise(TokenError("Invalid element for cell vector c.", &
                        "In DL_POLY file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                    exit body
                end if

                call input_file%read_line(line, stat=istat)
                if (istat /= 0) then
                    ! SYNTAX ERROR: the keys line is incosistent with both the CONFIG and HISTORY formats

                end if

                ibox = matrix_inverse(box)
            end if

!------ Allocate data arrays
            allocate(tags(num_atoms))
            tags = " "
            allocate(positions(3,num_atoms))
            positions = 0
            if (present(velocities)) then
                allocate(velocities(3,num_atoms))
                velocities = 0
            end if
            if (present(forces)) then
                allocate(forces(3,num_atoms))
                forces = 0
            end if

!------- Read atomic data
            read_atoms = 0
            do i=1, num_atoms

                ! Read the tag of the current atom
                start = 1
                do while(start < len(line) .and. is_whitespace(line(start:start)))
                    start = start + 1
                end do
                end = start + 1
                do while (end < len(line) .and. .not. is_whitespace(line(end:end)))
                    end = end + 1
                end do

                if (is_whitespace(line(start:end-1))) then
                    call istat%raise(DlpolyFileError(input_file, &
                        "Expected atom tag.", line, line(start:end-1)))
                    exit body
                end if

                ! Read the identifier of the current atom
                read(line(end:),*,iostat=iostat) j
                if (iostat /= 0) then
                    call istat%raise(DlpolyFileError(input_file, &
                        "Invalid atom identifier.", line, trim(adjustl(line(9:)))))
                    exit body
                end if
                tags(j) = line(start:end-1)
                read_atoms = read_atoms + 1

                call input_file%read_line(line, stat=istat)
                if (istat /= 0) then
                    call istat%discard
                    call istat%raise(DlpolyFileError(input_file, "Unexpected end of file."))
                    exit body
                end if

                ! Read the position of the current atom
                call read_array(line, positions(:,j), n)
                if (n < 3) then
                    call get_token(line, n+1, error_token)
                    error_token%line_number = input_file%current_line
                    call istat%raise(TokenError("Invalid atomic position.", &
                        "In DL_POLY file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                    exit body
                end if
                positions(:,j) = matmul(ibox, positions(:,j))

                ! Read the velocity of the current atom
                if (num_fields > 0) then
                    call input_file%read_line(line, stat=istat)
                    if (istat /= 0) then
                    call istat%discard
                        call istat%raise(DlpolyFileError(input_file, "Unexpected end of file."))
                        exit body
                    end if

                    if (present(velocities)) then
                        call read_array(line, velocities(:,j), n)
                        if (n < 3) then
                            call get_token(line, n+1, error_token)
                            error_token%line_number = input_file%current_line
                            call istat%raise(TokenError("Invalid atomic velocity.", &
                                "In DL_POLY file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                            exit body
                        end if
                    end if
                end if

                ! Read the force of the current atom
                if (num_fields > 1) then
                    call input_file%read_line(line, stat=istat)
                    if (istat /= 0) then
                    call istat%discard
                        call istat%raise(DlpolyFileError(input_file, "Unexpected end of file."))
                        exit body
                    end if

                    if (present(forces)) then
                        call read_array(line, forces(:,j), n)
                        if (n < 3) then
                            call get_token(line, n+1, error_token)
                            error_token%line_number = input_file%current_line
                            call istat%raise(TokenError("Invalid atomic force.", &
                                "In DL_POLY file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                            exit body
                        end if
                    end if
                end if

                if (i==num_atoms) exit

                ! Read the next line
                call input_file%read_line(line, stat=istat)
                if (istat /= 0) then
                    call istat%discard
                    call istat%raise(DlpolyFileError(input_file, "Unexpected end of file."))
                    exit body
                end if
            end do

            if (present(info)) then
                info%data_format = trim(STRUCTURE_DLPOLY)
                info%num_frames = 1
                info%compression_format = input_file%compression_format()
            end if
        end block body

        if (present(stat)) call stat%transfer(istat)
    end subroutine
!******************************************************************************!
!> Write a structure to a DL_POLY CONFIG file.
!******************************************************************************!
    module subroutine write_dlpoly(filename, box, positions, tags, include, elements, velocities, forces, stat)
        character(*), intent(in) :: filename
        real(e), dimension(3,3), intent(in) :: box
        real(e), dimension(:,:), intent(in) :: positions
        character(TAG_LENGTH), dimension(:), intent(in) :: tags
        logical, dimension(:), intent(in), optional :: include
        type(Element), dimension(:), intent(in), optional :: elements
        real(e), dimension(:,:), intent(in), optional :: velocities
        real(e), dimension(:,:), intent(in), optional :: forces
        type(FException), intent(out), optional :: stat
!------
        integer :: i, j, atoms_count, n, fields_count
        logical, dimension(:), allocatable :: mask
        character(10000) :: buffer
        real(c_double), dimension(:), allocatable :: tmpdat
        type(FileObject) :: output_file
        type(FException) :: istat
!------

        body: block

!------ Prepare the file
            call output_file%init(filename, "w", stat=istat)

            ! Return if the file could not be opened
            if (istat /= 0) then
                exit body
            end if

            if (size(tags) /= size(positions, 2)) then
                call istat%raise(FExceptionDescription("The tags array has the wrong size.", &
                    "When writing the DL_POLY file " // bold(filename) // ":"))
                exit body
            end if

            ! Set whether we need to save atomic velocities
            if (present(velocities)) then
                if (size(velocities, 2) /= size(positions, 2) .or. size(velocities, 1) /= 3) then
                    call istat%raise(FExceptionDescription("The velocities array has the wrong size.", &
                        "When writing the DL_POLY file " // bold(filename) // ":"))
                    exit body
                end if
            end if

            if (present(forces)) then
                if (size(forces, 2) /= size(positions, 2) .or. size(forces, 1) /= 3) then
                    call istat%raise(FExceptionDescription("The forces array has the wrong size.", &
                        "When writing the DL_POLY file " // bold(filename) // ":"))
                    exit body
                end if
            end if

!------ Build the inclusion mask
            if (present(include)) then
                if (size(include) /= size(positions, 2)) then
                    call istat%raise(FExceptionDescription("The atomic mask has the wrong size.", &
                        "When writing the DL_POLY file " // bold(filename) // ":"))
                    exit body
                end if
                allocate(mask, source=include)
                mask = include
            else
                allocate(mask(size(positions, 2)))
                mask = .true.
            end if
            atoms_count = count(mask)

            fields_count = 0
            if (present(velocities)) fields_count = 1
            if (present(forces)) fields_count = 2

!------ Print header
            call output_file%write_line(get_mark())
            call output_file%write_line(fields_count // " 3 " // atoms_count)
            call output_file%write_line(box(1,1) // " " // box(2,1) // " " // box(3,1))
            call output_file%write_line(box(1,2) // " " // box(2,2) // " " // box(3,2))
            call output_file%write_line(box(1,3) // " " // box(2,3) // " " // box(3,3))

!------ Special case: empty configuration
            if (atoms_count == 0) then
                call output_file%close
                exit body
            end if

            allocate(tmpdat(3))

!------ Print atomic data
            j = 0
            do i=1, size(positions, 2)
                if (.not.mask(i)) cycle

                j = j + 1
                call output_file%write_line(tags(i) // " " // j)

                tmpdat(:) = matmul(box, positions(:,i))
                n = writeArray(buffer, tmpdat, 3_c_size_t, len(buffer))
                call output_file%write_line(trim(buffer))

                if (fields_count > 0) then
                    if (present(velocities))  then
                        tmpdat(:) = velocities(:,i)
                    else
                        tmpdat(:) = 0.0
                    end if
                    n = writeArray(buffer, tmpdat, 3_c_size_t, len(buffer))
                    call output_file%write_line(trim(buffer))
                end if

                if (fields_count > 1) then
                    tmpdat(:) = forces(:,i)
                    n = writeArray(buffer, tmpdat, 3_c_size_t, len(buffer))
                    call output_file%write_line(trim(buffer))
                end if
            end do

            call output_file%close
        end block body

        if (present(stat)) call stat%transfer(istat)
    end subroutine
end submodule
