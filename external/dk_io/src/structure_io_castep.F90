!******************************************************************************!
!                    dk_structure_io_castep submodule
!------------------------------------------------------------------------------!
!> Specific procedures to read and write CASTEP structure files.
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
submodule (dk_structure_io) dk_structure_io_castep
    use dk_math,  only: matrix_inverse
    use dk_alloc, only: reallocate

    use dk_token,      only: Token, TokenError, get_token
    use dk_exception,  only: FException, FExceptionDescription, EVENT_LEVEL_LOG

    use ext_character

    implicit none(external, type)

contains
!******************************************************************************!
!> Read a structure from a CASTEP cell file.
!******************************************************************************!
    module subroutine read_castep(input_file, box, positions, tags, velocities, info, stat)
        type(FileObject), intent(inout) :: input_file
        real(e), dimension(3,3), intent(out) :: box
        real(e), dimension(:,:), allocatable, intent(out) :: positions
        character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
        real(e), dimension(:,:), allocatable, intent(out), optional :: velocities
        type(FileInfo), intent(out), optional :: info
        type(FException), intent(out), optional :: stat
!------
        type(FException) :: istat
        integer :: num_atoms, space_group_number, i
        real(e), dimension(6) :: parameters
        logical :: cartesian_mode, fractional_mode
        character(:), allocatable :: keyword
        character(:), allocatable :: line
        real(e), dimension(3,3) :: ibox
!------

        parameters = 0
        box = 0
        space_group_number = 0
        cartesian_mode = .false.
        fractional_mode = .false.
        num_atoms = 0

        body: block

            ! Process each line
            call input_file%read_nonempty_line(line, ["#", "!"], stat=istat)

            do while(istat == 0)

                ! Ignore comments
                i = index(line, "#")
                if (i > 0) line(i:) = " "

                ! Extract the keyword
                call get_keyword(line, keyword)

                if (allocated(keyword)) then

                    ! Proces each block depending on the keyword
                    select case(lower_case(keyword))

                    case ("lattice_abc")
                       call parse_cell(istat)

                    case ("lattice_cart")
                       call parse_vectors(istat)

                    case ("positions_frac")
                       fractional_mode = .true.
                       call parse_positions(istat)

                    case ("positions_abs")
                       cartesian_mode = .true.
                       call parse_positions(istat)

                    case ("ionic_velocities")
                       call parse_velocities(istat)

                    case default
                    end select
                end if
                if (istat /= 0) exit body

                call input_file%read_nonempty_line(line, ["#", "!"], stat=istat)
            end do

             if (all(box == 0)) then
                if (all(parameters == 0)) then
                    call istat%raise(FExceptionDescription("Missing lattice parameters.", &
                        "In CASTEP file "//bold(input_file%name)//":"))
                    exit body
                else
                    box = box_with_parameters(parameters)
                end if
             end if

             if (.not.allocated(positions)) then
                call istat%raise(FExceptionDescription("Missing atomic positions.", &
                    "In CASTEP file "//bold(input_file%name)//":"))
                exit body
             end if

            ! Convert positions to fractional if needed
             if (cartesian_mode) then
                ibox = matrix_inverse(box)
                do i=1, num_atoms
                    positions(:,i) = matmul(ibox, positions(:,i))
                end do
             end if

            if (present(info)) then
                info%data_format = trim(STRUCTURE_CASTEP)
                info%num_frames = 1
                info%compression_format = input_file%compression_format()
            end if
        end block body

        if (istat == "Trying to read past end of file.") call istat%discard

        if (present(stat)) call stat%transfer(istat)
    contains
        subroutine get_keyword(line, keyword)
            character(*), intent(in) :: line
            character(:), allocatable, intent(out) :: keyword
!------
            integer :: start, end
!------
            if (len_trim(line) == 0) return

            start = 1
            do while (start<len(line) .and. is_whitespace(line(start:start)))
                start = start + 1
            end do

            if (line(start:start) /= "%") return

            if (lower_case(line(start:start+6)) == "%block") then
                start = start+6
                do while(is_whitespace(line(start:start)))
                    start = start + 1
                end do
            end if

            end = start + 1
            do while (end < len(line) .and. .not.is_whitespace(line(end:end)))
                end = end + 1
            end do

            keyword = line(start:end-1)
        end subroutine

        subroutine parse_cell(pstat)
            type(FException), intent(out) :: pstat
!------
            integer :: n
            real(e) :: length_scale
            type(Token) :: error_token
!------
            call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            if (pstat /= 0) return

            select case(lower_case(trim(adjustl(line))))
            case ("bohr", "a0")
                length_scale = 0.529177_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case ("m")
                length_scale = 1.0e10_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case ("cm")
                length_scale = 1.0e8_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case ("nm")
                length_scale = 10.0_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case ("ang")
                length_scale = 1.0_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case default
                length_scale = 1.0_e
            end select
            if (pstat /= 0) return

            call read_array(line, parameters(:3), n)
            if (n < 3) then
                call get_token(line, n+1, error_token)
                error_token%line_number = input_file%current_line
                call pstat%raise(TokenError("Invalid lattice parameter value.", &
                    "In CASTEP file "//bold(input_file%name)//":", error_token))
                return
            end if
            parameters(:3) = parameters(:3) * length_scale

            call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            if (pstat /= 0) return
            call read_array(line, parameters(4:), n)
            if (n < 3) then
                call get_token(line, n+1, error_token)
                error_token%line_number = input_file%current_line
                call pstat%raise(TokenError("Invalid lattice angle value.", &
                    "In CASTEP file "//bold(input_file%name)//":", error_token))
                return
            end if

            ! Get the closing line for the lattice_abc block
            call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            if (pstat /= 0) return

            if (lower_case(line) /= "%endblock lattice_abc") then
                call get_token(line, 1, error_token)
                error_token%line_number = input_file%current_line
                call pstat%raise(TokenError("Expected end of LATTICE_ABC block.", &
                    "In CASTEP file "//bold(input_file%name)//":", error_token))
                return
            end if
        end subroutine

        subroutine parse_vectors(pstat)
            type(FException), intent(out) :: pstat
!------
            integer :: n
            type(Token) :: error_token
            real(e) :: length_scale
!------
            call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            if (pstat /= 0) return

            ! Ignore comments
            i = index(line, "#")
            if (i > 0) line(i:) = " "

            ! Get next line if the current line is a units definition
            select case(lower_case(trim(adjustl(line))))
            case ("bohr", "a0")
                length_scale = 0.529177_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case ("m")
                length_scale = 1.0e10_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case ("cm")
                length_scale = 1.0e8_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case ("nm")
                length_scale = 10.0_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case ("ang")
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
                length_scale = 1.0_e
            case default
                length_scale = 1.0_e
            end select
            if (pstat /= 0) return

            call read_array(line, box(:,1), n)
            if (n < 3) then
                call get_token(line, n+1, error_token)
                error_token%line_number = input_file%current_line
                call pstat%raise(TokenError("Invalid element for lattice vector a.", &
                    "In CASTEP file "//bold(input_file%name)//":", error_token))
                return
            end if

            call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            if (pstat /= 0) return

            ! Ignore comments
            i = index(line, "#")
            if (i > 0) line(i:) = " "

            call read_array(line, box(:,2), n)
            if (n < 3) then
                call get_token(line, n+1, error_token)
                error_token%line_number = input_file%current_line
                call pstat%raise(TokenError("Invalid element for lattice vector b.", &
                    "In CASTEP file "//bold(input_file%name)//":", error_token))
                return
            end if

            call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            if (pstat /= 0) return

            ! Ignore comments
            i = index(line, "#")
            if (i > 0) line(i:) = " "

            call read_array(line, box(:,3), n)
            if (n < 3) then
                call get_token(line, n+1, error_token)
                error_token%line_number = input_file%current_line
                call pstat%raise(TokenError("Invalid element for lattice vector c.", &
                    "In CASTEP file "//bold(input_file%name)//":", error_token))
                return
            end if

            ! Get the closing line for the lattice_cart block
            call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            if (pstat /= 0) return

            if (lower_case(line) /= "%endblock lattice_cart") then
                call get_token(line, 1, error_token)
                error_token%line_number = input_file%current_line
                call pstat%raise(TokenError("Expected end of LATTICE_ABC block.", &
                    "In CASTEP file "//bold(input_file%name)//":", error_token))
                return
            end if

            ! Rescale the box
            box(:,:) = box(:,:) * length_scale
        end subroutine

        subroutine parse_positions(pstat)
            type(FException), intent(out) :: pstat
!------
            integer :: n, start, end
!------
            allocate(tags(10))
            allocate(positions(3,10))

            call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)

            do while(pstat == 0)

                ! Get the atomic tag
                start = 1
                do while(start < len(line) .and. is_whitespace(line(start:start)))
                    start = start + 1
                end do
                end = start + 1
                do while (end < len(line) .and. .not. is_whitespace(line(end:end)))
                    end = end + 1
                end do

                ! Stop here if we reached the end of the block
                if (lower_case(line(start:end-1)) == "%endblock") then
                    ! We should check that the right kind of block is being ended
                    exit
                end if

                ! Add an atom
                num_atoms = num_atoms + 1

                ! Make sure the data arrays have enough space
                if (num_atoms>size(tags)) call reallocate(tags, 2*num_atoms)
                if (num_atoms>size(positions, 2)) call reallocate(positions, [3, 2*num_atoms])

                ! Set the current atom's tag
                tags(num_atoms) = line(start:end-1)

                ! Set the current atom's position
                call read_array(line(end:), positions(:,num_atoms), n)
                if (n < 3) then
                    call pstat%raise("Syntax error when reading atomic position")
                    return
                end if

                ! Read the next line
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            end do

            ! Resize the data arrays
            if (num_atoms<size(tags)) call reallocate(tags, num_atoms)
            if (num_atoms<size(positions, 2)) call reallocate(positions, [3, num_atoms])
        end subroutine

        subroutine parse_velocities(pstat)
            type(FException), intent(out) :: pstat
!------
            integer :: n, i, start, end
            type(Token) :: error_token
            real(e) :: scale
!------
            call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)

            allocate(velocities(3,10))
            i = 0

            ! Get next line if the current line is a units definition
            select case(lower_case(trim(line)))
            case ("auv")
                scale = 21876.912622_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case ("ang/ps")
                scale = 1.0_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case ("ang/fs")
                scale = 1.0e3_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case ("bohr/ps")
                scale = 0.529177211_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case ("bohr/fs")
                scale = 529.177211_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case ("m/s")
                scale = 1.0e-2_e
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            case default
                scale = 1.0_e
            end select
            if (pstat /= 0) return

            do while(pstat == 0)

                ! Get the atomic tag
                start = 1
                do while(start < len(line) .and. is_whitespace(line(start:start)))
                    start = start + 1
                end do
                end = start + 1
                do while (end < len(line) .and. .not. is_whitespace(line(end:end)))
                    end = end + 1
                end do

                ! Stop here if we reached the end of the block
                if (lower_case(line(start:end-1)) == "%endblock") then
                    ! We should check that the right kind of block is being ended
                    exit
                    return
                end if

                ! Make sure the velocities array is large enough
                i = i + 1
                if (i>size(velocities, 2)) then
                    call reallocate(velocities, [3, 2*i])
                end if

                ! Set the current atom's velocity
                call read_array(line, velocities(:,i), n)
                if (n < 3) then
                    call get_token(line, n+1, error_token)
                    error_token%line_number = input_file%current_line
                    call pstat%raise(TokenError("Invalid velocity element.", &
                        "In CASTEP file "//bold(input_file%name)//":", error_token))
                    return
                end if

                ! Read the next line
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            end do

            ! Resize the data array
            if (i<size(velocities, 2)) call reallocate(velocities, [3, i])
        end subroutine

        subroutine parse_symmetry(pstat)
            type(FException), intent(out) :: pstat
!------
!------
            call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)

            do while(pstat == 0)


                ! Read the next line
                call input_file%read_nonempty_line(line, ["#", "!"], stat=pstat)
            end do

        end subroutine
    end subroutine

    module subroutine write_castep(filename, box, positions, tags, masses, velocities, include, elements, stat)
        character(*), intent(in) :: filename
        real(e), dimension(3,3), intent(in) :: box
        real(e), dimension(:,:), intent(in) :: positions
        character(TAG_LENGTH), dimension(:), intent(in) :: tags
        real(e), dimension(:), intent(in), optional :: masses
            real(e), dimension(:,:), intent(in), optional :: velocities
        logical, dimension(:), intent(in), optional :: include
        type(Element), dimension(:), intent(in), optional :: elements
        type(FException), intent(out), optional :: stat
!------
        integer :: i, atoms_count, n
        logical, dimension(:), allocatable :: mask
        character(1000) :: buffer
        type(FileObject) :: output_file
        type(FException) :: istat
!------

        body: block

!------ Prepare the file
            call output_file%init(filename, "write", istat)

            ! Return if the file could not be opened
            if (istat /= 0) then
                exit body
            end if

            if (size(tags) /= size(positions, 2)) then
                call istat%raise(FExceptionDescription("The tags array has the wrong size.", &
                    "When writing the CASTEP file " // bold(filename) // ":"))
                exit body
            end if

!------ Build the inclusion mask
            if (present(include)) then
                if (size(include) /= size(positions, 2)) then
                    call istat%raise(FExceptionDescription("The atomic mask has the wrong size.", &
                        "When writing the CASTEP file " // bold(filename) // ":"))
                    exit body
                end if
                allocate(mask, source=include)
            else
                allocate(mask(size(positions, 2)))
                mask = .true.
            end if
            atoms_count = count(mask)

!------ Print header
            call output_file%write_line(get_mark())

!------ Write the unit cell
            call output_file%write_line("")
            call output_file%write_line("%BLOCK LATTICE_CART")
            call output_file%write_line("ang")
            call output_file%write_line("" // box(1,1) // " " // box(2,1) // " " // box(3,1))
            call output_file%write_line("" // box(1,2) // " " // box(2,2) // " " // box(3,2))
            call output_file%write_line("" // box(1,3) // " " // box(2,3) // " " // box(3,3))
            call output_file%write_line("%ENDBLOCK LATTICE_CART")

!------ Special case: empty configuration
            if (atoms_count == 0) then
                call output_file%close
                exit body
            end if

!------ Print atomic positions
            call output_file%write_line("")
            call output_file%write_line("%BLOCK POSITIONS_FRAC")
            do i=1, size(positions, 2)
                if (.not.mask(i)) cycle
                buffer = " "
                n = writeArray(buffer, positions(:,i), 3_c_size_t, len(buffer))
                call output_file%write_line(tags(i) // " " // trim(buffer))
            end do
            call output_file%write_line("%ENDBLOCK POSITIONS_FRAC")

!------ Print atomic velocities
            if (present(velocities)) then
                call output_file%write_line("")
                call output_file%write_line("%BLOCK IONIC_VELOCITIES")
                call output_file%write_line("ang/ps")
                do i=1, size(velocities, 2)
                    buffer = " "
                    n = writearray(buffer, velocities(:,i), 3_c_size_t, len(buffer))
                    call output_file%write_line(trim(buffer))
                end do
                call output_file%write_line("%ENDBLOCK IONIC_VELOCITIES")
            end if

            call output_file%close
        end block body

        if (present(stat)) call stat%transfer(istat)
    end subroutine
end submodule
