!******************************************************************************!
!                    dk_structure_io_gulp submodule
!------------------------------------------------------------------------------!
!> Specific procedures to read and write GULP structure files.
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
submodule (dk_structure_io) dk_structure_io_gulp
    use dk_math,       only: matrix_inverse
    use dk_alloc,      only: reallocate

    use dk_token,      only: Token, TokenError, get_token
    use dk_exception,  only: FException, FExceptionDescription
    use dk_spacegroup, only: hall_number
    use ext_character

    implicit none(external, type)

contains
!******************************************************************************!
!> Read a structure from a GULP file.
!******************************************************************************!
    module subroutine read_gulp(input_file, box, positions, tags, info, stat)
        type(FileObject), intent(inout) :: input_file
        real(e), dimension(3,3), intent(out) :: box
        real(e), dimension(:,:), allocatable, intent(out) :: positions
        character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
        type(FileInfo), intent(out), optional :: info
        type(FException), intent(out), optional :: stat
!------
        type(FException) :: istat, pstat
        integer :: num_atoms, space_group_number, i
        real(e), dimension(6) :: parameters
        logical :: cartesian_mode, fractional_mode
        type(Token) :: keyword, space_group_token
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
            call input_file%read_nonempty_line(line, ["#"], stat=istat)

            do while(istat == 0 .and. pstat == 0)

                ! Ignore comments
                i = index(line, "#")
                if (i > 0) line(i:) = " "

                ! Combine continued lines

                ! Extract the keyword
                keyword%start = 0
                call get_next_word(line, keyword)

                ! Proces each block depending on the keyword
                select case(lower_case(keyword%value))
            
                case ("cell")
                    call parse_cell(istat)
                    if (istat /= 0) exit body
                    call input_file%read_nonempty_line(line, ["#"], stat=istat)

                case ("vec", "vectors")
                    call parse_vectors(istat)
                    if (istat /= 0) exit body
                    call input_file%read_nonempty_line(line, ["#"], stat=istat)

                case ("frac", "fractional")
                    fractional_mode = .true.
                    call parse_positions(istat)
                
                case ("cart", "cartesian")
                    cartesian_mode = .true.
                    call parse_positions(istat)

                case ("space", "spacegroup")
                    call parse_space_group(istat)
                    if (istat /= 0) exit
                    call input_file%read_nonempty_line(line, ["#"], stat=istat)

                case default
                    call input_file%read_nonempty_line(line, ["#"], stat=istat)
                end select

                if (istat == "Trying to read past end of file.") exit
            end do

            ! Report an error if the cell parameters are missing
            if (all(box == 0) .and. all(parameters == 0)) then
                call istat%raise(FExceptionDescription("Missing cell parameters.", &
                    "When reading the GULP file " // bold(input_file%name) // ":"))
            else if (all(box == 0)) then
                ! Otherwise, set the cell matrix
                box = box_with_parameters(parameters)
            end if

            ! Report an error if the atomic positions are missing
            if (.not.allocated(positions)) then
                call istat%raise(FExceptionDescription("Missing atomic positions.", &
                    "When reading the GULP file " // bold(input_file%name) // ":"))
            end if

            ! Convert positions to fractional if needed
             if (cartesian_mode) then
                ibox = matrix_inverse(box)
                do i=1, num_atoms
                    positions(:,i) = matmul(ibox, positions(:,i))
                end do
             end if

            ! Apply the space group
            if (allocated(space_group_token%value)) then
                    block
                        use dk_spacegroup, only: get_equivalent_positions_multiple, hall_number

                        ! Check that the space group is valid
                        if (hall_number(space_group_token%value) == 0) then
                            call istat%raise(TokenError("Unknown space group.", &
                                "In GULP file "//bold(input_file%name//":"//input_file%current_line)//":", space_group_token))
                            exit body
                        end if

                        call get_equivalent_positions_multiple(space_group_token%value, positions, tags)
                    end block
            end if

            if (present(info)) then
                info%data_format = trim(STRUCTURE_GULP)
                info%num_frames = 1
                info%compression_format = input_file%compression_format()
            end if
        end block body

        if (istat == "Trying to read past end of file.") call istat%discard

        if (present(stat)) call stat%transfer(istat)
    contains
        subroutine get_next_word(line, word)
            character(*), intent(in) :: line
            type(Token), intent(inout) :: word
!------
            integer :: i
!------
            ! Check whether we are starting
            if (word%start <= 0) then
                i = 1
            else
                i = word%start + len(word%value) + 1
            end if

            ! Ignore leading whitespace
            do while (i<len(line))
                if (.not. is_whitespace(line(i:i))) exit
                i = i + 1
            end do
            word%start = i

            ! Scan to the next whitespace character or the end of the line
            do while (i<len(line))
                if (is_whitespace(line(i:i))) exit
                i = i + 1
            end do

            ! Set the token
            word%value = line(word%start:i-1)
            word%line = trim(line)
            word%line_number = input_file%current_line
        end subroutine
!******************************************************************************!
!> Parse a `vectors` block in a GULP file.
!******************************************************************************!
        subroutine parse_vectors(pstat)
            type(FException), intent(out) :: pstat
!------
            integer :: i, n
            type(Token) :: error_token
!------            
            call input_file%read_nonempty_line(line, ["#"], stat=pstat)
            if (pstat /= 0) return

            ! Ignore comments
            i = index(line, "#")
            if (i > 0) line(i:) = " "

            call read_array(line, box(:,1), n)
            if (n < 3) then
                call get_token(line, n+1, error_token)
                error_token%line_number = input_file%current_line
                call pstat%raise(TokenError("Invalid element for cell vector a.", &
                    "In GULP file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                return
            end if

            call input_file%read_nonempty_line(line, ["#"], stat=pstat)
            if (pstat /= 0) return

            ! Ignore comments
            i = index(line, "#")
            if (i > 0) line(i:) = " "

            call read_array(line, box(:,2), n)
            if (n < 3) then
                call get_token(line, n+1, error_token)
                error_token%line_number = input_file%current_line
                call pstat%raise(TokenError("Invalid element for cell vector b.", &
                    "In GULP file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                return
            end if

            call input_file%read_nonempty_line(line, ["#"], stat=pstat)
            if (pstat /= 0) return

            ! Ignore comments
            i = index(line, "#")
            if (i > 0) line(i:) = " "

            call read_array(line, box(:,3), n)
            if (n < 3) then
                call get_token(line, n+1, error_token)
                error_token%line_number = input_file%current_line
                call pstat%raise(TokenError("Invalid element for cell vector c.", &
                    "In GULP file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                return
            end if
        end subroutine
!******************************************************************************!
!> Parse a `cell` block in a GULP file.
!******************************************************************************!
        subroutine parse_cell(pstat)
            type(FException), intent(out) :: pstat
!------
            type(Token) :: error_token
            integer :: i, n
!------
            call input_file%read_nonempty_line(line, ["#"], stat=pstat)
            if (pstat /= 0) return

            ! Ignore comments
            i = index(line, "#")
            if (i > 0) line(i:) = " "

            call read_array(line, parameters, n)
            if (n < 6) then
                call get_token(line, n+1, error_token)
                error_token%line_number = input_file%current_line
                call pstat%raise(TokenError("Invalid lattice parameter value.", &
                    "In GULP file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
            end if
        end subroutine
!******************************************************************************!
!> Parse a `fractional` or `cartesian` block in a GULP file.
!******************************************************************************!
        subroutine parse_positions(pstat)
            type(FException), intent(out) :: pstat
!------
            type(Token) :: word
            integer :: i, n, iostat, j
            real(e) :: a, b
!------
            allocate(tags(10))
            allocate(positions(3,10))

            call input_file%read_line(line, stat=pstat)
            i = 0
            parse: do while(pstat == 0)

                if (line == "") exit

                i = i + 1

                ! Make sure the data arrays have enough space
                if (i>size(tags)) call reallocate(tags, 2*i)
                if (i>size(positions, 2)) call reallocate(positions, [3, 2*i])
                
                ! Get the atomic tag of the current atom
                word%start = -1
                call get_next_word(line, word)
                tags(i) = word%value

                ! See whether the particle is a core or a shell
                call get_next_word(line, word)
                if (lower_case(word%value)=="c" .or. lower_case(word%value)=="core") then
                    call get_next_word(line, word)
                else if (lower_case(word%value)=="s" .or. lower_case(word%value)=="shell") then
                    tags(i) = trim(tags(i)) // "_Shell"
                    call get_next_word(line, word)
                end if

                ! Get the position
                do n=1, 3
                    j = index(word%value, "/")
                    if (j == 0) then
                        ! Try to read the parameter
                        call word%read_real(positions(n,i), iostat)

                        ! If it failed, it probably means we left the block of atomic positions
                        if (iostat /= 0) then
                            ! The first line *must* be an atomic position, otherwise no need to raise an exception
                            if (i == 1) then
                                call pstat%raise(TokenError("Invalid atomic coordinate.", &
                                    "In GULP file "//bold(input_file%name//":"//input_file%current_line)//":", word))
                            end if
                            exit parse
                        end if
                    else
                        read(word%value(:j-1),*,iostat=iostat) a
                        if (iostat /= 0) then
                            call pstat%raise(TokenError("Invalid atomic coordinate.", &
                                "In GULP file "//bold(input_file%name//":"//input_file%current_line)//":", word))
                            exit parse
                        end if

                        read(word%value(j+1:),*,iostat=iostat) b
                        if (iostat /= 0) then
                            call pstat%raise(TokenError("Invalid atomic coordinate.", &
                                "In GULP file "//bold(input_file%name//":"//input_file%current_line)//":", word))
                            exit parse
                        end if

                        positions(n,i) = a/b
                        
                    end if

                    ! Get the next token
                    call get_next_word(line, word)
                end do

                ! Proceed to the next line
                call input_file%read_line(line, stat=pstat)
            end do parse

            ! Resize the data arrays
            num_atoms = i
            call reallocate(tags, num_atoms)
            call reallocate(positions, [3,num_atoms])
        end subroutine
!******************************************************************************!
!> Parse a `space` block in a GULP file.
!******************************************************************************!
        subroutine parse_space_group(pstat)
            type(FException), intent(out) :: pstat
!------
            call input_file%read_nonempty_line(line, ["#"], stat=pstat)
            if (pstat /= 0) return

            space_group_token%value = trim(line)
            space_group_token%line = trim(line)
            space_group_token%line_number = input_file%current_line
        end subroutine
    end subroutine
!******************************************************************************!
!> Write a structure to a GULP file.
!******************************************************************************!
    module subroutine write_gulp(filename, box, positions, tags, masses, include, elements, stat)
        character(*), intent(in) :: filename
        real(e), dimension(3,3), intent(in) :: box
        real(e), dimension(:,:), intent(in) :: positions
        character(TAG_LENGTH), dimension(:), intent(in) :: tags
        real(e), dimension(:), intent(in), optional :: masses
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
                    "When writing the GULP file " // bold(filename) // ":"))
                exit body
            end if

!------ Build the inclusion mask
            if (present(include)) then
                if (size(include) /= size(positions, 2)) then
                    call istat%raise(FExceptionDescription("The atomic mask has the wrong size.", &
                        "When writing the GULP file " // bold(filename) // ":"))
                    exit body
                end if
                allocate(mask, source=include)
            else
                allocate(mask(size(positions, 2)))
                mask = .true.
            end if
            atoms_count = count(mask)

!------ Print header
            call output_file%write_line("single")
            call output_file%write_line(get_mark())

            call output_file%write_line("")
            call output_file%write_line("vectors")
            call output_file%write_line("" // box(1,1) // " " // box(2,1) // " " // box(3,1))
            call output_file%write_line("" // box(1,2) // " " // box(2,2) // " " // box(3,2))
            call output_file%write_line("" // box(1,3) // " " // box(2,3) // " " // box(3,3))

!------ Special case: empty configuration
            if (atoms_count == 0) then
                call output_file%close
                exit body
            end if

!------ Print first atom data
            call output_file%write_line("")
            call output_file%write_line("fractional")
            do i=1, size(positions, 2)
                if (.not.mask(i)) cycle
                n = writeArray(buffer, positions(:,i), 3_c_size_t, len(buffer))
                call output_file%write_line(tags(i) // " core " // trim(buffer))
            end do

            call output_file%close
        end block body

        if (present(stat)) call stat%transfer(istat)
    end subroutine
end submodule
