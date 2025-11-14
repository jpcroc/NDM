!******************************************************************************!
!                    dk_structure_io_xsf submodule
!------------------------------------------------------------------------------!
!> Specific procedures to read and write XCrySDen XSF structure files.
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
submodule (dk_structure_io) dk_structure_io_xsf
    use dk_math,       only: matrix_inverse
    use dk_string,     only: FString
    use ext_character, only: operator(//), split, bold, is_whitespace

    use dk_token,      only: Token, TokenError, get_token
    use dk_exception,  only: FExceptionDescription, FException

    implicit none(external, type)

contains
!******************************************************************************!
!> Read a structure from a XSF file.
!******************************************************************************!
    module subroutine read_xsf(input_file, box, positions, tags, velocities, info, stat)
        type(FileObject), intent(inout) :: input_file
        real(e), dimension(3,3), intent(out) :: box
        real(e), dimension(:,:), allocatable, intent(out) :: positions
        character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
        real(e), dimension(:,:), allocatable, intent(out), optional :: velocities
        type(FileInfo), intent(out), optional :: info
        type(FException), intent(out), optional :: stat
!------
        integer :: i, iostat, read, a, n, atoms_count
        real(e), dimension(3,3) :: ibox
        character(:), allocatable :: line
        real(e) :: x
        integer, dimension(:), allocatable :: composition
        character(:), dimension(:), allocatable :: element_names
        type(Token) :: word, error_token
        type(FException) :: istat
        logical :: start
        real(e), dimension(3,3) :: convbox, primbox
!------

        body: block

            start = .true.

            ! Process each line in the input file
            do
                call input_file%read_nonempty_line(line, ['#'], stat=istat)

                if (istat /= 0) exit

                ! Get the file format
                if (start) then
                    select case(lower_case(adjustl(line)))
                    case ("molecule", "polymer", "slab")
                        call get_token(line, 1, error_token)
                        error_token%line_number = input_file%current_line
                        call istat%raise(TokenError("Only CRYSTAL files are supported.", &
                            "In XSF file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                    case ("crystal")

                    case default
                        call get_token(line, 1, error_token)
                        error_token%line_number = input_file%current_line
                        call istat%raise(TokenError("Invalid keyword.", &
                            "In XSF file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                    end select

                    start = .false.
                    if (istat /= 0) exit body
                    cycle
                end if

                ! Parse blocks
                select case(lower_case(adjustl(line)))
                case ("primvec")
                    call parse_box(box, istat)
                case ("convvec")
                    call parse_box(convbox, istat)
                case ("primcoord")
                    call parse_coords(positions, istat)
                case ("convcoord")
                    call parse_coords(positions, istat)
                case default
                    call get_token(line, 1, error_token)
                    error_token%line_number = input_file%current_line
                    call istat%raise(TokenError("Invalid keyword.", &
                        "In XSF file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                end select

                if (istat /= 0) exit body
            end do

            ! Set fractional coordinates
            ibox = matrix_inverse(box)
            do i=1, atoms_count
                positions(:,i) = matmul(ibox, positions(:,i))
            end do

            if (present(info)) then
                info%data_format = trim(STRUCTURE_XSF)
                info%num_frames = 1
                info%compression_format = input_file%compression_format()
            end if
        end block body

        if (istat == "Trying to read past end of file.") then
            call istat%discard
        end if

        ! Cleanup
        if (present(stat)) call stat%transfer(istat)
    contains
        subroutine parse_box(H, pstat)
            real(e), dimension(3,3), intent(out) :: H
            type(FException), intent(out) :: pstat
!------
            integer :: n, pos
            type(Token) :: error_token
            character(1), dimension(*), parameter :: name = ["a", "b", "c"]
!------
            H = 0

            do i=1, 3

                ! Read the vector
                call input_file%read_nonempty_line(line, ['#'], stat=istat)
                call read_array(line, H(:,i), n, pos)

                ! Check whether the whole vector could be read
                if (n < 3) then
                    call get_token(line, 1, error_token, offset=pos)
                    error_token%line_number = input_file%current_line
                    call pstat%raise(TokenError("Invalid element for cell vector "//name(i)//".", &
                        "In XSF file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                    return
                end if

                ! Check whether the line contains something after the a vector
                if (trim(line(pos:)) /= " " .and. trim(line(pos:)) /= C_NULL_CHAR) then
                    call get_token(line, 1, error_token, offset=pos)
                    error_token%line_number = input_file%current_line
                    call pstat%raise(TokenError("Junk after cell vector "//name(i)//".", &
                        "In XSF file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                    return
                end if
            end do
        end subroutine

        subroutine parse_coords(coord, pstat)
            real(e), dimension(:,:), allocatable, intent(out) :: coord
            type(FException), intent(out) :: pstat
!------
            real(e), dimension(6) :: val
            integer :: i, n, m
!------
            call input_file%read_nonempty_line(line, ['#'], stat=istat)

            ! Read the number of atoms
            call read_array(line, val, n)
            if (n < 2) then
                call get_token(line, n+1, error_token)
                error_token%line_number = input_file%current_line
                call pstat%raise(TokenError("Syntax error.", &
                    "In XSF file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                return
            end if
            if (nint(val(1)) /= val(1) .or. val(1)<0) then
                call get_token(line, 1, error_token)
                error_token%line_number = input_file%current_line
                call pstat%raise(TokenError("Invalid number of atoms.", &
                    "In XSF file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                return
            end if
            atoms_count = nint(val(1))

            allocate(tags(atoms_count))
            allocate(coord(3,atoms_count))

            do i=1, atoms_count
                call input_file%read_line(line, stat=istat)
                if (istat == "Trying to read past end of file.") then
                    call istat%discard
                    call pstat%raise(FileError(input_file, "Missing atoms."))
                end if
                if (istat /= 0) return

                ! Get the atomic tag
                n=1
                do while(line(n:n)==" " .and. n<len(line))
                    n = n + 1
                end do
                m = n
                do while(line(m:m)/=" " .and. m<len(line))
                    m = m + 1
                end do
                tags(i) = line(n:m-1)

                call read_array(line(m:), val, n)
                if (n < 3) then
                    call get_token(line, n+2, error_token)
                    error_token%line_number = input_file%current_line
                    call pstat%raise(TokenError("Invalid atomic coordinate.", &
                        "In XSF file "//bold(input_file%name//":"//input_file%current_line)//":", error_token))
                    return
                end if
                coord(:,i) = val(:3)
            end do
        end subroutine
    end subroutine
!******************************************************************************!
!> Write a structure to a XSF file.
!******************************************************************************!
    module subroutine write_xsf(filename, box, positions, tags, masses, include, elements, stat)
        character(*), intent(in) :: filename
        real(e), dimension(3,3), intent(in) :: box
        real(e), dimension(:,:), intent(in) :: positions
        character(TAG_LENGTH), dimension(:), intent(in) :: tags
        real(e), dimension(:), intent(in), optional :: masses
        logical, dimension(:), intent(in), optional :: include
        type(Element), dimension(:), intent(in), optional :: elements
        type(FException), intent(out), optional :: stat
!------
        type(FileObject) :: output_file
        character(:), allocatable :: str
        integer :: a, i, num_elements, n, num_elements
        integer, dimension(:), allocatable :: label
        character(TAG_LENGTH), dimension(:), allocatable :: element_names
        type(FException) :: istat
        character(1000,c_char) :: buffer
        type(Element), dimension(:), allocatable :: elements_
        integer, dimension(:), allocatable :: element_id
!------

        body: block

            if (size(tags) /= size(positions, 2)) then
                call istat%raise(FExceptionDescription("The tags array has the wrong size.", &
                    "When writing the XSF structure file " // bold(filename) // ":"))
                exit body
            end if

            call output_file%init(filename, "w", istat)

            ! Return if the file could not be opened
            if (istat /= 0) then
                exit body
            end if

            ! Write the comment line
            call output_file%write_line(get_mark())

            call output_file%write_line("")
            call output_file%write_line("CRYSTAL")
            call output_file%write_line("")

!------ Write the shape matrix
            call output_file%write_line("PRIMVEC")
            n = writeArray(buffer, box(:,1), 3_c_size_t, len(buffer))
            call output_file%write_line(buffer(:n))
            n = writeArray(buffer, box(:,2), 3_c_size_t, len(buffer))
            call output_file%write_line(buffer(:n))
            n = writeArray(buffer, box(:,3), 3_c_size_t, len(buffer))
            call output_file%write_line(buffer(:n))
            call output_file%write_line("")

!------ Get the atoms labels
            allocate(element_id(size(tags)))
            element_id = -42
            allocate(elements_(size(positions, 2)))
            num_elements = 0
            atom_loop: do i=1, size(tags)
                do a=1, num_elements
                    if (tags(i) == elements_(a)%tag) then
                        element_id(i) = a
                        cycle atom_loop
                    end if
                end do
                num_elements = num_elements + 1
                elements_(num_elements)%tag = tags(i)
                element_id(i) = num_elements
            end do atom_loop

!------ Write the atomic positions
            call output_file%write_line("PRIMCOORD")
            call output_file%write_line(size(positions, 2) // " 1")
            do i=1, size(positions, 2)
                n = writeArray(buffer,  matmul(box, positions(:,i)), 3_c_size_t, len(buffer))
                call output_file%write_line(element_id(i) // " " // buffer(:n))
            end do

            call output_file%close

        end block body

        if (present(stat)) call stat%transfer(istat)
    end subroutine
end submodule
