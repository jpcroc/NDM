!******************************************************************************!
!                    dk_structure_io_vasp submodule
!------------------------------------------------------------------------------!
!> Specific procedures to read and write VASP structure files.
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
submodule (dk_structure_io) dk_structure_io_vasp
    use dk_math,       only: matrix_inverse
    use dk_string,     only: FString
    use ext_character, only: operator(//), split, bold, is_whitespace

    use dk_token,      only: Token, TokenError, get_token
    use dk_exception,  only: FExceptionDescription, FException

    implicit none(external, type)

contains
!******************************************************************************!
!> Read a structure from a VASP POSCAR file.
!******************************************************************************!
    module subroutine read_vasp(input_file, box, positions, tags, velocities, info, stat)
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
        character(:), allocatable :: buffer
        real(e) :: x
        integer, dimension(:), allocatable :: composition
        logical :: cartesianMode, cartesian_velocities
        character(:), dimension(:), allocatable :: element_names
        type(Token) :: word, error_token
        type(FException) :: stat_
!------

        body: block

!------ Process the header

            ! First line can be ignored
            call input_file%read_line(buffer, stat=stat_)

            ! Second line contains the length scale
            call input_file%read_line(buffer, stat=stat_)
            call get_words(buffer, word, 1)
            if (buffer(word%start+len(word%value):) /= " ") then
                call stat_%raise(FExceptionDescription( &
                    "Expected end of line."))
                exit body
            end if
            read(word%value,*,iostat=iostat) x
            if (iostat /= 0) then
                call get_token(buffer, 1, error_token)
                error_token%line_number = input_file%current_line
                call stat_%raise(TokenError("Second line is not a valid length scale.", &
                    "In VASP file "//bold(input_file%name)//":", error_token))
                exit body
            end if

            ! Third line contains the first box vector
            call input_file%read_line(buffer, stat= stat_)
            call read_array(buffer, box(:,1), n)
            if (n < 3) then
                call get_token(buffer, n+1, error_token)
                error_token%line_number = input_file%current_line
                call stat_%raise(TokenError("Invalid element for lattice vector a.", &
                    "In VASP file "//bold(input_file%name)//":", error_token))
                return
            end if

            ! Fourth line contains the second box vector
            call input_file%read_line(buffer, stat= stat_)
            call read_array(buffer, box(:,2), n)
            if (n < 3) then
                call get_token(buffer, n+1, error_token)
                error_token%line_number = input_file%current_line
                call stat_%raise(TokenError("Invalid element for lattice vector b.", &
                    "In VASP file "//bold(input_file%name)//":", error_token))
                return
            end if

            ! Fifth line contains the third box vector
            call input_file%read_line(buffer, stat= stat_)
            call read_array(buffer, box(:,3), n)
            if (n < 3) then
                call get_token(buffer, n+1, error_token)
                error_token%line_number = input_file%current_line
                call stat_%raise(TokenError("Invalid element for lattice vector c.", &
                    "In VASP file "//bold(input_file%name)//":", error_token))
                return
            end if

            box = box * x

!------ Ion species section

            ! Sixth line is a list of species or the number of ions of each species
            call input_file%read_line(buffer, stat=stat_)
            if (stat_ /= 0) exit body

            ! Check whether the next line contains species names or population
            call split(trim(adjustl(buffer)), " ", element_names, group=.true.)
            allocate(composition(size(element_names)))

            ! Try to read the number of atoms
            iostat = 0
            do i=1, size(element_names)
                read(element_names(i),*,iostat=iostat) composition(i)
                if (iostat /= 0) exit
            end do

            ! If it failed, then the line should contain tags
            if (iostat /= 0) then
                call input_file%read_line(buffer, stat=stat_)
                read(buffer,*,iostat=iostat) composition
                if (iostat /= 0) then
                    ! Wrong ions per species line
                    call stat_%raise(FExceptionDescription( &
                        "Expected number of atoms per species."))
                    exit body
                end if
            end if
            atoms_count = sum(composition)
            allocate(positions(3,atoms_count))
            allocate(tags(atoms_count))

!------ Ion positions section

            ! Seventh line might or might not start with s
            call input_file%read_line(buffer, stat= stat_)
            buffer = adjustl(buffer)
            if (buffer(1:1) == "s" .or. buffer(1:1) == "S") then
                ! If it does, just ignore it
                call input_file%read_line(buffer, stat=stat_)
            end if

            ! Next line might or might not start with c or k
            buffer = adjustl(buffer)
            if (buffer(1:1) == "c" .or. buffer(1:1) == "C" .or. &
                buffer(1:1) == "k" .or. buffer(1:1) == "K") then

                ! If it does, then the atoms coordinates are cartesian
                cartesianMode = .true.
            else

                ! Otherwise, the atoms coordinates are fractional
                cartesianMode = .false.
            end if

            ! Process the atoms' descriptions
            iostat = 0
            do i=1, sum(composition)
                call input_file%read_line(buffer, stat=stat_)
                if (stat_ /= 0 .and. stat_ /= "End of file.") then
                    call stat_%raise(FExceptionDescription( &
                        "Error when reading atomic positions."))
                    exit body
                end if

                call read_array(buffer, positions(:,i), n)
                if (n < 3) then
                    call get_token(buffer, n+1, error_token)
                    error_token%line_number = input_file%current_line
                    call stat_%raise(TokenError("Invalid ionic position.", &
                        "In VASP file "//bold(input_file%name)//":", error_token))
                    return
                end if
            end do

            ! Set chemical species
            i = 0
            do a=1, size(composition)
                do n=1, composition(a)
                    i = i + 1
                    if (allocated(element_names)) then
                        tags(i) = element_names(a)
                    else
                        tags(i) = "" // a
                    end if
                end do
            end do

            ! Set fractional coordinates
            if (cartesianMode) then
                ibox = matrix_inverse(box)

                do i=1, atoms_count
                    positions(:,i) = matmul(ibox, positions(:,i))
                end do
            end if

            if (present(info)) then
                info%data_format = trim(STRUCTURE_VASP)
                info%num_frames = 1
                info%compression_format = input_file%compression_format()
            end if

!------ Lattice velocities section
            exit body
            ! See whether there might be more lines
            call input_file%read_line(buffer, stat=stat_)
            if (stat_ == "Trying to read past end of file.") then
                call stat_%discard
                exit body
            end if
            if (buffer(1:1) == "l" .or. buffer(1:1) == "L") then
                call input_file%read_line(buffer, stat=stat_)
                call input_file%read_line(buffer, stat=stat_)
                call input_file%read_line(buffer, stat=stat_)
                call input_file%read_line(buffer, stat=stat_)
                call input_file%read_line(buffer, stat=stat_)
                call input_file%read_line(buffer, stat=stat_)
                call input_file%read_line(buffer, stat=stat_)
            end if

!------ Ionic velocities section
            call input_file%read_line(buffer, stat=stat_)
            if (buffer(1:1) == "C" .or. buffer(1:1) == "c" .or. &
                buffer(1:1) == "K" .or. buffer(1:1) == "k" .or. buffer == " ") then

                cartesian_velocities = .true.
            else if (buffer(1:1) == "D" .or. buffer(1:1) == "d") then
                cartesian_velocities = .false.
            else
                error_token%value = trim(buffer)
                error_token%line = trim(buffer)
                error_token%start = 1
                error_token%line_number = input_file%current_line
                call stat_%raise(TokenError("Expected ionic velocities block.", &
                    "In VASP file "//bold(input_file%name)//":", error_token))
                exit body
            end if

            if (present(velocities)) allocate(velocities, mold=positions)

            ! Process the atoms' descriptions
            iostat = 0
            do i=1, sum(composition)
                call input_file%read_line(buffer, stat=stat_)
                if (stat_ /= 0 .and. stat_ /= "End of file.") then
                    call stat_%raise(FExceptionDescription( &
                        "Error when reading atomic velocities."))
                    exit body
                end if

                call read_array(buffer, velocities(:,i), n)
                if (n < 3) then
                    call get_token(buffer, n+1, error_token)
                    error_token%line_number = input_file%current_line
                    call stat_%raise(TokenError("Invalid ionic velocity.", &
                        "In VASP file "//bold(input_file%name)//":", error_token))
                    return
                end if
            end do

        end block body

        if (present(stat)) call stat%transfer(stat_)
    contains
        subroutine get_words(line, word, n)
            character(*), intent(in) :: line
            type(Token), intent(out) :: word
            integer, intent(in), optional :: n
!------
            integer :: i, m, read
!------
            ! Set the number of words to read
            if (present(n)) then
                m = n
            else
                m = 1
            end if
            read = 0

            ! Ignore leading whitespace
            i = 1
            do while (i<len(line))
                if (.not. is_whitespace(line(i:i))) exit
                i = i + 1
            end do
            word%start = i

            ! Scan to the next whitespace character or the end of the line
            do while (i<len(line))
                if (is_whitespace(line(i:i))) then
                    read = read + 1
                    if (read == m) exit
                end if
                i = i + 1
            end do

            ! Set the token
            word%value = line(word%start:i-1)
            word%line = trim(line)
        end subroutine
    end subroutine
!******************************************************************************!
!> Write a structure to a VASP POSCAR file.
!******************************************************************************!
    module subroutine write_vasp(filename, box, positions, tags, masses, include, elements, stat)
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
        integer :: a, i, num_elements
        character(TAG_LENGTH), dimension(:), allocatable :: element_names
        type(FException) :: istat
!------

        body: block

            if (size(tags) /= size(positions, 2)) then
                call istat%raise(FExceptionDescription("The tags array has the wrong size.", &
                    "When writing the VASP structure file " // bold(filename) // ":"))
                exit body
            end if

            call output_file%init(filename, "w", istat)

            ! Return if the file could not be opened
            if (istat /= 0) then
                exit body
            end if

            ! Write the comment line
            call output_file%write_line(get_mark())

            ! Write the length scale
            call output_file%write_line("1")

            ! Write the simulation box
            call output_file%write_line("" // box(1,1) // " " // box(1,2) // " " // box(1,3))
            call output_file%write_line("" // box(2,1) // " " // box(2,2) // " " // box(2,3))
            call output_file%write_line("" // box(3,1) // " " // box(3,2) // " " // box(3,3))

            ! Write the chemical elements
            allocate(element_names(size(tags)))
            num_elements = 0
            element_loop: do i=1, size(tags)
                do a=1, num_elements
                    if (tags(i) == element_names(a)) cycle element_loop
                end do
                num_elements = num_elements + 1
                element_names(num_elements) = tags(i)
            end do element_loop
            str = element_names(1)
            do a=2, num_elements
                str = str // " " // element_names(a)
            end do
            call output_file%write_line(str)

            ! Write the atoms population by element
            str = ""
            do a=1, num_elements
                str = str // " " // count(tags == element_names(a))
            end do
            call output_file%write_line(str)

            call output_file%write_line("Direct")

            ! Write the atoms positions
            do a=1, num_elements
                do i=1, size(positions, 2)
                    if (tags(i) == element_names(a)) then
                        call output_file%write_line( &
                            ""  // positions(1,i) // &
                            " " // positions(2,i) // &
                            " " // positions(3,i))
                    end if
                end do
            end do
            call output_file%close

        end block body

        if (present(stat)) call stat%transfer(istat)
    end subroutine
end submodule
