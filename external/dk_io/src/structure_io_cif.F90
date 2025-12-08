!******************************************************************************!
!                    dk_structure_io_cif submodule
!------------------------------------------------------------------------------!
!> Specific procedures to read and write CIF structure files.
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
submodule (dk_structure_io) dk_structure_io_cif
    use dk_alloc,      only: reallocate
    use dk_string,     only: FString
    use dk_datatable,  only: DataTable
    use dk_datacolumn, only: DataColumn

    use dk_exception,  only: FException, FExceptionDescription
    use dk_token, only: Token, TokenError
    use ext_character

    implicit none(external, type)

    integer, parameter :: TOKEN_UNDEF = 0
    integer, parameter :: TOKEN_DATA_BLOCK  = 1
    integer, parameter :: TOKEN_DATA_NAME = 2
    integer, parameter :: TOKEN_DATA_VALUE = 3
    integer, parameter :: TOKEN_KEYWORD = 4
    integer, parameter :: TOKEN_CIF_VERSION = 5
    integer, parameter :: TOKEN_WORD = 6
    integer, parameter :: TOKEN_END = 7

contains
!******************************************************************************!
!> Write a structure to a CIF file.
!******************************************************************************!
    module subroutine write_cif(filename, box, positions, tags, include, scale, stat)
        character(*), intent(in) :: filename
        real(e), dimension(3,3), intent(in) :: box
        real(e), dimension(:,:), intent(in) :: positions
        character(TAG_LENGTH), dimension(:), intent(in) :: tags
        logical, dimension(:), intent(in), optional :: include
        real(e), intent(in), optional :: scale
        type(FException), intent(out), optional :: stat
!------
        logical, dimension(:), allocatable :: mask
        type(FException) :: istat
        integer :: i
        real(e), dimension(6) :: parameters
        type(FileObject) :: output_file
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
                    "When writing the CIF file " // bold(filename) // ":"))
                exit body
            end if

!------ Build the inclusion mask
            if (present(include)) then
                if (size(include) /= size(positions, 2)) then
                    call istat%raise(FExceptionDescription("The atomic mask has the wrong size.", &
                        "When writing the CIF file " // bold(filename) // ":"))
                    exit body
                end if
                allocate(mask, source=include)
            else
                allocate(mask(size(positions, 2)))
                mask = .true.
            end if

            parameters = parameters_with_box(box)

            call output_file%write_line("#\#CIF_1.1")
            call output_file%write_line("data_")
            call output_file%write_line("_symmetry_int_tables_number 1")
            call output_file%write_line("_symmetry_space_group_name_H-M 'P 1'")
            call output_file%write_line("_cell_length_a " // parameters(1))
            call output_file%write_line("_cell_length_b " // parameters(2))
            call output_file%write_line("_cell_length_c " // parameters(3))
            call output_file%write_line("_cell_angle_alpha " // parameters(4))
            call output_file%write_line("_cell_angle_beta " // parameters(5))
            call output_file%write_line("_cell_angle_gamma " // parameters(6))
            call output_file%write_line("loop_")
            call output_file%write_line(" _atom_site_label")
            call output_file%write_line(" _atom_site_fract_x")
            call output_file%write_line(" _atom_site_fract_y")
            call output_file%write_line(" _atom_site_fract_z")
            do i=1, size(tags)
                if (.not.mask(i)) cycle
                call output_file%write_line(" " // trim(tags(i)) // (i-1) // " " // &
                    positions(1,i) // " " // positions(2,i) // " " // positions(3,i))
            end do

        end block body
        call output_file%close

        if (present(stat)) call stat%transfer(istat)
    end subroutine
!******************************************************************************!
!> Read a structure from a CIF file.
!******************************************************************************!
    module subroutine read_cif(input_file, box, positions, tags, info, stat)
        type(FileObject), intent(inout) :: input_file
        real(e), dimension(3,3), intent(out) :: box
        real(e), dimension(:,:), allocatable, intent(out) :: positions
        character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
        type(FileInfo), intent(out), optional :: info
        type(FException), intent(out), optional :: stat
!------
        type(FException) :: istat
        integer :: cif_version, space_group_number
        character(:), allocatable :: space_group
        type(Token), target :: tokens
        type(Token), pointer :: current_token, space_group_token
        real(e), dimension(6) :: parameters
!------

        parameters = -1
        box = 0
        space_group_number = 0

        body: block
            call get_tokens(input_file, tokens, istat)
            if (istat /= 0) exit body

            current_token => tokens
            do while(associated(current_token))

                call parse_document(current_token, istat)
                if (istat /= 0) exit body

                if (allocated(current_token%next)) then
                    current_token => current_token%next
                else
                    current_token => null()
                end if
            end do

            ! Check whether any information is missing
            if (parameters(1) < 0) then
                call istat%raise(FExceptionDescription("Missing definition for "//bold("_cell_length_a")//".", &
                    "In CIF file " // bold(input_file%name) // ":"))
                exit body
            end if
            if (parameters(2) < 0) then
                call istat%raise(FExceptionDescription("Missing definition for "//bold("_cell_length_b")//".", &
                    "In CIF file " // bold(input_file%name) // ":"))
                exit body
            end if
            if (parameters(3) < 0) then
                call istat%raise(FExceptionDescription("Missing definition for "//bold("_cell_length_c")//".", &
                    "In CIF file " // bold(input_file%name) // ":"))
                exit body
            end if
            if (parameters(4) < 0) then
                call istat%raise(FExceptionDescription("Missing definition for "//bold("_cell_angle_alpha")//".", &
                    "In CIF file " // bold(input_file%name) // ":"))
                exit body
            end if
            if (parameters(5) < 0) then
                call istat%raise(FExceptionDescription("Missing definition for "//bold("_cell_angle_beta")//".", &
                    "In CIF file " // bold(input_file%name) // ":"))
                exit body
            end if
            if (parameters(6) < 0) then
                call istat%raise(FExceptionDescription("Missing definition for "//bold("_cell_angle_gamma")//".", &
                    "In CIF file " // bold(input_file%name) // ":"))
                exit body
            end if
            if (.not.allocated(tags)) then
                call istat%raise(FExceptionDescription("Missing atomic tags.", "In CIF file " // bold(input_file%name) // ":"))
                exit body
            end if
            if (.not.allocated(positions)) then
                call istat%raise(FExceptionDescription("Missing atomic positions.", "In CIF file " // bold(input_file%name) // ":"))
                exit body
            end if

            box = box_with_parameters(parameters)

            if (allocated(space_group)) then
#ifdef HAVE_SPGLIB
                if (space_group /= "P 1") then
                    block
                        use dk_spacegroup, only: get_equivalent_positions_multiple, hall_number

                        ! Check that the space group is valid
                        if (hall_number(space_group) == 0) then
                            call istat%raise(TokenError("Unknown space group.", &
                                "In CIF file "//bold(input_file%name)//":", space_group_token))
                            exit body
                        end if

                        call get_equivalent_positions_multiple(space_group, positions, tags)
                    end block
                end if
#else
                if (space_group /= "P 1") then
                            call istat%raise(TokenError("Only the P1 space group is supported " // &
                                "because spglib was not enabled.", &
                                "In CIF file "//bold(input_file%name)//":", space_group_token))
                    exit body
                end if
#endif
            end if

            if (present(info)) then
                info%data_format = trim(STRUCTURE_CIF)
                info%num_frames = 1
                info%compression_format = input_file%compression_format()
            end if
        end block body

        if (present(stat)) call stat%transfer(istat)
    contains
        subroutine parse_document(current_token, stat)
            type(Token), intent(inout), pointer :: current_token
            type(FException), intent(out) :: stat
!------
            if (current_token%type == TOKEN_CIF_VERSION) then
                select case(current_token%value)
                case ("#\#CIF_2.0")
                    cif_version = 20
                case ("#\#CIF_1.1")
                    cif_version = 11
                case ("#\#CIF_1.0")
                    cif_version = 10
                case default
                    call stat%raise(TokenError("Invalid CIF version number.", &
                        "In CIF file "//bold(input_file%name)//":", current_token))
                    return
                end select

                if (.not.allocated(current_token%next)) return
                current_token => current_token%next
            end if

            if (len(current_token%value) < 5) then
                call stat%raise(TokenError("Expected data_ block.", &
                    "In CIF file "//bold(input_file%name)//":", current_token))
            else
                if (current_token%value(:5) == "data_") then
                    if (.not.allocated(current_token%next)) return
                    current_token => current_token%next
                    call parse_data_block(current_token, stat)
                else
                    call stat%raise(TokenError("Expected data_ block.", &
                        "In CIF file "//bold(input_file%name)//":", current_token))
                end if
            end if
        end subroutine

        recursive subroutine parse_data_block(current_token, stat)
            type(Token), intent(inout), pointer :: current_token
            type(FException), intent(out) :: stat
!------
            do while(associated(current_token))

                if (current_token%type /= TOKEN_WORD) then
                    call istat%raise(TokenError("Expected loop_ or data_ block, or property name.", &
                        "In CIF file "//bold(input_file%name)//":", current_token))
                end if

                if (current_token%type == TOKEN_WORD) then
                    if (current_token%value(1:1) == "_") then
                        call parse_property(current_token, stat)
                    else if (current_token%value(:) == "loop_") then
                        call parse_loop(current_token, stat)
                    else if (current_token%value(:5) == "data_") then
                        if (.not.allocated(current_token%next)) return
                        current_token => current_token%next
                        call parse_data_block(current_token, stat)
                    else
                        call istat%raise(TokenError("Invalid property name.", &
                            "In CIF file "//bold(input_file%name)//":", current_token))
                    end if
                end if

                if (stat /= 0) return
                if (.not.allocated(current_token%next)) return
                current_token => current_token%next
            end do
        end subroutine

        subroutine parse_property(current_token, stat)
            type(Token), intent(inout), pointer :: current_token
            type(FException), intent(out) :: stat
!------
            character(:), allocatable :: name, value
            integer :: iostat
!------
            name = current_token%value
            if (.not.allocated(current_token%next)) then
                call stat%raise(TokenError("Missing property value.", &
                    "In CIF file "//bold(input_file%name)//":", current_token))
                return
            end if
            if (current_token%next%type == TOKEN_END) then
                call stat%raise(TokenError("Missing property value.", &
                    "In CIF file "//bold(input_file%name)//":", current_token))
                return
            end if

            current_token => current_token%next

            value = current_token%value

            select case(name)
            case ("_cell_length_a")
                read(value,*,iostat=iostat) parameters(1)
                if (iostat > 0) then
                    call stat%raise(TokenError("Invalid value for _cell_length_a.", &
                        "In CIF file "//bold(input_file%name)//":", current_token))
                    return
                end if
            case ("_cell_length_b")
                read(value,*,iostat=iostat) parameters(2)
                if (iostat > 0) then
                    call stat%raise(TokenError("Invalid value for _cell_length_b.", &
                        "In CIF file "//bold(input_file%name)//":", current_token))
                    return
                end if
            case ("_cell_length_c")
                read(value,*,iostat=iostat) parameters(3)
                if (iostat > 0) then
                    call stat%raise(TokenError("Invalid value for _cell_length_c.", &
                        "In CIF file "//bold(input_file%name)//":", current_token))
                    return
                end if
            case ("_cell_angle_alpha")
                read(value,*,iostat=iostat) parameters(4)
                if (iostat > 0) then
                    call stat%raise(TokenError("Invalid value for _cell_angle_alpha.", &
                        "In CIF file "//bold(input_file%name)//":", current_token))
                    return
                end if
            case ("_cell_angle_beta")
                read(value,*,iostat=iostat) parameters(5)
                if (iostat > 0) then
                    call stat%raise(TokenError("Invalid value for _cell_angle_beta.", &
                        "In CIF file "//bold(input_file%name)//":", current_token))
                    return
                end if
            case ("_cell_angle_gamma")
                read(value,*,iostat=iostat) parameters(6)
                if (iostat > 0) then
                    call stat%raise(TokenError("Invalid value for _cell_angle_gamma.", &
                        "In CIF file "//bold(input_file%name)//":", current_token))
                    return
                end if
            case ("_symmetry_int_tables_number")
                read(value,*,iostat=iostat) space_group_number
                if (iostat > 0) then
                    call stat%raise(TokenError("Invalid value for _symmetry_int_tables_number.", &
                        "In CIF file "//bold(input_file%name)//":", current_token))
                    return
                end if
            case ("_symmetry_space_group_name_H-M")
                space_group = value
                space_group_token => current_token
            end select
        end subroutine

        subroutine parse_loop(current_token, stat)
            type(Token), intent(inout), pointer :: current_token
            type(FException), intent(out) :: stat
!------
            type(DataTable) :: table
            integer :: row, column, iostat, i
            type(DataColumn), pointer :: tag_column, x_column, y_column, z_column
            real(e) :: x
!------
            if (.not.allocated(current_token%next)) return
            current_token => current_token%next

            ! Read the properties
            do while(associated(current_token))
                if (current_token%type == TOKEN_WORD .and. current_token%value(1:1) == "_") then
                    ! Add the property to the list of loop properties
                    call table%add_column(current_token%value)
                else
                    exit
                end if
                if (.not.allocated(current_token%next)) return
                current_token => current_token%next
            end do

            ! Read the data
            row = 0
            row_loop: do
                row = row + 1
                do column=1, table%num_columns()

                    if (current_token%type == TOKEN_END) exit row_loop

                    read(current_token%value,*,iostat=iostat) x
                    if (iostat == 0) then
                        call table%append(column, x)
                    else
                        call table%append(column, current_token%value)
                    end if

                    if (.not.allocated(current_token%next)) then
                        if (column == 1) then
                            exit row_loop
                        else
                            call stat%raise("Syntax error: missing loop data")
                            return
                        end if
                    end if

                    ! Get out of the loop if the next token is a property name
                    if (allocated(current_token%next)) then
                        if (allocated(current_token%next%value) .and. current_token%next%type==TOKEN_WORD) then
                            if (current_token%next%value(1:1) == "_") exit row_loop
                            if (current_token%next%value == "loop_") exit row_loop
                        end if
                    end if

                    ! Otherwise, just process the next token
                    current_token => current_token%next
                end do
            end do row_loop

            ! Set atomic positions if the loop has atoms
            tag_column => table%column("_atom_site_label")
            x_column => table%column("_atom_site_fract_x")
            y_column => table%column("_atom_site_fract_y")
            z_column => table%column("_atom_site_fract_z")
            if (associated(tag_column) .and. associated(x_column) .and. associated(y_column) .and. associated(z_column)) then
                allocate(tags(table%num_rows()))
                allocate(positions(3,table%num_rows()))
                do row=1, table%num_rows()
                    tags(row) = tag_column%strings(row)%string()
                    do i=len_trim(tags(row)), 1, -1
                        if (iachar(tags(row)(i:i)) <= iachar('9') .and. iachar(tags(row)(i:i)) >= iachar('0')) tags(row)(i:i) = " "
                    end do
                    positions(1,row) = x_column%value(row)
                    positions(2,row) = y_column%value(row)
                    positions(3,row) = z_column%value(row)
                end do
            end if
        end subroutine
    end subroutine
!******************************************************************************!
!> Tokenise a CIF file.
!******************************************************************************!
    subroutine get_tokens(input_file, tokens, stat)
        type(FileObject), intent(inout) :: input_file
        type(Token), intent(out), target :: tokens
        type(FException), intent(out) :: stat
!------
        type(FException) :: istat
        character(:), allocatable :: line
        integer :: i, j
        type(Token), pointer :: current_token
!------
        body: block

            call input_file%read_line(line, stat=istat)
            current_token => tokens

            ! The first line of the document might be a CIF format version number
            if (len(line) >= 7) then
                if (line(1:7) == "#\#CIF_") then
                    current_token%type = TOKEN_CIF_VERSION
                    current_token%value = trim(line)
                    current_token%line = trim(line)
                    current_token%line_number = input_file%line_number()
                    current_token%start = 1
                    allocate(current_token%next)
                    current_token => current_token%next

                    call input_file%read_line(line, stat=istat)
                end if
            end if

            do while(istat == 0)

                i = 0
                line_loop: do while(i<len(line))
                    i = i + 1

                    ! Process the next character if the current one is whitespace
                    if (is_whitespace(line(i:i))) cycle line_loop

                    if (line(i:i) == "#") exit line_loop ! Process the next line if we reached a comment

                    if (line(i:i) == "'") then ! Process single-quoted strings
                        do j=i+1, len(line)
                            if (line(j:j) == "'") then
                                current_token%type = TOKEN_DATA_VALUE
                                current_token%value = line(i+1:j-1)
                                current_token%line = trim(line)
                                current_token%line_number = input_file%line_number()
                                current_token%start = i+1
                                allocate(current_token%next)
                                current_token => current_token%next

                                i = j
                                cycle line_loop
                            end if
                        end do
                        current_token%value = line(i:)
                        current_token%line = trim(line)
                        current_token%line_number = input_file%line_number()
                        current_token%start = i
                        call stat%raise(TokenError("Unterminated string.", &
                            "In CIF file "//bold(input_file%name)//":", current_token))
                        return
                    end if

                    if (line(i:i) == """") then ! Process double-quoted strings
                        do j=i+1, len(line)
                            if (line(j:j) == """") then
                                current_token%type = TOKEN_DATA_VALUE
                                current_token%value = line(i+1:j-1)
                                current_token%line = trim(line)
                                current_token%line_number = input_file%line_number()
                                current_token%start = i+1
                                allocate(current_token%next)
                                current_token => current_token%next

                                i = j
                                cycle line_loop
                            end if
                        end do
                        current_token%value = line(i:)
                        current_token%line = trim(line)
                        current_token%line_number = input_file%line_number()
                        current_token%start = i
                        call stat%raise(TokenError("Unterminated string.", &
                            "In CIF file "//bold(input_file%name)//":", current_token))
                        return
                    end if

                    if (line(i:i) == ";") then ! Process multi-line values
                        current_token%type = TOKEN_DATA_VALUE
                        current_token%value = trim(line(i+1:))
                        current_token%line = trim(line)
                        current_token%line_number = input_file%line_number()
                        current_token%start = i

                        call input_file%read_line(line, stat=istat)
                        do while(istat == 0)
                            if (line == ";") then
                                allocate(current_token%next)
                                current_token => current_token%next

                                cycle line_loop
                            end if

                            current_token%value = current_token%value // new_line("") // trim(line)

                            call input_file%read_line(line, stat=istat)
                        end do

                        if (istat == "Trying to read past end of file.") then
                            ! We reached the end of file without finding the end of the value
                            call istat%discard

                            i = index(current_token%value, new_line(""))
                            if (i > 0) then
                                current_token%value = current_token%value(:i)
                            end if
                            call istat%raise(TokenError("Unterminated multi-line value.", &
                                "In CIF file "//bold(input_file%name)//":", current_token))
                        end if

                        ! An unknown error happened
                        call stat%transfer(istat)
                        return
                    end if

                    ! Process barewords
                    do j=i+1, len(line)
                        if (is_whitespace(line(j:j))) then
                            current_token%type = TOKEN_WORD ! 6
                            current_token%value = line(i:j-1)
                            current_token%line = trim(line)
                            current_token%line_number = input_file%line_number()
                            current_token%start = i
                            allocate(current_token%next)
                            current_token => current_token%next

                            i = j
                            cycle line_loop
                        end if
                    end do

                    current_token%type = TOKEN_WORD ! 6
                    current_token%value = line(i:)
                    current_token%line = trim(line)
                    current_token%line_number = input_file%line_number()
                    current_token%start = 1
                    allocate(current_token%next)
                    current_token => current_token%next

                    exit line_loop
                end do line_loop

                call input_file%read_line(line, stat=istat)
            end do
        end block body

        if (istat == "Trying to read past end of file.") call istat%discard
        if (istat == "Unexpected end of file.") call istat%discard

        call stat%transfer(istat)
    end subroutine
end submodule
