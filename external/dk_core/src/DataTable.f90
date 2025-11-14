module dk_datatable

    use dk_datacolumn, only: DataColumn, dp
    use ext_character
    use dk_string, only: FString

    implicit none(external, type)
    private
    save

    type, public :: DataTable
        private
        character(:), allocatable :: separator
        integer :: length = 0
        type(DataColumn), allocatable :: data
    contains
        procedure, public :: get_header !< Get a header with the column names

        ! Output
        procedure, public :: print_header !< Print a header with the column names
        procedure, public :: print_row    !< Print the data corresponding to a given row
        procedure, public :: print        !< Print the whole table

        !
        procedure, public :: num_rows     !< Get the number of rows in the data table
        procedure, public :: num_columns  !< Get the number of columns in the data table

        !
        procedure, public :: add_column   !< Add a column to the data table

        ! Append data
        procedure, public :: append_real_name       !< Append a value to a named column of the data table
        procedure, public :: append_real_index      !< Append a value to the column of the data table at a given index
        procedure, public :: append_character_name  !< Append a value to a named column of the data table
        procedure, public :: append_character_index !< Append a value to the column of the data table at a given index
        procedure, public :: append_string_name     !< Append a value to a named column of the data table
        procedure, public :: append_string_index    !< Append a value to the column of the data table at a given index
        generic, public :: append => append_real_name, append_real_index, &
            append_character_name, append_character_index, &
            append_string_name, append_string_index    !< Append a value to a column of the data table

        ! Getter functions
        procedure, public :: column_with_name  !< Get a pointer to a column object with a given name
        procedure, public :: column_with_index !< Get a pointer to the column object at a given position in the data table
        generic, public :: column => column_with_name, column_with_index !< Get a pointer to a column in the data table
        procedure, public :: column_index !< Get the index of a column

        procedure, public :: column_name

        procedure, public :: value_name   !< Get the value corresponding to given row and column name
        procedure, public :: value_index  !< Get the value corresponding to given row and column index
        generic, public :: value => value_name, value_index !< Get the value corresponding to given row and column

        procedure, public :: get_row      !< Get the values in a row
    end type

contains
!******************************************************************************!
!> Return a numerical value from a data table corresponding to given row and 
!> column name.
!******************************************************************************!
    function value_name(this, name, row) result(value)
        class(DataTable), target, intent(in) :: this
        character(*), intent(in) :: name
        integer, intent(in), optional :: row
        real(dp), pointer :: value
!------
        type(DataColumn), pointer :: current
!------
        current => this%data
        value => null()

        do while(associated(current))
            if (current%name() == name) then
                if (present(row)) then
                    value => current%value(row)
                else
                    value => current%value(current%num_rows)
                end if
                return
            end if
            current => current%next
        end do
    end function
!******************************************************************************!
!> Return a numerical value from a data table corresponding to given row and 
!> column index.
!******************************************************************************!
    function value_index(this, index, row) result(value)
        class(DataTable), target, intent(in) :: this
        integer, intent(in) :: index
        integer, intent(in), optional :: row
        real(dp), pointer :: value
!------
        type(DataColumn), pointer :: current
        integer :: i
!------
        current => this%data
        value => null()
        i = 1

        do while(associated(current))
            if (i == index) then
                if (present(row)) then
                    value => current%value(row)
                else
                    value => current%value(current%num_rows)
                end if
                return
            end if
            current => current%next
            i = i + 1
        end do
    end function
!******************************************************************************!
!> Get the index of a given column.
!******************************************************************************!
    function column_index(this, name)
        class(DataTable), target, intent(in) :: this
        character(*), intent(in) :: name
        integer :: column_index
!------
        type(DataColumn), pointer :: current
!------
        current => this%data
        column_index = 1

        do while(associated(current))
            if (current%name() == name) return
            current => current%next
            column_index = column_index + 1
        end do
        column_index = 0
    end function
!******************************************************************************!
!> Return a pointer to a column object.
!******************************************************************************!
    function column_with_name(this, name) result(column)
        class(DataTable), target, intent(in) :: this
        character(*), intent(in) :: name
        type(DataColumn), pointer :: column
!------
        column => this%data

        do while(associated(column))
            if (column%name() == name) return
            column => column%next
        end do
    end function
!******************************************************************************!
!> Return a pointer to a column object.
!******************************************************************************!
    function column_with_index(this, index) result(column)
        class(DataTable), target, intent(in) :: this
        integer, intent(in) :: index
        type(DataColumn), pointer :: column
!------
        integer :: i
!------
        column => this%data
        i = 1

        do while(associated(column))
            if (i == index) return
            column => column%next
            i = i + 1
        end do
    end function
!******************************************************************************!
!> Return the name of a column.
!******************************************************************************!
    function column_name(this, index) result(name)
        class(DataTable), target, intent(in) :: this
        integer, intent(in) :: index
        character(:), allocatable :: name
!------
        integer :: i
        type(DataColumn), pointer :: column
!------
        column => this%data
        i = 1

        do while(associated(column))
            if (i == index) name = column%name()
            column => column%next
            i = i + 1
        end do
    end function
!******************************************************************************!
!> Append a value to one of the columns of a data table.
!******************************************************************************!
    subroutine append_real_name(this, name, value)
        class(DataTable), target, intent(inout) :: this
        character(*), intent(in) :: name
        real(dp), intent(in) :: value
!------
        type(DataColumn), pointer :: c
!------
        c => column_with_name(this, name)

        if (associated(c)) call c%append(value)
    end subroutine
!******************************************************************************!
!> Append a value to one of the columns of a data table.
!******************************************************************************!
    subroutine append_real_index(this, index, value)
        class(DataTable), target, intent(inout) :: this
        integer, intent(in) :: index
        real(dp), intent(in) :: value
!------
        type(DataColumn), pointer :: c
!------
        c => column_with_index(this, index)

        if (associated(c)) call c%append(value)
    end subroutine
!******************************************************************************!
!> Append a value to one of the columns of a data table.
!******************************************************************************!
    subroutine append_character_name(this, name, value)
        class(DataTable), target, intent(inout) :: this
        character(*), intent(in) :: name
        character(*), intent(in) :: value
!------
        type(DataColumn), pointer :: c
!------
        c => column_with_name(this, name)

        if (associated(c)) call c%append(value)
    end subroutine
!******************************************************************************!
!> Append a value to one of the columns of a data table.
!******************************************************************************!
    subroutine append_character_index(this, index, value)
        class(DataTable), target, intent(inout) :: this
        integer, intent(in) :: index
        character(*), intent(in) :: value
!------
        type(DataColumn), pointer :: c
!------
        c => column_with_index(this, index)

        if (associated(c)) call c%append(value)
    end subroutine
!******************************************************************************!
!> Append a value to one of the columns of a data table.
!******************************************************************************!
    subroutine append_string_name(this, name, value)
        class(DataTable), target, intent(inout) :: this
        character(*), intent(in) :: name
        type(FString), intent(in) :: value
!------
        type(DataColumn), pointer :: c
!------
        c => column_with_name(this, name)

        if (associated(c)) call c%append(value)
    end subroutine
!******************************************************************************!
!> Append a value to one of the columns of a data table.
!******************************************************************************!
    subroutine append_string_index(this, index, value)
        class(DataTable), target, intent(inout) :: this
        integer, intent(in) :: index
        type(FString), intent(in) :: value
!------
        type(DataColumn), pointer :: c
!------
        c => column_with_index(this, index)

        if (associated(c)) call c%append(value)
    end subroutine

    subroutine get_header(this, header)
        use iso_c_binding, only: c_char, C_NULL_CHAR, c_size_t, c_sizeof
        class(DataTable), target, intent(in) :: this
        character(:), allocatable, intent(out) :: header
!------
        character(2000) :: buffer
        type(DataColumn), pointer :: current
        integer :: offset, increment, i, j
!------
        if (.not.allocated(this%data)) return

        buffer = ""

        offset = 2
        current => this%data
        increment = current%string_length()
        i = offset + (current%string_length()-apparentLength(current%string()))/2
        j = offset+increment-1 - (current%string_length()-apparentLength(current%string()))/2
        buffer(i:j) = current%string()
        offset = offset + increment + 3

        do while(allocated(current%next))
            current => current%next

            increment = current%string_length()
            i = offset + (current%string_length()-apparentLength(current%string()))/2
            j = offset+increment-1 - (current%string_length()-apparentLength(current%string()))/2
            buffer(i:j) = current%string()
            offset = offset + increment + 3
        end do
       buffer(offset:) = ""

       header = trim(buffer)

    end subroutine
!******************************************************************************!
!> Print a header with the names and units of the columns in a data table.
!******************************************************************************!
    subroutine print_header(this, show_line)
        use iso_c_binding, only: c_char, C_NULL_CHAR, c_size_t, c_sizeof
        class(DataTable), target, intent(in) :: this
        logical, intent(in), optional :: show_line !< Indicates whether some space needs to be left to align the columns with a table printed with row numbers
!------
        character(:), allocatable :: header
        logical :: show_line_
!------

       call this%get_header(header)

       if (.not.allocated(header)) return

        if (present(show_line)) then
            show_line_ = show_line
        else
            show_line_ = .false.
        end if

        if (show_line_) then
            write(*,'('//ceiling(log10(max(this%num_rows(),2)*1.0_dp))//'x,1x,a)') header
        else
            write(*,'(a)') header
        end if
    end subroutine
!******************************************************************************!
!> Print a line with the data from a given row of a data table. Print the last
!> row if no row is given.
!******************************************************************************!
    subroutine print_row(this, row, show_line)
        class(DataTable), target, intent(in) :: this !< The data table
        integer, intent(in), optional :: row !< The row to print (the last row is printed if not present)
        logical, intent(in), optional :: show_line !< Indicates whether the row number should be printed
!------
        character(2000) :: buffer
        type(DataColumn), pointer :: current
        integer :: offset, row_, c
        character(:), allocatable :: format
        logical :: show_line_
!------
        buffer = ""

        if (present(row)) then
            row_ = row
        else
            row_ = this%num_rows()
        end if

        if (present(show_line)) then
            show_line_ = show_line
        else
            show_line_ = .false.
        end if

        ! Prepare the string with the data
        offset = 2
        current => this%data
        c = 0
        do while(associated(current))
            c = c + 1

            ! Scan the format string to see whether we are printing an integer
            format = current%format()
            if (associated(current%value(row_))) then
                call current%get_value_string(row_, buffer(offset:offset+current%string_length()-1))
            end if

            ! Move the cursor and get the next column
            offset = offset + current%string_length() + 3
            current => current%next
        end do
        buffer(offset:) = ""

        ! Write the string
        if (show_line_) then
            write(*,'(i'//ceiling(log10(this%num_rows()*1.0_dp))//',1x,a)') row_, trim(buffer)
        else
            write(*,'(a)') trim(buffer)
        end if
    end subroutine
!******************************************************************************!
!> Print a representation of a data table, including header and each row.
!******************************************************************************!
    subroutine print(this)
        class(DataTable), target, intent(in) :: this
!------
        integer :: i
!------
        call this%print_header(show_line=.true.)
        do i=1, this%num_rows()
            call this%print_row(i, show_line=.true.)
        end do
    end subroutine
!******************************************************************************!
!> Get a copy of the data in a row of a data table.
!******************************************************************************!
    subroutine get_row(this, row, data)
        class(DataTable), target, intent(in) :: this !< The data table
        real(dp), dimension(:), allocatable, intent(inout) :: data !< The data array containing the row
        integer, intent(in), optional :: row !< The row to extract
!------
        type(DataColumn), pointer :: current
        integer :: i, row_
!------
        ! Use the last row if no explicit row number is given
        if (present(row)) then
            row_ = row
        else
            row_ = this%num_rows()
        end if

        ! Stop here if the data table does not have enough rows
        if (row_ < 1 .or. row_ > this%num_rows()) return

        ! Make sure the data array is allocated with the right size
        if (allocated(data)) then
            if (size(data) /= this%num_rows()) deallocate(data)
        end if
        if (.not.allocated(data)) allocate(data(this%num_rows()))

        ! Prepare the string with the data
        i = 0
        current => this%data
        do while(associated(current))
            i = i + 1
            data(i) = current%value(row_)
            current => current%next
        end do
    end subroutine

    subroutine set_row(this, row, data)
        class(DataTable), target, intent(inout) :: this
        integer, intent(in) :: row
        real(dp), dimension(:), intent(in) :: data
!------
        type(DataColumn), pointer :: column
        integer :: i
!------
        column => this%data
        do i=1, size(data)
            if (.not.associated(column)) then
                ! Error: mismatch between the number of columns and the size of the data array
                return
            end if

            ! Make sure the column has enough rows
            call column%resize(row)

            ! Set the value
            column%value(row) = data(i)
            column => column%next
        end do
    end subroutine
!******************************************************************************!
!> Get the number of rows in a data table. Returns the size of the longest
!> column if not all columns in the data table have the same size.
!******************************************************************************!
    function num_rows(this)
        class(DataTable), target, intent(in) :: this !< The data table
        integer :: num_rows !< The number of rows
!------
        type(DataColumn), pointer :: column
!------
        column => this%data
        num_rows = 0

        do while(associated(column))
            num_rows = max(num_rows, column%num_rows)
            column => column%next
        end do
    end function
!******************************************************************************!
!> Get the number of columns in a data table.
!******************************************************************************!
    function num_columns(this)
        class(DataTable), target, intent(in) :: this !< THe data table
        integer :: num_columns !< The number of columns
!------
        type(DataColumn), pointer :: column
!------
        num_columns = 0

        if (.not.allocated(this%data)) return
        
        column => this%data

        do while(associated(column))
            num_columns = num_columns + 1
            column => column%next
        end do
    end function
!******************************************************************************!
!> Add a column to a data table. No column is added if there is already a
!> column with the same combination of name and unit.
!******************************************************************************!
    subroutine add_column(this, name, unit, length, format)
        class(DataTable), intent(inout) :: this
        character(*), intent(in) :: name
        character(*), intent(in), optional :: unit
        integer, intent(in), optional :: length
        character(*), intent(in), optional :: format
!------
        if (allocated(this%data)) then
            call this%data%add_column(name, unit=unit, length=length, format=format)
        else
            allocate(this%data)
            call this%data%init(name, unit=unit, length=length, format=format)
        end if
    end subroutine
end module
