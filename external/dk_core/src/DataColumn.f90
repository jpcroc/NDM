module dk_datacolumn

    use dk_alloc,     only: reallocate
    use dk_string,    only: FString, reallocate
    use dk_exception, only: FException

    implicit none(external, type)
    private
    save

    integer, parameter, public :: dp = kind(0.0d0)
    integer, parameter :: BUFFER_LENGTH = 1000

    integer, parameter :: COLUMN_UNDEF = 0
    integer, parameter :: COLUMN_REAL = 1
    integer, parameter :: COLUMN_CHARACTER = 2

    type, public :: DataColumn
        private

        integer, public :: num_rows = 0
        character(:), allocatable :: name_
        character(:), allocatable :: unit_
        character(:), allocatable :: format_
        integer :: type = COLUMN_UNDEF
        integer :: string_length_ = -1
        real(dp), dimension(:), allocatable, public :: data
        type(FString), dimension(:), allocatable, public :: strings
        type(DataColumn), allocatable, public :: next
    contains
        procedure, public :: add_column
        procedure, public :: init

        procedure, public :: append_character
        procedure, public :: append_real
        procedure, public :: append_string
        generic, public :: append => append_character, append_real, append_string

        procedure, public :: name
        procedure, public :: unit
        procedure, public :: string
        procedure, public :: format
        procedure, public :: string_length
        procedure, public :: value
        procedure, public :: get_value_string
        procedure, public :: resize

        procedure, public :: print_summary
    end type

contains

    subroutine print_summary(this)
        class(DataColumn), intent(in) :: this !< The data column

        real(dp) :: average, std_deviation, sum, y, t, c, variance
        integer :: i

        write(*,*) this%name()
        if (allocated(this%data)) then

            write(*,*) "Minimum:", minval(this%data)
            write(*,*) "Maximum:", maxval(this%data)

            ! Kahan summation algorithm for the average
            sum = 0.0_dp
            c = 0.0_dp
            do i = 1, size(this%data)
                y = this%data(i) - c
                t = sum + y
                c = (t - sum) - y
                sum = t
            end do
            average = sum / size(this%data)
            write(*,*) "Average:", average

            ! Kahan summation algorithm for the variance
            sum = 0.0d0
            c = 0.0d0
            do i = 1, size(this%data)
                y = (this%data(i) - average)**2 - c
                t = sum + y
                c = (t - sum) - y
                sum = t
            end do
            variance = sum / size(this%data)
            std_deviation = sqrt(variance)
            write(*,*) "Standard deviation:", std_deviation
        else

        end if


    end subroutine
!******************************************************************************!
!> Setup a data column. If the size is specified, the data array is
!> pre-allocated, which can minimise re-allocations when adding frows to the
!> data column later,
!******************************************************************************!
    subroutine init(this, name, unit, length, format)
        class(DataColumn), intent(out) :: this !< The data column
        character(*), intent(in) :: name !< The name of the data column
        character(*), intent(in), optional :: unit !< The physical units of the data column
        integer, intent(in), optional :: length !< The size to allocate for the data column
        character(*), intent(in), optional :: format !< The format string to print data from the column
!------
        this%name_ = name
        if (present(unit)) this%unit_ = unit
        if (present(length)) allocate(this%data(length))
        if (present(format)) this%format_ = format
        this%string_length_ = string_length(this)
    end subroutine

    subroutine resize(this, length)
        class(DataColumn), intent(out) :: this  !< The data column
        integer, intent(in), optional :: length !< The size to allocate for the data column
!------
        real(dp), dimension(:), allocatable :: tmp
        type(FString), dimension(:), allocatable :: tmp_string
!------
        if (this%type == COLUMN_REAL) then
            if (.not.allocated(this%data)) then
                allocate(this%data(length))
            else if (size(this%data) < length) then
                allocate(tmp(length))
                tmp(:size(this%data)) = this%data
                deallocate(this%data)
                call move_alloc(from=tmp, to=this%data)
            end if
        else if (this%type == COLUMN_CHARACTER) then
            if (.not.allocated(this%strings)) then
                allocate(this%strings(length))
            else if (size(this%strings) < length) then
                allocate(tmp_string(length))
                tmp_string(:size(this%strings)) = this%strings
                deallocate(this%strings)
                call move_alloc(from=tmp_string, to=this%strings)
            end if
        end if
    end subroutine
!******************************************************************************!
!> Append a data point to a data column.
!******************************************************************************!
    subroutine append_real(this, data)
        class(DataColumn), intent(inout) :: this
        real(dp), intent(in) :: data
!------
        if (this%type == COLUMN_UNDEF) call set_type(this, COLUMN_REAL)

        if (this%type /= COLUMN_REAL) then

            return
        end if

        ! Allocate the data array if it was not initialised
        if (.not.allocated(this%data)) then
            allocate(this%data(100))
            this%data = 0
        end if

        ! Make some room in the data array if needed
        if (this%num_rows == size(this%data)) then
            call reallocate(this%data, size(this%data)*2)
        end if

        ! Set the value
        this%num_rows = this%num_rows + 1
        this%data(this%num_rows) = data
    end subroutine
!******************************************************************************!
!> Append a characer string to a data column.
!******************************************************************************!
    subroutine append_character(this, string)
        class(DataColumn), intent(inout) :: this
        character(*), intent(in) :: string
!------
        if (this%type == COLUMN_UNDEF) call set_type(this, COLUMN_CHARACTER)

        if (this%type /= COLUMN_CHARACTER) then

            return
        end if

        ! Allocate the data array if it was not initialised
        if (.not.allocated(this%strings)) then
            allocate(this%strings(100))
        end if

        ! Make some room in the data array if needed
        if (this%num_rows == size(this%strings)) then
            call reallocate(this%strings, size(this%strings)*2)
        end if

        ! Set the value
        this%num_rows = this%num_rows + 1
        this%strings(this%num_rows) = string
    end subroutine
!******************************************************************************!
!> Append a characer string to a data column.
!******************************************************************************!
    subroutine append_string(this, string)
        class(DataColumn), intent(inout) :: this
        type(FString), intent(in) :: string
!------
        if (this%type == COLUMN_UNDEF) call set_type(this, COLUMN_CHARACTER)

        if (this%type /= COLUMN_CHARACTER) then

            return
        end if

        ! Allocate the data array if it was not initialised
        if (.not.allocated(this%strings)) then
            allocate(this%strings(100))
        end if

        ! Make some room in the data array if needed
        if (this%num_rows == size(this%strings)) then
            call reallocate(this%strings, size(this%strings)*2)
        end if

        ! Set the value
        this%num_rows = this%num_rows + 1
        this%strings(this%num_rows) = string
    end subroutine
!******************************************************************************!
!> Add a column to a column list. Does nothing if there is already a column
!> with the same name in the chain list.
!******************************************************************************!
    recursive subroutine add_column(this, name, unit, length, format)
        class(DataColumn), intent(inout) :: this
        character(*), intent(in) :: name
        character(*), intent(in), optional :: unit
        integer, intent(in), optional :: length
        character(*), intent(in), optional :: format
!------
        if (this%name_ == name) return

        if (allocated(this%next)) then
            call this%next%add_column(name, unit, length, format=format)
        else
            allocate(this%next)
            call this%next%init(name, unit, length=length, format=format)
        end if
    end subroutine
!******************************************************************************!
!> Get the number of rows in the column.
!******************************************************************************!
    pure function column_length(this)
        class(DataColumn), intent(in) :: this
        integer :: column_length
!------
        column_length = this%num_rows
    end function
!******************************************************************************!
!> Get the name of the column.
!******************************************************************************!
    pure function name(this)
        class(DataColumn), intent(in) :: this
        character(:), allocatable :: name
!------
        if (allocated(this%name_)) name = this%name_
    end function
!******************************************************************************!
!> Get a string representation of the physical units of the data.
!******************************************************************************!
    pure function unit(this)
        class(DataColumn), intent(in) :: this
        character(:), allocatable :: unit
!------
        if (allocated(this%unit_)) unit = this%unit_
    end function
!******************************************************************************!
!> Get a string representation of the column header.
!******************************************************************************!
    pure function string(this)
        class(DataColumn), intent(in) :: this
        character(:), allocatable :: string
!------
        if (allocated(this%name_) .and. allocated(this%unit_)) then
            if (this%unit_ /= "") then
                string = this%name_ // " [" // this%unit_ // "]"
            else
                string = this%name_
            end if
        else if (allocated(this%name_)) then
            string = this%name_
        else if (allocated(this%unit_)) then
            string = "[" // this%unit_ // "]"
        end if
    end function
!******************************************************************************!
!> Get the format string to use for data output.
!******************************************************************************!
    pure function format(this)
        class(DataColumn), intent(in) :: this
        character(:), allocatable :: format
!------
        if (allocated(this%format_)) then
            format =  "(" // this%format_ // ")"
        else
            if (this%type == COLUMN_REAL) then
                format = "(es12.4e3)"
            else if (this%type == COLUMN_CHARACTER) then
                format = "(a12)"
            end if
        end if
    end function
!******************************************************************************!
!> Get the length of a numerical value printed with the column's format.
!******************************************************************************!
    function string_length(this)
        class(DataColumn), intent(in) :: this
        integer :: string_length
!------
        character(BUFFER_LENGTH) :: buffer
        integer :: i, value_length, header_length
        character(:), allocatable :: format
!------
        if (this%string_length_ > 0) then
            string_length = this%string_length_
        else
            ! Scan the format string to see whether we are printing an integer
            format = this%format()

            if (format /= "") then
                i = 1
                do while(format(i:i) == " " .or. format(i:i) == "(")
                    i = i + 1
                end do
                if (format(i:i) == "a") then
                    write(buffer,this%format()) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
                else if (format(i:i) /= "i") then
                    write(buffer,this%format()) -huge(1.0d0)
                else
                    write(buffer,this%format()) -huge(1)
                end if

                i = 1
                do while(buffer(i:i) == " ")
                    i = i + 1
                end do

                value_length = index(buffer(i:), " ") + i - 1

                header_length = len(this%string())

                string_length = max(value_length, header_length)
            else
                string_length = len(this%string())
            end if
        end if
    end function
!******************************************************************************!
!> Get the value for a specified row in a data column.
!******************************************************************************!
    function value(this, row)
        class(DataColumn), target, intent(in) :: this !< The data column
        integer, intent(in) :: row !< The row
        real(dp), pointer :: value
!------
        value => null()
        if (row <= this%num_rows .and. row > 0) value => this%data(row)
    end function
!******************************************************************************!
!> Get the value for a specified row in a data column.
!******************************************************************************!
    subroutine get_value_string(this, row, value_string)
        class(DataColumn), target, intent(in) :: this !< The data column
        integer, intent(in) :: row !< The row
        character(*), intent(out) :: value_string
!------
        if (row <= this%num_rows .and. row > 0) then
            select case(this%type)
            case (COLUMN_REAL)
                write(value_string,this%format()) this%data(row)
            case (COLUMN_CHARACTER)
                value_string(:) = this%strings(row)%string()
            end select
        end if
    end subroutine
!******************************************************************************!
!> Assign a value to a row of the data column.
!******************************************************************************!
    subroutine set_real(this, row, data)
        class(DataColumn), intent(inout) :: this
        integer, intent(in) :: row
        real(dp), intent(in) :: data
!------
        if (this%type == COLUMN_UNDEF) call set_type(this, COLUMN_REAL)

        if (this%type /= COLUMN_REAL) return

        ! Allocate the data array if it was not initialised
        if (.not.allocated(this%data)) then
            allocate(this%data(row))
            this%data = 0
        end if

        ! Make some room in the data array if needed
        if (row > size(this%data)) then
            call reallocate(this%data, max(row, size(this%data)*2))
        end if

        ! Set the value
        this%num_rows = max(row, this%num_rows)
        this%data(row) = data
    end subroutine
!******************************************************************************!
!> Assign a character value to a row of the data column.
!******************************************************************************!
    subroutine set_character(this, row, string)
        class(DataColumn), intent(inout) :: this
        integer, intent(in) :: row
        character(*), intent(in) :: string
!------
        if (this%type == COLUMN_UNDEF) call set_type(this, COLUMN_CHARACTER)

        if (this%type /= COLUMN_CHARACTER) return
        
        ! Allocate the data array if it was not initialised
        if (.not.allocated(this%strings)) then
            allocate(this%strings(row))
        end if

        ! Make some room in the data array if needed
        if (row > size(this%strings)) then
            call reallocate(this%strings, max(row, size(this%strings)*2))
        end if

        ! Set the value
        this%num_rows = max(row, this%num_rows)
        this%strings(row) = string
    end subroutine
!******************************************************************************!
!> Assign a string value to a row of the data column.
!******************************************************************************!
    subroutine set_string(this, row, string)
        class(DataColumn), intent(inout) :: this
        integer, intent(in) :: row
        type(FString), intent(in) :: string
!------
        if (this%type == COLUMN_UNDEF) call set_type(this, COLUMN_CHARACTER)

        if (this%type /= COLUMN_CHARACTER) return
        
        ! Allocate the data array if it was not initialised
        if (.not.allocated(this%strings)) then
            allocate(this%strings(row))
        end if

        ! Make some room in the data array if needed
        if (row > size(this%strings)) then
            call reallocate(this%strings, max(row, size(this%strings)*2))
        end if

        ! Set the value
        this%num_rows = max(row, this%num_rows)
        this%strings(row) = string
    end subroutine
!******************************************************************************!
!> Set the type of the data in a data column.
!******************************************************************************!
    subroutine set_type(this, type, stat)
        class(DataColumn), intent(inout) :: this
        integer, intent(in) :: type
        type(FException), intent(out), optional :: stat
!------
        type(FException) :: istat
!------
        if (this%type == COLUMN_UNDEF) then
            select case(type)
            case (COLUMN_CHARACTER)
                this%type = COLUMN_CHARACTER
                this%string_length_ = string_length(this)
            case(COLUMN_REAL)
                this%type = COLUMN_REAL
                this%string_length_ = string_length(this)
            case default
                call istat%raise("Trying to set the column to an invalid type.")
            end select
        else
            call istat%raise("The column already has a type.")
        end if

        if (present(stat)) call stat%transfer(istat)
    end subroutine
end module
