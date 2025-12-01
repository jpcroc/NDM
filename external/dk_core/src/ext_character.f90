!******************************************************************************!
!                           ext_character module
!------------------------------------------------------------------------------!
!> Operators for character variables.
!
!  Part of the Core library version 2.0.
!  Written by Paul Fossati, <paul.fossati@gmail.com>
!  Copyright (c) 2009-2018 Paul Fossati
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
module ext_character
    use iso_c_binding

    use core, only: e, sp, dp, no, yes

    implicit none (external, type)
    private
    save

!******************************************************************************!
!> Overloaded assignment operator to create character representation of
!> numerical variables
!******************************************************************************!
    interface assignment(=)
        module procedure characterWithInteger, characterWithIntegerArray, &
            characterWithSingle, characterWithSingleArray, &
            characterWithDouble, characterWithDoubleArray, &
            characterWithLogical
    end interface
    public :: assignment(=)

!******************************************************************************!
!> Overloaded concatenation for character and numerical variables
!******************************************************************************!
    interface operator(//)
        module procedure character_appendInteger, integer_appendCharacter, &
            character_appendSingle,  single_appendCharacter, &
            character_appendDouble,  double_appendCharacter, &
            character_appendLogical, logical_appendCharacter, &
            character_appendIntegerArray, integerArray_appendCharacter, &
            character_appendSingleArray, singleArray_appendCharacter, &
            character_appendDoubleArray, doubleArray_appendCharacter
    end interface
    public :: operator(//)

    ! ANSI colour codes
    character(*), parameter, public :: ANSI_COLOUR_DEFAULT = achar(27)//"[0m"
    character(*), parameter, public :: ANSI_COLOUR_RED     = achar(27)//"[31m"
    character(*), parameter, public :: ANSI_COLOUR_GREEN   = achar(27)//"[32m"
    character(*), parameter, public :: ANSI_COLOUR_YELLOW  = achar(27)//"[33m"
    character(*), parameter, public :: ANSI_COLOUR_BLUE    = achar(27)//"[34m"
    character(*), parameter, public :: ANSI_COLOUR_CYAN    = achar(27)//"[35m"
    character(*), parameter, public :: ANSI_COLOUR_MAGENTA = achar(27)//"[36m"

    ! Set the colour of a character variable
    public :: hex
    public :: colour
    public :: bold
    public :: apparentLength
    public :: c_f_string
    public :: c_f_buffer
    public :: FTimeToString
    public :: substitute_characters
    public :: to_upper
    public :: to_lower
    public :: upper_case
    public :: lower_case
    public :: hex_string
    public :: is_whitespace
    public :: split

    ! Interface to grisu3
    interface
        ! dtoa_grisu3(double v, char *dst);
        pure function grisu3(v, dst) bind(C, name="dtoa_grisu3")
            use iso_c_binding
            real(c_double), value, intent(in) :: v ! number to write
            type(c_ptr), value, intent(in) :: dst ! Character string; length should be 25
            integer(c_int) :: grisu3 ! Number of characters written
        end function
    end interface
    public :: grisu3

    ! Interface to POSIX's isatty
    interface
        function isatty(id) bind(C, name="isatty")
            use iso_c_binding
            integer(c_int), intent(in), value :: id
            integer(c_int) :: isatty
        end function
    end interface

contains
    subroutine split(string, separator, strings, group)
        character(*), intent(in) :: string
        character(1), intent(in) :: separator
        character(:), dimension(:), allocatable, intent(out) :: strings
        logical, intent(in), optional :: group
!------
        integer :: i, j, n, l
        logical :: group_
!------
        if (present(group)) then
            group_ = group
        else
            group_ = .false.
        end if

        ! Count the number of strings and their maximum length
        n = 0
        l = 0
        j = 0
        do i=1, len(string)
            j = j + 1
            if (string(i:i) == separator) then
                if (group_ .and. i>1) then
                    if (string(i-1:i-1) == separator) then
                        j = 0
                        cycle
                    end if
                end if
                n = n + 1
                l = max(l, j-1)
                j = 0
            end if
        end do
        if (string(len(string):len(string)) /= separator) then
            n = n + 1
            l = max(l, j)
        end if

        ! Allocate the substrings array
        allocate(character(l) :: strings(n))

        ! Copy the substrings
        n = 0
        l = 0
        j = 0
        do i=1, len(string)
            j = j + 1
            if (string(i:i) == separator) then
                if (group_ .and. i>1) then
                    if (string(i-1:i-1) == separator) then
                    j = 0
                        cycle
                    end if
                end if
                n = n + 1
                strings(n) = string(i-j+1:i-1)
                l = max(l, j-1)
                j = 0
            end if
        end do
        if (string(len(string):len(string)) /= separator) strings(n+1) = string(len(string)-j+1:)
    end subroutine
!******************************************************************************!
!> Determine whether a string is white space.
!******************************************************************************!
    function is_whitespace(string)
        character(*), intent(in) :: string
        logical :: is_whitespace
!------
        integer :: i
!------
        is_whitespace = .true.
        do i=1, len(string)
            is_whitespace = is_whitespace .and. ( &
                iachar(string(i:i)) == 32 .or. &
                iachar(string(i:i)) == 9 .or. &
                iachar(string(i:i)) == 10 .or. &
                iachar(string(i:i)) == 13)
        end do
    end function
!******************************************************************************!
!> Get a hexadecimal representation of a character string.
!******************************************************************************!
    function hex_string(string) result(out)
        character(*), intent(in) :: string
        character(:), allocatable :: out
!------
        integer :: i
!------
        if (len(string) <= 0) return
        allocate(character(len(string)*3-1) :: out)
        do i = 1, len(string)
            write(out((i-1)*3+1:(i-1)*3+3), '(Z2.2,1x)') ichar(string(i:i))
        end do
    end function
!******************************************************************************!
!> Convert a character string to upper case.
!>
!> This is the function version for convenience. The more efficient, in-place
!> version is the `[[to_upper(subroutine)]]` subroutine.
!******************************************************************************!
    pure function upper_case(string) result(out)
        character(*), intent(in) :: string
        character(:), allocatable :: out
!------
        out = string
        call to_upper(out)
    end function
!******************************************************************************!
!> Convert a character string to lower case.
!>
!> This is the function version for convenience. The more efficient, in-place
!> version is the `[[to_lower(subroutine)]]` subroutine.
!******************************************************************************!
    pure function lower_case(string) result(out)
        character(*), intent(in) :: string
        character(:), allocatable :: out
!------
        out = string
        call to_lower(out)
    end function
!******************************************************************************!
!> Convert a character string to upper case.
!******************************************************************************!
    pure subroutine to_upper(string)
        character(*), intent(inout) :: string
!------
        integer :: shift, i
!------
        shift = iachar("A") - iachar("a")
        do i=1, len(string)
            if (string(i:i) >= "a" .and. string(i:i) <= "z") string(i:i) = char(iachar(string(i:i)) + shift)
        end do
    end subroutine
!******************************************************************************!
!> Convert a character string to lower case.
!******************************************************************************!
    pure subroutine to_lower(string)
        character(*), intent(inout) :: string
!------
        integer :: shift, i
!------
        shift = iachar("A") - iachar("a")
        do i=1, len(string)
            if (string(i:i) >= "A" .and. string(i:i) <= "Z") string(i:i) = char(iachar(string(i:i)) - shift)
        end do
    end subroutine

    pure function substitute_characters(string, from, to) result(str)
        character(*), intent(in) :: string
        character(1), intent(in), optional :: from
        character(1), intent(in), optional :: to
        character(:), allocatable :: str
!------
        integer :: i
        character(1) :: a, b
!------
        if (present(from)) then
            a = from
        else
            a = " "
        end if
        if (present(to)) then
            b = to
        else
            b = "_"
        end if

        str = string

        do i=1, len(str)
            if (str(i:i) == a) str(i:i) = b
        end do
    end function
!******************************************************************************!
!> Convert a time in seconds to a human-friendly character variable. Format is
!> hh:mm:ss.ss s
!******************************************************************************!
    function FTimeToString(time) result(string)
        real(e), intent(in) :: time
        character(:), allocatable :: string
!------
        real(e) :: time_
        integer :: hourCount, minuteCount, secondCount, milliSecondCount
!------

        time_ = time

        ! The time is larger than one hour
        if (time > 3600) then
            hourCount = floor(time_ / 3600)
            time_ = time_ - (3600*hourCount)
        else
            hourCount = 0
        end if

        ! The time is larger than one minute
        if (time > 60) then
            minuteCount = floor(time_ / 60)
            time_ = time_ - (60*minuteCount)
        else
            minuteCount = 0
        end if

        ! The time is larger than one second
        if (time > 1) then
            secondCount = floor(time_)
            time_ = time_ - secondCount
        else
            secondCount = 0
        end if

        millisecondCount = nint(time_ * 1000)

        if (hourCount > 0) then
            string = hourCount // ":" // minuteCount // ":" // secondCount // "." // millisecondCount // " s"
        else
            string = minuteCount // ":" // secondCount // "." // millisecondCount // " s"
        end if
    end function
!******************************************************************************!
!> Conversion to hexadecimal representation.
!******************************************************************************!
    function hex(input)
        integer, intent(in) :: input
        character(:), allocatable :: hex
!------
        integer :: i, s
        integer, dimension(:), allocatable :: int
!------
        s = storage_size(input) / storage_size(int) ! s is the size in units of the size of int

        allocate(int(s))
        int(:) = transfer(input, int)
        allocate(character(s*8) :: hex)
        hex(:) = ""
        do i=1, s
            write(hex((i-1)*8+1:i*8),'(4z8.8)') int(i)
        end do
    end function
!******************************************************************************!
!> Conversion from C to Fortran strings.
!******************************************************************************!
    subroutine c_f_string(stringptr, fstring)
        type(c_ptr), intent(in) :: stringptr
        character(:), allocatable, intent(out) :: fstring
!------
        character(1,c_char), dimension(:), pointer :: buffer
        integer(c_size_t) :: i
        integer(c_size_t), dimension(1) :: l

        interface
            function strlen(string) bind(C)
                import
                type(c_ptr), value :: string
                integer(c_size_t) :: strlen
            end function
        end interface
!------
        if (.not.c_associated(stringptr)) return
        l(1) = strlen(stringptr)
        call c_f_pointer(stringptr, buffer, l)
        allocate(character(l(1)) :: fstring)
        do i=1, l(1)
            fstring(i:i) = buffer(i)
        end do
    end subroutine

    subroutine c_f_buffer(buffptr, buffer, length)
        type(c_ptr), intent(in) :: buffptr
        character(*), intent(out) :: buffer
        integer(c_size_t), intent(in) :: length
!------
        character(1,c_char), dimension(:), pointer :: ptr
        integer(c_size_t) :: i
!------
        if (.not.c_associated(buffptr)) return

        call c_f_pointer(buffptr, ptr, [length])
        do i=1, length
            buffer(i:i) = ptr(i)
        end do
    end subroutine
!******************************************************************************!
!> Format a string by turning it into bold characters using an ANSI escape
!> sequence.
!******************************************************************************!
    function bold(string)
        character(*), intent(in) :: string
        character(:), allocatable :: bold
!------
        if (isatty(1) == 1) then
            bold = char(27)//"[1m"//string//char(27)//"[0m"
        else
            bold = string
        end if
    end function
!******************************************************************************!
!> Set the colour of a character variable for terminal output
!******************************************************************************!
    pure function colour(string, colourName) result(this)
        character(*), intent(in)  :: string
        character(*), intent(in)  :: colourName
        character(:), allocatable :: this
!------
        select case(colourName)
            case ("red")
                this = ANSI_COLOUR_RED//string//ANSI_COLOUR_DEFAULT
            case ("green")
                this = ANSI_COLOUR_GREEN//string//ANSI_COLOUR_DEFAULT
            case ("yellow")
                this = ANSI_COLOUR_YELLOW//string//ANSI_COLOUR_DEFAULT
            case ("blue")
                this = ANSI_COLOUR_BLUE//string//ANSI_COLOUR_DEFAULT
            case ("cyan")
                this = ANSI_COLOUR_CYAN//string//ANSI_COLOUR_DEFAULT
            case ("magenta")
                this = ANSI_COLOUR_MAGENTA//string//ANSI_COLOUR_DEFAULT
        end select
    end function
!******************************************************************************!
!> Get the length of a string object. The length returned is the number of
!> unicode glyphs, which might be different from the length if the internal
!> representation. It also ignore ANSI escape codes.
!******************************************************************************!
    function apparentLength(this) result(length)
        character(*), intent(in) :: this
        integer :: length !< Length of the string
!------
        integer :: i, code
!------
        length = 0

        i = 1

        ! Process each character in the string
        do while (i <= len(this))

            ! Get the code for the current character
            code = ichar(this(i:i))

            if (code == 27) then
                ! The current character is ESC

                ! Make sure that it is not the last character
                ! Exit if it is (the effective length does not need to be incremented)
                if (i >= len(this)-1) exit

                ! Get the code of the character following ESC
                code = ichar(this(i+1:i+1))

                if (this(i+1:i+1) == "[") then
                    ! This is a multiple-character sequence

                    ! Look for the end of the sequence
                    i = i + 1
                    do while(i <= len(this))
                        i = i + 1
                        code = ichar(this(i:i))

                        if (code >= 64 .and. code <= 126) exit
                    end do
                else
                    ! This is not the begining of an ANSI sequence, or
                    ! an unsupported two-character ANSI sequence
                    i = i + 1
                    cycle
                end if

            else if (code < int(z'80') .or. code > int(z'bf')) then
                ! The current character is not a continuation
                ! The string effective length needs to be incremented
                length = length + 1
            end if

            i = i + 1
        end do
    end function
!******************************************************************************!
!> Concatenate a character variable and a representation of an integer array
!******************************************************************************!
    pure function character_appendIntegerArray(string, x) result(out)
        character(*), intent(in) :: string
        integer, dimension(:), intent(in) :: x
        character(:), allocatable :: out
!------
        out = x
        out = string // out
    end function
!******************************************************************************!
!> Concatenate a representation of an integer array and a character variable
!******************************************************************************!
    pure function integerArray_appendCharacter(x, string) result(out)
        character(*), intent(in) :: string
        integer, dimension(:), intent(in) :: x
        character(:), allocatable :: out
!------
        out = x
        out = out // string
    end function
!******************************************************************************!
!> Concatenate a representation of a single-precision real array and a
!> character variable
!******************************************************************************!
    pure function character_appendSingleArray(string, x) result(out)
        character(*), intent(in) :: string
        real(sp), dimension(:), intent(in) :: x
        character(:), allocatable :: out
!------
        out = x
        out = string // out
    end function
!******************************************************************************!
!> Concatenate a character variable with a representation of a single-precision
!> real array
!******************************************************************************!
    pure function singleArray_appendCharacter(x, string) result(out)
        character(*), intent(in) :: string
        real(sp), dimension(:), intent(in) :: x
        character(:), allocatable :: out
!------
        out = x
        out = out // string
    end function
!******************************************************************************!
!> Concatenate a representation of a double-precision real array and a
!> character variable
!******************************************************************************!
    pure function character_appendDoubleArray(string, x) result(out)
        character(*), intent(in) :: string
        real(dp), dimension(:), intent(in) :: x
        character(:), allocatable :: out
!------
        out = x
        out = string // out
    end function
!******************************************************************************!
!> Concatenate a character variable with a representation of a double-precision
!> real array
!******************************************************************************!
    pure function doubleArray_appendCharacter(x, string) result(out)
        character(*), intent(in) :: string
        real(dp), dimension(:), intent(in) :: x
        character(:), allocatable :: out
!------
        out = x
        out = out // string
    end function
!******************************************************************************!
!> Append an integer value to a character variable
!******************************************************************************!
    pure function character_appendInteger(string, x) result(out)
        character(*), intent(in) :: string
        integer, intent(in) :: x
        character(:), allocatable :: out
!------
        out = x
        out = string // out
    end function
!******************************************************************************!
!> Append an integer value to a character variable
!******************************************************************************!
    pure function integer_appendCharacter(x, string) result(out)
        integer, intent(in) :: x
        character(*), intent(in) :: string
        character(:), allocatable :: out
!------
        out = x
        out = out // string
    end function
!******************************************************************************!
!> Append a single-precision value to a character variable
!******************************************************************************!
    pure function character_appendSingle(string, x) result(out)
        character(*), intent(in) :: string
        real(sp), intent(in) :: x
        character(:), allocatable :: out
!------
        out = x
        out = string // out
    end function
!******************************************************************************!
!> Append a character variable to a single-precision value
!******************************************************************************!
    pure function single_appendCharacter(x, string) result(out)
        real(sp), intent(in) :: x
        character(*), intent(in) :: string
        character(:), allocatable :: out
!------
        out = x
        out = out // string
    end function
!******************************************************************************!
!> Append a double-precision value to a character variable
!******************************************************************************!
    pure function character_appendDouble(string, x) result(out)
        character(*), intent(in) :: string
        real(dp), intent(in) :: x
        character(:), allocatable :: out
        character(:), allocatable :: tmp
!------
        out = x
        allocate(character(1000) :: tmp)
        tmp = string // out
        out = trim(tmp)
    end function
!******************************************************************************!
!> Append a character variable to a doube-precision value
!******************************************************************************!
    pure function double_appendCharacter(x, string) result(out)
        real(dp), intent(in) :: x
        character(*), intent(in) :: string
        character(:), allocatable :: out
!------
        out = x
        out = out // string
    end function
!******************************************************************************!
!> Append a logical value to a character variable
!******************************************************************************!
    pure function character_appendLogical(string, x) result(out)
        character(*), intent(in) :: string
        logical, intent(in) :: x
        character(:), allocatable :: out
!------
        out = x
        out = string // out
    end function
!******************************************************************************!
!> Append a character variable to a logical value
!******************************************************************************!
    pure function logical_appendCharacter(x, string) result(out)
        logical, intent(in) :: x
        character(*), intent(in) :: string
        character(:), allocatable :: out
!------
        out = x
        out = out // string
    end function
!******************************************************************************!
!> Create a character variable representing an integer number
!******************************************************************************!
    pure subroutine characterWithInteger(this, i)
        character(:), allocatable, intent(out) :: this
        integer, intent(in) :: i
!------
        integer :: n, iostat
        character(10) :: format
!------
        if (i == 0) then
            this = "0"
        else if (i > 0) then
            n = floor(log10(1.0_e*i)) + 1

            allocate(character(n) :: this)
            write(format,'(a2,i2,a1)') "(i", n, ")"
            write(this,format,iostat=iostat) i

            if (iostat /= 0) this = "<format error>"

        else if (i < 0) then
            n = floor(log10(-1.0_e*i)) + 2

            allocate(character(n) :: this)
            write(format,'(a2,i2,a1)') "(i", n, ")"
            write(this,format,iostat=iostat) i

            if (iostat /= 0) this = "<format error>"
        end if
    end subroutine
!******************************************************************************!
!> Create a character variable representing an integer array
!******************************************************************************!
    pure subroutine characterWithIntegerArray(this, x)
        character(:), allocatable, intent(out) :: this
        integer, dimension(:), intent(in) :: x
!------
        integer :: n
!------
        if (size(x) > 0) then
            ! Get the number list
            this = x(1)
            do n=2, size(x)
                this = this // ", " // x(n)
            end do

            ! Add the brackets
            if (size(x) <= 4) then
                this = "[" // this // "]"
            else
                this = "(" // this // ")"
            end if
        else
            this = "[]"
        end if
    end subroutine
!******************************************************************************!
!> Create a string object representing a single-precision number
!******************************************************************************!
    pure subroutine characterWithSingle(this, number)
        character(:), allocatable, intent(out) :: this
        real(sp), intent(in) :: number
!------
        include 'characterWithReal.inc'
    end subroutine

    pure subroutine characterWithSingleArray(this, x)
        character(:), allocatable, intent(out) :: this
        real(sp), dimension(:), intent(in) :: x
!------
        integer :: n
!------
        if (size(x) > 0) then
            ! Get the number list
            this = x(1)
            do n=2, size(x)
                this = this // ", " // x(n)
            end do

            ! Add the brackets
            if (size(x) <= 4) then
                this = "[" // this // "]"
            else
                this = "(" // this // ")"
            end if
        else
            this = "[]"
        end if
    end subroutine
!******************************************************************************!
!> Create a string object representing a double-precision real number
!******************************************************************************!
    pure subroutine characterWithDouble(this, number)
        use iso_c_binding

        character(:), allocatable, intent(out) :: this
        real(dp), intent(in) :: number
!------
        character(25,c_char), target :: buf
        integer :: i
!------
        i = grisu3(number, c_loc(buf))
        this = buf(:i)
    end subroutine
!******************************************************************************!
!> Create a string object representing a double-precision real array
!******************************************************************************!
    pure subroutine characterWithDoubleArray(this, x)
        character(:), allocatable, intent(out) :: this
        real(dp), dimension(:), intent(in) :: x
!------
        integer :: n
!------
        if (size(x) > 0) then
            ! Get the number list
            this = x(1)
            do n=2, size(x)
                this = this // ", " // x(n)
            end do

            ! Add the brackets
            if (size(x) <= 4) then
                this = "[" // this // "]"
            else
                this = "(" // this // ")"
            end if
        else
            this = "[]"
        end if
    end subroutine
!******************************************************************************!
!> Create a string object representing a logical variable
!******************************************************************************!
    pure subroutine characterWithLogical(this, value)
        character(:), allocatable, intent(out) :: this
        logical, intent(in) :: value
!------
        if (value) then
            this = YES
        else
            this = NO
        end if
    end subroutine
end module
