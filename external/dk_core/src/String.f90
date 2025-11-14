!******************************************************************************!
!                             mod_string module
!______________________________________________________________________________!
!> `[[FString(type)]]` derived type and TBPs, to handle variable-length
!> character strings.
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
module dk_string
    implicit none (external, type)
    private
    save

!******************************************************************************!
!> Derived type containing a variable-length character scalar and related
!> type-bound procedures.
!******************************************************************************!
    type, public :: FString
        private
        character(:), allocatable :: string_ !< Internal representation
        procedure(FStringComparisonFunction), pointer :: compare => null() !< Function to use when comparing with another string object
    contains
        procedure, public :: string         !< Get the internal representation of the string
        procedure, public :: trim => trim_string !< Get the trimmed string
        procedure, public :: length         !< Get the (logical) length of the string
        procedure, public :: internalLength !< Get the length of the internal representation
        procedure, public :: substring      !< Get a substring

        ! Test the content of the string
        procedure, public :: isBlank        !< Check whether the string is empty
        procedure, public :: startsWith     !< Check whether the string starts with a given substring
        procedure, public :: endsWith       !< Check whether the string ends with a given substring

        ! Overloaded concatenation operator
        procedure, private :: appendToChar  !< Append the string to a character variable
        procedure, private :: appendString  !< Concatenate two strings
        procedure, private, pass(this) :: appendCharTo !< Append a character variable to the string
        generic, public :: operator(//) => &
            appendToChar, appendCharTo, appendString   !< Concatenate the string with
                                                       !< other strings or character
                                                       !< variables

        ! Overloaded == comparison operator
        procedure, private :: equalsChar   !< Check whether the string is equivalent to a character variable
        procedure, private :: equalsString !< Check whether two strings are equivalent
        procedure, private, pass(this) :: charEquals   !< Check whether a character variable is equivalent to the string
        !> Test whether the string is equal to another string or a character variable
        generic, public :: operator(==) => equalsChar, equalsString, charEquals

        ! Overloaded /= comparison operator
        procedure, private :: differsFromChar   !< Check whether the string differs from a character variable
        procedure, private :: differsFromString !< Check whether two strings differ
        procedure, private, pass(this) :: charDiffers !< Check whether a character variable differs from the string
        !> Test whether the string is equal to another string or a character variable
        generic, public :: operator(/=) => differsFromChar, charDiffers, differsFromString

        ! Overloaded > comparison operator
        procedure, private :: greaterThanChar   !< Check whether a string is greater than a character variable
        procedure, private :: greaterThanString !< Check whether the string is greater than another string
        !> Test whether the string is greater than another string or a character variable
        generic, public :: operator(>) => greaterThanChar, greaterThanString

        ! Overloaded >= comparison operator
        procedure, private :: greaterOrEqual_Char   !<
        procedure, private :: greaterOrEqual_String !<
        !>
        generic, public :: operator(>=) => greaterOrEqual_Char, greaterOrEqual_String

        ! Overloaded < comparison operator
        procedure, private :: lessThanChar
        procedure, private :: lessThanString
        !> Test whether the string is smaller than another string or a character variable
        generic, public :: operator(<) => lessThanChar, lessThanString

        ! Overloaded <= comparison operator
        procedure, private :: lessOrEqual_Char   !<
        procedure, private :: lessOrEqual_String !<
        !>
        generic, public :: operator(<=) => lessOrEqual_Char, lessOrEqual_String

        ! Overloaded assignment for scalars
        procedure, private :: string_assignString                  !<
        procedure, private :: string_assignCharacter               !< Conversion from character to string
        procedure, pass(this), private :: string_assignToCharacter !< Conversion from string to character
        !> Assign a character variable or another string to the string
        generic, public :: assignment(=) => string_assignCharacter, string_assignToCharacter, &
            string_assignString
    end type

!******************************************************************************!
!> Constructor interface.
!******************************************************************************!
    interface FString
        module procedure string_initWithChar, string_initWithString
    end interface

!******************************************************************************!
!> Interface for comparison functions. The purpose of this interface is to
!> allow user-defined comparison functions to implement different ways of
!> comparing strings.
!******************************************************************************!
    abstract interface
       pure function FStringComparisonFunction(from, to) result(comp)
            import
            class(FString), intent(in) :: from !< String to compare
            character(*),   intent(in) :: to   !< String to compare to
            integer :: comp !< `0` if both strings are equal, `1` if `from` is
                            !< greater than `to`, and `-1` otherwise
        end function
    end interface
    public :: FStringComparisonFunction

    public :: naturalComparisonFunction
    public :: fortranComparisonFunction

    public :: sortStringArray

    interface reallocate
        module procedure reallocate_string_1d
    end interface
    public :: reallocate
contains
!******************************************************************************!
!> Get a trimmed copy of a string.
!******************************************************************************!
    pure function trim_string(this) result(string)
        class(FString), intent(in) :: this !< Subject string
        type(FString) :: string !< Trimmed copy of the subjects tring
!------
        if (allocated(this%string_)) then
            string = trim(this%string_)
        else
            string = ""
        end if
    end function
!******************************************************************************!
!> Check whether a string begins with a given substring.
!******************************************************************************!
    pure function isBlank(this) result(blank)
        class(FString), intent(in) :: this !< Subject string
        logical :: blank !< `.true.` if the string is empty, `.false.` otherwise
!------
        if (allocated(this%string_)) then
            blank = (this%string_ == "")
        else
            blank = .true.
        end if
    end function
!******************************************************************************!
!> Check whether a string begins with a given substring.
!******************************************************************************!
    pure function startsWith(this, pattern) result(match)
        class(FString), intent(in) :: this    !< Subject string
        character(*),   intent(in) :: pattern !< Pattern to test
        logical :: match !< `.true.` if the string `this` starts with `pattern`,
                         !< `.false.` otherwise
!------
        match = .false.
        if (len(pattern) > len(this%string_)) return

        if (allocated(this%string_)) then
            match = (this%string_(1:len(pattern)) == pattern)
        end if
    end function
!******************************************************************************!
!> Check whether a string begins with a given substring.
!******************************************************************************!
    pure function endsWith(this, pattern) result(match)
        class(FString), intent(in) :: this    !< Subject string
        character(*),   intent(in) :: pattern !< Pattern to test
        logical :: match !< `.true.` if the string `this` ends with `pattern`,
                         !< `.false.` otherwise
!------
        match = .false.
        if (len(pattern) > len(this%string_)) return

        if (allocated(this%string_)) then
            match = (this%string_(len(this%string_)-len(pattern)+1:) == pattern)
        end if
    end function
!******************************************************************************!
!> Initialise a string object with a character variable.
!******************************************************************************!
    pure function string_initWithChar(char, compare) result(this)
        character(*), intent(in) :: char !< Character variable with the string's content
        character(*), intent(in), optional :: compare !< Pointer to the comparison
                                                      !< function to be used with this string
        type(FString) :: this !< Initialised string
!------
        this%string_ = char

        ! Set the comparison function
        if (present(compare)) then
            ! Use a standard comparison function
            select case(compare)
                case ("fortran")
                    this%compare => fortranComparisonFunction
                case("natural")
                    this%compare => naturalComparisonFunction
                case default
                    this%compare => naturalComparisonFunction
            end select
        else
            ! Default comparison function: intrinsic comparisons
            this%compare => naturalComparisonFunction
        end if
    end function
!******************************************************************************!
!> Initialise a string object with another string.
!******************************************************************************!
    pure function string_initWithString(string, compare) result(this)
        type(FString), intent(in) :: string !< String object with the string's content
        character(*), intent(in), optional :: compare !< Pointer to the comparison
                                                      !< function to be used with this string
        type(FString) :: this !< Initialised string
!------
        if (present(compare)) then
            ! Use the provided comparison function if there is one
            this = string_initWithChar(string%string_, compare)
        else
            ! Otherwise, use the same comparison function as the model
            this = string_initWithChar(string%string_)
        end if
    end function
!******************************************************************************!
!> Get a character variable containing the representation of a string object.
!> This is mainly useful in contexts where character variables are accepted but
!> not `[[FString(type)]]` objects.
!******************************************************************************!
    pure function string(this) result(char)
        class(FString), intent(in) :: this
        character(:), allocatable :: char  !< Content of the string
!------

        if (.not.allocated(this%string_)) then
            char = ""
        else
            char = this%string_
        end if
    end function
!******************************************************************************!
!> Get the length of a string object. The length returned is the number of
!> unicode glyphs, which might be different from the length if the internal
!> representation. It also ignore ANSI escape codes.
!******************************************************************************!
    elemental function length(this)
        class(FString), intent(in) :: this
        integer :: length !< Length of the string
!------
        integer :: i, code
!------
        length = 0
        if (allocated(this%string_)) then
            i = 1

            ! Process each character in the string
            do while (i <= len(this%string_))

                ! Get the code for the current character
                code = ichar(this%string_(i:i))

                if (code == 27) then
                    ! The current character is ESC

                    ! Make sure that it is not the last character
                    ! Exit if it is (the effective length does not need to be incremented)
                    if (i >= len(this%string_)-1) exit

                    ! Get the code of the character following ESC
                    code = ichar(this%string_(i+1:i+1))

                    if (this%string_(i+1:i+1) == "[") then
                        ! This is a multiple-character sequence

                        ! Look for the end of the sequence
                        i = i + 1
                        sequenceLoop: do while(i <= len(this%string_))
                            i = i + 1
                            code = ichar(this%string_(i:i))

                            if (code >= 64 .and. code <= 126) exit sequenceLoop
                        end do sequenceLoop
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
        end if
    end function
!******************************************************************************!
!> Get the length of the internal representation of a string object. The length
!> does not take into account unicode glyphs or ANSI codes, and thus can be
!> greater than the number of glyphs (e.g. the length of the string when printed
!> to a terminal).
!******************************************************************************!
    pure function internalLength(this)
        class(FString), intent(in) :: this
        integer :: internalLength !< Length of the internal representation of the string
!------
        internalLength = 0
        if (allocated(this%string_)) internalLength = len(this%string_)
    end function
!******************************************************************************!
!> Comparison function using standard Fortran order.
!******************************************************************************!
    pure function fortranComparisonFunction(from, to) result(comp)
        class(FString), intent(in) :: from
        character(*),   intent(in) :: to
        integer :: comp !< `0` if both strings are equal, `1` if `from` is
                        !< greater than `to`, and `-1` otherwise
!------
        if (from%string_ < to) then
            comp = -1
        else if (from%string_ > to) then
            comp = 1
        else
            comp = 0
        end if
    end function
!******************************************************************************!
!> Internal comparison function.
!******************************************************************************!
    pure recursive function naturalComparison(from, to) result(comp)
        character(*), intent(in) :: from
        character(*),   intent(in) :: to
        integer :: comp !< `0` if both strings are equal, `1` if `from` is
                        !< greater than `to`, and `-1` otherwise
!------
        integer :: start_from, end_from, start_to, end_to, i
        integer :: s, ifrom, ito, n
        logical :: void_from, void_to
        character(1) :: b, c, d
        logical :: b_number, c_number, d_number
!------

!------ Handle simple cases
        ! Determine whether any of the strings is empty
        void_to = (to == "")
        void_from = (from == "")

        ! Both strings are empty
        if (void_from .and. void_to) then
            comp = 0
            return
        end if

        ! The reference string is empty
        if (void_to.and..not.void_from) then
            comp = +1
            return
        end if

        ! The other string is empty
        if (void_from.and..not.void_to) then
            comp = -1
            return
        end if

        ! Both strings are identical
        if (from == to) then
            comp = 0
            return
        end if

        ! One string is a substring of the other
        s = min(len(from), len(to))
        if (from(:s) == to(:s)) then
            ! The greatest string is the longer one
            if (len(from) > len(to)) then
                comp = +1
            else
                comp = -1
            end if
            return
        end if

!------ General case
        ! Look for the first different character
        s = 0
        do i=1, min(len(from), len(to))
            if (from(i:i) == to(i:i)) cycle
            s = i
            exit
        end do

        ! Check whether at least one of the different characters is not a number
        if (s > 1) then
            b = from(s-1:s-1)
        else
            b = " "
        end if
        c = from(s:s)
        d = to(s:s)
        b_number = (iachar(b) >= iachar("0") .and. iachar(b) <= iachar("9"))
        c_number = (iachar(c) >= iachar("0") .and. iachar(c) <= iachar("9"))
        d_number = (iachar(d) >= iachar("0") .and. iachar(d) <= iachar("9"))
        if (.not.c_number .and. .not.d_number) then
            if (c == d) then
                comp =  0
            else if (c > d) then
                comp =  1
            else
                comp = -1
            end if
            return
        else if (c_number .and. .not.d_number .and. .not.b_number) then
            comp = -1
            return
        else if (.not.c_number .and. d_number .and. .not.b_number) then
            comp =  1
            return
        end if

        ! Try to backtrack to get the index of the first digit
        if (b_number) then
            do i=s-1, 1, -1
                if (iachar(from(i:i)) >= iachar("0") .and. iachar(from(i:i)) <= iachar("9")) then
                    start_from = i
                    start_to = i
                else
                    exit
                end if
            end do
        else
            start_from = s
            start_to = s
        end if

        ! Get the last digit in the first string
        end_from = start_from
        do i=start_from, len(from)
            if (iachar(from(i:i)) >= iachar("0") .and. iachar(from(i:i)) <= iachar("9")) then
                end_from = i
            else
                exit
            end if
        end do

        ! Get the last digit in the second string
        end_to = start_to
        do i=start_to, len(to)
            if (iachar(to(i:i)) >= iachar("0") .and. iachar(to(i:i)) <= iachar("9")) then
                end_to = i
            else
                exit
            end if
        end do

        ifrom = 0
        n = -1
        do i=end_from, start_from, -1
            n = n + 1
            ifrom = ifrom + (iachar(from(i:i)) - 48) * 10**n
        end do

        ito = 0
        n = -1
        do i=end_to, start_to, -1
            n = n + 1
            ito = ito + (iachar(to(i:i)) - 48) * 10**n
        end do

        if (ifrom == ito) then
            comp = naturalComparison(from(end_from+1:), to(end_to+1:))
        else if (ifrom > ito) then
            comp =  1
        else if (ifrom < ito) then
            comp = -1
        end if
    end function
!******************************************************************************!
!> Comparison function using "natural" order.
!******************************************************************************!
    pure function naturalComparisonFunction(from, to) result(comp)
        class(FString), intent(in) :: from
        character(*),   intent(in) :: to
        integer :: comp !< `0` if both strings are equal, `1` if `from` is
                        !< greater than `to`, and `-1` otherwise
!------
        comp = 0
        if (allocated(from%string_)) then
            comp = naturalComparison(from%string_, to)
        else
            comp = naturalComparison("", to) !LCOV_EXCL_LINE
        end if
    end function
!******************************************************************************!
!> Assign a string object to a character variable. In contrast with the
!> constructor `[[FString(interface)]]`, the comparison function cannot be
!> provided. The resulting `[[FString(type)]]` object will use the default
!> comparison function.
!******************************************************************************!
    pure subroutine string_assignToCharacter(char, this)
        character(:), allocatable, intent(inout) :: char
        class(FString), intent(in) :: this !< Initialised string
!------
        if (allocated(char)) then
            if (len(char) /= len(this%string_)) deallocate(char)
        end if
        if (allocated(this%string_) .and. .not.allocated(char)) then
            allocate(char, source=this%string_)
        else if (allocated(this%string_)) then
            char(:) = this%string_(:)
        end if
    end subroutine
!******************************************************************************!
!> Assignment to a string object. In contrast with the constructor
!> `[[FString(interface)]]`, the comparison function cannot be provided. The
!> resulting `[[FString(type)]]` object will use the default comparison
!> function.
!******************************************************************************!
    pure subroutine string_assignCharacter(this, char)
        class(FString), intent(inout) :: this
        character(*), intent(in) :: char !< Initialised string
!------
        this%string_ = char
        this%compare => naturalComparisonFunction
    end subroutine
    pure subroutine string_assignString(this, string)
        class(FString), intent(inout) :: this
        type(FString), intent(in) :: string !< Initialised string
!------
        this%string_ = string%string_
        if (associated(string%compare)) &
            this%compare => string%compare
    end subroutine
!******************************************************************************!
!> Concatenate a string object and a character scalar.
!******************************************************************************!
    pure function appendToChar(this, string) result(out)
        class(FString), intent(in) :: this
        character(*), intent(in) :: string
        character(:), allocatable :: out
!------
        out = this%string_ // string
    end function
!******************************************************************************!
!> Concatenate a character scalar and a string object.
!******************************************************************************!
    pure function appendCharTo(string, this) result(out)
        class(FString), intent(in) :: this
        character(*), intent(in) :: string
        character(:), allocatable :: out
!------
        out = string // this%string_
    end function
!******************************************************************************!
!> Concatenate two string objects.
!******************************************************************************!
    pure function appendString(this, string) result(out)
        class(FString), intent(in) :: this
        class(FString), intent(in) :: string
        character(:), allocatable :: out
!------
        out = this%string_ // string%string_
    end function
!******************************************************************************!
!> Compare a string object and a character variable.
!******************************************************************************!
    elemental function equalsChar(this, char) result(comp)
        class(FString), intent(in) :: this
        character(*), intent(in) :: char
        logical :: comp
!------
        integer :: i
!------
        ! Handle variables with a 0 length
        if (.not.allocated(this%string_)) then
            comp = (char == "")
        else
            ! Call the comparison function
            if (associated(this%compare)) then
                i = this%compare(char)
            else
                i = fortranComparisonFunction(this, char)
            end if
            if (i == 0) then
                comp = .true.
            else
                comp = .false.
            end if
        end if
    end function
!******************************************************************************!
!> Compare a string object and a character variable.
!******************************************************************************!
    elemental function charEquals(char, this) result(comp)
        character(*),   intent(in) :: char
        class(FString), intent(in) :: this
        logical :: comp
!------
        comp = equalsChar(this, char)
    end function
!******************************************************************************!
!> Compare two string objects.
!******************************************************************************!
    elemental function equalsString(this, string) result(comp)
        class(FString), intent(in) :: this
        class(FString), intent(in) :: string
        logical :: comp
!------
        if (allocated(string%string_)) then
            comp = equalsChar(this, string%string_)
        else
            comp = equalsChar(this, "")
        end if
    end function
!******************************************************************************!
!> Compare a string object and a character variable.
!******************************************************************************!
    elemental function differsFromChar(this, char) result(comp)
        class(FString), intent(in) :: this
        character(*), intent(in) :: char
        logical :: comp
!------
        comp = .not.equalsChar(this, char)
    end function
!******************************************************************************!
!> Compare a string object and a character variable.
!******************************************************************************!
    elemental function charDiffers(char, this) result(comp)
        character(*), intent(in) :: char
        class(FString), intent(in) :: this
        logical :: comp
!------
        comp = .not.equalsChar(this, char)
    end function
!******************************************************************************!
!> Compare two string objects.
!******************************************************************************!
    elemental function differsFromString(this, string) result(comp)
        class(FString), intent(in) :: this
        class(FString), intent(in) :: string
        logical :: comp
!------
        if (allocated(string%string_)) then
            comp = .not.equalsChar(this, string%string_)
        else
            comp = .not.equalsChar(this, "")
        end if
    end function
!******************************************************************************!
!> Compare a string object and a character variable.
!******************************************************************************!
    elemental function greaterThanChar(this, char) result(comp)
        class(FString), intent(in) :: this
        character(*), intent(in) :: char
        logical :: comp
!------
        integer :: i
!------

        ! Handle variables with a 0 length
        if (.not.allocated(this%string_)) then
            comp = .false.
        else
            ! Call the comparison function
            if (associated(this%compare)) then
                i = this%compare(char)
            else
                i = fortranComparisonFunction(this, char)
            end if
            if (i > 0) then
                comp = .true.
            else
                comp = .false.
            end if
        end if
    end function
!******************************************************************************!
!> Compare two string objects.
!******************************************************************************!
    elemental function greaterThanString(this, string) result(comp)
        class(FString), intent(in) :: this
        class(FString), intent(in) :: string
        logical :: comp
!------
        if (allocated(string%string_)) then
            comp = greaterThanChar(this, string%string_)
        else
            comp = greaterThanChar(this, "")
        end if
    end function
!******************************************************************************!
!> Compare a string object and a character variable.
!******************************************************************************!
    elemental function greaterOrEqual_Char(this, char) result(comp)
        class(FString), intent(in) :: this
        character(*), intent(in) :: char
        logical :: comp
!------
        integer :: i
!------

        ! Handle variables with a 0 length
        if (.not.allocated(this%string_)) then
            comp = .false.
        else
            ! Call the comparison function
            if (associated(this%compare)) then
                i = this%compare(char)
            else
                i = fortranComparisonFunction(this, char)
            end if
            if (i >= 0) then
                comp = .true.
            else
                comp = .false.
            end if
        end if
    end function
!******************************************************************************!
!> Compare two string objects.
!******************************************************************************!
    elemental function greaterOrEqual_String(this, string) result(comp)
        class(FString), intent(in) :: this
        class(FString), intent(in) :: string
        logical :: comp
!------
        if (allocated(string%string_)) then
            comp = greaterOrEqual_Char(this, string%string_)
        else
            comp = greaterOrEqual_Char(this, "")
        end if
    end function
!******************************************************************************!
!> Compare a string object and a character variable.
!******************************************************************************!
    elemental function lessThanChar(this, char) result(comp)
        class(FString), intent(in) :: this
        character(*), intent(in) :: char
        logical :: comp
!------
        integer :: i
!------

        ! Handle variables with a 0 length
        if (.not.allocated(this%string_) .and. char == "") then
            comp = .false.
        else if (.not.allocated(this%string_) .and. char /= "") then
            comp = .true.
        else
            ! Call the comparison function
            if (associated(this%compare)) then
                i = this%compare(char)
            else
                i = fortranComparisonFunction(this, char)
            end if
            if (i < 0) then
                comp = .true.
            else
                comp = .false.
            end if
        end if
    end function
!******************************************************************************!
!> Compare two string objects.
!******************************************************************************!
    elemental function lessThanString(this, string) result(comp)
        class(FString), intent(in) :: this
        class(FString), intent(in) :: string
        logical :: comp
!------
        if (allocated(string%string_)) then
            comp = lessThanChar(this, string%string_)
        else
            comp = lessThanChar(this, "")
        end if
    end function
!******************************************************************************!
!> Compare a string object and a character variable.
!******************************************************************************!
    elemental function lessOrEqual_Char(this, char) result(comp)
        class(FString), intent(in) :: this
        character(*), intent(in) :: char
        logical :: comp
!------
        integer :: i
!------

        ! Handle variables with a 0 length
        if (.not.allocated(this%string_) .and. char == "") then
            comp = .false.
        else if (.not.allocated(this%string_) .and. char /= "") then
            comp = .true.
        else
            ! Call the comparison function
            if (associated(this%compare)) then
                i = this%compare(char)
            else
                i = fortranComparisonFunction(this, char)
            end if
            if (i <= 0) then
                comp = .true.
            else
                comp = .false.
            end if
        end if
    end function
!******************************************************************************!
!> Compare two string objects.
!******************************************************************************!
    elemental function lessOrEqual_String(this, string) result(comp)
        class(FString), intent(in) :: this
        class(FString), intent(in) :: string
        logical :: comp
!------
        if (allocated(string%string_)) then
            comp = lessOrEqual_Char(this, string%string_)
        else
            comp = lessOrEqual_Char(this, "")
        end if
    end function
!******************************************************************************!
!> Get a substring defined by a beginning and an end indices. The substring
!> bounds corresponds to unicode glyphs, ANSI code are not supported yet.
!> If `from` is not present, the substring will start from the beginning of the
!> string. If `to` is not present, the substring will end at the end of the
!> string.
!******************************************************************************!
    pure function substring(this, from, to)
        type(FString) :: substring  !< Extracted substring; `""` if `from` is greater than `to`.
        class(FString), intent(in) :: this
        integer, intent(in), optional :: from !< Index of the first character of the substring
        integer, intent(in), optional :: to   !< Index of the last character of the substring
!------
        integer :: i, j, c, start, end, n, code
!------
        ! Get the lower bound of the substring
        if (present(from)) then
            i = from
        else
            i = 1
        end if

        ! Get the upper bound of the substring
        if (present(to)) then
            j = to
        else
            j = this%length()
        end if

        ! Special case: empty string
        if (j<i) then
            substring = ""
            return
        end if

        ! Calculate the indices corresponding to the bounds
        start = 0
        end = 0
        n = 0
        do c=1, len(this%string())
            ! Check whether the current character is a continuation
            code = ichar(this%string_(c:c))

            if (code < int(z'80') .or. code > int(z'bf')) then
                ! If not, increment the string length
                n = n + 1
                if (n == i) then
                    start = c
                end if
                if (n-1 == j) then
                    end = c-1
                end if
            end if
        end do
        if (j==this%length()) end = len(this%string_)

        ! Copy the substring
        substring = this%string_(start:end)
    end function
!******************************************************************************!
!> Sort a string array using a quicksort algorithm.
!******************************************************************************!
    pure recursive subroutine sortStringArray(A)
        type(FString), dimension(:), intent(inout) :: A
!------
        integer :: iq
!------
        if (size(A) > 1) then
           call Partition_string(A, iq)
           if (iq >= 50) then
               call sortStringArray(A(:iq-1))
           else
               call insertionSort_string(A(:iq-1))
           end if
           if ((size(A)-iq) >= 50) then
               call sortStringArray(A(iq:))
           else
               call insertionSort_string(A(iq:))
           end if
        endif
    end subroutine

    pure subroutine Partition_string(A, marker)
        type(FString), dimension(:), intent(inout) :: A
        integer, intent(out) :: marker
!------
        integer :: i, j
        type(FString) :: temp, x, small, large
!------
      ! Find the median between A(1), A(size(A)) and A(size(A)/2)
      j = size(A)
      i = size(A)/2
      if (A(1) < A(j)) then
          small = A(1)
          large = A(j)
      else
          small = A(j)
          large = A(1)
      end if
      if (A(i) > large) then
          x = large
      else if (A(i) < small) then
          x = small
      else
          x = A(i)
      end if

      i= 0
      j= size(A) + 1

      do
         j = j-1
         do
            if (A(j) <= x) exit
            j = j-1
         end do
         i = i+1
         do
            if (A(i) >= x) exit
            i = i+1
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = A(i)
            A(i) = A(j)
            A(j) = temp
         elseif (i == j) then
            marker = i+1
            return
         else
            marker = i
            return
         endif
      end do
    end subroutine

    pure subroutine insertionSort_string(A)
      type(FString), dimension(:), intent(inout) :: A
!------
        type(FString) :: temp
        integer :: i, j
!------
        do i = 2, size(A)
            j = i - 1
            temp = A(i)
            do
                if (j == 0) exit
                if (a(j) <= temp) exit
                A(j+1) = A(j)
                j = j - 1
            end do
            a(j+1) = temp
        end do
    end subroutine

    subroutine reallocate_string_1d(array, newsize, minsize)
        type(FString), dimension(:), allocatable, intent(inout) :: array
        integer, intent(in) :: newsize
        integer, intent(in), optional :: minsize
!------
        type(FString), dimension(:), allocatable :: t
!------
        if (size(array) > newsize) return
        if (present(minsize)) then
            if (size(array) > minsize) return
        end if
        allocate(t(newsize))
        if (allocated(array)) then
            t(:size(array)) = array
            deallocate(array)
        end if
        call move_alloc(from=t, to=array)
    end subroutine
end module
