!******************************************************************************!
!                          dk_token module
!------------------------------------------------------------------------------!
!> `[[Token(type)]]` derived type and type-bound procedures. This module
!> provide facilities to help writing parsers.
!
!  Part of the dk_io library version 0.1.
!  Written by Paul Fossati, <paul.fossati@cea.fr>
!  Copyright (c) 2025 CEA
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
module dk_token
    use iso_c_binding, only: C_NULL_CHAR

    use dk_parameters, only: e

    use dk_exception, only: FexceptionDescription, FError
    use ext_character, only: colour

    implicit none(external, type)
    private

    integer, parameter, public :: TOKEN_UNDEF = 0
    integer, parameter, public :: TOKEN_DATA_BLOCK  = 1
    integer, parameter, public :: TOKEN_DATA_NAME = 2
    integer, parameter, public :: TOKEN_DATA_VALUE = 3
    integer, parameter, public :: TOKEN_KEYWORD = 4
    integer, parameter, public :: TOKEN_CIF_VERSION = 5
    integer, parameter, public :: TOKEN_WORD = 6
    integer, parameter, public :: TOKEN_END = 7
    integer, parameter, public :: TOKEN_NEW_LINE = 8

!******************************************************************************!
!> Derived type containing a lexical token and some of its context to help
!> parsers processing it.
!******************************************************************************!
    type, public :: Token
        character(:), allocatable :: value
        character(:), allocatable :: line
        type(Token), allocatable :: next
        integer :: line_number = -1
        integer :: start = -1
        integer :: type = TOKEN_END
    contains
        procedure, public :: get_error_context
        procedure, public :: read_real
        procedure, public :: read_integer
    end type

    public :: TokenError
    public :: get_token
    public :: tokenise
contains

    subroutine tokenise(string, separator, tokens, group)
        character(*), intent(in) :: string
        character(1), intent(in) :: separator
        type(Token), intent(out), target :: tokens
        logical, intent(in), optional :: group
!------
        integer :: i, j, n, l
        logical :: group_
        type(Token), pointer :: current_token, previous_token
!------
        if (present(group)) then
            group_ = group
        else
            group_ = .false.
        end if

        ! Copy the substrings
        n = 0
        l = 0
        j = 0
        current_token => tokens
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

                current_token%value = string(i-j+1:i-1)
                current_token%line = string
                current_token%start = i-j+1

                allocate(current_token%next)
                previous_token => current_token
                current_token => current_token%next

                l = max(l, j-1)
                j = 0
            end if
        end do

        if (string(len(string):len(string)) /= separator) then
            current_token%value = string(len(string)-j+1:)
            current_token%line = string
            current_token%start = len(string)-j+1
        else
            if (associated(previous_token)) deallocate(previous_token%next)
        end if
    end subroutine

    subroutine get_token(string, n, t, offset)
        character(*), intent(in) :: string
        integer, intent(in) :: n
        type(Token), intent(out) :: t
        integer, intent(in), optional :: offset
!------
        integer :: i, start, end, c
!------
        if (present(offset)) then
            start = offset-1
        else
            start = 0
        end if

        end = 0
        c = 0
        i = start
        do while(i < len(string))
            i = i + 1
            if (string(i:i) == " ") cycle

            c = c + 1
            if (c == n) then
                start = i
                end = i
            end if

            do while(i < len(string))
                i = i + 1
                if (string(i:i) /= " " .and. string(i:i) /= C_NULL_CHAR) then
                    if (c==n) then
                        end = i
                    end if
                    cycle
                end if
                exit
            end do

            if (c == n) exit
        end do

        ! Set the token object
        if (start > 0) then
            t%value = string(start:end)
            t%line = string
            t%start = start
            t%type = TOKEN_DATA_VALUE
        end if
    end subroutine
!******************************************************************************!
!> Get an error message to explain why a token is responsible for a parse
!> error.
!******************************************************************************!
    function TokenError(message, context, word) result(error)
        character(*), intent(in) :: message
        character(*), intent(in) :: context
        type(Token), intent(in) :: word
        type(FExceptionDescription) :: error
!------
        character(:), allocatable :: error_context
!------
        if (allocated(word%value)) then
            call word%get_error_context(error_context)
            error = FError(message, &
                context//new_line("")//new_line("")//error_context)
        else
            error = FError(message, context)
        end if
    end function
!******************************************************************************!
!> Try to read a real number from a token.
!******************************************************************************!
    subroutine read_real(this, r, iostat)
        class(Token), intent(in) :: this
        real(e), intent(out) :: r
        integer, intent(out) :: iostat
!------
        if (allocated(this%value)) then
            read(this%value,*, iostat=iostat) r
        else
            r = 0
            iostat = -1
        end if
    end subroutine
!******************************************************************************!
!> Try to read an integer number from a token.
!******************************************************************************!
    subroutine read_integer(this, r, iostat)
        class(Token), intent(in) :: this
        integer, intent(out) :: r
        integer, intent(out) :: iostat
!------
        if (allocated(this%value)) then
            read(this%value,*, iostat=iostat) r
        else
            r = 0
            iostat = -1
        end if
    end subroutine
!******************************************************************************!
!> Format a token's context for user-facing errors.
!******************************************************************************!
    subroutine get_error_context(error_token, context)
        class(Token), intent(in) :: error_token
        character(:), allocatable, intent(out) :: context
!------
        character(1000) :: line_number_buffer, empty, caret
        integer :: line_number_len, line_len, total_len, i
!------
        ! Get the size of the terminal

        ! If the output is not a terminal, then use a size of 100 characters
        total_len = 80

        ! Get the number of positions for the line number
        line_number_buffer = " "
        empty = " "
        write(line_number_buffer,'(i0)') error_token%line_number
        line_number_len = len_trim(line_number_buffer)

        ! Get the length of the line to print
        line_len = total_len - 4 - line_number_len

        ! Get the caret to indicate the error location
        caret = " "
        if (len_trim(error_token%value)==1) then
            caret = colour("^", "red")
        else if (len_trim(error_token%value) > 1) then
            do i=1, len_trim(error_token%value)
                caret(i:i) = "~"
            end do
            caret = colour(trim(caret), "red")
        end if

        ! Get the context
        context = " "//trim(line_number_buffer)//" | "//trim(error_token%line)//new_line("")// &
            empty(:1+line_number_len)//" | "//empty(:error_token%start-1)//trim(caret)
    end subroutine
end module
