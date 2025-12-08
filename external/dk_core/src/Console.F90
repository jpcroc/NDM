!******************************************************************************!
!                            mod_console module
!______________________________________________________________________________!
!> `[[FConsole(type)]]` derived type and TBPs, to handle output to a terminal.
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
module dk_console
    use iso_fortran_env, only: OUTPUT_UNIT, ERROR_UNIT

    use core
    use ext_character
    use event_levels

    use dk_event,  only: FEvent
    use dk_string, only: FString

    implicit none (type)
    private
    save

    type, public :: FConsole
        private
        integer, public :: print_level = EVENT_LEVEL_VERBOSE !< Level above which the messages are printed
        integer, public :: abort_level = EVENT_LEVEL_ERROR   !< Level above which the messages are printed
        integer :: errorUnit = ERROR_UNIT            !< Unit number of the error stream
        integer :: outputUnit = OUTPUT_UNIT          !< Unit number of the output stream
        integer :: messageWidth = 30                 !< Width of the message column of the console output
        logical :: colour_mode = .true.              !< Use color when printing messages
        integer :: lineLength = 100
    contains

        procedure, public :: level

        procedure :: printOption                  !<
        procedure :: traceback                    !<
        procedure :: printError                   !< Print an error to the error stream
        procedure :: printWarning                 !< Print a warning to the error stream
        procedure :: debug                        !< Print a debug message to the error stream
        procedure :: log                          !< Print an logging message to the output stream

        procedure, private :: console_init
        generic, public :: init => console_init

        procedure, private :: print_char
        procedure, private :: print_string
        procedure, private :: print_real
        procedure, private :: print_real_vector
        procedure, private :: print_real_matrix
        procedure, private :: print_integer
        procedure, private :: print_integer_vector
        procedure, private :: print_integer_matrix
        procedure, private :: print_logical
        generic, public :: print => print_char, print_string, &
            print_real, print_real_vector, print_real_matrix, &
            print_integer, print_integer_vector, print_integer_matrix, &
            print_logical

        procedure, public :: printheader
        procedure, public :: setVerbosity
    end type

    character(*), dimension(EVENT_LEVEL_NONE:*), parameter :: head_mono = [ &
        "                ", &
        "Debug:          ", &
        "Log:            ", &
        "                ", &
        "                ", &
        "WARNING:        ", &
        "ERROR:          ", &
        "ERROR:          ", &
        "INTERNAL ERROR: "  &
    ]

    character(*), dimension(EVENT_LEVEL_NONE:*), parameter :: head_colour = [ &
        "                         ", &
        ANSI_COLOUR_YELLOW // "Debug: " // ANSI_COLOUR_DEFAULT // "         ", &
        ANSI_COLOUR_GREEN  // "Log: " // ANSI_COLOUR_DEFAULT // "           ", &
        "                         ", &
        "                         ", &
        ANSI_COLOUR_RED    // "Warning: " // ANSI_COLOUR_DEFAULT // "       ", &
        ANSI_COLOUR_RED    // "Error: " // ANSI_COLOUR_DEFAULT // "         ", &
        ANSI_COLOUR_RED    // "Error: " // ANSI_COLOUR_DEFAULT // "         ", &
        ANSI_COLOUR_RED    // "Internal error: " // ANSI_COLOUR_DEFAULT  &
    ]

    interface
        function getTerminalWidth() bind(C, name="getTerminalWidth")
            use iso_c_binding
            integer(c_int) :: getTerminalWidth
        end function
    end interface

    interface
        function isatty(id) bind(C, name="isatty")
            use iso_c_binding
            integer(c_int), intent(in), value :: id
            integer(c_int) :: isatty
        end function
    end interface

    !> Default console object
    type(FConsole), public :: console
contains
!******************************************************************************!
!> Get the level above which messages will be printed bt a console object.
!******************************************************************************!
    pure function level(this)
        class(FConsole), intent(in) :: this
        integer :: level
!------
        level = this%print_level
    end function
!******************************************************************************!
!> Print an traceback message to the error stream of a console object
!******************************************************************************!
    subroutine traceback(this, string)
        class(FConsole), intent(in) :: this
        character(*), intent(in) :: string
!------
        call print_char(this, &
            message=string, &
            level=EVENT_LEVEL_TRACEBACK &
        )
        stop 1
    end subroutine
!******************************************************************************!
!> Print the description of a command-line option.
!******************************************************************************!
    subroutine printOption(this, option, help)
        class(FConsole), intent(in) :: this
        character(*),    intent(in) :: option
        character(*),    intent(in) :: help
!------
        call this%print_char("  "//bold(option), help)
    end subroutine
!******************************************************************************!
!> Print a warning message.
!******************************************************************************!
    subroutine printWarning(this, message, showPrefix)
        class(FConsole), intent(in) :: this
        character(*), intent(in) :: message
        logical, intent(in), optional :: showPrefix

        if (present(showPrefix)) then
            call print_char(this, level=EVENT_LEVEL_WARNING, message=message, showPrefix=showPrefix)
        else
            call print_char(this, level=EVENT_LEVEL_WARNING, message=message)
        end if
    end subroutine
!******************************************************************************!
!> Print an error message.
!******************************************************************************!
    subroutine printError(this, message, showPrefix)
        class(FConsole), intent(in) :: this
        character(*), intent(in) :: message
        logical, intent(in), optional :: showPrefix
!------
        if (present(showPrefix)) then
            call print_char(this, level=EVENT_LEVEL_ERROR, message=message, showPrefix=showPrefix)
        else
            call print_char(this, level=EVENT_LEVEL_ERROR, message=message)
        end if
        stop 42
    end subroutine
!******************************************************************************!
!> Print a log message.
!******************************************************************************!
    subroutine log(this, message, showPrefix)
        class(FConsole), intent(in) :: this
        character(*), intent(in) :: message
        logical, intent(in), optional :: showPrefix
!------
        if (present(showPrefix)) then
            call print_char(this, level=EVENT_LEVEL_LOG, message=message, showPrefix=showPrefix)
        else
            call print_char(this, level=EVENT_LEVEL_LOG, message=message)
        end if
    end subroutine
!******************************************************************************!
!> Print a debug message.
!******************************************************************************!
    subroutine debug(this, message, showPrefix)
        class(FConsole), intent(in) :: this
        character(*), intent(in) :: message
        logical, intent(in), optional :: showPrefix
!------
        if (present(showPrefix)) then
            call print_char(this, level=EVENT_LEVEL_DEBUG, message=message, showPrefix=showPrefix)
        else
            call print_char(this, level=EVENT_LEVEL_DEBUG, message=message)
        end if
    end subroutine
!******************************************************************************!
!> Print a message accompanied by an standard precision real number to
!>  the output stream of a console object.
!******************************************************************************!
    subroutine print_real(this, message, x)
        class(FConsole), intent(in) :: this
        character(*),    intent(in) :: message
        real(wp),         intent(in) :: x
!------
        call print_char(this, message, "" // x)
    end subroutine
!******************************************************************************!
!> Print a message accompanied by an standard precision real vector to the
!>  output stream of a console object.
!******************************************************************************!
    subroutine print_real_vector(this, message, x)
        class(FConsole), intent(in) :: this
        character(*), intent(in) :: message
        real(wp), dimension(:), intent(in) :: x
!------
        character(200) :: string
        character(27) :: format
        integer :: stat
!------
        format = ""
        write(format,'(a,i1,a)') "('[',f0.6,", size(x)-1,"(', ',f0.6),']')"
        write(string,format,iostat=stat) x
        if (stat == 0) then
            call print_char(this, message, trim(string))
        else
            call print_char(this, message, "<format error>")
        end if
    end subroutine
!******************************************************************************!
!> Print a message accompanied by an standard precision real matrix to the
!>  output stream of a console object.
!******************************************************************************!
    subroutine print_real_matrix(this, message, x)
        class(FConsole), intent(in) :: this
        character(*), intent(in) :: message
        real(wp), dimension(:,:), intent(in) :: x
!------
        character(200) :: string
        character(27) :: format
        integer :: linePosition, i, pos, iostat
        character(:), allocatable :: message_
!------
        format = ""
        pos = max(ceiling(log10(maxval(abs(x)))), 1) + 6
        if (minval(x) < 0) pos = pos + 1 ! Add a position if there are negative numbers
        if (pos < 0) pos = 6
        linePosition = (size(x,1)+1) / 2
        do i=1, size(x,1)
            if (i==linePosition) then
                message_ = message
            else
                message_ = ""
            end if
            if (i==1) then
                write(format,'(a,i1,a,i0,a)') "('⎛',", size(x,2),"(f",pos,".4,1x),'⎞')"
            else if (i==size(x,1)) then
                write(format,'(a,i1,a,i0,a)') "('⎝',", size(x,2),"(f",pos,".4,1x),'⎠')"
            else
                write(format,'(a,i1,a,i0,a)') "('⎜',", size(x,2),"(f",pos,".4,1x),'⎟')"
            end if
            write(string,format,iostat=iostat) x(i,:)
            if (iostat==0) then
                call print_char(this, message_, trim(string))
            else
                call print_char(this, message_, "****")
            end if
        end do
    end subroutine
!******************************************************************************!
!> Print a message accompanied by a logical value to the output stream of
!>  a console object.
!******************************************************************************!
    subroutine print_logical(this, message, x)
        class(FConsole), intent(in) :: this
        character(*), intent(in) :: message
        logical, intent(in) :: x
!------
        call print_char(this, message, "" // x)
    end subroutine
!******************************************************************************!
!> Print a message accompanied by an integer number to the output stream
!>  of a console object.
!******************************************************************************!
    subroutine print_integer(this, message, x)
        class(FConsole), intent(in) :: this
        character(*), intent(in) :: message
        integer, intent(in) :: x
!------
        call print_char(this, message, "" // x)
    end subroutine
!******************************************************************************!
!> Print a message accompanied by an integer vector to the output stream
!>  of a console object.
!******************************************************************************!
    subroutine print_integer_vector(this, message, x)
        class(FConsole), intent(in) :: this
        character(*), intent(in) :: message
        integer, dimension(:), intent(in) :: x
!------
        character(200) :: string
        character(23) :: format
!------
        format = ""
        write(format,'(a,i1,a)') "('[',i0,", size(x)-1,"(', ',i0),']')"
        write(string,format) x
        call print_char(this, message, trim(string))
    end subroutine
!******************************************************************************!
!> Print a message accompanied by an integer matrix to the output stream of a
!> console object.
!******************************************************************************!
    subroutine print_integer_matrix(this, message, x)
        class(FConsole), intent(in) :: this
        character(*), intent(in) :: message
        integer, dimension(:,:), intent(in) :: x
!------
        character(200) :: string
        character(27) :: format
        integer :: linePosition, i, pos, iostat
        character(:), allocatable :: message_
!------
        format = ""
        pos = max(ceiling(log10(maxval(abs(x*1.0_e)))), 1) + 6
        if (minval(x) < 0) pos = pos + 1 ! Add a position if there are negative numbers
        if (pos < 0) pos = 6
        linePosition = (size(x,1)+1) / 2
        do i=1, size(x,1)
            if (i==linePosition) then
                message_ = message
            else
                message_ = ""
            end if
            if (i==1) then
                write(format,'(a,i1,a,i0,a)') "('⎛',", size(x,2),"(i",pos,",1x),'⎞')"
            else if (i==size(x,1)) then
                write(format,'(a,i1,a,i0,a)') "('⎝',", size(x,2),"(i",pos,",1x),'⎠')"
            else
                write(format,'(a,i1,a,i0,a)') "('⎜',", size(x,2),"(i",pos,",1x),'⎟')"
            end if
            write(string,format,iostat=iostat) x(i,:)
            if (iostat==0) then
                call print_char(this, message_, trim(string))
            else
                call print_char(this, message_, "****")
            end if
        end do
    end subroutine
!******************************************************************************!
!> Print a message accompanied by an character string to the output
!>  stream of a console object.
!******************************************************************************!
    subroutine print_string(this, message, string)
        class(FConsole), intent(in) :: this
        class(FString), intent(in) :: message
        class(FString), intent(in), optional :: string
!------
        call print_char(this, message%string())
    end subroutine
!******************************************************************************!
!> Print an section header to the output stream of a console object.
!******************************************************************************!
    subroutine printHeader(this, name)
        class(FConsole), intent(in) :: this
        character(*), intent(in) :: name
!------
        integer :: padlen
        character(:), allocatable :: padding
!------
        call print_char(this)
        call print_char(this, " ┌──────────────────────────────" // &
            "───────────────────────────" // &
            "────────────────┐")
        padlen = 80 - 2 - 6 - 5 - apparentLength(adjustl(name))
        allocate(character(padlen) :: padding)
        padding(:) = " "
        call print_char(this, " │      " // adjustl(name) // padding // "│")
        call print_char(this, " └──────────────────────────────" // &
            "───────────────────────────" // &
            "────────────────┘")
        call this%print
    end subroutine
!******************************************************************************!
!> Print a message (for internal use only).
!******************************************************************************!
   subroutine print_char(this, message, column2, context, level, showPrefix)
        class(FConsole), intent(in)           :: this
        character(*),    intent(in), optional :: message
        character(*),    intent(in), optional :: column2
        character(*),    intent(in), optional :: context
        integer,         intent(in), optional :: level
        logical,         intent(in), optional :: showPrefix
!------
        integer :: unit, level_, i, l, columnWidth
        character(:), allocatable :: prefix_, column1_, column2_, space, newColumn
!------

        ! Set the logging level
        if (present(level)) then
            level_ = level
        else
            level_ = EVENT_LEVEL_VERBOSE
        end if

        ! Return if the level of the message is lower than the minimum level
        if (level_ < this%print_level) return

        ! Set the unit to use to print the message
        if (level_ >= EVENT_LEVEL_WARNING) then
            unit = this%errorUnit
        else
            unit = this%outputUnit
        end if

        ! Print an empty line if none of the optional parameters is set
        if (.not.present(message) .and. .not.present(column2) .and. .not.present(context)) then
            write(unit,*)
            return
        end if

        ! Print an empty line if none of the optional parameters is set
        if (.not.present(message) .and. .not.present(column2) .and. .not.present(context)) then
            write(unit,*)
            return
        end if

        ! Prevent simultaneous output to the console from concurrent threads
        !$omp critical

        ! Set the prefix
        if (this%colour_mode) then
            if (level_ <= EVENT_MAX_LEVEL) then
                prefix_ = trim(head_colour(level_))
            else
                prefix_ = ""
            end if
        else
            if (level_ <= EVENT_MAX_LEVEL) then
                prefix_ = trim(head_mono(level_))
            else
                prefix_ = ""
            end if
        end if

        if (present(showPrefix)) then
            if (.not.showPrefix) then
                prefix_ = ""
            end if
        end if

        ! Add the first column to the prefix
        if (present(message)) then
            column1_ = prefix_ // message
        else
            column1_ = prefix_
        end if

        ! Print the context if any
        if (present(context)) write(unit,'(a)') context

        ! Print the first column
        write(unit,'(a)',advance="no") column1_

        ! Print the second column if needed
        if (present(column2)) then
            column2_ = column2

            ! Complete the first column with padding
            if (apparentLength(column1_) < this%messageWidth) then
                i = this%messageWidth - apparentLength(column1_)
                allocate(character(i) :: space)
                write(space,'(a)')
                write(unit,'(a)',advance="no") space
                deallocate(space)
                columnWidth = this%lineLength - this%messageWidth
            else
                write(unit,'(a)',advance="no") " "
                columnWidth = this%lineLength - mod(apparentLength(column1_), this%lineLength) - 1
            end if

            ! Complete the first line
            l = min(columnWidth, len(column2_))

            ! Get the position of the last space in the string
            if (l < len(column2_)) then
                ! The limiting length is that of the column

                i = index(column2_(:l), " ", back=.true.)
                if (i /= 0) l = i
            end if

            write(unit,'(a)') column2_(:l)
            newColumn = column2_(l+1:)
            column2_ = newColumn

            ! Get the second column size for the other lines
            columnWidth = this%lineLength - this%messageWidth

            ! Get the spaces for the first column
            allocate(character(this%messageWidth) :: space)
            write(space,'(a)')

            l = min(columnWidth, len(column2_))
            do while(l > 0)
                ! Get the position of the last space in the string
                if (l < len(column2_)) then
                    i = index(column2_(:l), " ", back=.true.)
                    if (i /= 0) l = i
                end if

                write(unit,'(a)') space // column2_(:l)
                column2_ = column2_(l+1:)
                l = min(columnWidth, len(column2_))
            end do
        else
            write(unit,'(a)')
        end if

        ! Print a stack trace if needed
        if (level_ >= EVENT_LEVEL_IPE) then
#ifdef HAS_BACKTRACE
            write(unit,'(a)') "Backtrace:"
            call backtrace
#endif
        end if

        ! Abort if the level is high enough
        if (level_ >= this%abort_level) then
            stop 42
        end if

        !$omp end critical
    end subroutine
!******************************************************************************!
!> Initialise a `[[FConsole(type)]]` object.
!******************************************************************************!
    subroutine console_init(this, level, abort_level, colour_mode)
        class(FConsole), intent(out) :: this
        integer, intent(in) :: level
        integer, intent(in), optional :: abort_level
        logical, intent(in), optional :: colour_mode
!------
        this%print_level = level
        if (present(abort_level)) this%abort_level = abort_level
        if (present(colour_mode)) then
            this%colour_mode = colour_mode
        else
            this%colour_mode = (isatty(1) == 1)
        end if
    end subroutine
!******************************************************************************!
!> Set the verbosity level of a console object.
!******************************************************************************!
    subroutine setVerbosity(this, level)
        class(FConsole), intent(inout) :: this
        integer, intent(in) :: level
!------
        select case(level)
            case (EVENT_LEVEL_DEBUG, &
                EVENT_LEVEL_LOG, &
                EVENT_LEVEL_VERBOSE, &
                EVENT_LEVEL_QUIET, &
                EVENT_LEVEL_WARNING, &
                EVENT_LEVEL_TRACEBACK, &
                EVENT_LEVEL_ERROR)

                ! Set the level in the console object
                this%print_level = level
            case default

                ! Invalid level code
        end select
    end subroutine
end module
