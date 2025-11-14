!set******************************************************************************!
!                             dk_exception module
!------------------------------------------------------------------------------!
!> `[[FException(type)]]` derived type and type-bound procedures, for basic
!> error reporting.
!
!  Part of the dk_core library version 0.1.
!  Written by Paul Fossati, <paul.fossati@cea.fr>
!  Copyright (c) 2023-2024 CEA
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
module dk_exception

    use event_levels

    use dk_event, only: FEvent
    use dk_console, only: FConsole, default_console => console

    implicit none (external, type)
    private
    save

!******************************************************************************!
! Event codes from mod_event
!******************************************************************************!
    public :: EVENT_LEVEL_NONE
    public :: EVENT_LEVEL_DEBUG
    public :: EVENT_LEVEL_LOG
    public :: EVENT_LEVEL_VERBOSE
    public :: EVENT_LEVEL_QUIET
    public :: EVENT_LEVEL_WARNING
    public :: EVENT_LEVEL_TRACEBACK
    public :: EVENT_LEVEL_ERROR
    public :: EVENT_LEVEL_IPE
    public :: EVENT_MAX_LEVEL

!******************************************************************************!
! Message prefixes
!******************************************************************************!
    character(*), parameter :: WARNING_PREFIX = char(27)//"[1;31mWarning:"//char(27)//"[0m "
    character(*), parameter :: WARNING_PREFIX_PLAIN = "Warning: "
!    "Warning: "
    character(*), parameter :: ERROR_PREFIX   = char(27)//"[1;31mError:"//char(27)//"[0m "
    character(*), parameter :: ERROR_PREFIX_PLAIN   = "Error: "
!    "Error: "
    character(*), parameter :: LOG_PREFIX     = char(27)//"[1;32mLog"//char(27)//"[0m "
    character(*), parameter :: LOG_PREFIX_PLAIN     = "Log: "
!    "Log: "

    type, public :: FExceptionDescription
        private
        character(:), allocatable :: message
        character(:), allocatable :: context
        integer, public :: level_ = EVENT_LEVEL_ERROR
    contains
        procedure, public :: string => exceptionDescription_string
        procedure, public :: report => exceptionDescription_report
    end type
    interface FExceptionDescription
        module procedure exceptionDescription_init
    end interface

!******************************************************************************!
!> Runtime error description
!******************************************************************************!
    type, public :: FException
        private
        class(FExceptionDescription), allocatable, public :: description
        logical :: show = .false.             !< Whether the error has to be shown
    contains
        ! Exception setup
        procedure, private :: exception_raise                !<
        procedure, private :: exception_raiseWithMessage     !<
        procedure, private :: exception_raiseWithDescription !<
        generic, public :: raise => exception_raise, exception_raiseWithMessage, &
            exception_raiseWithDescription

        ! Exception reports
        procedure, public :: report           !< Show the error message
        procedure, public :: discard          !< Note that the error has been handled and should not be automatically shown
        procedure, public :: transfer         !< Transfer an error status from an exception object
        procedure, public :: string           !< Get a character string describing an exception
        procedure, public :: setLevel         !< Change the priority level of an exception
        procedure, public :: level            !< Get the priority level of an exception

        final :: exception_finalise           !< Finalise the eror objects

        ! Exceptions comparison
        procedure, private :: exception_integerDiffer
        procedure, private :: exception_errorDiffer
        procedure, private :: exception_characterDiffer
        generic, public :: operator(/=) => exception_integerDiffer, exception_errorDiffer, exception_characterDiffer
        procedure, private :: exception_integerEqual
        procedure, private :: exception_errorEqual
        procedure, private :: exception_characterEqual
        generic, public :: operator(==) => exception_integerEqual, exception_errorEqual, exception_characterEqual
    end type

    public :: FError
    public :: FWarning
contains
!******************************************************************************!
!> Change the level of an error description object.
!******************************************************************************!
    subroutine setLevel(this, level)
        class(FException), intent(inout) :: this
        integer, intent(in) :: level
!------ 
        if (allocated(this%description)) then
            this%description%level_ = level
        end if
    end subroutine
!******************************************************************************!
!> Initialise an error description object.
!******************************************************************************!
    pure function exceptionDescription_init(message, context, level) result(this)
        character(*), intent(in), optional :: message
        character(*), intent(in), optional :: context
        integer, intent(in), optional :: level
        type(FExceptionDescription) :: this
!------
        if (present(message)) this%message = message
        if (present(context)) this%context = context
        if (present(level)) this%level_ = level
    end function
!******************************************************************************!
!> Get the message from an error description object.
!******************************************************************************!
    pure function exceptionDescription_string(this) result(string)
        class(FExceptionDescription), intent(in) :: this
        character(:), allocatable :: string
!------
        if (allocated(this%message)) then
            string = this%message
        else
            string = ""
        end if
    end function
!******************************************************************************!
!> Get the priority level of an exception.
!******************************************************************************!
    pure function level(this)
        class(FException), intent(in) :: this
        integer :: level
!------
        if (allocated(this%description)) then
            level = this%description%level_
        else
            ! If the description is not set, assume it is an error
            level = EVENT_LEVEL_ERROR
        end if
    end function
!******************************************************************************!
!> Create an exception description object representing an error condition.
!******************************************************************************!
    pure function FError(message, context) result(desc)
        type(FExceptionDescription) :: desc
        character(*), intent(in) :: message
        character(*), intent(in), optional :: context
!------
        desc%message = message
        desc%level_ = EVENT_LEVEL_ERROR
        if (present(context)) desc%context = context
    end function
!******************************************************************************!
!> Create an exception description object representing a warning condition.
!******************************************************************************!
    pure function FWarning(message, context) result(desc)
        type(FExceptionDescription) :: desc
        character(*), intent(in) :: message
        character(*), intent(in), optional :: context
!------
        desc%message = message
        desc%level_ = EVENT_LEVEL_WARNING
        if (present(context)) desc%context = context
    end function
!******************************************************************************!
!> Raise an exception with a description from an already existing exception
!> description object.
!******************************************************************************!
    subroutine exception_raiseWithDescription(this, description)
        class(FException), intent(inout) :: this
        class(FExceptionDescription), intent(in) :: description
!------
        ! Use the provided description object
        if (allocated(this%description)) deallocate(this%description)
        allocate(this%description, source=description)
        this%show = .true.
    end subroutine
!******************************************************************************!
!> Raise an exception with a default message and a user-provided error level.
!******************************************************************************!
    subroutine exception_raise(this, level)
        class(FException), intent(inout) :: this
        integer, intent(in), optional :: level
!------
        if (allocated(this%description)) deallocate(this%description)
        allocate(FExceptionDescription :: this%description)

        ! Use a default message
        this%description%message = "Unknown error."

        ! Set the loging level if a level is provided
        if (present(level)) then
            this%description%level_ = level
        else
            ! Otherwise, asssume hard error
            this%description%level_ = EVENT_LEVEL_ERROR
        end if

        ! By default, the error has to be signalled
        this%show = .true.
    end subroutine
!******************************************************************************!
!> Raise an exception with use-provided message and a error level.
!******************************************************************************!
    subroutine exception_raiseWithMessage(this, message, level)
        class(FException), intent(inout) :: this
        character(*), intent(in) :: message
        integer, intent(in), optional :: level
!------
        if (allocated(this%description)) deallocate(this%description)
        allocate(FExceptionDescription :: this%description)

        ! Set the error message
        this%description%message = message

        ! Set the loging level if a level is provided
        if (present(level)) then
            this%description%level_ = level
        else
            ! Otherwise, asssume hard error
            this%description%level_ = EVENT_LEVEL_ERROR
        end if

        ! By default, the error has to be signalled
        this%show = .true.
    end subroutine
!******************************************************************************!
!> Report an exception to the user. By default, an error is reported when it is
!> finalised, unless it has been discarded explicitly.
!******************************************************************************!
!LCOV_IGNORE_START
    subroutine exceptionDescription_report(this, console)
        class(FExceptionDescription), intent(in) :: this
        type(FConsole), intent(in) :: console
!------
        ! Do nothing if there is no message
        if (.not.allocated(this%message)) return

        ! Print the message
        if (allocated(this%context)) then
            call console%print( &
                message=this%message, &
                context=this%context, &
                level=this%level_, &
                showprefix=.true. &
            )
        else
            call console%print( &
                message=this%message, &
                level=this%level_, &
                showprefix=.true. &
            )
        end if
    end subroutine
!LCOV_IGNORE_STOP
!******************************************************************************!
!> Report an exception to the user. Simply delegates to the exception
!> description object.
!******************************************************************************!
    subroutine report(this, console)
        class(FException), intent(in) :: this
        type(FConsole), intent(in), optional :: console
!------
        if (allocated(this%description)) then
            if (present(console)) then
                call this%description%report(console)
            else
                call this%description%report(default_console)
            end if
        end if
    end subroutine
!******************************************************************************!
!> Discard an exception. It will not be shown during finalisation. This is the
!> way of ignoring exception once they have been dealt with.
!******************************************************************************!
    subroutine discard(this)
        class(FException), intent(inout) :: this
!------
        if (allocated(this%description)) deallocate(this%description)
        this%show = .false.
    end subroutine
!******************************************************************************!
!> Finalizer for exception objects. The exception will be reported if it has
!> not been discarded.
!******************************************************************************!
    elemental impure subroutine exception_finalise(this)
        type(FException), intent(inout) :: this
!------
        ! Show the error if it has not been taken care of already
        if (allocated(this%description)) then
            if (this%show .and. this%description%level_ /= EVENT_LEVEL_NONE) then
                call this%report
            end if
        end if
    end subroutine
!******************************************************************************!
!> Transfer an error status from an exception object. The original exception is
!> discarded
!******************************************************************************!
    subroutine transfer(this, from)
        class(FException), intent(inout) :: this
        class(FException), intent(inout) :: from
!------
        if (allocated(from%description)) &
            call move_alloc(from=from%description, to=this%description)
        this%show = from%show
        call from%discard
    end subroutine
!******************************************************************************!
!> Compare two exception objects. They are identical if they have the same
!> message, or if both of them are not set.
!******************************************************************************!
    elemental function exception_errorEqual(err1, err2) result(identical)
        class(FException), intent(in) :: err1
        class(FException), intent(in) :: err2
        logical :: identical
!------
        if (allocated(err1%description) .and. allocated(err2%description)) then
            ! Both messages have been allocated
            identical = (err1%description%message == err2%description%message)
        else if (.not.allocated(err1%description) .and. .not.allocated(err2%description)) then
            ! Both messages are unallocated
            identical = .true.
        else
            ! One of the messages is not allocated
            identical = .false.
        end if
    end function
!******************************************************************************!
!> Compare two exception objects. They are different if they have the different
!> messages, or if one of them is not set.
!******************************************************************************!
    elemental function exception_errorDiffer(err1, err2) result(different)
        class(FException), intent(in) :: err1
        class(FException), intent(in) :: err2
        logical :: different
!------
        if (allocated(err1%description) .and. allocated(err2%description)) then
            ! Both messages have been allocated
            different = (err1%description%message /= err2%description%message)
        else if (.not.allocated(err1%description) .and. .not.allocated(err2%description)) then
            ! Both messages are unallocated
            different = .false.
        else
            ! One of the messages is not allocated
            different = .true.
        end if
    end function
!******************************************************************************!
!> Compare an exception object with a message in a character scalar.
!******************************************************************************!
    elemental function exception_characterEqual(err1, message) result(identical)
        class(FException), intent(in) :: err1
        character(*), intent(in) :: message
        logical :: identical
!------
        if (allocated(err1%description)) then
            identical = (err1%description%message == message)
        else
            identical = .false.
        end if
    end function
!******************************************************************************!
!> Compare an exception object with a message in a character scalar.
!******************************************************************************!
    elemental function exception_characterDiffer(err1, message) result(different)
        class(FException), intent(in) :: err1
        character(*), intent(in) :: message
        logical :: different
!------
        if (allocated(err1%description)) then
            different = (err1%description%message /= message)
        else
            different = .true.
        end if
    end function
!******************************************************************************!
!> Check the error status of an exception object. An exception is equals to
!> zero if `[[raise(subroutine)]]` has not been called. This is a way of
!> supporting classic error checks using a sytax similar to that used for
!> integer error codes.
!******************************************************************************!
    elemental function exception_integerEqual(err, stat) result(identical)
        class(FException), intent(in) :: err
        integer, intent(in) :: stat
        logical :: identical
!------
        identical = (allocated(err%description) .eqv. (stat /= 0))
    end function
!******************************************************************************!
!> Check the error status of an exception object. An exception is different
!> from zero if `[[raise(subroutine)]]` has been called. This is a way of
!> supporting classic error checks using a sytax similar to that used for
!> integer error codes.
!******************************************************************************!
    elemental function exception_integerDiffer(err, stat) result(different)
        class(FException), intent(in) :: err
        integer, intent(in) :: stat
        logical :: different
!------
        different = (allocated(err%description) .eqv. (stat == 0))
    end function
!******************************************************************************!
!> Returns the message from an error object. The string is empty if the message
!> in the error is not allocated (meaning the error has not been set using the
!> raise subroutine).
!******************************************************************************!
    pure function string(this)
        class(FException), intent(in) :: this
        character(:), allocatable :: string
!------
        if (allocated(this%description)) then
            string = this%description%string()
        else
            string = ""
        end if
    end function
end module
