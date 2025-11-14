
module event_levels
    implicit none(external, type)
    public
!******************************************************************************!
! Importance levels for an event
!******************************************************************************!
    integer, parameter :: EVENT_LEVEL_NONE      = 0   !< Non-event
    integer, parameter :: EVENT_LEVEL_DEBUG     = 1   !< Event relevant for debugging purposes only
    integer, parameter :: EVENT_LEVEL_LOG       = 2   !< Event to be logged
    integer, parameter :: EVENT_LEVEL_VERBOSE   = 3   !< Common event, should be printed
                                                      !< most of the time
    integer, parameter :: EVENT_LEVEL_QUIET     = 4   !< Important event, should be
                                                      !< printed in quiet contexts
    integer, parameter :: EVENT_LEVEL_WARNING   = 5   !< Event justifying printing a warning
    integer, parameter :: EVENT_LEVEL_TRACEBACK = 6   !< Tracebacks, used for internal errors
    integer, parameter :: EVENT_LEVEL_ERROR     = 7   !< Event justifying printing an error message and aborting
    integer, parameter :: EVENT_LEVEL_IPE       = 8   !< Event justifying printing an error message,
                                                      !< additional debugging information, and aborting
    integer, parameter :: EVENT_MAX_LEVEL       = 8   !< Highest possible level
end module

module dk_event
    use iso_c_binding, only: c_size_t

    use core, only: e, long

    use event_levels

    implicit none (external, type)
    private

    integer, parameter :: lkind = c_size_t
!******************************************************************************!
!> Event object
!******************************************************************************!
    type, public :: FEvent
        class(FEvent), pointer :: next => null()
        class(FEvent), pointer :: last => null()
        class(FEvent), pointer :: previous => null()
        integer(long) :: timestamp = 0
        real(e) :: cpu_time = 0
        integer :: level = EVENT_LEVEL_NONE
    contains
        procedure, public :: setup
        procedure, public :: getLast
        procedure, public :: getMessage
    end type
!******************************************************************************!
!> Event object describing a PROCEDURE-START event.
!******************************************************************************!
    type, public, extends(FEvent) :: FProcedureStartedEvent
        character(:), allocatable :: name
        character(:), allocatable :: file
        integer :: line = 0
        class(FProcedureStartedEvent), pointer :: caller => null()
    contains
        procedure, public :: getMessage => procedureStartedEvent_getMessage
        procedure, public :: setup => procedureStartedEvent_setup
    end type
    interface FProcedureStartedEvent
        module procedure procedureStartedEvent_init
    end interface
!******************************************************************************!
!> Event object describing a PROCEDURE-END event.
!******************************************************************************!
    type, public, extends(FEvent) :: FProcedureEndedEvent
        character(:), allocatable :: name
        character(:), allocatable :: file
        integer :: line = 0
        class(FProcedureStartedEvent), pointer :: caller => null()
        class(FProcedureStartedEvent), pointer :: starter => null()
    contains
        procedure, public :: getMessage => procedureEndedEvent_getMessage
        procedure, public :: setup => procedureEndedEvent_setup
    end type
    interface FProcedureEndedEvent
        module procedure procedureEndedEvent_init
    end interface

    type, public, extends(FEvent) :: FSignpostBeginEvent
        character(:), allocatable :: key
    end type
    interface FSignpostBeginEvent
        module procedure signpostBeginEvent_init
    end interface

    type, public, extends(FEvent) :: FSignpostEndEvent
        character(:), allocatable :: key
    end type
    interface FSignpostEndEvent
        module procedure signpostEndEvent_init
    end interface

    type, public, extends(FEvent) :: FMemoryAllocEvent
        character(:), allocatable :: name
        integer(lkind) :: size = 0_lkind
        character(:), allocatable :: procname
        character(:), allocatable :: file
        integer :: line = 0
    end type
    interface FMemoryAllocEvent
        module procedure FMemoryAllocEvent_init
    end interface

    type, public, extends(FEvent) :: FMemoryDeallocEvent
        character(:), allocatable :: name
        integer(lkind) :: size = 0_lkind
        character(:), allocatable :: procname
        character(:), allocatable :: file
        integer :: line = 0
    end type
    interface FMemoryDeallocEvent
        module procedure FMemoryDeallocEvent_init
    end interface
contains
!******************************************************************************!
!> Setup a SIGNPOST-BEGIN event.
!******************************************************************************!
    function signpostBeginEvent_init(key) result(event)
        character(*), intent(in) :: key
        type(FSignpostBeginEvent) :: event
!------
        event%key = key
        event%level = EVENT_LEVEL_LOG
    end function
!******************************************************************************!
!> Setup a SIGNPOST-END event.
!******************************************************************************!
    function signpostEndEvent_init(key) result(event)
        character(*), intent(in) :: key
        type(FSignpostEndEvent) :: event
!------
        event%key = key
        event%level = EVENT_LEVEL_LOG
    end function
!******************************************************************************!
!> Setup a MEMORY-ALLOC event.
!******************************************************************************!
    function FMemoryAllocEvent_init(name, size, procname, file, line) result(event)
        character(*), intent(in) :: name
        integer(lkind), intent(in) :: size
        character(*), intent(in) :: procname
        character(*), intent(in) :: file
        integer, intent(in) :: line
        type(FMemoryAllocEvent) :: event
!------
        event%name = name
        event%size = size
        event%procname = procname
        event%file = file
        event%line = line
        event%level = EVENT_LEVEL_DEBUG
    end function
!******************************************************************************!
!> Setup a MEMORY-DEALLOC event.
!******************************************************************************!
    function FMemoryDeallocEvent_init(name, size, procname, file, line) result(event)
        character(*), intent(in) :: name
        integer(lkind), intent(in) :: size
        character(*), intent(in) :: procname
        character(*), intent(in) :: file
        integer, intent(in) :: line
        type(FMemoryDeallocEvent) :: event
!------
        event%name = name
        event%size = size
        event%procname = procname
        event%file = file
        event%line = line
        event%level = EVENT_LEVEL_DEBUG
    end function
!******************************************************************************!
!> Setup a PROCEDURE-START event.
!******************************************************************************!
    function procedureStartedEvent_init(name, file, line) result(event)
        character(*), intent(in) :: name
        character(*), intent(in) :: file
        integer, intent(in) :: line
        type(FProcedureStartedEvent) :: event
!------
        event%name = name
        event%file = file
        event%line = line
        event%level = EVENT_LEVEL_DEBUG
    end function
!******************************************************************************!
!>
!******************************************************************************!
    subroutine procedureStartedEvent_setup(this)
        class(FProcedureStartedEvent), intent(inout) :: this
!------
        class(FEvent), pointer :: v
!------
        ! look for the ProcedureStarted event for the caller
        this%caller => null()
        if (associated(this%previous)) then
            v => this%previous
            callerLoop: do while(associated(v))
                select type(p => v)
                    class is (FProcedureStartedEvent)
                        this%caller => p
                        exit callerLoop
                    class is (FProcedureEndedEvent)
                        v => p%starter
                end select
                if (associated(v)) v => v%previous
            end do callerLoop
        end if
    end subroutine
!******************************************************************************!
!> Get a string representation of a procedureStarted event.
!******************************************************************************!
    subroutine procedureStartedEvent_getMessage(this, message)
        class(FProcedureStartedEvent), intent(in) :: this
        character(:), allocatable, intent(out) :: message
!------
        integer :: s
!------
        ! Print the message
        if (associated(this%caller)) then
            s = len(this%name) + len(this%file) + 84 + nint(log10(this%line*1.0_e))+1 + len(this%caller%name)
            allocate(character(s) :: message)
            write(message, '(a,i0,a)')  '{"event" : "ProcedureStarted", "name" : "'//this%name//'", "file" : "'//this%file//'", "line" : "', &
                this%line, '", "caller" : "'//this%caller%name//'"}'
        else
            s = len(this%name) + len(this%file) + 86 + nint(log10(this%line*1.0_e))+1
            allocate(character(s) :: message)
            write(message, '(a,i0,a)')  '{"event" : "ProcedureStarted", "name" : "'//this%name//'", "file" : "'//this%file//'", "line" : "', &
                this%line, '", "caller" : null}'
        end if
    end subroutine
!******************************************************************************!
!> Setup a procedureEnded event.
!******************************************************************************!
    function procedureEndedEvent_init(name, file, line) result(event)
        character(*), intent(in) :: name
        character(*), intent(in) :: file
        integer, intent(in) :: line
        type(FProcedureEndedEvent) :: event
!------
        event%name = name
        event%file = file
        event%line = line
        event%level = EVENT_LEVEL_DEBUG
    end function
!******************************************************************************!
!>
!******************************************************************************!
    subroutine procedureEndedEvent_setup(this)
        class(FProcedureEndedEvent), intent(inout) :: this
!------
        class(FEvent), pointer :: v
!------
        ! look for the ProcedureStarted event for the caller
        this%caller => null()
        if (associated(this%previous)) then
            v => this%previous
            callerLoop: do while(associated(v))
                select type(p => v)
                    class is (FProcedureStartedEvent)
                        if (p%name==this%name .and. p%file==this%file) then
                            this%starter => p
                        else
                            this%caller => p
                            exit callerLoop
                        end if
                    class is (FProcedureEndedEvent)
                        v => p%starter
                end select
                if (associated(v)) v => v%previous
            end do callerLoop
        end if
    end subroutine
!******************************************************************************!
!> Get a string representation of a procedureEnded event.
!******************************************************************************!
    subroutine procedureEndedEvent_getMessage(this, message)
        class(FProcedureEndedEvent), intent(in) :: this
        character(:), allocatable, intent(out) :: message
!------
        integer :: s
!------
        ! Print the message
        if (associated(this%caller)) then
            s = len(this%name) + len(this%file) + 84 + nint(log10(this%line*1.0_e))+1 + len(this%caller%name)
            allocate(character(s) :: message)
            write(message, '(a,i0,a)')  '{"event" : "ProcedureEnded",   "name" : "'//this%name//'", "file" : "'//this%file//'", "line" : "', &
                this%line, '", "caller" : "'//this%caller%name//'"}'
        else
            s = len(this%name) + len(this%file) + 86 + nint(log10(this%line*1.0_e))+1
            allocate(character(s) :: message)
            write(message, '(a,i0,a)')  '{"event" : "ProcedureEnded",   "name" : "'//this%name//'", "file" : "'//this%file//'", "line" : "', &
                this%line, '", "caller" : null}'
        end if
    end subroutine
!******************************************************************************!
!> Get a description of an event.
!******************************************************************************!
    subroutine getMessage(this, message)
        class(FEvent), intent(in) :: this
        character(:), allocatable, intent(out) :: message
!-----
        message = ""
    end subroutine
!******************************************************************************!
!> Finalise the setup of an event object once it is integrated in an event
!> chain.
!******************************************************************************!
    subroutine setup(this)
        class(FEvent), intent(inout) :: this
    end subroutine
!******************************************************************************!
!> Get the last event of an event chain.
!******************************************************************************!
    subroutine getLast(this, last)
        class(FEvent), target, intent(in) :: this
        class(FEvent), pointer, intent(out) :: last
!------
        last => null()
        if (associated(this%last)) then
            last => this%last
        else
            last => this
        end if
    end subroutine
end module
