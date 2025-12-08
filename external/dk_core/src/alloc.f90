
module dk_alloc
    implicit none(external, type)
    private

    integer, parameter :: dp = kind(0.0d0)

    interface reallocate
        module procedure reallocate_character_1d
        module procedure reallocate_double_1d
        module procedure reallocate_double_2d
        module procedure reallocate_real_1d
        module procedure reallocate_real_2d
        module procedure reallocate_integer_1d
        module procedure reallocate_integer_2d
    end interface
    public :: reallocate, reallocate_character_1d

contains
    subroutine reallocate_character_1d(array, newsize, minsize)
        character(*), dimension(:), allocatable, intent(inout) :: array
        integer, intent(in) :: newsize
        integer, intent(in), optional :: minsize
!------
        character(:), dimension(:), allocatable :: t
        integer :: i
!------
        if (.not.allocated(array)) then
            allocate(array(newsize))
            return
        end if

        if (size(array) > newsize) then
            allocate(character(len(array)) :: t(newsize))
            t(:) = array(:newsize)
            deallocate(array)
            call move_alloc(from=t, to=array)
            return
        end if

        if (present(minsize)) then
            if (size(array) > minsize) return
        end if
        allocate(character(len(array)) :: t(newsize))
        if (allocated(array)) then
            do i=1, size(array)
                t(i) = array(i)
            end do
            do i=size(array)+1, size(t)
                t(i) = " "
            end do
            deallocate(array)
        else
            do i=1, size(t)
                t(i) = " "
            end do
        end if
        call move_alloc(from=t, to=array)
    end subroutine

    subroutine reallocate_double_1d(array, newsize, minsize)
        real(dp), dimension(:), allocatable, intent(inout) :: array
        integer, intent(in) :: newsize
        integer, intent(in), optional :: minsize
!------
        real(dp), dimension(:), allocatable :: t
!------
        if (size(array) > newsize) then
            allocate(t(newsize))
            t(:) = array(:newsize)
            deallocate(array)
            call move_alloc(from=t, to=array)
            return
        end if

        if (present(minsize)) then
            if (size(array) > minsize) return
        end if
        allocate(t(newsize))
        t = 0
        if (allocated(array)) then
            t(:size(array)) = array
            deallocate(array)
        end if
        call move_alloc(from=t, to=array)
    end subroutine

    subroutine reallocate_real_1d(array, newsize, minsize)
        real, dimension(:), allocatable, intent(inout) :: array
        integer, intent(in) :: newsize
        integer, intent(in), optional :: minsize
!------
        real, dimension(:), allocatable :: t
!------
        if (size(array) > newsize) then
            allocate(t(newsize))
            t(:) = array(:newsize)
            deallocate(array)
            call move_alloc(from=t, to=array)
            return
        end if

        if (present(minsize)) then
            if (size(array) > minsize) return
        end if
        allocate(t(newsize))
        t = 0
        if (allocated(array)) then
            t(:size(array)) = array
            deallocate(array)
        end if
        call move_alloc(from=t, to=array)
    end subroutine

    subroutine reallocate_integer_1d(array, newsize, minsize)
        integer, dimension(:), allocatable, intent(inout) :: array
        integer, intent(in) :: newsize
        integer, intent(in), optional :: minsize
!------
        integer, dimension(:), allocatable :: t
!------
        if (size(array) > newsize) then
            allocate(t(newsize))
            t(:) = array(:newsize)
            deallocate(array)
            call move_alloc(from=t, to=array)
            return
        end if

        if (present(minsize)) then
            if (size(array) > minsize) return
        end if
        allocate(t(newsize))
        t = 0
        if (allocated(array)) then
            t(:size(array)) = array
            deallocate(array)
        end if
        call move_alloc(from=t, to=array)
    end subroutine

    subroutine reallocate_double_2d(array, newsize, minsize)
        real(dp), dimension(:,:), allocatable, intent(inout) :: array
        integer, intent(in), dimension(2) :: newsize
        integer, intent(in), dimension(2), optional :: minsize
!------
        real(dp), dimension(:,:), allocatable :: t
        integer :: a, b
!------
        if (.not.allocated(array)) then
            allocate(array(newsize(1),newsize(2)))
            array = 0
            return
        end if

        if (present(minsize)) then
            if (size(array, 1) > minsize(1) .or. size(array, 2) > minsize(2)) return
        end if
        allocate(t(newsize(1), newsize(2)))
        t = 0
        if (allocated(array)) then
            a = min(size(array, 1), newsize(1))
            b = min(size(array, 2), newsize(2))
            t(:a,:b) = array(:a,:b)
            deallocate(array)
        end if
        call move_alloc(from=t, to=array)
    end subroutine

    subroutine reallocate_real_2d(array, newsize, minsize)
        real, dimension(:,:), allocatable, intent(inout) :: array
        integer, intent(in), dimension(2) :: newsize
        integer, intent(in), dimension(2), optional :: minsize
!------
        real, dimension(:,:), allocatable :: t
        integer :: a, b
!------
        if (.not.allocated(array)) then
            allocate(array(newsize(1),newsize(2)))
            array = 0
            return
        end if

        if (present(minsize)) then
            if (size(array, 1) > minsize(1) .or. size(array, 2) > minsize(2)) return
        end if
        allocate(t(newsize(1), newsize(2)))
        t = 0
        if (allocated(array)) then
            a = min(size(array, 1), newsize(1))
            b = min(size(array, 2), newsize(2))
            t(:a,:b) = array(:a,:b)
            deallocate(array)
        end if
        call move_alloc(from=t, to=array)
    end subroutine

    subroutine reallocate_integer_2d(array, newsize, minsize)
        integer, dimension(:,:), allocatable, intent(inout) :: array
        integer, intent(in), dimension(2) :: newsize
        integer, intent(in), dimension(2), optional :: minsize
!------
        integer, dimension(:,:), allocatable :: t
        integer :: a, b
!------
        if (.not.allocated(array)) then
            allocate(array(newsize(1),newsize(2)))
            array = 0
            return
        end if

        if (present(minsize)) then
            if (size(array, 1) > minsize(1) .or. size(array, 2) > minsize(2)) return
        end if
        allocate(t(newsize(1), newsize(2)))
        t = 0
        if (allocated(array)) then
            a = min(size(array, 1), newsize(1))
            b = min(size(array, 2), newsize(2))
            t(:a,:b) = array(:a,:b)
            deallocate(array)
        end if
        call move_alloc(from=t, to=array)
    end subroutine
end module
