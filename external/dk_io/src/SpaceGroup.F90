module dk_spacegroup
    use iso_c_binding, only: c_int, c_double
    use dk_parameters, only: e
    use ext_character, only: to_upper, to_lower
#ifdef HAVE_SPGLIB
    use spglib_f08
#endif
    implicit none(external, type)
    private

    public :: hall_number
    public :: get_equivalent_positions
    public :: get_equivalent_positions_multiple

    enum, bind(C)
    enumerator :: CENTERING_ERROR = 0, &
        PRIMITIVE, &
        BODY, &
        FACE, &
        A_FACE, &
        B_FACE, &
        C_FACE, &
        BASE, &
        R_CENTER
    end enum

    type, private :: SpaceGroup
        integer(c_int) :: number
        character(len=7) :: schoenflies
        character(len=17) :: hall_symbol
        character(len=32) :: international
        character(len=20) :: international_full
        character(len=11) :: international_short
        character(len=6) :: choice
        integer(c_int) :: centering
        integer(c_int) :: pointgroup_number
    end type

    include 'space_group_list.inc'

    integer, parameter :: TAG_LENGTH = 8
contains

    function hall_number(space_group)
        character(*), intent(in) :: space_group
        integer :: hall_number
!------
        character(:), allocatable :: space_group_
        integer :: i, j
!------
        ! Remove any space in the space group name
        if (index(trim(space_group), " ") == 0) then
            space_group_ = space_group
        else
            allocate(character(len_trim(space_group)) :: space_group_)
            space_group_(:) = " "
            j = 0
            do i=1, len(space_group)
                if (space_group(i:i) /= " ") then
                    j = j + 1
                    space_group_(j:j) = space_group(i:i)
                end if
            end do
        end if
        call to_upper(space_group_(1:1))
        call to_lower(space_group_(2:))

        ! Get the Hall number for the space group
        hall_number = 0
        do i=1, size(SPACE_GROUPS)
            if (SPACE_GROUPS(i)%international_short == space_group_) then
                hall_number = i
                exit
            end if
        end do
    end function

    subroutine get_equivalent_positions_multiple(space_group, positions, tags)
#ifdef HAVE_SPGLIB
        use spglib_f08, only: spg_get_symmetry_from_database
#endif
        character(*), intent(in) :: space_group
        real(e), dimension(:,:), allocatable, intent(inout) :: positions
        character(TAG_LENGTH), dimension(:), allocatable, intent(inout) :: tags
#ifdef HAVE_SPGLIB
        integer :: number, i, j, n, num_symmetries, num_sites
        integer(c_int), dimension(:,:,:), allocatable :: rotation
        real(c_double), dimension(:,:), allocatable :: translation
        logical, dimension(:), allocatable :: mask
        real(e), dimension(:,:), allocatable :: equivalent_positions
        character(TAG_LENGTH), dimension(:), allocatable :: equivalent_tags
!------
        number = hall_number(space_group)

        if (number == 0) return

        num_sites = size(positions, 2)

        ! Get the symmetry operations for the space group
        allocate(rotation(3,3,192))
        allocate(translation(3,192))
        num_symmetries = spg_get_symmetry_from_database(rotation, translation, number)

        ! Apply all the symmetry operations
        allocate(equivalent_positions(3,num_symmetries*num_sites))
        allocate(equivalent_tags(num_symmetries*num_sites))
        n = 0
        do j=1, num_sites
            do i=1, num_symmetries
                n = n + 1
                equivalent_positions(:,n) = matmul(transpose(rotation(:,:,i)), positions(:,j)) + translation(:,i)
                equivalent_tags(n) = tags(j)
            end do
        end do

        ! Enforce periodic boundary conditions
        where(equivalent_positions < 0)
            equivalent_positions = equivalent_positions + 1
        end where
        where(equivalent_positions >= 1)
            equivalent_positions = equivalent_positions - 1
        end where

        ! Detect identical positions
        allocate(mask(size(equivalent_positions, 2)))
        mask = .true.
        do i=1, size(equivalent_positions, 2)
            if (.not.mask(i)) cycle

            do j=i+1, size(equivalent_positions, 2)
                if (.not.mask(j)) cycle
                if (norm2(equivalent_positions(:,i)-equivalent_positions(:,j)) < 1.0e-6_e) then
                    mask(j) = .false.
                end if
            end do
        end do

        ! Remove identical positions
        deallocate(positions)
        allocate(positions(3,count(mask)))
        deallocate(tags)
        allocate(tags(count(mask)))
        n = 0
        do i=1, size(equivalent_positions, 2)
            if (.not.mask(i)) cycle
            n = n + 1
            positions(:,n) = equivalent_positions(:,i)
            tags(n) = equivalent_tags(i)
        end do
#endif
    end subroutine

    subroutine get_equivalent_positions(space_group, position, equivalent)
#ifdef HAVE_SPGLIB
        use spglib_f08, only: spg_get_symmetry_from_database
#endif
        character(*), intent(in) :: space_group
        real(e), dimension(3), intent(in) :: position
        real(e), dimension(:,:), allocatable, intent(out) :: equivalent
#ifdef HAVE_SPGLIB
        integer :: number, i, j, n, num_symmetries
        integer(c_int), dimension(:,:,:), allocatable :: rotation
        real(c_double), dimension(:,:), allocatable :: translation
        logical, dimension(:), allocatable :: mask
        real(e), dimension(:,:), allocatable :: tmp
!------
        number = hall_number(space_group)

        if (number == 0) return

        ! Get the symmetry operations for the space group
        allocate(rotation(3,3,192))
        allocate(translation(3,192))
        num_symmetries = spg_get_symmetry_from_database(rotation, translation, number)

        ! Apply all the symmetry operations
        allocate(equivalent(3,num_symmetries))
        do i=1, num_symmetries
            equivalent(:,i) = matmul(transpose(rotation(:,:,i)), position) + translation(:,i)
        end do

        ! Enforce periodic boundary conditions
        where(equivalent < 0)
            equivalent = equivalent + 1
        end where
        where(equivalent >= 1)
            equivalent = equivalent - 1
        end where

        ! Detect identical positions
        allocate(mask(size(equivalent, 2)))
        mask = .true.
        do i=1, size(equivalent, 2)
            if (.not.mask(i)) cycle

            do j=i+1, size(equivalent, 2)
                if (.not.mask(j)) cycle
                if (norm2(equivalent(:,i)-equivalent(:,j)) < 1.0e-6_e) then
                    mask(j) = .false.
                end if
            end do
        end do

        ! Remove identical positions
        allocate(tmp(3,count(mask)))
        n = 0
        do i=1, size(equivalent, 2)
            if (.not.mask(i)) cycle
            n = n + 1
            tmp(:,n) = equivalent(:,i)
        end do

        deallocate(equivalent)
        call move_alloc(from=tmp, to=equivalent)
#endif
    end subroutine
end module
