module dk_structure

    use dk_datatable,  only: DataTable
    use dk_parameters, only: e

    implicit none(external, type)
    private

    integer, parameter, public :: TAG_LENGTH = 8

    type, public :: Structure
        character(:), allocatable :: name
        character(:), allocatable :: description
        real(e), dimension(3,3) :: box = 0.0_e
        real(e), dimension(:,:), allocatable :: r
        real(e), dimension(:,:), allocatable :: v
        real(e), dimension(:,:), allocatable :: f
        character(TAG_LENGTH), dimension(:), allocatable :: tag
        type(DataTable) :: aux
    end type

end module
