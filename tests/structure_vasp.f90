program structure_vasp
    use iso_fortran_env, only: int64
    use dk_structure_io
    use dk_math, only: approx_equal
    use dk_fileobject, only: FileObject
    use dk_exception, only: FException
    use dk_test

    implicit none(external, type)

    integer, parameter :: e = kind(0.0d0)
    character(*), parameter :: filename = "POSCAR"
    character(*), parameter :: format = "vasp"
    character(*), parameter :: format_check = "vasp"

    call test_basic_uncompressed
    call test_basic_gzip
    call test_basic_bzip2
    call test_basic_zstandard

    call test_errors
    call test_features

    call tests_finished
contains
!******************************************************************************!
!> Test error conditions in the specific implementation of write_structure.
!******************************************************************************!
    subroutine test_errors
        real(e), dimension(3,3) :: box
        real(e), dimension(:,:), allocatable :: positions, velocities
        real(e), dimension(:), allocatable :: masses
        character(TAG_LENGTH), dimension(:), allocatable :: tags
        type(FException) :: stat
!------
        call get_test_structure(box, positions, tags)

        allocate(masses(size(tags)))
        allocate(velocities(3,size(tags)))

        call start_test_group("Error conditions")

        call write_structure("structure.poscar", box, positions, tags(:4), stat=stat)
        call check(stat=="The tags array has the wrong size.", "size of the tags array")
        call stat%discard
    end subroutine

    subroutine test_features
        real(e), dimension(3,3) :: box!, box2
        real(e), dimension(:,:), allocatable :: positions, velocities! positions2
        real(e), dimension(:), allocatable :: masses
        character(TAG_LENGTH), dimension(:), allocatable :: tags!, tags2
        type(FException) :: stat

        call get_test_structure(box, positions, tags)

        allocate(masses(size(tags)))
        masses = 42
        allocate(velocities(3,size(tags)))
        velocities = 0.42

        call start_test_group("Features")

        call write_structure("structure.poscar", box, positions, tags, velocities=velocities, stat=stat)
        call check(stat==0, "saving velocities")
        call stat%discard

        ! call write_structure("structure.poscar", box, positions, tags, masses=masses, stat=stat)
        ! write(*,*) "Saved masses"
        ! call stat%report
        ! call check(stat==0, "saving masses")
        ! call stat%discard

        call write_structure("structure.poscar", box, positions, tags, forces=velocities, stat=stat)
        call check(stat==0, "saving forces")
        call stat%discard
    end subroutine

include 'structure_basic_checks.inc'
end program

