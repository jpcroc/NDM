program structure_castep
    use iso_fortran_env, only: int64
    use dk_structure_io
    use dk_math, only: approx_equal
    use dk_fileobject, only: FileObject
    use dk_exception, only: FException
    use dk_test

    implicit none(external, type)

    integer, parameter :: e = kind(0.0d0)
    character(*), parameter :: filename = "structure.cell"
    character(*), parameter :: format = "castep"
    character(*), parameter :: format_check = "castep"

    call test_basic_uncompressed
    call test_basic_gzip
    call test_basic_bzip2
    call test_basic_zstandard

!     call test_errors
    call test_castep_format

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
!         call get_test_structure(box, positions, tags)
!
!         allocate(masses(size(tags)))
!         allocate(velocities(3,size(tags)))
!
!         call start_test_group("Error conditions")
!
!         call write_structure("structure.xyz", box, positions, tags(:4), stat=stat)
!         call check(stat=="The tags array has the wrong size.", "size of the tags array")
!         call stat%discard
!
!         call write_structure("structure.xyz", box, positions, tags, velocities=velocities(:,:4), stat=stat)
!         call check(stat=="The velocities array has the wrong size.", "size of the velocities array")
!         call stat%discard
!
!         call write_structure("structure.xyz", box, positions, tags, masses=masses(:4), stat=stat)
!         call check(stat=="The masses array has the wrong size.", "size of the masses array")
!         call stat%discard
!
!         call write_structure("structure.xyz", box, positions, tags, forces=velocities(:,:4), stat=stat)
!         call check(stat=="The forces array has the wrong size.", "size of the forces array")
!         call stat%discard
    end subroutine

    subroutine test_castep_format
        real(e), dimension(3,3) :: box, box2
        real(e), dimension(:,:), allocatable :: positions, velocities, positions2, velocities2
        real(e), dimension(:), allocatable :: masses
        character(TAG_LENGTH), dimension(:), allocatable :: tags, tags2
        character(:), allocatable :: msg
        integer :: stat
        type(FException) :: ex

        call get_test_structure(box2, positions2, tags2, velocities2)

        call start_test_group("Features")

        call write_structure("structure.cell", box2, positions2, tags2, velocities=velocities2, format="castep", stat=ex)
        call check(ex, "Saving velocities")

        call read_structure("structure.cell", box, positions, tags, velocities=velocities, format="castep", stat=stat, msg=msg)
        call check(stat==0, "Reading velocities", "Message: "//msg)
        call check(allocated(velocities), "Velocities array allocation status")
        if (allocated(velocities)) then
            call check(size(velocities)==size(velocities2), "Velocities array size")
            if (size(velocities)==size(velocities2)) &
                call check(all(approx_equal(velocities, velocities2)), "Velocities array values")
        end if

        call read_structure("valid_1.cell", box, positions, tags, format="castep", stat=stat, msg=msg)
        call check(stat==0, "Cell parameters", "Message: "//msg)
        call check(all(approx_equal(box,box2)), "simulation box")

        call read_structure("valid_2.cell", box, positions, tags, format="castep", stat=stat, msg=msg)
        call check(stat==0, "Cartesian positions", "Message: "//msg)
        call check(allocated(positions), "positions array allocation status")
        if (allocated(positions)) then
            call check(size(positions)==size(positions2), "positions array size")
            if (size(positions) == size(positions2)) then
                call check(all(approx_equal(positions, positions2)), "positions array values")
            end if
        end if

        call read_structure("valid_3.cell", box, positions, tags, format="castep", stat=stat, msg=msg)
        call check(stat==0, "Cell parameters in nm", "Message: "//msg)
        call check(all(approx_equal(box,box2)), "simulation box")

        call read_structure("valid_4.cell", box, positions, tags, format="castep", stat=stat, msg=msg)
        call check(stat==0, "Cell parameters in cm", "Message: "//msg)
        call check(all(approx_equal(box,box2)), "simulation box")

        call read_structure("valid_5.cell", box, positions, tags, format="castep", stat=stat, msg=msg)
        call check(stat==0, "Cell parameters in m", "Message: "//msg)
        call check(all(approx_equal(box,box2)), "simulation box")

        call read_structure("bad_1.cell", box, positions, tags, format="castep", stat=stat, msg=msg)
        call check(msg=="Invalid lattice parameter value.", "Wrong lattice parameter", "Message: "//msg)
        call read_structure("bad_2.cell", box, positions, tags, format="castep", stat=stat, msg=msg)
        call check(msg=="Invalid lattice angle value.", "Wrong lattice angle", "Message: "//msg)
        call read_structure("bad_3.cell", box, positions, tags, format="castep", stat=stat, msg=msg)
        call check(msg=="Expected end of LATTICE_ABC block.", "Missing end of block", "Message: "//msg)
        call read_structure("bad_4.cell", box, positions, tags, format="castep", stat=stat, msg=msg)
        call check(msg=="Missing lattice parameters.", "Missing lattice parameters", "Message: "//msg)
        call read_structure("bad_5.cell", box, positions, tags, format="castep", stat=stat, msg=msg)
        call check(msg=="Missing atomic positions.", "Missing atomic positions", "Message: "//msg)

    end subroutine

include 'structure_basic_checks.inc'
end program

