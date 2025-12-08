program structure_dlpoly
    use iso_fortran_env, only: int64
    use dk_structure_io
    use dk_math, only: approx_equal
    use dk_fileobject, only: FileObject
    use dk_exception, only: FException, EVENT_LEVEL_QUIET
    use dk_test

    implicit none(external, type)

    integer, parameter :: e = kind(0.0d0)
    character(*), parameter :: filename = "CONFIG"
    character(*), parameter :: format = "dlpoly"
    character(*), parameter :: format_check = "dlpoly"

    call test_basic_uncompressed
    call test_basic_gzip
!     call test_basic_bzip2
    call test_basic_zstandard

    call test_errors
    call test_features

    call tests_finished
contains
    subroutine test_features
        real(e), dimension(3,3) :: box
        real(e), dimension(:,:), allocatable :: positions, velocities, forces
        real(e), dimension(3,3) :: box2
        real(e), dimension(:,:), allocatable :: positions2, velocities2, forces2
        character(TAG_LENGTH), dimension(:), allocatable :: tags
        character(TAG_LENGTH), dimension(:), allocatable :: tags2
        type(FException) :: stat
        integer :: ret
        character(:), allocatable :: msg

        call start_test_group("CONFIG file features")

        call get_test_structure(box, positions, tags, velocities, forces)

!------ Positions only
        call write_structure(filename, box, positions, tags, format=format, stat=stat)
        call check(stat==0, "file output with positions only", "Message: "//stat%string())
        call stat%discard

        call read_structure(filename, box2, positions2, tags2, stat=ret)
        call check(all(approx_equal(box, box2)), "simulation box")
        call check(allocated(positions2), "positions array allocation status")
        call check(size(positions)==size(positions2), "positions array size")
        call check(all(approx_equal(positions, positions2)), "positions array values")

!------ Positions and velocities
        call write_structure(filename, box, positions, tags, velocities=velocities, format=format, stat=stat)
        call check(stat==0, "file output with positions and velocities", "Message: "//stat%string())
        call stat%discard

        call read_structure(filename, box2, positions2, tags2, velocities=velocities2, format=format, stat=ret)
        call check(all(approx_equal(box, box2)), "simulation box")
        call check(allocated(positions2), "positions array allocation status")
        call check(size(positions)==size(positions2), "positions array size")
        call check(all(approx_equal(positions, positions2)), "positions array values")

        call check(allocated(velocities2), "velocities array allocation status")
        call check(size(velocities)==size(velocities2), "velocities array size")
        call check(all(approx_equal(velocities, velocities2)), "velocities array values")

!------ Positions and forces
        call write_structure(filename, box, positions, tags, forces=forces, format=format, stat=stat)
        call check(stat==0, "file output with positions and forces", "Message: "//stat%string())
        call stat%discard

        call read_structure(filename, box2, positions2, tags2, forces=forces2, format=format, stat=ret)
        call check(all(approx_equal(box, box2)), "simulation box")
        call check(allocated(positions2), "positions array allocation status")
        call check(size(positions)==size(positions2), "positions array size")
        call check(all(approx_equal(positions, positions2)), "positions array values")

        call check(allocated(forces2), "forces array allocation status")
        call check(size(forces)==size(forces2), "forces array size")
        call check(all(approx_equal(forces, forces2)), "forces array values")

!------ Positions, velocities, and forces
        call write_structure(filename, box, positions, tags, velocities=velocities, forces=forces, format=format, stat=stat)
        call check(stat==0, "file output with positions, velocities, and forces", "Message: "//stat%string())
        call stat%discard

        call read_structure(filename, box2, positions2, tags2,  velocities=velocities2,  forces=forces2, format=format, stat=ret)
        call check(all(approx_equal(box, box2)), "simulation box")
        call check(allocated(positions2), "positions array allocation status")
        call check(size(positions)==size(positions2), "positions array size")
        call check(all(approx_equal(positions, positions2)), "positions array values")

        call check(allocated(velocities2), "velocities array allocation status")
        call check(size(velocities)==size(velocities2), "velocities array size")
        call check(all(approx_equal(velocities, velocities2)), "velocities array values")

        call check(allocated(forces2), "forces array allocation status")
        call check(size(forces)==size(forces2), "forces array size")
        call check(all(approx_equal(forces, forces2)), "forces array values")

!------ Reading HISTORY files
        call read_structure("HISTORY", box, positions, tags, msg=msg)
        call check(msg=="", "reading HISTORY files", "Message: "//msg)

    end subroutine
!******************************************************************************!
!> Test error conditions in the specific implementation of write_structure.
!******************************************************************************!
    subroutine test_errors
        real(e), dimension(3,3) :: box
        real(e), dimension(:,:), allocatable :: positions, velocities, forces
        logical, dimension(:), allocatable :: include
        character(TAG_LENGTH), dimension(:), allocatable :: tags
        type(FException) :: stat
        integer :: ret
        character(:), allocatable :: msg
!------
        call start_test_group("error conditions")

        call get_test_structure(box, positions, tags)
        deallocate(positions)
        allocate(positions(3,0))
        deallocate(tags)
        allocate(tags(0))

        call write_structure(filename, box, positions, tags, format="dlpoly", stat=stat)
        call check(stat==0, "writing empty structures", "Message: "//stat%string())
        call stat%discard

        call get_test_structure(box, positions, tags, velocities)
        allocate(include(size(tags)))

        call write_structure(filename, box, positions, tags(:4), format="dlpoly", stat=stat)
        call check(stat=="The tags array has the wrong size.", "size of the tags array", "Message: "//stat%string())
        call stat%discard

        call write_structure(filename, box, positions, tags, velocities=velocities(:,:4), format="dlpoly", stat=stat)
        call check(stat=="The velocities array has the wrong size.", "size of the velocities array", "Message: "//stat%string())
        call stat%discard

        call write_structure(filename, box, positions, tags, forces=velocities(:,:4), format="dlpoly", stat=stat)
        call check(stat=="The forces array has the wrong size.", "size of the forces array", "Message: "//stat%string())
        call stat%discard

!         call write_structure(filename, box, positions, tags, include=include(:4), format="dlpoly", stat=stat)
!         call check(stat=="The include mask has the wrong size.", "size of the include mask", "Message: "//msg)
!         call stat%discard

        call read_structure("bad_1.dlpoly", box, positions, tags, format="dlpoly", except=stat)
        call check(stat=="Invalid element for cell vector a.", "Invalid box vector a", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_2.dlpoly", box, positions, tags, format="dlpoly", except=stat)
        call check(stat=="Invalid element for cell vector b.", "Invalid box vector b", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_3.dlpoly", box, positions, tags, format="dlpoly", except=stat)
        call check(stat=="Invalid element for cell vector c.", "Invalid box vector c", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_4.dlpoly", box, positions, tags, format="dlpoly", except=stat)
        call check(stat=="Invalid atom identifier.", "Invalid atom identifier", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_5.dlpoly", box, positions, tags, format="dlpoly", except=stat)
        call check(stat=="Invalid atomic position.", "Invalid atom position", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_6.dlpoly", box, positions, tags, velocities=velocities, format="dlpoly", except=stat)
        call check(stat=="Invalid atomic velocity.", "Invalid atom velocity", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_7.dlpoly", box, positions, tags, forces=forces, format="dlpoly", except=stat)
        call check(stat=="Invalid atomic force.", "Invalid atom force", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_8.dlpoly", box, positions, tags, forces=forces, format="dlpoly", except=stat)
        call check(stat=="Unexpected end of file.", "Truncated file 1", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_9.dlpoly", box, positions, tags, forces=forces, format="dlpoly", except=stat)
        call check(stat=="Unexpected end of file.", "Truncated file 2", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_10.dlpoly", box, positions, tags, forces=forces, format="dlpoly", except=stat)
        call check(stat=="Unexpected end of file.", "Truncated file 3", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_11.dlpoly", box, positions, tags, forces=forces, format="dlpoly", except=stat)
        call check(stat=="Unexpected end of file.", "Truncated file 4", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_12.dlpoly", box, positions, tags, forces=forces, format="dlpoly", except=stat)
        call check(stat=="Expected atom tag.", "Truncated file 5", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

    end subroutine

include 'structure_basic_checks.inc'
end program
