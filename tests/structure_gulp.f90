program structure_gulp
    use iso_fortran_env, only: int64
    use dk_structure_io, only: read_structure, write_structure, TAG_LENGTH, get_file_format, get_test_structure, get_triclinic_test_structure
    use dk_math, only: approx_equal
    use dk_fileobject, only: FileObject
    use dk_exception, only: FException, EVENT_LEVEL_QUIET
    use dk_test

    implicit none(external, type)

    integer, parameter :: e = kind(0.0d0)
    character(*), parameter :: filename = "structure.gin"
    character(*), parameter :: format = "gulp"
    character(*), parameter :: format_check = "gulp"

    call test_basic_uncompressed
    call test_basic_gzip
    call test_basic_bzip2
    call test_basic_zstandard

    call test_errors
    call test_gulp_format
    call test_triclinic

    call tests_finished
contains
    subroutine test_triclinic
        real(e), dimension(3,3) :: box, box2
        real(e), dimension(:,:), allocatable :: positions, positions2
        type(FException) :: stat
        integer :: istat
        character(:), allocatable :: msg
        character(TAG_LENGTH), dimension(:), allocatable :: tags, tags2

        call start_test_group("Triclinic boxes")

        call read_structure("triclinic.gin", box, positions, tags, format="gulp", stat=istat, msg=msg)
        call check(istat==0, "Reading reference triclinic structure", "Message: "//msg)

        call write_structure("triclinic_test.gin", box, positions, tags, format="gulp", stat=stat)
        call check(stat==0, "Writing triclinic structure", stat%string())

        call read_structure("triclinic_test.gin", box2, positions2, tags2, format="gulp", stat=istat, msg=msg)
        call check(istat==0, "Reading triclinic structure", "Message: "//msg)

        call check(all(box==box2), "triclinic cell shape")
        write(*,*) box(:,1)
        write(*,*) box(:,2)
        write(*,*) box(:,3)
        write(*,*)
        write(*,*) box2(:,1)
        write(*,*) box2(:,2)
        write(*,*) box2(:,3)
    end subroutine
!******************************************************************************!
!> Test the detection of varions syntax error in structure files.
!******************************************************************************!
    subroutine test_gulp_format
        real(e), dimension(3,3) :: box, box2
        real(e), dimension(:,:), allocatable :: positions, positions2
        character(TAG_LENGTH), dimension(:), allocatable :: tags, tags2
        integer :: stat
        character(:), allocatable :: msg
        type(FException) :: istat

        call get_test_structure(box2, positions2, tags2)

        call start_test_group("GULP format")

        call read_structure("symmetry.gin", box, positions, tags, format="gulp", stat=stat, msg=msg)
        call check(stat==0, "Symmetrised fluorite structure", "Message: "//msg)
        call write_structure("symmetry.xfg", box, positions, tags)

        call read_structure("valid_1.gin", box, positions, tags, format="gulp", stat=stat, msg=msg)
        call check(stat==0, "No empty line after atomic positions", "Message: "//msg)

        call read_structure("valid_2.gin", box, positions, tags, format="gulp")
        call check(stat==0, "Cartesian coordinates", "Message: "//msg)
        call check(all(approx_equal(box,box2)), "simulation box")
        call check(allocated(positions), "positions array allocation status")
        if (allocated(positions)) then
            call check(size(positions)==size(positions2), "positions array size")
            call check(all(approx_equal(positions, positions2)), "positions array values")
        end if

        call read_structure("bad_2.gin", box, positions, tags, format="gulp", except=istat)
        call check(istat=="Missing cell parameters.", "Missing cell parameters", "Message: "//istat%string())
        call istat%setlevel(EVENT_LEVEL_QUIET)
        call istat%report
        call istat%discard

        call read_structure("bad_3.gin", box, positions, tags, format="gulp", except=istat)
        call check(istat=="Missing atomic positions.", "Missing atomic positions", "Message: "//istat%string())
        call istat%setlevel(EVENT_LEVEL_QUIET)
        call istat%report
        call istat%discard

        call read_structure("bad_value_1.gin", box, positions, tags, format="gulp", except=istat)
        call check(istat=="Invalid lattice parameter value.", "Invalid lattice parameter", "Message: "//istat%string())
        call istat%setlevel(EVENT_LEVEL_QUIET)
        call istat%report
        call istat%discard

        call read_structure("bad_value_2.gin", box, positions, tags, format="gulp", except=istat)
        call check(istat=="Invalid atomic coordinate.", "Invalid atomic position", "Message: "//istat%string())
        call istat%setlevel(EVENT_LEVEL_QUIET)
        call istat%report
        call istat%discard

        call read_structure("bad_value_3.gin", box, positions, tags, format="gulp", except=istat)
        call check(istat=="Unknown space group.", "Invalid space group", "Message: "//istat%string())
        call istat%setlevel(EVENT_LEVEL_QUIET)
        call istat%report
        call istat%discard

        call read_structure("bad_value_4.gin", box, positions, tags, format="gulp", except=istat)
        call check(istat=="Invalid element for cell vector a.", "Invalid cell vector a", "Message: "//istat%string())
        call istat%setlevel(EVENT_LEVEL_QUIET)
        call istat%report
        call istat%discard

        call read_structure("bad_value_5.gin", box, positions, tags, format="gulp", except=istat)
        call check(istat=="Invalid element for cell vector b.", "Invalid cell vector b", "Message: "//istat%string())
        call istat%setlevel(EVENT_LEVEL_QUIET)
        call istat%report
        call istat%discard

        call read_structure("bad_value_6.gin", box, positions, tags, format="gulp", except=istat)
        call check(istat=="Invalid element for cell vector c.", "Invalid cell vector c", "Message: "//istat%string())
        call istat%setlevel(EVENT_LEVEL_QUIET)
        call istat%report
        call istat%discard

        call read_structure("bad_value_7.gin", box, positions, tags, format="gulp", except=istat)
        call check(istat=="Invalid atomic coordinate.", "Invalid atomic coordinate numerator", "Message: "//istat%string())
        call istat%setlevel(EVENT_LEVEL_QUIET)
        call istat%report
        call istat%discard

        call read_structure("bad_value_8.gin", box, positions, tags, format="gulp", except=istat)
        call check(istat=="Invalid atomic coordinate.", "Invalid atomic coordinate denominator", "Message: "//istat%string())
        call istat%setlevel(EVENT_LEVEL_QUIET)
        call istat%report
        call istat%discard

    end subroutine
!******************************************************************************!
!> Test error conditions in the specific implementation of write_structure.
!******************************************************************************!
    subroutine test_errors
        real(e), dimension(3,3) :: box
        real(e), dimension(:,:), allocatable :: positions, velocities
        character(TAG_LENGTH), dimension(:), allocatable :: tags
        type(FException) :: stat
!------

        call start_test_group("Error conditions")

        allocate(positions(3,0))
        allocate(tags(0))
        call write_structure(filename, box, positions, tags, format="gulp", stat=stat)
        call check(stat==0, "Empty structure", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call get_test_structure(box, positions, tags, velocities)

        call write_structure(filename, box, positions, tags(:4), format="gulp", stat=stat)
        call check(stat=="The tags array has the wrong size.", "size of the tags array", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard
    end subroutine
include 'structure_basic_checks.inc'
end program
