program structure_xsf
    use iso_fortran_env, only: int64
    use dk_structure_io
    use dk_math, only: approx_equal
    use dk_fileobject, only: FileObject
    use dk_exception, only: FException, EVENT_LEVEL_QUIET
    use dk_test

    implicit none(external, type)

    integer, parameter :: e = kind(0.0d0)
    character(*), parameter :: filename = "structure.xsf"
    character(*), parameter :: format = "xsf"
    character(*), parameter :: format_check = "xsf"

    call test_errors
    call test_features

    call test_basic_uncompressed
    call test_basic_gzip
    call test_basic_bzip2
    call test_basic_zstandard

    call tests_finished
contains
!******************************************************************************!
!> Test error conditions in the specific implementation of write_structure.
!******************************************************************************!
    subroutine test_errors
        real(e), dimension(3,3) :: box
        real(e), dimension(:,:), allocatable :: positions
        character(TAG_LENGTH), dimension(:), allocatable :: tags
        type(FException) :: stat
!------
        call get_test_structure(box, positions, tags)

        call start_test_group("Error conditions")

        call write_structure("structure.xsf", box, positions, tags(:4), stat=stat)
        call check(stat=="The tags array has the wrong size.", "size of the tags array")
        call stat%discard

        call read_structure("bad_1.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Invalid atomic coordinate.", "comment in atoms block")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_2.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Missing atoms.", "missing atom")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_3.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Invalid atomic coordinate.", "invalid atom coordinate")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_4.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Invalid element for cell vector c.", "missing cell vector")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_5.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Invalid element for cell vector c.", "invalid primitive cell vector c")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_6.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Invalid element for cell vector a.", "invalid primitive cell vector a")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_7.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Invalid element for cell vector b.", "invalid primitive cell vector b")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_8.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Invalid keyword.", "invalid structure type")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_9.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Only CRYSTAL files are supported.", "unsupported structure type")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_10.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Invalid keyword.", "invalid keyword")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_11.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Invalid element for cell vector a.", "invalid conventional cell vector a")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_12.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Invalid element for cell vector b.", "invalid conventional cell vector b")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_13.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Junk after cell vector c.", "invalid conventional cell vector c")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_14.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Syntax error.", "invalid conventional cell vector c")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_15.xsf", box, positions, tags, format="xsf", except=stat)
        call check(stat=="Invalid number of atoms.", "invalid number of atoms")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard
    end subroutine

    subroutine test_features
!         real(e), dimension(3,3) :: box, box2
!         real(e), dimension(:,:), allocatable :: positions, velocities, positions2, velocities2, forces, forces2
!         real(e), dimension(:), allocatable :: masses
!         character(TAG_LENGTH), dimension(:), allocatable :: tags, tags2
!         type(FException) :: stat
!
!         call get_test_structure(box, positions, tags, forces=forces, velocities=velocities, masses=masses)
!
!         call start_test_group("Features")
!
!         call write_structure("structure.xyz", box, positions, tags, velocities=velocities, stat=stat)
!         call check(stat==0, "Saving velocities")
!         call stat%discard
!         call read_structure("structure.xyz", box2, positions2, tags2, velocities=velocities2, format="xyz", except=stat)
!         call check(stat==0, "Reading velocities", "Message: "//stat%string())
!         call check(all(approx_equal(box,box2)), "Simulation box")
!         call check(allocated(positions2), "Positions array allocation status")
!         if (allocated(positions2)) then
!             call check(size(positions)==size(positions2), "Positions array size")
!             call check(all(approx_equal(positions, positions2)), "Positions array values")
!         end if
!         call check(allocated(velocities2), "Velocities array allocation status")
!         if (allocated(velocities2)) then
!             call check(size(velocities)==size(velocities2), "Velocities array size")
!             call check(all(approx_equal(velocities, velocities2)), "Velocities array values")
!         end if
!
!         ! call write_structure("structure.xyz", box, positions, tags, masses=masses, stat=stat)
!         ! write(*,*) "Saved masses"
!         ! call stat%report
!         ! call check(stat==0, "saving masses")
!         ! call stat%discard
!
!         call write_structure("structure.xyz", box, positions, tags, forces=forces, stat=stat)
!         call check(stat==0, "Saving forces")
!         call stat%discard
!         call read_structure("structure.xyz", box2, positions2, tags2, forces=forces2, format="xyz", except=stat)
!         call check(stat==0, "Reading forces", "Message: "//stat%string())
!         call check(all(approx_equal(box,box2)), "Simulation box")
!         call check(allocated(positions2), "Positions array allocation status")
!         if (allocated(positions2)) then
!             call check(size(positions)==size(positions2), "Positions array size")
!             call check(all(approx_equal(positions, positions2)), "Positions array values")
!         end if
!         call check(allocated(forces2), "Forces array allocation status")
!         if (allocated(forces2)) then
!             call check(size(forces)==size(forces2), "Forces array size")
!             call check(all(approx_equal(forces, forces2)), "Forces array values")
!         end if
!
!         call read_structure("valid_1.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat==0, "Cell matrix comment", "Message: "//stat%string())
!         call check(all(approx_equal(box,box2)), "simulation box")
!         call check(allocated(positions), "positions array allocation status")
!         if (allocated(positions)) then
!             call check(size(positions)==size(positions2), "positions array size")
!             call check(all(approx_equal(positions, positions2)), "positions array values")
!         end if
!
!         call read_structure("valid_2.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat==0, "Cell parameters comment", "Message: "//stat%string())
!         call check(all(approx_equal(box,box2)), "simulation box")
!         call check(allocated(positions), "positions array allocation status")
!         if (allocated(positions)) then
!             call check(size(positions)==size(positions2), "positions array size")
!             call check(all(approx_equal(positions, positions2)), "positions array values")
!         end if
!
!         call read_structure("valid_3.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat==0, "Spaces around = sign", "Message: "//stat%string())
!         call check(all(approx_equal(box,box2)), "simulation box")
!         call check(allocated(positions), "positions array allocation status")
!         if (allocated(positions)) then
!             call check(size(positions)==size(positions2), "positions array size")
!             call check(all(approx_equal(positions, positions2)), "positions array values")
!         end if
!
!         call read_structure("valid_4.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat==0, "Bareword value in the middle of the property line", "Message: "//stat%string())
!         call check(all(approx_equal(box,box2)), "simulation box")
!         call check(allocated(positions), "positions array allocation status")
!         if (allocated(positions)) then
!             call check(size(positions)==size(positions2), "positions array size")
!             call check(all(approx_equal(positions, positions2)), "positions array values")
!         end if
!
!         call read_structure("bad_1.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat=="Unexpected end of file.", "Truncated file", "Message: "//stat%string())
!         call stat%setlevel(EVENT_LEVEL_QUIET)
!         call stat%report
!         call stat%discard
!
!         call read_structure("bad_2.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat=="Invalid number value.", "Wrong atomic position", "Message: "//stat%string())
!         call stat%setlevel(EVENT_LEVEL_QUIET)
!         call stat%report
!         call stat%discard
!
!         call read_structure("bad_3.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat=="Missing data on line.", "Truncated line", "Message: "//stat%string())
!         call stat%setlevel(EVENT_LEVEL_QUIET)
!         call stat%report
!         call stat%discard
!
!         call read_structure("bad_4.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat=="Extended XYZ files are supported only if the first column contains the species names.", &
!             "Wrong auxiliary field order", "Message: "//stat%string())
!         call stat%setlevel(EVENT_LEVEL_QUIET)
!         call stat%report
!         call stat%discard
!
!         call read_structure("bad_5.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat=="Expected '=' sign after property name.", "Missing = character", "Message: "//stat%string())
!         call stat%setlevel(EVENT_LEVEL_QUIET)
!         call stat%report
!         call stat%discard
!
!         call read_structure("bad_6.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat=="Unterminated double-quoted string.", "Unterminated quoted value", "Message: "//stat%string())
!         call stat%setlevel(EVENT_LEVEL_QUIET)
!         call stat%report
!         call stat%discard
!
!         call read_structure("bad_7.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat=="Expected value after '=' sign.", "Missing property value", "Message: "//stat%string())
!         call stat%setlevel(EVENT_LEVEL_QUIET)
!         call stat%report
!         call stat%discard
!
!         call read_structure("bad_8.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat=="Expected '=' sign after property name.", "Whitespace in property name", "Message: "//stat%string())
!         call stat%setlevel(EVENT_LEVEL_QUIET)
!         call stat%report
!         call stat%discard
!
!         call read_structure("bad_9.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat=="Invalid cell shape.", "Invalid cell box", "Message: "//stat%string())
!         call stat%setlevel(EVENT_LEVEL_QUIET)
!         call stat%report
!         call stat%discard
!
!         call read_structure("bad_10.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat=="Only real auxiliary fields are supported in extended XYZ files.", "Integer auxiliary value", "Message: "//stat%string())
!         call stat%setlevel(EVENT_LEVEL_QUIET)
!         call stat%report
!         call stat%discard
!
!         call read_structure("bad_11.xyz", box2, positions2, tags2, format="xyz", except=stat)
!         call check(stat=="Only scalar auxiliary fields are supported in extended XYZ files.", "Vector auxiliary value", "Message: "//stat%string())
!         call stat%setlevel(EVENT_LEVEL_QUIET)
!         call stat%report
!         call stat%discard
    end subroutine

include 'structure_basic_checks.inc'
end program

