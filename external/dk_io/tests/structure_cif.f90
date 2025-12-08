program structure_cif
    use iso_fortran_env, only: int64
    use dk_structure_io
    use dk_math, only: approx_equal
    use dk_fileobject, only: FileObject
    use dk_exception, only: FException
    use dk_test

    implicit none(external, type)

    integer, parameter :: e = kind(0.0d0)
    character(*), parameter :: filename = "structure.cif"
    character(*), parameter :: format = "cif"
    character(*), parameter :: format_check = "cif"

    call test_triclinic

    call test_basic_uncompressed
    call test_basic_gzip
    call test_basic_bzip2
    call test_basic_zstandard

    call test_errors
!     call test_features
    call test_cif_format

    call tests_finished
contains
!******************************************************************************!
!> Test the detection of varions syntax error in structure files.
!******************************************************************************!
    subroutine test_cif_format
        real(e), dimension(3,3) :: box
        real(e), dimension(:,:), allocatable :: positions
        character(TAG_LENGTH), dimension(:), allocatable :: tags
        integer :: stat
        character(:), allocatable :: msg

        call start_test_group("CIF format")

        call read_structure("symmetry.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(stat==0, "Symmetrised fluorite structure", "Message: "//msg)
        call write_structure("symmetry.xfg", box, positions, tags)

        call read_structure("symmetry2.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(stat==0, "Symmetrised perovskite structure", "Message: "//msg)
        call write_structure("symmetry2.xfg", box, positions, tags)

        call read_structure("bad_1.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(msg=="Invalid property name.", "Invalid property name", "Message: "//msg)

        call read_structure("bad_2.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(msg=="Expected data_ block.", "Missing data_ block", "Message: "//msg)

        call read_structure("bad_3.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(msg=="Unterminated string.", "Invalid property value", "Message: "//msg)
        ! call read_structure("bad_3.cif", box, positions, tags, format="cif")

        call read_structure("bad_4.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(msg=="Missing atomic tags.", "Missing atomic tags", "Message: "//msg)

        call read_structure("bad_5.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(msg=="Unterminated string.", "Unterminated quoted value", "Message: "//msg)

        call read_structure("bad_6.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(msg=="Unterminated string.", "Unterminated double-quoted value", "Message: "//msg)

        call read_structure("bad_7.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(msg=="Unterminated multi-line value.", "Unterminated multi-line value", "Message: "//msg)

        call read_structure("bad_8.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(msg=="Unknown space group.", "Unknown space group", "Message: "//msg)

        call read_structure("bad_value_1.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(msg=="Invalid value for _cell_length_a.", "Invalid lattice parameter a", "Message: "//msg)

        call read_structure("bad_value_2.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(msg=="Invalid value for _cell_length_b.", "Invalid lattice parameter b", "Message: "//msg)

        call read_structure("bad_value_3.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(msg=="Invalid value for _cell_length_c.", "Invalid lattice parameter c", "Message: "//msg)

        call read_structure("bad_value_4.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(msg=="Invalid value for _cell_angle_alpha.", "Invalid lattice angle alpha", "Message: "//msg)

        call read_structure("bad_value_5.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(msg=="Invalid value for _cell_angle_beta.", "Invalid lattice angle beta", "Message: "//msg)

        call read_structure("bad_value_6.cif", box, positions, tags, format="cif", stat=stat, msg=msg)
        call check(msg=="Invalid value for _cell_angle_gamma.", "Invalid lattice angle gamma", "Message: "//msg)

    end subroutine
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

        call write_structure("structure.cif", box, positions, tags(:4), format="cif", stat=stat)
        call check(stat=="The tags array has the wrong size.", "size of the tags array", "Message: "//stat%string())
        call stat%discard
    end subroutine

    subroutine test_triclinic
        real(e), dimension(3,3) :: box, box2
        real(e), dimension(:,:), allocatable :: positions, positions2
        type(FException) :: stat
        integer :: istat
        character(:), allocatable :: msg
        character(TAG_LENGTH), dimension(:), allocatable :: tags, tags2

        call start_test_group("Triclinic boxes")

        call read_structure("triclinic.cif", box, positions, tags, format="cif", stat=istat, msg=msg)

        call get_triclinic_test_structure(box2, positions2, tags2)
        call write_structure("triclinic_test.cif", box2, positions2, tags2, format="cif", stat=stat)
    end subroutine

    ! subroutine test_features
    !     real(e), dimension(3,3) :: box, box2
    !     real(e), dimension(:,:), allocatable :: positions, positions2, velocities
    !     real(e), dimension(:), allocatable :: masses
    !     character(TAG_LENGTH), dimension(:), allocatable :: tags, tags2
    !     character(:), allocatable :: format
    !     character(:), allocatable :: compression
    !     integer :: i
    !     integer(int64) :: start_time, end_time, count_rate
    !     type(FileObject) :: file_read
    !     character(:), allocatable :: line
    !     type(FException) :: stat

    !     call get_structure(box, positions, tags)

    !     allocate(masses(size(tags)))
    !     masses = 42
    !     allocate(velocities(3,size(tags)))
    !     velocities = 0.42

    !     call start_test_group("Features")

    !     call write_structure("structure.xfg", box, positions, tags, velocities=velocities, stat=stat)
    !     call check(stat==0, "saving velocities")
    !     call stat%discard

    !     call write_structure("structure.xfg", box, positions, tags, masses=masses, stat=stat)
    !     call check(stat==0, "saving masses")
    !     call stat%discard

    !     call write_structure("structure.xfg", box, positions, tags, forces=velocities, stat=stat)
    !     call check(stat==0, "saving forces")
    !     call stat%discard
    ! end subroutine
include 'structure_basic_checks.inc'
end program
