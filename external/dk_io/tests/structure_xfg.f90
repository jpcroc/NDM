program structure_xfg
    use iso_fortran_env, only: int64
    use dk_structure_io
    use dk_math, only: approx_equal
    use dk_fileobject, only: FileObject
    use dk_exception, only: FException, EVENT_LEVEL_QUIET
    use dk_test

    implicit none(external, type)

    integer, parameter :: e = kind(0.0d0)
    character(*), parameter :: filename = "structure.xfg"
    character(*), parameter :: format = "xfg"
    character(*), parameter :: format_check = "xfg"

    call test_basic_uncompressed
    call test_basic_gzip
    call test_basic_bzip2
    call test_basic_zstandard

    call test_errors
    call test_xfg_format
    call test_features

!     block
!         use iso_fortran_env, only: real64, int64
!
!         real(e), dimension(3,3) :: box
!         real(e), dimension(:,:), allocatable :: positions
!         character(TAG_LENGTH), dimension(:), allocatable :: tags
!         real(real64) :: count_rate
!         integer(int64) :: start_time, end_time
!
!         write(*,*) "Read test (large)"
!         call system_clock(start_time, count_rate)
!         call read_structure("large.xfg", box, positions, tags, format="xfg")
!         call system_clock(end_time)
!         write(*,*) "Timing:", (end_time-start_time)*1.0d0 / count_rate
!         write(*,*)
!     end block

    ! ! Zstd-compressed read test
    ! write(*,*) "Zstd-compresed read test (large)"
    ! call system_clock(start_time, count_rate)
    ! call read_structure("large.xfg.zst", box2, positions2, tags2)
    ! call system_clock(end_time)
    ! write(*,*) "Timing:", (end_time-start_time)*1.0d0 / count_rate
    ! write(*,*) " Box:", all(box==box2)
    ! write(*,*) " Positions alloc:", allocated(positions2)
    ! write(*,*) " Positions size:", size(positions)==size(positions2)
    ! write(*,*) " Positions content:", all(positions==positions2)
    ! write(*,*)

    ! ! Gzip-compressed read test
    ! write(*,*) "Gzip-compresed read test (large)"
    ! call system_clock(start_time, count_rate)
    ! call read_structure("large.xfg.gz", box2, positions2, tags2)
    ! call system_clock(end_time)
    ! write(*,*) "Timing:", (end_time-start_time)*1.0d0 / count_rate
    ! write(*,*) " Box:", all(box==box2)
    ! write(*,*) " Positions alloc:", allocated(positions2)
    ! write(*,*) " Positions size:", size(positions)==size(positions2)
    ! write(*,*) " Positions content:", all(positions==positions2)
    ! write(*,*)

!    ! Large file performance test
!    write(*,*) "Testing large compressed file"
!    call system_clock(start_time, count_rate)
!    call read_structure("imag_Cu_100_300k.7500000.cfg.gz", box2, positions2, tags2)
!    call system_clock(end_time)
!    write(*,*) "Timing:", (end_time-start_time)*1.0d0 / count_rate
!    write(*,*)

!     ! Read the data from the file
!     write(*,*) "Testing simple reading"
!     call system_clock(start_time, count_rate)
!     call file_read%init("large_uncompressed.cfg", "r")
!     i = -1
!     do while(stat == 0)
!         call file_read%read_line(line, trim=.false., stat=stat)
! !        write(*,'(i6,a)') len_trim(line), ":" // trim(line)
!         i = i + 1
!     end do
!     write(*,*) "lines read:", i, len_trim(line)
!     call system_clock(end_time)
!     call stat%discard
!     ! write(*,*) trim(line)
!     write(*,*) "Timing:", (end_time-start_time)*1.0d0 / count_rate
!     write(*,*)
! !    stop

    ! write(*,*) "Testing Fortran simple reading"
    ! call system_clock(start_time, count_rate)
    ! call file_read%init_fortran("large_uncompressed.cfg", "r")
    ! i = -1
    ! do while(stat == 0)
    !     call file_read%read_line(line, trim=.true., stat=stat)
    !     i = i + 1
    ! end do
    ! write(*,*) "lines read:", i, len_trim(line)
    ! call system_clock(end_time)
    ! ! write(*,*) trim(line)
    ! write(*,*) "Timing:", (end_time-start_time)*1.0d0 / count_rate
    ! write(*,*)


    ! ! Large file performance test
    ! write(*,*) "Testing large uncompressed file"
    ! call system_clock(start_time, count_rate)
    ! call read_structure("large_uncompressed.cfg", box2, positions2, tags2)
    ! call system_clock(end_time)
    ! write(*,*) "Timing:", (end_time-start_time)*1.0d0 / count_rate
    ! write(*,*)
    ! write(*,*) "Positions:", size(positions2, 1), size(positions2, 2)
    ! write(*,*) positions2(:,1)
    ! write(*,*) positions2(:,2)
    ! write(*,*) positions2(:,3)
    ! write(*,*) "..."
    ! write(*,*) positions2(:,size(positions2, 2)-2)
    ! write(*,*) positions2(:,size(positions2, 2)-1)
    ! write(*,*) positions2(:,size(positions2, 2))

    call tests_finished
contains
!******************************************************************************!
!> Test the detection of varions syntax error in structure files.
!******************************************************************************!
    subroutine test_xfg_format
        real(e), dimension(3,3) :: box
        real(e), dimension(:,:), allocatable :: positions
        character(TAG_LENGTH), dimension(:), allocatable :: tags
        type(FException) :: stat


        call start_test_group("XFG format")

        call read_structure("empty.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat==0, "empty configuration", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_1.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Missing entry_count line.", "missing number of properties", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_2.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="The first non-comment line does not have the number of particles.", &
            "missing number of particles", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_3.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Invalid value for A.", "invalid A", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_4.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Invalid auxiliary index.", "invalid auxiliary index", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_5.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Auxiliary index out of bounds.", "auxiliary index out of bounds", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_6.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Invalid value for the number of particles.", "invalid number of particles", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_7.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Invalid value for entry_count.", "invalid entry_count value", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_8.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="At least 6 fields are necessary for atomic positions and velocities.", "missing fields 1", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_9.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="At least 3 fields are necessary for atomic positions.", "missing fields 2", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_10.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Invalid keyword.", "invalid keyword", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_11.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Invalid value for H0(1,1).", "invalid value for H0(1,1)", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_12.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Invalid value for H0(1,2).", "invalid value for H0(1,2)", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_13.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Invalid value for H0(1,3).", "invalid value for H0(1,3)", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_14.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Invalid value for H0(2,1).", "invalid value for H0(2,1)", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_15.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Invalid value for H0(2,2).", "invalid value for H0(2,2", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_16.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Invalid value for H0(2,3).", "invalid value for H0(2,3)", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_17.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Invalid value for H0(3,1).", "invalid value for H0(3,1)", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_18.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Invalid value for H0(3,2).", "invalid value for H0(3,2)", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call read_structure("bad_19.xfg", box, positions, tags, format="xfg", except=stat)
        call check(stat=="Invalid value for H0(3,3).", "invalid value for H0(3,3)", "Message: "//stat%string())
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

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

        call write_structure("structure.xfg", box, positions, tags(:4), stat=stat)
        call check(stat=="The tags array has the wrong size.", "size of the tags array")
        call stat%discard

        call write_structure("structure.xfg", box, positions, tags, velocities=velocities(:,:4), stat=stat)
        call check(stat=="The velocities array has the wrong size.", "size of the velocities array")
        call stat%discard

        call write_structure("structure.xfg", box, positions, tags, masses=masses(:4), stat=stat)
        call check(stat=="The masses array has the wrong size.", "size of the masses array")
        call stat%discard

        call write_structure("structure.xfg", box, positions, tags, forces=velocities(:,:4), stat=stat)
        call check(stat=="The forces array has the wrong size.", "size of the forces array")
        call stat%discard
    end subroutine

    subroutine test_features
        real(e), dimension(3,3) :: box!, box2
        real(e), dimension(:,:), allocatable :: positions, velocities!, positions2
        real(e), dimension(:), allocatable :: masses
        character(TAG_LENGTH), dimension(:), allocatable :: tags!, tags2
!         character(:), allocatable :: format
!         character(:), allocatable :: compression
!         integer :: i
!         integer(int64) :: start_time, end_time, count_rate
!         type(FileObject) :: file_read
!         character(:), allocatable :: line
        type(FException) :: stat

        call get_test_structure(box, positions, tags)

        allocate(masses(size(tags)))
        masses = 42
        allocate(velocities(3,size(tags)))
        velocities = 0.42

        call start_test_group("Features")

        call write_structure("structure.xfg", box, positions, tags, velocities=velocities, stat=stat)
        call check(stat==0, "saving velocities")
        call stat%discard

        call write_structure("structure.xfg", box, positions, tags, masses=masses, stat=stat)
        call check(stat==0, "saving masses")
        call stat%discard

        call write_structure("structure.xfg", box, positions, tags, forces=velocities, stat=stat)
        call check(stat==0, "saving forces")
        call stat%discard
    end subroutine

include 'structure_basic_checks.inc'
end program

