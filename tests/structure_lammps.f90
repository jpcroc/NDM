program structure_lammps
    use iso_fortran_env, only: int64
    use dk_structure_io
    use dk_math, only: approx_equal
    use dk_fileobject, only: FileObject
    use dk_exception, only: FException
    use dk_test

    implicit none(external, type)

    integer, parameter :: e = kind(0.0d0)
    character(*), parameter :: filename = "structure.lconfig"
    character(*), parameter :: format = "lammps"
    character(*), parameter :: format_check = "lammps"

    call test_basic_uncompressed
    call test_basic_gzip
    call test_basic_bzip2
    call test_basic_zstandard

!     call test_errors
!     call test_features


!     write(*,*) "Read test (large)"
!     call system_clock(start_time, count_rate)
!     call read_structure("large.lconfig", box, positions, tags)
!     call system_clock(end_time)
!     write(*,*) "Timing:", (end_time-start_time)*1.0d0 / count_rate
!     write(*,*)

    ! ! Zstd-compressed read test
    ! write(*,*) "Zstd-compresed read test (large)"
    ! call system_clock(start_time, count_rate)
    ! call read_structure("large.lconfig.zst", box2, positions2, tags2)
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
    ! call read_structure("large.lconfig.gz", box2, positions2, tags2)
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
include 'structure_basic_checks.inc'
end program

