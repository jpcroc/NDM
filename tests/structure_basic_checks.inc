!******************************************************************************!
!> Generic i/o tests with bzip2-compressed files.
!******************************************************************************!
    subroutine test_basic_bzip2
        real(e), dimension(3,3) :: box, box2
        real(e), dimension(:,:), allocatable :: positions, positions2
        character(TAG_LENGTH), dimension(:), allocatable :: tags, tags2
        type(FException) :: stat
        character(:), allocatable :: iformat
!------
        call get_test_structure(box, positions, tags)

        ! Writing test
        call start_test_group("Basic i/o, bzip2")
        call write_structure(filename//".bz2", box, positions, tags, format=format, stat=stat)
        call check(stat==0, "writing structure file")

        call get_file_format(filename, iformat, stat=stat)
        call check(stat == 0, "getting file format", "Message: "//stat%string())
        if (stat == 0) then
            call check(iformat==format_check, "file format", "Format: '"//iformat//"' instead of '"//format_check//"'")
        end if
        call stat%discard

        ! Reading test
        call read_structure(filename//".bz2", box2, positions2, tags2, format=format)
        call check(all(approx_equal(box, box2)), "simulation box")
        if (.not.all(approx_equal(box,box2))) then
            write(*,*) box(1,:)
            write(*,*) box(2,:)
            write(*,*) box(3,:)
            write(*,*)
            write(*,*) box2(1,:)
            write(*,*) box2(2,:)
            write(*,*) box2(3,:)
        end if
        call check(allocated(positions2), "positions array allocation status")
        if (allocated(positions2)) then
            call check(size(positions)==size(positions2), "positions array size")
            if (size(positions)==size(positions2)) &
                call check(all(approx_equal(positions, positions2)), "positions array values")
        end if
    end subroutine
!******************************************************************************!
!> Generic i/o tests with gzip-compressed files.
!******************************************************************************!
    subroutine test_basic_gzip
        integer, parameter :: e = kind(0.0d0)
        real(e), dimension(3,3) :: box, box2
        real(e), dimension(:,:), allocatable :: positions, positions2
        character(TAG_LENGTH), dimension(:), allocatable :: tags, tags2
        type(FException) :: stat
        character(:), allocatable :: iformat
!------
        call get_test_structure(box, positions, tags)

        ! Writing test
        call start_test_group("Basic i/o, gzip")
        call write_structure(filename//".gz", box, positions, tags, format=format, stat=stat)
        call check(stat==0, "writing structure file")

        call get_file_format(filename, iformat, stat=stat)
        call check(stat == 0, "getting file format", "Message: "//stat%string())
        if (stat == 0) then
            call check(iformat==format_check, "file format", "Format: '"//iformat//"' instead of '"//format_check//"'")
        end if
        call stat%discard

        ! Reading test
        call read_structure(filename//".gz", box2, positions2, tags2, format=format)
        call check(all(approx_equal(box, box2)), "simulation box")
        if (.not.all(approx_equal(box,box2))) then
            write(*,*) box(1,:)
            write(*,*) box(2,:)
            write(*,*) box(3,:)
            write(*,*)
            write(*,*) box2(1,:)
            write(*,*) box2(2,:)
            write(*,*) box2(3,:)
        end if
        call check(allocated(positions2), "positions array allocation status")
        if (allocated(positions2)) then
            call check(size(positions)==size(positions2), "positions array size")
            if (size(positions)==size(positions2)) &
                call check(all(approx_equal(positions, positions2)), "positions array values")
        end if
    end subroutine
!******************************************************************************!
!> Generic i/o tests with uncompressed files.
!******************************************************************************!
    subroutine test_basic_uncompressed
        integer, parameter :: e = kind(0.0d0)
        real(e), dimension(3,3) :: box, box2
        real(e), dimension(:,:), allocatable :: positions, positions2
        character(TAG_LENGTH), dimension(:), allocatable :: tags, tags2
        type(FException) :: stat
        character(:), allocatable :: iformat
!------
        call get_test_structure(box, positions, tags)

        ! Writing test
        call start_test_group("Basic i/o uncompressed")

        call write_structure("wrong/"//filename, box, positions, tags, format=format, stat=stat)
        call check(stat=="No such file or directory.", "wrong file name", "Message: "//stat%string())
        call stat%discard

        call write_structure(filename, box, positions, tags, format=format, stat=stat)
        call check(stat==0, "writing structure file")

        call get_file_format(filename, iformat, stat=stat)
        call check(stat == 0, "getting file format", "Message: "//stat%string())
        if (stat == 0) then
            call check(iformat==format_check, "file format", "Format: '"//iformat//"' instead of '"//format_check//"'")
        end if
        call stat%discard

        ! Reading tests
        call read_structure(filename, box2, positions2, tags2, format=format)
        call check(all(approx_equal(box,box2)), "simulation box")
        if (.not.all(approx_equal(box,box2))) then
            write(*,*) box(1,:)
            write(*,*) box(2,:)
            write(*,*) box(3,:)
            write(*,*)
            write(*,*) box2(1,:)
            write(*,*) box2(2,:)
            write(*,*) box2(3,:)
        end if
        call check(allocated(positions2), "positions array allocation status")
        if (allocated(positions2)) then
            call check(size(positions)==size(positions2), "positions array size")
            if (size(positions)==size(positions2)) &
                call check(all(approx_equal(positions, positions2)), "positions array values")
        end if
    end subroutine
!******************************************************************************!
!> Generic i/o tests with zstandard-compressed files.
!******************************************************************************!
    subroutine test_basic_zstandard
        integer, parameter :: e = kind(0.0d0)
        real(e), dimension(3,3) :: box, box2
        real(e), dimension(:,:), allocatable :: positions, positions2
        character(TAG_LENGTH), dimension(:), allocatable :: tags, tags2
        type(FException) :: stat
        character(:), allocatable :: iformat
!------
        call get_test_structure(box, positions, tags)

        ! Writing test
        call start_test_group("Basic i/o, zstandard")
        call write_structure(filename//".zst", box, positions, tags, format=format, stat=stat)
        call check(stat==0, "writing structure file")

        call get_file_format(filename, iformat, stat=stat)
        call check(stat == 0, "getting file format", "Message: "//stat%string())
        if (stat == 0) then
            call check(iformat==format_check, "file format", "Format: '"//iformat//"' instead of '"//format_check//"'")
        end if
        call stat%discard

        ! Reading test
        call read_structure(filename//".zst", box2, positions2, tags2, format=format)
        call check(all(approx_equal(box, box2)), "simulation box")
        if (.not.all(approx_equal(box,box2))) then
            write(*,*) box(1,:)
            write(*,*) box(2,:)
            write(*,*) box(3,:)
            write(*,*)
            write(*,*) box2(1,:)
            write(*,*) box2(2,:)
            write(*,*) box2(3,:)
        end if
        call check(allocated(positions2), "positions array allocation status")
        if (allocated(positions2)) then
            call check(size(positions)==size(positions2), "positions array size")
            if (size(positions)==size(positions2)) &
                call check(all(approx_equal(positions, positions2)), "positions array values")
        end if
    end subroutine
