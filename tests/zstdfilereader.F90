program test
    use iso_fortran_env, only: int64
    use dk_fileobject, only: FileObject
    use dk_exception, only: FException
    use dk_test

    implicit none(external, type)

    integer :: i
    character(1000) :: buffer

    do i=1, command_argument_count()
        call get_command_argument(i, buffer)
        if (buffer=="-p" .or. buffer=="--test-performances") call test_performances
    end do

    call test_append

    call test_errors

    call test_bad_files

    call tests_finished
contains
    subroutine test_performances
        integer, parameter :: e = kind(0.0d0)
        type(FileObject) :: file_write, file_read
        character(:), allocatable :: line
        character(100), dimension(:), allocatable :: data
        integer :: i, j, n, lines_count
        real(e) :: x, data_size
        integer, dimension(:), allocatable :: seed
        integer(int64) :: start_time, end_time, count_rate
        type(FException) :: stat
!------
        ! Set the data array
        lines_count = 250000
        allocate(data(lines_count))
        call random_seed(size=n)
        allocate(seed(n))
        seed(:) = 42
        call random_seed(put=seed)
        allocate(character(100) :: line)
        do i=1, lines_count
            do j=1, 100
                call random_number(x)
                n = nint(x)*25+1
                line(j:j) = achar(nint(x*90)+33)
            end do
            data(i) = line
        end do
        data_size = size(data)*storage_size(data)/8.0d6
        write(*,*) "File size:", data_size, "MB"
        write(*,*)

        ! Write the data to a file
        write(*,*) "Testing zstd-compressed writing"
        call system_clock(start_time, count_rate)
        call file_write%init("test_zstd.txt.zst", "w")
        do i=1, lines_count
            call file_write%write_line(data(i))
        end do
        call file_write%close
        call system_clock(end_time)
        write(*,*) "Timing:", (end_time-start_time)*1.0d0 / count_rate
        write(*,*) "Writing speed:", data_size / ((end_time-start_time)*1.0d0 / count_rate), "MB/s"
        write(*,*)

        ! Read the data from the file
        write(*,*) "Testing zstd-compressed reading"
        call system_clock(start_time, count_rate)
        call file_read%init("test_zstd.txt.zst", "r")
        do i=1, lines_count
            call file_read%read_line(line)
            if (line /= data(i)) then
            !LCOV_EXCL_START
                write(*,*) "Error on line ", i
                write(*,'(a)') line
                write(*,'(a)') data(i)
                error stop
            !LCOV_EXCL_STOP
            end if
        end do
        call system_clock(end_time)
        write(*,*) "Timing:", (end_time-start_time)*1.0d0 / count_rate
        write(*,*) "Reading speed:", data_size / ((end_time-start_time)*1.0d0 / count_rate), "MB/s"
        write(*,*)

        call start_test_group("file rewinding")
        call file_read%rewind(stat=stat)
        call check(stat==0, "rewinding file in read mode", stat%string())
        call stat%discard
        call file_read%read_line(line, stat=stat)
        call check(stat==0, "reading from file beginning", stat%string())
        call stat%discard
        call check(line==data(1), "checking file content", line)
        call file_read%close

        call file_write%init("rewind.txt.zst", "w")
        call file_write%write_line("LINE 1")
        call file_write%write_line("LINE 2")
        call file_write%rewind(stat=stat)
        call check(stat==0, "rewinding file in write mode", stat%string())
        call stat%discard
        call file_write%write_line("LINE 3", stat=stat)
        call check(stat==0, "writing to file beginning", stat%string())
        call stat%discard
        call file_write%close
        call file_read%init("rewind.txt.zst", "r")
        call file_read%read_line(line)
        call check(line=="LINE 3", "checking file content", line)
    end subroutine
    
    subroutine test_append
        type(FileObject) :: file
        type(FException) :: stat
        character(:), allocatable :: line

        call start_test_group("APPEND mode")
        call file%init("append.txt.zst", "w", stat=stat)
        call check(stat==0, "opening file in WRITE mode", stat%string())
        call stat%discard
        call file%write_line("qwertyuiop", stat=stat)
        call check(stat==0, "writing line 1", stat%string())
        call stat%discard
        call file%close(stat=stat)
        call check(stat==0, "closing file", stat%string())
        call stat%discard

        call file%init("append.txt.zst", "a", stat=stat)
        call check(stat==0, "opening file in APPEND mode", stat%string())
        call stat%discard
        call file%write_line("1234567890", stat=stat)
        call check(stat==0, "writing line 2", stat%string())
        call stat%discard
        call file%close(stat=stat)
        call check(stat==0, "closing file", stat%string())
        call stat%discard

        call file%init("append.txt.zst", "r", stat=stat)
        call check(stat==0, "opening file in READ mode", stat%string())
        call stat%discard
        call file%read_line(line, stat=stat)
        call check(stat==0, "reading line 1", stat%string())
        call stat%discard
        call check(line=="qwertyuiop", "checking line 1", line)
        call file%read_line(line, stat=stat)
        call check(stat==0, "reading line 2", stat%string())
        call stat%discard
        call check(line=="1234567890", "checking line 2", line)
    end subroutine

    subroutine test_errors
#ifdef HAVE_LIBZSTD
        use dk_zstdfilereader, only: ZstdFileReader

        type(ZstdFileReader) :: reader
        type(FException) :: stat
        character(:), allocatable :: line
        integer :: n
!------
        ! Make sure the test file is not already there
        call execute_command_line("rm -f newfile.txt")
        allocate(character(100) :: line)

        call start_test_group("error conditions")

        call reader%open_file("newfile.txt", "invalid", stat)
        call check(stat=="Invalid mode", "opening file with invalid mode", stat%string())
        call stat%discard

        call reader%open_file("newfile.txt", "read", stat)
        call check(stat=="No such file or directory", "opening non-existing file for reading", stat%string())
        call stat%discard

        call reader%open_file("invalid/newfile.txt", "write", stat)
        call check(stat=="No such file or directory", "opening non-existing file for writing", stat%string())
        call stat%discard

        call reader%open_file("invalid/newfile.txt", "append", stat)
        call check(stat=="No such file or directory", "opening non-existing file for appending", stat%string())
        call stat%discard

        call reader%read(line, n, stat)
        call check(stat=="File not open.", "reading from un-opened file", stat%string())
        call stat%discard

        call reader%write(line, n, stat)
        call check(stat=="File not open.", "writing to un-opened file", stat%string())
        call stat%discard

        call reader%open_file("newfile.txt", "write", stat)
        call check(stat==0, "opening file for writing", stat%string())
        call stat%discard

        call reader%open_file("newfile.txt", "write", stat)
        call check(stat=="File reader already in use", "opening a file twice", stat%string())
        call stat%discard

        call reader%read(line, n, stat)
        call check(stat=="File not open for reading", "reading from a file opened for writing", stat%string())
        call stat%discard

        call reader%close_file
        call reader%open_file("newfile.txt", "read", stat)
        call reader%write(line, n, stat)
        call check(stat=="File not open for writing", "writing to a file opened for reading", stat%string())
        call stat%discard

        call reader%close_file
        call reader%open_file("newfile.txt", "write", stat)
        call check(stat==0, "opening file for writing", stat%string())
        call stat%discard
        call execute_command_line("chmod -w newfile.txt")
        call reader%write("qwertyuiop", 0, stat)
        call check(stat==0, "writing to a read-only file", stat%string())
#endif
    end subroutine

    subroutine test_bad_files
#ifdef HAVE_LIBZSTD
        use dk_zstdfilereader, only: ZstdFileReader

        type(ZstdFileReader), allocatable :: reader
        type(FException) :: stat
        character(:), allocatable :: line
        integer :: n
!------
        allocate(character(100) :: line)

        call start_test_group("malformed files")

        allocate(reader)
        call reader%open_file("magic.zst", "r", stat)
        call reader%read(line, n, stat)
        deallocate(reader)
        call check(stat=="Unknown frame descriptor", "wrong magic number", stat%string())
        call stat%discard

        allocate(reader)
        call reader%open_file("corrupted.zst", "r", stat)
        call reader%read(line, n, stat)
        deallocate(reader)
        call check(stat=="Restored data doesn't match checksum", "corrupted file", stat%string())
        call stat%discard

        allocate(reader)
        call reader%open_file("incomplete.zst", "r", stat)
        call reader%read(line, n, stat)
        deallocate(reader)
        call check(stat=="End of file.", "incomplete file", stat%string())
        call stat%discard
#endif
    end subroutine
end program
