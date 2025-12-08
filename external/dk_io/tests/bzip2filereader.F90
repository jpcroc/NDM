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
        write(*,*) "Testing bzip2-compressed writing"
        call system_clock(start_time, count_rate)
        call file_write%init("test_bzip2.txt.bz2", "w")
        do i=1, lines_count
            call file_write%write_line(data(i))
        end do
        call file_write%close
        call system_clock(end_time)
        write(*,*) "Timing:", (end_time-start_time)*1.0d0 / count_rate
        write(*,*) "Writing speed:", data_size / ((end_time-start_time)*1.0d0 / count_rate), "MB/s"
        write(*,*)

        ! Read the data from the file
        write(*,*) "Testing bzip2-compressed reading"
        call system_clock(start_time, count_rate)
        call file_read%init("test_bzip2.txt.bz2", "r")
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
        call check(stat==0, "rewinding file in read mode")
        call file_read%read_line(line, stat=stat)
        call check(stat==0, "reading from file beginning")
        call check(line==data(1), "checking file content")
        call file_read%close

        call file_write%init("rewind.txt.bz2", "w")
        call file_write%write_line("LINE 1")
        call file_write%write_line("LINE 2")
        call file_write%rewind(stat=stat)
        call check(stat==0, "rewinding file in write mode")
        call file_write%write_line("LINE 3", stat=stat)
        call check(stat==0, "writing to file beginning")
        call file_write%close
        call file_read%init("rewind.txt.bz2", "r")
        call file_read%read_line(line)
        call check(line=="LINE 3", "checking file content", line)
    end subroutine
    
    subroutine test_append
        type(FileObject) :: file
        type(FException) :: stat
        character(:), allocatable :: line

        call start_test_group("APPEND mode")
        call file%init("append.txt.bz2", "w", stat=stat)
        call check(stat==0, "opening file in WRITE mode", stat%string())
        call stat%discard
        call file%write_line("qwertyuiop", stat=stat)
        call check(stat==0, "writing line 1", stat%string())
        call stat%discard
        call file%close(stat=stat)
        call check(stat==0, "closing file", stat%string())
        call stat%discard

        call file%init("append.txt.bz2", "a", stat=stat)
        call check(stat==0, "opening file in APPEND mode", stat%string())
        call stat%discard
        call file%write_line("1234567890", stat=stat)
        call check(stat==0, "writing line 2", stat%string())
        call stat%discard
        call file%close(stat=stat)
        call check(stat==0, "closing file", stat%string())
        call stat%report
        call stat%discard

        call file%init("append.txt.bz2", "r", stat=stat)
        call check(stat==0, "opening file in READ mode", stat%string())
        call file%read_line(line, stat=stat)
        call check(stat==0, "reading line 1", stat%string())
        call check(line=="qwertyuiop", "checking line 1", line)
        call file%read_line(line, stat=stat)
        call check(stat==0, "reading line 2", stat%string())
        call check(line=="1234567890", "checking line 2", line)
    end subroutine

    subroutine test_errors
#ifdef HAVE_LIBBZ2
        use dk_bzip2filereader, only: Bzip2FileReader
!------
        type(Bzip2FileReader) :: reader
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
        call check(stat=="The stream was opened for compression", "reading from a file opened for writing", stat%string())
        call stat%discard

        call reader%close_file
        call reader%open_file("newfile.txt", "read", stat)
        call reader%write(line, n, stat)
        call check(stat=="The stream was opened for decompression", "writing to a file opened for reading", stat%string())
        call stat%discard
#endif
    end subroutine

    subroutine test_bad_files
#ifdef HAVE_LIBBZ2
        use dk_bzip2filereader, only: Bzip2FileReader
!------
        type(Bzip2FileReader), allocatable :: reader
        type(FException) :: stat
        character(:), allocatable :: line
        integer :: n
!------
        allocate(character(100) :: line)

        call start_test_group("malformed files")

        allocate(reader)
        call reader%open_file("magic.bz2", "r", stat)
        call reader%read(line, n, stat)
        deallocate(reader)
        call check(stat=="The file does not contain valid bzip2-compressed data", "wrong magic number", stat%string())
        call stat%discard

        allocate(reader)
        call reader%open_file("corrupted.bz2", "r", stat)
        call reader%read(line, n, stat)
        deallocate(reader)
        call check(stat=="Data integrity error", "corrupted file", stat%string())
        call stat%discard

        allocate(reader)
        call reader%open_file("incomplete.bz2", "r", stat)
        call reader%read(line, n, stat)
        deallocate(reader)
        call check(stat=="Trying to read past end of file.", "incomplete file", stat%string())
        call stat%discard
#endif
    end subroutine
end program
