program test
    use iso_fortran_env, only: int64
    use dk_fileobject, only: FileObject
    use dk_exception, only: FException, EVENT_LEVEL_QUIET
    use dk_test

    implicit none(external, type)

    call test_append
    call test_read
    call test_write

    call tests_finished
contains
    subroutine test_read
        type(FileObject) :: file
        type(FException) :: stat
        character(:), allocatable :: line
!------
        call start_test_group("Read mode")

        call file%read_line(line, stat=stat)
        call check(stat=="File not open.", "uninitialised file")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call file%init("./test", "r", stat)
        call check(stat==0, "opening a file with a path")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call file%init("./test", "w", stat)
        call file%read_line(line, stat=stat)
        call check(stat=="File not open for reading.", "file open for writing")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call file%init("missing/missing", "r", stat)
        call check(stat=="No such file or directory.", "reading from a file in a non-existing directory")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call file%init("missing", "r", stat)
        call check(stat=="No such file or directory.", "opening a non-existing file")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call execute_command_line("echo 'qwertyuiop' > unreadable")
        call execute_command_line("chmod -r unreadable")
        call file%init("unreadable", "r", stat)
        call check(stat=="Cannot open file 'unreadable': Permission denied", "opening a file without permissions")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard
        call execute_command_line("rm unreadable")

        call execute_command_line("touch empty")
        call file%init("empty", "r", stat)
        call check(stat==0, "opening an empty file")
        call execute_command_line("rm empty")

        call execute_command_line("echo 'qwertyuiop' > unfinished")
        call execute_command_line("echo -n '1234567' >> unfinished")
        call file%init("unfinished", "r", stat)

        call file%read_line(line, trim=.true., stat=stat)
        call check(line=="qwertyuiop", "first line")
        write(*,*) line
        call check(file%line_number()==1, "line number (1)")
        call check(len(line)==10, "length of the first line")

        call file%read_line(line, trim=.true., stat=stat)
        call check(stat=="Trying to read past end of file.", "file without a final new line")
        write(*,*) line
        call check(line=="1234567", "second line")
        call check(file%line_number()==2, "line number (2)")
        call check(len(line)==7, "length of the second line")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard
    end subroutine

    subroutine test_write
        type(FileObject) :: file
        type(FException) :: stat
!------
        call start_test_group("Write mode")

        call file%write_line("qwertyuiop", stat=stat)
        call check(stat=="File not open.", "uninitialised file")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call file%init("./test", "w", stat)
        call check(stat==0, "opening a file with a path")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call file%init("./test", "r", stat)
        call file%write_line("qwertyuiop", stat=stat)
        call check(stat=="File not open for writing.", "file open for reading")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call file%init("missing/missing", "w", stat)
        call check(stat=="No such file or directory.", "writing to a file in a non-existing directory")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard
    end subroutine

    subroutine test_append
        type(FileObject) :: file
        type(FException) :: stat
!------
        call start_test_group("Append mode")

        call file%init("./test", "a", stat)
        call check(stat==0, "opening a file with a path")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call file%init("missing/missing", "a", stat)
        call check(stat=="No such file or directory.", "appending to a file in a non-existing directory")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard

        call execute_command_line("echo 'qwertyuiop' > unreadable")
        call execute_command_line("chmod -r unreadable")
        call file%init("unreadable", "a", stat)
        call check(stat=="Cannot open file 'unreadable': Permission denied", "opening a file without permissions")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call stat%discard
        call execute_command_line("rm unreadable")

        call execute_command_line("touch empty")
        call file%init("empty", "a", stat)
        call check(stat==0, "opening an empty file")
        call stat%setlevel(EVENT_LEVEL_QUIET)
        call stat%report
        call execute_command_line("rm empty")
    end subroutine
end program
