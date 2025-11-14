program error_conditions
    use dk_exception, only: FException
    use dk_fileobject, only: FileObject
    
    implicit none(external, type)

    type(FileObject) :: file
    character(:), allocatable :: line
    type(FException) :: stat
    
    ! Read from an uninitialise file
    call start_test("Testing reading from an unopened file: ")
    call file%read_line(line, stat=stat)
    call validate_test(stat, "File not open.")
    
    ! Write to an uninitialised file
    call start_test("Testing writing to an unopened file: ")
    line = ""
    call file%write_line(line, stat=stat)
    call validate_test(stat, "File not open.")
    
    ! Open a non-existing file for reading
    call start_test("Testing reading from a non-existing file: ")
    call file%init("invalid", "r", stat=stat)
    call validate_test(stat, "No such file or directory.")

    ! Open a file in a non-existing directory for reading
    call start_test("Testing reading from a file in a non-existing directory: ")
    call file%init("invalid/test.txt", "r", stat=stat)
    call validate_test(stat, "No such file or directory.")

    ! Open a file in a non-existing directory for writing
    call start_test("Testing writing to a file in a non-existing directory: ")
    call file%init("invalid/test.txt", "w", stat=stat)
    call validate_test(stat, "No such file or directory.")

    call start_test("Testing reading to a file opened for writing: ")
    call file%init("test.txt", "w", stat=stat)
    call assert_success(stat)
    call file%read_line(line, stat=stat)
    call validate_test(stat, "File not open for reading.")
    
!    call file%write_line("42")
contains
    subroutine start_test(message)
        character(*), intent(in) :: message
        write(*,'(a)',advance="no") message
    end subroutine

    subroutine assert_success(stat)
        type(FException), intent(inout) :: stat
!------
        if (stat /= 0) then
            call stat%report
            error stop 'Assertion failed'
        end if
        call stat%discard
    end subroutine

    subroutine validate_test(stat, message)
        use ext_character, only: colour
        use dk_exception, only: FException, EVENT_LEVEL_VERBOSE
!------        
        type(FException), intent(inout) :: stat
        character(*), intent(in) :: message
!------
        if (stat == message) then
            write(*,'(a)') colour("ok", "green")
            call stat%setLevel(EVENT_LEVEL_VERBOSE)
            call stat%report
            call stat%discard
            write(*,*)
        else
            write(*,'(a)') colour("failed", "red")
            call stat%report
            error stop
        end if
    end subroutine
end program
