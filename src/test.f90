module dk_test
    use ext_character

    use dk_exception, only: FException, EVENT_LEVEL_WARNING

    implicit none
    private
    save

    integer, private :: expected_passed = 0
    integer, private :: unexpected_passed = 0
    integer, private :: unexpected_failed = 0
    integer, private :: expected_failed = 0

    interface check
        module procedure check_condition, check_status
    end interface

    public :: tests_finished
    public :: check
    public :: xfail
    public :: start_test_group
contains
!******************************************************************************!
!> Print the header for a test group
!******************************************************************************!
    subroutine start_test_group(group_name, description)
        character(*), intent(in) :: group_name
        character(*), intent(in), optional :: description
!------
        write(*,*)
        write(*,*) "Testing "//group_name
        if (present(description)) then
            write(*,*) "  "//description
            write(*,*)
        end if
    end subroutine
!******************************************************************************!
!> Test a condition, and print a message indicating the test result
!******************************************************************************!
    subroutine check_condition(test, name, description)
        logical, intent(in) :: test
        character(*), intent(in) :: name
        character(*), intent(in), optional :: description
!------
        if (test) then
            write(*,*) "["//colour("passed", "green")//"] "//name
            expected_passed = expected_passed + 1
            return
        end if

        unexpected_failed = unexpected_failed + 1
        write(*,*) "["//colour("failed", "red")//"] "//name
        if (present(description)) write(*,*) description
    end subroutine
!******************************************************************************!
!> Test a condition, and print a message indicating the test result and the
!> error message if needed.
!******************************************************************************!
    subroutine check_status(test, name, description)
        class(FException), intent(inout) :: test
        character(*), intent(in) :: name
        character(*), intent(in), optional :: description
!------
        if (test==0) then
            write(*,*) "["//colour("passed", "green")//"] "//name
            expected_passed = expected_passed + 1
            return
        end if

        unexpected_failed = unexpected_failed + 1
        write(*,*) "["//colour("failed", "red")//"] "//name
        if (present(description)) write(*,*) description
        call test%setlevel(EVENT_LEVEL_WARNING)
        call test%report
        call test%discard
    end subroutine
!******************************************************************************!
!> Test a condition, and print a message indicating the test result, for tests
!> that are expected to fail
!******************************************************************************!
    subroutine xfail(test, name, description)
        logical, intent(in) :: test
        character(*), intent(in) :: name
        character(*), intent(in), optional :: description
!------
        if (.not.test) then
            expected_failed = expected_failed + 1
            write(*,*) "["//colour("xfail", "yellow")//"] "//name
            return
        end if

        unexpected_passed = unexpected_passed + 1
        write(*,*) "["//colour("xpass", "red")//"] "//name
        if (present(description)) write(*,*) description
    end subroutine
!******************************************************************************!
!> Print a summary of the tests
!******************************************************************************!
    subroutine tests_finished
        write(*,*)
        write(*,'(a,i0)') "Tests summary:"
        write(*,'(a,i0)') "Passed:            ", expected_passed
        write(*,'(a,i0)') "Expected failures: ", expected_failed
        write(*,'(a,i0)') "Failed:            ", unexpected_failed
        write(*,'(a,i0)') "Unexpected passes: ", unexpected_passed
        write(*,'(a,i0)') "Total:             ", (expected_passed+expected_failed+unexpected_passed+unexpected_failed)

        ! Return a non-zero error code if at least one test has failed
        if (unexpected_failed /= 0) error stop
    end subroutine
end module

