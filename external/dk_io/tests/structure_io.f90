program test_structure_io
    use dk_test

    implicit none(external, type)

    integer, parameter :: e = kind(0.0d0)

    call test_array_io

    call tests_finished
contains
    subroutine test_array_io
        use dk_structure_io, only: read_array
        use dk_token, only: Token, get_token

        real(e), dimension(5) :: data
        character(:), allocatable :: line
        type(Token) :: t


        integer :: n

        line = "1.0394 0.30132 5.2304"
        call read_array(line, data(:3), n)
        write(*,*)
        write(*,*) "Test 1"
        write(*,*) "Read: ", n
        write(*,*) "Input: ", line
        write(*,*) "Result:", data(:3)

        line = "1.0394 0.30132 5.2304"
        call read_array(line, data(:4), n)
        write(*,*)
        write(*,*) "Test 2"
        write(*,*) "Read: ", n
        write(*,*) "Input: ", line
        write(*,*) "Result:", data(:4)

        line = "1.0394 wrong 5.2304"
        call read_array(line, data(:3), n)
        write(*,*)
        write(*,*) "Test 3"
        write(*,*) "Read: ", n
        write(*,*) "Input: ", line
        write(*,*) "Result:", data(:3)

        call get_token("1.0394 wrong 5.2304", n+1, t)
        write(*,*) "Wrong token: ", t%value

        line = "1.0394 0.30132 wrong"
        call read_array(line, data(:3), n)
        write(*,*)
        write(*,*) "Test 4"
        write(*,*) "Read: ", n
        write(*,*) "Input: ", line
        write(*,*) "Result:", data(:3)

        line = "1.0394 0.30132 5.2304 junk"
        call read_array(line, data(:3), n)
        write(*,*)
        write(*,*) "Test 5"
        write(*,*) "Read: ", n
        write(*,*) "Input: ", line
        write(*,*) "Result:", data(:3)

        line = "1.0394 0.30132 5.2304junk"
        call read_array(line, data(:3), n)
        write(*,*)
        write(*,*) "Test 6"
        write(*,*) "Read: ", n
        write(*,*) "Input: ", line
        write(*,*) "Result:", data(:3)

    end subroutine
end program
