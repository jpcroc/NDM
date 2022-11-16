! afile=test_010_hello_world ; gfortran -o ${afile}.exe ${afile}.f90 && ${afile}.exe
! just for test...

program main
  implicit none

  write (*, "(a)") "hello world !"
end program
