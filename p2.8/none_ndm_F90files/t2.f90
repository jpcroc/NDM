program grep
  integer ::i
  character*30 sub (1000) 
  open(unit=10, file='ls4')

  do i=1,1000
	write(6,*)'echo'
	write(6,*)'echo'
     read(10,*)sub(i)
     write(6,*)sub(i)
!     write(6,*)'grep ', sub, ' *.F90 *.f |grep call' 
  end do
end program grep
