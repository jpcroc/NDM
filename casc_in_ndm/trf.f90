program grep
  character sub*30
  open(unit=10, file='ls4')
  do
     write(6,*)'echo'
     write(6,*)'echo'
     read(10,*)sub
     write(6,*)'grep ', sub, ' *.F90 *.f |grep call' 
  end do
end program grep
