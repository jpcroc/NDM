program transf

  do i=1,10000
     read(5,*) x,rep
     if( x<0.5) rep=480
     write(6,*)x,rep
  end do
end program transf
