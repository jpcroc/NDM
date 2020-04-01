subroutine Pextension(it,extension,lenfn)
  implicit none
  integer::lenfn,it
  character*9::extension


    lenfn = 9
    write(extension,'(i9.9)') it

134 format(i6)
135 format(A,3f10.4,I7)
136 format(A,3f10.4,D14.5,I7)



200 format(i2)
300 format(i3)
400 format(i4)
500 format(i5)
600 format(i6)
700 format(i7)
800 format(i8)
900 format(i9)
101 format(a1)
201 format(a2)
301 format(a3)
401 format(a4)
501 format(a5)
601 format(a6)
701 format(a7)
801 format(a8)
901 format(a9)
 
  end subroutine Pextension
