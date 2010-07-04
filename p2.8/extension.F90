subroutine Pextension(it,extension,lenfn)
  implicit none
  integer::lenfn,it
  character*9::extension
  write(6,*)'Pextension it',it
    open(unit=17, file='tampon', form='formatted', status='unknown')
    if (it<=9) then
       write (17, '(I1)') it
       rewind 17
       read (17, 101) extension
       write(6,*)extension
    end if
    if (it<=99.and.it>9) then
       write (17, 200) it
       rewind 17
       read (17, 201) extension
       write(6,*)extension
    end if
    if (it<=999.and.it>99) then
       write (17, 300) it
       rewind 17
       read (17, 301) extension
    end if
    if (it<=9999.and.it>999) then
       write (17, 400) it
       rewind 17
       read (17, 401) extension
    end if

    if (it<=99999.and.it>9999) then
       write (17, 500) it
       rewind 17
       read (17, 501) extension
    end if
    if (it<=999999.and.it>99999)  then
       write (17, 600) it
       rewind 17
       read (17, 601) extension
    end if
    if (it<=9999999.and.it>999999)  then
       write (17, 700) it
       rewind 17
       read (17, 701) extension
    end if
    if (it<=99999999.and.it>9999999)  then
       write (17, 800) it
       rewind 17
       read (17, 801) extension
    end if
    if (it<=999999999.and.it>99999999)  then
       write (17, 900) it
       rewind 17
       read (17, 901) extension
    end if
    if (it>999999999) then
       write (6, *) 'probleme de format dans rasmol.f'
       stop
    endif


    lenfn=index(extension,' ')-1
    write(6,*)'Pextension = ',extension,lenfn
134 format(i6)
135 format(A,3f10.4,I7)
136 format(A,3f10.4,D14.5,I7)

    close(17)

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
