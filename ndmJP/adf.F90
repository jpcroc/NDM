
subroutine adf
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m
  use tab_imm_m


  !******************************************************************
  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i1, i2, i3, i4, koo, ko1, ko2, i, &
       j, k, m,m1,lutriplet,ka,ma,m2

  real(double) :: thetaijk, c11, c21, c31, c12, c22, c32, &
       incre
  real(double) :: rspace2, invincre, &
       costheta
  real(double) :: aaa,bbb,ccc
  real(double), dimension(ntyp,ntyp,ntyp) :: intfda,angle
  character :: triplet1*20,triplet2*20,triplet3*20
  integer :: lenftriplet1,lenftriplet2,lenftriplet3
  character :: ftriplet1*80, ftriplet2*80, ftriplet3*80
  character :: ftripletangle*80,ftrip*80
  integer :: lusauvfda
  character :: charsauvfda*8
  !-----------------------------------------------

  incre = (thetamax-thetamin)/cont
  invincre = 1/incre

  lutriplet = 13
  if (rang==0) then
     write (6, *)
     write(6,*) '-----------------------------------------'
     write (6, *) '------- Calcul des Distributions angulaires --------'
     write(6,*)

     write (6, '(A,I5,A,D10.3)') '*  ITERATION  = ', it, '  time = ', timel
  end if

  ! On calcule maintenant la distribution angulaire moyennee

  if(.not.linstantfda) then
     do i3=1,ntyp
        if(nad(i3)==0) cycle

        do i1=1,ntyp
           if(nad(i1)==0) cycle

           do i2=i3,ntyp
              if(nad(i2)==0) cycle
              triplet1= ty(i1)
              triplet2= ty(i2)
              triplet3= ty(i3)
              ftriplet1=triplet1
              ftriplet2=triplet2
              ftriplet3=triplet3
              lenftriplet1=index(ftriplet1,' ')-1
              lenftriplet2=index(ftriplet2,' ')-1
              lenftriplet3=index(ftriplet3,' ')-1
              ftripletangle = ftriplet2(1:lenftriplet2)// &
                   '_'//ftriplet1(1:lenftriplet1)// &
                   '_'//ftriplet3(1:lenftriplet3)//'.angle'

              ftrip = ftriplet2(1:lenftriplet2)// &
                   '_'//ftriplet1(1:lenftriplet1)// &
                   '_'//ftriplet3(1:lenftriplet3)


              do m1=int(thetamin*invincre),cont
                 m2 = m1+1
                 rspace2 = (m2*incre)**2
                 aaa = fda(i2,i1,i3,m2)/nfda
                 intfda(i2,i1,i3)=intfda(i2,i1,i3)+aaa
                 angle(i2,i1,i3)=angle(i2,i1,i3)+aaa*m1*incre

                 if(rang==0) then
                    write (13,*) m2*incre*180/pi, aaa
                 endif
              end do


              if(angle(i2,i1,i3).gt.0.1) then
                 if (rang==0) open(lutriplet, file = ftripletangle, status = 'unknown')
                 angle(i2,i1,i3)=angle(i2,i1,i3)/intfda(i2,i1,i3)
                 if(rang.eq.0) write(6,'("Angle moyen du triplet",1X,A8,"=",1x,f7.3)')& 
                      ftrip,angle(i2,i1,i3)*180/pi
                 do m1=int(thetamin*invincre),cont
                    m2 = m1+1
                    rspace2 = (m2*incre)**2
                    aaa = fda(i2,i1,i3,m2)/nfda
                    if(rang==0) then
                       write (13,*) m2*incre*180/pi, aaa
                    endif
                 end do

              endif
           end do    !i2
        end do
     end do

     if(rang==0)         close(lutriplet)
     if(rang==0) then
        write(6,*)'--------------------------------------'
        write(6,*) '--------------------------------------'
     end if
  else
     if (rang==0) then
        open(unit=32,file='tampfda',form='formatted',status='unknown')

        if (it<=9) write (32, '(I1)') it
        if (it<=99.and.it>9) write (32, 200) it
        if (it<=999.and.it>99) write (32, 300) it
        if (it<=9999.and.it>999) write (32, 400) it
        if (it<=99999.and.it>9999) write (32, 500) it
        if (it<=999999.and.it>99999) write(32, 600) it
        if (it<=9999999.and.it>999999) write(32, 700) it
        if (it<=99999999.and.it>9999999) write(32, 800) it
        if (it<=999999999.and.it>99999999) write(32, 900) it
        if  (it>999999999) then
           write (6, *) 'probleme de format dans calcangle.f90'
           stop
        endif
        rewind 32

        read (32, 1000) charsauvfda

        lusauvfda=31
     end if
     if (rang==0) then
        write(6,*) ' sauvegarde fda it=',it
        write(6,*)
     end if
     do i3=1,ntyp
        if(nad(i3)==0) cycle
        do i1=1,ntyp
           if(nad(i1)==0) cycle
           do i2=1,ntyp
              if(nad(i2)==0) cycle

              triplet1= ty(i1)
              triplet2= ty(i2)
              triplet3= ty(i3)
              ftriplet1=triplet1
              ftriplet2=triplet2
              ftriplet3=triplet3
              lenftriplet1=index(ftriplet1,' ')-1
              lenftriplet2=index(ftriplet2,' ')-1
              lenftriplet3=index(ftriplet3,' ')-1
              ftripletangle = ftriplet2(1:lenftriplet2)// &
                   '_'//ftriplet1(1:lenftriplet1)// &
                   '_'//ftriplet3(1:lenftriplet3)//'.'//charsauvfda
              ftrip = ftriplet2(1:lenftriplet2)// &
                   '_'//ftriplet1(1:lenftriplet1)// &
                   '_'//ftriplet3(1:lenftriplet3)

              if (rang==0) open(unit=lusauvfda, file = ftripletangle, status ='unknown')

              do m1=int(thetamin*invincre),cont
                 m2 = m1+1
                 rspace2 = (m2*incre)**2
                 intfda(i2,i1,i3)=intfda(i2,i1,i3)+fda(i2,i1,i3,m2)
                 angle(i2,i1,i3)=angle(i2,i1,i3)+fda(i2,i1,i3,m2)*&
                      m1*incre
                 if(rang==0) then
                    write (lusauvfda,*) m2*incre*180/pi, fda(i2,i1,i3,m2)
                 endif
              end do

              angle(i2,i1,i3)=angle(i2,i1,i3)/intfda(i2,i1,i3)

              if(rang==0) then
                 write(6,'("Angle instantane moyen du triplet",1x,A8,"=",1x,f7.3)')&
                      ftrip,angle(i2,i1,i3)*180/pi
              endif
           enddo
        enddo
     enddo

200  format(i2)
300  format(i3)
400  format(i4)
500  format(i5)
600  format(i6)
700  format(i7)
800  format(i8)
900  format(i9)
1000 format(a6)
     if (rang==0)then
        write(6,*)'--------------------------------------'
        write(6,*) '--------------------------------------'
        close(lusauvfda)
        close(32,status='DELETE')
     endif
  endif

  fda(:,:,:,:)=0.D0
  intfda(:,:,:)=0.D0
  angle(:,:,:)=0.D0
  return
end subroutine adf

