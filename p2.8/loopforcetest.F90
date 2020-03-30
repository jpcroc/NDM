module loopforcetest_mod
        use calfo_mod
        implicit none
        contains
! ************************************************
!           Sous-programme dmloop.f
!          Version MPI du 21 fevrier 2001
! ************************************************

subroutine loopforcetest(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer  :: ielat(imm)
  integer  :: iwmax(imm)
  integer  :: ityp(imm)
  real(double)  :: xp(3,imm)
  real(double)  :: xpp(3,imm)
  real(double)  :: vp(3,imm)
  real(double)  :: ax(3,imm)
  real(double)  :: fp(3,imm)
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, iti,ic,is,test_force
  real(double) :: epot0,deltaE,deltaf1,fps(3,imm)
  !-----------------------------------------------
  !

  unitE=1.0
  cunitE=' erg'

  unitP=1.0
  cunitP='d/cm2'


  ! MPI
  if (rang==0) write (6, *) '***** test des forces  ****'

  imd = im
  nad(:ntyp) = na(:ntyp)




test_force=2

  select case (test_force)
  case(1)


     !      write(6,*)
     !      write(6,*)'***** ITERATION  ****', it

     ! appel de la routine generale des forces
     !position de départ
     vp=0.0
     do i=1,im
        write(6,*)i,xp(1,i),xp(2,i),xp(3,i)
     end do
     fp=0.
     call calfo
     write (6, '(A,D21.12)') '*Epot = ', potist
     epot0=potist
     write(6,*)'forces'
     do i=1,im
        write(6,'(I2,3D21.12)')i,fp(1,i),fp(2,i),fp(3,i)
     end do
     fps=fp


     do i=1,im
        do ic=1,3
           !ic =1
           do is=1,-1,-2
              xp(ic,i)=xp(ic,i)+is*deltax
              write(6,*)
              write(6,*) 'i,X is', i, ic,is

              call calfo
              deltaE=potist-epot0
              write (6, '(A,D21.12,A,D21.12)') '*Epot = ', potist,' deltaE= ',deltaE
              !                  deltaf1= (-1.*is*deltaE/deltax-fps(ic,i))/fps(ic,i)
              deltaf1= (-1.*is*deltaE/deltax)/fps(ic,i)
              write (6, '(A,3D18.10)') '*-deltaE/deltax ', -1.*is*deltaE/deltax,fps(ic,i),deltaf1
              !write (1, '(A,3D18.10)') '',xp(ic,i), deltaE/(is*deltax),fps(ic,i)
              xp(ic,i)=xp(ic,i)-is*deltax
           end do
        end do
     end do

     stop
     !close(1)
     return
  case(2)
     do i=1,im
        write(6,*)i,xp(1,i),xp(2,i),xp(3,i)
     end do
     fp=0.
     call calfo
     write(789,*)(xp(1,2)-xp(1,1))*1.d8,potist


     do while (xp(1,1).lt.xp(1,2))
        xp(1,1)=xp(1,1)+deltax
        call calfo
        write(789,*)(xp(1,2)-xp(1,1))*1.d8,potist*erg2ev
     end do
  end select

end subroutine loopforcetest
end module
