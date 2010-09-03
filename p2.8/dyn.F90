! *************************************************************
subroutine dyn
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use jqmod
   use tab_imm_m

#if(PARA)
  use mod_mpi
#endif

  implicit none
  !----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, iti, ic
  real(double), dimension(ntyp) :: aux
  real(double) :: xprov ,vv
  real(double), save :: tmoyinst, imesureT
  real(double), external :: tempinst
  real(double) :: tempavant,tmoy
  real(double) :: unitE,deltaE
  character*5 :: cunitE
  real(double):: eatommoy
  if(lEev) then
     unitE=erg2eV
     cunitE='  eV '
  else
     unitE=1.0
     cunitE=' erg '
  end if
  if (associated(eatom))  eatom(:)=eatom(:)+0.5*cm(ityp(:))*(vp(1,:)**2+vp(2,:)**2+vp(3,:)**2)

  if (lnemd) then
     eatommoy=0.
     do i=1,imd
!        eatom(i)=eatom(i)+0.5*cm(ityp(i))*(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)
        eatommoy=eatommoy+eatom(i)/float(imd)
     end do
     do i=1,imd

        fp(1,i)=fp(1,i)+(eatom(i)-eatommoy)*Fnemd
     end do
!     do i=1,imd
!        eatom(i)=eatom(i)-0.5*cm(ityp(i))*(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)
!     end do
  end if



  timel = timel+tstep
  aux(:ntyp) = tstep**2/cm(:ntyp)
  !      write(6,*)'aux ',aux


  IF (lFrozen) THEN     ! Some atoms are frozen
          do i = 1, imd
             IF (Free(i)) THEN  ! This atom is free to move
                     do ic = 1, 3
                        xprov = (xp(ic,i)-xpp(ic,i))+xp(ic,i)+aux(ityp(i))*fp(ic,i)
                        vp(ic,i) = (xprov-xpp(ic,i))*usdh
                        xpp(ic,i) = xp(ic,i)
                        xp(ic,i) = xprov
                     end do
             ELSE       ! This atom is frozen
                     vp(1:3,i) = aux(ityp(i))*fp(1:3,i)*usdh
                     xpp(1:3,i) = xp(1:3,i)
             END IF
          end do
  ELSE                  ! All atoms can move
          do i = 1, imd
             do ic = 1, 3
                xprov = (xp(ic,i)-xpp(ic,i))+xp(ic,i)+aux(ityp(i))*fp(ic,i)
                vp(ic,i) = (xprov-xpp(ic,i))*usdh
                xpp(ic,i) = xp(ic,i)
                xp(ic,i) = xprov
             end do
          end do
  END IF

  if (lcalcjq) then
     eatommoy=0.
     jqp=jq ; jqk=0.0 !; expvect(:)=0.0
     do i=1,imd
!        eatom(i)=eatom(i)+0.5*cm(ityp(i))*(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)
        if(lnemd) eatommoy=eatommoy+eatom(i)/float(imd)
        expvect(:)=expvect(:)+eatom(i)*xpp(:,i)
        do ic=1,3              
           jqk(ic)=jqk(ic)+eatom(i)*vp(ic,i)
           !               jq(ic)=jq(ic)+eatom(i)*vp(ic,i)
        end do
     end do
     jq=jqp+jqk



     write(65,'(I8,3D16.8)')it-1,jq(1),jq(2),jq(3)
     !         write(66,'(I8,3D15.6)')it-1,expvect(1),expvect(2),expvect(3)
  end if


  ! temp const moyenne depuis le dernier rescale
  if(lTcon) then
     tempavant=tempinst(vp,ityp)
     tmoyinst = tmoyinst +tempavant 
     imesureT = imesureT + 1
     if(mod(it,iteTconst).eq.0) then
        write(6,*)
        tmoy=tmoyinst/imesureT 
        write(6,'(A,3F12.2)')'> Temperature Constante = Text, Tmoy, Tinst ', Text,tmoy,tempavant


        deltaE=(tmoy-Text)*3*bk*unitE
        write(6,'(A,I0,D21.12,A)')'IT   modification d_energie par atome ',it,deltaE,&
     &cunitE
!       write(6,'(A,I,D21.12,A)')'IT , modification d_energie par atome ',IT,deltaE, cunitE
        write(6,*)

        vv = sqrt(Text/(tmoy))
        xpp(:,:im) = xp(:,:im)-(xp(:,:im)-xpp(:,:im))*vv
        vp(:,:im) = vp(:,:im)*vv
        tmoyinst=0.
        imesureT=0
     endif
  endif

  if (lperiod)       call period 





  return
end subroutine dyn

