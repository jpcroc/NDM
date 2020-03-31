module tabcr

  integer,pointer :: lastcr (:,:),natocr(:),ielatcr(:)

contains
! ******************************************************************
subroutine caltabtcr (natperc,nox,noy,noz,xpcr,im,imm,bg,at)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use tab_imm_m
  !          Version du 01 fevrier 2001
  ! ******************************************************************

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------

  integer,intent(in)::natperc,nox,noy,noz,im,imm
  real(double),intent(in)  :: xpcr(3,imm),bg(3,3),at(3,3)
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, ic, icell, kx, ky, kz, koo,noxyz
  real(double) :: aux, auy, auz
  !  real(double)::xpnp(3,imm)
  !-----------------------------------------------
  !
  ! --------- Initialisation --------------
 !     write(6,*)'entree caltabtcr'


  noxyz=nox*noy*noz
  !       write (*,*) 'entree caltabt >>>>>>>>>>>>>>>>>'

  !  if (lperiod) then
  !     xpnp(:,:)=xpcr(:,:)
  !  else 
  !     call notperiod(xp,xpnp)
  !  end if

  !      write(6,*)'entree caltabt noxyz',noxyz
  !  -------- cas sans cellule  -----------
  if (noxyz==1) then
     natocr(1) = im
     do i = 1, im
        ielatcr(i) = 1
        lastcr(i,1) = i
        !      write(6,*) 'last,i,im=',last(i,1),i,im
     end do
  else

     !  -------- Initialisations  -----------
     natocr(0:noxyz) = 0
     lastcr(natperc,:noxyz) = 0

     ! -------------------------------------------
     !   1. loop: lattice-coordinates of all atoms
     ! - - - - - - - - - - - - - - - - - - - - - -

     !debug       write (*,*) 'sub caltabt 1',it,xp(1,1)
     call cryst_to_cart (imm, xpcr, bg, -1) !cart vers cryst
     !debug       write (*,*) 'sub caltabt 2',it,xp(1,1)


     do i = 1, im
        aux = xpcr(1,i)*nox
        auy = xpcr(2,i)*noy
        auz = xpcr(3,i)*noz
        kx = int(aux)
        ky = int(auy)
        kz = int(auz)
        !               write(6,*)
        !               write(6,*)xp(1,i),xp(2,i),xp(3,i)
        !               write (6,*)'kxyzA',i,kx,ky,kz
        !==============================================================
        ! Ajout Emmanuel au cas ou les coordonnÃÂ©es rÃÂ©duites ne sont pas
        ! comprises entre 0 et 1 (en particulier pour lperio=.FALSE.)
        kx = Modulo(kx,nox)
        ky = Modulo(ky,noy)
        kz = Modulo(kz,noz)
        !==============================================================
        koo = 1+kx+nox*(ky+noy*kz)
        !	       if ((i.lt.5).or.(i.gt.(im-5))) then
        !	       write(*,*) bg
        !	       write(*,*) i, xp(1,i),xp(2,i),xp(3,i)
        !	       write(*,'(i5,3f12.5,6i5,i7)') i,aux,auy,auz,kx,ky,kz,nox,noy,noz,koo
        !               end if
        IF ( (koo.GT.noxyz).OR.(koo.LT.0) ) THEN
           WRITE(0,'(a,i0,a,3g20.12)') &
                'Problem with atom ', i, ', x,y,z = ', xp(1:3,i)
           WRITE(0,'(2(a,i0))') ' koo = ', koo, ' - noxyz = ', noxyz
           STOP
        END IF
        ielatcr(i) = koo
        natocr(koo) = natocr(koo)+1
        !debug	       write(*,*) MAXVAL(nato(:)),koo,i
        lastcr(natocr(koo),koo) = i
     end do
     !debug            call cryst_to_cart (imm, xpnp, at, 1)  !cryst vers cart

     call cryst_to_cart (imm, xpcr, at, 1)  !cryst vers cart


  endif
  !       write(6,*)'sortie caltabt'

  return
end subroutine caltabtcr
end module tabcr
