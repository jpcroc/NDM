module caltabt_mod
        use notperiod_mod
        use cryst_to_cart_mod
#ifdef PARA
        use mod_para
#endif      
        implicit none
        contains
! ******************************************************************
subroutine caltabt
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
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
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, ic, icell, kx, ky, kz, koo
  real(double) :: aux, auy, auz
  
  real(double), dimension(:,:), allocatable :: xpnp !
  !
  ! --------- Initialisation --------------
  !
!   write(6,*)'caltabt',it
   nato(:noxyz) = 0
  last(natperc,:noxyz) = 0

!     do i = 1, im
!     if ((it.ge.1000).and.(i.lt.20)) write(6,'(I5,3G15.7)')i, xp(1,i),xp(2,i),xp(3,i)
!  end do
  !  -------- cas sans cellule  -----------
  if (noxyz==1) then
          nato(1) = im
          do i = 1, im
             ielat(i) = 1
             last(i,1) = i
          end do
  else

          ALLOCATE(xpnp(3,imm))        
          if (lperiod) then            
             xpnp(:,:)=xp(:,:)         
          else                         
             call notperiod(xp,xpnp)   
          end if                       
     !  -------- Initialisations  -----------
     nato(0:noxyz) = 0
     last(natperc,:noxyz) = 0

     ! -------------------------------------------
     !   1. loop: lattice-coordinates of all atoms
     ! - - - - - - - - - - - - - - - - - - - - - -

     !debug       write (*,*) 'sub caltabt 1',it,xp(1,1)
     call cryst_to_cart (imm, xpnp, bg, -1) !cart vers cryst
     !debug       write (*,*) 'sub caltabt 2',it,xp(1,1)

!     if (it.gt.1000) write(6,*)'CALTABT',it
     do i = 1, im
!     if  ((it.ge.1000).and.(i.lt.20)) write(6,'(I5,3G15.7)')i, xpnp(1,i),xpnp(2,i),xpnp(3,i)
        aux = xpnp(1,i)*nox
        auy = xpnp(2,i)*noy
        auz = xpnp(3,i)*noz
        kx = int(aux)
        ky = int(auy)
        kz = int(auz)
        !write(*,*) i, nox,noy,noz, kx,ky,kz
        kx = Modulo(kx,nox)
        ky = Modulo(ky,noy)
        kz = Modulo(kz,noz)
!          if  ((it.ge.1000).and.(i.lt.20))  write(6,'(I5,3G15.7)')i, kx,ky,kz
        !==============================================================
        koo = 1+kx+nox*(ky+noy*kz)

        IF ( (koo.GT.noxyz).OR.(koo.LT.0) ) THEN
           WRITE(0,'(a,i0,a,3g20.12)') &
                'Problem with atom ', i, ', x,y,z = ', xp(1:3,i)
           WRITE(0,'(2(a,i0))') ' koo = ', koo, ' - noxyz = ', noxyz
           STOP
        END IF
        ielat(i) = koo
        nato(koo) = nato(koo)+1
        !!$write(*,*) MAXVAL(nato(:)),koo,i        ! DEBUG
        ! ==== MODIF Clouet =====================
        IF (nato(koo).GT.natperc) THEN
                WRITE(0,'(a)') 'You need to increase the maximal number of atoms per cell'
                WRITE(0,'(a,i0)') 'current value: natperc=', natperc
                STOP '< Caltabt >'
        END IF
        ! ==== Fin MODIF Clouet =================
        last(nato(koo),koo) = i
     end do
     !debug            call cryst_to_cart (imm, xpnp, at, 1)  !cryst vers cart


!        open(unit=809, file='cell.csv', form='formatted', &
!             status='unknown')
!    do i=1,im_glob
!       write(809,'(A,3I6,3G22.13)')'Cel ', it,i,ielat(i),xp(:,i)
!       write(809,'(A,3I6,3G22.13)')'Cel ', it,i,ielat(i),xpnp(:,i)
!    end do
 !             do i=1,noxyz
 !                write(6,*) i, nato(i) 
 !             end do
     DEALLOCATE(xpnp)   ! MODIF CLOUET
  endif

  !!$write(6,*)'sortie caltabt'     ! DEBUG

!  write(6,*)'maxnato', maxval(nato)

  return
end subroutine caltabt
end module
