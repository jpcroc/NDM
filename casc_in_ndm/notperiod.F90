

! *****************************************************************
subroutine notperiod(xp, xpnp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m

  !       version du 09 decembre 2003

  ! *****************************************************************
  ! Cette routine applique les conditions périodiques par décalage
  !  +/-1 (triclinique).
  !Le décalage est fait sur xp xp ET ax. Ce décalage sur ax 
  !permet de mesurer correctement le déplacements des atomes
  ! depuis leur position de départ

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(double)  :: xp(3,imm)
  real(double)  :: xpnp(3,imm)

  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, ic,icp
  real(double)::dz,trav,ecav,ecap
  real(double):: cpp,xpici,ctest
  integer,save  :: iperiod


  !      iperiod=iperiod+1
  !      
  !      write(*,*) 'PBC PBC PBC capitala tarii e ....          notperiod',iperiod
  !
  xpnp(:,:)=xp(:,:)


     call cryst_to_cart (imm, xpnp,  bg,  -1) !cart vers cryst

     do i=1,im
        do ic=1,3
           xpici=xpnp(ic,i)

           if ( (xpici < 0.d0 ).OR.( xpici >= 1.d0 ) ) then
              if ( (xpici > -low_limit).and.(xpici<0.d0) ) then
                 !	       write(*,*) 'low',low_limit,xpici
                 xpnp(ic,i)=zero
              else
                 xpnp (ic,i) = xpici - Dble(Floor(xpici))
              end if
           end if

        end do
     end do

     call cryst_to_cart (imm, xpnp, at, 1)  !cryst vers cart





end subroutine notperiod

