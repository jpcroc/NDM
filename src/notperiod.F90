module notperiod_mod
  USE cryst_to_cart_mod,only: cryst_to_cart
  implicit none
contains

  ! *****************************************************************
  subroutine notperiod(im,xp, xpnp,at,bg,lperiod)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:low_limit,zero

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
    integer,intent(in)::im
    real(double),intent(in)  :: xp(3,im),at(3,3),bg(3,3)
    real(double)  :: xpnp(3,im)
    logical, intent(in)::lperiod

    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: i, ic
    real(double):: xpici


    !      iperiod=iperiod+1
    !      
    !      write(*,*) 'PBC PBC PBC capitala tarii e ....          notperiod',iperiod
    !

    xpnp(:,:)=xp(:,:)

    if (lperiod) then
       return
    else 
       call cryst_to_cart (im, xpnp,  bg,  -1) !cart vers cryst

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

       call cryst_to_cart (im, xpnp, at, 1)  !cryst vers cart

    end if



  end subroutine notperiod
end module notperiod_mod
