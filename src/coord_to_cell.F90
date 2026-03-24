module coord_to_cell_mod
  USE gen_com_m, only:uwrt,lwrt, low_limit,zero,lperiod
  use cryst_to_cart_mod,only:cryst_to_cart
  use boxconfig,only: box_config
  use notperiod_mod,only:notperiod
   
  implicit none
contains
  subroutine coord_to_cell(tab_coord, cell,box,nox,noy,noz)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double


    implicit none

    !
    ! Cette routine retourne dans cell le numero de cellule
    ! contenant les coordonnees tab_coord
    !
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double)  :: tab_coord(3)
    class(box_config)::box
    integer       :: cell,nox,noy,noz
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: kx, ky, kz,k
    real(double) :: aux, auy, auz,cpp,xpici
!    real(double) :: coord_loc(3)   ! permet de ne pas ecraser coord
    ! lors de l'appel a cryst_to_cart                       
    !-----------------------------------------------
    real(double) :: xpnp(3,1) !
    call notperiod(1,tab_coord,xpnp,box%at,box%bg,lperiod)




    call cryst_to_cart (1,  xpnp(:,1), box%bg, -1) !cart vers cryst
    do k=1,3
       xpici= xpnp(k,1)
       if ( (xpici < 0.d0 ).OR.( xpici >= 1.d0 ) ) then
          if ( (xpici > -low_limit).and.(xpici<0.d0) ) then
              xpnp(k,1)=zero
          else
             cpp  = Dble(Floor(xpnp(k,1)))
             xpnp(k,1) = xpici     - cpp
          end if
       end if
    end do



    aux = xpnp(1,1)*nox
    auy = xpnp(2,1)*noy
    auz = xpnp(3,1)*noz

    kx = int(aux)
    ky = int(auy)
    kz = int(auz)

    cell = 1+kx+nox*(ky+noy*kz)

    return
  end subroutine coord_to_cell
end module coord_to_cell_mod
