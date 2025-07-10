module coord_to_cell_mod
  USE gen_com_m, ONLY: low_limit,zero
  use cryst_to_cart_mod,only:cryst_to_cart
  implicit none
contains
  subroutine coord_to_cell(tab_coord, cell,bg,nox,noy,noz)
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
    real(double)  :: tab_coord(3),bg(3,3)
    integer       :: cell,nox,noy,noz
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: kx, ky, kz,k
    real(double) :: aux, auy, auz,cpp,xpici
    real(double) :: coord_loc(3)   ! permet de ne pas ecraser coord
    ! lors de l'appel a cryst_to_cart                       
    !-----------------------------------------------

    coord_loc(:) = tab_coord(:)


    call cryst_to_cart (1, coord_loc, bg, -1) !cart vers cryst
    do k=1,3
       xpici=coord_loc(k)
       if ( (xpici < 0.d0 ).OR.( xpici >= 1.d0 ) ) then
          if ( (xpici > -low_limit).and.(xpici<0.d0) ) then
             coord_loc(k)=zero
          else
             cpp  = Dble(Floor(coord_loc(k)))
             coord_loc(k) = xpici     - cpp
          end if
       end if
    end do



    aux = coord_loc(1)*nox
    auy = coord_loc(2)*noy
    auz = coord_loc(3)*noz

    kx = int(aux)
    ky = int(auy)
    kz = int(auz)

    cell = 1+kx+nox*(ky+noy*kz)

    return
  end subroutine coord_to_cell
end module coord_to_cell_mod
