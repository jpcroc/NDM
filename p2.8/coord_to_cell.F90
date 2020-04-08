module coord_to_cell_mod
        USE cryst_to_cart_mod
        implicit none
        contains
subroutine coord_to_cell(tab_coord, cell)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:bg,nox,noy,noz

  implicit none

  !       version du 10 janvier 2007
  !
  ! Cette routine retourne dans cell le numero de cellule
  ! contenant les coordonnees tab_coord
  !
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(double)  :: tab_coord(3)
  integer       :: cell
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: kx, ky, kz
  real(double) :: aux, auy, auz
  real(double) :: coord_loc(3)   ! permet de ne pas ecraser coord
  ! lors de l'appel a cryst_to_cart                       
  !-----------------------------------------------

  coord_loc(:) = tab_coord(:)
  call cryst_to_cart (1, coord_loc, bg, -1) !cart vers cryst

  aux = coord_loc(1)*nox
  auy = coord_loc(2)*noy
  auz = coord_loc(3)*noz

  kx = int(aux)
  ky = int(auy)
  kz = int(auz)

  cell = 1+kx+nox*(ky+noy*kz)

  return
end subroutine coord_to_cell
end module
