integer function calcul_cellule(coordx,coordy,coordz)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m

  real(double), dimension(3,1) :: coord_tab
  real(double) :: aux, auy, auz
  integer      :: kx,ky,kz


     coord_tab(1,1)=coordx
     coord_tab(2,1)=coordy
     coord_tab(3,1)=coordz

     call cryst_to_cart (1, coord_tab, bg, -1) !cart vers cryst

     aux = coord_tab(1,1)*nox
     auy = coord_tab(2,1)*noy
     auz = coord_tab(3,1)*noz
     kx = int(aux)
     ky = int(auy)
     kz = int(auz)
     calcul_cellule = 1+kx+nox*(ky+noy*kz)




  return
end function calcul_cellule
