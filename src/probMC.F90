module probMC
  USE T_kind_param_m, ONLY:  double
  USE cryst_to_cart_mod, ONLY: cryst_to_cart
  implicit none
contains
!-----------------------------------------------------------------------
  subroutine probMC1(coord,proba,bg,at,fdmc_1,fdmc_2)
    real(double), dimension(3,1),intent(in) :: coord
    real(double), intent(in) :: fdmc_1,fdmc_2,bg(3,3),at(3,3)
    real(double), intent(out) :: proba

     
    real(double) :: dist_tot, alpha
    dist_tot = 0.0
    alpha = 0.0
    call calcul_dist(coord, dist_tot,bg,at) !calcul des distances
    !write(*,*) i, dist_tot
    call fct_alpha(dist_tot, alpha,fdmc_1,fdmc_2) !passage dans la fct alpha
    proba = dist_tot * alpha
  end subroutine probMC1

  subroutine calcul_dist(coord_atom, dist,bg,at)
    implicit none
    real(double), intent(in) :: bg(3,3),at(3,3)
    real(double) :: dist, m, n
    real(double), dimension(3,1) :: coord_atom, ref, distance
    integer :: j, k, l
    integer, dimension(8) :: liste_entier
    !attention, routine exacte uniquement pour les boites 4x4x4 de UO2

    liste_entier = (/3,7,11,15,19,23,27,31/)
    !write(*,*) 'liste_entier' , liste_entier(1)
    !write(*,*) 'coord atom avant', coord_atom(:,1)
    call cryst_to_cart(1,coord_atom,bg,-1)
    !write(*,*) 'coord atom apres', coord_atom(:,1)
    DO j=1,3 !calculer la distance au site interstitiel 'parfait' le plus proche
       m = 2
       do k=1, 8
          l = liste_entier(k) 
          n = abs(l*0.03125 - coord_atom(j,1))
          if (n .lt. m) then
             ref(j,1) = l*0.03125
             m = n
          end if
       end do
       !write(*,*) ref(j,1), coord_atom(j,1)
       distance(j,1) = ref(j,1) - coord_atom(j,1)
       if (distance(j,1) .gt. 0.5) then
          distance(j,1) = distance(j,1) -1
       end if !CP si >0.5
       if (distance(j,1) .lt. -0.5) then
          distance(j,1) = distance(j,1) +1
       end if !CP si <-0.5
       dist = dist + (distance(j,1)*at(j,j))**2
    END DO !boucle sur les coord
    dist = dsqrt(dist)*1E8

  end subroutine calcul_dist

  subroutine fct_alpha(dist, dist_alpha,fdmc_1,fdmc_2)
    implicit none
    real(double), intent(in) :: fdmc_1,fdmc_2
    real(double) :: dist, dist_alpha

    dist_alpha = 1.0-(1.0/( (dexp( (dist-fdmc_1) / fdmc_2)) +1.0) )
  end subroutine fct_alpha

end module probMC
