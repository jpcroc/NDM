module defcdp
  use T_kind_param_m


  integer :: &
       itecdp, &     ! introduction de DP tout les itecdp pas
       nintrodp, &   ! nombre de DP intrduit à chaque fois
       imin, &       ! indice minimal possible pour les atomes déplacés
       imax, &       ! indice maximal possible pour les atomes déplacés
       nposI,&        ! nombre de positions interstitielles
       iseed, &      ! racine des nombres aléatoires
       typint        ! type d'introduction des Intestitiels : 0 dans les sites prédéfinis, 1 aléatoirement

  real(double), dimension(:,:), pointer :: xposint ! positions des interstitiels
  real(double), dimension(20) :: Ed ! Energies de seuils
  real (double) :: dminins
  integer:: ioxdef
end module defcdp
