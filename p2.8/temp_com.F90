module temp_com
  USE T_kind_param_m
  implicit none


  integer,allocatable::na(:),nas(:),nad(:)
  integer :: imm,im
!  integer :: im_glob,imm_glob ! taille complète des tableaux 
  integer:: imm_loc ! taille des conf atomique par proc 

  real(double), dimension(3,3) :: at, bg ! at : vecteurs de base de la boite (BOND en cm) bg: vecteur du reseau reciproque
  real(double) :: volu      ! volume
  real(double), dimension(3) :: zl, zls2,nzl    ! largeur de la boite et largeur sur 2
  real(double), dimension(3,3) :: h0     ! Vecteurs de base de la boite de reference en A (Parrinello, Rahman)
  real(double), dimension(3) :: normat ! norme de at
  real(double), dimension(:,:,:),allocatable :: sigc ! contrainte par cel
  real(double), dimension(:,:,:),allocatable :: sigat ! contrainte par atome
  real(double), dimension (:),allocatable ::tempc,tempcm,celpm1,tm1,celpp,tcp,pmc,patcel,patcelmax,pm1

  logical :: ltabvois
  
  integer :: nox, noy, noz, noxy, noxyz	     !nb de cel suivant x y z et total (DOIT REMPLACER nce)
  integer :: natperc
  real(double), dimension(3) :: celsize	     ! taille des cel
  integer, dimension(:,:), allocatable :: ncel  ! ncel(i,j) indice de la jeme cel voisines de la cel i
  integer, dimension(:), allocatable :: nato    ! nb d'atome dans la ieme cel
  integer, dimension(:,:),allocatable :: atincel	! last (i,j) numero du ieme atome de la jeme cel
  integer, dimension(:,:,:),allocatable :: deltadist ! decalage a appliquer sur la cel

  integer :: cell_debx, cell_deby, cell_debz     !numero de la premiere cellule locale suivant x, y et z
  integer :: cell_finx, cell_finy, cell_finz     !numero de la derniere cellule locale  suivant x, y et z
  integer :: nb_cell_x, nb_cell_y, nb_cell_z     !nb de cel locales suivant x y z


  real(double),allocatable::eatom(:) ! energie par atome
  real(double),allocatable::eatomtotm(:) ! energie par atome
  integer, allocatable,dimension(:) :: indi ! table des voisins
  integer, allocatable, dimension(:) :: indi2 ! table de voision pour les constantes de force
  integer :: nvois   ! nb de voisins max dans toute la boite = nb d'atome * nb de voisins (/2)
end module temp_com
