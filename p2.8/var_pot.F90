module var_pot
  USE T_kind_param_m
  implicit none

  integer,parameter::npotmax=100
  integer, parameter :: nkmax = 10000
  integer, parameter :: contmax = 2000


  logical ::lprtpot

  integer :: ntyp, ntyp_buffer					! nb de type
  integer :: npair						! = ntyp*(ntyp+1)/2
  integer :: ntrip						! = ntyp*ntyp *(ntyp+1)/2


  integer, dimension(:), allocatable  :: na			! nb d'atomes par type
  integer, dimension(:,:), allocatable  :: ipo			! indice des paires d'atomes
  real(double), dimension(:), allocatable :: cm, cm_buffer, catom, q, rc	! masse, numero atomique, charge ionique, rayon de coup.
  character , dimension(:), allocatable  :: ty*3
  character , dimension(:), allocatable  :: ty_buffer*3
  real(double),dimension(:),allocatable::gamlt(:)

  real(double)::rclu(20), eatref(20)   ! rayon et energie des types d'atomes

  integer :: ipotentiel,npotentiel    ! type du potentiel COURANT 1=BMH, 2=buckingham 3=watanabe,4=UO2; etc...
  logical :: lpotentiel (0:npotmax)
  logical, allocatable :: typ_and_pot(:,:) ! typ_and_pot(iti,ipot)=.true. si le type iti interagit (en autres) par le potentiel ipot
  !  integer::lupotin=95
  integer, allocatable:: typ_pot_pair(:) ! donne le type d'interaction de la paire
  logical, allocatable::lu_roff_pair(:)
  logical,allocatable::lue_typ(:),lue_trip(:)
  integer :: ipotrep ! repulsion courte distance 0 = rien, 1 =polynom ; 2 =ziegler
  logical, allocatable, dimension (:) :: lue_paire

  integer,allocatable:: ipo_2_pair_tab(:)
  real(double),allocatable::pot_pair_tab(:,:,:)
  integer::ngr


  logical :: l3c ! somme d'Ewald terme a trois corps
  integer :: iewald

  integer:: ngrid  ! taille de la grille des pot de paire : lue dans .din


  real(double), dimension(:,:,:), allocatable :: pot ! table des pot splines
  real(double), dimension(:,:,:), allocatable :: pot_d ! table des pot splines
  real(double), dimension(:,:), allocatable :: potw ! table des pot a spliner
  real(double), dimension(:), allocatable :: ray, shel, bm !parametres du pot
  real(double), dimension(:), allocatable :: Dmorse, amorse, remorse !parametres du pot Morse
  real(double), dimension(:), allocatable :: zz ! qi*qj
  real(double) :: eta ! rayon de coupure et amortissement d'Ewald
  real(double) ::  rumax,csive ! rayonde coupure ; pas de la grille d'interpolation du potentiel
  real(double),allocatable::rue_pair(:)
  real(double),allocatable::rue_pot(:)
  integer :: ncouc3 ! nombre de couche dans la sommation d'Ewald
  integer :: n2max  ! valeur de ncouc3 au carre
  integer :: ncoucx, ncoucy, ncoucz,nvecttot !couches en x y et z de la sommation d'Ewald
  real(double) :: precis ! precision du calcul de la sommation d'Ewald

  ! 2 corps watanabe
  real (kind=double) :: & ! 2 corps
       epswat, & ! unite reduite d'energie
       sigmawat,&  ! unite reduite de longueur
       gm1,gm2,gm3,gm4,gm5,gR,gD,&       ! parametres de fonct g watamabe
       
       csive_g  ! taille grille pour discretiser la fonction g


  real(double),allocatable, dimension (:) :: gz,fcr! tab. des valeurs des fonctions
  real(double), dimension (:), allocatable ::            Awat,    &! tab. parametre de Stilliger-Weber a 2 corps
       Bwat,    &! tab. parametre de Stilliger-Weber a 2 corps
       pwat,   & ! tab. exposants de SW a 2 corps
       qwat,   & ! tab. exposants de SW a 2 corps
       rawat,  & ! tab. rayons de coupure des inter. a 2 corps
       rawat2    ! carrÃ£Â© des rayons de coupure



  !3 corps a la sauce JDT
  real(double), dimension (:), allocatable ::lamb, cangle,C3C
  real(double), dimension (:,:), allocatable :: gam,coup3c,coup3c2
  integer, dimension(:,:,:), allocatable :: ipo3c
  logical, dimension(:), allocatable :: l3ctyp
  logical, dimension(:), allocatable :: l3cpair
  real(double) :: r3cm,r3cm2

  !4 Potentiel UO2
  real(double), dimension (6) :: poly5
  real(double), dimension (4) :: poly3
  real(double) :: rbp5, rp5p3, rp3c

  real(double),allocatable, dimension (:) :: bspg,cspg,dspg,bspf,cspf,dspf ! spline de watanabe

  real(double),parameter  :: evA62ergcm6=1.6021892D-60   ! conversion eV.A^6 --> erg.cm^6
  real(double), dimension(:), allocatable :: ro, dip, pm, roff1, roff2, a_factor,r8p ! potentiel
  real(double), dimension(:,:), allocatable :: bspw, cspw, dspw ! spline
  real(double) :: alpha


  !6 Stillinger Weber Vashista JAP 101, 103515 (07)
  real(double):: lambda,xsi
  real(double),allocatable::capHij(:),capDij(:),capWij(:)
  integer,allocatable :: ietaij(:)


  ! EAM
  logical :: lforcetabulate ! if the first derivative is tabulate.
  real(double) :: potisrep, potisglue,potiseam ! energie potentielle EAM
  real(double),dimension(:,:,:),allocatable :: eamrep,eamrho,eamglue       ! tableaux des splines du pot EAM 
  real(double),dimension(:,:,:),allocatable :: eamrep_d,eamrho_d,eamglue_d ! tableaux des splines du pot_d EAM 
  real(double) :: rhomin=1d30,rhomax=0

  real(double), allocatable, dimension(:,:,:) :: digr, coord
  integer, dimension(:), allocatable :: nad, nas, nai  ! fracture

  real(double),allocatable, dimension(:,:,:,:) :: fda

  real(double), dimension(nkmax) :: gdertot,strucfactot,strucfactneu

  !PME
  integer :: kpmex, kpmey, kpmez   !taille de grille de PME
  integer :: kpme                  ! max des precedants
  integer :: maxorder, iorder  !ordre de la PME (bspline)
  integer :: npoint != kpmex*kpmey*kpmez
  integer :: nfft1, nfft2, nfft3  ! ~kpmex
  integer :: nff, nf1, nf2, nf3
  integer :: ntable      ! pour fftfront
  real(double) :: pterm, volterm, auxe ! constantes pour PME
  real(double), dimension(:),allocatable :: bsmod1   !bspline
  real(double), dimension(:),allocatable :: bsmod2
  real(double), dimension(:),allocatable :: bsmod3
  real(double), dimension(:,:),allocatable :: table  !pour fftfront
  integer, dimension(:,:),allocatable :: iiim,ijim,ikim  !calcul de qgrid
  real(double), dimension(:),allocatable :: fr1,fr2,fr3  !calcul de qgrid
  real(double), dimension(:),allocatable :: de1,de2,de3  !calcul de fp
  !jm       real(double), dimension(:),allocatable :: w1pme,w2pme,w3pme




end module var_pot
