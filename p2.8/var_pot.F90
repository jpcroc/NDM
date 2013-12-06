module var_pot
  use T_kind_param_m
  implicit none

  integer,parameter::npotmax=100
  integer, parameter :: nkmax = 4000
  integer, parameter :: contmax = 2000


  logical ::lprtpot

  integer :: ntyp, ntyp_buffer					! nb de type
  integer :: npair						! = ntyp*(ntyp+1)/2
  integer :: ntrip						! = ntyp*ntyp *(ntyp+1)/2


  integer, dimension(:), pointer  :: na			! nb d'atomes par type
  integer, dimension(:,:), pointer  :: ipo			! indice des paires d'atomes
  real(double), dimension(:), pointer :: cm, cm_buffer, catom, q, rc	! masse, numero atomique, charge ionique, rayon de coup.
  character , dimension(:), pointer  :: ty*3
  character , dimension(:), pointer  :: ty_buffer*3

  real(double)::rclu(20), eatref(20)   ! rayon et energie des types d'atomes

  integer :: ipotentiel,npotentiel    ! type du potentiel COURANT 1=BMH, 2=buckingham 3=watanabe,4=UO2; etc...
  logical :: lpotentiel (0:npotmax)
  logical, pointer :: typ_and_pot(:,:) ! typ_and_pot(iti,ipot)=.true. si le type iti interagit (en autres) par le potentiel ipot
  !  integer::lupotin=95
  integer, pointer:: typ_pot_pair(:) ! donne le type d'interaction de la paire
  logical, pointer::lu_roff_pair(:)
  logical,pointer::lue_typ(:),lue_trip(:)
  integer :: ipotrep ! repulsion courte distance 0 = rien, 1 =polynom ; 2 =ziegler
  logical, pointer, dimension (:) :: lue_paire

  integer,pointer:: ipo_2_pair_tab(:)
  real(double),pointer::pot_pair_tab(:,:,:)
  integer::ngr


  logical :: l3c ! somme d'Ewald terme a trois corps
  integer :: iewald

  integer:: ngrid  ! taille de la grille des pot de paire : lue dans .din


  real(double), dimension(:,:,:), pointer :: pot ! table des pot splines
  real(double), dimension(:,:,:), pointer :: pot_d ! table des pot splines
  real(double), dimension(:,:), pointer :: potw ! table des pot a spliner
  real(double), dimension(:), pointer :: ray, shel, bm !parametres du pot
  real(double), dimension(:), pointer :: Dmorse, amorse, remorse !parametres du pot Morse
  real(double), dimension(:), pointer :: zz ! qi*qj
  real(double) :: eta ! rayon de coupure et amortissement d'Ewald
  real(double) ::  rumax,csive ! rayonde coupure ; pas de la grille d'interpolation du potentiel
  real(double),pointer::rue_pair(:)
  real(double),pointer::rue_pot(:)
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


  real(double),pointer, dimension (:) :: gz,fcr! tab. des valeurs des fonctions
  real(double), dimension (:), pointer ::            Awat,    &! tab. parametre de Stilliger-Weber a 2 corps
       Bwat,    &! tab. parametre de Stilliger-Weber a 2 corps
       pwat,   & ! tab. exposants de SW a 2 corps
       qwat,   & ! tab. exposants de SW a 2 corps
       rawat,  & ! tab. rayons de coupure des inter. a 2 corps
       rawat2    ! carrÃ£Â© des rayons de coupure



  !3 corps a la sauce JDT
  real(double), dimension (:), pointer ::lamb, cangle,C3C
  real(double), dimension (:,:), pointer :: gam,coup3c,coup3c2
  integer, dimension(:,:,:), pointer :: ipo3c
  logical, dimension(:), pointer :: l3ctyp
  logical, dimension(:), pointer :: l3cpair
  real(double) :: r3cm,r3cm2

  !4 Potentiel UO2
  real(double), dimension (6) :: poly5
  real(double), dimension (4) :: poly3
  real(double) :: rbp5, rp5p3, rp3c

  real(double),pointer, dimension (:) :: bspg,cspg,dspg,bspf,cspf,dspf ! spline de watanabe

  real(double),parameter  :: evA62ergcm6=1.6021892D-60   ! conversion eV.A^6 --> erg.cm^6
  real(double), dimension(:), pointer :: ro, dip, pm, roff1, roff2, a_factor,r8p ! potentiel
  real(double), dimension(:,:), pointer :: bspw, cspw, dspw ! spline
  real(double) :: alpha


  !6 Stillinger Weber Vashista JAP 101, 103515 (07)
  real(double):: lambda,xsi
  real(double),pointer::capHij(:),capDij(:),capWij(:)
  integer,pointer :: ietaij(:)


  ! EAM
  logical :: lforcetabulate ! if the first derivative is tabulate.
  real(double) :: potisrep, potisglue,potiseam ! energie potentielle EAM
  real(double),dimension(:,:,:),pointer :: eamrep,eamrho,eamglue       ! tableaux des splines du pot EAM 
  real(double),dimension(:,:,:),pointer :: eamrep_d,eamrho_d,eamglue_d ! tableaux des splines du pot_d EAM 
  real(double) :: rhomin=1d30,rhomax=0

  real(double), pointer, dimension(:,:,:) :: digr, coord
  integer, dimension(:), pointer :: nad, nas, nai  ! fracture

  real(double),pointer, dimension(:,:,:,:) :: fda

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
  real(double), dimension(:),pointer :: bsmod1   !bspline
  real(double), dimension(:),pointer :: bsmod2
  real(double), dimension(:),pointer :: bsmod3
  real(double), dimension(:,:),pointer :: table  !pour fftfront
  integer, dimension(:,:),pointer :: iiim,ijim,ikim  !calcul de qgrid
  real(double), dimension(:),pointer :: fr1,fr2,fr3  !calcul de qgrid
  real(double), dimension(:),pointer :: de1,de2,de3  !calcul de fp
  !jm       real(double), dimension(:),pointer :: w1pme,w2pme,w3pme




end module var_pot
