module gen_com_m
  use T_kind_param_m
  implicit none


#ifdef para

  ! Declarations MPI
  !      include 'mpif.h'
  ! option pour realiser une trace Vampir
  !     include 'VT.inc'
  !      integer rang,code,nb_procs
  !      common/parampi/rang,code,nb_procs

#else

  !      integer :: rang, tranche
#endif



  integer :: rang
  logical :: parallele

  integer :: natperc                        ! nb d'atome par cel
  integer :: nvperat    ! Nombre moyen de voisins par atome

  integer :: imm, num_paire                 !imm taille des tableaux dependant du nombre d'atome
  integer :: imm_glob


  integer,parameter::npotmax=100
  real(double),parameter :: pi=3.141592654D0, bk= 1.380622D-16, &
       ecgs=1.6021892D-12, utemps= 1.0D-15, angst= 1.0D08, umass= 1.660056D-24, &
       inv_angst=1.d0/angst
  real(double),parameter :: zero=0d0, one=1.0d0, two=2.0d0, thr=3.0d0, five=5.0d0&
       , six=6.0d0, half=0.5d0
  real(double), parameter :: precexp =0.004   ! induit une precision de exp 10^-100
  real(double), parameter :: ev2erg=1.602d-12, erg2eV=1.d0/eV2erg   
  real(double), parameter :: erg2joule=1.d-7, joule2erg=1.d7        
  real(double), parameter :: low_limit=10.d0*epsilon(1.d0)
  real(double), parameter :: ang2cm=1d-8
  real(double), parameter :: e2on4pieps0= 23.06134575D-20


  integer :: im						! nb local d'atomes (=global en sequentiel)
  integer :: im_glob					! nb global d'atomes
  integer :: ntyp						! nb de type
  integer :: npair						! = ntyp*(ntyp+1)/2
  integer :: ntrip						! = ntyp*ntyp *(ntyp+1)/2
  integer, dimension(:), pointer  :: na			! nb d'atomes par type
  integer, dimension(:,:), pointer  :: ipo			! indice des paires d'atomes
  real(double), dimension(:), pointer :: cm, catom, q, rc	! masse, numero atomique, charge ionique, rayon de coup.
  real(double)::rclu(20), eatref(20)   ! rayon et energie des types d'atomes
  character , dimension(:), pointer  :: ty*3

  real(double), dimension(3) :: zl, zls2,nzl    ! largeur de la boite et largeur sur 2
  real(double) :: volu      ! volume
  integer, dimension(3) :: lat      ! generation: nb de repetition de cel unite
  real(double), dimension(3,3) :: at, bg ! at : vecteurs de base de la boite (BOND en cm) bg: vecteur du reseau reciproque
  real(double), dimension(3,3) :: h0     ! Vecteurs de base de la boite de reference en A (Parrinello, Rahman)
  logical :: lUcell                 ! affiche l'energie potentielle de la boite
  ! (cela suppose que h0 corresponde a l'etat de reference pour lequelle la contrainte est nulle)     

  logical :: lfrozen    ! .true.: certains atomes sont bloques (pas de dynamique)
  logical :: lbulle    ! .true.: bulle
  logical :: ldesinteg    ! .true.: insertion appelle init_insert
  integer:: nstepdes,ides, pm1des,itdes,imdesup,imdesdeb,imdesdn
  real(double)::lambdades,deltaF
  
  real(double),pointer:: xpchup(:,:),vpchup(:,:),xpchdn(:,:),vpchdn(:,:),xpchdeb(:,:),vpchdeb(:,:)
  integer,pointer :: itichup(:),itichdn(:),itichdeb(:),num_at_globdesup(:),num_at_globdesdn(:),num_at_globdesdeb(:)
  real(double):: kspr,xpspr(3),Espr,deltaEspr,xpspr0(3),tempdes
  integer:: typspr


  logical, dimension(:), pointer :: free ! free(i)=.true. si l'atome i compte dans l'energie 
  logical, dimension(:,:), pointer :: frozen ! frozen(ix,i)=.true. si la coordonnee ix de l'atome i est libre de relaxer
  integer::imfree ! nb d'atoems libres

  real(double), dimension(3) :: normat ! norme de at
  integer :: ipotentiel,npotentiel    ! type du potentiel COURANT 1=BMH, 2=buckingham 3=watanabe,4=UO2; etc...
  logical :: lpotentiel (0:npotmax)
  logical, pointer :: typ_and_pot(:,:) ! typ_and_pot(iti,ipot)=.true. si le type iti interagit (en autres) par le potentiel ipot
  !  integer::lupotin=95
  integer, pointer:: typ_pot_pair(:) ! donne le type d'interaction de la paire
  logical, pointer::lu_roff_pair(:)
  logical,pointer::lue_typ(:),lue_trip(:)
  logical :: lpotrep ! repulsion courte distance
  logical, pointer, dimension (:) :: lue_paire

  integer,pointer:: ipo_2_pair_tab(:)
  real(double),pointer::pot_pair_tab(:,:,:)
  integer::ngr

  integer, dimension(:,:), pointer :: ncel  ! ncel(i,j) indice de la jeme cel voisines de la cel i

  integer, dimension(:), pointer :: nato    ! nb d'atome dans la ieme cel

  integer, dimension(:,:),pointer :: last	! last (i,j) numero du ieme atome de la jeme cel

  integer, dimension(:,:,:),pointer :: deltadist ! decalage a appliquer sur la cel
  integer :: nox, noy, noz, noxy, noxyz	     !nb de cel suivant x y z et total (DOIT REMPLACER nce)
  integer :: cell_debx, cell_deby, cell_debz     !numero de la premiere cellule locale suivant x, y et z
  integer :: cell_finx, cell_finy, cell_finz     !numero de la derniere cellule locale  suivant x, y et z
  integer :: nb_cell_x, nb_cell_y, nb_cell_z     !nb de cel locales suivant x y z
  real(double), dimension(3) :: celsize	     ! taille des cel


  integer :: it, itmax, igen ! iteration courante, finale , type de generation
  integer :: lenfnam
  integer :: fmt_cin
  integer :: iewald
  character :: fnam*80, fnamout*80, fnamcout*80, fnamcoutxp*80,fnamcoutnonpbcxp*80
  logical :: ltranche ! surface
  integer:: iteplz,nplz ! distribution suivant des tranches en z
  real(double):: rulayer
  integer, dimension(:), pointer :: nad, nas, nai  ! fracture
  integer :: imgs, imgi, imd, itefrac !fracture IMD nombre d'atomes sur lesquels on fait la dynamique normale
  real(double) :: cougel, zincr !fracture

  logical :: l3c ! somme d'Ewald terme a trois corps


  integer:: ngrid  ! taille de la grille des pot de paire : lue dans .din


  real(double), dimension(:,:,:), pointer :: pot ! table des pot splines
  real(double), dimension(:,:), pointer :: potw ! table des pot a spliner
  real(double), dimension(:), pointer :: ray, shel, bm !parametres du pot
  real(double), dimension(:), pointer :: Dmorse, amorse, remorse !parametres du pot Morse
  real(double), dimension(:), pointer :: zz ! qi*qj
  real(double) :: eta ! rayon de coupure et amortissement d'Ewald
  real(double) ::  rumax,csive ! rayonde coupure ; pas de la grille d'interpolation du potentiel
  real(double),pointer::rue_pair(:)
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


  real(double), dimension(:), pointer :: ro, dip, pm, roff1, roff2, a_factor,r8p ! potentiel
  real(double), dimension(:,:), pointer :: bspw, cspw, dspw ! spline
  real(double) :: alpha
  real(double) :: potist ! energie potentielle totale
  real(double):: potisP,potis1, potis2, potis3, potis0, potcp ! energie potentielle de paire
  real(double) :: potisTersoff ! energie potentielle de tersoff


  !6 Stillinger Weber Vashista JAP 101, 103515 (07)
  real(double):: lambda,xsi
  real(double),pointer::capHij(:),capDij(:),capWij(:)
  integer,pointer :: ietaij(:)


  ! EAM
  real(double) :: potisrep, potisglue,potiseam ! energie potentielle EAM
  real(double),dimension(:,:,:),pointer :: eamrep,eamrho,eamglue ! tableaux des splines du pot EAM 
  real(double) :: rhomin,rhomax

  real(double), dimension(:), pointer :: h2sm ! delta t carre sur 2 m
  real(double) :: tstep, oldtstep, usdh, timel  
  integer :: itetemp, itesigma, itedepla, itecoordo, iterdf, nrdf, & 
       iterasmol, iteangle,nfda,itetemp2,iteanapos, itefcc,itecfg
  real(double)::rcangle,rcrdf

  real(double), dimension(3,3) :: sig ! contrainte
  real(double), dimension(3,3) :: sigtot
  real(double), dimension(3,3) :: sigkine

  real(double), dimension(:,:,:),pointer :: sigc ! contrainte par cel
  real(double), dimension(:,:,:),pointer :: sigat,sigtyp,sigtyp_loc ! contrainte par atome
  real(double), dimension(:,:,:,:),pointer :: sigtyptyp,sigtyptyp_loc ! contrainte par atome
  logical :: lsigat,lEparat,lsigtyp  ! calcul et affichage dans rasmol de la contrainte atomique; affichage Ã£Â£Ã¢Â£Ã£Â¢Ã¢Â£Ã£Â£Ã¢Â¢Ã£Â¢Ã¢Â©nergie par atome,calcul bond valence
  integer:: itebdv ! frequence de calcul des bond valence
  logical :: ljqbh ! calcul de la conductivitÃ© thermique par la mÃ©thode directe
  logical :: lnemd  ! calcul de la conductivitÃ© thermique par NEMD
  integer ::njqbh,ittherm,ntr
  real(double) :: epsil,epcoud,kthg


  real(double)::fnemd
  real(double):: fpstop ! critere de conv. sur la force par atome max  pour les trempes UNITE = EV/ANG
  real(double):: sigstop ! critere de conv. sur les composantes de contraintes  pour les trempes UNITE = kbar
  real(double):: fsumstop ! critere de conv. sur la force sqrt ( sum_f F_i^2 )  pour les trempes UNITE = EV/ANG
  real(double),pointer::eatom(:) ! energie par atome
  real(double),pointer::eatomtotm(:) ! energie par atome
  logical :: lprteat, lprtfat,lprteattotm  ! calcul et ecriture de l'energie et force par atome, de l'energie par atome totale (pot+cin) moyenne
  logical :: lposmoy ! ecrit à la fin la position moyenne des atomes
  real(double) :: tdepla, tdepla2 ! seuils de deplacement
  logical :: lfilm, linstantrdf,linstantfda, lrestart, ltpcel, lfilmext !film, RDF, restart, moyenne par cel

  !Correlations et Cie
  logical :: lcorrelvp ! ecriture de l'autocorrelation des vitesses*Masses
  logical :: lcalcjq

  logical :: lcdp ! algorithme d'accumulation de defauts ponctuels
  logical :: lheat ! algorithme de chauffage local
  real(double)::rheat,theat,Eheat
  integer :: iteheat


  real(double) :: tinit !temp initiale
  logical :: lvpread  ! vitesse lue dans le fichier .cin
  integer:: iseed ! graine du gerateur aleatoire des vitesses
  integer :: dmtype, itab, itetabvois, itetimestep, itederive ! type dynamique, periode de repartition entre cel, periode de calc. tab des voisins, periode de chgt du pas en temps, poeriode de correction de la derive
  real(double) :: tempstop, ttol, tfroi, tcooling, tcou, tfcou, epcou, &! temperature d'arret, max, visee si max, taux de refroidissement, temp de la couche externe et epaisseur
       tsfact, vmax, tgc, dfpred ! gestion du pas en temps
  real(double) :: deltaestop ! decroissance de la temperature moyenne
  integer :: nbmoye
  integer :: ibordcou
  integer :: itesauv, formatsauv, itesauvposition ! periode de sauvegarde format de sauvegarde periode de d'ecriture des positions en formatted
  real(double), dimension(3) :: vh ! vitesse de la boite
  real(double) :: pext, wbox, tbox ! pext poids de la boite temps d'amortissment de la boite
  logical ::  lpcon2,lprtzlm ! pression constante sans et avec amortissement
  logical :: lTcon, lTberendsen,lTandersen,lTNose,lTHoover,landerscou ! temp constante (3 algorithmes differents)
  real(double) :: Text ! T exterieure
  logical :: lLangevin ! Langevin MD
  real(double) :: gamlang  ! gamma de Langevin
  integer :: ilangevin ! 1=std ; 2=Athenes
  real(double) :: tauTcon ! Temps berendsen
  real(double) :: nuandersen ! frequence andersen     
  integer :: iteTconst! periode de temp const
  integer :: iko ! indide du PAF
  real(double) :: eko, xko, yko, zko ! energie et direction du PAF
  real(double) :: xx0, yy0, zz0 ! position initiale du projectile
  logical :: lcasca,lderive ! cascade,correction derive ?
  real(double) :: pist, temp, pmean, tmean, kine, kinemean ! pression temp et moyennes associees


  integer :: nvois   ! nb de voisins max dans toute la boite = nb d'atome * nb de voisins (/2)
  integer, pointer,dimension(:) :: indi ! table des voisins
  real(double) :: rvois ! rayon de la table des voisins
  logical :: ltabvois ! table des voisins ?
  logical :: lconstrtot ! construction par double boucle (T) ou par cel (F)
  logical :: ldemitab ! construction d'une demi-table (T) ou d'une table complete (F)
  character :: nature*6 ! element chimique
  integer :: nvat
  !EWALD
  real(double), dimension(:,:,:),pointer :: tabv3
  real(double), dimension(:,:,:,:),pointer :: tabf3

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

  logical :: lalea  ! preparation d'une configuration aleatoire
  logical :: lopt   ! optimisation de Ewald par PME si TRUE
  real(double) :: rsep  !Distance de separation pour le tirage aleatoire

  ! NVT, NPT ensembles
  !      logical :: lnose, lnosepar,lpr ! lnose =Tcst  la Nose ; lnosepar=T&P cst a la Nose Parinello Rahman
  logical :: lpr,lprtrp ! l Parinello Rahman
  !      real(double) :: tomega, tbomega ! mass fictive du thermostat et du piston
  real(double), dimension(3,3) :: att, ati    !vitesse de la forme de la boite ; ati=(at^-1)
  real(double), dimension(3,3) :: sigext, pext_hydro  !contraintes externes appliques; contraintes calculees

  integer :: lurdfout
  real(double) :: lastcool


  real(double), parameter :: rmin = 0.5d-8
  integer, parameter :: kmax = 3000
  integer, parameter :: nkmax = 4000
  real(double), parameter :: qmax=12
  real(double), parameter :: qmin=0.7
  real(double), parameter :: increq=0.1d+8
  !      real(double), dimension(ntyp,ntyp,nkmax) :: digr, coord,strucfact
  real(double), pointer, dimension(:,:,:) :: digr, coord,strucfact
  real(double), dimension(nkmax) :: gdertot,strucfactot,strucfactneu
  integer, dimension(:,:), pointer :: voisins
  real(double), parameter :: rcut=12e-8    !cutoff pour le calcul de S(q)

  integer, parameter :: cont = 1000
  integer, parameter :: contmax = 2000
  real(double),pointer, dimension(:,:,:,:) :: fda
  real(double), parameter :: thetamin = 1.0D-7
  real(double), parameter :: thetamax = 6.2

  logical lEev,lPkbar   !unite
  ! energies potentielle, cinetique et totale de la boite en Parrinello-Rahman
  real(double):: Ecell, Kcell, Ucell      


  !Variables Nose
  REAL(double) :: ENose, KNose, UNose
  REAL(double) :: fNose
  real(double) :: wNose     ! Poids associe au thermostat de Nose
  INTEGER :: nHoover        ! Nombre de chaines de Hoover
  REAL(double), dimension(:), allocatable :: zHoover   ! Viscosite 

  real(double) :: deltax 
  real(double) :: dilat(3)



  integer :: imf     ! forces normales depuis 1 jusqu'Ã£Â£Ã¢Â£Ã£Â¢Ã¢Â£Ã£Â£Ã¢Â¢Ã£Â¢Ã¢Â£Ã£Â£Ã¢Â£Ã£Â¢Ã¢Â¢Ã£Â£Ã¢Â¢Ã£Â¢Ã¢Â  imf
  integer :: imana     ! configurations analysÃ£Â£Ã¢Â£Ã£Â¢Ã¢Â£Ã£Â£Ã¢Â¢Ã£Â¢Ã¢Â£Ã£Â£Ã¢Â£Ã£Â¢Ã¢Â¢Ã£Â£Ã¢Â¢Ã£Â¢Ã¢Â©es d

  logical :: ldislo  ! calcul de dislocation
  real(double) :: epcoudis,& !epaisseur de la couche avec ajout de force pour dislo
       &fdislo ! force appliquÃ£Â£Ã¢Â£Ã£Â¢Ã¢Â£Ã£Â£Ã¢Â¢Ã£Â¢Ã¢Â£Ã£Â£Ã¢Â£Ã£Â¢Ã¢Â¢Ã£Â£Ã¢Â¢Ã£Â¢Ã¢Â©e aux atomes de bords 
  integer, pointer :: latdebord(:)

  logical :: lcontr    ! dynamique contrainte (routine contrainte)
  logical :: lperiod   ! conditions periodiques
  logical :: lsuivinonpbc


  !---inNEB
  integer  :: ipath, npath,nebtype,nebrelaxation,maxneb,iteanaposneb, &
              neb_noise
  REAL(double) :: kspring,deltaRmax,neb_noise_scale
  !...inNEB


  ! chauffage cylindre
  logical :: lHcyl ! variable de type logique representant le chauffage du cylindre
  real(double) :: Ecyl ! energie totale des atomes dans le cylindre
  real(double), dimension (3) :: pc ! position du centre du cylindre
  real(double), dimension (3) :: vdc ! vecteur direction du cylindre
  real(double) :: rayonc ! rayon du cylindre
  real(double) :: lgc ! longueur du cylindre
  integer :: ncyl ! nombre d'atomes dans le cylindre
  logical, dimension(:), pointer :: cyl ! tableau pour savoir si atome dans cylindre

  ! selection des atomes distordus
  integer :: natdistordusvraiment

end module gen_com_m
