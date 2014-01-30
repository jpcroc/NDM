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



  integer :: rang, rangph
  logical :: parallele

  integer :: natperc                        ! nb d'atome par cel
  integer :: nvperat    ! Nombre moyen de voisins par atome

  integer :: imm                 !imm taille des tableaux dependant du nombre d'atome
  integer :: imm_glob


  real(double),parameter :: pi=3.141592654D0, bk= 1.380622D-16, &
       ecgs=1.6021764631580d-12, &  !debugCOS 1.6021892D-12, &
       utemps= 1.0D-15, angst= 1.0D08, umass= 1.660056D-24, &
       inv_angst=1.d0/angst 
       
  real(double),parameter :: zero=0d0, one=1.0d0, two=2.0d0, thr=3.0d0, five=5.0d0&
       , six=6.0d0, half=0.5d0
  real(double), parameter :: precexp =0.004   ! induit une precision de exp 10^-100
  real(double), parameter :: ev2erg=1.6021764631580d-12, erg2eV=1.d0/eV2erg   !eV -> erg and in inverse
  real(double), parameter :: evA2dyn=1.6021764631580d-4  ! eV/A -> dyn conversion
  real(double), parameter :: hbar=1.05457266d-27  ! hbar


  real(double), parameter :: erg2joule=1.d-7, joule2erg=1.d7        
  real(double), parameter :: low_limit=10.d0*epsilon(1.d0)/10000000000000000000000000000000000000000000.d0
  real(double), parameter :: e2on4pieps0= 23.06134575D-20
  real(double) :: A2cm =1.0d-8     !conversion A->cm


  integer :: im						! nb local d'atomes (=global en sequentiel)
  integer :: im_glob					! nb global d'atomes


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


  logical, dimension(:), pointer :: Free ! free(i)=.true. si l'atome i compte dans l'energie 
  logical, dimension(:,:), pointer :: Frozen ! Frozen(ix,i)=.true. si la coordonnee ix de l'atome i est gelee
  integer::imFree,imFirstFrozen ! nb d'atoems libres

  real(double), dimension(3) :: normat ! norme de at



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

  character :: fnam*80, fnamout*80, fnamcout*80, fnamcoutxp*80,fnamcoutnonpbcxp*80
  logical :: ltranche ! surface
  integer:: iteplz,nplz ! distribution suivant des tranches en z
  real(double):: rulayer

  integer :: imgs, imgi, imd, itefrac !fracture IMD nombre d'atomes sur lesquels on fait la dynamique normale
  real(double) :: cougel, zincr !fracture



  real(double) :: potist ! energie potentielle totale
  real(double):: potisP,potis1, potis2, potis3, potis0, potcp ! energie potentielle de paire
  real(double) :: potisTersoff ! energie potentielle de tersoff



  real(double) :: oldtstep  
  real(double) :: tstep, usdh, timel  
  integer :: itetemp, itesigma, itedepla, itecoordo, iterdf, nrdf, & 
       iterasmol, iteangle,nfda,itetemp2,iteanapos, itefcc,itecfg
  integer::ivisu     ! format de sortie dans rasmol.f90 : ivisu=1=.mol, ivisu=2=vsim mal codé, ivisu=2=xred
  real(double)::rcangle,rcrdf

  real(double), dimension(3,3) :: sig ! contrainte
  real(double), dimension(3,3) :: sigtot
  real(double), dimension(3,3) :: sigkine

  real(double), dimension(:,:,:),pointer :: sigc ! contrainte par cel
  real(double), dimension(:,:,:),pointer :: sigat,sigtyp,sigtyp_loc ! contrainte par atome
  real(double), dimension(:,:,:,:),pointer :: sigtyptyp,sigtyptyp_loc ! contrainte par atome
  logical :: lsigat,lEparat,lsigtyp  ! calcul et affichage dans rasmol de la contrainte atomique; affichage ©nergie par atome,calcul bond valence
  integer:: itebdv ! frequence de calcul des bond valence
  logical :: ljqbh ! calcul de la conductivitÃ© thermique par la mÃ©thode directe
  logical :: lnemd  ! calcul de la conductivitÃ© thermique par NEMD
  integer ::njqbh,ittherm,ntr
  real(double) :: epsil,epcoud,kthg


  logical :: lFire      ! If true (default), fire algorithm is used for quenching
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
  logical::ldeplainit
  real(double)::tempdeplainit,debyetemp
  logical :: lvpread  ! vitesse lue dans le fichier .cin
  integer:: iseed ! graine du gerateur aleatoire des vitesses
  integer :: dmtype, itab, itetabvois, itetimestep, itederive ! type dynamique, periode de repartition entre cel, periode de calc. tab des voisins, periode de chgt du pas en temps, poeriode de correction de la derive
  real(double) :: tempstop, tempstopcel,ttol, tfroi, tcooling, tcou, tfcou, epcou, &! temperature d'arret, max, visee si max, taux de refroidissement, temp de la couche externe et epaisseur
       tsfact, vmax, tgc, dfpred ! gestion du pas en temps
  real(double)::maxtcel
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
  integer::ibrake   ! electronic slowing in cascades
  integer::ngrdel
  real(double):: elosselec,elosselec1 ! electronic losses for all atoms ; the PKA
  real(double),pointer::elstopforce(:,:,:)



  real(double) :: pist, temp, pmean, tmean, kine, kinemean ! pression temp et moyennes associees


  integer :: nvois   ! nb de voisins max dans toute la boite = nb d'atome * nb de voisins (/2)
  integer, pointer,dimension(:) :: indi ! table des voisins
  integer, pointer, dimension(:) :: indi2 ! table de voision pour les constantes de force
  real(double) :: rvois ! rayon de la table des voisins
  logical :: ltabvois      ! table des voisins ?
  logical :: lconstrtot ! construction par double boucle (T) ou par cel (F)
  logical :: ldemitab ! construction d'une demi-table (T) ou d'une table complete (F)
  character :: nature*6 ! element chimique
  integer :: nvat
  !EWALD
  real(double), dimension(:,:,:),pointer :: tabv3
  real(double), dimension(:,:,:,:),pointer :: tabf3


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

  real(double), parameter :: qmax=12
  real(double), parameter :: qmin=0.7
  real(double), parameter :: increq=0.1d+8

  real(double)::strucfact

  integer, dimension(:,:), pointer :: voisins
  real(double), parameter :: rcut=12e-8    !cutoff pour le calcul de S(q)

  integer, parameter :: cont = 1000


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



  integer :: imf     ! 
  integer :: imana     ! 

  logical :: ldislo  ! calcul de dislocation
  real(double) :: epcoudis,& !epaisseur de la couche avec ajout de force pour dislo
       &fdislo ! force appliqu
  integer, pointer :: latdebord(:)

  logical :: lcontr    ! dynamique contrainte (routine contrainte)
  logical :: lperiod   ! conditions periodiques
  logical :: lsuivinonpbc


  !---inNEB
  integer  :: ipath, npath,nebtype,nebrelaxation,maxneb,iteanaposneb, &
              neb_noise,mdcg_noise
  REAL(double) :: kspring,deltaRmax,neb_noise_scale,mdcg_noise_scale
  LOGICAL :: lPathFromGin      !if T : read initial path in gin files *.1.gin, *.2.gin, ... (NEB calculaion)
  !...inNEB

  !...inPHONDY
  integer  :: HessianOrder
  !...inPHONDY
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

  real(double)  :: kappa,text_teledyn,lanczos_step
  integer       :: niteration,nchemin_teledyn


  ! ------------- c o n d.   a u x   l i m i t e s   s p e c i f i q u e s ----------
  integer :: ibound         !*utilisat0 of specific bound condit0 (free or rigid)  !*!
  real(double) :: thickness	    ! epaisseur de la surface (ibound = 1, 2 ou 3)
  real(double) :: gap		    ! gap entre surface sup et surface inf
  real(double) :: thick_cryst       ! epaisseur de la surface (en coord. cristallines)
  real(double) :: layer_surf        ! aire de la surface XZ
  real(double) :: Lx_cm
  real(double) :: Lz_cm
  integer      :: i_surfINF
  integer      :: i_surfSUP
  integer      :: i_surfMAX
  integer      :: i_surfMIN
  real(double) :: y_max
  real(double) :: y_min
  real(double) :: y_2nd_max
  integer, dimension(:), pointer  :: b2sINF ! appartenance à la surface inférieure !*!
  integer, dimension(:), pointer  :: b2sSUP ! appartenance à la surface supérieure !*!
  logical      :: flag_fin
  real(double) :: ef_strain
  logical      :: ldecal_bc
  real(double) :: decal_bc		! pour les dislocations vis - decalage selon X
				! pour des potentiels EAM (implementé pour calfoeamtabvois)
  integer      :: itespebcout   ! frequence a laquelle on genere des .cfg (films mvt de dislo)
  real(double) :: inXMdis1
  real(double) :: inXMdis2
  logical      :: ldyn2D        ! .true. ->  dynamique 2D   ;    .false. ->  bords libres (par défaut)

! cas ibound = 1 :
  real(double):: user_strainrate   !*strain rate choosen by the user	     !*!
  real(double):: speed_user 	   ! surface atom speed (depends on user_strainrate)!*!

! cas ibound = 2 :
  real(double):: user_stress_yz	   !*stress applied on the cryst. surface	     !*!

! cas ibound = 3 : 
 !real(double):: user_strainrate   !*est aussi nécessaire		     !*!
 !real(double):: user_strainrate   !*est aussi nécessaire		     !*!
  real(double):: currentstress     ! stress applied - corrected at each time step !*!
  real(double):: fdbkcoef	   !*coef de la boucle de feedback de correct0 de currentstress
  real(double):: forceatsup        ! force on sup. surface atom (stress controlled) !*!
  real(double):: forceatinf        ! force on inf. surface atom (stress controlled) !*!

!      RQ:  les parametres comportant une * devant leur description doivent etre
!      definis dans le fichier .din
! ----------------------------------------------------------------------------------
  


end module gen_com_m
