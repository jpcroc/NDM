module gen_com_m
  USE T_kind_param_m
  implicit none


#ifdef PARA

  ! Declarations MPI
  !      include 'mpif.h'
  ! option pour realiser une trace Vampir
  !     include 'VT.inc'
  !      integer rang,code,nb_procs
  !      common/parampi/rang,code,nb_procs

#else

  !      integer :: rang, tranche
#endif

  integer::nvperat


  integer :: rang, rangph, rangml, rangmab, ja_phondy, ja_ml
  logical :: parallele

!  integer :: natperc                        ! nb d'atome par cel

  integer :: ivoismax

!  integer :: imm                 !imm taille des tableaux dependant du nombre d'atome mélange par proc et pas par proc, à remplacer par les deux suivants
!  integer :: imm_glob ! taille complète des tableaux 
  integer:: imm_loc ! taille des conf atomique par proc 


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


!  integer :: im						! nb local d'atomes (=global en sequentiel)
  integer :: im_glob,imm_glob					! nb global d'atomes


!
!  real(double) :: volu      ! volume
!  integer, dimension(3) :: lat      ! generation: nb de repetition de cel unite
!  real(double), dimension(3,3) :: at, bg ! at : vecteurs de base de la boite (BOND en cm) bg: vecteur du reseau reciproque
!  real(double), dimension(3,3) :: h0     ! Vecteurs de base de la boite de reference en A (Parrinello, Rahman)
  logical :: lUcell                 ! affiche l'energie potentielle de la boite
  ! (cela suppose que h0 corresponde a l'etat de reference pour lequelle la contrainte est nulle)     

  logical :: lfrozen    ! .true.: certains atomes sont bloques (pas de dynamique)
  logical :: lbulle    ! .true.: bulle
  logical :: ldesinteg    ! .true.: insertion appelle init_insert
  logical :: lrctest    ! .true.: test sur rc ; false pas de test
  integer:: nstepdes,ides, pm1des,itdes,imdesup,imdesdeb,imdesdn
  real(double)::lambdades,deltaF
  
  real(double),allocatable:: xpchup(:,:),vpchup(:,:),xpchdn(:,:),vpchdn(:,:),xpchdeb(:,:),vpchdeb(:,:)
  integer,allocatable :: itichup(:),itichdn(:),itichdeb(:),num_at_globdesup(:),num_at_globdesdn(:),num_at_globdesdeb(:)
  real(double):: kspr,xpspr(3),Espr,deltaEspr,xpspr0(3),tempdes
  integer:: typspr


!  logical, dimension(:), allocatable :: Free ! free(i)=.true. si l'atome i compte dans l'energie 
!  logical, dimension(:,:), allocatable :: Frozen ! Frozen(ix,i)=.true. si la coordonnee ix de l'atome i est gelee
  integer::imFree,imFirstFrozen ! nb d'atoems libres

!  real(double), dimension(3) :: normat ! norme de at



!!$  integer, dimension(:,:), allocatable :: ncel  ! ncel(i,j) indice de la jeme cel voisines de la cel i
!!$  integer, dimension(:), allocatable :: nato    ! nb d'atome dans la ieme cel
!!$  integer, dimension(:,:),allocatable :: atincel	! last (i,j) numero du ieme atome de la jeme cel
!!$  integer, dimension(:,:,:),allocatable :: deltadist ! decalage a appliquer sur la cel
!  integer :: nox, noy, noz, noxy, noxyz	     !nb de cel suivant x y z et total (DOIT REMPLACER nce)
!  real(double), dimension(3) :: celsize	     ! taille des cel


  integer :: it, itmax, nitmax,igen ! iteration courante, finale , type de generation
  real(double)::timemax ! temps max simul
  integer :: lenfnam
  integer :: fmt_cin

  character :: fnam*80, fnamout*80, fnamcout*80, fnamcoutxp*80,fnamcoutfp*80, fnamcoutnonpbcxp*80
  logical :: ltranche ! surface
  integer:: iteplz,nplz ! distribution suivant des tranches en z
  real(double):: rulayer

  integer :: imgs, imgi, itefrac !fracture IMD nombre d'atomes sur lesquels on fait la dynamique normale
  real(double) :: cougel, zincr !fracture



  real(double) :: potist ! energie potentielle totale
  real(double):: potisP,potis1, potis2, potis3, potis0, potcp ! energie potentielle de paire
  real(double) :: potisTersoff ! energie potentielle de tersoff

!  INTEGER::imd ! HISTORIQUE A DEGAGER LE PLUS TOT POSSIBLE

  real(double) :: oldtstep  
  real(double) :: tstep, usdh, timel  
  integer :: itetemp, itesigma, itedepla, itecoordo, iterdf, nrdf, & 
       iterasmol, iteangle,nfda,itetemp2,iteanapos, itefcc,itecfg
  integer::ivisu     ! format de sortie dans rasmol.f90 : ivisu=1=.mol, ivisu=2=vsim mal codﾃｩ, ivisu=2=xred
  real(double)::rcangle,rcrdf

  real(double), dimension(3,3) :: sig ! contrainte
  real(double), dimension(3,3) :: sigtot
  real(double), dimension(3,3) :: sigkine

!  real(double), dimension(:,:,:),allocatable :: sigc ! contrainte par cel
!  real(double), dimension(:,:,:),allocatable :: sigat ! contrainte par atome
  logical :: lEparat  ! calcul et affichage dans rasmol de la contrainte atomique; affichage ﾂｩnergie par atome,calcul bond valence
  integer:: itebdv ! frequence de calcul des bond valence
  logical :: ljqbh ! calcul de la conductivitﾃδｩ thermique par la mﾃδｩthode directe
  logical :: lnemd  ! calcul de la conductivitﾃδｩ thermique par NEMD
  integer ::njqbh,ittherm,ntr
  real(double) :: epsil,epcoud,kthg


  logical :: lFire      ! If true (default), fire algorithm is USEd for quenching
  real(double)::fnemd
  real(double):: fpstop ! critere de conv. sur la force par atome max  pour les trempes UNITE = EV/ANG
  real(double):: sigstop ! critere de conv. sur les composantes de contraintes  pour les trempes UNITE = kbar
  real(double):: fsumstop ! critere de conv. sur la force sqrt ( sum_f F_i^2 )  pour les trempes UNITE = EV/ANG
!  real(double),allocatable::eatom(:) ! energie par atome
!  real(double),allocatable::eatomtotm(:) ! energie par atome
  logical :: lPrtSigat, lprteat, lprtfat,lprteattotm  ! calcul et ecriture de la contrainte, l'energie et force par atome, de l'energie par atome totale (pot+cin) moyenne
  logical :: lsigatcel !ecriture de la contrainte atomique moyenne sur cellule
  logical :: lsigat ! la contrainte atomique est calcul馥 (rendu vrai par lprtsigat ou lsigatcel)
  logical :: lax ! stockage positions initiales
  logical :: lposmoy ! ecrit a la fin la position moyenne des atomes
  real(double) :: tdepla, tdepla2 ! seuils de deplacement
  logical :: lfilm, linstantrdf,linstantfda, lrestart, ltpcel, lfilmext !film, RDF, restart, moyenne par cel
  logical:: ldecoup !if T: cherche les nombres de procs optimums, voir decoup3D (ne marche su'en séquentiel (évidemment)) 
  real*8,dimension(4)::tpseuils ! 1:Tmin; 2:abs(T') ; ; 3:abs(P); 4:abs(P')
  !Correlations et Cie
  logical :: lcorrelvp ! ecriture de l'autocorrelation des vitesses*Masses
  logical :: lcalcjq

  logical :: lcdp ! algorithme d'accumulation de defauts ponctuels
  logical :: lheat ! algorithme de chauffage local
  real(double)::rheat,theat,Eheat
  integer :: iteheat


  real(double) :: tinit !temp initiale
  real(double)::tempdeplainit,debyetemp
  logical :: lvpread  ! vitesse lue dans le fichier .cin
  integer:: iseed ! graine du gerateur aleatoire des vitesses
  integer :: dmtype, itab, itetabvois, itetimestep, itederive ! type dynamique, periode de repartition entre cel, periode de calc. tab des voisins, periode de chgt du pas en temps, poeriode de correction de la derive
  real(double):: depmaxts,tsmin
  real(double) :: tempstop, tempstopcel,ttol, tfroi, tcooling, tcou, tfcou, epcou, &! temperature d'arret, max, visee si max, taux de refroidissement, temp de la couche externe et epaisseur
       tsfact, vmax, tgc, dfpred ! gestion du pas en temps
  real(double)::maxtcel
  real(double) :: deltaestop ! decroissance de la temperature moyenne
  integer :: nbmoye
  integer :: ibordcou
  integer :: itesauv, formatsauv, itesauvposition, itesauvforce,itesauvinter  ! periode de sauvegarde format de sauvegarde periode 
                                                                              ! de d'ecriture des positions et/ou forces en formatted ; 
  !itesauvinter=sauvegarde reguliere .cout.it qui n'efface pas les fichiers .cout precedent
  logical::lWgin ! ecriture finale de .newgin
  real(double), dimension(3) :: vh ! vitesse de la boite
  real(double) :: pext, wbox, tbox ! pext poids de la boite temps d'amortissment de la boite
  logical ::  lpcon2,lprtzlm ! pression constante sans et avec amortissement
  logical ::  lpconxyz      ! pression constante - buy only the diagonal term of box matrix can  change. 
  logical :: lTcon, lTberendsen,lTandersen,lTNose,lTHoover,landerscou ! temp constante (3 algorithmes differents)
  real(double) :: Text ! T exterieure
  logical :: lLangevin ! Langevin MD
  real(double) :: gamlg  ! gamma de Langevin
  integer :: ilangevin ! 1=std ; 2=Athenes
  real(double) :: tauTcon ! Temps berendsen
  real(double) :: nuandersen ! frequence andersen     
  integer :: iteTconst! periode de temp const
  integer :: iko ! indide du PAF
  real(double) :: eko, xko, yko, zko ! energie et direction du PAF
  real(double) :: xx0, yy0, zz0 ! position initiale du projectile
  logical :: lcasca,lderive ! cascade,correction derive ?



  real(double) :: pist, temp, pmean, tmean, kine, kinemean ! pression temp et moyennes associees
  real(double) :: tempEP ! temperature for slow moving atoms (EP=elec-phon)

!  integer :: nvois   ! nb de voisins max dans toute la boite = nb d'atome * nb de voisins (/2)
!  integer, allocatable,dimension(:) :: indi ! table des voisins
!  integer, allocatable, dimension(:) :: indi2 ! table de voision pour les constantes de force
!  real(double) :: rvois ! rayon de la table des voisins
!  logical :: ltabvois      ! table des voisins ?
  logical :: lconstrtot ! construction par double boucle (T) ou par cel (F)
  logical :: ldemitab ! construction d'une demi-table (T) ou d'une table complete (F)
  character :: nature*6 ! element chimique
  integer :: nvat
  !EWALD


  logical :: lalea  ! preparation d'une configuration aleatoire
  logical :: lopt   ! optimisation de Ewald par PME si TRUE
  real(double) :: rsep  !Distance de separation pour le tirage aleatoire

  ! NVT, NPT ensembles
  !      logical :: lnose, lnosepar,lpr ! lnose =Tcst  la Nose ; lnosepar=T&P cst a la Nose Parinello Rahman
  logical :: lpr,lprtrp ! l Parinello Rahman
  !      real(double) :: tomega, tbomega ! mass fictive du thermostat et du piston
  real(double), dimension(3,3) :: att, ati    !vitesse de la forme de la boite ; ati=(at^-1)
  real(double), dimension(3,3) :: ihbox0      ! the degree of freebom of the box. If is 1 everywhere all the shape  can change.
                                              ! If you put on diagonal 1 and the rest is 0 you can chage only anlong x,y and z.  

  real(double), dimension(3,3) :: sigext, pext_hydro  !contraintes externes appliques; contraintes calculees

  integer :: lurdfout
  real(double) :: lastcool


  real(double), parameter :: rmin = 0.5d-8

  real(double), parameter :: qmax=12
  real(double), parameter :: qmin=0.7
  real(double), parameter :: increq=0.1d+8

  real(double)::strucfact

  integer, dimension(:,:), allocatable :: voisins
  real(double), parameter :: rcut=12e-8    !cutoff pour le calcul de S(q)

  integer, parameter :: cont888 = 1000


  real(double), parameter :: thetamin = 1.0D-7
  real(double), parameter :: thetamax = 6.2

  logical lEev,lPkbar   !unite
  real(double) :: unitE,unitP
  character*5 :: cunitE, cunitP


  ! energies potentielle, cinetique et totale de la boite en Parrinello-Rahman
  real(double):: EcellPR, Kcell, Ucell      


  !Variables Nose
  REAL(double) :: ENose, KNose, UNose
  REAL(double) :: fNose
  real(double) :: wNose     ! Poids associe au thermostat de Nose
  INTEGER :: nHoover        ! Nombre de chaines de Hoover
  REAL(double), dimension(:), allocatable :: zHoover   ! Viscosite 

  real(double) :: deltax 
  real(double) :: dilat(3)




  logical :: ldislo  ! calcul de dislocation
  real(double) :: epcoudis,& !epaisseur de la couche avec ajout de force pour dislo
       &fdislo ! force appliqu
  integer, allocatable :: latdebord(:)

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
  logical, dimension(:), allocatable :: cyl ! tableau pour savoir si atome dans cylindre

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
  integer, dimension(:), allocatable  :: b2sINF ! appartenance a la surface infﾃｩrieure !*!
  integer, dimension(:), allocatable  :: b2sSUP ! appartenance a la surface supﾃｩrieure !*!
  logical      :: flag_fin
  real(double) :: ef_strain
  logical      :: ldecal_bc
  real(double) :: decal_bc		! pour les dislocations vis - decalage selon X
				! pour des potentiels EAM (implemente pour calfoeamtabvois)
  integer      :: itespebcout   ! frequence a laquelle on genere des .cfg (films mvt de dislo)
  real(double) :: inXMdis1
  real(double) :: inXMdis2
  logical      :: ldyn2D        ! .true. ->  dynamique 2D   ;    .false. ->  bords libres (par defaut)

! cas ibound = 1 :
  real(double):: USEr_strainrate   !*strain rate choosen by the user	     !*!
  real(double):: speed_USEr 	   ! surface atom speed (depends on user_strainrate)!*!

! cas ibound = 2 :
  real(double):: USEr_stress_yz	   !*stress applied on the cryst. surface	     !*!

! cas ibound = 3 : 
 !real(double):: USEr_strainrate   !*est aussi necessaire		     !*!
 !real(double):: USEr_strainrate   !*est aussi necessaire		     !*!
  real(double):: currentstress     ! stress applied - corrected at each time step !*!
  real(double):: fdbkcoef	   !*coef de la boucle de feedback de correct0 de currentstress
  real(double):: forceatsup        ! force on sup. surface atom (stress controlled) !*!
  real(double):: forceatinf        ! force on inf. surface atom (stress controlled) !*!

!      RQ:  les parametres comportant une * devant leur description doivent etre
!      definis dans le fichier .din
! ----------------------------------------------------------------------------------
  
!  real(double), dimension (:),allocatable ::tempc,tempcm,celpm1,tm1,celpp,tcp,pmc,patcel,patcelmax
  real(double), dimension (:,:,:),allocatable ::sigatcel
  integer, dimension (:),allocatable ::natchk
  logical, dimension (:),allocatable ::lprtcel(:)


!ZBL 
  real(double)::potiszbl ! energie pot de ZBl quand ajoute independement 
  
  logical :: l2T
  real(double), dimension (:),allocatable ::elossCel

  logical ::lspaceNDM ! TRUE= para space NDM/ false= paraspace LAMMPS
  character (len=15):: units_lammps
  real(double)::rskin,position_conversion_lammps, energy_conversion_lammps, pressure_conversion_lammps ! epaisseur pour lammps (equivalent rvois-rue)
 real(kind=8) , allocatable, dimension(:)  ::  posa, forca
 logical :: firsttime_lammps
 integer:: iverbose ! verbosity (0 = pas de détails, défaut, 1 = détails)


 logical :: latcomp ! masters (myidsp=0) have the complete positions (for sauvegardeT), rasmolT
 
end module gen_com_m
