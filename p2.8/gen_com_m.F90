module gen_com_m
  USE T_kind_param_m
  implicit none


  integer,target :: rang, rangph, rangml, rangmab, ja_phondy, ja_ml
  logical :: parallele

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


  integer :: imm_glob 

  real(double), dimension(3,3) :: h0     ! Vecteurs de base de la boite de reference en A (Parrinello, Rahman)
  logical :: lUcell                 ! affiche l'energie potentielle de la boite
  ! (cela suppose que h0 corresponde a l'etat de reference pour lequelle la contrainte est nulle)     
  logical :: lrctest    ! .true.: test sur rc ; false pas de test

  integer,target :: it ! iteration courante, finale , type de generation
  integer:: itloopmax, itmax, nitmax,igen ! iteration fin de boucle DM, finale , type de generation
  real(double)::timemax ! temps max simul
  integer :: lenfnam
  integer :: fmt_cin

  character :: fnam*80, fnamout*80, fnamcout*80, fnamcoutxp*80,fnamcoutfp*80, fnamcoutnonpbcxp*80
  real(double):: rulayer

  integer::idirectionmcgc

  real(double),target :: potist ! energie potentielle totale
  real(double):: potisP,potis1, potis2, potis3, potis0, potcp ! energie potentielle de paire
  real(double) :: potisTersoff ! energie potentielle de tersoff

!  INTEGER::imd ! HISTORIQUE A DEGAGER LE PLUS TOT POSSIBLE

  real(double) :: oldtstep  
  real(double) :: tstep, usdh, timel  
  integer :: itetemp, itesigma, iteprtsigma,itedepla, itecoordo, iterdf, nrdf, & 
       iterasmol, iteangle,nfda,itetemp2,iteanapos, itefcc
  integer::ivisu     ! format de sortie dans rasmol.f90 : ivisu=1=.mol, ivisu=2=vsim mal codﾃｩ, ivisu=2=xred
  real(double)::rcangle,rcrdf

  real(double), dimension(3,3),target :: sig ! contrainte
  real(double), dimension(3,3) :: sigtot
  real(double), dimension(3,3) :: sigkine

  logical :: lEparat  ! calcul et affichage dans rasmol de la contrainte atomique; affichage ﾂｩnergie par atome,calcul bond valence
  integer:: itebdv ! frequence de calcul des bond valence
  logical :: ljqbh ! calcul de la conductivitﾃδｩ thermique par la mﾃδｩthode directe
  logical :: lnemd  ! calcul de la conductivitﾃδｩ thermique par NEMD
  integer ::njqbh,ntr
  real(double) :: epsil,epcoud


  logical :: lFire      ! If true (default), fire algorithm is USEd for quenching
  real(double)::fnemd
  real(double):: fpstop ! critere de conv. sur la force par atome max  pour les trempes UNITE = EV/ANG
  real(double):: sigstop ! critere de conv. sur les composantes de contraintes  pour les trempes UNITE = kbar
  real(double):: fsumstop ! critere de conv. sur la force sqrt ( sum_f F_i^2 )  pour les trempes UNITE = EV/ANG
  logical :: lPrtSigat, lprteat, lprtfat,lprteattotm  ! calcul et ecriture de la contrainte, l'energie et force par atome, de l'energie par atome totale (pot+cin) moyenne
  logical :: lsigatcel !ecriture de la contrainte atomique moyenne sur cellule
  logical :: lsigat ! la contrainte atomique est calcul馥 (rendu vrai par lprtsigat ou lsigatcel)
  logical :: lax ! stockage positions initiales
  logical :: lposmoy ! ecrit a la fin la position moyenne des atomes
  real(double) :: tdepla, tdepla2 ! seuils de deplacement
  logical :: lfilm, linstantrdf,linstantfda, lrestart, ltpcel, lfilmext !film, RDF, restart, moyenne par cel
  logical :: lrestartmcgc
  logical:: ldecoup !if T: cherche les nombres de procs optimums, voir decoup3D (ne marche su'en séquentiel (évidemment)) 
  real*8,dimension(4)::tpseuils ! 1:Tmin; 2:abs(T') ; ; 3:abs(P); 4:abs(P')
  !Correlations et Cie
  logical :: lcorrelvp ! ecriture de l'autocorrelation des vitesses*Masses
  logical :: lcalcjq

  logical :: lcdp ! algorithme d'accumulation de defauts ponctuels
  real(double) :: tinit !temp initiale
  real(double)::tempdeplainit,debyetemp
  logical :: lvpread  ! vitesse lue dans le fichier .cin
  integer:: iseed ! graine du gerateur aleatoire des vitesses
  integer :: dmtype, itab, itetimestep, itederive ! type dynamique, periode de repartition entre cel, periode de calc. tab des voisins, periode de chgt du pas en temps, poeriode de correction de la derive
  integer,target :: itetabvois
  real(double):: depmaxts,tsmin
  real(double) :: tempstop, tempstopcel, tcou, tfcou, epcou, &! temperature d'arret, max, visee si max, taux de refroidissement, temp de la couche externe et epaisseur
       tsfact, vmax, tgc, dfpred ! gestion du pas en temps
  real(double)::maxtcel
  integer :: ibordcou
  integer :: itesauv,  itesauvposition, itesauvforce,itesauvinter  ! periode de sauvegarde periode 
                                                                              ! de d'ecriture des positions et/ou forces en formatted ; 
  !itesauvinter=sauvegarde reguliere .cout.it qui n'efface pas les fichiers .cout precedent
  logical::lWgin ! ecriture finale de .newgin
  real(double), dimension(3) :: vh ! vitesse de la boite
  real(double) :: pext, wboxf, tbox ! pext poids de la boite temps d'amortissment de la boite
  logical ::  lpcon2 ! pression constante sans et avec amortissement

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
  logical :: lconstrtot ! construction par double boucle (T) ou par cel (F)
  logical :: ldemitab ! construction d'une demi-table (T) ou d'une table complete (F)
  integer :: nvat
  !EWALD


  logical :: lopt   ! optimisation de Ewald par PME si TRUE

  ! NVT, NPT ensembles
  logical :: lprahman,lprtrp ! l Parinello Rahman
  real(double), dimension(3,3) :: att, ati    !vitesse de la forme de la boite ; ati=(at^-1)
  integer, dimension(3,3) :: ihbox0 ! integer pour bétonner les tests      ! the degree of freebom of the box. If is 1 everywhere all the shape  can change.

  real(double), dimension(3,3) :: sigext, pext_hydro  !contraintes externes appliques; contraintes calculees

!ci-dessous choses à modulariser
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

  logical :: lcontr    ! dynamique contrainte (routine contrainte)
  logical,target :: lperiod   ! conditions periodiques
  logical :: lsuivinonpbc


  !---inNEB
  integer  :: ipath, npath,nebtype,nebrelaxation,maxneb,iteanaposneb, &
              neb_noise,mdcg_noise
  REAL(double) :: kspring,deltaRmax,neb_noise_scale,mdcg_noise_scale
  LOGICAL :: lPathFromGin      !if T : read initial path in gin files *.1.gin, *.2.gin, ... (NEB calculaion)

  ! selection des atomes distordus

!  real(double), dimension (:),allocatable ::tempc,tempcm,celpm1,tm1,celpp,tcp,pmc,patcel,patcelmax
  real(double), dimension (:,:,:),allocatable ::sigatcel
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
