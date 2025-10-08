! *****************************************************************
module readdm_mod
  use read_val
     USE arret_ndm_mod,only:arret_ndm
  implicit none
contains


  subroutine readdm
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    use Tpara,only:nprocs,mpi_world
    USE gen_com_m, ONLY:a2cm,debyetemp,deltax,depmaxts,dfpred,gamprfact,&
         &epcou,ev2erg,fmt_cin,fpstop,fsumstop,gamlg,couxyz,&
         &igen,ilangevin,iseed,itab,itederive,&
         &itesauvforce,itesauvposition,itetabvois,itetconst,itetimestep,&
         &landerscou,lcdp,lconstrtot,lcorrelvp,lderive,lfire,&
         &ljqbh,lpcon2,lPcube,lrctest,lrestart,ltandersen,&
         &ltcon,lvpread,mdcg_noise,&
         &nhoover,nitmax,nuandersen,pext,&
         &rskin,rulayer,sigext,sigstop,tbox,tempdeplainit,tempstop,tempstopcel,tgc,&
         &timemax,tinit,tsfact,tsmin,two,units_lammps,usdh,utemps,wboxf,wnose,&
         &ihbox0,cunite,cunitp,dmtype,erg2ev,fnemd,&
         &iteanapos,iteangle,itebdv,itecoordo,itedepla,&
         &iterasmol,iterdf,itesauv,itesauvinter,itesigma,iteprtsigma,itetemp,itetemp2,itmax,ivisu,l2t,lcalcjq,&
         &lcasca,lcontr,ldemitab,leev,leparat,lfilm,linstantfda,linstantrdf,&
         &llangevin,lnemd,lperiod,lpkbar,lposmoy,lprahman,lprteat,lprteattotm,lprtfat,lprtsigat,lsigat,lsigatcel,&
         &lsuivinonpbc,ltberendsen,lthoover,ltnose,ltpcel,lucell,lwgin,nfda,h0,&
         &nrdf,rang,rcangle,rcrdf,tautcon,tdepla,tdepla2,text,tfcou,iteprtkin&
         &,tpseuils,tstep,unite,unitp,lenfnam,fnam,lanaposart,lmultin,lbabar&
         &, lax,ldecoup,lspaceNDM,latcomp,dilat,lrestartmcgc,lspecialinit,lmaxvp,vplim
#ifdef LAMMPS_VERSION
     USE gen_com_m, ONLY: energy_conversion_lammps, position_conversion_lammps, pressure_conversion_lammps
#endif
    use read_val
    use WGC_mod,only:ndir,nstep,betaguess,ncgtry,lvarstop,fstpdecr,beta35,gammas,gammav
    USE var_pot, ONLY:lforcetabulate,lprtpot,maxorder,ngrid,npotentiel,eatref,ipotentiel,npotmax,ntyp,lpotentiel       
    USE jqmod
    USE eloss, ONLY : tcelec,ecelec,ibrake,ngrdel
    USE arret_ndm_mod,only: arret_ndm
    use neb_module,only: lvzeroneb,kspring,lclimb,nwclimb,nebrelaxation,npath,lpathfromgin,iteanaposneb,&
         &neb_noise,neb_noise_scale,mdcg_noise_scale,maxneb,deltarmax,nebtype,i_neb_drag
    USE montecarlo_mod, ONLY: pas_lambda_mc,distminat,n_path,lparapath, nparapath,idirectionmcgc, &
         &lbiais_retrait,lbiais_inser, fdmc_1, fdmc_2,nbatplus,itypcalc,R0mcgc,fdfactmcgc,ins_typ,bublcenter,&
         &typswitch1,typswitch2,izlins,zlcenter,lspring,k_spring,protocol_mcc
    use ForceMatrix_mod,only: ndecal,decal,lparafm,nparafm,lwritefreq,lwfm
    use Parrinello_Rahman,only:TinitBox
    use constrconf_mod,only: ldecalcor
    use arps_mod,only:kmin,kmax,noxyzkmin,noxyzkmax,lpartarps!,lxyz
    use babar_mod,only:ntempbabar,nbabarprocs,bbtempmin,bbtempmax



    ! *****************************************************************

    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: ludin, lufilm, lufilmpaf,  i,itean, ic,ic2, iThermo,itecompcr,ipotcont
    character :: fnamdin*80
    logical :: lginread,ltriclin,lpcon,lfissure,tpot,lpr,lseedcom
    integer::itecfg,np2
    logical :: lpconx,lpcony,lpconz,lpconxyz,ltest,lseecom
    !-----------------------------------------------
    !
    !
    !
    !

    namelist /input/itab, itetabvois, itetemp, itesigma,iteprtsigma,  itedepla, tdepla, lfilm, &
         tempstop, tempstopcel,dmtype, lFire,  itecoordo, tstep, itetimestep, tsfact, &
         tinit,  tfcou, epcou, couxyz,lcasca, lfissure, itmax,nitmax, itean, kspring,  &
         itederive, igen, linstantrdf, iterdf, nrdf,nfda, linstantfda, itesauv,  &
         lrestart, lPathFromGin, tgc, ltabvois, rvois, rskin,ltpcel, nox, noy, noz, imm, dfpred, &
         rulayer,iterasmol, lpcon, pext, wboxf, wNose, lpcon2, lpconxyz,lpconx,lpcony,lpconz, tbox, &
         iteangle,  itesauvposition, itesauvforce,  tdepla2, lpcube,&
         lTcon,Text,iteTconst, lTberendsen, lTNose, lTHoover, nHoover, tauTcon, &
         maxorder, ipotentiel,lpotentiel,beta35,R0mcgc,izlins,zlcenter,fdfactmcgc,ins_typ,bublcenter,&
         h0, sigext,lconstrtot,lEev,lPkbar,deltax,lcorrelvp,lvpread,&
         lcalcjq,dilat,lderive,lTandersen,nuandersen,landerscou,Llangevin,gamlg,ilangevin,&
         lcdp, ljqbh,lEparat,itebdv,itetemp2,itecompcr,iteanapos,&
         lnemd,fnemd,fpstop,iseed,fsumstop,sigstop,lcontr,lpr,lUcell,ngrid,lperiod,&
         lprteat,lprteattotm,lprtfat,lprtsigat,lsigatcel,itecfg,npath,nebtype,nebrelaxation,maxneb,deltaRmax,&
         rcangle,rcrdf,fmt_cin,lginread,ltriclin,iteanaposneb,ntyp,ihbox0,&
         neb_noise,neb_noise_scale,lsuivinonpbc,lposmoy,gammas,gammav,lanaposart,&
         eatref,mdcg_noise_scale, mdcg_noise, lforcetabulate,ivisu,idirectionmcgc,nbatplus,&
         tempdeplainit,debyetemp,ibrake,lprtpot,ngrdel,timemax,tpseuils,lrctest,tcelec,Ecelec,l2T,depmaxts,tsmin,&
         itesauvinter,units_lammps,lWgin,lvzeroneb,lclimb,nwclimb,pas_lambda_mc,n_path,lax,ldecoup,distminat,&
         ndir,nstep,betaguess,ncgtry,lvarstop,fstpdecr,itypcalc,gamprfact,TinitBox,&
         &nparapath,lparapath,lrestartmcgc, lbiais_retrait,lbiais_inser,fdmc_1,lmultin,&
         &fdmc_2,ndecal,decal,lparafm,nparafm,lwritefreq,lwfm,ldecalcor,kmin,kmax,iteprtkin,lspecialinit,&
         &noxyzkmin,noxyzkmax,lpartarps,lspring,k_spring,i_neb_drag,protocol_mcc,lmaxvp,vplim,ntempbabar,&
         &nbabarprocs,bbtempmin,bbtempmax


    !
    !   set default values for variables in namelist
    !

    fnamdin = fnam(1:lenfnam)//'.din'
    ! variables de dynamique
    kspring=1

    lspaceNDM=.true.
    imm = 0                     !dimensionnement des tableaux atomiques
    itab = 10                   !period of cell repartition
    itetabvois = 10             !periode de calcul de la table des voisins
    tempstop = -1.0             !temperature of run stop
    tempstopcel = -1.0             !temperature of run stop
    dmtype = 0  
    idirectionmcgc=-2           !direction pour le montecarlo 0 ou 1 a designer par l'utilisateur
    lbiais_retrait = .false.    !biais ou non sur les retraits dans le montecarlo
    lbiais_inser = .false.    !biais ou non sur les retraits dans le montecarlo
    fdmc_1 = -1000.0            !param de fermi dirac A DEF PAR UTILISATEUR pour la fct discriminante du biais dans MC
    fdmc_2 = -1000.0            !valeur devant etre changee    
    !dmtype = type of calculation : 1 -> MD
    !                               2 -> quench (trempe) 
    !                               21 -> quench (trempe)Vcst
    !                               22 -> quench (trempe)PCst
    !                               23 ->  fire quench Vcst
    !                               24 ->  fire quench Pcst
    !                               3 -> gradient conjugue générique pointe vers 31 par défaut 
    !                              30 -> VIEUX gradient conjugue sur les coordonnes reduites
    !                              31 -> VIEUX gradient conjugue sur les coordonnes cartesiennes
    !                              32 -> steepest descent
    !                              33 -> gradient conjugue
    !                              34 -> gradient conjugue modifié Fletcher-Reeves
    !                              35 -> relaxation ADAMD. :Kingma and J. Ba, “Adam: A Method for Stochastic Optimization,” in International Conference on Learning Representations (ICLR), 2015.
    !                               4 -> Velocity Verlet
    !                               41 -> Velocity Verlet ARPS with partial forces
    !                               42 -> Velocity Verlet ARPS with regular forces
    !                               5 -> test des forces 
    !                               6 -> analyse des positions en fin de cascade 
    !                               7 -> calcul des phonons
    !                               8 -> PR
    !                               9 -> NEB
    !                              11 -> UN SEUL CALCUL DE FORCES
    !                              12 -> ART
    !                              16 -> BABAR
    !                              17 -> MAB
    !                              18 -> ML
    !                              19 -> matrice de forces
    !                              15 -> montecarlo_mCC
    !                              151 -> montecarlo_mcGC
    !                              112 -> histogramme des distances entre atomes
    !                              113 -> recherche de la position la plus éloigné des atomes
    lFire = .false.              ! Fire algorithm is USEd for quenching (cf tr_fire.F90)
    tstep = 1.0                 !timestep in 10^-15 sec unit
    itetimestep = -1            !period of check in timestep
    tsfact = 10.0               !change in time step factor
    tinit = -1.0                !initial temperature
    tinitbox = -1.0                !initial temperature
    tfcou = -1.0                !temperature of the border of the box
    epcou = -1.0                !width of the border of the box
    couxyz(:)=1
    lcasca = .FALSE.            !cascade Y/N
    lfissure = .FALSE.          !crack Y/N
    itmax = -1                  !maximum number of iterations
    nitmax = -1                 !maximum number of new iterations after restart
    itederive = -1              !"derive" correction
    igen = -2                 !type de generation :0 a partir de.gin, +1 a partir de .cin; -1 de gin vers cin puis stop +2 cintogin ; +3 modification de cin puis stop
    lrestart = .FALSE.          !if T : restarting from an interrupt job
    ldecoup = .FALSE.          !if T: cherche les nombres de procs optimums, voir decoup3D (ne marche su'en séquentiel (évidemment)) 
    lPathFromGin = .FALSE.      !if T : read initial path in gin files *.1.gin, *.2.gin, ... (NEB calculaion)
    tgc = 0.0                   ! threshold for CG calculation
    ltabvois = .FALSE.          ! methode de la table des voisins
    lconstrtot=.FALSE.           !!construction de la table des voisins T=double boucle F=via cel.
    rvois = 0.0                 ! rayon de la table des voisins
    rskin=0.4
    ltpcel = .FALSE.            ! output of temperature and stress in each cell
    lforcetabulate = .FALSE.    ! The derivative of the energy is NOT tabulated. TRUE if it is.


    itesauv = 1000               !period for saving
    itesauvposition = 0         !periode pour sauvegarde des positions en binaire
    itesauvforce = 0            !periode pour sauvegarde des forces en binaire
    itesauvinter=0
    fmt_cin=1                  !format des fichiers .cin 0 : initiale, 1 = para
    dfpred = 0.1            ! eguess for GC calculations and quenching
    nox = -1
    noy = -1
    noz = -1
    lprahman=.false.                 ! parinnelo rahman �あ contrainte constante
    lpr=lprahman
    ihbox0(:,:) = 1   ! all the dimension of the box can change
    sigext = 0.0                ! Symetric tensor related to the external stress
    !=== Modif Emmanuel Clouet ================
    h0(1:3,1:3) = 0.d0          ! Vecteurs de base de la bite de reference en A (Parrinello, Rahman)
    lUcell=.false.              ! affiche l'energie potentielle de la boite
    ! (cela suppose que h0 correspond a l'etat de
    ! reference, ie etat pour laquelle la
    ! contrainte est nulle)
    !=== Fin des modifications ================
    lpcon = .FALSE.             !algorithm a pression constante a la hache
    lpcon2 = .FALSE.            !amortissement de la deformation de la boite
    lpconxyz = .FALSE.          !the relaxation are allowed only along the X, Y and Z axis
    lpconx = .FALSE.          !the relaxation are allowed only along the X axis 
    lpcony = .FALSE.          !the relaxation are allowed only along the  Y  axis
    lpconz = .FALSE.          !the relaxation are allowed only along the  Z axis
    lpcube=.false.

    pext = 0.0                  !pression  par defaut
    wboxf = 1.0                  ! facteur masse de la boite pour Parrinello-Rahman (par defaut egale a 0.5*masse totale
    wNose = 0.0                 ! masse de la boite pour thermostat de Nose (par defaut egale a wbox)
    tbox = 1000.0               !"temps" de la boite
    lTcon=.false.               !algorithme a temperature constante
    lTberendsen=.false.         ! algorithme a temperature constante
    lTNose=.false.              ! algorithme a temperature constante de Nose
    lTHoover=.false.            ! algorithme a temperature constante de Hoover
    nHoover=1
    tauTcon=200.0               ! The rescales "time" for the Berendsen algorithm
    Text=-1.

    ipotentiel = -1              ! definit type potentiel : 0=Born-Mayer-Huggins, 1=Buckingham, 2=watanabe,3=buck8,4=UO2, 5 terme Morse, 6=SW �πｴﾎｵ縺､� la Vashista ; 7 pot paire tabule ; 10 EAM; 12 ZrC JuLi(+Tersoff Doan)  ; 13 Tersoff coupure COS; 14 Tersoff coupure FD ; 15 tersoff coupure SIN (original) ; 11 Ercollesi ;; -10=LAMMPS atom style atomic; -11 LAMMPS atom style charge (changes only simple.potin) ! 8 bandura 2017= Bukingham +Morse+Fermi-Dirac+Inverse gaussian ! 16= CRG ! 20=MILADY
    npotentiel = 1              ! nb de potentiels
    lpotentiel(:)=.false.
    ntyp=-1                    ! le nombre de type DOIT etre specifie si le nombre de potentiel est superieur �πｴﾎｵ縺､� 1
    ! PME
    maxorder=10                                ! Ordre du developpement maximal de la PME
    lvpread=.true.
    dilat(:)=0.0
    lderive=.false.
    lTandersen=.false.       ! temperature constante a la Andersen
    nuandersen=1.0d14        !frequence de tirage aleatoire des vitesses en Hz (valeur elevee = pour cascades)
    landerscou=.false.       ! Andersen seulement sur les bords
    lLangevin=.false.        ! Langevin MD
    gamlg =5d12            ! Gamma deLangevin (= 0.005/1d-15 fera vp*0.995 pour tstep=1d-15)
    gamprfact=0.1
    ilangevin=1

    lcdp=.false.             ! algorithme d'accumulation de defauts ponctuels
    lnemd=.false.  ! Kth par la methode NEMD Evans, P7229
    fnemd=0
    fpstop =-0.05 ! critere de conv. sur la force par atome max  pour les trempes UNITE = EV/ANG
    sigstop =-0.05 ! critere de conv. sur les contraintes par direction UNITE = kbar
    fsumstop =-0.1 ! critere de conv. sur la force sqrt ( sum_f F_i^2 )  pour les trempes UNITE = EV/ANG
    lcontr=.false. ! dynamique contrainte (routine contrainte)
    iseed=0 ! si <>0 controle le tirage aleatoire des vitesses
    !variables d'analyse

    itetemp = 100                !period of temperature calculation
    iteTconst =itetemp
    itesigma = -1               !period of stress calculation
    iteprtsigma = -1               !period of stress calculation
    itedepla = -100              !period of displacement cal.
    tdepla = 1.0                !threshold for displacement
    tdepla2 = -1.0              !second seuil pour calcul des atomes deplaces
    lfilm = .FALSE.             !film making of displaced atoms
    itecoordo = -100             !period of coordination calculation
    itean = 0                  !general control for analysis
    iterdf = -1                 !period of RDF calc. : -1 never ; 0 : nrdf last iterations; +iterdf every iterdf iterations
    linstantrdf = .FALSE.       !F: calculates and prints the average of the RDF, T : calculent RDFrdf iterations
    linstantfda = .FALSE.       !F: calculates and prints the average of the RDF, T : calcuatent RDFrdf iterations
    nrdf = 0
    nfda=0

    lEev=.true.
    lPkbar=.true.
    ! definition des rayons de coupure pour le calcul des coordinences autour de chaque type atomique
    deltax=0.0
    !    rclu=2.0


    iterasmol = -1                             ! <0 --> genere aucun fichier positions pour logiciel rasmol
    iteangle = -1                              ! pilote creation de fichier positions pour
    ! >=0 debut et fin d'execution

    lufilm = 89
    lufilmpaf = 79
    ludin = 94
    lcorrelvp=.false.
    lcalcjq=.false.
    lEparat=.false.             ! calul et affichage de l'energie par atom
    itebdv=-1  ! frequence de calcul des bond valence
    itetemp2=-1   ! frequence d'ecriture de la temperature dans fichier separe
    itecompcr=-1  ! remplacee par iteanapos
    iteanapos=-1  ! frequence de comparaison avec cristal de reference

    lperiod=.true.    ! si conditions periodiques (ipbc=1) et lperiod: coordonnées réduites entre 0 et 1. Si ipbc=1 et .not.lperiod coordonées peuvent dépasser
    lprteat=.false.   ! if you want to print the energy on atom
    lprtsigat=.false. ! calul et affichage de la contrainte sur chaque atome
    lsigatcel=.false. ! calul et affichage de la contrainte atomique moyenne sur la cellule
    lax=.false. ! stockage positions initiales
    lprteattotm=.false.   ! energie par atome totale (pot+cin) moyenne
    ngrid = 20000  ! taille de la grille des potentiels
    itecfg=-1   ! ecriture de fichiers .cfg pour AtomEye

    !---inNEB
    nebtype=2        ! NEB is the default
    nebrelaxation=2  ! We relax all the atoms; if nebrelaxation==1 only the most "deplaced" atoms      
    npath = 15       ! 15 images of the neb is the default
    maxneb = 700     ! the MAX of NEB steps
    deltaRmax=1.d-2
    neb_noise_scale=0.001      ! this will affect the 4th digit
    mdcg_noise_scale=0.001      ! this will affect the 4th digit
    ! x + x*neb_noise_scale*random,
    ! where "random" is a random number between 
    ! 0 and 1  
    neb_noise=0                 ! 0 without noise, 1 with noise
    mdcg_noise=0                ! 0 without noise, 1 with noise
    !      	lperiod=.false.  ! pas de conditions periodiques

    !...inNEB
    rcangle=3.0
    rcrdf=5.0
    ltriclin=.true.
    lprtfat=.false.

    iteanaposneb=0
    lsuivinonpbc=.false.  ! enable or disable a copy of non folded positions (by the pbc conditions)  in binary form each itetimestep. 
    lposmoy=.false.       ! writes the average position and energy of the atoms in a .mol file
    eatref(:)=0.
    lanaposart=.false.  ! anapos a la ART : decalage + defauts en WS, concu pour le cas des I dans UO2
    ivisu=4    ! format de sortie dans rasmol.f90 : ivisu=1=.mol, ivisu=2=vsim mal code supprime, ivisu=3=xred , ivisu=4 CFG, ivisu=6 xfg ; 7=xyz type à la Babel
    !4==> 40= pas de vitesses; 41 vitesses
    !6==> 60= pas de vitesses; 61 vitesses


    !   ----------------------------------------------------------------------------------------    !*!


    tempdeplainit=-1
    debyetemp=-1
    lprtpot=.false.
    ibrake =0   ! if =1 electronic slowing for cascades (acting on all atoms)
    ngrdel=500

    timemax=1d25
    
    tpseuils(:)=0 ! 1:Tmin; 2:abs(T') ; 3: abs(T'') ; 1:abs(P); 2:abs(P') ; 3: abs(P'')

    lrctest=.true.
    tcelec=0 ! temp鬧ｻature de coupure pour les pertes 鬧脇ctroniques
    Ecelec=0 ! temp鬧ｻature de coupure pour les pertes 鬧脇ctroniques


    l2T=.false. ! 2T model
    depmaxts=0.02
    tsmin=2.0

    units_lammps='TO_BE_SPECIFIED'
    lWgin=.false. ! =true écrit un fichier .newgin à la fin
    lvzeroneb=.false. ! si true , met vp à 0 ente chaque iteration neb (comportement pre ndm2020), defaut = false==> calcul plus rapide
    lparapath=.false.
    nparapath=1
    itypcalc=1
    pas_lambda_mc = -100 !valeur negative par defaut pour que l'utilisateur la change
    n_path = -100 !valeur negative par defaut pour que l'utilisateur la change
    distminat=1 ! distance minimale en Angstrom de l'atome inséré aux autres atomes en Monte-Carlo (défaut = pas de distance min=n'importe où)
    nbatplus=1
    ipbc(1:3)=1 ! 1=PBC; 2=wall... dimension 3 =plans bc; ac;ab

    ndecal=2 ! nombre de décalage dans le calcul de la matrice de force (dmtype=19)
    decal=0.1 ! décalage dans le calcul de la matrice de force (dmtype=19) (Angstroms)
    lparafm=.true.
    nparafm=nprocs
    lwfm=.false.
    lwritefreq=.true.


    ndir=50   !nombre de direction dans steepest descent
    nstep=50  ! nombre de pas dans la minimisation sur une ligne en steepes descent
    ncgtry=10  ! nombre de relaxations positions/celulle/positions/celluel,etc.
    betaguess=1d-6
    lvarstop=.false.
    fstpdecr=10.
    beta35=1d-10
    gammas=0.999
    gammav=0.9

    protocol_mcc="MCP"
    fdfactmcgc=18.0
    R0mcgc=-1.0
    bublcenter(:)=0.5
    izlins=-1
    zlcenter(:)=0.
    ins_typ=0   ! 0 everywhere, with bias) : 1 in a sphere, 3 in a slice ; with spring : 11 in a site, 33 in a plane , 44 in a line , 55 in a bubble
    typswitch1=0
    typswitch2=0
    ldecalcor=.true.
    kmin=0. ! min kinetic energy for arps
    kmax=0. ! max kinetic energy for arps
    lpartarps=.false.
    noxyzkmin(1:3)=-1
    noxyzkmax(1:3)=1000000
    iteprtkin=-1

    lspecialinit=.false. ! driver for specail initialization : cascade, press or heat burst etc.
    lspring=.false. ! MCC calculation with a slowly vanishing vabishing string
    k_spring=1.0 ! spring strenght (1 eV/Ang**2)

    lclimb=.false. ! set to true for climbin NEB
    nwclimb=3 ! starts the climbing at the second evaluation of forces (in VASP =1, in Henkelmann is set to "a few iterations")
    i_neb_drag=5

    lmaxvp=.false. ! if true velocities are caped at vplim in pr2.F90 (very crude way of stabilizing dynamics)
    vplim=5d6 
    

    ntempbabar=0
    nbabarprocs=ntempbabar
    bbtempmin=0; bbtempmax=0

    
    lmultin=.false. ! T==> reads multiple condfiguration files

    if (rang == 0) write (6, *) 'nom fichier din=', fnamdin

    open(unit=ludin, file=fnamdin, status='unknown', err=456)
    
    read (ludin, nml=input)

    rcangle=rcangle*1d-8
    rcrdf=rcrdf*1d-8
    lprahman=lpr

    distminat=distminat*1d-8
    tsmin=tsmin*1d-15
    depmaxts=depmaxts*1d-8
    if ((l2T).and.(tsmin==2d-15)) then
       tsmin=2.d-16
       !     depmaxts=0.002
    end if

    timemax=timemax*1d-15

    if (itecfg.gt.0) then
       write(6,*)'ITECFG desactive, reactivez (in readdm )"at your own risks"'
       write(6,*) 'utilisez ivisu=4 pour sortir des .cfg'
       ivisu=40
    end if
    if (ivisu==4) ivisu=40
    if (ivisu==6) ivisu=60
    if (rang == 0) then
       if (imm <= 0) then
          write (6, *) rang,'nombre d''atomes nul-> stop'
          call arret_ndm
       endif
    endif                                      ! fin rang=0


    if(lpcon) then
       if (rang==0) write(6,*)
       if (rang==0) write(6,*)'!!!!!!!!!!!!!!!!!!!'
       if (rang==0) write(6,*)'!!!!!!!!!!!!!!!!!!!'
       if (rang==0) write(6,*)'!!!!!!!!!!!!!!!!!!!'
       if (rang==0) write(6,*)
       if (rang==0) write(6,*)'lpcon n_existe plUSEst historique utiliser plutot lpr pour un Parinnello Rahman propre '
       if (rang==0) write(6,*)
       if (rang==0) write(6,*)
       if (rang==0) write(6,*)
       if (rang==0) write(6,*)'!!!!!!!!!!!!!!!!!!!'
       if (rang==0) write(6,*)'!!!!!!!!!!!!!!!!!!!'
       if (rang==0) write(6,*)'!!!!!!!!!!!!!!!!!!!'
       call arret_ndm
    end if

    if (any(ipbc.ne.1)) then
       do ic=1,3
          select case (ipbc(ic))
          case(1)
             if (rang==0) write(6,*)'direction', ic,' : regular PBC'
          case(2)
             if (rang==0) write(6,*)'direction', ic,' : WALL Boundary Conditions'
          case default 
             if (rang==0) write(6,*)'wrong bound ary conditions, stop'
             call arret_ndm
          end select
       end do
       if (dmtype==1) then
          dmtype=4
          if (rang==0) write(6,*)'dmtype changed from 1 to 4 for energy conservation'
       end if
    end if

    if(.not.lperiod) then
       if (rang==0)then
          write(6,*)
          write(6,*)'coordinates can get out of  0-1'
          if (ipbc(1)==1) write(6,*)'along  a'
          if (ipbc(2)==1) write(6,*)'along  b'
          if (ipbc(3)==1) write(6,*)'along  c'
       end if
    end if

    tstep = tstep*utemps
    tauTcon=tauTcon*utemps
    epcou = epcou*1D-8
    tdepla = tdepla*1D-8
    tdepla2 = tdepla2*1D-8
    if (itean .ne. 0) then
       itetemp = itean
       itesigma = itean
       itedepla = itean
       itecoordo = itean
       lfilm=.true.
       iterasmol=itean
       iterdf=itean
       linstantrdf=.true.
       iteangle=itean
       linstantfda=.true.
    endif



    if (lpkbar.EQV..true.)then
       pext=pext*1.0d9
       sigext=sigext*1.0d9
    end if


    ! input check
    if (lnemd) then
       lcalcjq=.true.
       if (fnemd==0) then
          if (rang==0) write(6,*) rang,'fnemd = 0 stop'
          call arret_ndm
       end if
       if (ljqbh) then
          if (rang==0) write(6,*) rang,'ljqbh et lnemd  stop'
          call arret_ndm
       end if
    end if
    select case (dmtype)
    case (23)
       lfire=.true.
    case(24)
       lfire=.true.
       lprahman=.true.
    end select
    if (lfire) then
       if ((dmtype==2).or.(dmtype==23).or.(dmtype==24))then
          if (lprahman) then
             dmtype=24
          else
             dmtype=23
          end if
       else if (dmtype==9) then
          if (rang==0) then
             write(6,*)'lfire NEB'
          end if
       else
          if (rang==0) then
             write(6,*)'lfire and not dmtype=2 ? remove lfire and choose dmtype'
             write(6,*)'dmtype=23 fire Vcst'
             write(6,*)'dmtype=24 fire LPRahman'
             write(6,*)'dmtype=21 fast quenching Vcst'
             write(6,*)'dmtype=22 fast quenching LPRahman'
             write(6,*)'dmtype=9 Fire NEB Vcst'
          end if
          call arret_ndm
       end if
    else
       select case (dmtype)
       case(2)
          if (lprahman) then
             dmtype=22
          else
             dmtype=21
          end if
       case(4)
          if (lprahman) then
             dmtype=8
             if (llangevin) dmtype=88
          end if
          
       end select
    end if


    select case(ipotentiel)
    case(-10,-11)
       npotentiel=1
       lspaceNDM=.false. ; latcomp=.true.
       if (rang==0) write(6,*)'POTENTIELS LAMMPS ; PARA_SPACE VERSION=LAMMPS NOT NDM !!'
    case(20)
       npotentiel=1
       lspaceNDM=.false. ; latcomp=.true.
       ltabvois=.true.
       itetabvois=100000000
       if (rang.eq.0) write (6, *) '    MILADY POTENTIALS',rvois
       if (rang==0) write(6,*)' PARA_SPACE VERSION=MILADY NOT NDM !!'
    case default
       latcomp=.false.
    end select

    



    if (lrestart) then
       igen = 1
       if (rang == 0) write (6, *) '****** RESTART FROM FILE **'
    endif
    if ((igen<-1).or.igen>2) then                    !+1 from file -1 generate then stop 0 generate then run
       if (rang==0) write (6, *) rang,'wrong igen stop'
       call arret_ndm
    endif

    if (itab <= 0) then
       if (rang==0) write (6, *) rang,'wrong itab < 1 '
       call arret_ndm
    endif
    if ((dmtype==41).or.(dmtype==42)) then
       ltabvois=.false.
       kmin=kmin*ev2erg; kmax=ev2erg*kmax
       if ((kmin==0.).or.(kmax==0.)) then
          write(6,*)'kmin==0. or kmax==0 '
          call arret_ndm
       end if
    end if
#ifdef PARA
    select case(dmtype)
    case(21,22,4,3,1,30,31,32,33,34,35,23,24,8,88,41,42)
       if (ltabvois) then
          select case (ipotentiel)
          case(20)
          case default
             ltabvois=.false.
             rvois=0
             if (rang==0) write(6,*)'LTABVOIS MIS A FALSE en PARA'
          end select
       end if
    case (19)
       if (nparafm.ne.nprocs) then
          if (ltabvois) then
             select case (ipotentiel)
             case(20)
             case default
                ltabvois=.false.
                if (rang==0) write(6,*)'LTABVOIS MIS A FALSE en PARA'
             end select
          end if
       end if
        
    case(9)
       np2=npath-2
       if (np2.ne.nprocs) then
          if (ltabvois) then
             select case (ipotentiel)
             case(20)
             case default
                ltabvois=.false.
                if (rang==0) write(6,*)'LTABVOIS MIS A FALSE en PARA'
             end select
          end if
       end if
    case(12)
       if (lparapath) then 
          if (mod(nprocs,nparapath).ne.0) then
             write(6,*)'nprocs/nparapath <>0 STOP'
             call arret_ndm
          end if
       end if
    case(15,151)
!!$       if ((lpr).and.((nox==-1).or.(noy==-1).or.(noz==-1))) then
!!$          if (rang==0) then
!!$             write(6,*)'MONTECARLO constant pressure : nox/noy/noz must be set in the din file STOP'
!!$          end if
!!$          call arret_ndm
!!$       end if
!       if (igen.ne.0) then
!          write(6,*)'IGEN MUST BE ZERO (dont know why) stop'
!          call arret_ndm
!       end if
       if (itypcalc.lt.0) then
          write(6,*)'itypcalc<0'
          call arret_ndm
       end if
       if (n_path.lt.0) then
          write(6,*)'n_path<0'
          call arret_ndm
       end if
       np2=nparapath*2
       if (np2.ne.nprocs) then
          if (ltabvois) then
             select case (ipotentiel)
             case(20)
             case default
                ltabvois=.false.
                if (rang==0) write(6,*)'LTABVOIS MIS A FALSE en PARA'
             end select
          end if
       end if
    case(16)
       if (ntempbabar==0) then
          if (rang==0) write(6,*)'DMTYPE=16 and NTEMPBABAR=0 stop'
          call arret_ndm
       end if
       if ((bbtempmin==0).or.(bbtempmax==0)) then
          if (rang==0) write(6,*)'DMTYPE=16 and bbtempmin/max=0 stop'
          call arret_ndm
       end if
       lbabar=.true.
    case default
       write(6,*)'DMTYPE',dmtype
       if (rang==0) write(6,*) 'WARNING : VERSION PARALLELE seulement avec ',&
& 'dmtype=21,22,4,3,9,15,151,1,30,31,32,33,34,19,35,23,24,41'
    end select
#endif

#ifdef DECOUP
    ltabvois=.false.;rvois=0.
#endif

    if ((.not.ltabvois).and. itab/=1) then
       !     if (rang == 0) then
       !        write (6, *) ' '
       !        write (6, *) 'Modification obligatoire a itab = 1 '
       !        write (6, *) ' '
       !     endif
       itab = 1
    endif

    if ((itmax.GT.0).and.(nitmax.GT.0)) then
       write(6,*)' NITMAX PASSE DEVANT ITMAX'
    end if
    if ((.not.lrestart).and.(nitmax.GT.0)) then
       write(6,*)'Nitmax>0 et pas restart ?'
       call arret_ndm
    end if
    if ((timemax.gt.0).and.(itmax==-1)) itmax=1000000000

    if ((dmtype.EQ.11).or.(dmtype.eq.15).or.(dmtype.eq.151).or.(dmtype.eq.9)) itmax=1

!    if (itmax < 0) then
!       if (rang==0) write (6, *) rang,'wrong itmax < 0 '
!       call arret_ndm
!    endif

    if (itedepla>=1 .and. tdepla<0.0) then
       if (rang==0) write (6, *) rang,'wrong tdepla < 0 '
       call arret_ndm
    endif

    !     if ((itedepla.ge.1).and.(tdepla2.lt.0.0)) then
    !       write(6,*)'wrong tdepla2 < 0 '
    !       call arret_ndm
    !     endif


    if (lfilm ) then
       if (itedepla < 1) then
          if (rang==0) write (6, *) rang,'contradiction itedepla <-> lfilm'
          call arret_ndm
       endif
    endif



    if (tstep<1D-20 .or. tstep>1D-13) then
       if (rang==0) write (6, *) rang,'mauvais pas en temps = ', utemps
       call arret_ndm
    endif

    if (tfcou*epcou <= 0.0) then
       write (6, *) rang,'probleme TFCOU EPCOU = ', tfcou, epcou
       call arret_ndm
    endif

    if ((ltriclin.EQV..false.).and.(rang==0))then
       if (rang==0) write(6,*)'LTRICLIN=FALSE N_EXISTE PLUS'
       if (rang==0) write(6,*)'INPUT HISTORIQUE ??'
    end if

    if(dmtype==9) lprteat=.true.
!!$    if(dmtype==12) then
!!$       lprteat=.true.
!!$       ltabvois=.false.
!!$    end if
    if(dmtype==16) lprteat=.true.
    if (lposmoy.EQV..true.) then 
       lprteattotm=.true.
       write(6,*)'LPOSMOY, stocke les positions moyennes dans posmoyx et les ecrit a la fin avec les energies moyennes'
    end if
    if (lprteattotm.EQV..true.) then
       lprteat=.true.
       write(6,*)'LPRTEATTOTM calcule les energies moyenne de chaque atome et les ecrit en retranchant eatref en eV (=0 par defaut)'
    end if


    deltax=deltax*A2cm
    if (dmtype==5.and.deltax.le.0) then
       write(6,*) rang,'dmtype 5 deltax 0'
       call arret_ndm
    endif
    if (dmtype==7.and.deltax.le.0) then
       if (.not.lEev) then
          write(6,*)  'PHONDY: lEev should be set on .true.'
          write(6,*)  'PHONDY: Change accordingly and try again!'
          call arret_ndm
       end if
       write(6,*) rang,'For dmtype 7 deltax must be deltax > 0'
       write(6,*) rang,'Change deltax!'
       call arret_ndm
    endif

    if (dmtype==7) then
       write(6,*)'dmtype==7 is deprecated'
       call arret_ndm
    end if

    if (ipotentiel==20) then
       ldemitab=.FALSE.



    end if



    iThermo=0
    if(lTandersen) iThermo=iThermo+1
    IF ( lTBerendsen ) iThermo=iThermo+1
    IF ( lTcon ) iThermo=iThermo+1
    IF ( lTNose ) iThermo=iThermo+1
    IF ( lTHoover ) iThermo=iThermo+1
    IF ( iThermo .GE. 1) THEN
       if (text.le.0) then
          IF (RANG==0) WRITE(6,*) 'T Const et Text<=0 : stop '
          call arret_ndm
       end if

       if (rang==0)then
          if(lTcon)                WRITE(6,*) 'TEMPERATURE CONSTANTE A LA :','lTcon  '
          if (lTberendsen)         WRITE(6,*) 'TEMPERATURE CONSTANTE A LA :','lTBerendsen  '
          if(lTNose)               WRITE(6,*) 'TEMPERATURE CONSTANTE A LA :','lTNose '
          if(lTHoover)             WRITE(6,*) 'TEMPERATURE CONSTANTE A LA :','lTHoover '
          if(lTandersen)           WRITE(6,*) 'TEMPERATURE CONSTANTE A LA :','lTandersen '
       end IF
       IF ( iThermo .GT. 1) THEN
          IF (RANG==0) WRITE(6,*) 'Thermostat "normal", de Berendsen, de Nose, de Hoover ou de Andersen: en choisir un seul'
          STOP
       end IF
    END IF

    if (linstantrdf .and. iterdf==0) then
       if (rang==0) write (6, *) rang,'linstantrdf et iterdf incompatibles'
       call arret_ndm
    endif

    if (ltabvois .and. rvois==0.0) then
       if (rang==0) write (6, *) rang,'erreur rvois ', rvois
       call arret_ndm
    endif
    if ((.not.ltabvois) .and. rvois>0.0) then
       if (rang==0) write (6, *) 'ltabvois=false et rvois >0 ', rvois
       call arret_ndm
    endif
    rvois=rvois*1.0d-8





    if (lcasca) then
       lspecialinit=.true.
       lax=.true.
       lfilm=.true.
    end if

    if(lTcon.and.dmtype>1) then
       if (rang==0) write(6,*)rang,'lTcon dmtype>1'
       call arret_ndm
    end if

    if(lTcon.and.lcasca) then
       if (rang==0) write(6,*)rang,'lTcon lacasca'
       call arret_ndm
    end if

    if(lTcon.and.(Text<0.))then
       if (rang==0) write(6,*)rang,'Text <0'
       call arret_ndm
    endif


    if((lTberendsen).and.( (dmtype.EQ.21).OR.(dmtype.EQ.22).OR.(dmtype.EQ.3).OR.(dmtype.EQ.30)&
         &.OR.(dmtype.EQ.31).OR.(dmtype.EQ.32).OR.(dmtype.EQ.34).OR.(dmtype.EQ.33).OR.(dmtype.EQ.35) )) then
       write(6,*) 'Berendsen pas possible';stop
    end if

    if (ipotentiel==-1) then
       tpot=.false.
       lpt:     do i=1,npotmax
          if (lpotentiel(i).EQV..true.)then
             tpot=.true.
             exit lpt
          end if
       end do lpt
       if (.not.tpot) then
          if (rang==0) write(6,*)'probleme ipotentiel npotentiel',ipotentiel,lpotentiel
          call arret_ndm
       end if
    end if
    !  if ((all(lpotentiel)==.false.).and.(ipotentiel==-1)) then
    !  end if
    if ((ipotentiel.gt.0).and.(any(lpotentiel))) then
       write(6,*)'choose iptentiel or lpotentiel, not both'
       call arret_ndm
       stop
    end if
    if (ipotentiel.ge.0) lpotentiel(ipotentiel)=.true.
    npotentiel=0
    do ipotcont=0,npotmax
       if (lpotentiel(ipotcont).EQV..true.) npotentiel =npotentiel+1
    end do
    if ((lpotentiel(0).EQV..true.).and.(npotentiel.gt.1)) then
       if (rang==0) write(6,*)'npotentiel>1 et lpotentiel(0)=T'
       call arret_ndm
    end if

    if ((ntyp==-1).and.(npotentiel.gt.1))then
       if (rang==0) write(6,*)'npotentiel>1 et ntyp=-1'
       call arret_ndm
    end if
!    if ((npotentiel.gt.1).and.(lpotentiel(10).eqv..true.)) then
       !if (rang==0)write(6,*)'**** npotentiel >1 ET EAM ==> EAM TAB only!'
!    end if


    if ((lpotentiel(12).EQV..true.).and.(ltabvois.EQV..true.))ldemitab=.false.
    if (ipotentiel.le.-10)then
       if ((rang==0).and.(ltabvois)) write(6,*)'LAMMPS +ltabvois ; impossible pour l instant, ltabvois à false'
       ltabvois=.false.
    end if

    if(lrestart.and.lcorrelvp) then
       if (rang==0) write(6,*)rang,'pas de restart et de correlation'
       call arret_ndm
    end if

    if(lcalcjq) then
       !     fnamjqbis = fnam(1:lenfnam)//'.E_xp'
       fnamjq = fnam(1:lenfnam)//'.jq'
       if(rang==0)open(unit=65,file=fnamjq)
       !         open(unit=66,file=fnamjqbis)
    end if

    if (itetimestep>0)  then
       if ((dmtype.eq.1).or.(dmtype.eq.21).or.(dmtype.eq.22).or.(dmtype.eq.4).or.(dmtype.eq.41).or.(dmtype.eq.42)) then
          if (rang==0) write(6,*) 'The time step changed each', itetimestep,' steps'
       else 
          if (rang==0) write(6,*)'itetimestep seulement avec dmtype =1, 2  or 4 '
          if (rang==0) write(6,*) 'STOP in readdm'
          call arret_ndm
       end if
    end if

    if(lLangevin.and.(Text.le.0.0)) then
       if (rang==0) write(6,*)'Langevin avec Text pas defini : stop'
       call arret_ndm
    end if
    if(llangevin) then
       if (lpr) then
          dmtype=88
       else
          ltest=.false.
          if ((dmtype==4).or.(dmtype==41).or.(dmtype==42))ltest=.true.
          if (.not.ltest)then
             write(6,*)'llangevin only with dmtype =4,41, 42'
             call arret_ndm
          end if
       end if
    end if

    ! end check
    if (lpconxyz) then
       if (.NOT.lprahman) then
          if (rang==0) then
             write(6,*) 'lpconxyz can be USEd ONLY is with PR dynamics or lpr=.true'
             write(6,*) 'stop in <readdm>'
          end if
          call arret_ndm
       end if
    end if


    if (lprahman) then
       itesigma=1
       iteprtsigma=itetemp
       select case(dmtype)
       case(21)
          dmtype=22
       case(22)
          if (wboxf==1) wboxf=0.2
          if (sigstop.le.0) sigstop =0.05 ! critere de conv. sur les contraintes par direction UNITE = kbar

       case(24)
          if (wboxf==1) wboxf=0.2
          if (sigstop.le.0) sigstop =0.05 ! critere de conv. sur les contraintes par direction UNITE = kbar
       case (3,30,31,32,33,34,35)
          lEev=.true.
          if (sigstop.le.0) sigstop =0.05 ! critere de conv. sur les contraintes par direction UNITE = kbar
       end select

       if ((pext.ne.0.).or.(any(sigext.ne.0))) then
          if (rang==0) write(6,*)'SIGEXT en cgs', sigext
          if (rang==0) write(6,*)'Pext en cgs',pext
          do ic=1,3
             sigext(ic,ic)=sigext(ic,ic)+pext
          end do
          if (rang==0) write(6,*) 'transforme en '
          if (rang==0) write(6,*)'SIGEXT', sigext
       end if
       !=== Modif Emmanuel Clouet ================
       ! Verifie si un etat de reference a ete donne
       IF (Sum(h0(1:3,1:3)**2).GE.1.d-30) lUcell=.TRUE.
       ! Transformation A => cm pour le repere de reference
       Pext = (sigext(1,1)+sigext(2,2)+sigext(3,3))/3.d0
       !=== Fin des modifications ================
       h0(1:3,1:3) = 1e-8*h0(1:3,1:3)
       if (lpconxyz) then
          ihbox0(:,:)=0   ! ALL the dimension are blockef except ...
          ihbox0(1,1)=1   ! X ...
          ihbox0(2,2)=1   ! Y ...
          ihbox0(3,3)=1   ! and Z.
       end if
       if ((lpconx).or.(lpcony).or.(lpconz)) then
          ihbox0(:,:)=0   ! ALL the dimension are blockef except ...
          if (lpconx) ihbox0(1,1)=1   ! X ...
          if (lpcony)ihbox0(2,2)=1   ! Y ...
          if (lpconz)  ihbox0(3,3)=1   ! and Z.

       end if
       if (lpcube) then
          ihbox0(:,:)=0   ! ALL the dimension are blockef except ...
          ihbox0(1,1)=1   ! X ...
          ihbox0(2,2)=1   ! Y ...
          ihbox0(3,3)=1   ! and Z.
       end if
       
       do ic=1,3
          do ic2=1,3
             if ((ihbox0(ic,ic2).ne.0).and.(ihbox0(ic,ic2).ne.1))then
                if (rang==0) write(6,*)'non zero ihbox0(',ic,ic2,ihbox0(ic,ic2)
                call arret_ndm
             end if
          end do
       end do
       do ic=1,3
          if ((ihbox0(ic,ic)==1).and.(ipbc(ic).ne.1))then
             write(6,*)'direction ',ic,' lprahman and no pbc ipbc =',ipbc(ic)
             call arret_ndm
          end if
       end do
    else
       iteprtsigma=1
       itesigma=1
    end if
    ! read for cascade
    select case (ibrake)
    case(3)
       if(rang==0) then 
          write(6,*)'electronic stopping with constant coeff read in  elstop.in for Ec> Ecelec  >' , Ecelec
       end if
    case(0)
       if ((tcelec.gt.0).or.(ecelec.gt.0)) then
          write(6,*) 'tcelec > 0 et pas de pertes electroniques : stop'
          call arret_ndm
       end if
    case(1)
       if(rang==0) then 
          write(6,*)'electronic stopping according to elstop.in from SRIM POUR DES EC >' , Ecelec,'!!!!!!!!!!'
          write(6,*)'electronic stopping according to elstop.in from SRIM POUR DES Tempcel  >' , tcelec,'!!!!!!!!!!'

       endif
    case(2)
       if (ecelec.le.0) then
          write(6,*) 'Ecelec <  0 et  perte electronique Langevin : stop'
          call arret_ndm
       end if
       if (Text.le.0) then
          write(6,*) 'Text <  0 et  perte electronique Langevin : stop'
          call arret_ndm
       end if
       if(rang==0) then 
          write(6,*)'electronic stopping according to elstop.in from SRIM ET CONNECTION A LANGEVIN pour EC >' , Ecelec,'!!!!!!!!!!'
       endif
       llangevin=.true.
       gamlg=-1.
    end select
    if ((ibrake.gt.0).and.(l2t.eqv..true.)) then
       if (rang==0) write(6,*)'L2T AND ibrake = 0 STOP'
       call arret_ndm
    end if
    



    if (iterdf >= 0) then

       if (iterdf > 0) nrdf = 0
       if (iteangle > 0) nfda = 0
    endif

    !  if (itedepla.gt.0) lfilm=.true.
    if (lfilm) then
       if((rang==0).and.(lcasca))         open(unit=lufilmpaf, file='filmpaf', status='unknown')
    end if
    

    if(dilat(1).ne.0.0) then
       if(igen.lt.1) then
          if (rang==0) write (6,*) rang,'dilat<>0 et igen<>1 stop'
          call arret_ndm
       end if
       if (dilat(2)==0.0) dilat(2)=dilat(1)
       if (dilat(3)==0.0) dilat(3)=dilat(1)
       if (rang==0) write(6,*)'dilat 1 2 3 ', dilat(1),dilat(2),dilat(3)
    end if

    ! MPI
    if (itesauvinter.gt.0) then
       if (mod(itesauvinter,itesauv).ne.0) then
          write(6,*)'itesauvinter n est pas un multiple de intesauv : stop'
          call arret_ndm
       end if
    end if
    if (rang==0) write (6, *)
    if (rang==0) write (6, '(a,I2)') ' -------- caracteristiques du run DM--------', dmtype
    select case (dmtype)
    case(111)
       if (rang==0) write (6,'(a)') '|=========       ONE STEP           ===============|'

    case (1)
       if (rang==0) write (6, '(a)') '     DYNAMIQUE MOLECULAIRE VERLET STANDARD'
    case (21)
       if (rang==0) write (6, '(a)') '     FAST Quenching Vcst'
    case (22)
       if (rang==0) write (6, '(a)') '     FAST Quenching Pcst'
       lprahman=.true.
    case (23)
       if (rang==0) write (6, '(a)') '     FIRE Quenching Vcst'
    case (24)
       if (rang==0) write (6, '(a)') '     FIRE Quenching Pcst'
       lprahman=.true.
    case (31)
       if (rang==0) write (6, '(a)') '     VIEUX GRADIENT CONJUGUE par défaut = 31 sur les coordonnees cartésiennes '
    case (30)
       if (rang==0) write (6, '(a)') '     VIEUX GRADIENT CONJUGUE sur les coordonnees REDUITES'
    case (3)
       if (rang==0) write (6, '(a)') '     GRADIENT CONJUGUE STANDARD '
       dmtype=33
    case (32)
       if (rang==0) write (6, '(a)') '     STEEPEST DESCENT'
    case (33)
       if (rang==0) write (6, '(a)') '     GRADIENT CONJUGUE STANDARD '
    case (34)
       if (rang==0) write (6, '(a)') '     GRADIENT CONJUGUE FLETCHER-REEVES'
    case (35)
       if (rang==0) write (6, '(a)') '     ADAM relaxation D. Kingma and J. Ba, 2015'
       if(itmax==-1)itmax=300
    case (4)
       if (rang==0) write (6,'(a)') '      DYNAMIQUE MOLECULAIRE VELOCITY VERLET'
       if (llangevin.and.(rang==0)) write (6,'(a)') '      Tcst LANGEVIN'
    case (41)
       if (rang==0) write (6,'(a)') '      MOLECULAR DYNAMICS ADAPTATIVE RESTRAINED PARTICLE SIMULATION with restrained forces'
    case (42)
       if (rang==0) write (6,'(a)') '      DYNAMIQUE MOLECULAIRE ADAPTATIVE RESTRAINED PARTICLE SIMULATION with coplete forces'
    case (5)
       if (rang==0) write (6,'(a)') '      TEST DES FORCES '
    case (8)
       if (rang==0) write (6,'(a)') '      PARRINELLO RAHMAN AUTOCOHERENT '
       lprahman=.true.
    case (88)
       if (rang==0) write (6,'(a)') '      PCst en LANGEVIN'
       lprahman=.true.
    case (6)
       if (rang==0) write (6,'(a)') '      ANALYSE DES POSITIONS EN FIN DE CASCADE '
    case (9)
       if (rang==0) write (6,'(a)') '      DRAG OR NEB DYNAMICS ' 
       itesauvposition=-1
       itesauvforce=-1
       itetemp=-1;itesigma=-1
    case (11)
       if (rang==0) write (6,'(a)') '      UN CALCUL DE FORCES '
       if (rang==0) write (6,*)
    case (19)
       if (rang==0) then
          write (6,'(a)') '      FORCE Matrix calculation '
          if (ltabvois) write (6,'(a)') '  BE SURE THAT RVOIS>RUE+DECAL'
       end if
#ifndef MKL
       if (rang==0)then
          write(6,*)"dmtype=19 works with lapack or MKL"
          write(6,*)"these libraries are NOT linked by default"
          write(6,*)"the force matrix will be written to binary file"
       end if
       lwfm=.true.
!       call arret_ndm
#endif

       decal=decal*1d-8
       if((ndecal.le.0).or.(decal.le.0)) then 
          if (rang==0) write (6,*)' problem decal, ndecal:',decal,ndecal
          call arret_ndm
       end if
       if (rang==0) write (6,*)' decal, ndecal lparaFM, naparaFM:',decal,ndecal,lparafm,nparafm
       if (rang==0) write (6,*)
    case (15,151)
       if ((rang==0).and.(dmtype==151)) write (6,'(a)') '      CALCUL MONTE CARLO GRAND CANONIQUE '
       if ((rang==0).and.(dmtype==15)) write (6,'(a)') '      CALCUL MONTE CARLO DES CHEMINS '
       if (rang==0) write (6,*)'LPARAPATH NPARAPATH', lparapath, nparapath
!!$       if ((nparapath.gt.1).and.(.not.lparapath)) then
!!$          write(6,*)'nparapath >1, needs lparapath = TRUE'
!!$          call arret_ndm
!!$       end if

       if (rang==0) write (6,*)
       if ((lparapath).and.(nparapath.le.1)) then
          write(6,*)'lparapath ET nparapath=1 stop'
          call arret_ndm
       end if
!       if (dmtype.ne.12) then
          if ((.not.lrestartmcgc) .and. (.not.((idirectionmcgc==0).or.(idirectionmcgc==1)))) then
             write(6,*)'set idirectionmcgc to 0 or 1 '
             call arret_ndm
          end if
          if (rang==0) then
             write(6,*)'MCGC starts in direction, idirectionmcgc ', idirectionmcgc, " protocol=", protocol_mcc
          end if
!       end if
#ifdef PARA
#else
       if (lparapath) then
          write(6,*)'lparapath=true et sequentiel==> lparapath=.false.'
          lparapath=.false.
       end if
#endif

!#ifdef ART    
    case (12)
       if (rang==0) write (6,'(a)') '|=========NDM ENTERTAINMENTS presents:===============|'
       if (rang==0) write (6,'(a)') '|---------ART nouveau by N MOUSSEAU.---------------|'
       if (rang==0) write (6,'(a)') '|======== colored by Cosmin Marinica!==============|'
       if (rang==0) write (6,'(a)') '|======== updated by J-P Crocombette!==============|'
!#endif
#ifdef SUNDAE    
    case (16)
       if (rang==0) write (6,'(a)') '|=========       NDM + SUNDAE       ===============|'
       if (rang==0) write (6,'(a)') '|---------..........................---------------|'
       if (rang==0) write (6,'(a)') '|==================================================|'
#endif
#ifdef MAB    
    case (17)
       if (rang==0) write (6,'(a)') '|=========       NDM + MAB          ===============|'
       if (rang==0) write (6,'(a)') '|---------..........................---------------|'
       if (rang==0) write (6,'(a)') '|==================================================|'
#endif


    case(112)
       write(6,*)'simple test de distance entre atomes'
    case(113)
       write(6,*)'recvherche de la position la plus éloignée des atomes'
    case(16)
       if (rang==0) write (6, *) 'CALCUL BABAR,',dmtype
    case default
       if (rang==0) write (6, *) 'mauvais type de calcul dmtype=TTT',dmtype
       call arret_ndm
    end select

    if (any(ihbox0==0))then
       if (rang==0) then
          write(6,*)'incomplete cell relaxation'
          do ic=1,3
             do ic2=1,3
                if (ihbox0(ic,ic2)==1)then
                   write(6,'(A,2I2,A,I2)')' ihbox0(',ic,ic2,')=',ihbox0(ic,ic2)
                end if
             end do
          end do
       end if
    end if


    if (lcontr)  write (6, '(a)') '******************* CONTRAINTE !!! *****'


    if (lTcon) then
       if (rang==0) write (6, *) 'TEMPERATURE CONSTANTE a la main Text= ',text
    endif
    if (lTberendsen) then
       if (text.le.0) then
          write(6,*)'text<0' ;stop
       endif
       if (rang==0) write (6, *) 'TEMPERATURE CONSTANTE a la Berendsen Text= ',text
    endif
    if ((text.gt.0).and.(.not.((dmtype==15).or.(dmtype==151)))) then
              if (.not.(ltberendsen.or.llangevin.or.lThoover.or.lTnose)) then
          write(6,*)'text<0 mais pas dalgo',dmtype ;stop
       end if
    end if
    select case (igen)
    case (-1)
       if (rang==0) write (6, *) 'generation du crystal'
    case (0)
       if (rang==0) write (6, *) 'generation du crystal ; puis run'
    case (1)
       if (rang==0) write (6, *) 'run a partir du fichier .cin'
    case (2)
       if (rang==0) write (6, *) 'écriture de gin à partir du fichier .cin'
    case (3)
       if (rang==0) write (6, *) 'modification du fichier .cin'
    case default
       if (rang==0) write (6, *) 'mauvais igen=', igen
       call arret_ndm
    end select

    if (ltabvois) then
       if (npotentiel.gt.1) then
          if (rang==0) write(6,*)'ltabvois avec plusieurs potentiels= pas programmee (demi table ou table complete = prise de tete'
          call arret_ndm
       end if

       !     if (tempstopcel.gt.0) ltpcel=.true.
       if (ltpcel) write (6, *) '   -> -> pas de contrainte par celulles'

       if (rang==0) write(6,*)'IPOTENTIEL',ipotentiel
       select case (ipotentiel)
       case(:9)
          ldemitab=.TRUE.
          if (rang.eq.0) write (6, *) '    DEMI-TABLE DES VOISINS rvois ',rvois
       case(11:17)
          ldemitab=.false.
          if (rang.eq.0) write (6, *) '    DEMI-TABLE DES VOISINS rvois ',rvois

       case(10)  ! Potentiel EAM
          if(dmtype==7) then 
             ldemitab=.FALSE.
             if (rang.eq.0) write(6,*)'    TABLE DES VOISINS COMPLETE rvois ',rvois
          else
             ldemitab=.TRUE.
             if (rang.eq.0) write (6, *) '    DEMI-TABLE DES VOISINS rvois ',rvois
          end if
       case(20)
          ltabvois=.true.
          ldemitab=.false.
          !itetabvois=0

       end select

       !     if(ipotentiel.le.10) then 
       !        ldemitab=.TRUE.
       !        if (rang.eq.0) write (6, *) '    DEMI-TABLE DES VOISINS rvois ',rvois
       !     else
       !        ldemitab=.false.
       !        if (rang.eq.0) write(6,*)'    TABLE DES VOISINS COMPLETE rvois ',rvois
       !     end if

    endif


    if ( (dmtype==21).or.(dmtype==22).or.(dmtype==3).or.(dmtype==30).or.(dmtype==32)&
         &.or.(dmtype==33).or.(dmtype==31).or.(dmtype==34).or.(dmtype==35).or.(dmtype==9)&
          ) then    
       if ( (fpstop<0).and.(fsumstop<0)) then
          if (rang==0) write(6,*) 'One of fpstop and fsumstop must be positive for dmtype=',dmtype
          if (rang==0) write(6,*) 'STOP in readdm',fpstop,fsumstop
          call arret_ndm
       end if
       if ( (fpstop < 0) .and. (dmtype==9) ) then
          if (rang==0) write(6,*) 'NEB and DRAG implementation only for positive fpstop'
          if (rang==0) write(6,*) 'STOP in readdm'
          call arret_ndm 
       end if
       if  ( (fpstop>0).and.(fsumstop>0) ) then
          if (rang==0) write(6,*) 'DANGER - WARNING - ACHTUNG:  both fpstop and fsumstop are positive !!!'
       end if
       if ((dmtype==8).and.(sigstop.LT.0)) then
          write(6,*)'trempe PR + sigstop <0 ; stop'; stop
       end if
    end if

    if  (dmtype==9) then
       if (lperiod) then
          if (rang==0)  write(6,*) 'There is no NEB and DRAG implementation for lperiod true'
          if (rang==0)  write(6,*) 'put your lperiod to false in din file and restart.'
          call arret_ndm
       end if
    end if

    if ( (.not.lperiod).and.(itesauvposition>0).and.lsuivinonpbc) then  
       if (rang==0) write(6,*)' lsuivinonpbc will be turn to FALSE'
       lsuivinonpbc=.false.
    end if
    if (lsuivinonpbc) then

       if ((dmtype.ne.4).and.(dmtype.ne.41).and.(dmtype.ne.42)) then
          if (rang==0) write(6,*) 'lsuivinonpbc is implemented only with velocity verlet'
          if (rang==0) write(6,*) 'Or dmtype=4. Change and restart until there I will stop for you.'
          call arret_ndm 
       end if

       if (itesauvposition<=0) then
          if (rang==0) write(6,*) 'itesauvpostion MUST be positive if you want lsuivinonpbc TRUE'
          if (rang==0) write(6,*) 'STOP in readdm'
          call arret_ndm
       end if
    end if
    if (lforcetabulate) then
       if (ipotentiel/=10) then
          write(*,*) 'There is no implementation for lforcetabulate TRUE and ipotentiel ', ipotentiel
          write(*,*) 'Change lforcetabulate of FALSE or ipotential to EAM (10) '
          write(*,*) 'stop in readdm'
          call arret_ndm
       end if
       if (lcasca) then
          write(*,*) 'There is no implementation for lforcetabulate TRUE and lcasc TRUE'
          write(*,*) 'Change lforcetabulate of FALSE or lcasc on FALSE'
          write(*,*) 'stop in readdm'
          call arret_ndm
       end if
    end if

    if (itetemp2==-1) itetemp2=itetemp
    if (rang==0) write(6,*)
    if (rang==0) write (6, *) '     ANALYSES '
    if (rang==0) write (6, *) 'itetemp=', itetemp, ' iteprtsigma=', iteprtsigma
    if (itecoordo>0)  write(6,*)  ' itecoordo=', itecoordo
    if (itedepla>0) then
       lax=.true.
       if (rang==0) write (6, '(A,I3,A,D9.3,A,D9.3,A,I3,A,I3)') ' itedepla=', itedepla, &
            ' tdepla=', tdepla*1D+8, ' tdepla2=', tdepla2*1D+8, ' itesauv=', &
            itesauv, ' itesauvposition=', itesauvposition
    end if
    if (lfilm) then
       if (rang==0) write (6, '(A,D11.3)') ' film; seuil=', tdepla*1D+8
       if (rang==0) write (6, '(A,D11.3)') ' film; seuil2=', tdepla2*1D+8
    endif
    if (rang==0) write(6,*)
    if (rang==0) write (6, *) '     CONTROLES '
    if (rang==0) write (6, '(A,D11.3)') 'tstep=', tstep
    if (rang==0) write (6, *) 'itmax=', itmax, 'timemax= ',timemax,' itab=', itab, ' itetimestep=', &
         itetimestep
    if (itederive>0)  write(6,*) ' itederive=', itederive
    if (rang==0) write (6, '(A,F10.1,A,F10.1,A,F10.1,A,F10.1)') 'tinit=', tinit
    if (tempdeplainit.GT.0) then
       if (debyetemp==-1) then
          write (6,*)"debyetemp  non definie mais tempdeplainit> 0" ; stop 
       end if
    end if
    if (rang==0) write (6, '(A,F10.1)') 'tempdeplainit=', tempdeplainit
    if (rang==0) write (6, '(A,F10.1)') 'debyetemp=', debyetemp



    if ((tempstop>0.).and.(rang==0)) write(6,*)' tempstop=', tempstop
    if ((tfcou>0.).and.(rang==0))         write (6, '(A,F10.1,A,F10.1,A,F10.1)') 'tfcou=', tfcou, ' epcou=', &
         epcou*1D+8


    if (rang==0) write (6, *) ' ----------------------------------'
    if (rang==0) write (6, *)
    if (rang==0) write (6, *)



    if(lTandersen.and.rang==0) write(6,*)'Tandersen nuandersen = ',nuandersen

    if ((lprtsigat.eqv..true.).or.(lsigatcel.eqv..true.))then 
       lsigat=.true.
    else
       lsigat=.false.
    end if

    if(lPrtSigat.and.(.not.ltabvois)) then
       write(6,'(a)')rang,'contrainte atomique programme en table des voisins&
            & avec un potentiel EAM ou un terme a deux corps seulement'
       call arret_ndm
    end if


    !  if(itetemp2.gt.0)then
    !     fnamdin = fnam(1:lenfnam)//'.2.T'
    !     open(unit=112, file=fnamdin, status='unknown')
    !  end if


    if ( ( (dmtype==3).OR.(dmtype==30).or.(dmtype==32).or.(dmtype==34).or.(dmtype==35)&
         &.or.(dmtype==33).or.(dmtype==31) ) &
         .and.(fpstop.le.0.0).and.(fsumstop.le.0.0)) then
       if (rang==0) write(6,*) rang,'critere de conv. sur la force par atome max negative' 
       if (rang==0) write(6,*) rang,'fpstop', fpstop
       if (rang==0) write(6,*) rang,'critere de conv. sur la force sqrt ( sum_f F_i^2 ) negative' 
       if (rang==0) write(6,*) rang,'fsumstop', fsumstop
       if (rang==0) write(6,*) rang,'un de deux doit etre > 0. Exemple:'
       if (rang==0) write(6,*) rang,'fpstop = 0.05, fsumstop=0.1 les unites sont eV/A'      
       call arret_ndm
    end if


    ! Gestion des atomes bloques
    if (rulayer.gt.0.0)then
       rulayer=rulayer*1.0d-8

       IF (rang==0)write(6,*)'atomes immobiles fixes par rulayer ', rulayer*1d8
    end if


    !  if (rang==0) write(6,*) 'sortie readdm'



    if ((iteanapos.eq.-1).and.(itecompcr.ne.-1)) iteanapos=itecompcr
    usdh = 1/(two*tstep)
    if (rang==0)then 
       write(6,*)'nb de potentiels', npotentiel
       if (npotentiel==1) then
          write(6,*)'ipotentiel',ipotentiel
          !        write(6,*)lpotentiel
       else
!!$          write(6,*)'npotentiel buggué stop'
!!$          call arret_ndm

          do ipotcont=1,npotmax
             if (lpotentiel(ipotcont).EQV..true.)write(6,*)'potentiel actif', ipotcont
          end do
          if (lcasca.eqv..true.) then
             if (rang==0) write(6,*)'ATTENTION!!! npotentiel>1 et ziegler surement faux !!!!'
             call arret_ndm
          end if
       end if
    end if
    if (rang==0) write(6,*)'fmt_cin',fmt_cin

    if(lPkbar) then
       unitP=1.0d-9
       cunitP='kbar'
    else
       unitP=1.0
       cunitP='d/cm2'
    endif
    if(lEev) then
       unitE=erg2eV
       cunitE='  eV'
    else
       unitE=1.0
       cunitE=' erg'
    end if
#ifdef LAMMPS_VERSION
    if (ipotentiel.lt.0) then 
    if(trim(units_lammps)=='metal') then
       energy_conversion_lammps=1/erg2ev
       position_conversion_lammps=A2cm
       pressure_conversion_lammps=1d6
    elseif(trim(units_lammps)=='real') then
       energy_conversion_lammps=0.043/erg2ev
       position_conversion_lammps=A2cm
       pressure_conversion_lammps=1013250.0
    elseif(trim(units_lammps)=='si') then
       energy_conversion_lammps=1e7
       position_conversion_lammps=1d2
       pressure_conversion_lammps=10.0
    elseif(trim(units_lammps)=='cgs') then
       energy_conversion_lammps=1
       position_conversion_lammps=1
       pressure_conversion_lammps=1.0
    elseif(trim(units_lammps)=='electron') then
       energy_conversion_lammps=27.211399/erg2ev
       position_conversion_lammps=A2cm*0.529177249
       pressure_conversion_lammps=10.
    else
       write(6,*)'error in units_lammps',units_lammps
       call arret_ndm
    end if
 end if
#endif     

    !condition d'arret du prog si l'utilisateur veut utilise la methode mcgc mais n'a pas defini le pas lambda pour l'integration de la particule    
    if((dmtype == 15).or.(dmtype==151))then
       if (pas_lambda_mc.lt.0) then
          write(6,*)'To use  MONTE-CARLO, indicate pas_lambda (number of insertion steps'
          call arret_ndm
       end if
       !idem dans le cas ou l'utilisatuer utilise le biais sur les retraits sans avoir defini la fonction alpha (fermi dirac) pilotant celui ci
       if(((dmtype == 15).or.(dmtype==151)) .and. (lbiais_retrait).and. (fdmc_1 .eq. -1000.0) .and. (fdmc_2 .eq. -1000.0)) then
          write(6,*)'Pour utiliser la methode MCGC avec le biais sur les retraits: indiquer les param pour le fermidirac'
          call arret_ndm
       end if
       if(((dmtype == 15).or.(dmtype==151)) .and. (lbiais_retrait).and. (nbatplus.lt.1))then 
          write(6,*)' methode MCGC avec  le biais sur les retraits: pas possible aev nbatplus>1'
          call arret_ndm
       end if
       if(((lbiais_retrait).or.(lbiais_inser)).and.(lspring)) then
          write(6,*)'SPRING OR BIAS not both....'
       end if
       if (lspring) then
          if (rang==0) write(6,*) 'vanishing spring of strenght ', k_spring
          k_spring=k_spring*1d16/erg2ev
          if ((ins_typ.ne.55).and.(ins_typ.ne.11).and.(ins_typ.ne.33).and.(ins_typ.ne.44)) then
             if (rang==0) write(6,*) 'lspring==> ins_typ=11 (site) or 33 (plane)  or 44 (line)  or 55 (bubble) ', ins_typ
             call arret_ndm
          end if
          
       end if
       select case(ins_typ)
       case(0)
          if (rang==0) write(6,*)' MCC N-> N+1 in all the box'
       case(1,3,33,44,55,11)
          
          
             select case (ins_typ)
             case(1)
 if (rang==0)      write(6,*)'MCC N-> N+1 in a sphere',r0mcgc,bublcenter
             case(11)
                if (rang==0) write(6,*)'MCC N-> N+1 in a site with spring pos/spring ',bublcenter,kspring
             case(3)
                do ic=1,3
                   if (ic==izlins) cycle
                   if ((zlcenter(ic).ne.0).or.(izlins==0)) then
  if (rang==0)         write(6,*)'ins_typ,izlins zlcenter inconsitstency 3', ins_typ,izlins,zlcenter(:)
                      call arret_ndm
                   end if
                end do
  if (rang==0)    write(6,*)'MCC N-> N+1 in a slice ',r0mcgc,izlins,zlcenter(izlins)
             case(33)
                do ic=1,3
                   if (ic==izlins) cycle
                   if ((zlcenter(ic).ne.0).or.(izlins==0)) then
       if (rang==0)  write(6,*)'ins_typ,izlins zlcenter inconsitstency 33', ins_typ,izlins,zlcenter(:)
                      call arret_ndm
                   end if
                end do
 if (rang==0)          write(6,*)'MCC N-> N+1 in a plane with spring norm/pos/spring ',izlins,zlcenter(izlins),kspring
             case(44)
                if ((zlcenter(izlins)==0).or.(izlins==0)) then
      if (rang==0)     write(6,*)'ins_typ,izlins zlcenter inconsitstency 44 ', ins_typ,izlins,zlcenter(:)
                   call arret_ndm
                end if
 if (rang==0)   write(6,*)'MCC N-> N+1 in a plane with spring norm/pos/spring ',izlins,zlcenter(izlins),kspring
             case (55)
                if (rang==0)      write(6,*)'MCC N-> N+1 in a sphere with a spring',r0mcgc,bublcenter, kspring
             end select
          
          if ((R0mcgc.lt.0).and.(lspring.eqv..false.)) then
             if (rang==0) write(6,*)' R0mcgc.lt.0'
             call arret_ndm
          end if
          R0mcgc=R0mcgc*1d-8
          if ((ins_typ==3).or.(ins_typ==33)) then
             if (izlins==-1) then
                if (rang==0) write(6,*)' izlins 1,2 or 3 ?'
                call arret_ndm
             end if
          end if
       case(2)
          idirectionmcgc=0
          if (rang==0) write(6,*)' MCC semi grand canonique'
          if ((typswitch1==0).or.(typswitch2==0).or.(typswitch1==typswitch2))then
             if (rang==0) write(6,*)'problem with typswitch'
             call arret_ndm
          end if
       case default
          if (rang==0) write(6,*)' ins_typ =0 or 1, 2'
          call arret_ndm
       end select
    end if
    if (lcdp) then
       select case(dmtype)
       case(1,4,8)
          itetimestep=1
       case(2,3,32,33,34,35,21,22,23,24)
          itetimestep=-1
       case default
          if (rang==0) write(6,*)'dmtype inconsistent with creaDP', dmtype
          call arret_ndm
       end select
    end if

    if (fstpdecr.le.1) then
       write(6,*)'fstpdecr must be >1 ; stop'
       stop
    endif
    If (Tinitbox==-1) then
       if (Tinit.GT.0) then
          Tinitbox=Tinit
       else if (Text.GT.0) then
          Tinitbox=Text
       else
          Tinitbox=0.
       end if
    end If
    lseedcom=.false.

    if (iseed.le.0) then
       call system_clock (iseed)

       iseed =iseed +10*rang
    else
       lseecom=.true.
    end if
    if (lspacendm.eqv..false.) then
#ifdef PARA
       call mpi_world%bcast(0,iseed)
       lseedcom=.true.
#endif
    end if
    if (lseedcom   ) then
       if (rang==0)    write(6,*)'all ranks  readdm iseed ',iseed
    else
       write(6,*)'rang readdm iseed ',rang,iseed
    end if
    
    return
456 print *,'Erreur lors de la lecture du fichier .din, verifier l''ajout de fmt_cin'
    if ((dmtype.ge.41).and.(dmtype.le.42)) then
       if ((any(noxyzkmin(:).ge.1)).and.(.not.lpartarps)) then
          write(6,*)'NOT lpartarps and noyzkmin >0 STOP'
          call arret_ndm
       end if
    end if
  end subroutine readdm
  
end module readdm_mod
