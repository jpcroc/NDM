! *****************************************************************
subroutine readdm
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use jqmod
#if(PARA)
  use mod_mpi
#endif
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
  integer :: ludin, lufilm, lufilmpaf,  i,itean, ic, iThermo,itecompcr,ipotcont
  character :: fnamdin*80
  logical :: lginread,ltriclin,lpcon,lfissure,tpot
  logical :: lxFrozen,lyFrozen,lzFrozen, lxyFrozen, lxzFrozen, lyzFrozen, lxyzFrozen
  !  integer :: imFree     ! nb d'atomes libres
  !-----------------------------------------------
  !
  !
  !
  !

  namelist /input/itab, itetabvois, itetemp, itesigma, itefcc, itedepla, tdepla, lfilm, &
       tempstop, tempstopcel,dmtype, lFire, ttol, tfroi, itecoordo, tstep, itetimestep, tsfact, &
       tinit, tcooling, tfcou, epcou, lcasca, lfissure, itmax, itean, itespebcout,  &
       itederive, igen, linstantrdf, iterdf, nrdf,nfda, linstantfda,rclu, itesauv, formatsauv, &
       lrestart, lPathFromGin, tgc, ltabvois, rvois, ltpcel, nox, noy, noz, imm, dfpred, &
       ltranche, rulayer,iterasmol, lpcon, lprtzlm,pext, wbox, wNose, lpcon2, lpconxyz, tbox, &
       iteangle, ipotentiel, lpotentiel, itesauvposition, itesauvforce, lfilmext, tdepla2, &
       lTcon,Text,iteTconst, lTberendsen, lTNose, lTHoover, nHoover, tauTcon, ldecal_bc, ldyn2D, &
       maxorder,  lalea, rsep, &
       h0, sigext,lconstrtot,lEev,lPkbar,deltax,lcorrelvp,lvpread,&
       lcalcjq,dilat,lderive,lTandersen,nuandersen,landerscou,Llangevin,gamlang,ilangevin,&
       lcdp,lsigtyp, ljqbh,lEparat,itebdv,itetemp2,itecompcr,iteanapos,ldislo,epcoudis,&
       fdislo,lnemd,fnemd,fpstop,iseed,fsumstop,sigstop,lcontr,lpr,lUcell,ibordcou,iteplz,nplz,ngrid,lperiod,&
       lprteat,lprteattotm,lprtfat,lprtsigat,itecfg,npath,nebtype,nebrelaxation,maxneb,kspring,deltaRmax,&
       rcangle,rcrdf,deltaestop,nbmoye,lHcyl,fmt_cin,lginread,ltriclin,nvperat, &
       lFrozen,lxFrozen,lyFrozen,lzFrozen,lxyFrozen,lxzFrozen,lyzFrozen,lxyzFrozen,imFree,imFirstFrozen,&
       natperc,iteanaposneb,ntyp,&
       lbulle,ldesinteg,nstepdes,ides, kspr,xpspr,typspr,tempdes,neb_noise,neb_noise_scale,lsuivinonpbc,lposmoy,&
       eatref,lheat,rheat,iteheat,theat,Eheat,HessianOrder,kappa,niteration,lanczos_step,mdcg_noise_scale, &
       mdcg_noise, lforcetabulate,ivisu,ibound,user_strainrate,user_stress_yz,fdbkcoef, decal_bc,&
       tempdeplainit,debyetemp,ibrake,lprtpot,ngrdel,timemax,tpseuils,lrctest


  !
  !   set default values for variables in namelist
  !
  if (rang.eq.0) write(6,*) '>>>>>>>>>>> entree readdm'


  fnamdin = fnam(1:lenfnam)//'.din'
  ! variables de dynamique

  imm = 0                     !dimensionnement des tableaux atomiques
  itab = 10                   !period of cell repartition
  itetabvois = 10             !periode de calcul de la table des voisins
  tempstop = -1.0             !temperature of run stop
  tempstopcel = -1.0             !temperature of run stop
  dmtype = 0  
 !dmtype = type of calculation : 1 -> MD
  !                               2 -> quench (trempe) or fire quench
  !                               3 -> gradient conjugue sur les coordonnes cartesiennes
  !                              30 -> gradient conjugue sur les coordonnes reduites
  !                               4 -> Velocity Verlet 
  !                               5 -> test des forces 
  !                               6 -> analyse des positions en fin de cascade 
  !                               7 -> calcul des phonons
  !                               8 -> PR
  !                               9 -> NEB
  !                              10 -> PARIN RAHMAN 
  !                              11 -> UN SEUL CALCUL DE FORCES
  !                              12 -> ART
  !                              16 -> SUNDAE
  !                              17 -> MAB
  lFire = .true.              ! Fire algorithm is used for quenching (cf tr_fire.F90)
  ttol = 0.0                  !max tolerance for temperature in %
  tfroi = -1.0                !imposed temperature
  tstep = 1.0                 !timestep in 10^-15 sec unit
  itetimestep = -1            !period of check in timestep
  tsfact = 10.0               !change in time step factor
  tinit = -1.0                !initial temperature
  tcooling = -1.0             !cooling rate
  tfcou = -1.0                !temperature of the border of the box
  epcou = -1.0                !width of the border of the box
  lcasca = .FALSE.            !cascade Y/N
  lfissure = .FALSE.          !crack Y/N
  itmax = -1                  !maximum number of iterations
  itederive = -1              !"derive" correction
  igen = -2                 !type de generation :0 a partir de.gin, +1 a partir de .cin; -1 de gin vers cin puis stop +2 modification de cin puis stop
  lrestart = .FALSE.          !if T : restarting from an interrupt job
  lPathFromGin = .FALSE.      !if T : read initial path in gin files *.1.gin, *.2.gin, ... (NEB calculaion)
  tgc = 0.0                   ! threshold for CG calculation
  ltabvois = .FALSE.          ! methode de la table des voisins
  lconstrtot=.FALSE.           !!construction de la table des voisins T=double boucle F=via cel.
  rvois = 0.0                 ! rayon de la table des voisins
  ltpcel = .FALSE.            ! output of temperature and stress in each cell
  lforcetabulate = .FALSE.    ! The derivative of the energy is NOT tabulated. TRUE if it is.


  itesauv = 100               !period for saving
  itesauvposition = 0         !periode pour sauvegarde des positions en binaire
  itesauvforce = 0            !periode pour sauvegarde des forces en binaire
  formatsauv = 3              !format of saving  2 MC triclin; 1 DM triclin
  fmt_cin=1                  !format des fichiers .cin 0 : initiale, 1 = para
  dfpred = 0.1            ! eguess for GC calculations and quenching
  nox = -1
  noy = -1
  noz = -1
  ltranche = .FALSE.          ! existence d'une trache gelee
  rulayer=0.0                 ! largeur de la tranche gelee par 
  ibordcou=0                  !refroidissement sur 3 bords ou seuleument z
  lpr=.false.                 ! parinnelo rahman   contrainte constante
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


  pext = 0.0                  !pression  par defaut
  wbox = 0.0                  ! masse de la boite pour Parrinello-Rahman (par defaut egale a 0.5*masse totale
  wNose = 0.0                 ! masse de la boite pour thermostat de Nose (par defaut egale a wbox)
  tbox = 1000.0               !"temps" de la boite
  lTcon=.false.               !algorithme a temperature constante
  lTberendsen=.false.         ! algorithme a temperature constante
  lTNose=.false.              ! algorithme a temperature constante de Nose
  lTHoover=.false.            ! algorithme a temperature constante de Hoover
  nHoover=1
  tauTcon=200.0               ! The rescales "time" for the Berendsen algorithm
  Text=-1.
  iteTconst =itetemp
  lalea = .FALSE.             ! structure initiale aleatoire
  rsep = 1.0               !Distance de separation pour le tirage aleatoire
  ipotentiel = -1              ! definit type potentiel : 0=Born-Mayer-Huggins, 1=Buckingham, 2=watanabe,3=buck8,4=UO2, 5 terme Morse, 6=SW ÃÂ  la Vashista ; 7 pot paire tabule ; 10 EAM; 12 ZrC JuLi(+Tersoff Doan)  ; 13 Tersoff coupure COS; 14 Tersoff coupure FD ; 15 tersoff coupure SIN (original) ; 11 Ercollesi
  npotentiel = 1              ! nb de potentiels
  lpotentiel(:)=.false.
  ntyp=-1                    ! le nombre de type DOIT etre specifie si le nombre de potentiel est superieur ÃÂ  1
  ! PME
  maxorder=10                                ! Ordre du developpement maximal de la PME
  lvpread=.true.
  dilat(:)=0.0
  lderive=.false.
  lTandersen=.false.       ! temperature constante a la Andersen
  nuandersen=1.0d14        !frequence de tirage aleatoire des vitesses en Hz (valeur elevee = pour cascades)
  landerscou=.false.       ! Andersen seulement sur les bords
  lLangevin=.false.        ! Langevin MD
  gamlang =0.001            ! Gamma de Langevin
  ilangevin=1
  iko=-1
  lcdp=.false.             ! algorithme d'accumulation de defauts ponctuels
  ldislo=.false. ! calcul de dislocation
  epcoudis= 0.0 !epaisseur de la couche avec ajout de force pour dislo
  fdislo=0.0     ! force appliquee aux atomes de bords 
  lnemd=.false.  ! Kth par la methode NEMD Evans, P7229
  fnemd=0
  fpstop =-0.05 ! critere de conv. sur la force par atome max  pour les trempes UNITE = EV/ANG
  sigstop =-0.05 ! critere de conv. sur les contraintes par direction UNITE = kbar
  fsumstop =-0.1 ! critere de conv. sur la force sqrt ( sum_f F_i^2 )  pour les trempes UNITE = EV/ANG
  lcontr=.false. ! dynamique contrainte (routine contrainte)
  iseed=0 ! si <>0 controle le tirage aleatoire des vitesses
  !variables d'analyse

  itetemp = 20                !period of temperature calculation
  itesigma = -1               !period of stress calculation
  itefcc=-1                   ! period of fcc structure analysis
  itedepla = -100              !period of displacement cal.
  tdepla = 1.0                !threshold for displacement
  tdepla2 = -1.0              !second seuil pour calcul des atomes deplaces
  lfilm = .FALSE.             !film making of displaced atoms
  lfilmext = .FALSE.          !film par iteration des atomes deplaces
  itecoordo = -100             !period of coordination calculation
  itean = 0                  !general control for analysis
  iterdf = -1                 !period of RDF calc. : -1 never ; 0 : nrdf last iterations; +iterdf every iterdf iterations
  linstantrdf = .FALSE.       !F: calculates and prints the average of the RDF, T : calculent RDFrdf iterations
  linstantfda = .FALSE.       !F: calculates and prints the average of the RDF, T : calcuatent RDFrdf iterations
  nrdf = 0
  nfda=0

  lprtzlm = .FALSE.           ! plot du nombres d'atomes par tranche suivant z
  lEev=.false.
  lPkbar=.false.
  ! definition des rayons de coupure pour le calcul des coordinences autour de chaque type atomique
  deltax=0.0
  rclu(:)=2.0


  iterasmol = -1                             ! <0 --> genere aucun fichier positions pour logiciel rasmol
  iteangle = -1                              ! pilote creation de fichier positions pour
  ! >=0 debut et fin d'execution

  lufilm = 89
  lufilmpaf = 79
  ludin = 94
  lcorrelvp=.false.
  lcalcjq=.false.
  lsigtyp=.false.             ! calul et affichage de la contrainte atomique
  lEparat=.false.             ! calul et affichage de l'energie par atom
  itebdv=-1  ! frequence de calcul des bond valence
  iteplz=0
  nplz=1000
  itetemp2=-1   ! frequence d'ecriture de la temperature dans fichier separe
  itecompcr=-1  ! remplacee par iteanapos
  iteanapos=-1  ! frequence de comparaison avec cristal de reference

  lperiod=.true.    ! conditions periodiques
  lprteat=.false.   ! if you want to print the energy on atom
  lprtsigat=.false. ! calul et affichage de la contrainte sur chaque atome
  lprteattotm=.false.   ! energie par atome totale (pot+cin) moyenne
  ngrid = 20000  ! taille de la grille des potentiels
  itecfg=-1   ! ecriture de fichiers .cfg pour AtomEye

  !---inNEB
  nebtype=2        ! drag methos is the default     
  nebrelaxation=2  ! We relax all the atoms if nebrelaxation==1 only
  ! the most "deplaced" atoms      
  npath = 15       ! 15 images of the neb is the default
  maxneb = 700     ! the MAX of NEB steps
  kspring = 1.0    ! the default value for the spring
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
  deltaestop=0.0
  nbmoye=100
  lHcyl=.false.
  ltriclin=.true.
  lprtfat=.false.

  lFrozen=.FALSE.
  lxfrozen=.FALSE.             ! .true.: certains atomes sont bloqueÂ¸ (pas de dynamique)
  lyfrozen=.FALSE.             
  lzfrozen=.FALSE.             
  lxyfrozen=.FALSE.             
  lxzfrozen=.FALSE.             
  lyzfrozen=.FALSE.             
  lxyzfrozen=.FALSE.             
  imFree=-1                   ! The index from which all the atoms with the index i >  imFree   are frozen.
                              !                          or with                  i <= imFree   are free
  imFirstFrozen=0             ! The index from which all the atoms with the index i <= imFirstFrozen are frozen
                              !                                         the index i >  imFirstFrozen are free
                              ! imFirstFree can be used in the same time with imFree
 
  nvperat=-1                  ! nb moyen de voisins par atomes
  natperc=-1   
  iteanaposneb=0
  lbulle=.false.
  ldesinteg=.false.    ! calcul du delta F de la desintegration d'un atome
  nstepdes=-1
  ides=1
  kspr=10.0
  xpspr(:)=-1000.
  typspr=0
  tempdes=-1.0
  lsuivinonpbc=.false.  ! enable or disable a copy of non folded positions (by the pbc conditions)  in binary form each itetimestep. 
  lposmoy=.false.       ! writes the average position and energy of the atoms in a .mol file
  eatref(:)=0.

  lheat=.false.
  rheat=0.
  Theat=0.0
  Eheat=0.
  ivisu=1    ! format de sortie dans rasmol.f90 : ivisu=1=.mol, ivisu=2=vsim mal codé, ivisu=3=xred

 ! management of the specific boundary conditions (free or rigid)  ---------------------------   !*!
  ibound = 0	! ( ibound = 0 <=> no spe BoundC, ibound = 1 <=> strain controlled BoundC, ibound = 2 <=> stress controlled BoundC)		!*!
  user_strainrate = 0.	! crystal strainrate (ibound=1)			!*!
  user_stress_yz  = 0.	! stress on the surface (ibound=2)		!*!
  fdbkcoef	  = 0.  ! feedback coefficient for the correction of applied stress (ibound=3)   !*! 
  decal_bc = 0.
  ldecal_bc = .false.
  itespebcout = -1      ! on n'écrit pas de .cfg pour le film
  ldyn2D = .false.      ! par défaut : bords libres selon Y
 !   ----------------------------------------------------------------------------------------    !*!

!.... in SUNDAE
  kappa = 1e6
  niteration=10000
  lanczos_step=1.0d-3
!.... in SUNDAE 

  tempdeplainit=-1
  debyetemp=-1
  lprtpot=.false.
  ibrake =0   ! if =1 electronic slowing for cascades (acting on all atoms)
  ngrdel=500

  timemax=1d20
  tpseuils(:)=0 ! 1:Tmin; 2:abs(T') ; 3: abs(T'') ; 1:abs(P); 2:abs(P') ; 3: abs(P'')

  lrctest=.true.

  if (rang == 0) write (6, *) 'nom fichier din=', fnamdin

  open(unit=ludin, file=fnamdin, status='unknown', err=456)
  read (ludin, nml=input)



  imm_glob = imm
#if(PARA)
  ! En parallele, on initialise le nombre maximum d'atomes d'un
  ! processus au nombre d'atomes locaux. Plus tard ce nombre sera
  ! complete par le nombre maximal d'atomes fantomes
  ! On suppose que la concentration max ne depasse pas 20%  de 
  ! la concentration moyenne
  imm      = min( imm_glob, int(1.2 * imm_glob / nprocs) )
  if (rang==0) write(6,*)'IMM PARA = ',imm
#endif




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
     if (rang==0) write(6,*)'lpcon n_existe plusest historique utiliser plutot lpr pour un Parinnello Rahman propre '
     if (rang==0) write(6,*)
     if (rang==0) write(6,*)
     if (rang==0) write(6,*)
     if (rang==0) write(6,*)'!!!!!!!!!!!!!!!!!!!'
     if (rang==0) write(6,*)'!!!!!!!!!!!!!!!!!!!'
     if (rang==0) write(6,*)'!!!!!!!!!!!!!!!!!!!'
     call arret_ndm
  end if

  if(.not.lperiod) then
     if (rang==0) write(6,*)
     if (rang==0) write(6,*)'Pas de conditions periodiques'
  end if


  tstep = tstep*utemps
  tauTcon=tauTcon*utemps
  epcou = epcou*1D-8
  tdepla = tdepla*1D-8
  tdepla2 = tdepla2*1D-8
  rsep=rsep*1.0d-8
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

  if(ljqbh) then
     !     open(unit=59,file='temptranche.mol')
     dmtype=1
     read(94,*)njqbh, epsil, epcoud,ittherm,ntr
     if (parallele) then
        if (njqbh.ne.4) then
           write(6,*)'njqbh  doit etre =4 en para'
           call arret_ndm
        end if
     end if
     if (njqbh==5) read(94,*) kthg
  end if


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

  if (parallele) then
     if (ltabvois) then
        ltabvois=.false.
        if (rang==0) write(6,*)'LTABVOIS MIS A FALSE en PARA'
     end if
     select case(dmtype)
        case(2,4)
        case default 
        if (rang==0) write(*,*) 'FATAL: VERSION PARALLELE seulement avec dmtype=4'
        if (rang==0) write(*,*) 'Stop in readdm'
        call arret_ndm
     end select
  end if


#if(DECOUP)
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


  if (dmtype.EQ.11) itmax=1
  if (itmax < 0) then
     if (rang==0) write (6, *) rang,'wrong itmax < 0 '
     call arret_ndm
  endif

  if (itedepla>=1 .and. tdepla<0.0) then
     if (rang==0) write (6, *) rang,'wrong tdepla < 0 '
     call arret_ndm
  endif

  !     if ((itedepla.ge.1).and.(tdepla2.lt.0.0)) then
  !       write(6,*)'wrong tdepla2 < 0 '
  !       call arret_ndm
  !     endif

  if (tdepla2>0.0 .and. tdepla2<tdepla) then
     if (rang==0) write (6, *) rang,'choisissez tdepla2 > tdepla '
     call arret_ndm
  endif

  if (lfilm .and. lfilmext) then
     if (itedepla < 1) then
        if (rang==0) write (6, *) rang,'contradiction itedepla <-> lfilm et lfilmext'
        call arret_ndm
     endif
  endif

  if (itedepla<1 .and. lfilm .and. (.not.lfilmext)) then
     if (rang==0) write (6, *) rang,'contradiction itedepla <-> lfilm '
     call arret_ndm
  endif

  if (itedepla<1 .and. lfilmext .and. (.not.lfilm)) then
     if (rang==0) write (6, *) rang,'contradiction itedepla <-> lfilmext '
     call arret_ndm
  endif

  !      if (ttol==0.0 .and. tfroi<=0.0) then
  !         write (6, *) 'contradiction ttol <-> tfroi '
  !         stop
  !      endif

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
  if(dmtype==12) lprteat=.true.
  if(dmtype==16) lprteat=.true.
  if (lposmoy.EQV..true.) then 
     lprteattotm=.true.
     write(6,*)'LPOSMOY, stocke les positions moyennes dans posmoyx et les ecrit a la fin avec les energies moyennes'
  end if
  if (lprteattotm.EQV..true.) then
     lprteat=.true.
     write(6,*)'LPRTEATTOTM calcule les energies moyenne de chaque atome et les ecrit en retranchant eatref en eV (=0 par defaut)'
  end if
  if ((lprteattotm.EQV..true.).and.(parallele.EQV..true.))then
     write(6,*)'eattotm et PARA pas prog'
     stop
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
              stop
          end if    
     write(6,*) rang,'For dmtype 7 deltax must be deltax > 0'
     write(6,*) rang,'Change deltax!'
     call arret_ndm
  endif

  if (dmtype==7) then
        ldemitab=.false.
        if (.not.((HessianOrder.eq.1).or.(HessianOrder.eq.2).or.(HessianOrder.eq.4))) then
         write (6,*) ' PHONDY: HessianOrder can have only the values 1, 2 or 4  '
         write (6,*) ' PHONDY: which corresponds to a Hessian on 2,3 or 5 points' 
         write (6,*) ' PHONDY: HessianOrder.........: ', HessianOrder  
         write (6,*) ' PHONDY: stop'
         stop
        end if
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
        stop
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
  if (formatsauv>3 .or. formatsauv<0) then
     if (rang==0) write (6, *) rang,'mauvais formatsauv = ', formatsauv
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


  if((lTberendsen).and.( (dmtype.EQ.2).OR.(dmtype.EQ.3).OR.(dmtype.EQ.30) ) ) stop
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
        call endrun 
     end if
  end if
  !  if ((all(lpotentiel)==.false.).and.(ipotentiel==-1)) then
  !  end if

  if (ipotentiel.ge.0) lpotentiel(ipotentiel)=.true.
  npotentiel=0
  do ipotcont=0,npotmax
     if (lpotentiel(ipotcont).EQV..true.) npotentiel =npotentiel+1
  end do
  if ((lpotentiel(0).EQV..true.).and.(npotentiel.gt.1)) then
     if (rang==0) write(6,*)'npotentiel>1 et lpotentiel(0)=T'
     stop
  end if

  if ((ntyp==-1).and.(npotentiel.gt.1))then
     if (rang==0) write(6,*)'npotentiel>1 et ntyp=-1'
     stop
  end if
  if ((npotentiel.gt.1).and.(lpotentiel(10).eqv..true.)) then
     if (rang==0)write(6,*)'**** npotentiel >1 ET EAM ==> EAM TAB only!'
  end if


  if ((lpotentiel(12).EQV..true.).and.(ltabvois.EQV..true.))ldemitab=.false.

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
     if ((dmtype.eq.1).or.(dmtype.eq.2).or.(dmtype.eq.4)) then
        if (rang==0) write(6,*) 'The time step changed each', itetimestep,' steps'
     else 
        if (rang==0) write(6,*)'itetimestep seulement avec dmtype =1, 2  or 4 '
        if (rang==0) write(6,*) 'STOP in readdm'
        stop
     end if
  end if

  if(lLangevin) then
     dmtype=4
     write(6,*)'langevin bugguÃ© voir Cosmin fabien'
     stop
  end if
  if(lLangevin.and.(Text.le.0.0)) then
     if (rang==0) write(6,*)'Langevin avec Text pas defini : stop'
     stop
  end if

  ! end check
  if (lpconxyz) then
   if (.NOT.lpr) then
    if (rang==0) then
       write(6,*) 'lpconxyz can be used only is with PR dynamics or lpr=.true'
       write(6,*) 'stop in <readdm>'
    end if
     stop
 end if 
end if       
   
   
  if (lpr) then
     if (dmtype==2) lprtrp=.true.
     dmtype=8
     itesigma=1
     if (pext.ne.0.) then
        if (rang==0) write(6,*)'SIGEXT', sigext
        if (rang==0) write(6,*)'Pext ',pext
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
     ihbox0(:,:) = 1.d0   ! all the dimension of the box can change
     if (lpconxyz) then
      ihbox0(:,:)=0.d0   ! ALL the dimension are blockef except ...
      ihbox0(1,1)=1.d0   ! X ...
      ihbox0(2,2)=1.d0   ! Y ...
      ihbox0(3,3)=1.d0   ! and Z.
     end if 

  end if
  ! read for cascade
  if (lcasca) then
     read (ludin, *) iko, eko, xko, yko, zko , xx0, yy0, zz0
     if (lrestart) lcasca=.false.
     !     xx0=xx0*1.D-8
     !     yy0=yy0*1.D-8
     !     zz0=zz0*1.D-8
     if (ibrake.gt.0) then
if(rang==0) then 
       write(6,*)'electronic stopping according to elstop.in from SRIM POUR DES EC > 1!!!!!!!!!!'
       write(6,*)'electronic stopping according to elstop.in from SRIM POUR DES EC > 1!!!!!!!!!!'
       write(6,*)'electronic stopping according to elstop.in from SRIM POUR DES EC > 1!!!!!!!!!!'
       write(6,*)'electronic stopping according to elstop.in from SRIM POUR DES EC > 1!!!!!!!!!!'
       write(6,*)'electronic stopping according to elstop.in from SRIM POUR DES EC > 1!!!!!!!!!!'
       write(6,*)'electronic stopping according to elstop.in from SRIM POUR DES EC > 1!!!!!!!!!!'
endif
     end if
  endif



  if (lHcyl) then
     read (ludin,*) Ecyl, pc(1), pc(2), pc(3), vdc(1), vdc(2), vdc(3), rayonc, lgc 
     if (lrestart) then
        lHcyl=.false.
     end if
     pc(:)=pc(:)*1.D-8
     vdc(:)=vdc(:)*1.D-8
     rayonc=rayonc*1.D-8
     lgc=lgc*1.D-8
  end if


  ! Initialisation

  if (iterdf >= 0) then

     !         lurdfout = 88
     !         fnamrdfout = fnam(1:lenfnam)//'.rdfout'
     !         open(unit=lurdfout, file=fnamrdfout, status='unknown')
     gdertot=0.0
     if (iterdf > 0) nrdf = 0
     if (iteangle > 0) nfda = 0
  endif

!  if (itedepla.gt.0) lfilm=.true.
  if (lfilm) then
     if(rang==0) then
        if(rang==0)         open(unit=lufilm, file='film', status='unknown')
        if(rang==0)         open(unit=lufilmpaf, file='filmpaf', status='unknown')
     end if
  endif
  if (tcooling > 0) lastcool = 0.0

  if(dilat(1).ne.0.0) then
     if(igen.ne.1) then
        if (rang==0) write (6,*) rang,'dilat<>0 et igen<>1 stop'
        call arret_ndm
     end if
     if (dilat(2)==0.0) dilat(2)=dilat(1)
     if (dilat(3)==0.0) dilat(3)=dilat(1)
     if (rang==0) write(6,*)'dilat 1 2 3 ', dilat(1),dilat(2),dilat(3)
  end if

  ! MPI

  if (rang==0) write (6, *)
  if (rang==0) write (6, '(a)') ' -------- caracteristiques du run DM--------'

  select case (dmtype)
  case (1)
     if (rang==0) write (6, '(a)') '     DYNAMIQUE MOLECULAIRE VERLET STANDARD'
  case (2)
     if (rang==0) write (6, '(a)') '     TREMPE RAPIDE'
  case (3)
     if (rang==0) write (6, '(a)') '     GRADIENT CONJUGUE sur les coordonnees CARTESIENNES'
  case (30)
     if (rang==0) write (6, '(a)') '     GRADIENT CONJUGUE sur les coordonnees REDUITES'
  case (4)
     if (rang==0) write (6,'(a)') '      DYNAMIQUE MOLECULAIRE VELOCITY VERLET'
  case (5)
     if (rang==0) write (6,'(a)') '      TEST DES FORCES '
  case (7)
     if (rang==0) write (6,'(a)') '      CALCUL DES PHONONS A PARTIR DE POSITIONS DE FORCES NULLES '
  case (8)
     if (rang==0) write (6,'(a)') '      PARRINELLO RAHMAN AUTOCOHERENT '
  case (6)
     if (rang==0) write (6,'(a)') '      ANALYSE DES POSITIONS EN FIN DE CASCADE '
  case (9)
     if (rang==0) write (6,'(a)') '      DRAG OR NEB DYNAMICS ' 
     itesauvposition=-1
     itesauvforce=-1
     itetemp=-1;itesigma=-1
  case (10)
     if (rang==0) write (6,'(a)') '      TREMPE FIRE '
     if (rang==0) write (6,'(a)') '  !!!  Attention les masses atomiques sont toutes celle du type 1!!! '
     if (rang==0) write (6,*)
  case (11)
     if (rang==0) write (6,'(a)') '      UN CALCUL DE FORCES '
     if (rang==0) write (6,*)
#if(ART)    
  case (12)
     if (rang==0) write (6,'(a)') '|=========NDM ENTERTAINMENTS presents:===============|'
     if (rang==0) write (6,'(a)') '|---------ART nouveau by N MOUSSEAU.---------------|'
     if (rang==0) write (6,'(a)') '|======== colored by Cosmin Marinica!==============|'
#endif
#if(SUNDAE)    
  case (16)
     if (rang==0) write (6,'(a)') '|=========       NDM + SUNDAE       ===============|'
     if (rang==0) write (6,'(a)') '|---------..........................---------------|'
     if (rang==0) write (6,'(a)') '|==================================================|'
#endif
#if(MAB)    
  case (17)
     if (rang==0) write (6,'(a)') '|=========       NDM + MAB          ===============|'
     if (rang==0) write (6,'(a)') '|---------..........................---------------|'
     if (rang==0) write (6,'(a)') '|==================================================|'
#endif
  case default
     if (rang==0) write (6, '(a)') 'mauvais type de calcul dmtype=', dmtype
     stop
  end select

  if (ltranche) then
     if (rang==0) write (6, '(a)') '******************* TRANCHE GELEE !!! *****'
!     rulayer=rulayer*1.0d-8

     lfrozen=.true.
     if (lcdp.EQV..true.) then
        write(6,*)'TRANCHE +DP = PAS POSSIBLE' ; stop
     end if


  end if


  if (lcontr)  write (6, '(a)') '******************* CONTRAINTE !!! *****'


  if (lTcon) then
     if (rang==0) write (6, *) 'TEMPERATURE CONSTANTE a la main Text= ',text
  endif
  if (lTcon) then
     if (rang==0) write (6, *) 'TEMPERATURE CONSTANTE a la Berendsen Text= ',text
  endif

  select case (igen)
  case (-1)
     if (rang==0) write (6, *) 'generation du crystal'
  case (0)
     if (rang==0) write (6, *) 'generation du crystal ; puis run'
  case (1)
     if (rang==0) write (6, *) 'run a partir du fichier .cin'
  case (2)
     if (rang==0) write (6, *) 'modification du fichier .cin'
  case default
     if (rang==0) write (6, *) 'mauvais igen=', igen
     stop
  end select

  if (ltabvois) then
     if (npotentiel.gt.1) then
        if (rang==0) write(6,*)'ltabvois avec plusieurs potentiels= pas programmee (demi table ou table complete = prise de tete'
        stop
     end if

!     if (tempstopcel.gt.0) ltpcel=.true.
     if (ltpcel) write (6, *) '   -> -> pas de contrainte par celulles'

     if (rang==0) write(6,*)'IPOTENTIEL',ipotentiel
     select case (ipotentiel)
     case(:9)
         ldemitab=.TRUE.
         if (rang.eq.0) write (6, *) '    DEMI-TABLE DES VOISINS rvois ',rvois
      case(11:)
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
      end select

     !     if(ipotentiel.le.10) then 
     !        ldemitab=.TRUE.
     !        if (rang.eq.0) write (6, *) '    DEMI-TABLE DES VOISINS rvois ',rvois
     !     else
     !        ldemitab=.false.
     !        if (rang.eq.0) write(6,*)'    TABLE DES VOISINS COMPLETE rvois ',rvois
     !     end if

  endif

  if ( (dmtype==2).or.(dmtype==3).or.(dmtype==30).or.(dmtype==9).or.(dmtype==10) ) then    
     if ( (fpstop<0).and.(fsumstop<0)) then
        if (rang==0) write(6,*) 'One of fpstop and fsumstop must be positive for dmtype=',dmtype
        if (rang==0) write(6,*) 'STOP in readdm'
        stop
     end if
     if ( (fpstop < 0) .and. (dmtype==9) ) then
        if (rang==0) write(6,*) 'NEB and DRAG implementation only for positive fpstop'
        if (rang==0) write(6,*) 'STOP in readdm'
        stop 
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
        stop
     end if
  end if

  if ( (.not.lperiod).and.(itesauvposition>0).and.lsuivinonpbc) then  
     if (rang==0) write(6,*)' lsuivinonpbc will be turn to FALSE'
     lsuivinonpbc=.false.
  end if
  if (lsuivinonpbc) then

     if (dmtype.ne.4) then
        if (rang==0) write(6,*) 'lsuivinonpbc is implemented only with velocity verlet'
        if (rang==0) write(6,*) 'Or dmtype=4. Change and restart until there I will stop for you.'
        stop 
     end if

     if (itesauvposition<=0) then
        if (rang==0) write(6,*) 'itesauvpostion MUST be positive if you want lsuivinonpbc TRUE'
        if (rang==0) write(6,*) 'STOP in readdm'
        stop
     end if
  end if
  if (lforcetabulate) then
    if (ipotentiel/=10) then
     write(*,*) 'There is no implementation for lforcetabulate TRUE and ipotentiel ', ipotentiel
     write(*,*) 'Change lforcetabulate of FALSE or ipotential to EAM (10) '
     write(*,*) 'stop in readdm'
     stop
    end if
    if (lcasca) then
     write(*,*) 'There is no implementation for lforcetabulate TRUE and lcasc TRUE'
     write(*,*) 'Change lforcetabulate of FALSE or lcasc on FALSE'
     write(*,*) 'stop in readdm'
     stop
    end if
  end if

  if (itetemp2==-1) itetemp2=itetemp
  if (rang==0) write(6,*)
  if (rang==0) write (6, *) '     ANALYSES '
  if (rang==0) write (6, *) 'itetemp=', itetemp, ' itesigma=', itesigma
  if (itecoordo>0)  write(6,*)  ' itecoordo=', itecoordo
  if (itedepla>0) then
     if (rang==0) write (6, '(A,I3,A,D9.3,A,D9.3,A,I3,A,I3)') ' itedepla=', itedepla, &
          ' tdepla=', tdepla*1D+8, ' tdepla2=', tdepla2*1D+8, ' itesauv=', &
          itesauv, ' itesauvposition=', itesauvposition
  end if
  if (lfilm) then
     if (rang==0) write (6, '(A,D11.3)') ' film; seuil=', tdepla*1D+8
     if (rang==0) write (6, '(A,D11.3)') ' film; seuil2=', tdepla2*1D+8
  endif
  if (lfilmext) then
     if (rang==0) write (6, '(A,D11.3)') ' filmext; seuil=', tdepla*1D+8
     if (rang==0) write (6, '(A,D11.3)') 'filmext; seuil2=', tdepla2*1D+8
  endif
  if (rang==0) write(6,*)
  if (rang==0) write (6, *) '     CONTROLES '
  if (rang==0) write (6, '(A,D11.3)') 'tstep=', tstep
  if (rang==0) write (6, *) 'itmax=', itmax, ' itab=', itab, ' itetimestep=', &
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



  if (ttol>0.)  write(6,*) ' ttol=', ttol
  if (tfroi>0.) write(6,*) ' tfroi=' ,tfroi
  if ((tempstop>0.).and.(rang==0)) write(6,*)' tempstop=', tempstop
  if ((tfcou>0.).and.(rang==0))         write (6, '(A,F10.1,A,F10.1,A,F10.1)') 'tfcou=', tfcou, ' epcou=', &
       epcou*1D+8
  if (tcooling>0.) write(6,*)' tcooling=', tcooling

  if (lcasca) then
     if (rang==0) write (6, *) '----CASCADE-----'
     if (rang==0) write (6, *) 'projectile=', iko, ' energie=', eko
     if (rang==0) write (6, *) 'direction=', xko, yko, zko
     if (rang==0) write (6, *) 'position de depart : ', xx0, yy0, zz0
  endif

  if (rang==0) write (6, *) ' ----------------------------------'
  if (rang==0) write (6, *)
  if (rang==0) write (6, *)



  if(lTandersen.and.rang==0) write(6,*)'Tandersen nuandersen = ',nuandersen

  if ((lsigtyp).and.(rang==0)) then
     write(6,*)
     write(6,*)'SIGTYP programme en 2 corps cellule seulement pour le type 3!!!!!!!!!!!!!!!!!!!!!!'
     write(6,*)'PRESSION A DIVISER PAR LES VOLUMES !!!!!!!!!'
     write(6,*)
  end if
  if(lPrtSigat.and.(.not.ltabvois)) then
     write(6,'(a)')rang,'contrainte atomique programme en table des voisins&
                & avec un potentiel EAM ou un terme a deux corps seulement'
     call arret_ndm
  end if


  if(itetemp2.gt.0)then
     fnamdin = fnam(1:lenfnam)//'.2.T'
     open(unit=112, file=fnamdin, status='unknown')
  end if

  if (ldislo.and.((epcoudis==0.0).or.(fdislo==0.0))) then
     epcoudis=epcoudis*1.0d-8
     write(6,*)rang,'ldislo vs epcoudis ou fdislo stop'
     call arret_ndm
  end if


  if ( ( (dmtype==3).OR.(dmtype==30) ) &
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
  IF (lFrozen.OR.lxyzFrozen) THEN
          lxFrozen=.true. ; lyFrozen=.true. ; lzFrozen=.true.
  END IF
  IF (lxyFrozen) THEN
          lxFrozen=.true. ; lyFrozen=.true.
  END IF
  IF (lxzFrozen) THEN
          lxFrozen=.true. ; lzFrozen=.true.
  END IF
  IF (lyzFrozen) THEN
          lyFrozen=.true. ; lzFrozen=.true.
  END IF
  
  IF (lxFrozen.OR.lyFrozen.OR.lzFrozen) THEN

     IF ( ( (imFree.gt.0) .AND. (parallele) ).OR. ( (imFirstFrozen/=0) .AND. (parallele) ) ) THEN
           WRITE(0,'(a)') 'Initialisation du tableau free(:) pour&
                & determiner les atomes bloques non implementes en&
                & parallele'
           STOP '< ReadDm >'
     end IF

     IF (dmType.EQ.8) THEN
        IF (RANG==0) WRITE(0,'(a)') 'Vous ne pouvez pas utiliser&
             & Parrinello-Rahman tout en maintenant fixes certains&
             & atomes'
        STOP '< ReadDm >'
     END IF

     if ((imFree==-1).and.(rulayer==0.0).and.(imFirstFrozen==0))then
        if(rang==0) write(6,*)'LFROZEN+IMFREE=-1 et RULAYER=0 et imFirstFrozen==0 == stop'
        stop
     end if

     ! Le tableau free controle quels atomes participent a l'energie (utilise par JP a priori)
     ! Le tableau frozen controle quelles coordonnees de quels atomes sont libres de relaxer
     !    i.e. quelles forces doivent Ãªtre annulees
     Allocate(Free(1:imm))
     Free(:)=.true.
     Allocate(Frozen(1:3,1:imm))
     Frozen(:,:)=.false.   ! Tout le monde bouge ... ...

     if (imFree.ne.-1) then
             IF (lxFrozen) Frozen(1,1+imFree:imm)=.true.   !x of i>imFree is frozen  
             IF (lyFrozen) Frozen(2,1+imFree:imm)=.true.   !y of i>imFree is frozen  
             IF (lzFrozen) Frozen(3,1+imFree:imm)=.true.   !z of i>imFree is frozen  
             IF ((rang==0).and.(imFree.ne.imm)) THEN
                     WRITE(6,'(a)') 'Dynamique / relaxation avec des atomes bloques'
                     if(imFirstFrozen==0) then
                      WRITE(6,'(a,i0,a)') "Seuls les atomes d'indice inferieur ou egal a ", &
                        imFree, " bougent"
                     else 
                      WRITE(6,'(a,i0,a,i0,a)') "Seuls les atomes d'indice inferieur ou egal a ", &
                        imFree, " et superieur et egal a ", imFirstFrozen+1, "bougent"
                     end if 
             end IF
     end if

     if (imFirstFrozen > 0) then
             IF (lxFrozen) Frozen(1,1:imFirstFrozen)=.true. ! x of i<=imFirstFrozen is frozen
             IF (lyFrozen) Frozen(2,1:imFirstFrozen)=.true. ! y of i<=imFirstFrozen is frozen
             IF (lzFrozen) Frozen(3,1:imFirstFrozen)=.true. ! z of i<=imFirstFrozen is frozen
            if (imFree == -1) WRITE(6,'(a)') 'Dynamique / relaxation avec des atomes bloques'
            if (imFree == -1) WRITE(6,'(a,i0,a)') "Seuls les atomes d'indice superior ou egal a ", &
                        imFirstFrozen+1, " bougent"
     end if 




     lFrozen=.true.

  end IF   !lxFrozen,lyFrozen,lzFrozen
  if (rulayer.gt.0.0)then
     rulayer=rulayer*1.0d-8
     
     IF (rang==0)write(6,*)'atomes immobiles fixes par rulayer ', rulayer*1d8
  end if


  !  if (rang==0) write(6,*) 'sortie readdm'



  if ((iteanapos.eq.-1).and.(itecompcr.ne.-1)) iteanapos=itecompcr
  if(deltaestop.ne.0)  write(6,*)'<<arret apres chnegement de Epot moyenne ', deltaestop, nbmoye



  usdh = 1/(two*tstep)

  if (rang==0)then 
     write(6,*)'nb de potentiels', npotentiel
     if (npotentiel==1) then
        write(6,*)'ipotentiel',ipotentiel
        !        write(6,*)lpotentiel
     else
        do ipotcont=1,npotmax
           if (lpotentiel(ipotcont).EQV..true.)write(6,*)'potentiel actif', ipotcont
        end do
!        if (lcasca.eqv..true.) then
!           if (rang==0) write(6,*)'ATTENTION!!! npotentiel>1 et ziegler surement faux !!!!'
!           stop
!        end if
     end if
  end if
  if (rang==0) write(6,*)'fmt_cin',fmt_cin

  if (ldesinteg)then
     !     itmax=nstepdes
     itesauv=0.
     if (dmtype.ne.4) then
        write(6,*) 'dmtype <> 4 et linsert'
        stop
     end if
     if (nstepdes.le.0) then
        write(6,*) 'desinteg et nstepdes<1'
        stop
     end if

     if(tempdes==-1)tempdes=Text

     if (rang==0)write(6,*)
     if (rang==0)write(6,*)'desintegration de l atome ',ides, 'mis Ã  1'
     if (rang==0)write(6,*) 'desinteg NE FONCTIONNE QUE AVEC DES POT DE PAIRES !!!'
     if (rang==0)write(6,'(A,D12.5,A)')' kspr=',kspr,'eV/Ang**2'
     kspr=kspr*1d16/erg2eV
     if (rang==0)write(6,*)
     if (rang==0)write(6,*)'ATTENTION EN PARA LES ATOMES NE DOIVENT PAS TROP VOYAGER PENDANT LES CHEMINS'
     if (rang==0)write(6,*)'ATTENTION EN PARA un atome ne doit pas aller d-un proc. Ã  un proc non voisin'



     if (xpspr(1)==-1000) then 
        if (rang==0)  write(6,*) 'position du ressort sur la position de l_atome ides'
     else
        if (rang==0) write(6,*) 'position du ressort',xpspr,'Ang'
        xpspr=xpspr*1d-8
     end if

     if (rang==0)write(6,*)        
     lambdades=1.0 ; pm1des=-1

  end if

  rheat=rheat*1d-8
  if (lheat.eqv..true.)then
     Eheat=Eheat*ev2erg
     dmtype=4
     if (rheat.le.0)then
        write(6,*)'heat et rheat<=0 stop'
        stop
     end if
     if((Eheat==0).and.(Theat==0)) then
        write(6,*)'heat et Eheat=0 Theat=0 stop'
        stop
     end if
     if((Eheat.ne.0).and.(Theat.ne.0)) then
        write(6,*)'heat et Eheat<>0 Theat<>0 stop'
        stop
     end if
  end if

  return
456 print *,'Erreur lors de la lecture du fichier .din, verifier l''ajout de fmt_cin'
end subroutine readdm


