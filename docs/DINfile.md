# din file (e.g. name.din)
General driving file for the MD calculation
This file is a fortran namelist with a few mandatory and many optional variables  
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
         itesauvinter,units_lammps,lWgin,lvzeroneb,pas_lambda_mc,n_path,lax,ldecoup,distminat,&  
         ndir,nstep,betaguess,ncgtry,lvarstop,fstpdecr,itypcalc,gamprfact,TinitBox,&  
         &nparapath,lparapath,lrestartmcgc, lbiais_retrait,lbiais_inser,fdmc_1,&  
         &fdmc_2,ndecal,decal,lparafm,nparafm,lwritefreq,lwfm,ldecalcor,kmin,kmax,iteprtkin,lspecialinit,&  
         &noxyzkmin,noxyzkmax,lpartarps,kspring=1
	 
##Defaults values and explantaions Mandatory variables are denoted with MMM

**imm = 0**  MMM dimension of the atomic arrays (must be larger or equal to the maximum number of total atoms in the box  
**tstep = 1.0** timestep in 10^-15 sec unit  
**tinit** = -1.0 initial temperature  
** lvpread**=.true. if .cin is read, reads the velocities  
Velocities are rescaled to tinit if they are read ,  otherwise randomly inititaed  
**igen = -2  MMM generation of atomic positions :0 from gin file +1 from cin file ;  -1 from gin to cin then stop (cin to gin with lwgin=.true.) **
**itmax** = -1  MMM (maximum) number of iterations  
timemax=1d25 maximum simulated time in fs  
**lrestart = .FALSE.  if T : restarting from an interrupt job (reads a .cout file)**
**dmtype = 0  MMM dmtype = type of calculation :**  
   **1** -> Standard Verlet  
   **2** -> quench fast quenching algorithm (==21) For quenching I suggest dmtype = 23 or 24   
   21 -> quench constant volume  
   22 -> quench constant pressure   
   ***23*** ->  fire quench algorithm constant volume  
   ***24*** -> fire quench constant pressure  
   3 -> old style conjugate gradient ==31 , remains for historical reasons, use at your own risks  
  30 ->old style conjugate gradient on reduced coodinates  , remains for historical reasons, use at your own risks  
  31 ->old style conjugate gradient on cartesian  coodinates  , remains for historical reasons, use at your own risks  
  32 -> steepest descent  
  33 -> new conjugate gradient  
  34 -> new conjugate gradient with Fletcher-Reeves  modification
  35 -> relaxation with  ADAM D. :Kingma and J. Ba, “Adam: A Method for Stochastic Optimization,” in International Conference on Learning Representations (ICLR), 2015.  
   **4** -> Velocity Verlet  
   41 -> Velocity Verlet ARPS with partial forces  
   42 -> Velocity Verlet ARPS with regular forces  
   5 -> force test  
   6 -> position analysis  
   **9** -> NEB calculaation  uses deb_calc.gin and fin_neb.gin files are starting and end points  
  **12** -> ART calculation  <---  16  SUNDAE  
  11  UN SEUL CALCUL DE FORCES 
  17  MAB --->  
  18 -> Milady fit from NDM  
  19 -> thermodynamics from force matrix diagonalization  
  15 -> Path Monte-Carlo calculation (==151)  
  151 -> Path Monte-Carlo calculation    OTHER MCC are not documented  
  112 -> histogram of distances between atoms

*nox* = -1 default : the box is automatically divided in cells with as many cells as possible w.r. the potential cutoff  
*noy* = -1 nox/noy/noz are then output for information    
*noz* = -1 at the opposite they can be specified in din file. Their maximum value being the default calcuated value. This can be useful to avoid nasty value that would spoil paralelization (e.g. prime number)  

##potential 
**ipotentiel = -1 MMM  defines the potential**  
> 0=Born-Mayer-Huggins,  
1=Buckingham,  
2=watanabe,  
3=buck8,  
4=Morelon UO2,  
5  Morse term,  
6=Stillinger Weber a  la Vashista  
 **7 tabulated pair (with or without charge)   **    
 **10 EAM;**  
 12 ZrC JuLi(+Tersoff Doan)  ;  
 13 Tersoff coupure COS;  
 14 Tersoff coupure FD ;  
 15 tersoff coupure SIN (original) ;  
 11 Ercollesi ;;  
 **-10=LAMMPS atom style atomic;**  
 **-11 LAMMPS atom style charge (changes only simple.potin)**  
 8 bandura 2017= Bukingham +Morse+Fermi-Dirac+Inverse gaussian  
 **16= CRG**  
 **20=MILADY**  
  
npotentiel = 1   number of potentials if >1 ipotential is used instead of ipotential   
lpotentiel(:)=.false. if npotential >1 lpotential indicates the used potential    
 ntyp=-1 the number of types in the simulation must be specified ONLY in npotentiel >1  

##algorithm choice
lFire = .false. ! if true with dmtype=2 forces fire algorithm  

**Text=-1.  (MMM) the external temperature, mandatory for cst temperature algorithms**  
**lLangevin=.false. Langevin MD for cst temperature**
wboxf = 1.0    mass of the box for PR algo (by default 0.5*SUM(mass all atoms)   
wNose = 0.0  mass of the box for Nose algo (default=wbox)   
tbox = 1000.0  relaxation time for the box    
lTcon=.false. if true isokinetic algo   
lTberendsen=.false.  T: Berendsen cst temp    
lTNose=.false. T: Nose  cst temp   
lTHoover=.false.  T: NoseHoovzr  cst temp     
nHoover=1  number of Hoover chains
tauTcon=200.0 The rescales "time" for the Berendsen algorithm  
lTandersen=.false. temperature constante a la Andersen   
nuandersen=1.0d14 frEQUENCY OF RANDOM SAMPLING OF VELOCITIES FOR aNDERSEN   
dfpred = 0.1 eguess for GC calculations and quenching  
lprahman=.false.  parinnelo rahman at constant stress
lUcell=.false.  if true plots the potential energy of the box. This assumes that the reference state is for 0 stress  
fpstop =-0.05  stopping criterion for quenching : maximum force per atom  UNIT = EV/ANG   
sigstop =-0.05 stopping criterion for quenching : maximum stress component  UNIT =  kbar  
fsumstop =-0.1 stopping criterion for quenching : maximum SUM of forces   UNIT = EV/ANG  
lcontr=.false.  dynamique contrainte (routine contrainte)  

**lpr=.false.  ! forces parinello-Rahman calculation for dynamics or quenching**  
**pext = 0.0  extrernal pressure in kbar**
lpcon2 = .FALSE. if true a damping term is added on the dynamics of the box  
lpconxyz = .FALSE. T: the relaxation are allowed only along the X, Y and Z axis wonceived for a tetragonal box  
sig0stara (sig0starb or sig0starc) =.false. If sig0stara=.true. , the cell vectors b and c are constant,. cell vector a is changed so that the force normal to the (b,c) plane (along stara) goes to zero  
lpr=lprahman  
ihbox0(:,:) = 1 all the dimension of the box can change. Set to zero to fix a cell lenght or angle  
sigext = 0.0 Symetric tensor of  the external stress  
h0(1:3,1:3) = 0.d0  reference box for the Parinello Rahman algorithm, by default equals the inkitial box    


gamlg =5d12 Gamma Langevin (= 0.005/1d-15 will do  vp*0.995 for  tstep=1d-15)  
gamprfact=0.1 factor of langevin for variable cell Langevin  
ilangevin=1  std Langevin type


##analysis variables  
itetemp = 100 period of temperature calculation  
tempstop = -1.0  stopping temperature for the run  If tempstop > 0 the run will stop on the first step whre the calcuated temperature is lower than tempstop  
tempstopcel = -1.0  stops if temperature in ALL cells are below temstopcel  
itesigma = -1   period of stress calculation  
iteprtsigma = -1   period of stress calculation  
itedepla = -100  period of displacement cal.  
tdepla = 1.0 threshold for displacement  
tdepla2 = -1.0  second threshold for displacement    
lfilm = .FALSE. film making of displaced atoms  
itecoordo = -100 period of coordination calculation  
itean = 0  general control for analysis  
iterdf = -1 period of RDF calc. : -1 never ; 0 : nrdf last iterations; +iterdf every iterdf iterations  
linstantrdf = .FALSE.   F: calculates and prints the average of the RDF, T : calculent RDFrdf iterations  
linstantfda = .FALSE.   F: calculates and prints the average of the RDF, T : calcuatent RDFrdf iterations  
nrdf = 0 period of the averaging of the rdf calc   
nfda=0  idem angular analysis  
lEparat=.false.  calcul et affichage de l'energie par atom  
itebdv=-1   frequence de calcul des bond valence  
itetemp2=-1frequence d'ecriture de la temperature dans fichier separe  
itecompcr=-1   remplacee par iteanapos  
iteanapos=-1   frequence de comparaison avec cristal de reference  
lposmoy=.false.writes the average position and energy of the atoms in a .mol file  
eatref(:)=0.  
lanaposart=.false.   anapos a la ART : decalage + defauts en WS, concu pour le cas des I dans UO2  


##miscellanae
itab = 10 period of cell repartition , can be increaszd if atoms will not move in the simulations  
itetabvois = 10  period  of neighbour table calculation active only if ltabvois=.true.  
ltabvois = .FALSE. true : builds a neighbour table for force calculations works only in sequential   
lconstrtot=.FALSE. T: builds the neighbour table from a complete double loop ; F builds neighbour table from cell dispatching   
rvois = 0.0  MMM if ltabvois=.true. radius of neighbour table   
itetimestep = -1  period of check in timestep if >0 adapts timestep for a mximum displacement every itetimestep     
tsfact = 10.0   change in time step factor  
tinitbox = -1.0 initial temperature of the box  
tfcou = -1.0 temperature of the border of the box  
epcou = -1.0width of the border of the box  
couxyz(1:3)=1  if set tà 0 temperature is NOT controlled on these surfaces  
nitmax = -1 maximum number of new iterations after restart   
ldecoup = .FALSE.  if True poutputs an efficiency score for the parallel calculations up to some numbner of cores (works only in sequentiall   

ltpcel = .FALSE. output of temperature and stress in each cell  
itesauv = 1000   period for saving  
itesauvposition = 0 period of saving of positions in binary file   
itesauvinter=0  
lcdp=.false. point defect accumulation algo  
iseed=0  si <>0 controls  the random velocity initialization (for calculation reproduction)  
lEev=.true.  T: output in eV (F output in cgs)
lPkbar=.true.   T: output inkbar (F output in cgs)
**iterasmol = -1**  \>0 formatted output every iterasmol steps  
**ivisu=4** format de sortie dans rasmol.f90 : ivisu=1=.mol, ivisu=2=vsim mal code supprime, ivisu=3=xred , ivisu=4 CFG, ivisu=6 xfg ; 7=xyz type à la Babel
>4==> 40= no velocities; 41 with velocities
6==> 60= no velocities; 61 with velocitiesvitesses

**lwgin**=.false.  =true écrit un fichier .newgin à la fin  
iteangle = -1  formatted angular position  

lperiod=.true. si conditions periodiques (ipbc=1) et lperiod: coordonnées réduites entre 0 et 1. Si ipbc=1 et .not.lperiod coordonées peuvent dépasser
lprteat=.false.if you want to print the energy on atom
lprtsigat=.false.  calul et affichage de la contrainte sur chaque atome
lsigatcel=.false.  calul et affichage de la contrainte atomique moyenne sur la cellule
lax=.false.  stockage positions initiales
lprteattotm=.false.energie par atome totale (pot+cin) moyenne
ngrid = 20000   taille de la grille des potentiels
tempdeplainit=-1
debyetemp=-1
lprtpot=.false.
ibrake =0 if =1 electronic slowing for cascades (acting on all atoms)



##NEB
nebtype=2 NEB is the default  
fpstop=0.01 !stopping energy criterion (eV/ang)
lPathFromGin = .FALSE.  if T : read initial path in gin files *.1.gin, *.2.gin, ... (NEB calculaion)  
npath = 15 Number of  images of the neb is the default (including the two extremes that won't be moved see NEBHowto.md)  
maxneb = 700  the MAX of NEB steps  
neb_noise_scale=0.001   this will affect the 4th digit  
mdcg_noise_scale=0.001   this will affect the 4th digit  x + x*neb_noise_scale*random,  where "random" is a random number between  0 and 1  
neb_noise=0  0 without noise, 1 with noise  
mdcg_noise=0 0 without noise, 1 with noise  
lperiod=.false. pas de conditions periodiques  VERY IMPORTANT !  
**lclimb=.false.** ! set to true for climbin NEB  
nwclimb=3 ! starts the climbing at the second evaluation of forces (in VASP =1, in Henkelmann is set to "a few iterations")  
i_neb_drag =5 makes 5 in hyperplan (drag like) relaxation between succesive calculation of the spring forces. 1 maybe a better choice.  
nebrelaxation=2   We relax all the atoms; if nebrelaxation==1 only the most "deplaced" atoms  
deltaRmax=1.d-2  tthreshold for relaxation
iteanaposneb=0




ngrdel=500

tpseuils(:)=0  1:Tmin; 2:abs(T') ; 3: abs(T'') ; 1:abs(P); 2:abs(P') ; 3: abs(P'')

lrctest=.true.
tcelec=0  temp鬧ｻature de coupure pour les pertes 鬧脇ctroniques
Ecelec=0  temp鬧ｻature de coupure pour les pertes 鬧脇ctroniques


l2T=.false.  2T model
depmaxts=0.02
tsmin=2.0

units_lammps='TO_BE_SPECIFIED'

lvzeroneb=.false.  si true , met vp à 0 ente chaque iteration neb (comportement pre ndm2020), defaut = false==> calcul plus rapide
lcimb=.false. set to true for climbing image NEB

lparapath=.false.
nparapath=1
itypcalc=1
pas_lambda_mc = -100 valeur negative par defaut pour que l'utilisateur la change
n_path = -100 valeur negative par defaut pour que l'utilisateur la change
distminat=1  distance minimale en Angstrom de l'atome inséré aux autres atomes en Monte-Carlo (défaut = pas de distance min=n'importe où)
nbatplus=1
ipbc(1:3)=1  1=PBC; 2=wall... dimension 3 =plans bc; ac;ab

ndecal=2  nombre de décalage dans le calcul de la matrice de force (dmtype=19)
decal=0.1  décalage dans le calcul de la matrice de force (dmtype=19) (Angstroms)
lparafm=.true.
nparafm=nprocs
lwfm=.false.
lwritefreq=.true.


ndir=50   nombre de direction dans steepest descent
nstep=50   nombre de pas dans la minimisation sur une ligne en steepes descent
ncgtry=10   nombre de relaxations positions/celulle/positions/celluel,etc.
betaguess=1d-6
lvarstop=.false.
fstpdecr=10.
beta35=1d-10
gammas=0.999
gammav=0.9
idirectionmcgc=-2   direction pour le montecarlo 0 ou 1 a designer par l'utilisateur
lbiais_retrait = .false.biais ou non sur les retraits dans le montecarlo
lbiais_inser = .false.biais ou non sur les retraits dans le montecarlo
fdmc_1 = -1000.0param de fermi dirac A DEF PAR UTILISATEUR pour la fct discriminante du biais dans MC
fdmc_2 = -1000.0valeur devant etre changee

protocl_mcc= "mcp" ! mcp = linear icnrease of the insertion, "cos"= non linear (1-cos ) insertion; (1-sin) desinsertion  
fdfactmcgc=18.0
R0mcgc=-1.0
bublcenter(:)=0.5
izlins=-1
zlcenter=0.5
ins_typ=0
typswitch1=0
typswitch2=0
ldecalcor=.true.
kmin=0.  min kinetic energy for arps
kmax=0.  max kinetic energy for arps
lpartarps=.false.
noxyzkmin(1:3)=-1
noxyzkmax(1:3)=1000000
iteprtkin=-1

lspecialinit=.false.  driver for specail initialization : cascade, press or heat burst etc.

