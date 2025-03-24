module montecarlo_mod
  USE arret_ndm_mod,only:arret_ndm
  USE gen_com_m,only:  lperiod, tstep, timel, tstep,  itetabvois,lenfnam,&
       & iterasmol,itetemp, temp, kine, pi, bk, Text, gamlg,gamprfact,one,pi,text,tinit,&
       &lspaceNDM,rang,iteration,firsttime_lammps,erg2ev,fnam,fnamcout,unitP,dmtype,&
       &lrestartmcgc,imm_glob,iseed,sig,lprahman,sigext,h0,kcell,ucell,ihbox0,sigtot,sigkine
  USE atomconfig,only:atom_config,atom_config_d, switch_atom
  USE cellconfig, only:cell_config, caltabtC
  USE var_pot,only:ntyp,cm,gamlt
  USE calfo_mod,only: calfo
  USE T_kind_param_m, ONLY:  double
  USE cryst_to_cart_mod, ONLY: cryst_to_cart
  USE boxconfig,only:box_config,periodbox,box_config_lpr,updatebox
  USE rasmolT_mod,only: rasmolT
  use sigkinetot_mod,only:sigkinetotMC
  USE scalebox_mod,only: scalebox
  USE tempinstT_mod,only: tempinstT
  USE sauvegardeT_mod,only:sauvegardeT
  use paraconfig,only:para_config,commconstr,initparapuresp
  USE Parrinello_Rahman,only:initlpr
  USE init_simple_mod,only:init_simple
   USE calerf_mod

#ifdef PARA
  use Tpara,only:grp_world,nprocs,myidsp,MPI_COMM_space,nprocspace,ierr,mpi_comm_world,&
       &NDM_MPI_REAL_DOUBLE,para_space_config,status,comm_space,mpi_world
  USE init_vois_mod,only: init_voisinage
#else
  use Tpara,only:myidsp,nprocspace,para_space_config
#endif
#ifdef ML
  use NDM_ML,only:init_config_ml
#endif

  use read_val,only:rvois,ltabvois
  use var_pot,only:ipotentiel,rumax
  USE parautils,only:initloc,pointer_caltabt_calfo
  USE calctemp_mod,only:calctemp
  use vect_dist_mod,only:distat,closest_at
  USE recips_mod,only:distmin
#ifdef LAMMPS_VERSION
  use vars_lammps
  use lammps_util_mod,only:init_lammps
#endif  
  use config2data_mod,only:config2data
  USE constrconf_mod,only:read_cin,lprt
  use probMC,only:probMC1
  use Parrinello_Rahman,only:sp,sdot, sdot_new,trh0,invh0,invtrh0,epsi,tension,volu0,invvolu0
  implicit none

  type, extends (atom_config_d):: atom_config_mc
     real(double), allocatable :: proba_des(:) !defini pr chaque atome mais utile que pour O dans notre cas
     real(double) :: proba_ins !defini par POSITION, donc pas un tableau
  contains
     procedure, pass :: init => init_atom_config_mc
     procedure, pass :: copy_atom => copy_atom_mc
     procedure, pass :: copy_config => copy_config_mc
     procedure, pass :: switch_atom => switch_atom_mc
  end type atom_config_mc

  type(para_config)::parapath ! division de tous les procs en nparapath chemins calculés simultanément
  type(para_config)::paramcgc ! division de parapath en 2*espace
  type(para_space_config)::pscgc
  real(double)::distminat
  type(atom_config_mc),pointer::atconf_n !type derive atom_config du systeme a n atomes
  type(atom_config_mc),pointer::atconf_nplus1 !type derive atom_config du systeme a n+1 atomes
  type(cell_config),pointer:: cells_n !type derive cell_config du systeme a n atomes
  type(cell_config),pointer:: cells_nplus1 !type derive cell_config du systeme a n+1 atomes

  type(atom_config_mc),allocatable,target::config_atom_n(:) !type derive atom_config du systeme a n atomes
  type(atom_config_mc),allocatable,target::config_atom_nplus1(:) !type derive atom_config du systeme a n+1 atomes
  type(cell_config),allocatable,target:: config_cells_n(:) !type derive cell_config du systeme a n atomes
  type(cell_config),allocatable,target:: config_cells_nplus1(:) !type derive cell_config du systeme a n+1 atomes
  type(cell_config):: config_cells_old !type derive cell_config du systeme a n+1 atomes

  type(atom_config_mc):: config_atom_old_0, config_atom_old_1, config_atom_new_0, config_atom_new_1 !config intermediaire pour suivre l'evolution des systemes: 0 -> syst N, 1 -> syst N+1.


  type(box_config_lpr)::boxmcgc , box_new1, box_old1, box_new0, box_old0 !(les 2 dernurs définis uniquement pour bigmasters
  type(box_config_lpr),allocatable,target::boxmcgcpath(:)
  type(box_config_lpr),pointer::boxmcgc_p

  !LPR  
  !  real(double), allocatable :: sp(:,:), sdot(:,:), sdot_new(:,:),sfp(:,:),spp(:,:)
  !  real(double), dimension(3,3) ::trh0,invh0,invtrh0,epsi
  !  real(double),dimension(3,3)::grsig,tension

  integer::nbatplus
  integer::typswitch1,typswitch2
  integer::  pas_lambda_mc,idirectionmcgc
  real(double) :: lambda_mc !lambda compris entre 0 et 1

  integer :: n_path ! nb de chemin d'insertion, a definir dans .din, par defaut 10

  real(double), dimension(3,3) :: sig_n, sig_nplus1,grsig
  real(double) :: potist_n, potist_nplus1
  real(double) :: Weff, Work,tempcell
  logical::lmaster !(true= master du système N ou du système N+1)
  logical::lbigmaster !(true= master d'un/du calcul de chemin(s))
  ! il y a 2*plus de masters que de bigmaster(s)
  logical::lmegamaster !(true= master de pilotage des calculs de chemins).
  !Si PARA ET nparapath >1 ET lparapath alors il existe un megamatser et nparapath bigmasters. Dans les autres cas lmegamaster=lbigmasters)
  class(atom_config),pointer::atmcgcloc
  type(cell_config),pointer::cellmcgcloc
  type(cell_config),target:: cellcible ! ne sert qu'à faire pointer cellnebloc sur quelquechose
  type(atom_config_d),target::atcible
  integer::rgcib,rgem !(cible et emeteur e confN+1)
  integer, dimension(:), allocatable :: iseedt
  integer::seed_size
  integer:: nparapath
  logical ::  lparapath
  logical :: lbiais(0:1),lbiais_retrait,lbiais_inser
  real(double)::fdmc_1, fdmc_2 !paramtres pilotant la fct_alpha utilisee dans le biais des retraits:forme fermi dirac
  integer::itypcalc

  integer,parameter::nrins=10000
  real(double)::R0mcgc,fdfactmcgc,probaR(0:nrins),bublcenter(3),zlmin,Frad(1:nrins),fecalcprob
  real(double)::ZLcenter(3)
  integer::iZLins
  real(double)::epotnp1min=1d12,beta
  integer::ins_typ
  real(double),allocatable::rcpath(:)
  real(double)::k_spring,FEspring

  logical lspring

contains 

  subroutine init_montecarlo(boxndm,rv)
    class(box_config)::boxndm
    real(double),intent(in)::rv
    logical::linitpot,lcalc
    integer::ipp,imdm=0,ic
    real(double)::fe2,res,xerf

    beta = 1.0/(bk*Text)    
    call random_seed(size=seed_size)
    allocate(iseedt(seed_size))
    iseedt(:)=iseed
    call    random_seed (put=iseedt)


    allocate (config_atom_n(nparapath))
    allocate (config_atom_nplus1(nparapath))
    allocate (config_cells_n(nparapath))
    allocate (config_cells_nplus1(nparapath))
    allocate(boxmcgcpath(nparapath))
    !       if (nparapath==1) then
    select type (boxndm)
    type is (box_config_lpr)
       boxmcgc=boxndm
    end select
    select case (ins_typ)
    case(1,3,33,44,55,11)
       allocate(rcpath(nparapath))
       rcpath=0
    end select
    lprt=.true.
    do ipp=1,nparapath
       if (ipp.gt.1) lprt=.false.
       lcalc=.false.
       if (lparapath) then
          if (parapath%image+1==ipp) lcalc=.true.
       else
          lcalc=.true.
       end if
       lcalc=.true.
       if (lcalc) then

          boxmcgc_p=>boxmcgcpath(ipp)
          atconf_n=> config_atom_n(ipp)
          cells_n=>config_cells_n(ipp)
          atconf_nplus1=>config_atom_nplus1(ipp)
          cells_nplus1=>config_cells_nplus1(ipp)
          boxmcgc_p=boxmcgc
          if (idirectionmcgc==0) then
             call atconf_n%init(imdm,imm_glob,ltabvois,rvois=rv)
          else
             call atconf_nplus1%init(imdm,imm_glob,ltabvois,rvois=rv)
          end if
          ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs

          if (ipp==1) then
             linitpot=.true.
          else
             linitpot=.false.
          end if
          if (idirectionmcgc==0) then
             call init_simple(atconf_n%atom_config_d,cells_n,boxmcgc_p,psc=pscgc,linitpot=linitpot)
             select case(ins_typ)
             case(1,3,33,44,11,55)
                call init_instyp
             end select
                
             call initNP1(ipp) ! initialise la configuration N+1

          else
             call init_simple(atconf_nplus1%atom_config_d,cells_nplus1,boxmcgc_p,psc=pscgc,linitpot=linitpot)
             select case(ins_typ)
             case(1,3,33,44,11,55)
                call init_instyp
             end select
             call initN(ipp) ! initialise la configuration N+1
          end if
          if (lprahman) then
             call initlpr(atconf_nplus1,cells_nplus1,boxmcgc_p,pscgc)
             !          call initMClpr(atconf_nplus1%im)
          end if
       end if
    end do
    !    if (lbigmaster) write(10+rang,*)rcpath
    !END PARAPATH
    atconf_n=> config_atom_n(1)
    cells_n=>config_cells_n(1)
    atconf_nplus1=>config_atom_nplus1(1)
    cells_nplus1=>config_cells_nplus1(1)

    if (rang==0) then
       write(6,*)'***************PATH MONTE-CARLO*****************'
       write(6,'(A,I6,A,I6,A,I4)')'pas_lambda=',pas_lambda_mc,' npath=',n_path,' naparapath=',nparapath
       write(6,'(A,I3,A)')'ins_typ=',ins_typ, ' (0=random; 1=sph 2=switch type, 3=slice, 11 site+spring, 33 slice +spring, 44 line+spring), 55 sphere +spring'
       write(6,*)'lbiais_retrait , lbiais_inser ',lbiais_retrait,lbiais_inser
       if (lbiais_inser) then
          select case(ins_typ)
          case (1)
             if (rang==0) write(6,*)'R0mcgc bublcenter ',R0mcgc,bublcenter
          case(3)
             if (rang==0) write(6,'(A,G15.7,I3,G15.7)')'R0mcgc, izlins(X,Y,Z), zlcenter ',R0mcgc,izlins, zlcenter(izlins)
          end select
       end if
    end if



    if (lspring) then
       fe2=0
       select case(ins_typ)
       case(11)
          FEspring=(bk*text*log(boxmcgc_p%volu) -(3*bk*text/2)*log(2*pi*bk*text/k_spring))*erg2ev
          xerf=0.5*boxmcgc_p%at(1,1)/dsqrt(2*pi*bk*text/k_spring)
          call calerf(xerf,res,0)
          fe2= erg2ev*(-3*bk*text*log(dsqrt(2*pi*bk*text/k_spring)*res/boxmcgc_p%zl(1)))
       case(33)
          FEspring=(bk*text*log(boxmcgc_p%zl(izlins))-(bk*text/2)*log(2*pi*bk*text/k_spring))*erg2ev
       case(44)
          Fespring=0
          do ic=1,3
             if (ic.ne.izlins) then 
                FEspring=Fespring+(bk*text*log(boxmcgc_p%zl(ic)))
             end if
          end do
          FEspring=Fespring-(bk*text)*log(2*pi*bk*text/k_spring)
          FEspring=Fespring*erg2eV
       case(55)
          FEspring=Fecalcprob ! bk*text*log(boxmcgc_p%volu)*erg2ev
       end select
       if (rang==0) then
          write(6,*)'***SPRING CALCULATION***'
          write(6,'(A)')'Free energy of the spring to ADD to the calculated chamical potential at the very end (in eV) (second value is better if non zero)'
          write(6,*)'FEspring=',FEspring,fe2
       end if
    end if

          
  end subroutine init_montecarlo

  subroutine montecarlo

    implicit none
    logical :: lextend

    integer :: direction
    integer :: n_accepted, n_accepted_0, n_accepted_1
    integer :: n_gen, n_gen_0, n_gen_1
    real(double) :: acceptance_rate, acceptance_rate_0, acceptance_rate_1

    integer :: i_path,ipp,ipch
    real(double)::zr1,xp_np1(3)

    integer :: iloc
    integer :: acceptation, test_acc

    real(double) :: W, Wprec, Wprecedent !sauvegarde Wprec pour posttraitement
    real(double) :: mu_moy, mu_wrmc, mu_NC, mu_SC,tempf,pressf,kindum

    real(double),allocatable:: Weff_npp(:), pot_npp(:)

    logical :: lchange,ldistrib,lcalc
    CHARACTER(len=89) :: fnamread
    real(double), dimension(2,22):: pot_cumul ! pour le calcul du pot chimique


    real(double) :: travail_prec,qeff
    integer :: dir_prec
    !########################################################################################################################
    !                                             Initialisation
    !########################################################################################################################

    !initialisation variables 
    !pour le premier chemin: sens positif, d'ajout d'une particule et acceptation
    !    write(6,*)'RANG UISEED',rang,iseed
    lbiais(0)=lbiais_inser
    lbiais(1)=lbiais_retrait
    direction=idirectionmcgc
#ifdef PARA
    cellmcgcloc=>cellcible
    atmcgcloc=>atcible
    lmaster=paramcgc%lmaster
    do ipp=1,nparapath
       if ((paramcgc%mpi_image%nproc.gt.1).and.(lspaceNDM.eqv..true.)) then
          call init_voisinage(config_cells_n(ipp),pscgc) !  a faire seulement une fois ?
       end if
    end do

#else
    lmaster=.true.
#endif


    lextend = .true.
    lperiod = .true.

    n_accepted   = 0
    n_accepted_0 = 0
    n_accepted_1 = 0
    n_gen   = 0
    n_gen_0 = 0
    n_gen_1 = 0
    acceptance_rate   = 0.0
    acceptance_rate_0 = 0.0
    acceptance_rate_1 = 0.0

    mu_moy = 0.0 
    mu_wrmc = 0.0
    mu_NC = 0.0
    mu_SC = 0.0

    test_acc = 0    


    pot_cumul(:,:) = 0.0
    if (nparapath.gt.0) then
       allocate (Weff_npp(nparapath))
       allocate (pot_npp(nparapath))
       Weff_npp(:)=0
       pot_npp(:)=0
    end if
    if (lmegamaster) then
       if (dmtype==151)open(UNIT= 754, FILE="analyse_file.151")!, STATUS = 'new')
    end if
    !  if (rang==0) write(6,*)' IN MCGC nbatplus, idirection',nbatplus,idirectionmcgc
    if (lrestartmcgc) then
       if (lmegamaster) then
          write(*,*) 'imm_n, imm_nplus1', config_atom_n(1)%imm, config_atom_nplus1(1)%imm
          write(*,*) 'im_n, im_nplus1', config_atom_n(1)%im, config_atom_nplus1(1)%im
          call config_atom_old_0%init(config_atom_n(1)%im, config_atom_n(1)%imm, config_atom_n(1)%ltabvois,&
               &im_glob=config_atom_n(1)%im_glob,imm_glob=imm_glob)
          call config_atom_old_1%init(config_atom_nplus1(1)%im, config_atom_n(1)%imm,config_atom_nplus1(1)%ltabvois,&
               &im_glob=config_atom_nplus1(1)%im_glob,imm_glob=imm_glob)
          fnamread= fnam(1:lenfnam)//'.N.cout'
          !  write(6,*)'R1',config_atom_n(1)%imm
          call read_cin(box_old0,1,config_atom_old_0%atom_config_d,config_atom_n(1)%imm,fnamread) ! 1=complet
          fnamread= fnam(1:lenfnam)//'.NP1.cout'
          ! write(6,*)'R2',config_atom_n(1)%imm
          call read_cin(box_old1,1,config_atom_old_1%atom_config_d,config_atom_n(1)%imm,fnamread) ! 1=complet
          !seul MEGAMASTER A LES POSITIONS OLD
          call config_atom_old_0%copy_config(config_atom_n(1), lrescl=.true.)          
          call config_atom_old_1%copy_config(config_atom_nplus1(1), lrescl=.true.) 
          call caltabtC(config_cells_n(1),config_atom_n(1),lperiod,box_old0)
          call caltabtC(config_cells_nplus1(1),config_atom_nplus1(1),lperiod,box_old1)
          if(dmtype==151) then
             if(direction==0) then
                config_cells_old=config_cells_n(1)
             else
                config_cells_old=config_cells_nplus1(1)
             end if
          end if
          !recalculer les probas du systeme
          atconf_nplus1=>config_atom_nplus1(1)
          if (direction == 0) then
             boxmcgc=box_old0
             boxmcgcpath(1)=box_old0
          else
             boxmcgcpath(1)=box_old1
             boxmcgc=box_old1
          end if
          boxmcgc_p=>boxmcgcpath(1)
          call calcul_proba_des
          xp_np1(:)=atconf_nplus1%xp(:,atconf_nplus1%im)
          atconf_nplus1%proba_ins= calcul_proba_ins (xp_np1)

          !recopier les nouvelles configs dans old 1
          call config_atom_nplus1(1)%copy_config(config_atom_old_1,lrescl=.true.)
          call restart_chemin(travail_prec, dir_prec)
          Wprec = travail_prec
          direction = 1 - dir_prec !si on avait 0 on repart de 1 et si on avait 1 on repart de 0.
          if (direction == 0) then
             call analyse_montecarlo(config_atom_n(1),config_cells_n(1),boxmcgcpath(1),'SystN_init')
          else
             call analyse_montecarlo(config_atom_nplus1(1),config_cells_nplus1(1),boxmcgcpath(1),'SystNP1_init')
          end if !analyse_montecarlo
       end if ! if megamaster
       !#ifdef PARA
       !       if (lprahman) call boxmcgc%master2slave(0,mpi_world)
       !#endif

       !il s'agit d'envoyer a tous les procs la direction, le Wprec et les config old et new
       if (lbigmaster) then

          call parapath%mpi_master%bcast(0,Wprec)
          call parapath%mpi_master%bcast(0,direction)

#ifdef PARA

          if (lparapath) then
             if (lprahman) call boxmcgcpath(1)%master2slave(0,parapath%mpi_master)
             call config_atom_n(1)%atom_config_d%send2all(0,parapath%mpi_master)
             call config_atom_nplus1(1)%atom_config_d%send2all(0,parapath%mpi_master)

             call config_cells_n(1)%send2all(0,parapath%mpi_master)
             call config_cells_nplus1(1)%send2all(0,parapath%mpi_master)
          end if
#endif

          do ipp=2,nparapath
             config_cells_n(ipp)= config_cells_n(1)
             config_cells_nplus1(ipp)= config_cells_nplus1(1)
             config_atom_n(ipp)=config_atom_n(1)
             config_atom_nplus1(ipp)=config_atom_nplus1(1)
             boxmcgcpath(ipp)=boxmcgcpath(1)
          end do


          ipch = 1 !utile pour la suite de la boucle sur les chemins
       end if !bigmaster

       !recalculer les forces
       call lambda(direction, nstep = 0, protocol_name = 'MCP') !initialisation dulambda a 0 pour le premier melange des forces
       iloc=1;lchange=.false.;ldistrib=.true.
       call calfoMCGC(iloc,lchange,ldistrib)
       ! a la fin de lrestart, tous les procs ont N et N+1 courants pareil + old0 et old1 sont connus + Wprec + direction 
       ! de meme que les proba ont ete calculees

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!***********************!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!***********************!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!***********************!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else !cas où lrestart = .false. 

       !deplacé !
       if (lbigmaster) then
          ! sauvegarde du système
          if (idirectionmcgc == 0) then
             boxmcgc_p=>boxmcgcpath(1)
             call config_atom_n(1)%copy_config(config_atom_old_0, lrescl=.true.)
             box_old0=boxmcgcpath(1)

          else
             box_old1=boxmcgcpath(1)
             do ipp=1,nparapath
                atconf_nplus1=>config_atom_nplus1(ipp)
                boxmcgc_p=>boxmcgcpath(ipp)

                call calcul_proba_des ! on initialise une désintégration qui va être accepté (car c'est la première) : Il faut calculer les proba pour les mettre dans OLD etdans config_atom_nplus1
                call atconf_nplus1%copy_config(config_atom_old_1, lrescl=.true.)
             end do
          end if
       end if
       if (idirectionmcgc == 0) then

          call analyse_montecarlo(config_atom_n(1),config_cells_n(1),boxmcgcpath(1),'SystN_init')
       else
          call analyse_montecarlo(config_atom_nplus1(1),config_cells_nplus1(1),boxmcgcpath(1),'SystNP1_init')
       end if
       if(dmtype==151) then
          if(direction==0) then
             config_cells_old=config_cells_n(1)
          else
             config_cells_old=config_cells_nplus1(1)
          end if


       end if

       do ipp=1,nparapath
          lcalc=.false.
          if (lparapath) then
             if (parapath%image+1==ipp) then
                lcalc=.true.
             end if
          else
             lcalc=.true.
          end if
          if (lcalc) then
             !          write(6,*)'CALC1',rang,ipp
             atconf_n=> config_atom_n(ipp)
             cells_n=>config_cells_n(ipp)
             atconf_nplus1=>config_atom_nplus1(ipp)
             cells_nplus1=>config_cells_nplus1(ipp)
             boxmcgc_p=>boxmcgcpath(ipp)


             direction = idirectionmcgc ! direction = 0 on ajoute un atome, = 1 on retire un atome

             call lambda(direction, nstep = 0, protocol_name = 'MCP') !initialisation du lambda a 0 pour le premier melange des forces
             iloc=1;lchange=.false.;ldistrib=.true.
             call calfoMCGC(iloc,lchange,ldistrib)
!             write(6,*)'MCCDBG1 ', rang
             !pour le premier chemin: sens positif, d'ajout d'une particule et acceptation

             if (lprahman) then
                call langevinLPR(direction, protocol = 'MCP')
             else
                call langevin(direction, protocol = 'MCP',qeff=qeff,work=work)
             end if
             !if (lbigmaster) write(6,'(A,I2,3G15.7)')'potist', direction,Weff*erg2ev,Qeff*erg2ev,work*erg2ev
             ! if LPR call langevinPc
             if (idirectionmcgc == 0) then 
                weff_npp(ipp)= +Weff
                pot_npp(ipp)= potist_nplus1
             else
                weff_npp(ipp)= -Weff
                pot_npp(ipp)= potist_nplus1
             end if
             if (lbigmaster) then

                if (idirectionmcgc==0) then
                   tempf=tempinstt(atconf_nplus1,kindum,latcomp=.true.)
                else
                   tempf=tempinstt(atconf_n,kindum,latcomp=.true.)
                end if
                pressF= (sigtot(1,1)+sigtot(2,2)+sigtot(3,3))/3.0
                select case(ins_typ)
                   case(1,3,33,44,11,55)
                      
                      !                   call parapath%mpi_master%sum(rcpath)
                      write(6,*)'Weff eV dist', ipp,weff_npp(ipp)*erg2eV,rcpath(ipp)
                      if (dmtype==151) write(754,*)'Weff eV dist', ipp,weff_npp(ipp)*erg2eV,rcpath(ipp)
                      !                   write(6,*)'Weff eV dist', ipp,weff_npp(ipp)*erg2eV,rcpath(ipp)
                      
                   case default
                      
                      if (dmtype==151) write(754,*)'Weff eV dist', ipp,weff_npp(ipp)*erg2eV
                      write(6,*)'Weff eV ', ipp,weff_npp(ipp)*erg2eV
                   end select
                end if
             !if (lbigmaster)write(6,*)'potist', ipp,potist_n,potist_nplus1
          end if

       end do
       if (dmtype==15) then
          if ((lbigmaster).and.(lparapath)) then
             call parapath%mpi_master%sum(weff_npp)
             call parapath%mpi_master%sum(pot_npp)
          end if

          if (lbigmaster) then 
             if (lmegamaster) then
                if (nparapath.gt.1) then
                   call random_number(zr1)

                   ipch=1+int(nparapath*zr1) ! choix aléatoire débile
                   write(6,*)'chemin choisi aléatoirement',ipch
                else
                   ipch=1
                end if
             end if

             call parapath%mpi_master%bcast(0,ipch)

          Weff=weff_npp(ipch) !on a effectué le chnment de signe si retrait

#ifdef PARA
          if (lparapath) then 
             call config_atom_n(ipch)%send2all(ipch-1,parapath%mpi_master)
             call config_atom_nplus1(ipch)%send2all(ipch-1,parapath%mpi_master)
             call config_cells_n(ipch)%send2all(ipch-1,parapath%mpi_master)
             call config_cells_nplus1(ipch)%send2all(ipch-1,parapath%mpi_master)
             call boxmcgcpath(ipch)%master2slave(ipch-1,parapath%mpi_master)
          end if


#endif          

          atconf_nplus1=>config_atom_nplus1(ipch)
          boxmcgc_p=>boxmcgcpath(ipch)
          if (idirectionmcgc == 0) then
             if (lmegamaster) then
                epotnp1min=pot_npp(ipch)
                fnamcout = fnam(1:lenfnam)//'.NP1min.cout'
                write(6,*)'new epotnp1min ', epotnp1min*erg2ev
                call sauvegardeT(config_atom_nplus1(ipch),config_cells_nplus1(ipch),boxmcgcpath(ipch),3,fnamcout,latcomp=.true.)
             end if
             call calcul_proba_des ! on vient de choisir ipch qui est accepté. On calcule les proba pour : 1:choisir les atomes à désintégrer et mettre dans old_1 pour les calculs du biais
             call config_atom_nplus1(ipch)%copy_config(config_atom_old_1, lrescl=.true.)
             box_old1=boxmcgcpath(ipch)
          else
             xp_np1(:)=atconf_nplus1%xp(:,atconf_nplus1%im)
             atconf_nplus1%proba_ins= calcul_proba_ins (xp_np1)

             call config_atom_n(ipch)%copy_config(config_atom_old_0, lrescl=.true.)
             box_old0=boxmcgcpath(ipch)
          end if
          do ipp=1,nparapath
             config_cells_n(ipp)= config_cells_n(ipch)
             config_cells_nplus1(ipp)= config_cells_nplus1(ipch)
             config_atom_n(ipp)=config_atom_n(ipch)
             config_atom_nplus1(ipp)=config_atom_nplus1(ipch)
             boxmcgcpath(ipp)=boxmcgcpath(ipch)
          end do

          W = WEff !avec le bon signe
          !W = Work
          Wprec = + W
          Wprecedent = Wprec
          if (lmegamaster) write(*,*) 'W0',W*erg2eV

       end if! sur bigmaster

       if (lmaster) then ! on est dans l'un des 2 masters

          call boxmcgcpath(1)%master2slave(0,paramcgc%mpi_master)
          boxmcgcpath(:)=boxmcgcpath(1)
          rgcib=1;rgem=0
          if(paramcgc%image==0) then !on est dans le master général
             call atconf_Nplus1%send2proc(rgcib,paramcgc%mpi_master,'x')
          else !on est dans le master de N+1
             call atconf_Nplus1%recv(rgem,paramcgc%mpi_master,'x')
          end if

       end if


       direction = 1 - idirectionmcgc !=0 si le premier pas était un retrait, =1 si le premier pas
              else if (dmtype==151) then

          if (lbigmaster) then
             do ipp=1,nparapath
                if (direction==0) then 
                   config_cells_n(ipp)= config_cells_old
                   config_atom_n(ipp)=config_atom_old_0
                   boxmcgcpath(ipp)=box_old0
                else
                   config_cells_nplus1(ipp)= config_cells_old
                   config_atom_nplus1(ipp)=config_atom_old_1
                   boxmcgcpath(ipp)=box_old1
                   !             config_atom_nplus1(ipp)=config_atom_nplus1(ipch)
                end if
             end do
          end if
       end if

    end if !if sur lrestart

    if (lmegamaster) then
       if (dmtype==15) then
          open(UNIT= 752, FILE="analyse_file", STATUS = 'new')

          open(UNIT= 85, FILE="restart_file", STATUS = 'new')
          !      open(UNIT= 753, FILE="nrj_pot_systacc", STATUS = 'new')
          if (nparapath .gt. 1) then
             write(752,*) '#ACC/REF  direction  ipchemin  WeV(x nparapath)&
                  & Wprec XPROB(x nparapath +1)'
          else
             write(752,*) '#ACC/REF  direction  WeV  Wprec  XPROB   XALEA'
          end if! sur nparapath
       else if (dmtype==151) then 
          open(UNIT= 754, FILE="analyse_file.151")!, STATUS = 'new')
       end if
    end if!sur megamaster

    !########################################################################################################################
    !                                             boucle sur lambda le long d'un chemin
    !########################################################################################################################

    DO i_path = 1, n_path ! boucle à faire pour tous les procs
#ifdef PARA
       call mpi_world%barrier
#endif
       Weff_npp(:)=0
       pot_npp(:)=0
       if (lmegamaster)then
          write(6,*)
          write(6,*)'path ',i_path,' in direction', direction, ' to ',1-direction
       end if

       if ((lbigmaster).and.(dmtype==15)) then
          do ipp=1,nparapath
             config_atom_nplus1(ipp)%vp(:,:)   = - config_atom_nplus1(ipp)%vp(:,:) !à chaque retour dans la boucle, on change de direction
             config_atom_n(ipp)%vp(:,:)   = - config_atom_n(ipp)%vp(:,:)
             if (lprahman) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!hdotpath(ipp,:,:)=-hdotpath(ipp,:,:)
             endif
          end do
          if (direction == 0) then
             call config_atom_n(ipch)%copy_config(config_atom_new_0, lrescl=.true.)
             box_new0=boxmcgcpath(ipch)
             boxmcgcpath(:)=boxmcgcpath(ipch)
          endif
          if (direction == 1) then
             call config_atom_nplus1(ipch)%copy_config(config_atom_new_1, lrescl=.true.)
             box_new1=boxmcgcpath(ipch)
             boxmcgcpath(:)=boxmcgcpath(ipch)
          endif

       end if !if lbigmaster
       do ipp=1,nparapath
          !write(200+rang,*)'IPPP2',ipp
          !call boxmcgcpath(ipp)%print(unit=200+rang)
#ifdef PARA
          if (lparapath) then
             call boxmcgc_p%master2slave(0,paramcgc%mpi_master)
          end if
#endif
          lcalc=.false.
          if (lparapath) then
             if (parapath%image+1==ipp) lcalc=.true.
          else
             lcalc=.true.
          end if

          if (lcalc) then

             !             if (lbigmaster)write(6,*)'parapath',rang,i_path,ipp
             atconf_n=> config_atom_n(ipp)
             cells_n=>config_cells_n(ipp)
             atconf_nplus1=>config_atom_nplus1(ipp)
             cells_nplus1=>config_cells_nplus1(ipp)
             boxmcgc_p=>boxmcgcpath(ipp)

             !choisir l'at a retirer ou ajouter + preparation des syst N et N+1 pour etre prets pour le langevin (cad decoupage cellules + calcul forces + melange des forces - se fait dans cette sous routine)
             select case (ins_typ)
             case(0,1,3,33,44,55,11)
                call ajout_retrait(direction,ipp)
             case(2)
                call type_switch(direction)
             end select
!             write(6,*)'MCCDBG2 ', rang
             ! pas de langevin
             if (lprahman) then
                call langevinLPR(direction, protocol = 'MCP')
             else
                call langevin(direction, protocol = 'MCP',qeff=qeff,work=work)
             end if
             ! if LPR call langevinPc
             !if (lbigmaster) write(6,'(A,I2,3G15.7)')'potist', direction,Weff*erg2ev,Qeff*erg2ev,work*erg2ev
             if (direction == 0) then
                Weff_npp(ipp)= +weff
                pot_npp(ipp)= potist_nplus1
             else
                Weff_npp(ipp)= -weff
                pot_npp(ipp)= potist_nplus1
             end if !sur direction
             if (lbigmaster) then
                if (direction==0) then
                   tempf=tempinstt(atconf_nplus1,kindum,latcomp=.true.)
                else
                   tempf=tempinstt(atconf_n,kindum,latcomp=.true.)
                end if
                pressF= (sigtot(1,1)+sigtot(2,2)+sigtot(3,3))/3.0
                !                write(6,*)'Weff eV Tempf PressF',ipp, weff_npp(ipp)*erg2eV,tempf,pressf*unitP
             end if

          end if
!!$          call mpi_finalize(ierr)
!!$          call arret_ndm
       end do !boucle nparapath
       select case(ins_typ)
       case(1,3,33,44,11,55)
          if (lbigmaster)then
             call parapath%mpi_master%sum(rcpath)
          end if
       end select
       if ((lbigmaster).and.(lparapath)) then
          call parapath%mpi_master%sum(weff_npp)
          call parapath%mpi_master%sum(pot_npp)

       end if
       select case(ins_typ)
       case(1,3,33,44,11,55)
          if (lmegamaster) then
             do ipp=1,nparapath
                write(6,'(A,I2,I4,2F20.10)')'Weff eV dist ', direction,ipp,weff_npp(ipp)*erg2eV,rcpath(ipp)
                if (dmtype==151) write(754,*)'Weff eV dist', ipp,weff_npp(ipp)*erg2eV,rcpath(ipp)
             end do
          end if
          rcpath=0
       case default
          if (lmegamaster) then
             do ipp=1,nparapath
                write(6,'(A,I2,I4,F20.10)')'Weff eV',direction, ipp,weff_npp(ipp)*erg2eV
                if (dmtype==151) write(754,*)'Weff eV dist', ipp,weff_npp(ipp)*erg2eV
             end do
          end if
       end select
       n_gen = n_gen + 1
       if (direction == 0) then
          n_gen_0 = n_gen_0 + 1
       else
          n_gen_1 = n_gen_1 + 1
       end if
       if (dmtype==15) then
!!!!!!!!!!!!!!!!!!!!!! TEST D'ACCEPTATION !!!!!!!!!!!!!!!!!!!
          if (nparapath.gt.1) then 
             call multiproposal(Weff_npp, pot_npp, Wprec, Wprecedent, ipch, direction, acceptation, test_acc, n_gen, &
                  & lbiais, mu_moy, mu_wrmc, mu_NC, mu_SC, pot_cumul)
          else
             call monoproposal(Weff_npp, pot_npp, Wprec, Wprecedent, ipch, direction, acceptation, test_acc, n_gen, &
                  & lbiais, mu_moy, mu_wrmc, mu_NC, mu_SC, pot_cumul)
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

          if (acceptation == 1) then
             n_accepted = n_accepted + 1
             !choupi
             if (direction == 0) then
                n_accepted_0 = n_accepted_0 + 1
             else 
                n_accepted_1 = n_accepted_1 + 1
             end if
          end if

          Wprecedent = Wprec
          iteration = iteration + 1

          if (direction == 0) then
             direction = 1
          else
             direction = 0
          end if

          acceptance_rate   = (real(n_accepted)/real(n_gen))*1.0d2
          acceptance_rate_0 = (real(n_accepted_0)/real(n_gen_0))*1.0d2
          acceptance_rate_1 = (real(n_accepted_1)/real(n_gen_1))*1.0d2
       else if (dmtype==151) then
          if (lbigmaster) then 
             do ipp=1,nparapath
                if (direction==0) then 
                   config_cells_n(ipp)= config_cells_old
                   config_atom_n(ipp)=config_atom_old_0
                   boxmcgcpath(ipp)=box_old0
                else
                   config_cells_nplus1(ipp)= config_cells_old
                   config_atom_nplus1(ipp)=config_atom_old_1
                   boxmcgcpath(ipp)=box_old1
                   !             config_atom_nplus1(ipp)=config_atom_nplus1(ipch)
                end if
             end do
          end if

       end if
    END DO !end do sur la boucle des chemins

    ! write(6,*)'BARRIERE FINALE ',rang
#ifdef PARA
    call MPI_BARRIER(parapath%mpi_orig%comm,ierr)
#endif

    if (lmegamaster) then
       write(*,*) ' taux d acceptation final   : ', acceptance_rate,  ' %'
       write(*,*) ' taux d acceptation alpha 0 : ', acceptance_rate_0,' %'
       write(*,*) ' taux d acceptation alpha 1 : ', acceptance_rate_1,' %'
       write(*,'(A10, G15.7,A10, G15.7, A10, G15.7, A10, G15.7)') &
            &'mu_moy',mu_moy ,'mu_wrmc', mu_wrmc, 'mu_NC', mu_NC, 'mu_SC', mu_SC
       close(752)
       close(85)
       !      close(753)
    end if
    !stop
       if ((rang==0).and.(lspring)) then
          write(6,*)'***SPRING CALCULATION***'
          write(6,*)'Free energy of the spring to ADD to the calculated chamical potential at the very end (in eV)'
          write(6,*)'FEspring=',FEspring
       end if

    
  end subroutine montecarlo


  subroutine monoproposal(travail_npp, nrjpot_npp, Wprece, Wpreced, ipchemin, dir, accepta, premier_accept, &
       &ngen, lbiais, pot_moy, pot_wrmc, pot_NC, pot_SC,tab_cumul)

    implicit none 

    real(double), dimension(nparapath) :: travail_npp, nrjpot_npp
    real(double) :: Wprece ! Wprec
    real(double) :: Wpreced !sauvegarde Wprec pour posttraitement
    integer :: ipchemin, dir, accepta, ngen
    logical :: lbiais(0:1)
    integer ::  premier_accept
    real(double) :: pot_moy, pot_wrmc, pot_NC, pot_SC


    real(double) :: biais
    real(double) :: W, xprob, xalea
    real(double) :: ln_xalea, ln_Wprec, ln_W, ln_xprob
    real(double) :: theta, xp_np1(3)
    real(double), dimension(nparapath+1) :: xprob_i ! naparapath = 1 dans le cas du monoproposal
    logical :: lextend, lperiod 
    real(double), dimension(2,22) :: tab_cumul

    theta = 0.5

    lextend = .true.
    lperiod = .true.
    xprob_i(:) = 0.0
    !CRC
    if (lbigmaster) then !master general
       if (lmegamaster) then ! megamaster et bigmaster sont indetiques ici je pense. Ce if ne sert à rien
          ipchemin=1
       end if ! mega master   
       Weff=travail_npp(ipchemin)  
       W = Weff ! le signe a ete inversé precedemment

       if (dir == 0) then
          !W = +Work
          atconf_nplus1=>config_atom_nplus1(ipchemin)
          boxmcgc_p=>boxmcgcpath(ipchemin)
          call calcul_proba_des ! calcul_proba placé étrangement. Devrait etre après l'acceptation. Mais ça marche car si acceptation alors va devenir la nouvelle conf N+1 a tester et va devenir _old_1 (ici copié en _new_1. Après acceptation : _new_1 copié en _old_1). 
          call config_atom_nplus1(ipchemin)%copy_config(config_atom_new_1, lrescl=.true.)
          box_new1=boxmcgcpath(ipchemin)
       else
          !W = - Work

          xp_np1(:)=config_atom_old_1%xp(:,config_atom_old_1%im)
          config_atom_old_1%proba_ins= calcul_proba_ins (xp_np1)

          call config_atom_n(ipchemin)%copy_config(config_atom_new_0, lrescl=.true.)

          box_new0=boxmcgcpath(ipchemin)
       endif


       if (dir == 0) then
          if (lbiais(0)) then
             biais = config_atom_nplus1(ipchemin)%proba_ins&
                  &/config_atom_old_1%proba_ins
          else
             biais = 1.0
          end if
       else
          if (lbiais(1)) then
             biais = config_atom_nplus1(ipchemin)%proba_des(config_atom_nplus1(ipchemin)%im)&
                  &/config_atom_old_1%proba_des(config_atom_nplus1(ipchemin)%im)
          else
             biais=1
          end if
       end if

       ln_Wprec  = (+beta*(dir-theta)*Wprece)
       ln_W      = (+beta*(dir-theta)*W)
       ln_xprob  = - dlog(1 + (dexp(ln_Wprec-ln_W)*biais)) !avec biais
       !ln_xprob  = - dlog(1 + dexp(ln_Wprec-ln_W))   !sans biais
       xprob     = dexp(ln_xprob)
       xprob_i(1) = xprob   !proba d'accpeter W
       xprob_i(2) = 1.0 - xprob !proba d'accepter Wprec
       call random_number(xalea)
       ln_xalea  = log(xalea)

       if (lmegamaster) write(*,'(A10, G25.16E3,A10, G25.16E3,A10, G25.16E3,A10, G25.16E3)')&
            & 'WeV', W*erg2eV, 'WpreceV', Wprece*erg2eV,'XPROB', xprob, 'XALEA', xalea
       !if (lmegamaster) write(752,'(I3, 4G25.16E3)')  dir,  W*erg2eV, Wprece*erg2eV, xprob, xalea 

       if (ln_xprob > ln_xalea) then    

!!!!!!!!!!!!!!!!!!!!!!!!!! ACCEPTATION   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (lmegamaster) write(*,*) 'ACCEPTATION, direction=', dir,'W', W*erg2eV, &
               &'Wprec', Wprece*erg2eV, '  LN_XPROB ', ln_xprob, '  XPROB ', xprob, '  XALEA ', xalea
          if (lmegamaster) write(752,'(I3, I3, 4G25.16E3)') 1 , dir,  W*erg2eV,Wprece*erg2eV, xprob, xalea
          if (lmegamaster) write(85,'(I3,G25.16E3)') dir, W
          premier_accept = 1
          Wprece = + W 
          accepta = 1

!!!!! etape 1 pot chimique !!!!!!!!!
          call potentiel_chimique(tab_cumul, pot_moy, pot_wrmc, pot_NC, pot_SC, premier_accept, ngen, ipchemin,&
               & 1, travail_npp, xprob_i, Wpreced)

          call config_atom_new_0%copy_config(config_atom_old_0, lrescl=.true.)
          call config_atom_new_1%copy_config(config_atom_old_1, lrescl=.true.)
          box_old0=box_new0
          box_old1=box_new1
          if(lmegamaster) then
             fnamcout = fnam(1:lenfnam)//'.N.cout'
             call sauvegardeT(config_atom_new_0%atom_config_d,cells_n,box_new0,3,fnamcout,latcomp=.true.)
             fnamcout = fnam(1:lenfnam)//'.NP1.cout'
             call sauvegardeT(config_atom_new_1%atom_config_d,cells_nplus1,box_new1,3,fnamcout,latcomp=.true.)
          end if

          if (lmegamaster) then
             if (idirectionmcgc == 1) then 
                call analyse_montecarlo(config_atom_n(ipchemin),config_cells_n(ipchemin),boxmcgcpath(ipchemin),'SystN_accepte')
             else
                if (epotnp1min.gt.nrjpot_npp(ipchemin)) then
                   epotnp1min=nrjpot_npp(ipchemin)
                   fnamcout = fnam(1:lenfnam)//'.NP1min.cout'
                   write(6,*)'new epotnp1min ', epotnp1min*erg2ev
                   call sauvegardeT(config_atom_new_1%atom_config_d,config_cells_nplus1(ipchemin),&
                        &box_new1,3,fnamcout,latcomp=.true.)
                end if

                call analyse_montecarlo(config_atom_nplus1(ipchemin),&
                     &config_cells_nplus1(ipchemin),boxmcgcpath(ipchemin),'SystNP1_accepte')
             end if
             ! Ecriture de l'energie potentiel du chemin accepte
             !            write(753, '(2G25.16E3)') nrjpot_npp(ipchemin), kine
          end if
       else
!!!!!!!!!!!!!!!!!!!!!!!!!! REFUS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (lmegamaster) write(*,*) ' REJECTION, direction=', dir, 'W', W*erg2eV, &
               &'Wprec', Wprece*erg2eV, '  LN_XPROB ', ln_xprob, '  XPROB ', xprob, '  XALEA ', xalea
          if (lmegamaster) write(752,'(I3, I3, 4G25.16E3)') 0 ,  dir,  W*erg2eV, Wprece*erg2eV, xprob, xalea

          !on accepte le sens opposé - changer des signes des vitesses 
          config_atom_old_0%vp(:,:)   = - config_atom_old_0%vp(:,:)
          config_atom_old_1%vp(:,:)   = - config_atom_old_1%vp(:,:)
!!!!!!!!!!!!!!!!!!!!!!          hdotCHANGE

          accepta = 0
          if (dir == 0) then
             call config_atom_old_1%copy_config(config_atom_nplus1(ipchemin), lrescl=.true.)
             boxmcgcpath(ipchemin)=box_old1
             call caltabtC(config_cells_nplus1(ipchemin),config_atom_nplus1(ipchemin),lperiod,boxmcgcpath(ipchemin))
          else
             call config_atom_old_0%copy_config(config_atom_n(ipchemin), lrescl=.true.)
             boxmcgcpath(ipchemin)=box_old0
             call caltabtC(config_cells_n(ipchemin),config_atom_n(ipchemin),lperiod,boxmcgcpath(ipchemin))
          end if

!!!!! etape 2 pot chimique !!!!!!!!!
          call potentiel_chimique(tab_cumul, pot_moy, pot_wrmc, pot_NC, pot_SC, premier_accept, ngen, ipchemin,&
               & 2, travail_npp, xprob_i, Wpreced)

       end if !test sur xprob         

!!!!! etape 3 pot chimique !!!!!!!!!
       call potentiel_chimique(tab_cumul, pot_moy, pot_wrmc, pot_NC, pot_SC, premier_accept, ngen, ipchemin,&
            & 3, travail_npp, xprob_i, Wpreced)

    end if !fin master general

  end subroutine monoproposal


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multiproposal(travail_npp, nrjpot_npp, Wprece, Wpreced, ipchemin, dir, accepta, premier_accept, ngen,&
       & lbiais, pot_moy, pot_wrmc, pot_NC, pot_SC, tab_cumul)

    implicit none 

    real(double), dimension(nparapath) :: travail_npp, nrjpot_npp
    real(double) :: Wprece ! Wprec
    real(double) :: Wpreced !sauvegarde Wprec pour posttraitement
    integer :: ipchemin, dir, accepta, ngen,  premier_accept
    logical :: lbiais(0:1)
    real(double) :: pot_moy, pot_wrmc, pot_NC, pot_SC




    real(double) :: theta, beta,xp_np1(3)
    integer :: ipp
    real(double), dimension(nparapath+1) :: xprob_i
    logical :: lextend, lperiod
    real(double), dimension(2,22) :: tab_cumul

    theta = 0.5
    beta = 1.0/(bk*Text)
    lextend = .true.
    lperiod = .true.
    xprob_i(:) = 0.0
    if (lbigmaster) then !master general

       if (lmegamaster) then
          call calcul_chemin(travail_npp, xprob_i, Wprece, dir, theta, ipchemin, lbiais) 
          call choix_chemin(xprob_i, ipchemin)

          write(6,*)'chemin choisi',ipchemin
       end if !megamaster

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ACCEPTATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if (ipchemin.le.nparapath) then ! si ipchemin <= nparapath: on choisit une des configs generees
          call parapath%mpi_master%bcast(0,ipchemin)

          if (lmegamaster) write(*,'(A20, 1G25.16E3, A10, 1G25.16E3, A10, 1G25.16E3)') 'ACCEPTATION WeV',&
               & (travail_npp(ipchemin)*erg2eV), 'WpreceV', (Wprece*erg2eV),'XPROB', xprob_i(ipchemin)
          if (lmegamaster) write(752,'(3I3, 50G25.16E3)') 1 , &
               &dir, ipchemin, (travail_npp*erg2eV), (Wprece*erg2eV), xprob_i
          if (lmegamaster) write(85,'(1I3,G25.16E3)') dir, travail_npp(ipchemin)
#ifdef PARA
          if (lparapath) then 
             call config_atom_n(ipchemin)%send2all(ipchemin-1,parapath%mpi_master)
             call config_atom_nplus1(ipchemin)%send2all(ipchemin-1,parapath%mpi_master)
             call config_cells_n(ipchemin)%send2all(ipchemin-1,parapath%mpi_master)
             call config_cells_nplus1(ipchemin)%send2all(ipchemin-1,parapath%mpi_master)
             call boxmcgcpath(ipchemin)%master2slave(ipchemin-1,parapath%mpi_master)
          end if
#endif          


          if (dir == 0) then
             atconf_nplus1=>config_atom_nplus1(ipchemin)
             boxmcgc_p=>boxmcgcpath(ipchemin)
             call calcul_proba_des ! on sait que ce atconf_nplu1 est le bon et accepté. Calcul de proba_des pour choix future quand désitégration et biais futurs copié en _new_1
             call config_atom_nplus1(ipchemin)%copy_config(config_atom_new_1, lrescl=.true.)
             box_new1=boxmcgcpath(ipchemin)
          else
             xp_np1(:)=config_atom_nplus1(ipchemin)%xp(:,config_atom_nplus1(ipchemin)%im)
             config_atom_new_1%proba_ins= calcul_proba_ins (xp_np1)
             call config_atom_n(ipchemin)%copy_config(config_atom_new_0, lrescl=.true.)
             box_new0=boxmcgcpath(ipchemin)
          end if
          do ipp=1,nparapath
             config_cells_n(ipp)= config_cells_n(ipchemin)
             config_cells_nplus1(ipp)= config_cells_nplus1(ipchemin)
             config_atom_n(ipp)=config_atom_n(ipchemin)
             config_atom_nplus1(ipp)=config_atom_nplus1(ipchemin)
             boxmcgcpath(ipp)=boxmcgcpath(ipchemin)
          end do

          Weff=travail_npp(ipchemin) ! il y a deja eu le changement de signe
          Wprece = Weff
          accepta = 1
          premier_accept = 1

!!!!! etape 1 pot chimique !!!!!!!!!
          call potentiel_chimique(tab_cumul, pot_moy, pot_wrmc, pot_NC, pot_SC, premier_accept, ngen, ipchemin,&
               & 1, travail_npp, xprob_i, Wpreced)


          call config_atom_new_0%copy_config(config_atom_old_0, lrescl=.true.)
          call config_atom_new_1%copy_config(config_atom_old_1, lrescl=.true.)
          box_old0=box_new0
          box_old1=box_new1
          if(lmegamaster) then
             fnamcout = fnam(1:lenfnam)//'.N.cout'
             call sauvegardeT(config_atom_new_0%atom_config_d,config_cells_n(ipchemin),box_new0,3,fnamcout,latcomp=.true.)
             fnamcout = fnam(1:lenfnam)//'.NP1.cout'
             call sauvegardeT(config_atom_new_1%atom_config_d,config_cells_nplus1(ipchemin),box_new1,3,fnamcout,latcomp=.true.)
             !CRC CHECK LES cells...
          end if

          if (lmegamaster) then
             if (idirectionmcgc == 1) then
                call analyse_montecarlo(config_atom_n(ipchemin),config_cells_n(ipchemin),boxmcgcpath(ipchemin),'SystN_accepte')
             else
                if (epotnp1min.gt.nrjpot_npp(ipchemin)) then
                   epotnp1min=nrjpot_npp(ipchemin)
                   fnamcout = fnam(1:lenfnam)//'.NP1min.cout'
                   write(6,*)'new epotnp1min ', epotnp1min*erg2ev
                   call sauvegardeT(config_atom_new_1%atom_config_d&
                        &,config_cells_nplus1(ipchemin),box_new1,3,fnamcout,latcomp=.true.)
                end if
                call analyse_montecarlo(config_atom_nplus1(ipchemin),config_cells_nplus1(ipchemin)&
                     &,boxmcgcpath(ipchemin),'SystNP1_accepte')
             end if
             ! Ecriture de l'energie potentiel du chemin accepte
             !            write(753, '(2G25.16E3)') nrjpot_npp(ipchemin), kine 
          end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   REFUS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       else ! si ipchemin > nparapath: on choisit la config prec

          if (lmegamaster) write(*,'(A15, 1G25.16E3, A10, 1G25.16E3, A10, 1G25.16E3)') 'REFUS WeV', &
               &(Wprece*erg2eV), 'WpreceV', (Wprece*erg2eV),'XPROB', xprob_i(ipchemin)
          if (lmegamaster) write(752,'(3I3, 50G25.16E3)')  0 , &
               &dir, ipchemin, travail_npp*erg2eV, (Wprece*erg2eV), xprob_i
          config_atom_old_0%vp(:,:)   = - config_atom_old_0%vp(:,:)
          config_atom_old_1%vp(:,:)   = - config_atom_old_1%vp(:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          hdot_old=-hdot_old

          accepta = 0
          ipchemin = 1 ! là où on recopie la configuration precedente, par defaut on ecrase le chemin 1
          call parapath%mpi_master%bcast(0,ipchemin) !on envoie le nouveau ipchemin a tous les procs, là où chemin prec va etre mis
          if (dir == 0) then
              
             call config_atom_old_1%copy_config(config_atom_nplus1(ipchemin), lrescl=.true.)
             boxmcgcpath(ipchemin)=box_old1
             call caltabtC(config_cells_nplus1(ipchemin),config_atom_nplus1(ipchemin),lperiod,boxmcgcpath(ipchemin))
          else
             call config_atom_old_0%copy_config(config_atom_n(ipchemin), lrescl=.true.)
             boxmcgcpath(ipchemin)=box_old0
             call caltabtC(config_cells_n(ipchemin),config_atom_n(ipchemin),lperiod,boxmcgcpath(ipchemin))
          end if
          !on envoie l'ancienne conf a tous les procs
#ifdef PARA
          if (lparapath) then 
             call config_atom_n(ipchemin)%send2all(ipchemin-1,parapath%mpi_master)
             call config_atom_nplus1(ipchemin)%send2all(ipchemin-1,parapath%mpi_master)
             call config_cells_n(ipchemin)%send2all(ipchemin-1,parapath%mpi_master)
             call config_cells_nplus1(ipchemin)%send2all(ipchemin-1,parapath%mpi_master)
             call boxmcgcpath(ipchemin)%master2slave(ipchemin-1,parapath%mpi_master)
          end if
#endif          

          !pour chaque procs, on copie l'ancienne conf dans tous les chemins          
          do ipp=1,nparapath
             config_cells_n(ipp)= config_cells_n(ipchemin)
             config_cells_nplus1(ipp)= config_cells_nplus1(ipchemin)
             config_atom_n(ipp)=config_atom_n(ipchemin)
             config_atom_nplus1(ipp)=config_atom_nplus1(ipchemin)
             boxmcgcpath(ipp)=boxmcgcpath(ipchemin)
          end do

!!!!! etape 2 pot chimique !!!!!!!!!
          call potentiel_chimique(tab_cumul, pot_moy, pot_wrmc, pot_NC, pot_SC, premier_accept, ngen, ipchemin,&
               & 2, travail_npp, xprob_i, Wpreced)

       end if !test sur ipchemin. acceptation ou refus (W ou Wprec)

!!!!! etape 3 pot chimique !!!!!!!!!
       call potentiel_chimique(tab_cumul, pot_moy, pot_wrmc, pot_NC, pot_SC, premier_accept, ngen, ipchemin,&
            & 3, travail_npp, xprob_i, Wpreced)

    end if !fin master general

  end subroutine multiproposal


  subroutine potentiel_chimique(cumul, moy, wrmc, NC, SC, prem_accept, nb_gen, chemin, etape,&
       & liste_travail, liste_proba, Wprecedent)

    implicit none
    real(double), dimension(nparapath) :: liste_travail
    real(double), dimension(nparapath+1) :: liste_proba
    real(double), dimension(2,22) :: cumul
    real(double) :: moy, wrmc, NC, SC, Wprecedent
    integer :: prem_accept, etape, nb_gen, chemin
    !etape correspond a l'etape du pot a laquelle on est:
    !etape == 1: on est dans l'acceptation
    !etape == 2: on est dans le refus
    !etape == 3: on est apres le test 

    real(double) :: beta, beta_eV

    real(double),dimension(2) :: contribut_accepte, f2, fminusf
    real(double),dimension(2) :: f_cumul, f2_cumul, fminusf_cumul
    real(double),dimension(2) :: f_wr, f2_wr, fminusf_wr
    real(double),dimension(2) :: f_wr_cumul, f2_wr_cumul, fminusf_wr_cumul
    real(double),dimension(2) :: f_cumul_inte, f2_cumul_inte, fminusf_cumul_inte
    real(double),dimension(2) :: f_wr_cumul_inte, f2_wr_cumul_inte, fminusf_wr_cumul_inte

    real(double),dimension(2) :: b_opt, b_wr_opt, est_opt, est_opt_bwr
    integer :: direc, ip
    real(double), dimension(nparapath+1) :: liste_travaux ! contient les Weff des chemins + Wprec en dernier element

    if (lmegamaster) then
       beta = 1.0/(bk*Text)
       beta_eV = 1.0/((8.617333262145E-5)*Text)

       liste_travaux(1:nparapath) = liste_travail(1:nparapath)
       liste_travaux(nparapath+1) =  Wprecedent

       contribut_accepte(1:2) = cumul(:,1)
       f2(1:2) = cumul(:,2)
       fminusf(1:2) = cumul(:,3)
       f_wr(1:2) = cumul(:,4)
       f2_wr(1:2) = cumul(:,5)
       fminusf_wr(1:2) = cumul(:,6)
       f_cumul(1:2) = cumul(:,7)
       f2_cumul(1:2) = cumul(:,8)
       fminusf_cumul(1:2) = cumul(:,9)
       f_wr_cumul(1:2) = cumul(:,10)
       f2_wr_cumul(1:2) = cumul(:,11)
       fminusf_wr_cumul(1:2) = cumul(:,12)
       b_opt(1:2) = cumul(:,13)
       b_wr_opt(1:2) = cumul(:,14) 
       est_opt(1:2) = cumul(:,15)
       est_opt_bwr(1:2) = cumul(:,16)
       f_cumul_inte(1:2) = cumul(:,17)
       f2_cumul_inte(1:2) = cumul(:,18)
       fminusf_cumul_inte(1:2) = cumul(:,19)
       f_wr_cumul_inte(1:2) = cumul(:,20)
       f2_wr_cumul_inte(1:2) = cumul(:,21)
       fminusf_wr_cumul_inte(1:2) = cumul(:,22)

       if (etape == 1) then
          contribut_accepte(1) = dexp(beta*liste_travaux(chemin)*0.5)
          f2(1) = dexp(beta*liste_travaux(chemin))
          fminusf(1) =(dexp(beta*liste_travaux(chemin)*0.5) - dexp(beta*liste_travaux(nparapath+1)*0.5))**2 

          contribut_accepte(2) = dexp(-beta*liste_travaux(chemin)*0.5)
          f2(2) = dexp(-beta*liste_travaux(chemin))
          fminusf(2) =(dexp(-beta*liste_travaux(chemin)*0.5) - dexp(-beta*liste_travaux(nparapath+1)*0.5))**2
          !write(*,*) 'acceptation', contribut_accepte 
          cumul(:,1) = contribut_accepte(1:2) 
          cumul(:,2) = f2(1:2) 
          cumul(:,3) = fminusf(1:2) 
          cumul(:,4) = f_wr(1:2) 
          cumul(:,5) = f2_wr(1:2) 
          cumul(:,6) = fminusf_wr(1:2) 
          cumul(:,7) = f_cumul(1:2) 
          cumul(:,8) = f2_cumul(1:2) 
          cumul(:,9) = fminusf_cumul(1:2) 
          cumul(:,10) = f_wr_cumul(1:2) 
          cumul(:,11) = f2_wr_cumul(1:2) 
          cumul(:,12) = fminusf_wr_cumul(1:2) 
          cumul(:,13) = b_opt(1:2) 
          cumul(:,14) = b_wr_opt(1:2) 
          cumul(:,15) = est_opt(1:2) 
          cumul(:,16) = est_opt_bwr(1:2) 
          cumul(:,17) = f_cumul_inte(1:2) 
          cumul(:,18) = f2_cumul_inte(1:2) 
          cumul(:,19) = fminusf_cumul_inte(1:2) 
          cumul(:,20) = f_wr_cumul_inte(1:2) 
          cumul(:,21) = f2_wr_cumul_inte(1:2) 
          cumul(:,22) = fminusf_wr_cumul_inte(1:2)
       end if

       if (etape == 2) then
          contribut_accepte(1) = dexp(beta*liste_travaux(nparapath+1)*0.5)
          f2(1) = dexp(beta*liste_travaux(nparapath+1))
          fminusf(1) = 0.0

          contribut_accepte(2) = dexp(-beta*liste_travaux(nparapath+1)*0.5)
          f2(2) = dexp(-beta*liste_travaux(nparapath+1))
          fminusf(2) = 0.0

          !write(*,*) 'refus', contribut_accepte, cumul(1:2,1)
          cumul(:,1) = contribut_accepte(1:2) 
          cumul(:,2) = f2(1:2) 
          cumul(:,3) = fminusf(1:2) 
          cumul(:,4) = f_wr(1:2) 
          cumul(:,5) = f2_wr(1:2) 
          cumul(:,6) = fminusf_wr(1:2) 
          cumul(:,7) = f_cumul(1:2) 
          cumul(:,8) = f2_cumul(1:2) 
          cumul(:,9) = fminusf_cumul(1:2) 
          cumul(:,10) = f_wr_cumul(1:2) 
          cumul(:,11) = f2_wr_cumul(1:2) 
          cumul(:,12) = fminusf_wr_cumul(1:2) 
          cumul(:,13) = b_opt(1:2) 
          cumul(:,14) = b_wr_opt(1:2) 
          cumul(:,15) = est_opt(1:2) 
          cumul(:,16) = est_opt_bwr(1:2) 
          cumul(:,17) = f_cumul_inte(1:2) 
          cumul(:,18) = f2_cumul_inte(1:2) 
          cumul(:,19) = fminusf_cumul_inte(1:2) 
          cumul(:,20) = f_wr_cumul_inte(1:2) 
          cumul(:,21) = f2_wr_cumul_inte(1:2) 
          cumul(:,22) = fminusf_wr_cumul_inte(1:2)
       end if

       if (etape == 3) then
!!!!! Calcul estimateur WR/SC !!!!!!!!!
          f_wr(:) = 0.0
          f2_wr(:) = 0.0
          fminusf_wr(:) = 0.0

          DO ip = 1, nparapath+1
             if (liste_proba(ip) .ne. 0.0) then !permet d'eviter lorsqu'un W trop grand est genere qu'il fasse devenir NaN le potchi
                f_wr(1) = f_wr(1) + liste_proba(ip)*dexp(beta*liste_travaux(ip)*0.5) 
                f2_wr(1) = f2_wr(1) + liste_proba(ip)*dexp(beta*liste_travaux(ip))

                f_wr(2) = f_wr(2) + liste_proba(ip)*dexp(-beta*liste_travaux(ip)*0.5) 
                f2_wr(2) = f2_wr(2) + liste_proba(ip)*dexp(-beta*liste_travaux(ip))
                !if (lmegamaster) write(*,*) 'proba' ,liste_proba(ip), 'WeV', liste_travaux(ip)*erg2eV, 'W', liste_travaux(ip)

                fminusf_wr(1) = fminusf_wr(1) + liste_proba(ip)*(dexp(beta*liste_travaux(ip)*0.5)&
                     &-dexp(beta*liste_travaux(nparapath+1)*0.5))**2
                fminusf_wr(2) = fminusf_wr(2) + liste_proba(ip)*(dexp(-beta*liste_travaux(ip)*0.5)&
                     &-dexp(-beta*liste_travaux(nparapath+1)*0.5))**2
             end if
          END DO
          !if (lmegamaster) write(*,*) 'f_wr', f_wr, 'f2_wr', f2_wr, 'fminusf_wr', fminusf_wr

!!!!! formation des sommes !!!!
          DO direc = 1,2
             f_cumul(direc) = f_cumul(direc) + contribut_accepte(direc)
             f2_cumul(direc) = f2_cumul(direc) + f2(direc)
             fminusf_cumul(direc) = fminusf_cumul(direc) + fminusf(direc)

             f_wr_cumul(direc) = f_wr_cumul(direc) + f_wr(direc)
             f2_wr_cumul(direc) = f2_wr_cumul(direc) + f2_wr(direc)
             fminusf_wr_cumul(direc) = fminusf_wr_cumul(direc) + fminusf_wr(direc)


!!!!! diviser par le nombre de chemin!!!!!
             f_cumul_inte(direc) = f_cumul(direc) / real(nb_gen)
             f2_cumul_inte(direc) = f2_cumul(direc) / real(nb_gen)
             fminusf_cumul_inte(direc) = fminusf_cumul(direc) / (real(nb_gen*2.0))

             f_wr_cumul_inte(direc) = f_wr_cumul(direc) / real(nb_gen)
             f2_wr_cumul_inte(direc) = f2_wr_cumul(direc) / real(nb_gen)
             fminusf_wr_cumul_inte(direc) = fminusf_wr_cumul(direc) / real(nb_gen*2.0)

!!!!! calcul de la variable de contrôle !!!!
             if (prem_accept == 1) then
                b_opt(direc)=(f2_cumul_inte(direc)-f_cumul_inte(direc)**2)&
                     & / fminusf_cumul_inte(direc)
                b_wr_opt(direc)=(f2_wr_cumul_inte(direc)-f_wr_cumul_inte(direc)**2)&
                     & / fminusf_wr_cumul_inte(direc)


!!!!! calcul de l'estimateur NC et SC !!!!!!
                est_opt(direc) = b_opt(direc)*f_wr_cumul_inte(direc) &
                     & + (1.0-b_opt(direc))*f_cumul_inte(direc)
                est_opt_bwr(direc) = b_wr_opt(direc)*f_wr_cumul_inte(direc) &
                     & + (1.0-b_wr_opt(direc))*f_cumul_inte(direc)
             end if
          END DO

          !write(*,*) 'f_cumul', f_cumul, 'f2_cumul', f2_cumul, 'fminusf_cumul', fminusf_cumul
          !write(*,*) 'f_wr_cumul', f_wr_cumul, 'f2_wr_cumul', f2_wr_cumul, 'fminusf_wr_cumul', fminusf_wr_cumul
          !write(*,*) 'b_opt', b_opt, 'b_wr_opt', b_wr_opt, 'est_opt',   est_opt, 'est_wr_opt', est_opt_bwr
          !write(*,*) 'direc b_opt et b_wr_opt', dir, est_opt(1), est_opt_bwr(1)

!!!! calcul de mu!!!!!
          if (prem_accept == 1) then
             moy = -(1/beta_eV)*dlog(f_cumul_inte(2) / f_cumul_inte(1)) !estimateur simple Im(f)
             wrmc = -(1/beta_eV)*dlog(f_wr_cumul_inte(2) / f_wr_cumul_inte(1)) !estimateur WRMC (moyenne des Wgen pondérée par les proba)
             NC = -(1/beta_eV)*dlog(est_opt(2) / est_opt(1))  !estimateur NC Jnc,M(f)
             SC = -(1/beta_eV)*dlog(est_opt_bwr(2) / est_opt_bwr(1)) !estimateur SC Jsc,M(f)
          end if
!!!!!!!!!! fin calcul mu !!!!!!!!!!!
          cumul(:,1) = contribut_accepte(1:2) 
          cumul(:,2) = f2(1:2) 
          cumul(:,3) = fminusf(1:2) 
          cumul(:,4) = f_wr(1:2) 
          cumul(:,5) = f2_wr(1:2) 
          cumul(:,6) = fminusf_wr(1:2) 
          cumul(:,7) = f_cumul(1:2) 
          cumul(:,8) = f2_cumul(1:2) 
          cumul(:,9) = fminusf_cumul(1:2) 
          cumul(:,10) = f_wr_cumul(1:2) 
          cumul(:,11) = f2_wr_cumul(1:2) 
          cumul(:,12) = fminusf_wr_cumul(1:2) 
          cumul(:,13) = b_opt(1:2) 
          cumul(:,14) = b_wr_opt(1:2) 
          cumul(:,15) = est_opt(1:2) 
          cumul(:,16) = est_opt_bwr(1:2) 
          cumul(:,17) = f_cumul_inte(1:2) 
          cumul(:,18) = f2_cumul_inte(1:2) 
          cumul(:,19) = fminusf_cumul_inte(1:2) 
          cumul(:,20) = f_wr_cumul_inte(1:2) 
          cumul(:,21) = f2_wr_cumul_inte(1:2) 
          cumul(:,22) = fminusf_wr_cumul_inte(1:2) 
       end if

    end if !if lmegamaster
  end subroutine potentiel_chimique


  subroutine ajout_retrait(direc,ipp)

    implicit none

    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    integer,intent(in) :: direc ,ipp
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    real(double), allocatable, dimension(:,:) :: cart_vec_nplus1
    integer,allocatable :: indice(:)
    integer:: i,nag,rgcib,rgem,iloc,i1,i2,j,iplus
    logical ::ldistrib,lchange
    real(double)::pins,poscenter(3,1),postest(3)
    allocate (cart_vec_nplus1(3,nbatplus))
    allocate(indice(nbatplus))

    lperiod = .true.

    !######################################### Direction 0 vers 1 (ajout) ##########################################

    if (direc == 0) then ! ajout d'une particule en N+1
       if (lbigmaster) then
          !tirer des positions aleatoires pour les N+nbatplus eme atome
          call atom_supp(cart_vec_nplus1,pins,ipp)
          !          do i=1,nbatplus
          !             write(*,'(A25, 3G25.16E3,A,I4)') 'atome supplementaire', cart_vec_nplus1(:,i), 'image',parapath%image+1
          !          end do
          call boucle_copy_atom(atconf_N,atconf_Nplus1, sens= .false.)   
          !addition de la n+1eme particule
          do i=1,nbatplus
             iplus=i+atconf_n%im
             atconf_nplus1%xp(1:3,iplus) = cart_vec_nplus1(1:3,i)
             atconf_nplus1%fp(1:3,iplus) = 0
             atconf_nplus1%ityp(iplus) = itypcalc
             atconf_nplus1%ielat(iplus) = -1
             nag=maxval(atconf_Nplus1%num_at_glob(1:iplus-1))
             atconf_Nplus1%num_at_glob(iplus) = nag+1
             atconf_Nplus1%proba_ins = pins
             !             write(6,*)'pinsN',i,pins
          end do
          call init_vitesse(atconf_nplus1,param = 0)
       end if
       !repreparer les config pour le prochain langevin

#ifdef PARA
       rgcib=1;rgem=0
       if(paramcgc%image==0) then !procs N
          if (lmaster)   call atconf_Nplus1%send2proc(rgcib,paramcgc%mpi_master)    
       else !procs N+1
          if (lmaster) call atconf_Nplus1%recv(rgem,paramcgc%mpi_master)
       end if
#endif

       iloc=1;lchange=.false.;ldistrib=.true.
       call calfoMCGC(iloc,lchange,ldistrib)

    end if

    !######################################### Direction 1 vers 0 (retrait) ##########################################


    if (direc == 1) then ! retrait d'une particule alea, la placer en N+1eme position, copier le syst pour le syst à N
       if (lbigmaster) then

          if (.not.lbiais(1)) then
             !SANS BIAIS
             call indice_alea(atconf_Nplus1,indice,itypcalc,nbatplus)
          else
             !AVEC BIAIS
             call atom_biais(atconf_Nplus1,indice)
          end if
          do i=1,nbatplus
             if (indice(i).gt.atconf_N%im) cycle
             do j=atconf_Nplus1%im,atconf_N%im+1,-1
                if (any(indice(1:nbatplus)==j)) then
                   cycle
                else
                   i1=indice(i)
                   i2=j
                   indice(i)=i2
                   call atconf_Nplus1%switch_atom(i1,i2)
                   exit
                end if
             end do
          end do
          select case(ins_typ)
          case(1)
             poscenter(:,1)=bublcenter(:)
             call cryst_to_cart(1,poscenter,boxmcgc_p%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m
             postest(:)=-1*(poscenter(:,1)-atconf_Nplus1%xp(:,atconf_Nplus1%im))
             rcpath(ipp)=sqrt(postest(1)**2+postest(2)**2+postest(3)**2)*1d8
          case(3,33)
             rcpath(ipp)=abs(atconf_Nplus1%xp(izlins,atconf_Nplus1%im)-zlcenter(izlins)*boxmcgc%at(izlins,izlins))*1d8
             
          end select
          call calcul_proba_des ! sans doute inutile
          !          do i=1,nbatplus
          !             iplus=atconf_N%im+i
          !             write(*,'(A25, 3G25.16E3, A10, G25.16E3, A10, I4 )') 'coord atome a retirer',  &
          !                  &atconf_Nplus1%xp(:,iplus),'proba', atconf_Nplus1%proba(iplus),' image ',parapath%image+1
          !on copie les N nouveaux premiers atomes du syst N+1 dans le systeme N
          !          end do
          call boucle_copy_atom(atconf_N,atconf_Nplus1, sens = .true.)  

       end if
#ifdef PARA

       rgcib=1;rgem=0
       if(paramcgc%image==0) then !procs N
          if (lmaster)   call atconf_Nplus1%send2proc(rgcib,paramcgc%mpi_master)
       else !procs N+1
          if (lmaster) call atconf_Nplus1%recv(rgem,paramcgc%mpi_master)
       end if
#endif
       iloc=1;lchange=.false.;ldistrib=.true.
       call calfoMCGC(iloc,lchange,ldistrib)      


    end if


  end subroutine ajout_retrait

  function  calcul_proba_ins(xpt) result (pinser)
    implicit none

    real(double), dimension(3)::xpt
    real (double)::pinser

    real(double)::poscenter(3,1),postest(3),rd,distzl

    select case (ins_typ)
    case(1)
       
       poscenter(:,1)=bublcenter(:)
       call cryst_to_cart(1,poscenter,boxmcgc_p%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m
       postest(:)=-1*(poscenter(:,1)-xpt(:))
       rd=sqrt(postest(1)**2+postest(2)**2+postest(3)**2)
       pinser= rd**2/(1+exp(fdfactmcgc*(rd-R0mcgc)))
    case(3)
       distzl=abs(ZLcenter(izlins)-dot_product(xpt,boxmcgc_p%as(:,izlins))/norm2(boxmcgc_p%as(:,izlins)))*norm2(boxmcgc_p%as(:,izlins))
       pinser=1/(1+exp(fdfactmcgc*(distzl-R0mcgc)))
       
    case default
       pinser=1
    end select
    !    write(6,*)'probainser',pinser
  end function calcul_proba_ins
!!!!!!!!!!!!!!!!!!!!!!

  subroutine calcul_proba_des
    implicit none
    real(double), dimension(atconf_Nplus1%imm) :: proba
    integer :: i
    real(double) ::   sum_norm

    !       DO i=1, atconf_Nplus1%im
    !          write(6,*)'AV',iteration,atconf_Nplus1%proba(1),atconf_Nplus1%proba(atconf_Nplus1%im)
    !       end DO
    if (lbiais(1)) then
       sum_norm = 0.0
       DO i=1, atconf_Nplus1%im
          !         write(6,*)'AV',iteration,atconf_Nplus1%proba(i)
          if (atconf_Nplus1%ityp(i) == itypcalc) then ! si l'atome est un oxygene
             call probMC1(atconf_Nplus1%xp(:,i),proba(i),boxmcgc_p%bg,boxmcgc_p%at,fdmc_1,fdmc_2)
             sum_norm=sum_norm+proba(i)
             !write(*,*) i, dist_tot, alpha
          else
             proba(i) = 0.0
          end if !si oxygene ou uranium

       END DO ! boucle atomes

    else
       sum_norm = 0.0
       DO i=1, atconf_Nplus1%im
          if (atconf_Nplus1%ityp(i) == itypcalc) then ! si l'atome est un oxygene
             proba(i)=1.
             sum_norm=sum_norm+1.
          else
             proba(i) = 0.0
          end if !si oxygene ou uranium

       END DO ! boucle atomes
    end if
    proba(:) = proba(:) / sum_norm
    atconf_Nplus1%proba_des(:) = proba(:)

  end subroutine calcul_proba_des



  subroutine calcul_chemin(travail, proba,Wp, dir, the, ind_chemin, lbiais)
    real(double), dimension(nparapath+1) :: proba
    real(double), dimension(nparapath) :: travail, biais 
    real(double) :: Wp, sum_expW, beta, the, lower, upper
    integer :: i, dir, ind_chemin
    logical :: lbiais(0:1)

    beta = 1.0/(bk*Text)
    sum_expW = 0.0
    lower = 1.0D-308
    upper = 1.0D308

    if (dir == 0) then
       if (lbiais(0)) then
          DO i = 1, nparapath
             biais(i) = config_atom_nplus1(i)%proba_ins&
                  &/config_atom_old_1%proba_ins
          end DO
       else
          biais(:) = 1.0
       end if
    else
       if (lbiais(1)) then
          DO i = 1, nparapath
             biais(i) = config_atom_nplus1(i)%proba_des(config_atom_nplus1(i)%im)&
                  &/config_atom_old_1%proba_des(config_atom_nplus1(ind_chemin)%im)
          end DO
       else
          biais(:)=1
       end if
    end if




    if ((exp(beta*(dir-the)*Wp) .lt. upper) .and. (exp(beta*(dir-the)*Wp) .gt. lower)) then
       !calcul de la somme des expW
       DO i = 1, nparapath
          if ((exp(beta*(dir-the)*(travail(i)-Wp)) .lt. upper) .and. (exp(beta*(dir-the)*(travail(i)-Wp)) .gt. lower)) then
             sum_expW = sum_expW + exp(beta*(dir-the)*(travail(i)-Wp))/biais(i)
          end if
       END DO

       !attribution d'une proba pour chaque chemin
       DO i = 1, nparapath
          if ((exp(beta*(dir-the)*(travail(i)-Wp)) .lt. upper) .and. (exp(beta*(dir-the)*(travail(i)-Wp)) .gt. lower)) then
             proba(i) = ( exp(beta*(dir-the)*(travail(i)-Wp))/biais(i) ) / (1.0 + sum_expW)
          else
             proba(i) = 0.0
          end if
       END DO
       !proba de tirer l'ancien chemin
       proba(nparapath+1) = 1.0 / (1.0 + sum_expW)


    else !dans le cas rare ou le premier chemin est si defavorable que exp(beta/2 Wprec)= 0, alors il ne faut pas passer par les soustractions W - Wp pour calculer la proba


       !calcul de la somme des expW
       DO i = 1, nparapath
          if ((exp(beta*(dir-the)*(travail(i))) .lt. upper) .and. (exp(beta*(dir-the)*(travail(i))) .gt. lower)) then
             if (dir==1) then
                if (lbiais(1)) then
                   sum_expW = sum_expW + exp(beta*(dir-the)*(travail(i)))/&
                        &(config_atom_nplus1(i)%proba_des(config_atom_nplus1(i)%im)) !1/alpha,new *exp(...)
                else
                   sum_expW = sum_expW + exp(beta*(dir-the)*(travail(i)))
                end if
             else
                if (lbiais(0)) then
                   sum_expW = sum_expW + exp(beta*(dir-the)*(travail(i)))/(config_atom_nplus1(i)%proba_ins)
                else
                   sum_expW = sum_expW + exp(beta*(dir-the)*(travail(i)))
                end if
             end if
          end if
       END DO

       !attribution d'une proba pour chaque chemin
       DO i = 1, nparapath
          if ((exp(beta*(dir-the)*(travail(i))) .lt. upper) .and. (exp(beta*(dir-the)*(travail(i))) .gt. lower)) then

             if (dir==1) then
                if (lbiais(1)) then
                   proba(i) = ( exp(beta*(dir-the)*(travail(i))) / &
                        &config_atom_nplus1(i)%proba_des(config_atom_nplus1(i)%im)  ) / ( sum_expW)
                else
                   proba(i) = ( exp(beta*(dir-the)*(travail(i))) ) / ( sum_expW)
                end if
             else
                if (lbiais(0)) then
                   proba(i) = ( exp(beta*(dir-the)*(travail(i))) / &
                        &config_atom_nplus1(i)%proba_ins) / ( sum_expW)

                else
                   proba(i) = ( exp(beta*(dir-the)*(travail(i))) ) / ( sum_expW)
                end if
             end if
          else
             proba(i) = 0.0
          end if
       END DO
       !proba de tirer l'ancien chemin
       proba(nparapath+1) = 0.0

    end if

  end subroutine calcul_chemin

  subroutine choix_chemin(proba, indice)
    integer :: indice, i
    real(double), dimension(nparapath+1) :: proba
    real(double) :: rand, somme

    call random_number(rand)
    somme = 0.0
    i = 0

    DO WHILE (somme .lt. rand)
       i = i+1
       somme = somme + proba(i)
       indice = i  
    END DO

  end subroutine choix_chemin

  subroutine init_vitesse(config, param) !juste pour les "N+1eme" atomes
    implicit none
    type(atom_config_mc)::config
    integer,intent(in)::param
    real(double) :: v0, v1, z1, z2, z3, z4
    integer :: i, ilast ! if param =1 alors v0 défini avec tinit (juste pour l'initialisation), sinon v0 défini avec Text  

    if ((param==0).and.(tinit.gt.0)) then
       v0 = sqrt(2.D0*bk*tinit)
    else 
       v0 = sqrt(2.D0*bk*Text)
    end if

    ilast=config%im+1-nbatplus 
    do i=config%im,ilast,-1
       !    i = config%im

       call random_number(z1)
       call random_number(z2)
       call random_number(z3)
       call random_number(z4)
       if(z1.eq.0.d0) z1=0.000000001d0
       if(z2.eq.0.d0) z2=0.000000001d0
       if(z3.eq.0.d0) z3=0.000000001d0
       if(z4.eq.0.d0) z4=0.000000001d0

       v1 = one/sqrt(cm(config%ityp(i)))
       config%vp(1,i) = v1*v0*sqrt((-log(z1)))*cos(2.0*pi*z3)
       config%vp(2,i) = v1*v0*sqrt((-log(z1)))*sin(2.0*pi*z3)
       config%vp(3,i) = v1*v0*sqrt((-log(z2)))*cos(2.0*pi*z4)
    end do
  end subroutine init_vitesse



  subroutine boucle_copy_atom(config_n,config_nplus1, sens)
    implicit none

    type(atom_config_mc)::config_n, config_nplus1
    integer :: i
    logical :: sens

    if (sens) then

       DO i=1, config_n%im
          call config_nplus1%copy_atom(i,config_n,i,lextend=.false.)
       END DO

    else
       DO i=1, config_n%im
          call config_n%copy_atom(i,config_nplus1,i,lextend=.false.)
       END DO

    end if

  end subroutine boucle_copy_atom

  subroutine indice_alea(config, ind,itc,nbp)

    implicit none

    type(atom_config_mc),intent(in)::config
    integer,intent(in)::itc,nbp
    integer :: ind(:)
    logical,allocatable::lchosen(:)
    real(double) :: rand
    integer::natyp,iatyp,i,indT
    integer,allocatable::indatyp(:)

    natyp=0
    do i=1,config%im
       if (config%ityp(i)==itc) then
          natyp=natyp+1
       end if
    end do
    allocate(lchosen(natyp))
    allocate(indatyp(natyp))
    lchosen(:)=.false.
    iatyp=0
    do i=1,config%im
       if (config%ityp(i)==itc) then
          iatyp=iatyp+1
          indatyp(iatyp)=i
       end if
    end do
    if (iatyp.ne.natyp) then
       write(6,*)'WTF iatyp natyp',iatyp,natyp
       call arret_ndm
    end if

    do i=1,nbp
3      continue
       call random_number(rand)
       indT = ( (natyp - 1) * rand ) + 1
       !       write(6,*)'atom_del', rang,rand,indT
       if (lchosen(indT).eqv..true.) goto 3
       lchosen(indT)=.true.
       ind(i)=indatyp(indT)
    end do
    !  write(*,*) 'indice atome supprimé', ind

  end subroutine indice_alea


  subroutine atom_biais(config, ind)
    implicit none
    type(atom_config_mc)::config
    integer :: ind(:), i
    real(double) :: rand, somme
    logical,allocatable::lchosen(:)
    integer::natyp,iatyp
    integer,allocatable::indatyp(:)
    natyp=0
    do i=1,config%im
       if (config%ityp(i)==itypcalc) then
          natyp=natyp+1
       end if
    end do
    allocate(lchosen(natyp))
    allocate(indatyp(natyp))
    lchosen(:)=.false.
    iatyp=0
    do i=1,config%im
       if (config%ityp(i)==itypcalc) then
          iatyp=iatyp+1
          indatyp(iatyp)=i
       end if
    end do
    if (iatyp.ne.natyp) then
       write(6,*)'WTF iatyp natyp',iatyp,natyp
       call arret_ndm
    end if



    loopatplus:    do i=1,nbatplus
33     continue
       somme = 0.0
       call random_number(rand)

       DO iatyp=1,natyp
          somme = somme + config%proba_des(indatyp(iatyp))

          if (somme.gt.rand) then
             !             write(6,*)'CH',rand,somme,iatyp
             if (lchosen(iatyp).eqv..true.) then
                goto 33
             else
                lchosen(iatyp)=.true.
                ind(i)=indatyp(iatyp)

                cycle loopatplus
             end if
          end if ! if sur les sommes
       END DO !boucle sur les atomes
    end do loopatplus


  end subroutine atom_biais


  subroutine analyse_montecarlo(atdml,celndm,box,name_file)

    implicit none

    type(atom_config_mc)::atdml
    type(cell_config)::celndm
    class(box_config)::box

    character(len=*) :: name_file
    logical :: lperiod

    lperiod = .true.
    if (lbigmaster) then
       if (itetemp>0) then
          if (mod(iteration,itetemp)==0) then
             call caltabtC(celndm,atdml,lperiod,box)
             call calctemp (temp,kine,atdml,celndm,latcomp=.true.)
          end if
       end if
       if (iterasmol>0) then
          if (mod(iteration,iterasmol)==0) then
             call caltabtC(celndm,atdml,lperiod,box)
             call calctemp (temp,kine,atdml,celndm,latcomp=.true.)
             call rasmolT(atdml,box,namefr=name_file,latcomp=.true.,lappend=.true.)
          end if
       end if

    end if
  end subroutine analyse_montecarlo

  subroutine calcul_U(atdml,celndm,box,potist,U_ini)

    use tempinstT_mod,only:tempinstT
    implicit none

    type(atom_config_mc)::atdml
    type(cell_config)::celndm
    class(box_config)::box
    real(double) :: U_ini, potist,Ti

    if (lbigmaster) then
       Ti=tempinstT(atdml,kine,latcomp=.true.)

       U_ini = potist + kine
    end if
  end subroutine calcul_U

  subroutine restart_chemin(travail, dir)

    implicit none
    integer :: dir
    real(double) :: travail

    open(86, FILE = 'restart_file_old', STATUS = 'old')
    DO 
       read(86,'(I3,G25.16E3)', END=10) dir, travail
    END DO
10  close(86)

  end subroutine restart_chemin


  subroutine atom_supp(vecteur,pins,ipp)

    implicit none
    integer,intent(in)::ipp
    real(double), dimension(:,:) :: vecteur
    real(double),intent(out)::pins
    real(double) :: x_nplus1, y_nplus1, z_nplus1,dist,rd !position initiale aleatoire de la N+1eme particule
    integer::i,j,itry,iat,jat
    real(double)::vec(3,1)
    itry=0
    dist=0
    select case (ins_typ)

    case(1,11,55)
       call atom_supp_sph(vecteur,pins,rd)
       rcpath(ipp)=rd

    case(3,33)
       call atom_supp_sl(vecteur,pins,rd)
       rcpath(ipp)=rd

    case(44)
       call atom_supp_line(vecteur,pins,rd)
       rcpath(ipp)=rd


    case(0)
       pins=1
       if ((nparapath .gt. 1) .and. (lparapath)) then 
          do iat=1,nbatplus
1            continue
             itry=itry+1
             !if (itry.gt.1) write(6,*)'INSER',rang,itry,dist
             !write(*,*) 'rang', rang
             do i=1,rang+1
                call random_number(x_nplus1)
                call random_number(y_nplus1)
                call random_number(z_nplus1)
             end do
             vec(1,1) = x_nplus1
             vec(2,1) = y_nplus1
             vec(3,1) = z_nplus1
             !write(6,*)'PLUS1',rang,x_nplus1,y_nplus1,z_nplus1
             call cryst_to_cart(1,vec,boxmcgc_p%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m
             vecteur(:,iat) = vec(:,1)
             do i=1,atconf_n%im
                call distat(vecteur(:,iat),atconf_n%xp(:,i),boxmcgc_p,dist)
                if (dist.le.distminat) then
                   goto 1
                end if
             end do
             do jat=1,iat-1
                call distat(vecteur(:,iat),vecteur(:,jat),boxmcgc_p,dist)
                if (dist.le.distminat) then
                   do j=1, nparapath*2 !ajout pour eviter de tirer le meme atome meme si rang differents.
                      call random_number(x_nplus1)
                      call random_number(y_nplus1)
                      call random_number(z_nplus1)
                   end do
                   goto 1
                end if
             end do
          end do
       else ! cas monoproposal  ou multi en serie
          !call random_seed
          do iat=1,nbatplus
11           continue
             call random_number(x_nplus1)
             call random_number(y_nplus1)
             call random_number(z_nplus1)

             vecteur(1,iat) = x_nplus1
             vecteur(2,iat) = y_nplus1
             vecteur(3,iat) = z_nplus1
             call cryst_to_cart(1,vecteur(:,iat),boxmcgc_p%at,1)
             do i=1,atconf_n%im
                call distat(vecteur(:,iat),atconf_n%xp(:,i),boxmcgc_p,dist)
                if (dist.le.distminat) then
                   goto 11
                end if
             end do
             do jat=1,iat-1
                call distat(vecteur(:,iat),vecteur(:,jat),boxmcgc_p,dist)
                if (dist.le.distminat) then
                   do j=1, nparapath*2 !ajout pour eviter de tirer le meme atome meme si rang differents.
                      call random_number(x_nplus1)
                      call random_number(y_nplus1)
                      call random_number(z_nplus1)
                   end do
                   goto 11
                end if
             end do
          end do

       end if !lparapath = true
    end select
    !  write(*,*) 'atome supplementaire', vecteur
  end subroutine atom_supp


  subroutine noise(var,NS)

    implicit none
    integer,intent(in)::NS
    real(double)  :: var(3,NS)

    real(double)  :: u_1, u_2, v1, v2, r
    integer :: i, ic


    DO i=1, NS
       do ic=1,3
          r = 2.0
          do while (r >= 1.0)
             call random_number(u_1)
             call random_number(u_2)
             v1 = 2.*u_1-1.
             v2 = 2.*u_2-1.
             r = v1*v1 + v2*v2
          end do
          var(ic,i) = v1*sqrt(-2.*log(r)/(r))
       end do
    END DO

  end subroutine noise


  subroutine langevin( direc, protocol,qeff,work) !LANGEVIN
    implicit none

    character(len=3), intent(in) :: protocol
    integer :: direc 
    integer :: i,ic, ip
    real(double) :: Ek_n, Ek_n_plus1, Ek_n_1s4, Ek_n_3s4, dQeff,qeff,&
         &dWeff, dWork
    real(double) :: U_0, U_1, U_l_n_m1, U_l_n, H_l_n, H_l_n_m1, H_l_ini,work

    real(double) :: beta
    real(double)  :: Gl(3,atconf_Nplus1%im)
    real(double)::rga, rga_s4,tempN,tempNP1,kineN,kineNP1

    real(double), dimension(ntyp) :: aux  !pour les calculs d'acceleration
    integer::rgcib,rgem,iloc
    logical::lchange,ldistrib

    !initialisation des energies
    if (lbigmaster) then ! Master général

       U_0 = 0.0
       U_1 = 0.0
       H_l_n_m1 = 0.0
       dQEff  = 0.0
       QEff   = 0.0
       dWEff  = 0.0
       WEff   = 0.0
       Work=0
       Ek_n = 0.0
       H_l_ini = 0.0
       H_l_n = 0.0
       U_l_n = 0.0  
       U_l_n_m1 = 0.0  
       dWork = 0.0
       Work = 0.0
       ip = 0

       beta = 1.0/(bk*Text)

       DO i=1, atconf_Nplus1%im
          do ic=1,3
             Ek_n = Ek_n + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)*cm(atconf_Nplus1%ityp(i))
          enddo
       ENDDO

       call lambda(direc, ip, protocol)
       !       lambda_mc = 0
       !write(6,*)'DIRECI',direc, ip,protocol,lambda_mc
       U_l_n = (1.d0-lambda_mc)*potist_n + lambda_mc*potist_nplus1

       H_l_ini = Ek_n + U_l_n 
       H_l_n   = H_l_ini
       !write(15,*)
       !       write(15,'(A15, G25.16E3,A15, G25.16E3,A15, G25.16E3,A15, G25.16E3)') &
       !            &'lambda_mc ' ,lambda_mc,'Ek_n' , Ek_n, 'U_l_n', U_l_n, 'H_l_ini', H_l_ini
       !             write(6,'(A)') 'lambda_mc,             Ek_n_plus1*erg2ev,        &
       ! &        U_l_n*erg2ev,             H_l_n*erg2ev,            WEff*erg2ev,           dWEff*erg2ev,&
       !&         dWORK*erg2eV,        dQeff*erg2ev'

       !             write(6,'(8G25.16E3)') lambda_mc, Ek_n_plus1*erg2ev,&
       !                       & U_l_n*erg2ev, H_l_n*erg2ev, WEff*erg2ev, dWEff*erg2ev,dWORK*erg2eV, dQeff*erg2ev

    end if ! Master général

    aux(:ntyp) = tstep/(cm(:ntyp)*2.d0)


    !########################################################################################################################
    !                               Ajout/Retrait d'une particule N+1: système N vers N+1 - direction = 0
    !########################################################################################################################

    !    if (direc == 0) then
    DO ip = 1, pas_lambda_mc
       !incrémentation de lambda
       call lambda(direc,ip, protocol)
       !if (rang==0)write(6,*)'DIRECR',rang,direc, ip,protocol,lambda_mc
       !lambda_mc = dble(ip)/dble(pas_lambda_mc)

       if (lbigmaster) then ! Master général
          Ek_n = 0.0
          Ek_n_plus1 = 0.0  
          Ek_n_1s4 = 0.0
          Ek_n_3s4 = 0.0
!                    write(6,*) 'lambda_mc ' ,lambda_mc,rang
          ! faire le pas de langevin (velocity verlet) pour determiner les nouvelles forces et positions

          ! step 1 First half-step velocities update, v(t) -> v(t+dt/2)
          timel = timel+tstep
          rga=exp((-gamlg)*tstep/2)
          rga_s4 = exp((-gamlg)*tstep/4)
          call noise(Gl, atconf_nplus1%im)
          DO i=1, atconf_Nplus1%im
             do ic=1,3
                Ek_n = Ek_n + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)*cm(atconf_Nplus1%ityp(i))

                atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i)*rga &
                     &+ Gl(ic,i)*sqrt(cm(atconf_Nplus1%ityp(i))*gamlg*rga_s4*tstep/beta)&
                     &/cm(atconf_Nplus1%ityp(i))
                Ek_n_1s4 = Ek_n_1s4 + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)*&
                     &cm(atconf_Nplus1%ityp(i))
                atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i) + aux(atconf_Nplus1%ityp(i))&
                     &*atconf_Nplus1%fp(ic,i)
             end do
          END DO
          ! step 2  Coordinate update, x(t)-> x(t+dt)
          DO i=1, atconf_Nplus1%im
             atconf_Nplus1%xp(1:3,i) = atconf_Nplus1%xp(1:3,i) + tstep*atconf_Nplus1%vp(1:3,i)
          END DO

          !recopier les nouvelles positions dans le syst N
          DO i=1, atconf_N%im
             atconf_N%xp(1:3,i) = atconf_Nplus1%xp(1:3,i)
             atconf_N%vp(1:3,i) = atconf_Nplus1%vp(1:3,i)
          END DO
          !conditions periodiques 
          if (lperiod)    then
             call periodbox(boxmcgc_p,atconf_N)
             call periodbox(boxmcgc_p,atconf_Nplus1)
          end if
       end if !on sort du master général
       !En ce point on doit transférer le système N+1 du master 0 vers le master 1

#ifdef PARA
       if (lmaster) then ! on est dans l'un des 2 masters
          rgcib=1;rgem=0
          if(paramcgc%image==0) then !on est dans le master général
             call atconf_Nplus1%send2proc(rgcib,paramcgc%mpi_master,'x')
          else !on est dans le master de N+1
             call atconf_Nplus1%recv(rgem,paramcgc%mpi_master,'x')
          end if
       end if
#endif        
       iloc=0;ldistrib=.false.;lchange=.true.
!       write(6,*)'MCCDBG3 ', rang
       call calfoMCGC(iloc,lchange,ldistrib)
!       write(6,*)'MCCDBG4 ', rang
       if (lbigmaster) then
!          write(6,*)'MCCDBG401 ', rang
          call sigkinetotMC(atconf_n,atconf_nplus1,boxmcgc_p,lambda_mc,sig,sigkine,sigtot)
!          write(6,*)'MCCDBG402 ', rang
          !mise a jour de U_l_n = (1-lambda_mc)*U_0 + lambda_mc*U_1
          U_l_n = (1.0-lambda_mc)*potist_n + lambda_mc*potist_nplus1
          !write(*,*) potist_n, potist_nplus1
          !affichage temperature
          iteration = iteration +1
!          write(6,*)'MCCDBG41 ', rang
          call noise(Gl,atconf_nplus1%im)
          DO i=1, atconf_Nplus1%im
             do ic=1,3
                atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i) + aux(atconf_Nplus1%ityp(i))&
                     &*atconf_Nplus1%fp(ic,i)
                Ek_n_3s4 = Ek_n_3s4 + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)&
                     &*cm(atconf_Nplus1%ityp(i))
                atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i)*rga &
                     & + Gl(ic,i)*sqrt(cm(atconf_Nplus1%ityp(i))*gamlg*rga_s4*tstep/beta)&
                     &/cm(atconf_Nplus1%ityp(i))
                Ek_n_plus1 = Ek_n_plus1 + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)&
                     &*cm(atconf_Nplus1%ityp(i))
             end do
          END DO
          DO i=1, atconf_N%im
             atconf_N%vp(1:3,i) = atconf_Nplus1%vp(1:3,i) 
          END DO
!          write(6,*)'MCCDBG42 ', rang
          call caltabtC(cells_nplus1,atconf_nplus1,lperiod,boxmcgc_p)
          call caltabtC(cells_n,atconf_n,lperiod,boxmcgc_p)
!          write(6,*)'MCCDBG43 ', rang
          call calctemp (tempN,kineN,atconf_N,cells_n,latcomp=.true.)
          call calctemp (tempNP1,kineNP1,atconf_Nplus1,cells_nplus1,latcomp=.true.)
       end if !master general

       if (lbigmaster) then !master general
!          write(6,*)'MCCDBG44 ', rang
          !calcul des energies et travail et chaleur efficaces
          U_l_n_m1 = U_l_n
          H_l_n_m1 = H_l_n
          H_l_n    = Ek_n_plus1  + U_l_n
          dWork = H_l_n - H_l_n_m1
          !          write(6,*)'compHLN',Ek_n_plus1*erg2ev ,U_l_n*erg2ev
          dQEff  = (Ek_n_1s4-Ek_n) + (Ek_n_plus1-Ek_n_3s4)
          !          write(6,*)'compqeff', ((Ek_n_1s4-Ek_n) + (Ek_n_plus1-Ek_n_3s4))*erg2ev
          QEff   = QEff + dQEff
          dWEff  = H_l_n - H_l_n_m1 - dQEff
          WEff   = WEff + dWEff
          work=work+dwork
       !   if (rang==0)  write(6,'(A,3G15.7)')'WEFF dW dH dQ',dweff*erg2ev,dqeff*erg2ev,dwork*erg2ev
          if (protocol == 'MCP') then
             !write(15,'(6A15)') '#lambda_mc ', 'Ek_n_plus1', 'U_l_n',&
             !         &'H_l_n', 'WEff',  'dWEff'
             !             write(6,'(9G25.16E3)') lambda_mc, Ek_n_plus1*erg2ev,&
             !                       & U_l_n*erg2ev, H_l_n*erg2ev, WEff*erg2ev, work*erg2ev, dWEff*erg2ev,dWORK*erg2eV, dQeff*erg2ev
             !write(*,*) 'lambda_mc ' ,lambda_mc, 'Ek_n_plus1', Ek_n_plus1, 'U_l_n',&
             !            & U_l_n, 'H_l_n', H_l_n, 'Work', Work, 'dWork', dWork
          end if
       end if
!       write(6,*)'MCCDBG5 ', rang
    END DO



  end subroutine langevin



  subroutine lambda(dir, nstep, protocol_name)
    implicit none
    integer :: nstep, dir
    character(len=3), intent(in) :: protocol_name


    !alpha = 1.0 !1.5 !doit etre superieur a 1 pou avoir insertion lente au debut et rapide vers lambda =1

    !dans le cas d'un ajout ou d'un retrait, lambda varie de 0 à 1 ou l'inverse
    if (protocol_name == "MCP") then
       lambda_mc = abs(dble(dir)-(dble(nstep)/dble(pas_lambda_mc)))
       !lambda_mc = (abs(dble(dir)-(dble(nstep)/dble(pas_lambda_mc))))**alpha
       !dans le cas d'un équilibrage, on laisse lambda constant à 0
    elseif (protocol_name == "eql") then
       lambda_mc = 0.0
    end if

  end subroutine lambda



  subroutine init_mpi_MCGC

#ifdef PARA
    if (rang==0)write(6,*)'INITMPIMCGC'
    if (lparapath) then 
       if (mod(nprocs,2*nparapath).ne.0) then
          write(6,*)'nprocs/2*nparapath <>0 STOP'
          call MPI_FINALIZE(ierr)
          call arret_ndm
       end if
       parapath%mpi_orig%nproc=nprocs
       parapath%mpi_orig%rank=rang
       parapath%nimage=nparapath
       call MPI_COMM_DUP(MPI_COMM_WORLD,parapath%mpi_orig%comm,ierr)
       call MPI_COMM_GROUP(parapath%mpi_orig%comm,parapath%mpi_orig%group,ierr)
       call commconstr(parapath)

    else
       call initparapuresp(parapath,rang,mpi_WORLD)
    end if


    paramcgc%mpi_orig%nproc= parapath%mpi_image%nproc ! =parapath%mpi_orig%nproc/nparapath
    !paramcgc%mpi_orig%comm= parapath%mpi_image%comm ! =parapath%mpi_orig%nproc/nparapath 
    paramcgc%mpi_orig%rank=parapath%mpi_image%rank
    paramcgc%nimage=2


    call MPI_COMM_DUP(parapath%mpi_image%comm,paramcgc%mpi_orig%comm,ierr)
    call MPI_COMM_GROUP(paramcgc%mpi_orig%comm,paramcgc%mpi_orig%group,ierr)
    call commconstr(paramcgc)

    myidsp=paramcgc%mpi_image%rank
    call MPI_COMM_free(mpi_comm_space,ierr)
    MPI_COMM_space=paramcgc%mpi_image%comm
    call comm_space%init(MPI_COMM_SPACE)
    nprocspace=paramcgc%mpi_image%nproc
    lbigmaster=parapath%lmaster
    lmaster=paramcgc%lmaster
    lmegamaster=.false.
    if (parapath%mpi_orig%rank==0) lmegamaster=.true.

#else

    parapath%mpi_orig%nproc=1
    parapath%mpi_orig%rank=0
    parapath%mpi_image%nproc=1
    parapath%lmaster=.true.



    paramcgc%mpi_orig%nproc=1
    paramcgc%mpi_orig%rank=0
    paramcgc%mpi_image%nproc=1
    myidsp=0
    paramcgc%lmaster=.true.
    nprocspace=1
    lbigmaster=.true.
    lmaster=.true.
    lmegamaster=.true.

#endif

  end subroutine init_mpi_MCGC


  subroutine initNP1(ipp) !PARAPATH DEFINIR LES POINTEURS atconf_nplus1 et atconf_n

    integer,intent(in)::ipp
    real(double), allocatable, dimension(:,:) :: cart_vec_nplus1

    integer::i,iplus
    character :: extension*4
    logical ::lc2d

#ifdef LAMMPS_VERSION
     character*80::namef
#endif
 
    logical::lwrite
    real(double)::pins
    !definir le systeme a N+1 en tirant une position aleatoire pour le N+1eme atome

    allocate (cart_vec_nplus1(3,nbatplus))
    if (ins_typ==2) then
       if (lbigmaster) then
          call atconf_nplus1%init(atconf_n%im,atconf_n%imm,atconf_n%ltabvois,&
               &im_glob=atconf_n%im_glob,imm_glob=imm_glob)
          call type_switch(idirectionmcgc)
          call cells_nplus1%init(boxmcgc_P,cells_n%nox,cells_n%noy,cells_n%noz, cells_n%natperc)
          call cells_n%copy(cells_nplus1,boxmcgc_p)
       else
          call atconf_nplus1%init(atconf_n%im,atconf_n%imm,atconf_n%ltabvois,&
               &im_glob=atconf_n%im_glob,imm_glob=imm_glob)
          call cells_nplus1%init(boxmcgc_p,cells_n%nox,cells_n%noy,cells_n%noz, cells_n%natperc)
          call cells_n%copy(cells_nplus1,boxmcgc_p)
       end if

    else
       if (lbigmaster) then
          call atom_supp(cart_vec_nplus1,pins,ipp)  
          call atconf_nplus1%init(atconf_n%im+nbatplus,atconf_n%imm,atconf_n%ltabvois,&
               &im_glob=atconf_n%im_glob+nbatplus,imm_glob=imm_glob)
          atconf_nplus1%ltabvois=atconf_n%ltabvois
          !call atconf_n%copy_config(atconf_nplus1,lrescl=.false.)
          call boucle_copy_atom(atconf_n,atconf_nplus1, sens= .false.)
          !addition de la n+1eme particule
          do i=1,nbatplus
             iplus=i+atconf_n%im
             atconf_nplus1%xp(1:3,iplus) = cart_vec_nplus1(1:3,i)
             atconf_nplus1%fp(1:3,iplus) = 0
             atconf_nplus1%ityp(iplus) = itypcalc
             atconf_nplus1%num_at_glob(iplus) = iplus
             atconf_Nplus1%proba_ins = pins
             !          write(6,*)'iex pins',i,pins
             !copie de cell puis caltabtC pour redecouper avec la n+1eme particule
          end do
          call init_vitesse(atconf_nplus1,param = 0)
          call cells_nplus1%init(boxmcgc_P,cells_n%nox,cells_n%noy,cells_n%noz, cells_n%natperc)
          call cells_n%copy(cells_nplus1,boxmcgc_p)
       else
          call atconf_nplus1%init(atconf_n%im+nbatplus,atconf_n%imm,atconf_n%ltabvois,&
               &im_glob=atconf_n%im_glob+nbatplus,imm_glob=imm_glob)
          call cells_nplus1%init(boxmcgc_p,cells_n%nox,cells_n%noy,cells_n%noz, cells_n%natperc)
          call cells_n%copy(cells_nplus1,boxmcgc_p)
       end if

    end if
#ifdef PARA

    rgcib=1;rgem=0
    if (lmaster) then 
       if(paramcgc%image==0) then !procs N
          call  atconf_nplus1%send2proc(rgcib,paramcgc%mpi_master)
       else !procs N+1
          call  atconf_nplus1%recv(rgem,paramcgc%mpi_master)
       end if
    end if

    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
       if(paramcgc%image==0) then !procs N
          call init_voisinage(cells_n,pscgc)

       else !procs N+1
          call atconf_nplus1%send2all(0,paramcgc%mpi_image)
          call init_voisinage(cells_nplus1,pscgc)
          call caltabtC(cells_nplus1,atconf_nplus1,lperiod,boxmcgc_p)
       end if

    end if
#endif

    lc2d=.false.
    if((ipotentiel==-10).or.(ipotentiel==-11)) then
       if (lparapath) then
          if(parapath%image+1==ipp) lc2d=.true.
       else
          if (ipp==1) lc2d=.true.
       end if
       if (lc2d) then
          call atconf_nplus1%send2all(0,paramcgc%mpi_orig)
          if (lbigmaster) then
             !           write(6,*)'write configuration N+1  to confNP1.XXX.lmp'
             lwrite=.true.
          else
             lwrite=.false.
          end if
          write(extension,'(i4.4)') ipp

          call config2data (atconf_nplus1%imm,atconf_nplus1%im,&
               atconf_nplus1%xp,atconf_nplus1%ityp,boxmcgc_p%at,ntyp,lwrite,filename='confNP1.'//extension//'.lmp') ! PARAPATH CHANGER LE NOM AVEC INDICE DE LA BOITE


#ifdef PARA
#ifdef LAMMPS_VERSION
          if (lparapath) then
             !          if(parapath%image+1=ipp) then ! assuré par lc2d
             if(paramcgc%image==0) then !procs N
                firsttime_lammps=.true.
                write(extension,'(i4.4)') ipp
                namef='in.lammps.'//extension//'.N'
                write(6,*)'callinit_lammps ',rang,namef
                call init_lammps(namef)
                !              call init_lammps('in.lammps.N')
             else !procs N+1
                firsttime_lammps=.true.
                write(extension,'(i4.4)') ipp
                namef='in.lammps.'//extension//'.NP1'
                write(6,*)'callinit_lammps ',rang,namef
                call init_lammps(namef)
             end if
             !          end if
          else
             !           write(6,*)'COUCOU',rang
             if(paramcgc%image==0) then !procs N
                firsttime_lammps=.true.

                call init_lammps(iopt=1)
             else !procs N+1
                firsttime_lammps=.true.
                call init_lammps(iopt=2)
             end if
          end if


#else
          write(6,*)'Ipotentiel<0 (lammps) et NON LAMMPS_VERSION : stop'
          call MPI_FINALIZE(ierr)
          call arret_ndm
#endif

#else
          write(6,*)'Ipotentiel<0 (lammps) et NON para en MCGC : stop'
          call arret_ndm
#endif       

       end if
    end if


  end subroutine initNP1

  !*********************************************************

  subroutine initN(ipp) !PARAPATH DEFINIR LES POINTEURS atconf_nplus1 et atconf_n

    integer,intent(in)::ipp

    real(double)::xp_np1(3)
    integer::i,i1,i2,j,iplus
    integer,allocatable::indice(:)
    character :: extension*4
    logical ::lc2d
#ifdef LAMMPS_VERSION
    character*80::namef
#endif
    logical::lwrite
    real(double)::poscenter(3,1),postest(3)

    !definir le systeme a N+1 en tirant une position aleatoire pour le N+1eme atome
    allocate(indice(nbatplus))
    call atconf_n%init(atconf_nplus1%im-nbatplus,atconf_nplus1%imm,atconf_nplus1%ltabvois,&
         &im_glob=atconf_nplus1%im_glob-nbatplus,imm_glob=imm_glob)
    atconf_n%ltabvois=atconf_nplus1%ltabvois

    xp_np1(:)=atconf_nplus1%xp(:,atconf_nplus1%im)
    atconf_nplus1%proba_ins= calcul_proba_ins (xp_np1)


    if (lbigmaster) then
       !call atconf_n%copy_config(atconf_nplus1,lrescl=.false.)
       if (.not.lbiais(1)) then
          !SANS BIAIS
          call indice_alea(atconf_Nplus1,indice,itypcalc,nbatplus)
       else
          !AVEC BIAIS
          call atom_biais(atconf_Nplus1,indice)
       end if
       do i=1,nbatplus
          if (indice(i).gt.atconf_N%im) cycle
          do j=atconf_Nplus1%im,atconf_N%im+1,-1
             if (any(indice(1:nbatplus)==j)) then
                cycle
             else
                i1=indice(i)
                i2=j
                indice(i)=i2
                call atconf_Nplus1%switch_atom(i1,i2)
                exit
             end if
          end do
       end do
       !       call atconf_Nplus1%switch_atom(indice,atconf_Nplus1%im)
       call calcul_proba_des ! on calcule les proba_des car on va accepter la désintégration (et donc on aura besoin des proba pour initialiser old)
       do i=1,nbatplus
          iplus=atconf_N%im+i
          !          write(*,'(A25, 3G25.16E3,  A10, G15.6E3, A, I4 )') 'coord atome a retirer', &
          !               &atconf_Nplus1%xp(:,iplus),&
          !               &'proba', atconf_Nplus1%proba_des(iplus),' image ',parapath%image+1
          !write(*,*) 'coord atome a retirer',  atconf_Nplus1%xp(:,atconf_Nplus1%im)
          !on copie les N nouveaux premiers atomes du syst N+1 dans le systeme N
          call boucle_copy_atom(atconf_N,atconf_Nplus1, sens = .true.)
       end do
       if (ins_typ==1) then
          poscenter(:,1)=bublcenter(:)
          call cryst_to_cart(1,poscenter,boxmcgc_p%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m
          postest(:)=-1*(poscenter(:,1)-atconf_Nplus1%xp(:,atconf_Nplus1%im))
          rcpath(ipp)=sqrt(postest(1)**2+postest(2)**2+postest(3)**2)*1d8
       end if


    else
       call atconf_n%init(atconf_nplus1%im-nbatplus,atconf_nplus1%imm,atconf_nplus1%ltabvois,&
            &im_glob=atconf_nplus1%im_glob-nbatplus,imm_glob=imm_glob)
    end if

    call cells_n%init(boxmcgc_p,cells_nplus1%nox,cells_nplus1%noy,cells_nplus1%noz, cells_nplus1%natperc)
    call cells_nplus1%copy(cells_n,boxmcgc_p)


#ifdef PARA
    !Je pense que l'inversion des lignes en dessous est inutile car le bigmaster connait déja N
!!$  rgcib=1;rgem=0
!!$  if (lmaster) then 
!!$     if(paramcgc%image==0) then !procs N
!!$        call  atconf_nplus1%send2proc(rgcib,paramcgc%mpi_master)
!!$     else !procs N+1
!!$        call  atconf_nplus1%recv(rgem,paramcgc%mpi_master)
!!$     end if
!!$  end if

    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
       if(paramcgc%image==0) then !procs N
          call atconf_n%send2all(0,paramcgc%mpi_image)
          call init_voisinage(cells_n,pscgc)
          call caltabtC(cells_n,atconf_n,lperiod,boxmcgc_p)
       else !procs N+1
          !deja fait dans init_simple
          !        call atconf_nplus1%send2all(0,paramcgc%mpi_image)
          !        call caltabtC(cells_nplus1,atconf_nplus1,lperiod,boxmcgc)
          call init_voisinage(cells_nplus1,pscgc)
       end if

    end if
#endif

    lc2d=.false.
    if((ipotentiel==-10).or.(ipotentiel==-11)) then
       if (lparapath) then
          if(parapath%image+1==ipp) lc2d=.true.
       else
          if (ipp==1) lc2d=.true.
       end if
       if (lc2d) then
          call atconf_n%send2all(0,paramcgc%mpi_orig)
          if (lbigmaster) then
             !           write(6,*)'write configuration N+1  to confNP1.XXX.lmp'
             lwrite=.true.
          else
             lwrite=.false.
          end if
          write(extension,'(i4.4)') ipp

          call config2data (atconf_n%imm,atconf_n%im,&
               atconf_n%xp,atconf_n%ityp,boxmcgc_p%at,ntyp,lwrite,filename='confN.'//extension//'.lmp') ! PARAPATH CHANGER LE NOM AVEC INDICE DE LA BOITE


#ifdef PARA
#ifdef LAMMPS_VERSION
          if (lparapath) then
             if(paramcgc%image==0) then !procs N
                firsttime_lammps=.true.
                write(extension,'(i4.4)') ipp
                namef='in.lammps.'//extension//'.N'
                write(6,*)'callinit_lammps ',rang,namef
                call init_lammps(namef)
                !              call init_lammps('in.lammps.N')
             else !procs N+1
                firsttime_lammps=.true.
                write(extension,'(i4.4)') ipp
                namef='in.lammps.'//extension//'.NP1'
                write(6,*)'callinit_lammps ',rang,namef
                call init_lammps(namef)
             end if
          else

             if(paramcgc%image==0) then !procs N
                firsttime_lammps=.true.
                call init_lammps(iopt=1)
             else !procs N+1
                firsttime_lammps=.true.
                call init_lammps(iopt=2)
             end if
          end if


#else
          write(6,*)'Ipotentiel<0 (lammps) et NON LAMMPS_VERSION : stop'
          call MPI_FINALIZE(ierr)
          call arret_ndm
#endif

#else
          write(6,*)'Ipotentiel<0 (lammps) et NON para en MCGC : stop'
          call arret_ndm
#endif       

       end if
    end if

  end subroutine initN



  subroutine calfoMCGC(iloc,lchange,ldistrib)
    use  tempinstT_mod
    integer,intent(in)::iloc
    logical, intent(in)::lchange,ldistrib
    integer::rgcib,rgem,i,iplus
    integer,save::ncalls=0
    logical::lcalcvois
    real(double)::forcebias(3),potisbias



    ncalls=ncalls+1
    lcalcvois=.false.



    if (iloc==1) then
       if (atconf_n%ltabvois)then
          lcalcvois=.true.
       else
          lcalcvois=.false.
       end if
    else
       if (lchange) then 
          if ((mod(iteration,itetabvois)==0).and.(atconf_n%ltabvois))then
             lcalcvois=.true.
          else
             lcalcvois=.false.
          end if
       end if
    end if

#ifdef PARA
    if(paramcgc%image==0) then !procs N
       if (iloc==1) call initloc(atconf_n,cells_n,atmcgcloc,cellmcgcloc,boxmcgc_p,paramcgc,rumax,lperiod&
            &,psc=pscgc,ldistrib=ldistrib,lcalcvois=lcalcvois,lboxchange=lprahman) !initloc contient caltabtc sur atloc
       call pointer_caltabt_calfo(sig_n,potist_n,atconf_n,cells_n,boxmcgc_p,atmcgcloc,cellmcgcloc,paramcgc,&
            &lperiod,lupdate=lchange,psc=pscgc,lcalcvois=lcalcvois,lboxchange=lprahman)
    else !procs N+1
       if (iloc==1)call initloc(atconf_nplus1,cells_nplus1,atmcgcloc,cellmcgcloc,boxmcgc_p,paramcgc,rumax,&
            &lperiod,psc=pscgc,ldistrib=ldistrib,lcalcvois=lcalcvois,lboxchange=lprahman) !initloc contient caltabtc sur atloc
       call pointer_caltabt_calfo(sig_nplus1,potist_nplus1,atconf_nplus1,cells_nplus1,boxmcgc_p,atmcgcloc,cellmcgcloc,paramcgc,&
            &lperiod,lupdate=lchange,psc=pscgc,lcalcvois=lcalcvois,lboxchange=lprahman)

    end if

    call MPI_BARRIER(paramcgc%mpi_orig%comm,ierr)
    !en ce point chacun des deux masters a les forces de son paquet datomes
    if (lmaster) then ! on est dans l'un des 2 masters7
       rgcib=0;rgem=1
       if(paramcgc%image==1) then !on est dans le master de N+1
          call atconf_nplus1%send2proc(rgcib,paramcgc%mpi_master,'fp')
          call paramcgc%mpi_master%send(potist_nplus1,rgcib,1000)
          call paramcgc%mpi_master%send(sig_nplus1,rgcib,1001)
       else !on est dans le master de N qui est le master général
          call atconf_nplus1%recv(rgem,paramcgc%mpi_master,'fp')
          call paramcgc%mpi_master%recv(potist_nplus1,rgem,1000)
          call paramcgc%mpi_master%recv(sig_nplus1,rgem,1001)
       end if
    end if
    !en ce point le master général (rang_orig=0) a les forces de N et N+1    

#else
    if (iloc==1) call initloc(atconf_n,cells_n,atmcgcloc,cellmcgcloc,boxmcgc_p,paramcgc,rumax,lperiod&
         &,psc=pscgc,ldistrib=ldistrib,lcalcvois=lcalcvois) !initloc contient caltabtc sur atloc
    call pointer_caltabt_calfo(sig_n,potist_n,atconf_n,cells_n,boxmcgc_p,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,lupdate=lchange,psc=pscgc,lcalcvois=lcalcvois)

    if (iloc==1) call initloc(atconf_nplus1,cells_nplus1,atmcgcloc,cellmcgcloc,boxmcgc_p,paramcgc,rumax,lperiod&
         &,psc=pscgc,ldistrib=ldistrib,lcalcvois=lcalcvois) !initloc contient caltabtc sur atloc
    call pointer_caltabt_calfo(sig_nplus1,potist_nplus1,atconf_nplus1,cells_nplus1,boxmcgc_p,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,lupdate=lchange,psc=pscgc,lcalcvois=lcalcvois)

#endif

    if (lbigmaster) then

       if (ins_typ.ne.2) then
          !Egalisation des forces pour les deux systemes
          
          do i=1,nbatplus
             iplus=atconf_n%im+i
             if (lspring) then
                call calcforcespring(atconf_nplus1%xp(:,iplus),forcebias,potisbias)
             else
                forcebias(:)=0 ; potisbias=0.
             end if
!             write(666,*)'pot',potist_n*erg2ev/atconf_n%im,potisbias*erg2ev
!             write(666,*)'force',norm2(atconf_nplus1%fp(:,iplus))*erg2ev*1d-8,norm2(forcebias)*erg2ev*1d-8
!             write(6,*)'pot',potist_n*erg2ev/atconf_n%im,potisbias*erg2ev
!             write(6,*)'force',norm2(atconf_nplus1%fp(:,iplus))*erg2ev*1d-8,norm2(forcebias)*erg2ev*1d-8
             atconf_nplus1%fp(:,iplus) =(1-lambda_mc)*forcebias(:)+ lambda_mc*atconf_nplus1%fp(:,iplus)
             potist_n=potist_n+potisbias
          end do
       end if
       DO i=1,atconf_n%im
          atconf_nplus1%fp(:,i) = (1-lambda_mc)*atconf_n%fp(:,i) + lambda_mc*atconf_nplus1%fp(:,i)
       END DO
       
       DO i=1,atconf_n%im
          atconf_n%fp(:,i) = atconf_nplus1%fp(:,i)
       END DO
       
       sig(:,:)=(1-lambda_mc)*sig_n(:,:) + lambda_mc*sig_nplus1(:,:)
       
    end if!end master general
    !write(6,*)'outcfmc',rang,ncalls
  end subroutine calfoMCGC

  subroutine init_atom_config_mc(atconf,imin,immin,ltabvois,nvois,rvois,lreallocate,im_glob,imm_glob)
    class(atom_config_mc),intent(inout)::atconf
    !type(atom_config_mc),intent(inout)::atconf
    integer,intent(in):: imin
    logical,optional, intent(in)::ltabvois,lreallocate
    integer, optional::immin,im_glob,imm_glob,nvois
    real(double),optional::rvois
    logical :: lrealloc

    lrealloc=.false.
    if (present(lreallocate))then
       lrealloc=lreallocate
    end if
    !initialisation de la partie atom_config_d

    call atconf%atom_config_d%init(imin,immin,ltabvois,rvois=rvois,lreallocate=lreallocate,im_glob=im_glob,imm_glob=imm_glob) 
    !initialisation de la partie mc ajoutée
    if ((lrealloc).and.(allocated(atconf%proba_des)))then
       deallocate(atconf%proba_des)
    end if
    if (.not.allocated(atconf%proba_des))then
       allocate(atconf%proba_des(atconf%imm))
    end if
    atconf%proba_des=0
  end subroutine init_atom_config_mc


  subroutine copy_config_mc(atsource,atcible,lrescl)
    implicit none
    class(atom_config_mc),intent(in)::atsource
    class(atom_config)::atcible
    !type(atom_config_mc),intent(in)::atsource
    !type(atom_config_mc)::atcible
    logical,intent(in)::lrescl

    call atsource%atom_config_d%copy_config(atcible, lrescl)

    select type(atcible)
    class is (atom_config_mc)
       select type (atsource)
       class is (atom_config_mc)
          atcible%proba_des(1:atsource%imm)=atsource%proba_des(1:atsource%imm)
          atcible%proba_ins=atsource%proba_ins
       end select
    end select
  end subroutine copy_config_mc



  subroutine copy_atom_mc(atsource,i,atcible,j,lextend,caracT)
    implicit none
    class(atom_config_mc), intent(in)::atsource
    !type(atom_config_mc),intent(in)::atsource
    integer,intent(in):: i
    class(atom_config), intent(inout)::atcible
    !type(atom_config_mc), intent(inout)::atcible
    integer,intent(in):: j
    logical , optional, intent(in) :: lextend
    character(len=*),optional,intent(in)::caracT
    logical::let
    let=.false.
    if (present(lextend)) then
       let=lextend
    end if
    call atsource%atom_config_d%copy_atom(i,atcible,j,let,caracT='xfniewdlpvrugasm')
    select type(atcible)
    class is (atom_config_mc)
       select type (atsource)
       class is (atom_config_mc)
          atcible%proba_des(j) = atsource%proba_des(i)
       end select
    end select
  end subroutine copy_atom_mc



  subroutine switch_atom_mc(atsource,ind_switch_1, ind_switch_2)
    implicit none
    class(atom_config_mc)::atsource
    !type(atom_config_mc) :: atsource
    integer :: ind_switch_1, ind_switch_2
    real(double) :: intermediaire

    call atsource%atom_config_d%switch_atom(ind_switch_1, ind_switch_2)
    intermediaire = atsource%proba_des(ind_switch_1)
    atsource%proba_des(ind_switch_1) =  atsource%proba_des(ind_switch_2)
    atsource%proba_des(ind_switch_2) = intermediaire

  end subroutine switch_atom_mc

  subroutine  atom_supp_sph(vec,pins,rd) ! routine à écrire qui tire une position et fixe proba_ins

    real(double),intent(out)::vec(:,:),pins ,rd! at this point vec should always be (3,1)
    real(double)::poscenter(3,1),postest(3),xins(3)
    real(double)::zf,zt,zr,fhi,theta,rex,somP,somPm1,dist
    integer::itry,i,iex
    !choose vecteur
    itry=0
    block
      real(double)::distm2
      integer::icl
      poscenter(:,1)=bublcenter(:)
      call cryst_to_cart(1,poscenter,boxmcgc_p%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m
      call closest_at(poscenter(:,1),atconf_n,cells_nplus1,boxmcgc_p,.true.,dist=distm2,iclose=icl)
      if (rang==0) then
         write(6,*)
         write(6,*)'closest atom',distm2,icl,atconf_n%ityp(icl)
         write(6,*)'bubl cent',bublcenter
         write(6,*)'atclose',atconf_n%xp(1,icl)/boxmcgc_p%zl(1),atconf_n%xp(2,icl)/boxmcgc_p%zl(2),atconf_n%xp(3,icl)/boxmcgc_p%zl(3)
         write(6,*)
      end if
    end block

22  continue
    itry=itry+1
    if (itry.gt.10000) then
       write(6,*)'ITRY 10000'
       call arret_ndm
    end if
    call random_number(zf)
    fhi=2*pi*zf
    call random_number(zt)
    theta=acos(2*zt-1)

    call random_number(zr)
    !    write(6,*)'atom_supp_sph', rang,zf,zt,zr
    somP=0.
    somPm1=0
    loopi:do i=1,nrins
       somPm1=somP
       somP=somP+probaR(i)
       !       write(6,*)i,probaR(i),somP
       if (zr.le.somP) then
          iex=i-1
          rex=(float(iex)+(zr-somPm1)/probaR(i))*zlmin/(2*nrins)
          pins=probaR(i)
          exit loopi
       end if
    end do loopi
    !   write(6,*)'IEX',iex,rex
    xins(1)=rex*sin(theta)*cos(fhi)
    xins(2)=rex*sin(theta)*sin(fhi)
    xins(3)=rex*cos(theta)
    poscenter(:,1)=bublcenter(:)
    call cryst_to_cart(1,poscenter,boxmcgc_p%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m
    postest(:)=poscenter(:,1)+xins(:)

      do i=1,atconf_n%im
         call distat(postest(:),atconf_n%xp(:,i),boxmcgc_p,dist)
         if (dist.le.distminat) then
!            if (rang==0) write(6,*)'iex TOO close',dist,distminat
            goto 22
         end if
      end do
    vec(:,1)=postest(:)
!!$    select case (ins-typ)
!!$    case(1)
!!$       pins=1/(1+exp(fdfactmcgc*(rex-R0mcgc)))
!!$    case(11)
!!$       pins=exp(-0.5*beta*k_string*(r**2))
!!$    end select

    !    if (rang==0)then
    !       if (itry.gt.1) write(6,*)'image',parapath%image,'NTRY',itry
    !    end if
    !    write(6,'(A,I3,4G15.7)')'atom_supp_sph vec', rang,vec(:,1),rex*1d8
    rd=rex*1d8
    return
  end subroutine atom_supp_sph


  subroutine  atom_supp_sl(vec,pins,rd) ! r

    real(double),intent(out)::vec(:,:),pins ,rd! at this point vec should always be (3,1)
    real(double)::poscenter(3,1),postest(3),xins(3),zx(3)
    real(double)::zt,rex,somP,somPm1,dist,dex,zr,frac
    integer::itry,i,iex,ic
    !choose vecteur
    itry=0
22  continue
!       write(6,*)'ITRY',itry
    itry=itry+1
    if (itry.gt.1000) then

       call arret_ndm
    end if

    call random_number(zx(1))
    call random_number(zx(2))
    call random_number(zx(3))
    !    write(6,*)'atom_supp_sph', rang,zf,zt,zr
    zt=zx(iZLins)
    
    somP=0.
    somPm1=0
    loopi:do i=1,nrins
       somPm1=somP
       somP=somP+probaR(i)

       if (zt.le.somP) then
          iex=i-1
          dex= -boxmcgc_p%zls2(izlins)+(float(i)/nrins)*boxmcgc_p%zl(izlins)
  !     write(6,*)i,probaR(i),somP,zt
 !         write(6,*)'DEX',dex
          pins=probaR(i)
          exit loopi
       end if
    end do loopi


!    write(6,*)'IEX',iex,dex,zt,r0mcgc
    frac=zlcenter(izlins)+dex/norm2(boxmcgc_p%as(:,izlins))
    poscenter=0
    do i=1,3
       if (i.ne.izlins) then
          poscenter(:,1)=poscenter(:,1)+zx(i)*boxmcgc_p%at(:,i)
       else
          poscenter(:,1)=poscenter(:,1)+frac*boxmcgc_p%at(:,i)
       end if
    end do
!    call cryst_to_cart(1,poscenter,boxmcgc_p%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m
    postest(:)=poscenter(:,1)
!    write(6,*)'postest',poscenter(:,1)
    do i=1,atconf_n%im
       call distat(postest(:),atconf_n%xp(:,i),boxmcgc_p,dist)
!       write(6,*)'DIST',i,dist
       if (dist.le.distminat) then
!          write(6,*)'iex TOO close',dist
          goto 22
       end if
    end do
    vec(:,1)=postest(:)
!
    !    if (rang==0)then
    !       if (itry.gt.1) write(6,*)'image',parapath%image,'NTRY',itry
    !    end if
        write(6,'(A,I3,A,3G15.7,A,G15.7)')'atom_supp_slice', rang,' pos=',vec(:,1),' distance= ',dex*1d8
    rd=abs(dex*1d8)
    return
  end subroutine atom_supp_sl


  subroutine  atom_supp_line(vec,pins,rd) ! r

    real(double),intent(out)::vec(:,:),pins ,rd! at this point vec should always be (3,1)
    real(double)::poscenter(3,1),postest(3),xins(3),zx(3)
    real(double)::zt,rex,somP,somPm1,dist,dex,zr,frac,angle
    integer::itry,i,iex,ic,iloop
    !choose vecteur
    itry=0
222  continue
    itry=itry+1
    if (itry.gt.10000) then
       write(6,*)'ITRY 10000'
       call arret_ndm
    end if

    call random_number(zx(1))
    call random_number(zx(2))
    call random_number(zx(3))
    !    write(6,*)'atom_supp_sph', rang,zf,zt,zr

    !    zx(1)  to get the distance from line zx(2) to get the angle around the line zx(3) to get the position along the line
    angle=2*pi*zx(2)
    somP=0.
    somPm1=0
    loopi:do i=1,nrins
       somPm1=somP
       somP=somP+probaR(i)
!       write(6,*)i,probaR(i),somP,zt
       if (zx(1).le.somP) then
          iex=i-1
          dex= (float(iex)+(zx(1)-somPm1)/probaR(i))*0.5*zlmin/nrins
          exit loopi
       end if
    end do loopi
    iloop=0
    do ic=1,3
       if (ic==izlins)then 
          poscenter(ic,1)=zx(3)*boxmcgc_p%at(ic,ic)
       else
          if (iloop==0) then
             poscenter(ic,1)=zlcenter(ic)*boxmcgc_p%at(ic,ic)+dex*cos(angle)
             iloop=1
          else
             poscenter(ic,1)=zlcenter(ic)*boxmcgc_p%at(ic,ic)+dex*sin(angle)
          end if
       end if
    end do
    
    postest(:)=poscenter(:,1)
!    write(6,*)'postest',poscenter(:,1)
    do i=1,atconf_n%im
       call distat(postest(:),atconf_n%xp(:,i),boxmcgc_p,dist)
       if (dist.le.distminat) then

          !          write(6,*)'iex TOO close'
          goto 222
       end if
    end do
    vec(:,1)=postest(:)
    pins=probaR(iex)
    rd=abs(dex*1d8)
    return
  end subroutine atom_supp_line


  subroutine init_instyp
    real(double):: r,zlm2,somP,dist
    integer::i
    zlmin = distmin(boxmcgc_p%at(:,1),boxmcgc_p%at(:,2))
    zlm2 = distmin(boxmcgc_p%at(:,1),boxmcgc_p%at(:,3))
    zlmin = min(zlmin,zlm2)
    zlm2 = distmin(boxmcgc_p%at(:,2),boxmcgc_p%at(:,3))
    zlmin = min(zlmin,zlm2)
    fdfactmcgc=fdfactmcgc*1d8

    somP=0.
    select case(ins_typ)
    case(1)
       do i=0,nrins
          r=float(i)*zlmin/(2*nrins)
          
          probaR(i)=r*r/(1+exp(fdfactmcgc*(r-R0mcgc)))
          if (rang==0)                 write(136,*)i,r,probaR(i)
          !       write(6,*)i,r,fdfactmcgc,r-R0mcgc,probaR(i),probaR(i)/(r*r)
          somP=somP+probaR(i)
       end do


    case(55)
       do i=0,nrins
          r=float(i)*zlmin/(2*nrins)
          
          probaR(i)=4*pi*r*r*form55(fdfactmcgc,r,R0mcgc)! form55=1/(1+exp(fdfactmcgc*(r-R0mcgc)))
          if (rang==0)                 write(136,*)i,r,probaR(i)
          !       write(6,*)i,r,fdfactmcgc,r-R0mcgc,probaR(i),probaR(i)/(r*r)
          somP=somP+probaR(i)

       end do
       call probaUP(Fecalcprob,Frad,form55,fdfactmcgc,R0mcgc)

    case(11)

   !    somP=0.
       do i=0,nrins
          r=float(i)*zlmin/(2*nrins)
          
          probaR(i)=r*r*exp(-0.5*beta*k_spring*(r**2))
          if (rang==0) write(136,*)i,r,probaR(i)
          somP=somP+probaR(i)    
       end do
  !     probaR(:)=probaR(:)/somP
    
    case(3)
       do i=1,nrins
          dist=abs(-boxmcgc_p%zls2(izlins)+(float(i)/nrins)*boxmcgc_p%zl(izlins))
          probaR(i)=dist**2/(1+exp(fdfactmcgc*(dist-R0mcgc)))
          somP=somP+probaR(i)    
          if (rang==0)                 write(136,*)i,dist,probaR(i)

       end do
    case(33)
       do i=1,nrins
          dist=abs(-boxmcgc_p%zls2(izlins)+(float(i)/nrins)*boxmcgc_p%zl(izlins))
          probaR(i)=exp(-0.5*beta*k_spring*(dist**2))
!          if (rang==0)write(6,*)i,dist,-0.5*beta*k_spring*(dist**2),probaR(i)
          somP=somP+probaR(i)    
          if (rang==0)                 write(136,*)i,dist,probaR(i)
       end do
      
    case(44)
       do i=1,nrins
          dist=i*0.5*zlmin/nrins
          probaR(i)=dist*exp(-0.5*beta*k_spring*(dist**2))
          somP=somP+probaR(i)   
          if (rang==0)                 write(136,*)i,dist,probaR(i)
 
       end do

    end select


    probaR(:)=probaR(:)/somP    
    
       
    !    call arret_ndm    
  end subroutine init_instyp

  subroutine calcukcell
    integer::i
    Kcell = 0.5d0*boxmcgc_p%wbox*Sum( boxmcgc_p%hDot(1:3,1:3)**2 )
    Tempcell=Kcell*2./(sum(ihbox0)*bk)
    epsi=0.5d0*MatMul( MatMul( invtrh0, boxmcgc_p%Gmat ), invh0 )
    DO i=1, 3
       epsi(i,i) = epsi(i,i) - 1.d0
    END DO

    grsig = boxmcgc_p%volu * MatMul(boxmcgc_p%invh, MatMul( sigext, boxmcgc_p%invtrh) )
    tension = invVolu0*MatMul( MatMul( h0, grsig), trh0 )
    ! Énergie potentielle de la cellule (Eq. 2.25, Ref.2)
    Ucell = volu0*Sum( tension(1:3,1:3) * epsi(1:3,1:3) )
  end subroutine calcukcell

  subroutine langevinLPR( direc, protocol) !LANGEVIN
    !    use Parrinello_Rahman,only:h,hdo
    implicit none


    character(len=3), intent(in) :: protocol
    integer :: direc 
    integer :: i,ic, ip
    real(double) :: Ek_n, Ek_n_plus1, Ek_n_1s4, Ek_n_3s4, dQeff, Qeff,&
         &dWeff, dWork,tempN,tempNP1,kineN,kineNP1,tempx
    real(double) :: U_0, U_1, U_l_n_m1, U_l_n, H_l_n, H_l_n_m1, H_l_ini
    real(double)::EkP_n_1s4 ,EkP_n_3s4 ,EkP_n_p1 ,EkP_n 


    real(double) :: beta,u1,u2
    real(double),dimension(3,3)::glanh
    real(double)  :: Gl(3,atconf_Nplus1%im)
    real(double)::rga, rga_s4,rgah

    real(double), dimension(ntyp) :: aux  !pour les calculs d'acceleration
    integer::rgcib,rgem,iloc,ic2
    logical::lchange,ldistrib

    !initialisation des energies
    if (lbigmaster) then ! Master général

       U_0 = 0.0
       U_1 = 0.0
       H_l_n_m1 = 0.0
       dQEff  = 0.0
       QEff   = 0.0
       dWEff  = 0.0
       WEff   = 0.0
       Ek_n = 0.0
       H_l_ini = 0.0
       H_l_n = 0.0
       U_l_n = 0.0  
       U_l_n_m1 = 0.0  
       dWork = 0.0
       Work = 0.0
       ip = 0

       beta = 1.0/(bk*Text)

       DO i=1, atconf_Nplus1%im
          do ic=1,3
             Ek_n = Ek_n + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)*cm(atconf_Nplus1%ityp(i))
          enddo
       ENDDO
       call calcUKcell

       call lambda(direc, ip, protocol)
       U_l_n = (1.d0-lambda_mc)*potist_n + lambda_mc*potist_nplus1

       H_l_ini = Ek_n + U_l_n +kcell+ucell
       H_l_n   = H_l_ini

       !if (lmegamaster) write(*,'(A15, G25.16E3,A15, G25.16E3,A15, G25.16E3,A15, G25.16E3)') &
       !      &'lambda_mc ' ,lambda_mc,'Ek_n' , Ek_n, 'U_l_n', U_l_n, 'H_l_ini', H_l_ini
    end if ! Master général

    aux(:ntyp) = tstep/(cm(:ntyp)*2.d0)

    !########################################################################################################################
    !                               Ajout/Retrait d'une particule N+1: système N vers N+1 - direction = 0
    !########################################################################################################################

    !    if (direc == 0) then
    DO ip = 1, pas_lambda_mc
       !incrémentation de lambda
       call lambda(direc,ip, protocol)
       !       lambda_mc = 0

       if (lbigmaster) then ! Master général
       call sigkinetotMC(atconf_n,atconf_nplus1,boxmcgc_p,lambda_mc,sig,sigkine,sigtot)
          Ek_n = 0.0
          Ek_n_plus1 = 0.0  
          Ek_n_1s4 = 0.0
          Ek_n_3s4 = 0.0
          ! step 1 First half-step velocities update, v(t) -> v(t+dt/2)
          timel = timel+tstep
          rga=exp((-gamlg)*tstep/2)
          rga_s4 = exp((-gamlg)*tstep/4)
          call noise(Gl,atconf_nplus1%im)
          ! Dérivée des coordonnées réduites des atomes à l'instant t+dt/2
          DO i=1, atconf_Nplus1%im
             do ic=1,3
                Ek_n = Ek_n + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)*cm(atconf_Nplus1%ityp(i))

                atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i)*rga &
                     &+ Gl(ic,i)*sqrt(cm(atconf_Nplus1%ityp(i))*gamlg*rga_s4*tstep/beta)&
                     &/cm(atconf_Nplus1%ityp(i))

                Ek_n_1s4 = Ek_n_1s4 + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)*&
                     &cm(atconf_Nplus1%ityp(i))
                atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i) + aux(atconf_Nplus1%ityp(i))&
                     &*atconf_Nplus1%fp(ic,i)
             end do
          END DO
          tempx= tempinstT(atconf_n)
          ! De même pour les vitesses au cas où, par exemple, on utilise le thermostat
          sdot(:,1: atconf_Nplus1%im) = MatMul(boxmcgc_p%invh(:,:), atconf_Nplus1%vp(:,1: atconf_Nplus1%im) )
          sp(:,1: atconf_Nplus1%im) = MatMul(boxmcgc_p%invh(:,:), atconf_Nplus1%xp(:,1: atconf_Nplus1%im) )


          rgah=exp(-gamlg*gamprfact*tstep/2)
          do ic=1,3
             do ic2=1,3
                call random_number(u1)
                call random_number(u2)
                glanh(ic,ic2)=sqrt(-2.*log(u1))*cos(2.*pi*u2)   
             end do
          end do
          call calcUKcell
          EkP_n=kcell

          boxmcgc_p%hdot(:,:) = (  boxmcgc_p%hdot(:,:)*rgah  &
               + (glanh(:,:)/boxmcgc_p%wbox)*sqrt(boxmcgc_p%wbox*bk*text*(1-rgah))  )*ihbox0(:,:)
          call calcUKcell
          EkP_n_1s4=kcell
          boxmcgc_p%hdot(:,:) =(boxmcgc_p%hdot(:,:) +&
               &tstep/(2.d0*boxmcgc_p%wBox)*boxmcgc_p%volu*MatMul(sigtot(:,:)-sigext(:,:),boxmcgc_p%invtrh(:,:)))*ihbox0(:,:)
          sp(:,1:atconf_Nplus1%im) = sp(:,1:atconf_Nplus1%im) + sdot(:,1:atconf_Nplus1%im)*tstep

          ! Tenseur h à l'instant t+dt
          boxmcgc_p%h(:,:) = boxmcgc_p%h(:,:) + boxmcgc_p%hdot(:,:)*tstep*ihbox0(:,:)
          ! Coordonnées réelles à l'instant t+dt
          atconf_Nplus1%xp(:,1:atconf_Nplus1%im) = MatMul( boxmcgc_p%h, sp(:,1:atconf_Nplus1%im) )
          call updatebox(boxmcgc_p,boxmcgc_p%h)
          atconf_Nplus1%vp(:,1:atconf_Nplus1%im) = MatMul( boxmcgc_p%h(:,:), sdot(:,1:atconf_Nplus1%im) ) ! retour à vp car transfert d'atomes  dans scalebox en PARA
          CALL ScaleBox(atconf_Nplus1,cells_nplus1,boxmcgc_p,pscgc)

          !recopier les nouvelles positions dans le syst N
          DO i=1, atconf_N%im
             atconf_N%xp(1:3,i) = atconf_Nplus1%xp(1:3,i)
             atconf_N%vp(1:3,i) = atconf_Nplus1%vp(1:3,i)
          END DO
          call caltabtC(cells_nplus1,atconf_nplus1,lperiod,boxmcgc_p)
          call calctemp (tempNP1,kineNP1,atconf_Nplus1,cells_nplus1,latcomp=.true.)
          call caltabtC(cells_n,atconf_n,lperiod,boxmcgc_p)
          call calctemp (tempN,kineN,atconf_N,cells_n,latcomp=.true.)
          call calcUKcell
          Tempcell=Kcell*2./(sum(ihbox0)*bk)

       end if !on sort du master général (bigmaster)
       !En ce point on doit transférer le système N+1 du master 0 vers le master 1

#ifdef PARA
       if (lmaster) then ! on est dans l'un des 2 masters
!!$          rgcib=1;rgem=0
!!$          if(paramcgc%image==0) then !on est dans le master général
!!$             call boxmcgc_p%send2proc(rgcib,paramcgc%mpi_master)
!!$          else !on est dans le master de N+1
!!$             call boxmcgc_p%send2proc(rgem,paramcgc%mpi_master)
!!$          end if

          call boxmcgc_p%master2slave(0,paramcgc%mpi_master) 
          rgcib=1;rgem=0
          if(paramcgc%image==0) then !on est dans le master général
             call atconf_Nplus1%send2proc(rgcib,paramcgc%mpi_master,'x')
          else !on est dans le master de N+1
             call atconf_Nplus1%recv(rgem,paramcgc%mpi_master,'x')
          end if
       end if
#endif        
       iloc=0;ldistrib=.false.;lchange=.true.

       call calfoMCGC(iloc,lchange,ldistrib)

       if (lbigmaster) then
       call sigkinetotMC(atconf_n,atconf_nplus1,boxmcgc_p,lambda_mc,sig,sigkine,sigtot)
          tempx= tempinstT(atconf_n)
          Kcell = 0.5d0*boxmcgc_p%wbox*Sum( boxmcgc_p%hDot(1:3,1:3)**2 )
          Tempcell=Kcell*2./(sum(ihbox0)*bk)
          !          write(6,*)'PR35',boxmcgc_p%hdot(1,1),boxmcgc_p%h(1,1)*1d8,sigtot(1,1)*unitP,tempx,tempcell


          !mise a jour de U_l_n = (1-lambda_mc)*U_0 + lambda_mc*U_1
          U_l_n = (1.0-lambda_mc)*potist_n + lambda_mc*potist_nplus1
          !write(*,*) potist_n, potist_nplus1
          !affichage temperature
          iteration = iteration +1
          call noise(Gl,atconf_nplus1%im)
          ! Dérivée des coordonnées réduites des atomes à l'instant t+dt/2
          Ek_n_plus1=0  ;      Ek_n_3s4 =0
          DO i=1, atconf_Nplus1%im
             do ic=1,3
                atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i) + aux(atconf_Nplus1%ityp(i))&
                     &*atconf_Nplus1%fp(ic,i)
                Ek_n_3s4 = Ek_n_3s4 + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)&
                     &*cm(atconf_Nplus1%ityp(i))

                atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i)*rga &
                     &+ Gl(ic,i)*sqrt(cm(atconf_Nplus1%ityp(i))*gamlg*rga_s4*tstep/beta)&
                     &/cm(atconf_Nplus1%ityp(i))

                Ek_n_plus1 = Ek_n_plus1 + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)&
                     &*cm(atconf_Nplus1%ityp(i))
             end do
          END DO

          DO i=1, atconf_N%im
             atconf_N%vp(1:3,i) = atconf_Nplus1%vp(1:3,i) 
          END DO
          ! De même pour les vitesses au cas où, par exemple, on utilise le thermostat
          sdot(:,1: atconf_Nplus1%im) = MatMul(boxmcgc_p%invh(:,:), atconf_Nplus1%vp(:,1: atconf_Nplus1%im) )


          boxmcgc_p%hdot(:,:) =(boxmcgc_p%hdot(:,:) +&
               &tstep/(2.d0*boxmcgc_p%wBox)*boxmcgc_p%volu*MatMul(sigtot(:,:)-sigext(:,:),boxmcgc_p%invtrh(:,:)))*ihbox0(:,:)

          call calcUKcell
          EkP_n_3s4=kcell

          boxmcgc_p%hdot(:,:) = (  boxmcgc_p%hdot(:,:)*rgah  &
               + (glanh(:,:)/boxmcgc_p%wbox)*sqrt(boxmcgc_p%wbox*bk*text*(1-rgah))  )*ihbox0(:,:)


          call calcUKcell
          EkP_n_p1=kcell
          !calcul des energies et travail et chaleur efficaces
          U_l_n_m1 = U_l_n
          H_l_n_m1 = H_l_n
          H_l_n    = Ek_n_plus1  + U_l_n+kcell+ucell
          !          write(6,*)'compHLN',Ek_n_plus1*erg2ev ,U_l_n*erg2ev,kcell*erg2ev,ucell*erg2ev
          dWork = H_l_n - H_l_n_m1

          Work = Work + dWork
          dQEff  = (Ek_n_1s4-Ek_n) + (Ek_n_plus1-Ek_n_3s4) + (EkP_n_1s4-EkP_n) + (EkP_n_p1-EkP_n_3s4)
          !          write(6,*)'compqeff', ((Ek_n_1s4-Ek_n) + (Ek_n_plus1-Ek_n_3s4))*erg2ev, ((EkP_n_1s4-EkP_n) + (EkP_n_p1-EkP_n_3s4))*erg2ev
          QEff   = QEff + dQEff
          dWEff  = H_l_n - H_l_n_m1 - dQEff
          WEff   = WEff + dWEff

       end if

    END DO

    !    endif


!!$
!!$    !########################################################################################################################


  end subroutine langevinLPR

  subroutine type_switch(direction)
    integer,intent(in)::direction
    integer,allocatable::indice(:)
    integer::i,i1,i2,j

    allocate(indice(nbatplus))

    if (direction==0) then 
       call indice_alea (atconf_N, indice,typswitch1,nbatplus)
       do i=1,nbatplus
          if (atconf_N%ityp(indice(i)).ne.typswitch1) then
             write(6,*)'WTF type_switch'
             call arret_ndm
          end if
          atconf_Nplus1%ityp(indice(i))=typswitch2
       end do
    else
       call indice_alea (atconf_Nplus1, indice,typswitch2,nbatplus)
       do i=1,nbatplus
          if (atconf_Nplus1%ityp(indice(i)).ne.typswitch2) then
             write(6,*)'WTF type_switch2'
             call arret_ndm
          end if
          atconf_N%ityp(indice(i))=typswitch1
       end do
    end if
    do i=1,nbatplus
       if (indice(i).gt.atconf_N%im-nbatplus) cycle
       do j=atconf_Nplus1%im,atconf_N%im-nbatplus+1,-1
          if (any(indice(1:nbatplus)==j)) then
             cycle
          else
             i1=indice(i)
             i2=j
             indice(i)=i2
             call atconf_Nplus1%switch_atom(i1,i2)
             exit
          end if
       end do
    end do

  end subroutine type_switch


  subroutine calcdistspring(pos,dist,normout)
    real(double),intent(in)::pos(3)
    real(double),intent(out)::dist,normout(3)
    real(double)::poscenter(3,1),postest(3),posred(3,1)
    integer::ic
    posred(:,1)=pos(:)
    normout(:)=0
    call cryst_to_cart(1,posred,boxmcgc_p%bg,-1) 
!    write(6,*)'POS', pos
!    write(6,*)'POSres', posred
    
    select case (ins_typ)
    case(11,55)
       poscenter(:,1)=bublcenter(:)
       call cryst_to_cart(1,poscenter,boxmcgc_p%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m
       postest(:)=pos(:)-poscenter(:,1)
!       write(6,*)'POSTEST',postest
       dist=sqrt(postest(1)**2+postest(2)**2+postest(3)**2)
       normout(:)=postest(:)/dist
    case(33)
       
       dist=abs(boxmcgc_p%at(izlins,izlins)*(posred(izlins,1)-zlcenter(izlins)))
!       write(6,*)'decd ',boxmcgc_p%at(izlins,izlins),posred(izlins,1),zlcenter(izlins)
       normout(izlins)=sign(1.,posred(izlins,1)-zlcenter(izlins))
    case(44)
       dist=0
       do ic=1,3
          if (ic.ne.izlins) then
             dist=dist+((posred(ic,1)-zlcenter(ic))*boxmcgc_p%at(ic,ic))**2
             normout(ic)=(posred(ic,1)-zlcenter(ic))*boxmcgc_p%at(ic,ic)
          end if
       end do
       dist=dsqrt(dist)
       normout=normout/dist
    end select

  end subroutine calcdistspring
  subroutine calcforcespring(pos,forceb,potisb)
    real(double),intent(in)::pos(3)
    real(double),intent(out):: forceb(3),potisb
    real(double)::dist,normout(3)
    integer::ic
    call calcdistspring(pos,dist,normout)
    if (ins_typ.ne.55) then 
       potisb=0.5*k_spring*dist**2
       forceb(:)=-1*k_spring*dist*normout(:)
!       write(6,'(A,5G17.5)')'DIST',dist,potisb*erg2ev,forceb
    else
       potisb=(1/beta)*log(1+exp(fdfactmcgc*dist))
       forceb(:)=-1*normout(:)*(fdfactmcgc*exp(fdfactmcgc*dist))/((1+exp(fdfactmcgc*dist))*beta)

    end if
  end subroutine calcforcespring
  
  function form55(fd,r,r0) result(pu)
      USE T_kind_param_m, ONLY:  double
    real(double),intent(in)::fd,r,r0
    real(double)::pu
    pu=1/(1+exp(fd*(r-R0)))
  end function form55

  subroutine probaUP(fcp,Fradc,form,fdfact,r0mc)
    interface
       function myfunc(fd,r,r0) result(pu)
         USE T_kind_param_m, ONLY:  double
         real(double),intent(in)::fd,r,r0
         real(double)::pu
       end function myfunc
    end interface
    procedure (myfunc) :: form
    real(double)::fcp,fradc(:),fdfact,r0mc

    integer::i
    real(double)::r,pins_inter,dr,rp05,rm05,Urp05,Urm05
    pins_inter=0.
    dr=zlmin/(2*nrins)
    do i=1,nrins
       r=float(i)*dr
       pins_inter=pins_inter+4*pi*r**2*dr*form(fdfact,r,r0mc)
       rp05=r+dr*0.5;       rm05=r-dr*0.5
       Urp05=-bk*text*log(form(fdfact,rp05,r0mc))
       Urm05=-bk*text*log(form(fdfact,rm05,r0mc))
       Fradc(i)=-(urp05-urm05)/dr
    end do
    fcp=-bk*text*log(pins_inter/boxmcgc_p%volu)
  end subroutine probaUP
  
end module montecarlo_mod
