module montecarlo_mod
  USE gen_com_m,only:  lperiod, tstep, timel, tstep, sig, itetabvois,&
       & iterasmol,itetemp, temp, kine, pi, bk, Text, gamlg,one,pi,text,tinit,&
       &lspaceNDM,rang,it,firsttime_lammps,posa,forca,erg2ev,parallele
  USE atomconfig,only:atom_config,atom_config_d, config2ndm, switch_atom
  USE cellconfig, only:cell_config, cellconfig2ndm, caltabtC
  USE var_pot,only:ntyp,cm,gamlt
  USE calfo_mod,only: calfo
  USE T_kind_param_m, ONLY:  double
  USE cryst_to_cart_mod, ONLY: cryst_to_cart
  USE caltabi_mod,only: caltabi
  USE boxconfig,only:box_config,periodbox
  USE rasmolT_mod,only: rasmolT
  use paraconfig,only:para_config,commconstr,initparapuresp
#ifdef PARA
  use Tpara,only:grp_world,nprocs,myidsp,MPI_COMM_space,nprocspace,ierr,mpi_comm_world,&
       &NDM_MPI_REAL_DOUBLE,para_space_config,status,comm_space,mpi_world
  use mod_para,only:maj_atomes_frt_ftm
  USE init_vois_mod,only: init_voisinage
#else
  use Tpara,only:myidsp,nprocspace,para_space_config
#endif
  use read_val,only:rvois,ltabvois
  use var_pot,only:ipotentiel,rumax
  USE parautils,only:initloc,pointer_caltabt_calfo
  USE calctemp_mod,only:calctemp

#ifdef LAMMPS_VERSION
  use vars_lammps
  use lammps_util_mod,only:init_lammps
#endif  
  use config2data_mod,only:config2data

  implicit none

  type, extends (atom_config_d):: atom_config_mc
  real(double), allocatable :: proba(:) !defini pr chaque atome mais utile que pour O dans notre cas
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


type(box_config)::boxmcgc

integer::  pas_lambda_mc
real(double) :: lambda_mc !lambda compris entre 0 et 1

integer :: n_path ! nb de chemin d'insertion, a definir dans .din, par defaut 10

real(double), dimension(3,3) :: sig_n, sig_nplus1
real(double) :: potist_n, potist_nplus1
real(double) :: Weff, Work
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
integer, dimension(12) :: seed 
integer:: nparapath
logical lparapath

contains 



  subroutine montecarlo


    implicit none

    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------

    type(atom_config_mc):: config_atom_old_0, config_atom_old_1, config_atom_new_0, config_atom_new_1 !config intermediaire pour suivre l'evolution des systemes: 0 -> syst N, 1 -> syst N+1.


    integer :: i,ic
    logical :: lextend

    integer :: direction, direc, dir
    integer :: n_accepted, n_accepted_0, n_accepted_1
    integer :: n_gen, n_gen_0, n_gen_1

    integer :: i_path,ipp,ipch
    real(double)::zr1


    real(double) :: biais
    real(double) :: theta
    real(double) :: beta, beta_eV

    real(double) :: W, Wprec, xprob, xalea
    real(double) :: ln_xalea, ln_Wprec, ln_W, ln_xprob
    real(double) :: acceptance_rate, acceptance_rate_0, acceptance_rate_1

    integer :: n,iloc
    integer :: acceptation, premier_accept

    real(double) :: Wprecedent !sauvegarde Wprec pour posttraitement
    real(double) :: mu_moy, mu_wrmc, mu_NC, mu_DC

    real(double) :: contribut_accepte(0:1), f2(0:1), fminusf(0:1)
    real(double) :: f_cumul(0:1), f2_cumul(0:1), fminusf_cumul(0:1)
    real(double) :: f_wr(0:1), f2_wr(0:1), fminusf_wr(0:1)
    real(double) :: f_wr_cumul(0:1), f2_wr_cumul(0:1), fminusf_wr_cumul(0:1)
    real(double) :: f_cumul_inte(0:1), f2_cumul_inte(0:1), fminusf_cumul_inte(0:1)
    real(double) :: f_wr_cumul_inte(0:1), f2_wr_cumul_inte(0:1), fminusf_wr_cumul_inte(0:1)

    real(double) :: b_opt(0:1), b_wr_opt(0:1), est_opt(0:1), est_opt_bwr(0:1)
    real(double),allocatable:: Weff_npp(:)
    real(double), dimension(nparapath+1) :: xprob_i

    logical :: lchange,ldistrib,lcalc

    !########################################################################################################################
    !                                             Initialisation
    !########################################################################################################################

    !initialisation variables 
    !pour le premier chemin: sens positif, d'ajout d'une particule et acceptation

!!$ do i=1,size(seed)
!!$    seed(i) = 152+rang*10*i*100
!!$ end do
 call random_seed!(PUT=seed(1:12))


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

    theta = 0.5
    beta = 1.0/(bk*Text)
    beta_eV = 1.0/((8.617333262145E-5)*Text)

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
    mu_DC = 0.0
    contribut_accepte(:) = 0.0
    f2(:) = 0.0
    fminusf(:) = 0.0
    f_wr(:) = 0.0
    f2_wr(:) = 0.0
    fminusf_wr(:) = 0.0
    f_cumul(:) = 0.0
    f2_cumul(:) = 0.0
    fminusf_cumul(:) = 0.0
    f_wr_cumul(:) = 0.0
    f2_wr_cumul(:) = 0.0
    fminusf_wr_cumul(:) = 0.0
    b_opt(:) = 0.0
    b_wr_opt(:) = 0.0 
    est_opt(:) = 0.0
    est_opt_bwr(:) = 0.0
    f_cumul_inte(:) = 0.0
    f2_cumul_inte(:) = 0.0
    fminusf_cumul_inte(:) = 0.0
    f_wr_cumul_inte(:) = 0.0
    f2_wr_cumul_inte(:) = 0.0
    fminusf_wr_cumul_inte(:) = 0.0

    xprob_i(:) = 0.0
    premier_accept = 0    

    !deplacé !
    if (lbigmaster) then
       ! sauvegarde du système
       call config_atom_n(1)%copy_config(config_atom_old_0, lrescl=.true.)
    end if


    if (nparapath.gt.0) then
       allocate (Weff_npp(nparapath))
       Weff_npp(:)=0
    end if

    do ipp=1,nparapath
       lcalc=.false.
       if (lparapath) then
          if (parapath%image+1==ipp) lcalc=.true.
       else
          lcalc=.true.
       end if
       if (lcalc) then 
          atconf_n=> config_atom_n(ipp)
          cells_n=>config_cells_n(ipp)
          atconf_nplus1=>config_atom_nplus1(ipp)
          cells_nplus1=>config_cells_nplus1(ipp)


          direction = 0 ! direction = 0 on ajoute un atome, = 1 on retire un atome

          call lambda(direction, nstep = 0, protocol_name = 'MCP') !initialisation du lambda a 0 pour le premier melange des forces
          iloc=1;lchange=.false.;ldistrib=.true.
          call calfoMCGC(iloc,lchange,ldistrib) 

          ! on relaxe le systeme initial 
          !call langevin(direction, protocol = 'eql') ! sinon deplace l'atome N+1

!!$ if (lbigmaster) then
!!$    ! sauvegarde du système
!!$    call atconf_n%copy_config(config_atom_old_0, lrescl=.true.)
!!$    !call analyse_montecarlo(atconf_n,cells_n,boxmcgc, 'syst_UO2n_postinit')
!!$ end if
        
          !pour le premier chemin: sens positif, d'ajout d'une particule et acceptation
          call langevin(direction, protocol = 'MCP')

          weff_npp(ipp)=Weff
          !if (lbigmaster)write(6,*)'potist', ipp,potist_n,potist_nplus1
       end if
    end do

    if ((lbigmaster).and.(lparapath)) then
       call parapath%mpi_master%sum(weff_npp)
    end if

    if (lbigmaster) then 
       if (lmegamaster) then
          if (nparapath.gt.1) then
             call random_number(zr1)

             ipch=1+int(nparapath*zr1) ! choix aléatoire débile
             write(6,*)'chemin choisi',ipch,zr1
          else
             ipch=1
          end if
       end if

       call parapath%mpi_master%bcast(0,ipch)

       Weff=weff_npp(ipch)

#ifdef PARA
       if (lparapath) then 
          call config_atom_n(ipch)%send2all(ipch-1,parapath%mpi_master)
          call config_atom_nplus1(ipch)%send2all(ipch-1,parapath%mpi_master)
!!$#ifdef PARA
!!$ call MPI_FINALIZE(ierr)
!!$#endif
!!$ stop

          call config_cells_n(ipch)%send2all(ipch-1,parapath%mpi_master)
          call config_cells_nplus1(ipch)%send2all(ipch-1,parapath%mpi_master)
       end if
#endif          

       atconf_nplus1=>config_atom_nplus1(ipch)
       call calcul_proba
       call config_atom_nplus1(ipch)%copy_config(config_atom_old_1, lrescl=.true.)      

       W = WEff
       !W = Work
       xprob = 1
       Wprec = + W
       Wprecedent = Wprec
       write(*,*) 'W0', W,W*erg2eV



    end if
    direction = 1

    !########################################################################################################################
    !                                             boucle sur lambda le long d'un chemin
    !########################################################################################################################

    DO i_path = 1, n_path ! boucle à faire pour tous les procs
       Weff_npp(:)=0
       if (lbigmaster) then
          config_atom_nplus1(ipch)%vp(:,:)   = - config_atom_nplus1(ipch)%vp(:,:) !à chaque retour dans la boucle, on change de direction
          config_atom_n(ipch)%vp(:,:)   = - config_atom_n(ipch)%vp(:,:)
          if (direction == 0) then
             call config_atom_n(ipch)%copy_config(config_atom_new_0, lrescl=.true.)
          endif
          if (direction == 1) then
             call config_atom_nplus1(ipch)%copy_config(config_atom_new_1, lrescl=.true.)
          endif
       end if
       do ipp=1,nparapath
          lcalc=.false.
          if (lparapath) then
             if (parapath%image+1==ipp) lcalc=.true.
          else
             lcalc=.true.
          end if

          if (lcalc) then 
             atconf_n=> config_atom_n(ipp)
             cells_n=>config_cells_n(ipp)
             atconf_nplus1=>config_atom_nplus1(ipp)
             cells_nplus1=>config_cells_nplus1(ipp)

             if (lbigmaster) then !!master general

                !!          write(*,*) ' '
                !!          write(*,*) ' '
                !!          write(*,*) ' '
                !!          write(*,*) 'Numéro de chemin', i_path,direction

                call random_number(xalea)
                ln_xalea  = log(xalea)
                
             end if !fin master general
             !choisir l'at a retirer ou ajouter + preparation des syst N et N+1 pour etre prets pour le langevin (cad decoupage cellules + calcul forces + melange des forces - se fait dans cette sous routine)
             call ajout_retrait(direction)
             
             ! pas de langevin
             call langevin(direction, protocol = 'MCP')
             !if (lbigmaster) write(6,*)'potist', ipp,potist_n,potist_nplus1

             Weff_npp(ipp)=weff
          end if
       end do
       if ((lbigmaster).and.(lparapath)) then
          call parapath%mpi_master%sum(weff_npp)
       end if
       if (lbigmaster) then !master general
          if (lmegamaster) then
             if (nparapath.gt.1) then
                !call calcul_chemin(Weff_npp, xprob_i, Wprec, direction, theta)
                !call choix_chemin(xprob_i, ipch)
                call random_number(zr1)
                ipch=1+int(nparapath*zr1) ! choix aléatoire débile
                !write(6,*)'chemin choisi',ipch
             else
                ipch=1
             end if
          end if
          call parapath%mpi_master%bcast(0,ipch)
          Weff=weff_npp(ipch)
#ifdef PARA
       if (lparapath) then 
          call config_atom_n(ipch)%send2all(ipch-1,parapath%mpi_master)
          call config_atom_nplus1(ipch)%send2all(ipch-1,parapath%mpi_master)
          call config_cells_n(ipch)%send2all(ipch-1,parapath%mpi_master)
          call config_cells_nplus1(ipch)%send2all(ipch-1,parapath%mpi_master)
       end if
#endif          
          
          do ipp=1,nparapath
             config_cells_n(ipp)= config_cells_n(ipch)
             config_cells_nplus1(ipp)= config_cells_nplus1(ipch)
             config_atom_n(ipp)=config_atom_n(ipch)
             config_atom_nplus1(ipp)=config_atom_nplus1(ipch)
          end do
      

      !call analyse_montecarlo(config_atom_n(ipch),config_cells_n(ipch),boxmcgc, 'UO2_syst_n_after_lang')
       if (direction == 0) then
          W = +WEff
          !W = +Work
          n_gen_0 = n_gen_0 + 1
          atconf_nplus1=>config_atom_nplus1(ipch)
          call calcul_proba
          call config_atom_nplus1(ipch)%copy_config(config_atom_new_1, lrescl=.true.)      
       else
          W = -WEff
          !W = - Work
          n_gen_1 = n_gen_1 + 1
          call config_atom_n(ipch)%copy_config(config_atom_new_0, lrescl=.true.)      
       endif
       n_gen = n_gen + 1
       ln_Wprec  = (+beta*(direction-theta)*Wprec)
       ln_W      = (+beta*(direction-theta)*W)
       if (direction == 0) then
          biais = 1
       else
          biais = config_atom_nplus1(ipch)%proba(config_atom_nplus1(ipch)%im)/config_atom_old_1%proba(atconf_nplus1%im)
       end if
       !ln_xprob  = - dlog(1 + (dexp(ln_Wprec-ln_W)/biais))
       ln_xprob  = - dlog(1 + dexp(ln_Wprec-ln_W)) 
       xprob     = dexp(ln_xprob)

       if (lmegamaster) write(*,*) 'WeV', W*erg2eV, 'WpreceV', Wprec*erg2eV,'XPROB', xprob, 'XALEA', xalea
       !call analyse_montecarlo(config_atom_n(ipch),config_cells_n(ipch),boxmcgc, 'UO2_syst_n_beforetest')
       if (ln_xprob > ln_xalea) then    
!!!!!!!!!!!!!!!!!!!!!!!!!! ACCEPTATION   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !if (lmegamaster) write(*,*) 'ACCEPTATION, direction=', direction,'W', W*erg2eV, &
          !     &'Wprec', Wprec*erg2eV, '  LN_XPROB ', ln_xprob, '  XPROB ', xprob, '  XALEA ', xalea

          premier_accept = 1

          if (direction == 0) then
             n_accepted_0 = n_accepted_0 + 1
          else
             n_accepted_1 = n_accepted_1 + 1
          endif

!!!!!!!Calcul de mu!!!!!!!!
!!!!! Calcul d'une moyenne classique / estimateur WR/NC !!!!!!!!!
          contribut_accepte(0) = dexp(beta*W*0.5)
          f2(0) = dexp(beta*W)
          fminusf(0) =(dexp(beta*W*0.5) - dexp(beta*Wprecedent*0.5))**2 

          contribut_accepte(1) = dexp(-beta*W*0.5)
          f2(1) = dexp(-beta*W)
          fminusf(1) =(dexp(-beta*W*0.5) - dexp(-beta*Wprecedent*0.5))**2
!!!!!!!!!! fin calcul mu !!!!!!!!!!!

          n_accepted = n_accepted + 1
          Wprec = + W 
          dir = direction
          acceptation = 1

          call config_atom_new_0%copy_config(config_atom_old_0, lrescl=.true.)
          call config_atom_new_1%copy_config(config_atom_old_1, lrescl=.true.)          

       else
!!!!!!!!!!!!!!!!!!!!!!!!!! REFUS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !if (lmegamaster) write(*,*) ' REJECTION, direction=', direction, 'W', W*erg2eV, &
          !     &'Wprec', Wprec*erg2eV, '  LN_XPROB ', ln_xprob, '  XPROB ', xprob, '  XALEA ', xalea

          !on accepte le sens opposé - changer des signes des vitesses 
          config_atom_old_0%vp(:,:)   = - config_atom_old_0%vp(:,:)
          config_atom_old_1%vp(:,:)   = - config_atom_old_1%vp(:,:)

          acceptation = 0
          if (direction == 0) then
             call config_atom_old_1%copy_config(config_atom_nplus1(ipch), lrescl=.true.)
             call caltabtC(config_cells_nplus1(ipch),config_atom_nplus1(ipch),lperiod,boxmcgc)
          end if

          if (direction == 1) then
             call config_atom_old_0%copy_config(config_atom_n(ipch), lrescl=.true.)
             call caltabtC(config_cells_n(ipch),config_atom_n(ipch),lperiod,boxmcgc)
          end if

!!!!!!!Calcul de mu!!!!!!!!
!!!!! Calcul d'une moyenne classique / estimateur WR/NC !!!!!!!!!
          contribut_accepte(0) = dexp(beta*Wprecedent*0.5)
          f2(0) = dexp(beta*Wprecedent)
          fminusf(0) = 0.0

          contribut_accepte(1) = dexp(-beta*Wprecedent*0.5)
          f2(1) = dexp(-beta*Wprecedent)
          fminusf(1) = 0.0
!!!!!!!!!! fin calcul mu !!!!!!!!!!!

       end if !test sur xprob


       if (nparapath.gt.1) then
          do ipp=1,nparapath
             config_cells_n(ipp)= config_cells_n(ipch)
             config_cells_nplus1(ipp)= config_cells_nplus1(ipch)
             config_atom_n(ipp)=config_atom_n(ipch)
             config_atom_nplus1(ipp)=config_atom_nplus1(ipch)
          end do
       end if



       acceptance_rate   = (real(n_accepted)/real(n_gen))*1.0d2
       acceptance_rate_0 = (real(n_accepted_0)/real(n_gen_0))*1.0d2
       acceptance_rate_1 = (real(n_accepted_1)/real(n_gen_1))*1.0d2

!!!!! Calcul de mu !!!!!!
!!!!! Calcul estimateur WR/DC !!!!!!!!!
       f_wr(0) = xprob*dexp(beta*W*0.5) + (1.0-xprob)*dexp(beta*Wprecedent*0.5)
       f2_wr(0) = xprob*dexp(beta*W)&
            & + (1.0-xprob)*dexp(beta*Wprecedent)
       fminusf_wr(0) = (1.0-xprob) * xprob * ( dexp(beta*W*0.5)&
            & - dexp(beta*Wprecedent*0.5) )**2

       f_wr(1) = xprob*dexp(-beta*W*0.5) + (1.0-xprob)*dexp(-beta*Wprecedent*0.5)
       f2_wr(1) = xprob*dexp(-beta*W)&
            & + (1.0-xprob)*dexp(-beta*Wprecedent)
       fminusf_wr(1) = (1.0-xprob) * xprob * ( dexp(-beta*W*0.5)&
            & - dexp(-beta*Wprecedent*0.5) )**2
       !write(*,*) 'contribut_accepte', contribut_accepte(0), contribut_accepte(1)
       !write(*,*) 'f_wr', f_wr(direction)

       !write(*,*) 'f2_wr', f2_wr(direction)
       !write(*,*) 'f2', f2(direction)

       !write(*,*) 'f_wr**2', f_wr(direction)**2
       !write(*,*) 'contribut_accepte**2', contribut_accepte(direction)**2

       !write(*,*) 'fminusf_wr', fminusf_wr(direction)
       !write(*,*) 'fminusf', fminusf(direction)

!!!!! formation des sommes !!!!
       DO direc = 0,1
          f_cumul(direc) = f_cumul(direc) + contribut_accepte(direc)
          f2_cumul(direc) = f2_cumul(direc) + f2(direc)
          fminusf_cumul(direc) = fminusf_cumul(direc) + fminusf(direc)

          f_wr_cumul(direc) = f_wr_cumul(direc) + f_wr(direc)
          f2_wr_cumul(direc) = f2_wr_cumul(direc) + f2_wr(direc)
          fminusf_wr_cumul(direc) = fminusf_wr_cumul(direc) + fminusf_wr(direc)


!!!!! diviser par le nombre de chemin!!!!!
          f_cumul_inte(direc) = f_cumul(direc) / real(n_gen)
          f2_cumul_inte(direc) = f2_cumul(direc) / real(n_gen)
          fminusf_cumul_inte(direc) = fminusf_cumul(direc) / (real(n_gen*2.0))
       
          f_wr_cumul_inte(direc) = f_wr_cumul(direc) / real(n_gen)
          f2_wr_cumul_inte(direc) = f2_wr_cumul(direc) / real(n_gen)
          fminusf_wr_cumul_inte(direc) = fminusf_wr_cumul(direc) / real(n_gen)


!!!!! calcul de la variable de contrôle !!!!
          if (premier_accept == 1) then
            b_opt(direc)=(f2_cumul_inte(direc)-f_cumul_inte(direc)**2)&
               & / fminusf_cumul_inte(direc)
            b_wr_opt(direc)=(f2_wr_cumul_inte(direc)-f_wr_cumul_inte(direc)**2)&
               & / fminusf_wr_cumul_inte(direc)


!!!!! calcul de l'estimateur NC et DC !!!!!!
            est_opt(direc) = b_opt(direc)*f_wr_cumul_inte(direc) &
               & + (1.0-b_opt(direc))*f_cumul_inte(direc)
            est_opt_bwr(direc) = b_wr_opt(direc)*f_wr_cumul_inte(direc) &
               & + (1.0-b_wr_opt(direc))*f_cumul_inte(direc)
          end if
       END DO
       !write(*,*) 'direc b_opt et b_wr_opt', direction, est_opt(0), est_opt_bwr(0)
       !write(*,*) 'direc b_opt et b_wr_opt', direction, est_opt(1), est_opt_bwr(1)

!!!! calcul de mu!!!!!
       if (premier_accept == 1) then
          mu_moy = -(1/beta_eV)*dlog(f_cumul_inte(1) / f_cumul_inte(0)) !estimateur simple Im(f)
          mu_wrmc = -(1/beta_eV)*dlog(f_wr_cumul_inte(1) / f_wr_cumul_inte(0)) !estimateur WRMC (moyenne des Wgen pondérée par les proba)
          mu_NC = -(1/beta_eV)*dlog(est_opt(1) / est_opt(0))  !estimateur NC Jnc,M(f)
          mu_DC = -(1/beta_eV)*dlog(est_opt_bwr(1) / est_opt_bwr(0)) !estimateur DC Jdc,M(f)
       end if
!!!!!!!!!! fin calcul mu !!!!!!!!!!!          

       if (i_path == 1) then
         if (lmegamaster) write(*,'(A)') 'acceptation  num_chemin  direction b_opt(0) b_opt(1) b_wr_opt(0) b_wr_opt(1)&
         & mu_moy   mu_WRMC  mu_NC  mu_DC   WeV   WprecedenteV  proba_acc proba_refus'
       end if

       if (lmegamaster) write(*,'(3I5, 12G20.13)') acceptation, i_path, direction, b_opt(0), b_opt(1), b_wr_opt(0), b_wr_opt(1)&
       &, mu_moy, mu_wrmc, mu_NC, mu_DC, W*erg2eV, Wprecedent*erg2eV, xprob, 1-xprob

       !write(*,*) i_path, direction, W*erg2eV, Wprecedent*erg2eV, xprob, 1-xprob, acceptation, Wprecedent,&
       !        &  W

       Wprecedent = Wprec

       !call analyse_montecarlo(atconf_n,cells_n,boxmcgc, 'UO2_syst_n_after_test')
       !call caltabtC(cells_nplus1,atconf_nplus1,lperiod,boxmcgc)
       !call analyse_montecarlo(atconf_nplus1,cells_nplus1,boxmcgc, 'UO2_syst_nplus1_after_test')
       !call analyse_montecarlo(atconf_nplus1,cells_nplus1,boxmcgc, 'syst_UO2nplus1_finboucle')
       !call analyse_montecarlo(config_atom_n(ipch),config_cells_n(ipch),boxmcgc, 'UO2_syst_n_aftertest')
    end if !fin master general

    it = it +1

    if (direction == 0) then
       direction = 1
    else
       direction = 0
    end if

 END DO !end do sur la boucle des chemins

 if (lmegamaster) then
    write(*,*) ' taux d acceptation final   : ', acceptance_rate,  ' %'
    write(*,*) ' taux d acceptation alpha 0 : ', acceptance_rate_0,' %'
    write(*,*) ' taux d acceptation alpha 1 : ', acceptance_rate_1,' %'
    !call analyse_montecarlo(atconf_nplus1,cells_nplus1,boxmcgc, 'syst_UO2nplus1_out')
    write(*,'(A, G15.7)') 'mu_moy',mu_moy ,'mu_wrmc', mu_wrmc, 'mu_NC', mu_NC, 'mu_DC', mu_DC
 end if
stop
end subroutine montecarlo


subroutine ajout_retrait(direc)

 implicit none

 !-----------------------------------------------
 !   G l o b a l   P a r a m e t e r s
 !-----------------------------------------------
 !  type(atom_config_d)::atconf_N, atconf_Nplus1 
 !  type(cell_config)::cel_N, cel_Nplus1
 !  type(box_config)::box
 integer :: direc 
 !-----------------------------------------------
 !   L o c a l   P a r a m e t e r s
 !-----------------------------------------------
 real(double), dimension(3,1) :: cart_vec_nplus1
 integer :: indice, i,nag,rgcib,rgem,iloc
 logical ::ldistrib,lchange


 lperiod = .true.

 !######################################### Direction 0 vers 1 (ajout) ##########################################

 if (direc == 0) then ! ajout d'une particule en N+1
    if (lbigmaster) then
       !tirer une position aleatoire pour le N+1eme atome
       call atom_supp(cart_vec_nplus1)
       !cart_vec_nplus1(1,1) = 0.03125 +0.125 !0.61999346353817297                
       !cart_vec_nplus1(2,1) = 0.03125 +0.125 !0.73351583558252054                 
       !cart_vec_nplus1(3,1) = 0.03125 +0.125 !0.18351953714045011
       call cryst_to_cart(1,cart_vec_nplus1,boxmcgc%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m
       !write(*,*) 'atome supplementaire', cart_vec_nplus1
       !copie du syst n dans n+1 
       !call analyse_montecarlo(atconf_n,cells_n,boxmcgc, 'UO2_syst_n_before')
       !call analyse_montecarlo(atconf_nplus1,cells_nplus1,boxmcgc, 'UO2_syst_nplus1_before')
       call boucle_copy_atom(atconf_N,atconf_Nplus1, sens= .false.)   
       !call analyse_montecarlo(atconf_n,cells_n,boxmcgc, 'UO2_syst_n_after')
       !call analyse_montecarlo(atconf_nplus1,cells_nplus1,boxmcgc, 'UO2_syst_nplus1_after')     
       !addition de la n+1eme particule
       atconf_Nplus1%xp(1:3,atconf_Nplus1%im) = cart_vec_nplus1(1:3,1)
       atconf_Nplus1%xpp(1:3,atconf_Nplus1%im) = cart_vec_nplus1(1:3,1)
       atconf_Nplus1%fp(1:3,atconf_Nplus1%im) = 0
       atconf_Nplus1%ityp(atconf_Nplus1%im) = 1
       atconf_Nplus1%ielat(atconf_Nplus1%im) = -1
       !atconf_Nplus1%vp(1:3,atconf_Nplus1%im) = 0
       call init_vitesse(atconf_Nplus1,param = 1)
       nag=maxval(atconf_Nplus1%num_at_glob(1:atconf_Nplus1%im-1))
       atconf_Nplus1%num_at_glob(atconf_Nplus1%im) = nag+1

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
       !call indice_alea(atconf_Nplus1,indice)
       ! write(*,*) 'indice et coord atome a retirer', indice,  atconf_Nplus1%xp(:,indice)
 
       call atom_biais(atconf_Nplus1,indice)
       call atconf_Nplus1%switch_atom(indice,atconf_Nplus1%im)
       call calcul_proba

       !on copie les N nouveaux premiers atomes du syst N+1 dans le systeme N
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


subroutine calcul_proba
 implicit none
 real(double), dimension(atconf_Nplus1%imm) :: proba
 integer :: i,j
 real(double) :: dist_tot, alpha, sum_norm, tot
 real(double), dimension(3,1) :: coord

 sum_norm = 0.0
 DO i=1, atconf_Nplus1%im
    if (atconf_Nplus1%ityp(i) == 1) then ! si l'atome est un oxygene
       dist_tot = 0.0
       alpha = 0.0
       coord(:,1) =  atconf_Nplus1%xp(:,i)
       call calcul_dist(coord, dist_tot) !calcul des distances
       !write(*,*)  dist_tot
       call fct_alpha(dist_tot, alpha) !passage dans la fct alpha
       sum_norm = sum_norm + dist_tot * alpha
       proba(i) = dist_tot * alpha
    else
       proba(i) = 0.0
    end if !si oxygene ou uranium

 END DO ! boucle atomes

 proba(:) = proba(:) / sum_norm
 atconf_Nplus1%proba(:) = proba(:)

!!$ !verification que la somme est bien = à 1
!!$ tot = 0.0
!!$ DO i=1, atconf_Nplus1%im
!!$    if (atconf_Nplus1%ityp(i) == 1) then
!!$       tot = tot + proba(i)
!!$       !write(*,*) proba(i)
!!$    end if
!!$ end do
 !write(*,*) 'sum proba', tot

end subroutine calcul_proba


subroutine calcul_dist(coord_atom, dist)
 implicit none
 real(double) :: dist, m, n
 real(double), dimension(3,1) :: coord_atom, ref, distance
 integer :: j, k
 !attention, routine exacte uniquement pour les boites 4x4x4

 !write(*,*) 'coord atom avant', coord_atom(:,1)
 call cryst_to_cart(1,coord_atom,boxmcgc%bg,-1)
 !write(*,*) 'coord atom apres', coord_atom(:,1)
 DO j=1,3 !calculer la distance au site interstitiel 'parfait' le plus proche
    m = 2
    do k=0, 7
       n = abs((2.0*k+1)*0.0625 - coord_atom(j,1))
       if (n .lt. m) then
          ref(j,1) = (2.0*k+1.0)*0.0625
          m = n
       end if
    end do
    !write(*,*) ref(j,1), coord_atom(j,1)
    distance(j,1) = ref(j,1) - coord_atom(j,1)
    if (distance(j,1) .gt. 0.5) then
       distance(j,1) = distance(j,1) -1
    end if !CP si >0.5
    if (distance(j,1) .lt. -0.5) then
       distance(j,1) = distance(j,1) +1
    end if !CP si <-0.5
    dist = dist + (distance(j,1)*boxmcgc%at(j,j))**2
 END DO !boucle sur les coord
 dist = dsqrt(dist)*1E8

end subroutine calcul_dist


subroutine fct_alpha(dist, dist_alpha)
 implicit none
 real(double) :: dist, dist_alpha

 dist_alpha = 1.0-(1.0/( (dexp( ((dist/0.8)-1.0) * 3.0)) +1.0) )
end subroutine fct_alpha

subroutine calcul_chemin(travail, proba,Wp, dir, the)
real(double), dimension(nparapath+1) :: proba
real(double), dimension(nparapath) :: travail
real(double) :: Wp, sum_expW, beta, the
integer :: i, dir

beta = 1.0/(bk*Text)
sum_expW = 0.0

!calcul de la somme des expW
DO i = 1, nparapath
  sum_expW = sum_expW + exp(beta*(dir-the)*(travail(i)-Wp))
END DO

!attribution d'une proba pour chaque chemin
DO i = 1, nparapath
  proba(i) = exp(beta*(dir-the)*(travail(i)-Wp)) / (1.0 + sum_expW)
END DO
!proba de tirer l'ancien chemin
proba(nparapath+1) = 1.0 / (1.0 + sum_expW)
 
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

subroutine init_vitesse(config, param) !juste pour le N+1eme atome
 implicit none
 type(atom_config_mc)::config
 real(double) :: v0, v1, z1, z2, z3, z4
 integer :: i, param ! if param =1 alors v0 défini avec tinit (juste pour l'initialisation), sinon v0 défini avec Text  

 if (param==0) then 
    v0 = sqrt(2.D0*bk*tinit)
 else 
    v0 = sqrt(2.D0*bk*Text)
 end if

 i = config%im

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

subroutine indice_alea(config, ind)

 implicit none

 type(atom_config_mc)::config
 integer :: ind

 real(double) :: rand


 call random_number(rand)
 ind = ( (config%im - 1) * rand ) + 1

 do while (config%ityp(ind) > 1)
    call random_number(rand)
    ind = ( (config%im - 1) * rand ) + 1 
 end do

 !  write(*,*) 'indice atome supprimé', ind

end subroutine indice_alea


subroutine atom_biais(config, ind)
   implicit none
   type(atom_config_mc)::config
   integer :: ind, i
   real(double) :: rand, somme

   somme = 0.0
   call random_number(rand)
   !write(*,*) 'rand',rand
   DO i=1, config%im
     if (config%ityp(i) == 1) then
       if (somme .lt. rand) then
         somme = somme + config%proba(i)
         !write(*,*) 'somme', somme, 'ind', i
         ind = i
       end if ! if sur les sommes
     end if ! si sur les oxygene
   END DO !boucle sur les atomes

end subroutine atom_biais


subroutine analyse_montecarlo(atdml,celndm,box,name_file)

 implicit none

 type(atom_config_mc)::atdml
 type(cell_config)::celndm
 type(box_config)::box

 character(len=*) :: name_file
 if (lbigmaster) then


    if (itetemp>0) then
       if (mod(it,itetemp)==0) then
          call calctemp (temp,kine,atdml,celndm,latcomp=.true.)
       end if
    end if
    if (iterasmol>0) then
       if (mod(it,iterasmol)==0) then

          call rasmolT(atdml,box,it,namefr=name_file,latcomp=.true.)
       end if
    end if

 end if
end subroutine analyse_montecarlo

subroutine calcul_U(atdml,celndm,box,potist,U_ini)

 use tempinstT_mod,only:tempinstT
 implicit none

 type(atom_config_mc)::atdml
 type(cell_config)::celndm
 type(box_config)::box
 real(double) :: U_ini, potist,Ti

 if (lbigmaster) then
    Ti=tempinstT(atdml,kine,latcomp=.true.)

    U_ini = potist + kine
 end if
end subroutine calcul_U




subroutine atom_supp(vecteur)

 implicit none

 real(double), dimension(3,1) :: vecteur
 real(double) :: x_nplus1, y_nplus1, z_nplus1 !position initiale aleatoire de la N+1eme particule
 integer::i
 do i=1,rang+1
    call random_number(x_nplus1)
    call random_number(y_nplus1)
    call random_number(z_nplus1)
 end do
! write(6,*)'PLUS1',rang,x_nplus1,y_nplus1,z_nplus1
 vecteur(1,1) = x_nplus1
 vecteur(2,1) = y_nplus1
 vecteur(3,1) = z_nplus1

 !  write(*,*) 'atome supplementaire', vecteur
end subroutine atom_supp


subroutine noise(var)

 implicit none

 real(double)  :: var(3,atconf_nplus1%im)
 real(double)  :: u_1, u_2, v1, v2, r
 integer :: i, ic


 DO i=1, atconf_nplus1%im
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



subroutine langevin( direc, protocol)
 implicit none

 character(len=3), intent(in) :: protocol
 integer :: direc 
 integer :: i,ic, ip, tot
 real(double) :: Ek_n, Ek_n_plus1, Ek_n_1s4, Ek_n_3s4, dQeff, Qeff,&
      &dWeff, dWork
 real(double) :: U_0, U_1, U_l_n_m1, U_l_n, H_l_n, H_l_n_m1, H_l_ini

 real(double) :: beta
 real(double)  :: Gl(3,atconf_N%im+1)
 real(double)::rga, rga_s4

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
    U_l_n = (1.d0-lambda_mc)*potist_n + lambda_mc*potist_nplus1

    H_l_ini = Ek_n + U_l_n 
    H_l_n   = H_l_ini

    !write(*,*) 'lambda_mc ' ,lambda_mc,'Ek_n' , Ek_n, 'U_l_n', U_l_n, 'H_l_ini', H_l_ini
 end if ! Master général

 aux(:ntyp) = tstep/(cm(:ntyp)*2.d0)


 !########################################################################################################################
 !                               Ajout d'une particule N+1: système N vers N+1 - direction = 0
 !########################################################################################################################
 if (direc == 0) then
    DO ip = 1, pas_lambda_mc

       !incrémentation de lambda
       call lambda(direc,ip, protocol)
       !lambda_mc = dble(ip)/dble(pas_lambda_mc)

       if (lbigmaster) then ! Master général
          !call analyse_montecarlo(atconf_nplus1,cells_nplus1,boxmcgc, 'syst_UO2nplus1_direc0')
          Ek_n = 0.0
          Ek_n_plus1 = 0.0  
          Ek_n_1s4 = 0.0
          Ek_n_3s4 = 0.0
          !write(*,*) 'lambda_mc ' ,lambda_mc
          ! faire le pas de langevin (velocity verlet) pour determiner les nouvelles forces et positions

          ! step 1 First half-step velocities update, v(t) -> v(t+dt/2)
          timel = timel+tstep
          rga=exp((-gamlg)*tstep/2)
          rga_s4 = exp((-gamlg)*tstep/4)
          call noise(Gl)
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
             atconf_Nplus1%xpp(1:3,i)=atconf_Nplus1%xp(1:3,i)
             atconf_Nplus1%xp(1:3,i) = atconf_Nplus1%xp(1:3,i) + tstep*atconf_Nplus1%vp(1:3,i)
          END DO

          !recopier les nouvelles positions dans le syst N
          DO i=1, atconf_N%im
             atconf_N%xpp(1:3,i) = atconf_Nplus1%xpp(1:3,i)
             atconf_N%xp(1:3,i) = atconf_Nplus1%xp(1:3,i)
          END DO

          !conditions periodiques 
          if (lperiod)    then
             call periodbox(boxmcgc,atconf_N)
             call periodbox(boxmcgc,atconf_Nplus1)
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
       call calfoMCGC(iloc,lchange,ldistrib)


       if (lbigmaster) then
          !mise a jour de U_l_n = (1-lambda_mc)*U_0 + lambda_mc*U_1
          U_l_n = (1.0-lambda_mc)*potist_n + lambda_mc*potist_nplus1
          !write(*,*) potist_n, potist_nplus1
          !affichage temperature
          it = it +1

          call noise(Gl)
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

       end if !master general

       if (lbigmaster) then !master general
          !calcul des energies et travail et chaleur efficaces
          U_l_n_m1 = U_l_n
          H_l_n_m1 = H_l_n
          H_l_n    = Ek_n_plus1  + U_l_n
          dWork = H_l_n - H_l_n_m1
          Work = Work + dWork
          dQEff  = (Ek_n_1s4-Ek_n) + (Ek_n_plus1-Ek_n_3s4)
          QEff   = QEff + dQEff
          dWEff  = H_l_n - H_l_n_m1 - dQEff
          WEff   = WEff + dWEff
          !call analyse_montecarlo(atconf_nplus1,cells_nplus1,boxmcgc, 'syst_UO2nplus1_direc0')
          if (protocol == 'MCP') then
             !call analyse_montecarlo(atconf_nplus1,cells_nplus1,boxmcgc, 'syst_UO2nplus1_direc0')
             !write(*,*) 'lambda_mc ' ,lambda_mc, 'Ek_n_plus1', Ek_n_plus1, 'U_l_n',&
             !          & U_l_n, 'H_l_n', H_l_n, 'WEff', WEff, 'dWEff', dWEff
             !write(*,*) 'lambda_mc ' ,lambda_mc, 'Ek_n_plus1', Ek_n_plus1, 'U_l_n',&
             !            & U_l_n, 'H_l_n', H_l_n, 'Work', Work, 'dWork', dWork
          end if
       end if

    END DO

 endif



 !########################################################################################################################
 !                            Déletion d'une particule N+1: système N+1 vers N - direction = 1
 !########################################################################################################################

 if (direc == 1) then
    DO ip = 1, pas_lambda_mc

       !incrémentation de lambda
       call lambda(direc,ip, protocol)
       !lambda_mc = 1.d0 - (dble(ip)/dble(pas_lambda_mc))

       if (lbigmaster) then
          !call analyse_montecarlo(atconf_nplus1,cells_nplus1,boxmcgc, 'syst_UO2nplus1_direc1_av_lang')
          Ek_n = 0.0
          Ek_n_plus1 = 0.0  
          Ek_n_1s4 = 0.0
          Ek_n_3s4 = 0.0
          ! faire le pas de langevin (velocity verlet) pour determiner les nouvelles forces et positions
          ! step 1 First half-step velocities update, v(t) -> v(t+dt/2)
          timel = timel+tstep
          !write(*,*) 'lambda_mc ' ,lambda_mc
          !write(6,*) 'avant langevin atconf_Nplus1%vp(:,23)=' ,atconf_Nplus1%vp(:,23)
          rga=exp((-gamlg)*tstep/2)
          rga_s4 = exp((-gamlg)*tstep/4)
          call noise(Gl)
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
             atconf_Nplus1%xpp(1:3,i)=atconf_Nplus1%xp(1:3,i)
             atconf_Nplus1%xp(1:3,i) = atconf_Nplus1%xp(1:3,i) + tstep*atconf_Nplus1%vp(1:3,i)
          END DO

          !recopier les nouvelles positions dans le syst N
          DO i=1, atconf_N%im
             atconf_N%xpp(1:3,i) = atconf_Nplus1%xpp(1:3,i)
             atconf_N%xp(1:3,i) = atconf_Nplus1%xp(1:3,i)
          END DO

          !conditions periodiques 
          if (lperiod)    then
             call periodbox(boxmcgc,atconf_N)
             call periodbox(boxmcgc,atconf_Nplus1)
          end if
       end if

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
       call calfoMCGC(iloc,lchange,ldistrib)

       if (lbigmaster) then

          !mise a jour de U_l_n = (1-lambda_mc)*U_0 + lambda_mc*U_1
          U_l_n = (1-lambda_mc)*potist_n + lambda_mc*potist_nplus1
          !write(*,*) potist_n, potist_nplus1
          !affichage temperature
          it = it +1

          ! Second half-step velocities update, v(t+1/2dt) -> v(t+dt)
          call noise(Gl)
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

          !repartir les nouvelles positions et forces dans les syst N et N+1
          !Systeme a N

          DO i=1, atconf_N%im
             atconf_N%vp(1:3,i) = atconf_Nplus1%vp(1:3,i) 
          END DO

       end if

       if (lbigmaster) then
          !calcul des energies et travail et chaleur efficaces
          U_l_n_m1 = U_l_n
          H_l_n_m1 = H_l_n
          H_l_n    = Ek_n_plus1  + U_l_n
          dWork = H_l_n - H_l_n_m1
          Work = Work + dWork           
          dQEff  = (Ek_n_1s4-Ek_n) + (Ek_n_plus1-Ek_n_3s4)
          QEff   = QEff + dQEff
          dWEff  = H_l_n - H_l_n_m1 - dQEff
          WEff   = WEff + dWEff
          if (protocol == 'MCP') then
             !call analyse_montecarlo(atconf_nplus1,cells_nplus1,boxmcgc, 'syst_UO2nplus1_direc1')
             !write(*,*) 'lambda_mc ' ,lambda_mc, 'Ek_n_plus1', Ek_n_plus1, 'U_l_n',&
             !        & U_l_n, 'H_l_n', H_l_n, 'WEff', WEff, 'dWEff', dWEff
             !write(*,*) 'lambda_mc ' ,lambda_mc, 'Ek_n_plus1', Ek_n_plus1, 'U_l_n',&
             !            & U_l_n, 'H_l_n', H_l_n, 'Work', Work, 'dWork', dWork
          end if

       end if
    END DO


 endif


end subroutine langevin



subroutine lambda(dir, nstep, protocol_name)
 implicit none
 integer :: nstep, dir
 character(len=3), intent(in) :: protocol_name
 real(double) :: alpha

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

if (lparapath) then 
   if (mod(nprocs,2*nparapath).ne.0) then
      write(6,*)'nprocs/2*nparapath <>0 STOP'
      call MPI_FINALIZE(ierr)
      stop
   end if
   parapath%mpi_orig%nproc=nprocs
   parapath%mpi_orig%rank=rang
   parapath%nimage=nparapath
   call MPI_COMM_DUP(MPI_COMM_WORLD,parapath%mpi_orig%comm,ierr)
   call MPI_COMM_GROUP(parapath%mpi_orig%comm,parapath%mpi_orig%group,ierr)
   call commconstr(parapath)
!   call parapath%print(unit=1000)
  
else
   call mpi_world%print(unit=50+rang)
   call initparapuresp(parapath,rang,mpi_WORLD)
end if


paramcgc%mpi_orig%nproc= parapath%mpi_image%nproc ! =parapath%mpi_orig%nproc/nparapath
!paramcgc%mpi_orig%comm= parapath%mpi_image%comm ! =parapath%mpi_orig%nproc/nparapath 
 paramcgc%mpi_orig%rank=parapath%mpi_image%rank
 paramcgc%nimage=2

!!$ paramcgc%mpi_orig%nproc=nprocs
!!$ lbigmaster=rang
!!$ !    paramcgc%mpi_orig%group=grp_world
!!$ paramcgc%nimage=2
!!$ !    paramcgc%mpi_orig%comm=MPI_COMM_WORLD

 call MPI_COMM_DUP(parapath%mpi_image%comm,paramcgc%mpi_orig%comm,ierr)
! call MPI_COMM_DUP(MPI_COMM_WORLD,paramcgc%mpi_orig%comm,ierr)
 call MPI_COMM_GROUP(paramcgc%mpi_orig%comm,paramcgc%mpi_orig%group,ierr)
  call parapath%print(rang)
  
 call commconstr(paramcgc)

 myidsp=paramcgc%mpi_image%rank
 call MPI_COMM_free(mpi_comm_space,ierr)
 MPI_COMM_space=paramcgc%mpi_image%comm
 call comm_space%init(MPI_COMM_SPACE)
 nprocspace=paramcgc%mpi_image%nproc
 if (nprocspace==1) parallele=.false.
 lbigmaster=parapath%lmaster
 lmaster=paramcgc%lmaster
 lmegamaster=.false.
 if (parapath%mpi_orig%rank==0) lmegamaster=.true.

 call paramcgc%print(rang+100)
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

subroutine distat(xi,x0,box,dist)
 type(box_config),intent(in)::box
 real(double), dimension(3),intent(in)::xi,x0
 real(double),dimension(3)::dx
 real(double)::dist
 real(double),dimension (3,2)::xat
 integer::ns=2
 xat(:,1)=xi(:)
 xat(:,2)=x0(:)
 call cryst_to_cart (ns,xat,box%bg,-1)
 dx(1:3)=xat(1:3,1)-xat(1:3,2)
 WHERE ( (dx.GT.0.5d0).OR.(dx.LT.-0.5d0) )
    dx(1:3) = dx(1:3) - Dble(Nint(dx(1:3)))
 END WHERE
 dx = MatMul(box%at,dx)
 dist = sqrt(Sum( dx(1:3)**2 ))
 return
end subroutine distat



subroutine initNP1(ipp) !PARAPATH DEFINIR LES POINTEURS atconf_nplus1 et atconf_n

  integer,intent(in)::ipp
  real(double), dimension(3,1) :: cart_vec_nplus1
  real(double)::distati
  integer::i
  character :: extension*4
  logical ::lc2d
  !definir le systeme a N+1 en tirant une position aleatoire pour le N+1eme atome

  if (lbigmaster) then

     call atom_supp(cart_vec_nplus1)
     !write(*,*) 'cart_vec_nplus1', cart_vec_nplus1(:,1)
     !cart_vec_nplus1(1,1) = 0.03125 +0.125
     !cart_vec_nplus1(2,1) = 0.03125 +0.125 
     !cart_vec_nplus1(3,1) = 0.03125 +0.125 
     call cryst_to_cart(1,cart_vec_nplus1,boxmcgc%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m
     !FAIT DANS atom_supp
     !    call cryst_to_cart(1,cart_vec_nplus1,boxmcgc%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m

     !copie du syst n dans n+1 
     call atconf_nplus1%init(atconf_n%im+1,atconf_n%imm,atconf_n%ltabvois)
     atconf_nplus1%ltabvois=atconf_n%ltabvois
     !call atconf_n%copy_config(atconf_nplus1,lrescl=.false.)
     call boucle_copy_atom(atconf_n,atconf_nplus1, sens= .false.)
     !addition de la n+1eme particule
     atconf_nplus1%xp(1:3,atconf_nplus1%im) = cart_vec_nplus1(1:3,1)
     atconf_nplus1%fp(1:3,atconf_nplus1%im) = 0
     atconf_nplus1%xpp(1:3,atconf_nplus1%im) =     atconf_nplus1%xp(1:3,atconf_nplus1%im) 
     atconf_nplus1%ityp(atconf_nplus1%im) = 1
     atconf_nplus1%num_at_glob(atconf_nplus1%im) = atconf_nplus1%im
     !atconf_Nplus1%vp(1:3,atconf_Nplus1%im) = 0
     call init_vitesse(atconf_nplus1,param = 0)
     !copie de cell puis caltabtC pour redecouper avec la n+1eme particule
     call cells_nplus1%init(cells_n%nox,cells_n%noy,cells_n%noz, cells_n%natperc)
     call cells_n%copy_cell(cells_nplus1)
  else
     call atconf_nplus1%init(atconf_n%im+1,atconf_n%imm,atconf_n%ltabvois)
     call cells_nplus1%init(cells_n%nox,cells_n%noy,cells_n%noz, cells_n%natperc)
     call cells_n%copy_cell(cells_nplus1)
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
        !?          call maj_atomes_frt_ftm(atconf_n,cells_n)

     else !procs N+1
        call atconf_nplus1%send2all(0,paramcgc%mpi_image)
        call caltabtC(cells_nplus1,atconf_nplus1,lperiod,boxmcgc)
        call init_voisinage(cells_nplus1,pscgc)
        !?          call maj_atomes_frt_ftm(atconf_nplus1,cells_nplus1)
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

        if (lbigmaster) then
           write(6,*)'write configuration N+1  to confNP1.lmp'
           write(extension,'(i4.4)') ipp

           call config2data (atconf_nplus1%imm,atconf_nplus1%im,&
                atconf_nplus1%xp,atconf_nplus1%ityp,boxmcgc%at,ntyp,filename='confNP1.'//extension//'.lmp') ! PARAPATH CHANGER LE NOM AVEC INDICE DE LA BOITE
        end if

#ifdef PARA
#ifdef LAMMPS_VERSION
        if (lparapath) then
           !          if(parapath%image+1=ipp) then ! assuré par lc2d
           if(paramcgc%image==0) then !procs N
              firsttime_lammps=.true.
              allocate (posa(3*atconf_n%im),  forca(3*atconf_n%im))

              call init_lammps('in.lammps.N')
           else !procs N+1
              firsttime_lammps=.true.
              allocate (posa(3*atconf_nplus1%im),  forca(3*atconf_nplus1%im))
              write(extension,'(i4.4)') ipp
              namef='in.lammps.'//extension//'.NP1'
              call init_lammps(namef)
           end if
           !          end if
        else

           if(paramcgc%image==0) then !procs N
              firsttime_lammps=.true.
              allocate (posa(3*atconf_n%im),  forca(3*atconf_n%im))

              call init_lammps('in.lammps.N')
           else !procs N+1
              firsttime_lammps=.true.
              allocate (posa(3*atconf_nplus1%im),  forca(3*atconf_nplus1%im))
              call init_lammps('in.lammps.NP1')
           end if
        end if


#else
     write(6,*)'Ipotentiel<0 (lammps) et NON LAMMPS_VERSION : stop'
     call MPI_FINALIZE(ierr)
     stop
#endif

#else
     write(6,*)'Ipotentiel<0 (lammps) et NON para en MCGC : stop'
     stop
#endif       

  end if
end if

end subroutine initNP1



subroutine calfoMCGC(iloc,lchange,ldistrib)
 integer,intent(in)::iloc
 logical, intent(in)::lchange,ldistrib
 integer::rgcib,rgem,i


#ifdef PARA
 if(paramcgc%image==0) then !procs N
    if (iloc==1) call initloc(atconf_n,cells_n,atmcgcloc,cellmcgcloc,boxmcgc,paramcgc,rumax,lperiod&
         &,psc=pscgc,ldistrib=ldistrib) !initloc contient caltabtc sur atloc
    call pointer_caltabt_calfo(sig,potist_n,atconf_n,cells_n,boxmcgc,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,atconf_n%ltabvois,it,itetabvois,lchg=lchange,psc=pscgc)
 else !procs N+1

    if (iloc==1)call initloc(atconf_nplus1,cells_nplus1,atmcgcloc,cellmcgcloc,boxmcgc,paramcgc,rumax,&
         &lperiod,psc=pscgc,ldistrib=ldistrib) !initloc contient caltabtc sur atloc
    call pointer_caltabt_calfo(sig,potist_nplus1,atconf_nplus1,cells_nplus1,boxmcgc,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,atconf_nplus1%ltabvois,it,itetabvois,lchg=lchange,psc=pscgc)
 end if
 !en ce point chacun des deux masters a les forces de son paquet datomes
 if (lmaster) then ! on est dans l'un des 2 masters7
    rgcib=0;rgem=1
    if(paramcgc%image==1) then !on est dans le master de N+1
       call atconf_nplus1%send2proc(rgcib,paramcgc%mpi_master,'f')
       !          call MPI_SEND(potist_nplus1, 1,NDM_MPI_REAL_DOUBLE,rgcib,1000,paramcgc%mpi_master%comm,ierr)
       call paramcgc%mpi_master%send(potist_nplus1,rgcib,1000)
    else !on est dans le master de N qui est le master général
       call atconf_nplus1%recv(rgem,paramcgc%mpi_master,'f')
       call paramcgc%mpi_master%recv(potist_nplus1,rgem,1000)
       !          call MPI_RECV(potist_nplus1, 1,NDM_MPI_REAL_DOUBLE,rgem,1000,paramcgc%mpi_master%comm,status,ierr)
    end if
 end if
 !en ce point le master général (rang_orig=0) a les forces de N et N+1    
#else
 call pointer_caltabt_calfo(sig,potist_n,atconf_n,cells_n,boxmcgc,atmcgcloc,cellmcgcloc,paramcgc,&
      &lperiod,atconf_n%ltabvois,it,itetabvois,lchg=lchange,psc=pscgc)
 call pointer_caltabt_calfo(sig,potist_nplus1,atconf_nplus1,cells_nplus1,boxmcgc,atmcgcloc,cellmcgcloc,paramcgc,&
      &lperiod,atconf_nplus1%ltabvois,it,itetabvois,lchg=lchange,psc=pscgc)


#endif

 if (lbigmaster) then
    DO i=1,atconf_n%im
       atconf_nplus1%fp(:,i) = (1-lambda_mc)*atconf_n%fp(:,i) + lambda_mc*atconf_nplus1%fp(:,i)
    END DO
    atconf_nplus1%fp(:,atconf_nplus1%im) = lambda_mc*atconf_nplus1%fp(:,atconf_nplus1%im)

    !Egalisation des forces pour les deux systemes
    DO i=1,atconf_n%im
       atconf_n%fp(:,i) = atconf_nplus1%fp(:,i)
    END DO
 end if!end master general

end subroutine calfoMCGC

subroutine init_atom_config_mc(atconf,imin,immin,ltabvois,nvois,rvois,lreallocate)
 class(atom_config_mc),intent(inout)::atconf
 !type(atom_config_mc),intent(inout)::atconf
 integer,intent(in):: imin
 logical,optional, intent(in)::ltabvois,lreallocate
 integer, optional::nvois,immin
 real(double),optional::rvois
 logical :: lrealloc

 lrealloc=.false.
 if (present(lreallocate))then
    lrealloc=lreallocate
 end if
 !initialisation de la partie atom_config_d
 call atconf%atom_config_d%init(imin,immin,ltabvois,nvois,rvois,lreallocate) 
 !initialisation de la partie mc ajoutée
 if ((lrealloc).and.(allocated(atconf%proba)))then
    deallocate(atconf%proba)
 end if
 if (.not.allocated(atconf%proba))then
    allocate(atconf%proba(atconf%imm))
 end if
 atconf%proba=0
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
       atcible%proba(1:atsource%imm)=atsource%proba(1:atsource%imm)
    end select
 end select
end subroutine copy_config_mc



subroutine copy_atom_mc(atsource,i,atcible,j,lextend)
 implicit none
 class(atom_config_mc), intent(in)::atsource
 !type(atom_config_mc),intent(in)::atsource
 integer,intent(in):: i
 class(atom_config), intent(inout)::atcible
 !type(atom_config_mc), intent(inout)::atcible
 integer,intent(in):: j
 logical , optional, intent(in) :: lextend
 logical::let
 let=.false.
 if (present(lextend)) then
    let=lextend
 end if
 call atsource%atom_config_d%copy_atom(i,atcible,j,let)
 select type(atcible)
    class is (atom_config_mc)
    select type (atsource)
       class is (atom_config_mc)
       atcible%proba(j) = atsource%proba(i)
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
 intermediaire = atsource%proba(ind_switch_1)
 atsource%proba(ind_switch_1) =  atsource%proba(ind_switch_2)
 atsource%proba(ind_switch_2) = intermediaire

end subroutine switch_atom_mc


end module montecarlo_mod
