module montecarlo_mod
  USE gen_com_m,only:  lperiod, tstep, timel, tstep, sig, potist,itetabvois,&
       & iterasmol,itetemp, temp, kine, pi, bk, Text, gamlg,one,pi,text,tinit,&
       &lspaceNDM,rang,it,firsttime_lammps,posa,forca
!  USE tab_imm_m, only:xp, xpp, fp, vp, num_at_glob, ityp, ielat, iwmax
  USE atomconfig,only:atom_config,atom_config_d, config2ndm, switch_atom
  USE period_mod,only: period 
  USE cellconfig, only:cell_config, cellconfig2ndm, caltabtC
  USE var_pot,only:ntyp,cm,gamlt
  USE calfo_mod,only: calfo
  USE T_kind_param_m, ONLY:  double
  USE cryst_to_cart_mod, ONLY: cryst_to_cart
  USE caltabi_mod,only: caltabi
  USE boxconfig,only:box_config,periodbox
  USE rasmolT_mod,only: rasmolT
  use paraconfig,only:para_config,commconstr
#ifdef PARA
  use Tpara,only:grp_world,nprocs,myidsp,MPI_COMM_space,nprocspace,ierr,mpi_comm_world,NDM_MPI_REAL_DOUBLE
  use mod_para,only:maj_atomes_frt_ftm
  USE init_vois_mod,only: init_voisinage
#else
  use Tpara,only:myidsp,nprocspace
#endif
  use read_val,only:rvois,ltabvois
  use var_pot,only:ipotentiel,rumax
  USE parautils,only:initloc,pointer_caltabt_calfo
  USE calctemp_mod,only:calctemp
  USE constrconf_mod,only:config2data
#ifdef LAMMPS_VERSION
  use lammps_util_mod
  use vars_lammps
#endif  
  implicit none

  type(para_config)::paramcgc

 real(double)::distminat
  type(atom_config_d)::config_atom_n !type derive atom_config du systeme a n atomes
  type(atom_config_d)::config_atom_nplus1 !type derive atom_config du systeme a n+1 atomes

  type(cell_config):: cells_n !type derive cell_config du systeme a n atomes
  type(cell_config):: cells_nplus1 !type derive cell_config du systeme a n+1 atomes

  type(box_config)::boxmcgc
  
  integer::  pas_lambda_mc
  real(double) :: lambda_mc !lambda compris entre 0 et 1

  real(double), dimension(3,3) :: sig_n, sig_nplus1
  real(double) :: potist_n, potist_nplus1
  real(double) :: Weff
  logical::lmaster !(true= master du système N ou du système N+1)
    class(atom_config),pointer::atmcgcloc
    type(cell_config),pointer::cellmcgcloc
    type(cell_config),target:: cellcible ! ne sert qu'à faire pointer cellnebloc sur quelquechose
    type(atom_config_d),target::atcible
    integer::rgcib,rgem !(cible et emeteur e confN+1)
contains 

  
  subroutine montecarlo


    implicit none

    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------

    type(atom_config_d):: config_atom_old_0, config_atom_old_1, config_atom_new_0, config_atom_new_1 !config intermediaire pour suivre l'evolution des systemes: 0 -> syst N, 1 -> syst N+1.



    integer :: i,ic
    logical :: lextend

    integer :: direction
    integer :: n_accepted, n_accepted_0, n_accepted_1
    integer :: n_gen, n_gen_0, n_gen_1

    integer :: i_path, n_path

    real(double) :: theta
    real(double) :: beta

    real(double) :: W, Wprec, xprob, xalea
    real(double) :: ln_xalea, ln_Wprec, ln_W, ln_xprob
    real(double) :: acceptance_rate, acceptance_rate_0, acceptance_rate_1

    integer, dimension(12) :: seed 
    integer :: n



    !########################################################################################################################
    !                                             Initialisation
    !########################################################################################################################

    !initialisation variables 
    !pour le premier chemin: sens positif, d'ajout d'une particule et acceptation     
#ifdef PARA
    cellmcgcloc=>cellcible
    atmcgcloc=>atcible
    lmaster=paramcgc%lmaster
    if ((paramcgc%npim.gt.1).and.(lspaceNDM.eqv..true.)) then
       call init_voisinage(cells_n)
    end if
#else
    lmaster=.true.
#endif    


    direction = 0 ! direction = 0 on ajoute un atome, = 1 on retire un atome
    lambda_mc = 0.0

    lextend = .true.
    lperiod = .true.

    n_path = 100
    theta = 0.5
    beta = 1.0/(bk*Text)

    n_accepted   = 0
    n_accepted_0 = 0
    n_accepted_1 = 0
    n_gen   = 0
    n_gen_0 = 0
    n_gen_1 = 0
    acceptance_rate   = 0.0
    acceptance_rate_0 = 0.0
    acceptance_rate_1 = 0.0




#ifdef PARA
    if(paramcgc%image==0) then !procs N
       call initloc(config_atom_n,cells_n,atmcgcloc,cellmcgcloc,boxmcgc,paramcgc,rumax,lperiod) !initloc contient caltabtc sur atloc
       call pointer_caltabt_calfo(sig,potist,config_atom_n,cells_n,boxmcgc,atmcgcloc,cellmcgcloc,paramcgc,&
            &lperiod,config_atom_n%ltabvois,it,itetabvois,lchg=.false.)
    else !procs N+1
       call initloc(config_atom_nplus1,cells_nplus1,atmcgcloc,cellmcgcloc,boxmcgc,paramcgc,rumax,lperiod) !initloc contient caltabtc sur atloc
       call pointer_caltabt_calfo(sig,potist,config_atom_nplus1,cells_nplus1,boxmcgc,atmcgcloc,cellmcgcloc,paramcgc,&
            &lperiod,config_atom_nplus1%ltabvois,it,itetabvois,lchg=.false.)
    end if
    call MPI_Barrier(MPI_COMM_SPACE,ierr)
    !en ce point chacun des deux masters a les forces de son paquet datomes
    if (lmaster) then ! on est dans l'un des 2 masters
       rgcib=0;rgem=1
       if(paramcgc%image==1) then !on est dans le master de N+1
          call config_atom_nplus1%send2proc(rgcib,paramcgc%comm_master,'f')
       else !on est dans le master de N qui est le master général
          call config_atom_nplus1%recv(rgem,paramcgc%comm_master,'f')
       end if
    end if
    !en ce point le master général (rang_orig=0) a les forces de N et N+1    
#else
    call pointer_caltabt_calfo(sig,potist,config_atom_n,cells_n,boxmcgc,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,config_atom_n%ltabvois,it,itetabvois,lchg=.false.)
    call pointer_caltabt_calfo(sig,potist,config_atom_nplus1,cells_nplus1,boxmcgc,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,config_atom_nplus1%ltabvois,it,itetabvois,lchg=.false.)
    !    call caltabtC(cells_n,config_atom_n,lperiod,boxmcgc)
    !    call caltabtC(cells_nplus1,config_atom_nplus1,lperiod,boxmcgc)
    !    CALL CalFo(sig_n,potist_n,config_atom_n,cells_n,boxmcgc) 
    !    CALL CalFo(sig_nplus1,potist_nplus1,config_atom_nplus1,cells_nplus1,boxmcgc)

#endif
    !initialisation pour le dyn_vverlet
    !calcul des forces des systemes N et N+1
    ! melange des forces des deux systemes N et N+1
    if (paramcgc%rang_orig==0) then
       DO i=1,config_atom_n%im
          config_atom_nplus1%fp(:,i) = (1-lambda_mc)*config_atom_n%fp(:,i) + lambda_mc*config_atom_nplus1%fp(:,i)
       END DO
       config_atom_nplus1%fp(:,config_atom_nplus1%im) = lambda_mc*config_atom_nplus1%fp(:,config_atom_nplus1%im)

       !Egalisation des forces pour les deux systemes
       DO i=1,config_atom_n%im
          config_atom_n%fp(:,i) = config_atom_nplus1%fp(:,i)
       END DO


       ! sauvegarde du système
       call config_atom_n%copy_config(config_atom_old_0, lrescl=.true.)
    end if

    !pour le premier chemin: sens positif, d'ajout d'une particule et acceptation

    call langevin(config_atom_n,config_atom_nplus1,cells_n,cells_nplus1,boxmcgc,direction)
    if (paramcgc%rang_orig==0) then
       call config_atom_nplus1%copy_config(config_atom_old_1, lrescl=.true.)

       W = WEff
       xprob = 1
       Wprec = + W



       seed(1) = 152533
       call random_seed(PUT=seed(1:12))
    end if
    direction = 1
    !########################################################################################################################
    !                                             boucle sur lambda le long d'un chemin
    !########################################################################################################################
    
    DO i_path = 1, n_path ! boucle à faire pour tous les procs
       if (paramcgc%rang_orig==0) then !!master general
          call random_number(xalea)
          ln_xalea  = log(xalea)
          config_atom_nplus1%vp(:,:)   = - config_atom_nplus1%vp(:,:) !à chaque retour dans la boucle, on change de direction
          !initialisation de lambda
          if (direction == 0) then
             lambda_mc = 0.0
             call config_atom_n%copy_config(config_atom_new_0, lrescl=.true.)
          endif
          if (direction == 1) then
             lambda_mc = 1.0 
             call config_atom_nplus1%copy_config(config_atom_new_1, lrescl=.true.)
          endif
          !choisir l'at a retirer ou ajouter + preparation des syst N et N+1 pour etre prets pour le langevin (cad decoupage cellules + calcul forces + melange des forces - se fait dans cette sous routine)
       end if !master general
       call ajout_retrait(config_atom_n,config_atom_nplus1,cells_n,cells_nplus1,boxmcgc,direction)

       call analyse_montecarlo(config_atom_n,cells_n,boxmcgc, 'UO2_syst_n_before_test')
       call analyse_montecarlo(config_atom_nplus1,cells_nplus1,boxmcgc, 'UO2_syst_nplus1_before_test')

          !       end if
       ! pas de langevin
       call langevin(config_atom_n,config_atom_nplus1,cells_n,cells_nplus1,boxmcgc,direction)

       if (paramcgc%rang_orig==0) then
          if (direction == 0) then
             W = +WEff
             n_gen_0 = n_gen_0 + 1
             call config_atom_nplus1%copy_config(config_atom_new_1, lrescl=.true.)
          else
             W = -WEff
             n_gen_1 = n_gen_1 + 1
             call config_atom_n%copy_config(config_atom_new_0, lrescl=.true.)
          endif
          n_gen = n_gen + 1
          ln_Wprec  = (+beta*(direction-theta)*Wprec)
          ln_W      = (+beta*(direction-theta)*W)
          ln_xprob  = - dlog(1 + dexp(ln_Wprec-ln_W))
          xprob     = dexp(ln_xprob)

          write(*,*) 'Wprec', Wprec, 'ln_Wprec', ln_Wprec
          write(*,*) 'W', W, 'ln_W', ln_W

          if (ln_xprob > ln_xalea) then    
!!!!!!!!!!!!!!!!!!!!!!!!!! ACCEPTATION   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             write(*,*) 'ACCEPTATION, direction=', direction,& 
                  &'  LN_XPROB ', ln_xprob, '  XPROB ', xprob, '  XALEA ', xalea

             if (direction == 0) then
                n_accepted_0 = n_accepted_0 + 1
             else
                n_accepted_1 = n_accepted_1 + 1
             endif

             n_accepted = n_accepted + 1
             Wprec = + W 

             call config_atom_new_0%copy_config(config_atom_old_0, lrescl=.true.)
             call config_atom_new_1%copy_config(config_atom_old_1, lrescl=.true.)          

          else
!!!!!!!!!!!!!!!!!!!!!!!!!! REFUS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             write(*,*) ' REJECTION, direction=', direction, &
                  & '  LN_XPROB ', ln_xprob, '  XPROB ', xprob, '  XALEA ', xalea

             !on accepte le sens opposé - changer des signes des vitesses 
             config_atom_old_0%vp(:,:)   = - config_atom_old_0%vp(:,:)
             config_atom_old_1%vp(:,:)   = - config_atom_old_1%vp(:,:)

             if (direction == 0) then
                call config_atom_old_1%copy_config(config_atom_nplus1, lrescl=.true.)
                call caltabtC(cells_nplus1,config_atom_nplus1,lperiod,boxmcgc)
             end if

             if (direction == 1) then
                call config_atom_old_0%copy_config(config_atom_n, lrescl=.true.)
                call caltabtC(cells_n,config_atom_n,lperiod,boxmcgc)
             end if


          end if

!!$          if (direction == 0) then
!!$             direction = 1
!!$          else
!!$             direction = 0
!!$          end if
          acceptance_rate   = (real(n_accepted)/real(n_gen))*1.0d2
          acceptance_rate_0 = (real(n_accepted_0)/real(n_gen_0))*1.0d2
          acceptance_rate_1 = (real(n_accepted_1)/real(n_gen_1))*1.0d2

       end if
       call analyse_montecarlo(config_atom_n,cells_n,boxmcgc, 'UO2_syst_n_after_test')
       call analyse_montecarlo(config_atom_nplus1,cells_nplus1,boxmcgc, 'UO2_syst_nplus1_after_test')
       it = it +1
       
       if (direction == 0) then
          direction = 1
       else
          direction = 0
       end if
    END DO

    write(*,*) ' taux d acceptation final   : ', acceptance_rate,  ' %'
    write(*,*) ' taux d acceptation alpha 0 : ', acceptance_rate_0,' %'
    write(*,*) ' taux d acceptation alpha 1 : ', acceptance_rate_1,' %'

  end subroutine montecarlo


subroutine ajout_retrait(atconf_N, atconf_Nplus1, cel_N, cel_Nplus1, box, direc)
  
  implicit none

    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
  type(atom_config_d)::atconf_N, atconf_Nplus1 
  type(cell_config)::cel_N, cel_Nplus1
  type(box_config)::box
  integer :: direc 
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
  real(double), dimension(3,1) :: cart_vec_nplus1
  integer :: indice, i,nag

  lperiod = .true.

!######################################### Direction 0 vers 1 (ajout) ##########################################

  if (direc == 0) then ! ajout d'une particule en N+1
     if (paramcgc%rang_orig==0) then
        !tirer une position aleatoire pour le N+1eme atome
        call atom_supp(cart_vec_nplus1,boxmcgc,atconf_N)
        
        call cryst_to_cart(1,cart_vec_nplus1,boxmcgc%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m
!POURQOI CA ? Le systm        
        !copie du syst n dans n+1 
        call atconf_N%copy_config(atconf_Nplus1,lrescl=.false.)        
    !addition de la n+1eme particule
        atconf_Nplus1%xp(1:3,atconf_Nplus1%im) = cart_vec_nplus1(1:3,1)
        atconf_Nplus1%xpp(1:3,atconf_Nplus1%im) = cart_vec_nplus1(1:3,1)
        atconf_Nplus1%fp(1:3,atconf_Nplus1%im) = 0
        atconf_Nplus1%ityp(atconf_Nplus1%im) = 1
        atconf_Nplus1%ielat(atconf_Nplus1%im) = -1
        call init_vitesse(atconf_Nplus1,param = 1)
        nag=maxval(atconf_Nplus1%num_at_glob(1:atconf_Nplus1%im-1))
        atconf_Nplus1%num_at_glob(atconf_Nplus1%im) = nag+1

     end if
!repreparer les config pour le prochain langevin

        ! redecouper les cellules N et N+1
#ifdef PARA
    if(paramcgc%image==0) then !procs N
       
       call initloc(atconf_n,cel_n,atmcgcloc,cellmcgcloc,box,paramcgc,rumax,lperiod) !initloc contient caltabtc sur atloc
       call pointer_caltabt_calfo(sig,potist,atconf_n,cel_n,box,atmcgcloc,cellmcgcloc,paramcgc,&
               &lperiod,atconf_n%ltabvois,it,itetabvois,lchg=.false.)
    else !procs N+1
       call initloc(atconf_nplus1,cel_nplus1,atmcgcloc,cellmcgcloc,box,paramcgc,rumax,lperiod) !initloc contient caltabtc sur atloc
       call pointer_caltabt_calfo(sig,potist,atconf_nplus1,cel_nplus1,box,atmcgcloc,cellmcgcloc,paramcgc,&
               &lperiod,atconf_nplus1%ltabvois,it,itetabvois,lchg=.false.)
    end if
    !en ce point chacun des deux masters a les forces de son paquet datomes
    if (lmaster) then ! on est dans l'un des 2 masters
       rgcib=0;rgem=1
       if(paramcgc%image==1) then !on est dans le master de N+1
         call atconf_nplus1%send2proc(rgcib,paramcgc%comm_master,'f')
       else !on est dans le master de N qui est le master général
          call atconf_nplus1%recv(rgem,paramcgc%comm_master,'f')
       end if
    end if
!en ce point le master général (rang_orig=0) a les forces de N et N+1    
#else
    call pointer_caltabt_calfo(sig,potist,atconf_n,cel_n,box,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,atconf_n%ltabvois,it,itetabvois,lchg=.false.)    
    call pointer_caltabt_calfo(sig,potist,atconf_nplus1,cel_nplus1,box,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,atconf_nplus1%ltabvois,it,itetabvois,lchg=.false.)
#endif

!!$        
!!$    call caltabtC(cel_Nplus1,atconf_Nplus1,lperiod,boxmcgc)
!!$    call caltabtC(cel_N,atconf_N,lperiod,boxmcgc)
!!$
!!$    !calcul des forces des systemes N et N+1
!!$    CALL CalFo(sig_n,potist_n,atconf_N,cel_N,box) 
!!$    CALL CalFo(sig_nplus1,potist_nplus1,atconf_Nplus1,cel_Nplus1,box)
    if (paramcgc%rang_orig==0) then
    
    ! melange des forces des deux systemes N et N+1
       DO i=1,atconf_N%im
          atconf_Nplus1%fp(:,i) = (1-lambda_mc)*atconf_N%fp(:,i) + lambda_mc*atconf_Nplus1%fp(:,i)
       END DO
       atconf_Nplus1%fp(:,atconf_Nplus1%im) = lambda_mc*atconf_Nplus1%fp(:,atconf_Nplus1%im)

    !Egalisation des forces pour les deux systemes
       DO i=1,atconf_N%im
          atconf_N%fp(:,i) = atconf_Nplus1%fp(:,i)
       END DO
    end if
 end if

!######################################### Direction 1 vers 0 (retrait) ##########################################


 if (direc == 1) then ! retrait d'une particule alea, la placer en N+1eme position, copier le syst pour le syst à N
        if (paramcgc%rang_orig==0) then
           call indice_alea(atconf_Nplus1,indice)
!          call atconf_Nplus1%switch_atom(indice,atconf_Nplus1%im)
           call boucle_copy_atom(atconf_N,atconf_Nplus1)           
        end if
#ifdef PARA
    if(paramcgc%image==0) then !procs N
       call initloc(atconf_n,cel_n,atmcgcloc,cellmcgcloc,box,paramcgc,rumax,lperiod) !initloc contient caltabtc sur atloc
       call pointer_caltabt_calfo(sig,potist,atconf_n,cel_n,box,atmcgcloc,cellmcgcloc,paramcgc,&
            &lperiod,atconf_n%ltabvois,it,itetabvois,lchg=.false.)
    else !procs N+1
       call initloc(atconf_nplus1,cel_nplus1,atmcgcloc,cellmcgcloc,box,paramcgc,rumax,lperiod) !initloc contient caltabtc sur atloc
       call pointer_caltabt_calfo(sig,potist,atconf_nplus1,cel_nplus1,box,atmcgcloc,cellmcgcloc,paramcgc,&
            &lperiod,atconf_nplus1%ltabvois,it,itetabvois,lchg=.false.)
    end if
    !en ce point chacun des deux masters a les forces de son paquet datomes
    if (lmaster) then ! on est dans l'un des 2 masters
       rgcib=0;rgem=1
       if(paramcgc%image==1) then !on est dans le master de N+1
         call atconf_nplus1%send2proc(rgcib,paramcgc%comm_master,'f')
       else !on est dans le master de N qui est le master général
          call atconf_nplus1%recv(rgem,paramcgc%comm_master,'f')
       end if
    end if
!en ce point le master général (rang_orig=0) a les forces de N et N+1    
#else
    call pointer_caltabt_calfo(sig,potist,atconf_n,cel_n,boxmcgc,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,atconf_n%ltabvois,it,itetabvois,lchg=.false.)    
    call pointer_caltabt_calfo(sig,potist,atconf_nplus1,cel_nplus1,boxmcgc,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,atconf_nplus1%ltabvois,it,itetabvois,lchg=.false.)
#endif
!!$
!!$
!!$!repreparer les config pour le prochain langevin
!!$
!!$    !redecouper les cellules N et N+1
!!$    call caltabtC(cel_N,atconf_N,lperiod,boxmcgc)
!!$    call caltabtC(cel_Nplus1,atconf_Nplus1,lperiod,boxmcgc)
!!$
!!$    !calcul des forces des systemes N et N+1
!!$    CALL CalFo(sig_n,potist_n,atconf_N,cel_N,box) 
!!$    CALL CalFo(sig_nplus1,potist_nplus1,atconf_Nplus1,cel_Nplus1,box)

    ! melange des forces des deux systemes N et N+1
    if (paramcgc%rang_orig==0) then
       DO i=1,atconf_N%im
          atconf_Nplus1%fp(:,i) = (1-lambda_mc)*atconf_N%fp(:,i) + lambda_mc*atconf_Nplus1%fp(:,i)
       END DO
       atconf_Nplus1%fp(:,atconf_Nplus1%im) = lambda_mc*atconf_Nplus1%fp(:,atconf_Nplus1%im)
       
       !Egalisation des forces pour les deux systemes
       DO i=1,atconf_N%im
          atconf_N%fp(:,i) = atconf_Nplus1%fp(:,i)
       END DO
    end if
 end if


end subroutine ajout_retrait


subroutine init_vitesse(config, param)
  implicit none
  type(atom_config_d)::config
  real(double) :: v0, v1, z1, z2, z3, z4
  integer :: i, param ! if param =1 alors v0 défini avec tinit (juste pour l'initialisation), sinon v0 défini avec Text  

  if (param==0) then 
    v0 = sqrt(2.D0*bk*tinit)
  else 
    v0 = sqrt(2.D0*bk*Text)
  end if

  do i = 1, config%im

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




subroutine boucle_copy_atom(config_n,config_nplus1)
  implicit none
 
  type(atom_config_d)::config_n, config_nplus1
  integer :: i

  DO i=1, config_n%im
    call config_nplus1%copy_atom(i,config_n,i,lextend=.false.)
  END DO


end subroutine boucle_copy_atom


subroutine indice_alea(config, ind)

  implicit none
 
  type(atom_config_d)::config
  integer :: ind
 
  real(double) :: rand
 

  call random_number(rand)
  ind = ( (config%im - 1) * rand ) + 1

  do while (config%ityp(ind) > 1)
    call random_number(rand)
    ind = ( (config%im - 1) * rand ) + 1 
  end do 
  
  write(*,*) 'atome supprimé', ind

end subroutine indice_alea





subroutine analyse_montecarlo(atdml,celndm,box,name_file)
  
  implicit none
 
  type(atom_config_d)::atdml
  type(cell_config)::celndm
  type(box_config)::box

  character(len=*) :: name_file
  if (paramcgc%rang_orig==0) then

  
  if (itetemp>0) then
     if (mod(it,itetemp)==0) then
        call calctemp (temp,kine,atdml,celndm,latcomp=.true.)
     end if
  end if

  !if (atdml%im == 12000) then
     !name_file = 'UO2_syst_n_____'
  !else
     !name_file = 'UO2_syst_nplus1'
  !end if

  if (iterasmol>0) then
     if (mod(it,iterasmol)==0) then
        
        call rasmolT(atdml,box,it,namefr=name_file,latcomp=.true.)
     end if
  end if

!  if(paramcgc%image==0)then
!     write(6,*)name_file
!     write (6, '(I10,G10.3,A,f0.3)') it,timel,'*Temp instantanee = ',temp  
!  end if
end if
end subroutine analyse_montecarlo









subroutine calcul_U(atdml,celndm,box,potist,U_ini)

  use tempinstT_mod,only:tempinstT
  implicit none

  type(atom_config_d)::atdml
  type(cell_config)::celndm
  type(box_config)::box
  real(double) :: U_ini, potist,Ti

  Ti=tempinstT(atdml,kine,latcomp=.true.)
  
  U_ini = potist + kine
end subroutine calcul_U







subroutine atom_supp(vecteur,box,atcf)

  implicit none
  type(box_config)::box
  class(atom_config)::atcf

  real(double), dimension(3,1) :: vecteur
  real(double) :: x_nplus1, y_nplus1, z_nplus1,distati,distatM !position initiale aleatoire de la N+1eme particule
  integer::itry,i

  if (rang==0) then 

     itry=1 
1    continue

     call random_number(x_nplus1)
     call random_number(y_nplus1)
     call random_number(z_nplus1)

     vecteur(1,1) = x_nplus1
     vecteur(2,1) = y_nplus1
     vecteur(3,1) = z_nplus1
     call cryst_to_cart(1, vecteur,box%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m

     if (distminat.gt.0) then
        distatm=100
        do i=1,config_atom_n%im
           call distat(atcf%xp(:,i),vecteur, boxmcgc,distati)
           if (distatm.gt.distati) distatm=distati

        end do
        !     write(66,*)itry,distatM
        if (distatM.lt.distminat) then
           itry=itry+1
           goto 1
        end if
        write(6,*)'ITRY',itry
     end if
  end if
#ifdef PARA
  call MPI_BCAST(vecteur, 3, NDM_MPI_REAL_DOUBLE, 0,MPI_COMM_WORLD,ierr) 
#endif

end subroutine atom_supp
   





subroutine noise(var)

  implicit none
  
  real(double)  :: var(3,config_atom_nplus1%im)
  real(double)  :: u_1, u_2
  integer :: i, ic
 
  DO i=1, config_atom_nplus1%im
    do ic=1,3     
      call random_number(u_1)
      call random_number(u_2)
      var(ic,i)=sqrt(-2.*log(u_1))*cos(2.*pi*u_2)
    enddo
  ENDDO

end subroutine noise


subroutine langevin(atconf_N, atconf_Nplus1, cel_N, cel_Nplus1, box, direc)


  implicit none

  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  type(atom_config_d)::atconf_N, atconf_Nplus1 
  type(cell_config)::cel_N, cel_Nplus1
  type(box_config)::box
  integer :: direc 

  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  integer :: i,ic

  real(double) :: Ek_n, Ek_n_plus1, Ek_n_1s4, Ek_n_3s4, dQeff, Qeff,&
       &dWeff
  real(double) :: U_0, U_1, U_l_n_m1, U_l_n, H_l_n, H_l_n_m1, H_l_ini

  real(double)  :: Gl(3,atconf_N%im+1)
  real(double)::rga

  real(double), dimension(ntyp) :: aux  !pour les calculs d'acceleration
  integer::rgcib,rgem

  !initialisation des energies
  
  if (paramcgc%rang_orig==0) then ! Master général

     U_0 = 0.0
     U_1 = 0.0
     H_l_n_m1 = 0.0
     dQEff  = 0.0
     QEff   = 0.0
     dWEff  = 0.0
     WEff   = 0.0
     Ek_n =0
     DO i=1, atconf_Nplus1%im
        do ic=1,3
           Ek_n = Ek_n + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)*cm(atconf_Nplus1%ityp(i))
        enddo
     ENDDO
     call calcul_U(atconf_N,cel_N,boxmcgc,potist_n, U_0)
     call calcul_U(atconf_Nplus1,cel_Nplus1,boxmcgc,potist_nplus1, U_1)


     U_l_n = (1-lambda_mc)*U_0 + lambda_mc*U_1

     H_l_ini = (Ek_n/2) + U_l_n 
     H_l_n     = H_l_ini

     !write(6,*) 'U_1 - nrj pot et Ek_n' , (U_1-potist_nplus1), Ek_n
  end if ! Master général
  
  aux(:ntyp) = tstep/cm(:ntyp)/2.d0

  !initialisation de lambda
  if (direc == 0) then
     lambda_mc = 1.0/pas_lambda_mc
  endif

  if (direc == 1) then
     lambda_mc = 1.0 - 1.0/pas_lambda_mc
  endif

  !########################################################################################################################
  !                               Ajout d'une particule N+1: système N vers N+1 - direction = 0
  !########################################################################################################################
  if (direc == 0) then
     DO WHILE (lambda_mc <= 1)

        if (paramcgc%rang_orig==0) then ! Master général

           Ek_n = 0.0
           Ek_n_plus1 = 0.0  
           Ek_n_1s4 = 0.0
           Ek_n_3s4 = 0.0
write(6,*) 'lambda_mc' ,lambda_mc
           ! faire le pas de langevin (velocity verlet) pour determiner les nouvelles forces et positions

           ! step 1 First half-step velocities update, v(t) -> v(t+dt/2)
           timel = timel+tstep
           rga=exp((-gamlg)*tstep/2)
           call noise(Gl)
           DO i=1, atconf_Nplus1%im
              do ic=1,3
                 Ek_n = Ek_n + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)*cm(atconf_Nplus1%ityp(i))
                 atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i)*rga &
                      &+ Gl(ic,i)*sqrt(cm(atconf_Nplus1%ityp(i))*bk*Text*&
                      &(1-rga))/cm(atconf_Nplus1%ityp(i))
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
              !        atconf_N%vp(1:3,i) = atconf_Nplus1%vp(1:3,i) 
           END DO

           !conditions periodiques 
           if (lperiod)    then
              call periodbox(box,atconf_N)
              call periodbox(box,atconf_Nplus1)
           end if

        end if !on sort du master général
        !En ce point on doit transférer le système N+1 du master 0 vers le master 1

#ifdef PARA
        call MPI_BARRIER(paramcgc%comm_orig,ierr)
        if (lmaster) then ! on est dans l'un des 2 masters
           rgcib=1;rgem=0
           if(paramcgc%image==0) then !on est dans le master général
              call atconf_Nplus1%send2proc(rgcib,paramcgc%comm_master)
           else !on est dans le master de N+1
              call atconf_Nplus1%recv(rgem,paramcgc%comm_master)
           end if
        end if
        if(paramcgc%image==0) then !procs N
           call pointer_caltabt_calfo(sig,potist,atconf_N,cel_N,box,atmcgcloc,cellmcgcloc,paramcgc,&
                &lperiod,atconf_N%ltabvois,it,itetabvois)
        else !procs N+1
           call pointer_caltabt_calfo(sig,potist,atconf_Nplus1,cel_Nplus1,box,atmcgcloc,cellmcgcloc,paramcgc,&
                &lperiod,atconf_nplus1%ltabvois,it,itetabvois)
        end if
        !en ce point chacun des deux masters a les forces de son paquet datomes
        if (lmaster) then ! on est dans l'un des 2 masters
           rgcib=0;rgem=1
           if(paramcgc%image==1) then !on est dans le master de N+1
              call atconf_nplus1%send2proc(rgcib,paramcgc%comm_master,'f')
           else !on est dans le master de N qui est le master général

              call atconf_nplus1%recv(rgem,paramcgc%comm_master,'f')
           end if
        end if

#else
        
    call pointer_caltabt_calfo(sig,potist,atconf_N,cel_N,box,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,atconf_N%ltabvois,it,itetabvois)
    call pointer_caltabt_calfo(sig,potist,atconf_Nplus1,cel_Nplus1,box,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,atconf_N%ltabvois,it,itetabvois)
!Nouvelles positions
!!$        call caltabtC(cel_N,atconf_N,lperiod,boxmcgc)
!!$        call caltabtC(cel_Nplus1,atconf_Nplus1,lperiod,boxmcgc)
!!$        ! Force calculation pour chacun des systèmes avec les nouvelles positions et forces melangées
!!$        CALL CalFo(sig_n,potist_n,atconf_N,cel_N,box)
!!$        CALL CalFo(sig_nplus1,potist_nplus1,atconf_Nplus1,cel_Nplus1,box)

#endif
        if (paramcgc%rang_orig==0) then
           !mise a jour de U_l_n = (1-lambda_mc)*U_0 + lambda_mc*U_1
           U_l_n = (1-lambda_mc)*U_0 + lambda_mc*U_1
           ! melange des forces des deux systemes N et N+1
           DO i=1,atconf_N%im
              atconf_Nplus1%fp(:,i) = (1-lambda_mc)*atconf_N%fp(:,i) + lambda_mc*atconf_Nplus1%fp(:,i)
           END DO
           atconf_Nplus1%fp(:,atconf_Nplus1%im) = lambda_mc*atconf_Nplus1%fp(:,atconf_Nplus1%im)
           !Egalisation des forces pour les deux systemes
           DO i=1,atconf_N%im
              atconf_N%fp(:,i) = atconf_Nplus1%fp(:,i)
           END DO
           ! Second half-step velocities update, v(t+1/2dt) -> v(t+dt)
           call noise(Gl)
           DO i=1, atconf_Nplus1%im
              do ic=1,3
                 atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i) + aux(atconf_Nplus1%ityp(i))&
                     &*atconf_Nplus1%fp(ic,i)
                 Ek_n_3s4 = Ek_n_3s4 + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)&
                      &*cm(atconf_Nplus1%ityp(i))
                 atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i)*rga + &
                      &Gl(ic,i)*sqrt(cm(atconf_Nplus1%ityp(i))*bk*Text*&
                     &(1-rga))/cm(atconf_Nplus1%ityp(i))
                 Ek_n_plus1 = Ek_n_plus1 + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)&
                     &*cm(atconf_Nplus1%ityp(i))
              end do
           END DO
           DO i=1, atconf_N%im
              atconf_N%vp(1:3,i) = atconf_Nplus1%vp(1:3,i) 
           END DO

           !incrementation de lambda
        end if !master general

        !A FAIRE POUR TOUS LES PROCS car pilote l'arrêt     
        lambda_mc = lambda_mc + (1./pas_lambda_mc)
        if (paramcgc%rang_orig==0) then !master general
           !calcul des energies et travail et chaleur efficaces
           U_l_n_m1 = U_l_n
           H_l_n_m1 = H_l_n
           H_l_n    = Ek_n_plus1/2  + U_l_n
           dQEff  = (Ek_n_1s4-Ek_n)/2 + (Ek_n_plus1-Ek_n_3s4)/2
           QEff   = QEff + dQEff
           dWEff  = H_l_n - H_l_n_m1 - dQEff
           WEff   = WEff + dWEff
        end if

     END DO
 
  endif



  !########################################################################################################################
  !                            Déletion d'une particule N+1: système N+1 vers N - direction = 1
  !########################################################################################################################

  if (direc == 1) then
     DO WHILE (lambda_mc > 0)

        if (paramcgc%rang_orig==0) then

           Ek_n = 0.0
           Ek_n_plus1 = 0.0  
           Ek_n_1s4 = 0.0
           Ek_n_3s4 = 0.0

           ! faire le pas de langevin (velocity verlet) pour determiner les nouvelles forces et positions

           ! step 1 First half-step velocities update, v(t) -> v(t+dt/2)
           timel = timel+tstep

           write(6,*) 'lambda_mc' ,lambda_mc
           !write(6,*) 'avant langevin atconf_Nplus1%vp(:,23)=' ,atconf_Nplus1%vp(:,23)

           rga=exp((-gamlg)*tstep/2)

           call noise(Gl)

           DO i=1, atconf_Nplus1%im

              do ic=1,3

                 Ek_n = Ek_n + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)*cm(atconf_Nplus1%ityp(i))


                 atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i)*rga &
                      &+ Gl(ic,i)*sqrt(cm(atconf_Nplus1%ityp(i))*bk*Text*&
                      &(1-rga))/cm(atconf_Nplus1%ityp(i))
                 !write(6,*) 'rga=' , rga
                 !write(6,*) 'bruit=' ,Gl(ic,i)*sqrt(cm(atconf_Nplus1%ityp(i))*bk*Text*&
                 !                                                  &(1-rga))/cm(atconf_Nplus1%ityp(i))

                 !write(6,*) 'sans fp atconf_Nplus1%vp(:,1)=' ,atconf_Nplus1%vp(:,1)

                 Ek_n_1s4 = Ek_n_1s4 + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)*&
                      &cm(atconf_Nplus1%ityp(i))

                 !write(6,*) 'fp*dt=' , aux(atconf_Nplus1%ityp(i))*atconf_Nplus1%fp(ic,i)

                 atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i) + aux(atconf_Nplus1%ityp(i))&
                      &*atconf_Nplus1%fp(ic,i)

                 !write(6,*) 'atconf_Nplus1%vp(:,1)=' ,atconf_Nplus1%vp(:,1)
                 !write(6,*) 'atconf_Nplus1%fp(:,1)=' ,atconf_Nplus1%fp(:,1)

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
              !        atconf_N%vp(1:3,i) = atconf_Nplus1%vp(1:3,i) 
           END DO

           !conditions periodiques 
           !write(6,*) 'atconf_Nplus1%im =' ,atconf_Nplus1%im
           if (lperiod)    then
              call periodbox(box,atconf_N)
              call periodbox(box,atconf_Nplus1)
           end if
        end if
#ifdef PARA
        if (lmaster) then ! on est dans l'un des 2 masters
           rgcib=1;rgem=0
           if(paramcgc%image==0) then !on est dans le master général
              call atconf_Nplus1%send2proc(rgcib,paramcgc%comm_master)
           else !on est dans le master de N+1
              call atconf_Nplus1%recv(rgem,paramcgc%comm_master)
           end if
        end if
        if(paramcgc%image==0) then !procs N
           call pointer_caltabt_calfo(sig,potist,atconf_N,cel_n,box,atmcgcloc,cellmcgcloc,paramcgc,&
                &lperiod,atconf_N%ltabvois,it,itetabvois)
        else !procs N+1
           call pointer_caltabt_calfo(sig,potist,atconf_Nplus1,cel_nplus1,box,atmcgcloc,cellmcgcloc,paramcgc,&
                &lperiod,atconf_Nplus1%ltabvois,it,itetabvois)
        end if
        !en ce point chacun des deux masters a les forces de son paquet datomes
        if (lmaster) then ! on est dans l'un des 2 masters
           rgcib=0;rgem=1
           if(paramcgc%image==1) then !on est dans le master de N+1
              call atconf_Nplus1%send2proc(rgcib,paramcgc%comm_master,'f')
           else !on est dans le master de N qui est le master général
              call atconf_Nplus1%recv(rgem,paramcgc%comm_master,'f')
           end if
        end if

#else
        call pointer_caltabt_calfo(sig,potist,atconf_N,cel_n,box,atmcgcloc,cellmcgcloc,paramcgc,&
             &lperiod,atconf_N%ltabvois,it,itetabvois)
        call pointer_caltabt_calfo(sig,potist,atconf_Nplus1,cel_nplus1,box,atmcgcloc,cellmcgcloc,paramcgc,&
             &lperiod,atconf_Nplus1%ltabvois,it,itetabvois)        ! repartition des atomes des syst N et N+1 avec les nouvelles positions

        ! repartition des atomes des syst N et N+1 avec les nouvelles positions
!!$        call caltabtC(cel_N,atconf_N,lperiod,boxmcgc)
!!$        call caltabtC(cel_Nplus1,atconf_Nplus1,lperiod,boxmcgc)
!!$        ! Force calculation pour chacun des systèmes avec les nouvelles positions et forces melangées
!!$        CALL CalFo(sig_n,potist_n,atconf_N,cel_N,box)
!!$        CALL CalFo(sig_nplus1,potist_nplus1,atconf_Nplus1,cel_Nplus1,box)

#endif
        if (paramcgc%rang_orig==0) then
           !mise a jour de U_l_n = (1-lambda_mc)*U_0 + lambda_mc*U_1
           U_l_n = (1-lambda_mc)*U_0 + lambda_mc*U_1

           ! melange des forces des deux systemes N et N+1
           DO i=1,atconf_N%im
              atconf_Nplus1%fp(:,i) = (1-lambda_mc)*atconf_N%fp(:,i) + lambda_mc*atconf_Nplus1%fp(:,i)
           END DO
           atconf_Nplus1%fp(:,atconf_Nplus1%im) = lambda_mc*atconf_Nplus1%fp(:,atconf_Nplus1%im)

           !Egalisation des forces pour les deux systemes
           DO i=1,atconf_N%im
              atconf_N%fp(:,i) = atconf_Nplus1%fp(:,i)
           END DO

           !affichage temperature
           !call analyse_montecarlo(atconf_Nplus1,cel_Nplus1,box)
           it = it +1

           ! Second half-step velocities update, v(t+1/2dt) -> v(t+dt)
           call noise(Gl)

           DO i=1, atconf_Nplus1%im
              do ic=1,3
                 atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i) + aux(atconf_Nplus1%ityp(i))&
                      &*atconf_Nplus1%fp(ic,i)

                 Ek_n_3s4 = Ek_n_3s4 + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)&
                      &*cm(atconf_Nplus1%ityp(i))


                 atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i)*rga + &
                      &Gl(ic,i)*sqrt(cm(atconf_Nplus1%ityp(i))*bk*Text*&
                      &(1-rga))/cm(atconf_Nplus1%ityp(i))


                 Ek_n_plus1 = Ek_n_plus1 + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)&
                      &*cm(atconf_Nplus1%ityp(i))

              end do
           END DO

           !write(6,*) 'apres langevin atconf_Nplus1%vp(:,23)=' ,atconf_Nplus1%vp(:,23)
           !write(6,*)

           !repartir les nouvelles positions et forces dans les syst N et N+1
           !Systeme a N

           DO i=1, atconf_N%im
              !atconf_N%xpp(1:3,i) = atconf_Nplus1%xpp(1:3,i)
              !atconf_N%xp(1:3,i) = atconf_Nplus1%xp(1:3,i)
              atconf_N%vp(1:3,i) = atconf_Nplus1%vp(1:3,i) 
           END DO

           !incrementation de lambda
        end if
        lambda_mc = lambda_mc - (1./pas_lambda_mc)
        if (paramcgc%rang_orig==0) then
           !calcul des energies et travail et chaleur efficaces
           U_l_n_m1 = U_l_n
           H_l_n_m1 = H_l_n
           H_l_n    = Ek_n_plus1/2  + U_l_n

           dQEff  = (Ek_n_1s4-Ek_n)/2 + (Ek_n_plus1-Ek_n_3s4)/2
           QEff   = QEff + dQEff
           dWEff  = H_l_n - H_l_n_m1 - dQEff
           WEff   = WEff + dWEff
           !write(6,*) 'dWEff WEff', dWEff, WEff
           !write(6,*) 'U_l_n_m1, U_l_n  H_l_n_m1, H_l_n', U_l_n_m1, U_l_n,  H_l_n_m1, H_l_n
           !write(6,*)
           !write(6,*) 'dQEff QEff dWEff WEff', dQEff, QEff, dWEff, WEff


        end if
     END DO


  endif


end subroutine langevin

  subroutine init_mpi_MCGC
    
!  subroutine init_mpi_neb()
    ! Routine d'initialisation de MPI pour la NEB
#ifdef PARA

    paramcgc%np_orig=nprocs
    paramcgc%rang_orig=rang
!    paramcgc%grp_orig=grp_world
    paramcgc%nimage=2
!    paramcgc%comm_orig=MPI_COMM_WORLD
    call MPI_COMM_DUP(MPI_COMM_WORLD,paramcgc%comm_orig,ierr)
    call MPI_COMM_GROUP(paramcgc%comm_orig,paramcgc%grp_orig,ierr)
    call commconstr(paramcgc)

    myidsp=paramcgc%rgim
    MPI_COMM_space=paramcgc%comm_image
    nprocspace=paramcgc%npim
#else
    paramcgc%np_orig=1
    paramcgc%rang_orig=0
    paramcgc%npim=1
    myidsp=0
    paramcgc%lmaster=.true.
    nprocspace=1
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

    
  
  subroutine initNP1

    real(double), dimension(3,1) :: cart_vec_nplus1
    real(double)::distati
    integer::i
    !definir le systeme a N+1 en tirant une position aleatoire pour le N+1eme atome
    call atom_supp(cart_vec_nplus1,boxmcgc,config_atom_n)

!FAIT DANS atom_supp
!    call cryst_to_cart(1,cart_vec_nplus1,boxmcgc%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m

    !copie du syst n dans n+1 
    call config_atom_nplus1%init(config_atom_n%im+1,config_atom_n%imm,config_atom_n%ltabvois)
    config_atom_nplus1%ltabvois=config_atom_n%ltabvois
    call config_atom_n%copy_config(config_atom_nplus1,lrescl=.false.)


    !addition de la n+1eme particule
    config_atom_nplus1%xp(1:3,config_atom_nplus1%im) = cart_vec_nplus1(1:3,1)
    config_atom_nplus1%fp(1:3,config_atom_nplus1%im) = 0
    config_atom_nplus1%xpp(1:3,config_atom_nplus1%im) =     config_atom_nplus1%xp(1:3,config_atom_nplus1%im) 
    config_atom_nplus1%ityp(config_atom_nplus1%im) = 1
    config_atom_nplus1%num_at_glob(config_atom_nplus1%im) = config_atom_nplus1%im
    call init_vitesse(config_atom_nplus1,param = 0)
    !copie de cell puis caltabtC pour redecouper avec la n+1eme particule
    call cells_nplus1%init(cells_n%nox,cells_n%noy,cells_n%noz, cells_n%natperc)
    call cells_n%copy_cell(cells_nplus1)

    call caltabtC(cells_nplus1,config_atom_nplus1,lperiod,boxmcgc)

#ifdef PARA
    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
       if(paramcgc%image==0) then !procs N
          call init_voisinage(cells_n)
          !?          call maj_atomes_frt_ftm(config_atom_n,cells_n)

       else !procs N+1
          call init_voisinage(cells_nplus1)
          !?          call maj_atomes_frt_ftm(config_atom_nplus1,cells_nplus1)
       end if

    end if
    !       call MPI_BARRIER(MPI_COMM_WORLD)
#endif


    if((ipotentiel==-10).or.(ipotentiel==-11)) then
       if (paramcgc%rang_orig==0) then
          write(6,*)'write configuration N+1  to confNP1.lmp'
          call config2data (config_atom_nplus1%imm,config_atom_nplus1%im,&
               config_atom_nplus1%xp,config_atom_nplus1%ityp,boxmcgc%at,ntyp,filename='conf.lmp.NP1')
       end if

#ifdef PARA
#ifdef LAMMPS_VERSION
    if(paramcgc%image==0) then !procs N
       firsttime_lammps=.true.
       allocate (posa(3*config_atom_n%im),  forca(3*config_atom_n%im))

       call init_lammps('in.lammps.N')
    else !procs N+1
       firsttime_lammps=.true.
       allocate (posa(3*config_atom_nplus1%im),  forca(3*config_atom_nplus1%im))
       call init_lammps('in.lammps.NP1')

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
  end subroutine initNP1

end module montecarlo_mod
