module montecarlo_mod
  USE gen_com_m,only:  lperiod, tstep, timel, tstep, sig, potist,itetabvois,&
       & it,iterasmol,itetemp, temp, kine, pi, bk, Text, gamlg,&
       &lspaceNDM,rang
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
  USE calctemp_mod,only: calctemp
  USE rasmolT_mod,only: rasmolT
  use paraconfig,only:para_config,commconstr
#ifdef PARA
  use mod_para,only:grp_world,nprocs,myidsp,MPI_COMM_space,nprocspace,ierr,mpi_comm_world
  USE init_vois_mod,only: init_voisinage
#else
  use mod_para,only:myidsp,nprocspace
#endif
  use read_val,only:rvois,ltabvois
  use var_pot,only:ipotentiel,rumax
  USE parautils,only:initloc,pointer_caltabt_calfo


  implicit none

  type(para_config)::paramcgc


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

  subroutine init_mpi_MCGC
    
!  subroutine init_mpi_neb()
    ! Routine d'initialisation de MPI pour la NEB
#ifdef PARA


    paramcgc%np_orig=nprocs
    paramcgc%rang_orig=rang
    paramcgc%grp_orig=grp_world

    paramcgc%nimage=2

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
  
subroutine montecarlo


    implicit none

    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    type(box_config)::boxndm
    
    type(atom_config_d):: config_atom_old_0, config_atom_old_1, config_atom_new_0, config_atom_new_1 !config intermediaire pour suivre l'evolution des systemes: 0 -> syst N, 1 -> syst N+1.

    real(double), dimension(3,1) :: cart_vec_nplus1
      
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


!defintion de la boite du syst a N atomes
!    call ndm2boxconfig(at,bg,zl,zls2,nzl,volu,normat,boxndm)
!    call caltabtC(cells_n,config_atom_n,lperiod,bg)
    !call cells_n%print

!definir le systeme a N+1 en tirant une position aleatoire pour le N+1eme atome
    call atom_supp(cart_vec_nplus1)

    call cryst_to_cart(1,cart_vec_nplus1,boxmcgc%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m

!copie du syst n dans n+1 
    call config_atom_nplus1%init(config_atom_n%im+1,config_atom_n%imm,config_atom_n%ltabvois)
    config_atom_nplus1%ltabvois=config_atom_n%ltabvois
    call config_atom_n%copy_config(config_atom_nplus1,lrescl=.false.)

    !write(6,*) 'config_atom_n%vp(:,56) =' ,config_atom_n%vp(:,56)
    !write(6,*) 'config_atom_nplus1%vp(:,56) =' ,config_atom_nplus1%vp(:,56)

!addition de la n+1eme particule
    config_atom_nplus1%xp(1:3,config_atom_nplus1%im) = cart_vec_nplus1(1:3,1)
    config_atom_nplus1%fp(1:3,config_atom_nplus1%im) = 0
    config_atom_nplus1%vp(1:3,config_atom_nplus1%im) = 0
    config_atom_nplus1%xpp(1:3,config_atom_nplus1%im) = 0
    config_atom_nplus1%ityp(config_atom_nplus1%im) = 1

!copie de cell puis caltabtC pour redecouper avec la n+1eme particule
    call cells_nplus1%init(cells_n%nox,cells_n%noy,cells_n%noz, cells_n%natperc)
    call cells_n%copy_cell(cells_nplus1)

!CRC
!syst N et N+1 definis et non repartis
!    call caltabtC(cells_nplus1,config_atom_nplus1,lperiod,bg)

!lorsque ltabvois = true, attention, il faut la recalculer pour le syst n+1
    !if (config_atom_nplus1%ltabvois)&
    !&call caltabi(config_atom_nplus1%atom_config,cells_nplus1)

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
!    WRITE(6,*)'PN'
    call pointer_caltabt_calfo(sig,potist,config_atom_n,cells_n,boxmcgc,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,config_atom_n%ltabvois,it,itetabvois,lchg=.false.)
!        WRITE(6,*)'PN+1'
    call pointer_caltabt_calfo(sig,potist,config_atom_nplus1,cells_nplus1,boxmcgc,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,config_atom_nplus1%ltabvois,it,itetabvois,lchg=.false.)
!    WRITE(6,*)'PN+2'
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
!    WRITE(6,*)'PL'
    
       call langevin(config_atom_n,config_atom_nplus1,cells_n,cells_nplus1,boxmcgc,direction)

       if (paramcgc%rang_orig==0) then
          call config_atom_nplus1%copy_config(config_atom_old_1, lrescl=.true.)
          
          W = WEff
          xprob = 1
          Wprec = + W
          
          direction = 1
          
          seed(1) = 152533
          call random_seed(PUT=seed(1:12))
       end if
          
!########################################################################################################################
!                                             boucle sur lambda le long d'un chemin
!########################################################################################################################

          DO i_path = 1, n_path 
             if (paramcgc%rang_orig==0) then
             
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
          end if
          call ajout_retrait(config_atom_n,config_atom_nplus1,cells_n,cells_nplus1,boxmcgc,direction)
          if (paramcgc%rang_orig==0) then
             call analyse_montecarlo(config_atom_n,cells_n,boxmcgc, 'UO2_syst_n_before_test')
             call analyse_montecarlo(config_atom_nplus1,cells_nplus1,boxmcgc, 'UO2_syst_nplus1_before_test')
          end if
             ! pas de langevin
!          WRITE(6,*)'P4N'
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
             
             if (direction == 0) then
                direction = 1
             else
                direction = 0
             end if
             
             
             acceptance_rate   = (real(n_accepted)/real(n_gen))*1.0d2
             acceptance_rate_0 = (real(n_accepted_0)/real(n_gen_0))*1.0d2
             acceptance_rate_1 = (real(n_accepted_1)/real(n_gen_1))*1.0d2
             
             call analyse_montecarlo(config_atom_n,cells_n,boxmcgc, 'UO2_syst_n_after_test')
             call analyse_montecarlo(config_atom_nplus1,cells_nplus1,boxmcgc, 'UO2_syst_nplus1_after_test')
             it = it +1
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
  integer :: indice, i

  lperiod = .true.

!######################################### Direction 0 vers 1 (ajout) ##########################################

  if (direc == 0) then ! ajout d'une particule en N+1
     if (paramcgc%rang_orig==0) then
        !tirer une position aleatoire pour le N+1eme atome
        call atom_supp(cart_vec_nplus1)
        
        call cryst_to_cart(1,cart_vec_nplus1,boxmcgc%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m
        
        !copie du syst n dans n+1 
        call copy_nplus1_xp(atconf_N,atconf_Nplus1)
        
    !addition de la n+1eme particule
        atconf_Nplus1%xp(1:3,atconf_Nplus1%im) = cart_vec_nplus1(1:3,1)


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
           call atconf_Nplus1%switch_atom(indice,atconf_Nplus1%im)
           call copy_n_xp(atconf_N,atconf_Nplus1)
           
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




subroutine copy_nplus1_xp(config_n,config_nplus1)

  implicit none
 
  type(atom_config_d)::config_n, config_nplus1
  integer :: i

  DO i=1, config_n%im
    config_nplus1%xp(:,i) = config_n%xp(:,i)
  END DO


end subroutine copy_nplus1_xp






subroutine copy_n_xp(config_n,config_nplus1)

  implicit none
 
  type(atom_config_d)::config_n, config_nplus1
  integer :: i

  DO i=1, config_n%im
    config_n%xp(:,i) = config_nplus1%xp(:,i)
  END DO

end subroutine copy_n_xp




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

  if (itetemp>0) then
     if (mod(it,itetemp)==0) then
        call calctemp (temp,kine,atdml,celndm)
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

  !write (6, '(I10,G10.3,A,f0.3)') it,timel,'*Temp instantanee = ',temp  


end subroutine analyse_montecarlo









subroutine calcul_U(atdml,celndm,box,potist,U_ini)

  implicit none

  type(atom_config_d)::atdml
  type(cell_config)::celndm
  type(box_config)::box
  real(double) :: U_ini, potist

  call calctemp (temp,kine,atdml,celndm)
  
  U_ini = potist + kine
  !write(6,*) 'kine' , kine
end subroutine calcul_U







subroutine atom_supp(vecteur)
  
  implicit none

  real(double), dimension(3,1) :: vecteur
  real(double) :: x_nplus1, y_nplus1, z_nplus1 !position initiale aleatoire de la N+1eme particule

  call random_number(x_nplus1)
  call random_number(y_nplus1)
  call random_number(z_nplus1)

  vecteur(1,1) = x_nplus1
  vecteur(2,1) = y_nplus1
  vecteur(3,1) = z_nplus1

  
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
!  WRITE(6,*)'IL',paramcgc%rang_orig
  
  if (paramcgc%rang_orig==0) then

     U_0 = 0.0
     U_1 = 0.0
     H_l_n_m1 = 0.0
     dQEff  = 0.0
     QEff   = 0.0
     dWEff  = 0.0
     WEff   = 0.0

     DO i=1, atconf_Nplus1%im
        do ic=1,3
           Ek_n = Ek_n + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)*cm(atconf_Nplus1%ityp(i))
        enddo
     ENDDO

     call calcul_U(atconf_N,cel_N,boxmcgc,potist_n, U_0)
     call calcul_U(atconf_Nplus1,cel_Nplus1,boxmcgc,potist_nplus1, U_1)
     !write(6,*) 'Ek_n' , Ek_n

     U_l_n = (1-lambda_mc)*U_0 + lambda_mc*U_1

     H_l_ini = (Ek_n/2) + U_l_n 
     H_l_n     = H_l_ini

     !write(6,*) 'U_1 - nrj pot et Ek_n' , (U_1-potist_nplus1), Ek_n


  end if
  aux(:ntyp) = tstep/cm(:ntyp)/2.d0
  !write(6,*) 'aux(:ntyp)' ,aux(:ntyp)

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
!  WRITE(6,*)'direc',direc
  if (direc == 0) then
     DO WHILE (lambda_mc <= 1)

        if (paramcgc%rang_orig==0) then

           Ek_n = 0.0
           Ek_n_plus1 = 0.0  
           Ek_n_1s4 = 0.0
           Ek_n_3s4 = 0.0

           ! faire le pas de langevin (velocity verlet) pour determiner les nouvelles forces et positions

           ! step 1 First half-step velocities update, v(t) -> v(t+dt/2)
           timel = timel+tstep

           !write(6,*) 'lambda_mc' ,lambda_mc
           !write(6,*)
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

        end if !on sort du master général
        !En ce point on doit transférer le système N+1 du master 0 vers le master 1

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
!        WRITE(6,*)'P2N'
        
    call pointer_caltabt_calfo(sig,potist,atconf_N,cel_N,box,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,atconf_N%ltabvois,it,itetabvois)
!    WRITE(6,*)'P2N1'
    call pointer_caltabt_calfo(sig,potist,atconf_Nplus1,cel_Nplus1,box,atmcgcloc,cellmcgcloc,paramcgc,&
         &lperiod,atconf_N%ltabvois,it,itetabvois)
!Nouvelles positions
!           WRITE(6,*)'P2N2'
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
        !A FAIRE POUR TOUS LES PROCS      
        lambda_mc = lambda_mc + (1./pas_lambda_mc)
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
           !write(6,*) 'U_l_n_m1, U_l_n  H_l_n_m1, H_l_n dQEff QEff dWEff WEff', U_l_n_m1, U_l_n,  H_l_n_m1,& 
           !&H_l_n, dQEff, QEff, dWEff, WEff
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

           !write(6,*) 'lambda_mc' ,lambda_mc
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
!        WRITE(6,*)'P3N'
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
           !it = it +1

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





end module montecarlo_mod
