module montecarlo_mod
  USE gen_com_m,only:  lperiod, tstep, timel, tstep, sig, itetabvois,&
       & iterasmol,itetemp, temp, kine, pi, bk, Text, gamlg,one,pi,text,tinit,&
       &lspaceNDM,rang,it,firsttime_lammps,posa,forca,erg2ev,parallele
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
  use Tpara,only:grp_world,nprocs,myidsp,MPI_COMM_space,nprocspace,ierr,mpi_comm_world,&
       &NDM_MPI_REAL_DOUBLE,para_space_config,status,comm_space
  use mod_para,only:maj_atomes_frt_ftm
  USE init_vois_mod,only: init_voisinage
#else
  use Tpara,only:myidsp,nprocspace,para_space_config
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
  type(para_space_config)::pscgc
 real(double)::distminat
  type(atom_config_d)::atconf_n !type derive atom_config du systeme a n atomes
  type(atom_config_d)::atconf_nplus1 !type derive atom_config du systeme a n+1 atomes

  type(cell_config):: cells_n !type derive cell_config du systeme a n atomes
  type(cell_config):: cells_nplus1 !type derive cell_config du systeme a n+1 atomes

  type(box_config)::boxmcgc
  
  integer::  pas_lambda_mc
  real(double) :: lambda_mc !lambda compris entre 0 et 1

  integer :: n_path ! nb de chemin d'insertion, a definir dans .din, par defaut 10

  real(double), dimension(3,3) :: sig_n, sig_nplus1
  real(double) :: potist_n, potist_nplus1
  real(double) :: Weff
  logical::lmaster !(true= master du système N ou du système N+1)
    class(atom_config),pointer::atmcgcloc
    type(cell_config),pointer::cellmcgcloc
    type(cell_config),target:: cellcible ! ne sert qu'à faire pointer cellnebloc sur quelquechose
    type(atom_config_d),target::atcible
    integer::rgcib,rgem !(cible et emeteur e confN+1)
    integer, dimension(12) :: seed 

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

    integer :: i_path

    real(double) :: theta
    real(double) :: beta

    real(double) :: W, Wprec, xprob, xalea
    real(double) :: ln_xalea, ln_Wprec, ln_W, ln_xprob
    real(double) :: acceptance_rate, acceptance_rate_0, acceptance_rate_1


    integer :: n,iloc

    logical :: lchange,ldistrib

    !########################################################################################################################
    !                                             Initialisation
    !########################################################################################################################

    !initialisation variables 
    !pour le premier chemin: sens positif, d'ajout d'une particule et acceptation     
#ifdef PARA
    cellmcgcloc=>cellcible
    atmcgcloc=>atcible
    lmaster=paramcgc%lmaster
    if ((paramcgc%mpi_image%nproc.gt.1).and.(lspaceNDM.eqv..true.)) then
       call init_voisinage(cells_n,pscgc)
    end if
#else
    lmaster=.true.
#endif    


    direction = 0 ! direction = 0 on ajoute un atome, = 1 on retire un atome
    lambda_mc = 0.0

    lextend = .true.
    lperiod = .true.

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


    iloc=1;lchange=.false.;ldistrib=.true.
    call calfoMCGC(iloc,lchange,ldistrib)

    !initialisation pour le dyn_vverlet
    !calcul des forces des systemes N et N+1
    ! melange des forces des deux systemes N et N+1

    if (paramcgc%mpi_orig%rank==0) then
       ! sauvegarde du système
       call atconf_n%copy_config(config_atom_old_0, lrescl=.true.)
    end if

    !pour le premier chemin: sens positif, d'ajout d'une particule et acceptation

    call langevin(direction)
    if (paramcgc%mpi_orig%rank==0) then
       call atconf_nplus1%copy_config(config_atom_old_1, lrescl=.true.)

       W = WEff
       xprob = 1
       Wprec = + W

       write(*,*) 'W0', W,W*erg2eV

       seed(1) = 152533
       call random_seed(PUT=seed(1:12))
    end if
    direction = 1
    !########################################################################################################################
    !                                             boucle sur lambda le long d'un chemin
    !########################################################################################################################
    
    DO i_path = 1, n_path ! boucle à faire pour tous les procs
       if (direction == 0) then
          lambda_mc = 0.0
       endif
       if (direction == 1) then
          lambda_mc = 1.0 
       endif

       if (paramcgc%mpi_orig%rank==0) then !!master general

          write(*,*) ' '
          write(*,*) ' '
          write(*,*) ' '
          write(*,*) 'Numéro de chemin', i_path,direction

          call random_number(xalea)
          ln_xalea  = log(xalea)
          atconf_nplus1%vp(:,:)   = - atconf_nplus1%vp(:,:) !à chaque retour dans la boucle, on change de direction
          !initialisation de lambda
          if (direction == 0) then
             lambda_mc = 0.0
             call atconf_n%copy_config(config_atom_new_0, lrescl=.true.)
          endif
          if (direction == 1) then
             lambda_mc = 1.0 
             call atconf_nplus1%copy_config(config_atom_new_1, lrescl=.true.)
          endif

       end if !master general

          !choisir l'at a retirer ou ajouter + preparation des syst N et N+1 pour etre prets pour le langevin (cad decoupage cellules + calcul forces + melange des forces - se fait dans cette sous routine)

       call ajout_retrait(direction)
if (paramcgc%mpi_orig%rank==0) then
       call analyse_montecarlo(atconf_n,cells_n,boxmcgc, 'UO2_syst_n_before_test')
       call caltabtC(cells_nplus1,atconf_nplus1,lperiod,boxmcgc)
       call analyse_montecarlo(atconf_nplus1,cells_nplus1,boxmcgc, 'UO2_syst_nplus1_before_test')
       ! pas de langevin
    end if
       call langevin(direction)
       if (paramcgc%mpi_orig%rank==0) then
          if (direction == 0) then
             W = +WEff
             n_gen_0 = n_gen_0 + 1
             call atconf_nplus1%copy_config(config_atom_new_1, lrescl=.true.)
          else
             W = -WEff
             n_gen_1 = n_gen_1 + 1
             call atconf_n%copy_config(config_atom_new_0, lrescl=.true.)
          endif
          n_gen = n_gen + 1
          ln_Wprec  = (+beta*(direction-theta)*Wprec)
          ln_W      = (+beta*(direction-theta)*W)
          ln_xprob  = - dlog(1 + dexp(ln_Wprec-ln_W))
          xprob     = dexp(ln_xprob)

          write(6,*) 'Wprec', Wprec, 'ln_Wprec', ln_Wprec,'Wprec eV',Wprec*erg2eV
          write(6,*) 'W', W, 'ln_W', ln_W,'W eV', W*erg2eV

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
                call config_atom_old_1%copy_config(atconf_nplus1, lrescl=.true.)
                call caltabtC(cells_nplus1,atconf_nplus1,lperiod,boxmcgc)
             end if

             if (direction == 1) then
                call config_atom_old_0%copy_config(atconf_n, lrescl=.true.)
                call caltabtC(cells_n,atconf_n,lperiod,boxmcgc)
             end if


          end if

          acceptance_rate   = (real(n_accepted)/real(n_gen))*1.0d2
          acceptance_rate_0 = (real(n_accepted_0)/real(n_gen_0))*1.0d2
          acceptance_rate_1 = (real(n_accepted_1)/real(n_gen_1))*1.0d2



          call analyse_montecarlo(atconf_n,cells_n,boxmcgc, 'UO2_syst_n_after_test')
          call caltabtC(cells_nplus1,atconf_nplus1,lperiod,boxmcgc)
          call analyse_montecarlo(atconf_nplus1,cells_nplus1,boxmcgc, 'UO2_syst_nplus1_after_test')
       end if
       it = it +1
       
       if (direction == 0) then
          direction = 1
       else
          direction = 0
       end if
    END DO

    if (paramcgc%mpi_orig%rank==0) then
    write(*,*) ' taux d acceptation final   : ', acceptance_rate,  ' %'
    write(*,*) ' taux d acceptation alpha 0 : ', acceptance_rate_0,' %'
    write(*,*) ' taux d acceptation alpha 1 : ', acceptance_rate_1,' %'
    end if

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
     if (paramcgc%mpi_orig%rank==0) then
        !tirer une position aleatoire pour le N+1eme atome
        call atom_supp(cart_vec_nplus1)
        
        call cryst_to_cart(1,cart_vec_nplus1,boxmcgc%at,1) !at vecteur de base de la boite en cm, defini dans gen_com_m
!POURQOI CA ? Le systm        
        !copie du syst n dans n+1 
        call boucle_copy_atom(atconf_N,atconf_Nplus1, sens= .false.)        
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
        if (paramcgc%mpi_orig%rank==0) then
           call indice_alea(atconf_Nplus1,indice)
           call atconf_Nplus1%switch_atom(indice,atconf_Nplus1%im)

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


subroutine init_vitesse(config, param) !juste pour le N+1eme atome
  implicit none
  type(atom_config_d)::config
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
 
  type(atom_config_d)::config_n, config_nplus1
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
 
  type(atom_config_d)::config
  integer :: ind
 
  real(double) :: rand
 

  call random_number(rand)
  ind = ( (config%im - 1) * rand ) + 1

  do while (config%ityp(ind) > 1)
    call random_number(rand)
    ind = ( (config%im - 1) * rand ) + 1 
  end do 

end subroutine indice_alea





subroutine analyse_montecarlo(atdml,celndm,box,name_file)
  
  implicit none
 
  type(atom_config_d)::atdml
  type(cell_config)::celndm
  type(box_config)::box

  character(len=*) :: name_file
  if (paramcgc%mpi_orig%rank==0) then

  
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

  type(atom_config_d)::atdml
  type(cell_config)::celndm
  type(box_config)::box
  real(double) :: U_ini, potist,Ti

  if (paramcgc%mpi_orig%rank==0) then
  Ti=tempinstT(atdml,kine,latcomp=.true.)
  
  U_ini = potist + kine
  end if
end subroutine calcul_U




subroutine atom_supp(vecteur)
  
  implicit none

  real(double), dimension(3,1) :: vecteur
  real(double) :: x_nplus1, y_nplus1, z_nplus1 !position initiale aleatoire de la N+1eme particule
  seed(1) = 152533
  call random_seed(PUT=seed(1:12))

  call random_number(x_nplus1)
  call random_number(y_nplus1)
  call random_number(z_nplus1)

  vecteur(1,1) = x_nplus1
  vecteur(2,1) = y_nplus1
  vecteur(3,1) = z_nplus1

  
end subroutine atom_supp

subroutine noise(var)

  implicit none
  
  real(double)  :: var(3,atconf_nplus1%im)
  real(double)  :: u_1, u_2
  integer :: i, ic
 
  DO i=1, atconf_nplus1%im
    do ic=1,3     
      call random_number(u_1)
      call random_number(u_2)
      var(ic,i)=sqrt(-2.*log(u_1))*cos(2.*pi*u_2)
    enddo
  ENDDO

end subroutine noise


subroutine langevin( direc)
  implicit none
  integer :: direc 
  integer :: i,ic
  real(double) :: Ek_n, Ek_n_plus1, Ek_n_1s4, Ek_n_3s4, dQeff, Qeff,&
       &dWeff
  real(double) :: U_0, U_1, U_l_n_m1, U_l_n, H_l_n, H_l_n_m1, H_l_ini

  real(double)  :: Gl(3,atconf_N%im+1)
  real(double)::rga

  real(double), dimension(ntyp) :: aux  !pour les calculs d'acceleration
  integer::rgcib,rgem,iloc
  logical::lchange,ldistrib

  !initialisation des energies
  if (paramcgc%mpi_orig%rank==0) then ! Master général

     U_0 = 0.0
     U_1 = 0.0
     H_l_n_m1 = 0.0
     dQEff  = 0.0
     QEff   = 0.0
     dWEff  = 0.0
     WEff   = 0.0
     Ek_n = 0.0

     DO i=1, atconf_Nplus1%im
        do ic=1,3
           Ek_n = Ek_n + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)*cm(atconf_Nplus1%ityp(i))
        enddo
     ENDDO

     U_l_n = (1-lambda_mc)*potist_n + lambda_mc*potist_nplus1

     H_l_ini = Ek_n + U_l_n 
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

        if (paramcgc%mpi_orig%rank==0) then ! Master général

           Ek_n = 0.0
           Ek_n_plus1 = 0.0  
           Ek_n_1s4 = 0.0
           Ek_n_3s4 = 0.0
           write(6,*) 'lambda_mc D0' ,rang,lambda_mc
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


    if (paramcgc%mpi_orig%rank==0) then
           !mise a jour de U_l_n = (1-lambda_mc)*U_0 + lambda_mc*U_1
           U_l_n = (1-lambda_mc)*potist_n + lambda_mc*potist_nplus1
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

        !incrementation de lambda    
        lambda_mc = lambda_mc + (1./pas_lambda_mc)
        if (paramcgc%mpi_orig%rank==0) then !master general
           !calcul des energies et travail et chaleur efficaces
           U_l_n_m1 = U_l_n
           H_l_n_m1 = H_l_n
           H_l_n    = Ek_n_plus1  + U_l_n
           dQEff  = (Ek_n_1s4-Ek_n) + (Ek_n_plus1-Ek_n_3s4)
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

        if (paramcgc%mpi_orig%rank==0) then
           Ek_n = 0.0
           Ek_n_plus1 = 0.0  
           Ek_n_1s4 = 0.0
           Ek_n_3s4 = 0.0
           ! faire le pas de langevin (velocity verlet) pour determiner les nouvelles forces et positions
           ! step 1 First half-step velocities update, v(t) -> v(t+dt/2)
           timel = timel+tstep
           write(6,*) 'lambda_mc D1' ,rang,lambda_mc
           !write(6,*) 'avant langevin atconf_Nplus1%vp(:,23)=' ,atconf_Nplus1%vp(:,23)
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

        if (paramcgc%mpi_orig%rank==0) then

           !mise a jour de U_l_n = (1-lambda_mc)*U_0 + lambda_mc*U_1
           U_l_n = (1-lambda_mc)*potist_n + lambda_mc*potist_nplus1

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


                 atconf_Nplus1%vp(ic,i) = atconf_Nplus1%vp(ic,i)*rga + &
                      &Gl(ic,i)*sqrt(cm(atconf_Nplus1%ityp(i))*bk*Text*&
                      &(1-rga))/cm(atconf_Nplus1%ityp(i))


                 Ek_n_plus1 = Ek_n_plus1 + 0.5*atconf_Nplus1%vp(ic,i)*atconf_Nplus1%vp(ic,i)&
                      &*cm(atconf_Nplus1%ityp(i))

              end do
           END DO

           !repartir les nouvelles positions et forces dans les syst N et N+1
           !Systeme a N

           DO i=1, atconf_N%im
              atconf_N%vp(1:3,i) = atconf_Nplus1%vp(1:3,i) 
           END DO

           !incrementation de lambda
        end if
        lambda_mc = lambda_mc - (1./pas_lambda_mc)
        if (paramcgc%mpi_orig%rank==0) then
           !calcul des energies et travail et chaleur efficaces
           U_l_n_m1 = U_l_n
           H_l_n_m1 = H_l_n
           H_l_n    = Ek_n_plus1  + U_l_n
           dQEff  = (Ek_n_1s4-Ek_n) + (Ek_n_plus1-Ek_n_3s4)
           QEff   = QEff + dQEff
           dWEff  = H_l_n - H_l_n_m1 - dQEff
           WEff   = WEff + dWEff

        end if
     END DO


  endif


end subroutine langevin

  subroutine init_mpi_MCGC
    
!  subroutine init_mpi_neb()
    ! Routine d'initialisation de MPI pour la NEB
#ifdef PARA

    paramcgc%mpi_orig%nproc=nprocs
    paramcgc%mpi_orig%rank=rang
!    paramcgc%mpi_orig%group=grp_world
    paramcgc%nimage=2
!    paramcgc%mpi_orig%comm=MPI_COMM_WORLD
    call MPI_COMM_DUP(MPI_COMM_WORLD,paramcgc%mpi_orig%comm,ierr)
    call MPI_COMM_GROUP(paramcgc%mpi_orig%comm,paramcgc%mpi_orig%group,ierr)
    call commconstr(paramcgc)

    myidsp=paramcgc%mpi_image%rank
    call MPI_COMM_free(mpi_comm_space,ierr)
    MPI_COMM_space=paramcgc%mpi_image%comm
    call comm_space%init(MPI_COMM_SPACE)
    nprocspace=paramcgc%mpi_image%nproc
    if (nprocspace==1) parallele=.false.
#else
    paramcgc%mpi_orig%nproc=1
    paramcgc%mpi_orig%rank=0
    paramcgc%mpi_image%nproc=1
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
    lmaster=paramcgc%lmaster
    if (paramcgc%mpi_orig%rank==0) then
       
       call atom_supp(cart_vec_nplus1)
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


    if((ipotentiel==-10).or.(ipotentiel==-11)) then
       if (paramcgc%mpi_orig%rank==0) then
          write(6,*)'write configuration N+1  to confNP1.lmp'
          call config2data (atconf_nplus1%imm,atconf_nplus1%im,&
               atconf_nplus1%xp,atconf_nplus1%ityp,boxmcgc%at,ntyp,filename='confNP1.lmp')
       end if

#ifdef PARA
#ifdef LAMMPS_VERSION
    if(paramcgc%image==0) then !procs N
       firsttime_lammps=.true.
       allocate (posa(3*atconf_n%im),  forca(3*atconf_n%im))

       call init_lammps('in.lammps.N')
    else !procs N+1
       firsttime_lammps=.true.
       allocate (posa(3*atconf_nplus1%im),  forca(3*atconf_nplus1%im))
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

    if (paramcgc%mpi_orig%rank==0) then
       DO i=1,atconf_n%im
          atconf_nplus1%fp(:,i) = (1-lambda_mc)*atconf_n%fp(:,i) + lambda_mc*atconf_nplus1%fp(:,i)
       END DO
       atconf_nplus1%fp(:,atconf_nplus1%im) = lambda_mc*atconf_nplus1%fp(:,atconf_nplus1%im)

       !Egalisation des forces pour les deux systemes
       DO i=1,atconf_n%im
          atconf_n%fp(:,i) = atconf_nplus1%fp(:,i)
       END DO
    end if

  end subroutine calfoMCGC
    
  
end module montecarlo_mod
