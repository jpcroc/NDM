module babar_mod
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m,only: rang,fnam,lenfnam,lmultin,imm_glob,dmtype,itloopmax,itmax,timemax,timeloopmax,iteration,lwrtb,unitwb,&
       &fnam,lenfnam,lmasterb,lspacendm,latcomp,ivisu,igen,text,bk,timel,tstep,potist
  USE atomconfig,only:atom_config_e,atom_config_d
  USE cellconfig, only:cell_config, caltabtC
  USE arret_ndm_mod,only:arret_ndm
  USE boxconfig,only:box_config,periodbox,box_config_lpr,updatebox
  use paraconfig,only:para_config,commconstr,initparapuresp
  USE sauvegardeT_mod,only:sauvegardeT
  use rasmolT_mod,only:rasmolT

#ifdef PARA
  USE mod_para,only:maj_atomes_frt_ftm
  use Tpara,only:grp_world,nprocs,myidsp,MPI_COMM_space,nprocspace,ierr,mpi_comm_world,&
       &NDM_MPI_REAL_DOUBLE,para_space_config,status,comm_space,mpi_world
  USE init_vois_mod,only: init_voisinage
#else
  use Tpara,only:myidsp,nprocspace,para_space_config
#endif
   use read_val,only:ltabvois
 USE init_simple_mod,only:init_simple
 USE init_pot_mod,only:init_pot
 USE sauvegardeT_mod,only:sauvegardeT
 USE dmloop_lpr_mod,only: dmloop_lpr
 use babar_def_mod
  use Tpara,only:mpi_communicator
  
  implicit none
  

  type(mpi_communicator),pointer::comm_master
  type(babar_config), allocatable,target:: babarloc(:)
  type(babar_config), pointer:: babarcur,babar1,babar2
  
  type(para_config),target::parababar ! division de parapath en 2*espace
  type(para_space_config)::pscbabar

  integer::ntempbabar,processes ! nombre de température
  integer::nbabarprocs! nombre de process gérant les températures (par défaut= processes)
  real(double)::bbtempmin,bbtempmax

  logical          :: ok

  integer::ntbbpp,itbbpp! nombre de température par process
  logical::lbigmaster !(true= master du calcul complet)
  
  type(atom_config_e),pointer::atconfb,atconfb1,atconfb2 !type derive atom_config du systeme a n atomes
  type(cell_config),pointer:: cellb !type derive cell_config du systeme a n atomes
  type(box_config_lpr),pointer::boxb
  type(atom_config_e),allocatable,target::config_atom_b(:) !type derive atom_config du systeme a n atomes
  type(cell_config),allocatable,target:: config_cell_b(:) !type derive cell_config du systeme a n atomes
  type(box_config_lpr),allocatable,target:: config_box_b(:) !type derive cell_config du systeme a n atomes

  type(para_space_config),pointer::pscbb
  type(para_space_config),allocatable,target::config_psc_b(:)

  integer:: itbtherm,itbprod
  integer::exchange_attempts,exchange_accepted_attempts
  
  bacontains

  subroutine init_mpi_babar

    integer::itbbtot,itbbpp
    character*80::namef,nameo
    character*6::extension
    real(double)::tempcur


    betamin=1/(bk*bbtempmin)
    betamax=1/(bk*bbtempmax)
    delta_beta=(betamax-betamin)/(processes-1)
    ALLOCATE(exchange_accepted(1:processes),id2babar(1:processes),rank2id(0:processes-1),exchange_attempted(1:processes))
    id2babar(:)%rank=0 ; id2babar%itbb=0
!    if (babartot(processes-1)%temp()-babartot(processes)%temp().lt.1.0) then
       write(6,*)' delta >=1K stop'
!       call arret_ndm
!    end if
       
    if (mod(nprocs,nbabarprocs).ne.0) then
       if (rang==0)  write(6,*)'nprocs/nbabarprocs <>0 STOP'
       call MPI_FINALIZE(ierr)
       call arret_ndm
    else

    end if
       
    if (mod(processes,nbabarprocs).ne.0) then
       if (rang==0)      write(6,*)'ntempbabar/nbabarprocs <>0 STOP'
       call MPI_FINALIZE(ierr)
       call arret_ndm
    else
       ntbbpp=processes/nbabarprocs
       allocate( babarloc(ntbbpp))
       allocate(config_atom_b(ntbbpp))
       config_atom_b%llangevin=.true.
       allocate(config_cell_b(ntbbpp))
       allocate(config_box_b(ntbbpp))
       allocate(config_psc_b(ntbbpp))
    end if

    
    parababar%mpi_orig%nproc=nprocs
    parababar%mpi_orig%rank=rang
    parababar%nimage=nbabarprocs
    call MPI_COMM_DUP(MPI_COMM_WORLD,parababar%mpi_orig%comm,ierr)
    call MPI_COMM_GROUP(parababar%mpi_orig%comm,parababar%mpi_orig%group,ierr)
    call commconstr(parababar)

    myidsp=parababar%mpi_image%rank
    call MPI_COMM_free(mpi_comm_space,ierr)
    MPI_COMM_space=parababar%mpi_image%comm
    call comm_space%init(MPI_COMM_SPACE)
    nprocspace=parababar%mpi_image%nproc
    lmasterb=parababar%lmaster

    comm_master=>parababar%mpi_master


    
    if (parababar%mpi_orig%rank==0) lbigmaster=.true.
    do itbbpp=1,ntbbpp
       idcur=itbbpp+(parababar%image)*ntbbpp
       if (lmasterb) then
          id2babar%rank=parababar%image
          id2babar%itbb=itbbpp
       end if
       babarcur=> babarloc(itbbpp)
       babarcur%indice=ntbbpp*parababar%image+itbbpp
       babarcur%beta=betamin+delta_beta*(itbbtot-1)*(betamax-betamin)/(processes-1)
       unitwb=1000+babarcur%temp()
       nameo=fnam(1:lenfnam)
       write(extension,'(i6.6)') int(babarcur%temp())
       namef=trim(nameo)//'.'//trim(extension)//'K.out'
       open(unit=unitwb, file=namef, status='unknown')
    end do

    call mpi_worl%sum(id2babar(:)%rank)
    call mpi_worl%sum(id2babar(:)%itbb)

    
    exchange_accepted  = 0
    exchange_attempted = 0
    self_rank=parababar%image
    self_id=self_rank+1
    replica_id= self_id
    
  end subroutine init_mpi_babar

  subroutine init_babar(rv)
    real(double),intent(in)::rv
    integer::itbbpp,itemp,i
    character*80::filename,nameo
    character*6:: extension
    real(double)::tinitb

    exchange_attempts=0
    call init_pot


    ALLOCATE(replica_betas(1:processes))
    FORALL (i=1:processes) replica_betas(i) = betamin + delta_beta * real(i-1,DOUBLE)/real(processes-1,DOUBLE)

    
    do itbbpp=1,ntbbpp
       babarcur=>babarloc(itbbpp)
       itemp=babarcur%indice
       atconfb=>config_atom_b(itbbpp)
       cellb=>config_cell_b(itbbpp)
       boxb=>config_box_b(itbbpp)
       pscbb=>config_psc_b(itbbpp)
       call atconfb%init(0,imm_glob,ltabvois,rvois=rv)
          nameo=fnam(1:lenfnam)
       if (lmultin.eqv..true.) then
          write(extension,'(i6.6)') int(babarcur%temp())
          filename=trim(nameo)//'.'//trim(extension)//'K'
       else
          filename=nameo
       end if

       tinitb=babarcur%temp()
       unitwb=1000+babarcur%temp()
       if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
          latcomp=.false.
       else
          latcomp=.true.
       end if
       write(6,*)int(babarcur%temp()), 'FILENAME', filename,igen
       call init_simple(atconfb,cellb,boxb,filename=trim(filename),psc=pscbb,linitpot=.false.,tinitr=tinitb)
    end do
  end subroutine init_babar
    
  subroutine babar
    use dmloop_pilot_mod,only:dmloop_pilot
    integer::itbbpp,iter,itapp
    character*80::fnamecout,nameo
    character*6:: extension
    integer::formatsauv=5



    dmtype=88
    if (itbtherm.gt.0) then
       
       iteration=0
       itloopmax=itbtherm
       timeloopmax=timemax
       do while (iteration.lt.itloopmax)
          iteration=iteration+1
          do itbbpp=1,ntbbpp
             iteration=1
             ite4=0
             babarcur=>babarloc(itbbpp)
             atconfb=>config_atom_b(itbbpp)
             cellb=>config_cell_b(itbbpp)
             boxb=>config_box_b(itbbpp)
             pscbb=>config_psc_b(itbbpp)
             if (lmasterb.eqv..true.)write(6,*)'BABAR RUN',babarcur%temp()
             if (lmasterb.eqv..true.)lwrtb=.true.
             unitwb=1000+babarcur%temp()
             text=babarcur%temp()
#ifdef PARA
             if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
                call maj_atomes_frt_ftm(atconfb,cellb,boxb,pscbb)
             end if
#endif
             !       do iteration=0,itloopmax
             
             call initbabarloop1 (atconfb,cellb,boxb,pscbb)
             
             do while ((iteration.lt.itloopmax).and.(MOD(iteration,4).ne.0))
                call babarloop1 (atconfb,cellb,boxb,pscbb)
             end do
          end do
          
          !      IF (MOD(iteration,4) == 0) THEN
          call mpi_world%barrier
          if (lmasterb) then
             do itbbpp=1,ntbbpp
                CALL rex_id_exchange(boxb%Htot, ok)
             !      end IF
             end do
          end if
          do itbbpp=1,ntbbpp
             if (lmasterb) then
                babacur%beta  = replica_betas(babacur%indice)
                text=babacur%temp()
             end if
             call comm_space%bcast(text,0)
          end do
          
       end do

       nameo=fnam(1:lenfnam)
       write(extension,'(i6.6)') int(babarcur%temp())
       fnamecout=trim(nameo)//'.'//trim(extension)//'K.cout'
       itapp=int(babarcur%temp()); write (6,*)'ITAPP',itapp,ivisu
       if ((lspacendm).and.(nprocspace.gt.1)) then
             call rasmolT(atconfb,boxb,itapp,'K',latcomp=.false.)
             call sauvegardeT(atconfb,cellb,boxb,formatsauv,fnamecout,latcomp=.false.)
          else
             call rasmolT(atconfb,boxb,itapp,'K',latcomp=.true.)
             call sauvegardeT(atconfb,cellb,boxb,formatsauv,fnamecout,latcomp=.true.)
          end if
       
       end if
 end if

 
end subroutine babar


  subroutine initbabarloop1(atpr,celndm,boxndm,psc)
    use Parrinello_Rahman, only:initlpr,unitw,lwrt,tinitbox
    type(para_space_config)::psc
    type(box_config_lpr)::boxndm
    class(atom_config_d)::atpr
    type(cell_config):: celndm

    tinitbox=babarcur%temp()
    if (lmasterb.eqv..true.) then
       lwrt=.true.
       unitw=unitwb
    else
       lwrt=.false.
    end if
    call initlpr(atpr,celndm,boxndm,psc)

  end subroutine initbabarloop1

  subroutine babarloop1(atpr,celndm,boxndm,psc)
    use Parrinello_Rahman,only:pr1
    USE analyseT_mod,only: analyseT
    USE controleT_mod,only: controleT
    type(para_space_config)::psc
    type(box_config_lpr)::boxndm
    class(atom_config_d)::atpr
    type(cell_config):: celndm
    integer,save::icall=0
    icall =icall+1
    do while (ite4.lT.4)
       write(6,*)'ITER ICAL',ite4,iteration,icall
       call pr1(atpr,celndm,boxndm,psc)
       timel=timel+tstep
       boxndm%Htot=potist+boxndm%ucell
       call analyseT (atpr,celndm,boxndm,psc,lwrtr=lwrtb,unitwr=unitwb)
    end do

    return
  end subroutine babarloop1


  subroutine rex_id_rank
    integer :: i , err
    rank2id = 0

    if(lmasterb) then
        rank2id(self_rank) = self_id
       call comm_master%sum(rank2id)
        do i = 0,processes-1
           id2rank(rank2id(i)) = i
        enddo
     end if
   end subroutine rex_id_rank
  


    subroutine rex_id_exchange(energy, ok)
        REAL(DOUBLE), intent(in) :: energy
        logical, intent(out) :: ok
        integer peer_rank , peer_id
        integer last_id,i
        call rex_id_rank
        self_id= babarcur%indice  ! C'est l'indice de la temperature
! IL FAUT CONSTRUIRE UN TABLEAU DES PEER AVEC LEUR RANK LEUR ITBB DE TAILLE NTBBPP ETUN TABLEAU DES ENERGIE TOTALE DE TAILLE PROCESSES
        if (mod(exchange_attempts, 2) == 0) then
            if (mod(self_id, 2) == 0) then
                peer_id = self_id + 1
            else
                peer_id = self_id - 1
            end if
        else
            if (mod(self_id, 2) == 0) then
                peer_id = self_id - 1
            else
                peer_id = self_id + 1
            end if
        end if

        last_id = self_id

        if (peer_id > 0 .and. peer_id <= processes ) then
           peer_rank = id2babar(peer_id)%rank
           peer_itbb=id2babar(peer_id)%itbb
        else
            peer_rank = -1
         endif
COMMENT ENVOYER la nature du peer ? Il faut envoyer le RANK ET le ITBB
        if (peer_rank >= 0 .and. peer_rank < processes) then
            if (self_rank > peer_rank) then
                call rex_exchange_as_leader(energy, peer_rank)
            else
                call rex_exchange_as_follower(energy, peer_rank)
            end if
        end if
        exchange_attempts = exchange_attempts + 1
        exchange_attempted(self_id) = exchange_attempted(self_id) + 1
        ok = self_id /= last_id
        babacur%indice=self_id
        if (ok) then
            exchange_accepted_attempts = exchange_accepted_attempts + 1
          !  exchange_accepted(peer_id) = exchange_accepted(peer_id) + 1
            exchange_accepted(self_id) = exchange_accepted(self_id) + 1
        endif

    end subroutine rex_id_exchange


  subroutine rex_exchange_as_leader(self_energy, peer_rank)
    REAL(DOUBLE), intent(in) :: self_energy
    integer, intent(in) :: peer_rank
    REAL(DOUBLE) self_beta, peer_beta, peer_energy
    integer peer_id, tmp_id
    REAL(DOUBLE) u_rand
    REAL(DOUBLE) metropolis
    REAL(DOUBLE) message(2)
    integer err
    integer::sendid
    call comm_master%recv(message,peer_rank,101)
!    call mpi_recv(message, 2, MPI_DOUBLE_PRECISION, peer_rank, TAG_EXCHANGE, MPI_COMM_WORLD, status_, err)

    !     write(*,*) ' status_ ', status_
    !     write(*,*) ' err     ', err
    !    write(*,*) ' message     ', message
    peer_id = int(message(1))
    peer_energy = message(2)
 !   if (lkappa) then
 !      self_beta = -replica_beta_kappa_epsilons(self_id)
 !      peer_beta = -replica_beta_kappa_epsilons(peer_id)
 !   else
       self_beta = replica_betas(self_id)
       peer_beta = replica_betas(peer_id)
 !   endif
!       call rex_rand_uniform(u_rand)
       call random_number(u_rand)
    metropolis = exp((self_beta - peer_beta) * (self_energy - peer_energy))

    !        write(*,*) 'metropolis ' , metropolis
    !        write(*,*) ' peer_beta  peer_energy err', peer_beta , peer_energy
    !        write(*,*) ' self_beta  self_energy ',self_beta , self_energy
    !        write(*,*) ' self_id peer_id' ,  self_id , peer_id

    if (u_rand < metropolis) then
       tmp_id = self_id
       self_id = peer_id
       peer_id = tmp_id
    end if
    sendid=peer_id
!    message(1) = REAL(peer_id,DOUBLE)
!    message(2) = peer_energy
    call comm_master%send(sendid,peer_rank,102)
!    call mpi_send(message, 2, MPI_DOUBLE_PRECISION, peer_rank, TAG_EXCHANGE, MPI_COMM_WORLD, err)

    return
  end subroutine rex_exchange_as_leader

  subroutine rex_exchange_as_follower(self_energy, peer_rank)
    REAL(DOUBLE), intent(in) :: self_energy
    integer, intent(in) :: peer_rank
    REAL(DOUBLE) message(2)
    integer err
    integer::recvid

    message(1) = REAL(self_id,DOUBLE)
    message(2) = self_energy

    !        write(*,*) 'message from follower ',message
    call comm_master%send(message, peer_rank,101)
    call comm_master%recv(recvid, peer_rank,102)
!    call mpi_send(message, 2, MPI_DOUBLE_PRECISION, peer_rank, TAG_EXCHANGE, MPI_COMM_WORLD, err)
!    call mpi_recv(message, 2, MPI_DOUBLE_PRECISION, peer_rank, TAG_EXCHANGE, MPI_COMM_WORLD, status_, err)
    self_id = int(message(1))

    return
  end subroutine rex_exchange_as_follower

    

  
end module babar_mod
  
