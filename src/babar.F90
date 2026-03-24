module babar_mod
  !NDM is © 2021, Jean-Paul Crocombette, CEA Saclay, SRMP
  !NDM is published and distributed under the Academic Software License v1.0 (ASL).
  !NDM is distributed in the hope that it will be useful for non-commercial academic research, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the ASL for more details.
  !You should have received a copy of the ASL along with this program; if not, write to jpcrocombette@cea.fr. It is also published at https://github.com/jpcroc/NDM/blob/ndm2025/LICENSE.md.
  !You may contact the original licensor at jpcrocombette@cea.fr.



  USE T_kind_param_m, ONLY:  double
  USE gen_com_m,only: uwrt,lwrt,rang,fnam,lenfnam,imm_glob,dmtype,itloopmax,itmax,timemax,timeloopmax,iteration,&
       &fnam,lenfnam,lspacendm,latcomp,ivisu,igen,text,bk,timel,tstep,potist,kine,lbabar,lregular,timel,tstep
  USE atomconfig,only:atom_config_e,atom_config_d
  USE cellconfig, only:cell_config, caltabtC
  USE arret_ndm_mod,only:arret_ndm
  USE boxconfig,only:box_config,periodbox,box_config_lpr,updatebox
  use paraconfig,only:para_config,commconstr,initparapuresp
  USE sauvegardeT_mod,only:sauvegardeT
  use rasmolT_mod,only:rasmolT
  USE analyseT_mod,only: analyseT
!  USE controleT_mod,only: controleT

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
  USE dmloop_lpr_mod,only: dmloop_lpr
  use babar_def_mod
  use Tpara,only:mpi_communicator

  implicit none

  logical::lmultin,lbetagrid
  type(mpi_communicator),pointer::comm_babar
  type(babar_config),allocatable,target:: babarloc(:)
  type(babar_config), pointer:: babarcur
  type (babar_config),allocatable,target::babartot(:)

  
  
  type(para_config),target::parababar ! division de parapath en 2*espace

  !NPROCS=nbabarprocs*nprospace
  !ntempbabar=nbabarprocs*ntbbpp
  integer::ntempbabar ! nombre de température
  integer::nbabarprocs! nombre de process gérant les températures (par défaut= processes)
  real(double)::bbtempmin,bbtempmax

  integer::ntbbpp ! nombre de température par process
  logical::lbigmaster,lmaster !(true= master du calcul complet)

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
  integer,allocatable::unitwb(:)
  logical::lwrtb=.false.

contains

  subroutine init_mpi_babar

    integer::icalc,itbbpp,iproc,idc,i
    character*80::nameo
    character*80,allocatable::namef(:)
    character*6::extension
    lregular=.false.

    betamin=1/(bk*bbtempmax)
    betamax=1/(bk*bbtempmin)
    delta_temp=(bbtempmax-bbtempmin)/(ntempbabar-1)
    delta_beta=(betamax -betamin)/(ntempbabar-1)
    
    write(uwrt,*)'TEMP0', bbtempmin,bbtempmax,delta_temp
    ALLOCATE(exchange_accepted(1:ntempbabar),exchange_attempted(1:ntempbabar))
    ALLOCATE(id_beta(1:ntempbabar))
    ALLOCATE(id_temp(1:ntempbabar))
    ALLOCATE(replica_energie(1:ntempbabar))
    ALLOCATE(id2calc(1:ntempbabar))

    do idc=1,ntempbabar
       if (lbetagrid) then
          id_beta(idc) = betamin  + delta_beta * (idc-1)
          id_temp(idc)=1/(bk*id_beta(idc))
       else
          id_temp(idc) = bbtempmin + delta_temp * (idc-1)
          id_beta(idc)=1/(bk*id_temp(idc))
       end if
       id2calc(idc)=idc
       if (lwrt)write(uwrt,*)'TEMPERATURE GRID', id_temp(idc)
    end do

    !write(uwrt,*)' delta >=1K stop'
    !       call arret_ndm

    if (mod(nprocs,nbabarprocs).ne.0) then
       if (rang==0)  write(uwrt,*)'nprocs/nbabarprocs <>0 STOP'
       call MPI_FINALIZE(ierr)
       call arret_ndm
    else

    end if

    if (mod(ntempbabar,nbabarprocs).ne.0) then
       if (rang==0)      write(uwrt,*)'ntempbabar/nbabarprocs <>0 STOP'
       call MPI_FINALIZE(ierr)
       call arret_ndm
    else
       ntbbpp=ntempbabar/nbabarprocs
       allocate( babarloc(ntbbpp))
       allocate(config_atom_b(ntbbpp))
       config_atom_b%llangevin=.true.
       allocate(config_cell_b(ntbbpp))
       allocate(config_box_b(ntbbpp))
       allocate(config_psc_b(ntbbpp))
       allocate(unitwb(ntbbpp))
       allocate(namef(ntbbpp))
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
    lmaster=parababar%lmaster

    comm_babar=>parababar%mpi_master

    if (lmaster) then 

       if (parababar%mpi_orig%rank==0)then
          lbigmaster=.true.
          allocate (babartot(ntempbabar))
          icalc=0
          do itbbpp=1,ntbbpp
             icalc=icalc+1
             call babartot(icalc)%setbetatemp(temp=id_temp(icalc))
             babartot(icalc)%rank=0
             babartot(icalc)%itbb=itbbpp
             babartot(icalc)%indice=icalc
             babartot(icalc)%ictot=icalc
          end do
          do iproc=1,nbabarprocs-1
             do itbbpp=1,ntbbpp
                icalc=icalc+1
                call babartot(icalc)%setbetatemp(temp=id_temp(icalc))
                babartot(icalc)%rank=iproc
                babartot(icalc)%itbb=itbbpp
                babartot(icalc)%indice=icalc
                babartot(icalc)%ictot=icalc
             end do
          end do
       end if
       call babardistrib
       do itbbpp=1,ntbbpp
          babarcur=>babarloc(itbbpp)
          nameo=fnam(1:lenfnam)
          unitwb(itbbpp)=1000+int(babarcur%temp)
          !             extension=trim(extension)
          write(extension, '(i0)') int(babarloc(itbbpp)%temp)  ! No leading spaces
          namef(itbbpp) = trim(nameo) // '.' // trim(adjustl(extension)) // 'K.out'
          do i = 1, len_trim(namef(itbbpp))
             if (namef(itbbpp)(i:i) == ' ') namef(itbbpp)(i:i) = '_'
          end do
          open(unit=unitwb(itbbpp), file=namef(itbbpp), status='replace')
       end do


    end if

    call comm_space%bcast(0,unitwb)
!    call comm_space%bcast(0,namef)
  end subroutine init_mpi_babar
  subroutine babardistrib
    integer::icalc,iproc,itbbpp
    if (lbigmaster)then
       icalc=0
       do itbbpp=1,ntbbpp
          icalc=icalc+1
          babarloc(itbbpp)=babartot(icalc)
       end do
       do iproc=1,nbabarprocs-1
          do itbbpp=1,ntbbpp
             icalc=icalc+1
             call babartot(icalc)%send2proc(iproc,comm_babar)
          end do
       end do
    else
       do itbbpp=1,ntbbpp
          call babarloc(itbbpp)%recv(0,comm_babar)
       end do
    end if
  end subroutine babardistrib
  
  subroutine init_babar(rv)
    real(double),intent(in)::rv
    integer::itbbpp,i,icalc,iproc
    character*80::filename,nameo
    character*6:: extension
    real(double)::tinitb

    exchange_accepted  = 0
    exchange_attempted = 0
    exchange_attempts=0
    call init_pot
    do itbbpp=1,ntbbpp
       call mpi_world%barrier
       do iproc=0,nbabarprocs-1
          call mpi_world%barrier
          if (iproc==parababar%image) then
             write(uwrt,*)'process iproc goes',iproc
             icalc=parababar%image*ntbbpp+itbbpp
             if (lmaster) then
                babarcur=>babarloc(itbbpp)
                tinitb=babarcur%temp
             end if
             call comm_space%bcast(0,tinitb)
             atconfb=>config_atom_b(itbbpp)
             cellb=>config_cell_b(itbbpp)
             boxb=>config_box_b(itbbpp)
             pscbb=>config_psc_b(itbbpp)
             !       call atconfb%init(0,imm_glob,ltabvois,rvois=rv)
             !       if (itbbpp==1) then
             nameo=fnam(1:lenfnam)
             if (lmultin.eqv..true.) then
                write(extension,'(i6.6)') int(babarcur%temp)
                filename=trim(nameo)//'.'//trim(extension)//'K'
             else
                filename=nameo
             end if
             if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
                latcomp=.false.
             else
                latcomp=.true.
             end if
             if (lmaster)write(uwrt,*)babarcur%temp, 'FILENAME', filename,igen,unitwb(itbbpp)
             block
               integer::uwstk
               uwstk=uwrt
               uwrt=unitwb(itbbpp)
               call init_simple(atconfb,cellb,boxb,filename=trim(filename),psc=pscbb,linitpot=.false.,t&
                    &initr=tinitb,lwrtiR=lmaster,lrepart=.true.)
               !         call atconfb%print(unit=unitwb(itbbpp))
               !         call mpi_world%barrier
               !         call arret_ndm
               uwrt=uwstk
             end block
             
             !       else
             !          atconfb=config_atom_b(1)
             !          cellb=config_cell_b(1)
             !          boxb=config_box_b(1)
             !          pscbb=config_psc_b(1)
             !       end if
          end if
       end do
       call mpi_world%barrier
    end do

           
  end subroutine init_babar

  subroutine babar
    use dmloop_pilot_mod,only:dmloop_pilot
    integer::itbbpp,iter,itapp,irank,icalc,index,tempwrt
    character*80::fnamecout,nameo
    character*6:: extension
    integer::formatsauv=5
    integer::iproc,itemp
    integer,allocatable::indice_matrice(:,:)

    allocate(indice_matrice(0:nbabarprocs-1,1:ntbbpp))
    indice_matrice(:,:)=0


    if (itbtherm.gt.0) then
       itloopmax=itbtherm
       timeloopmax=timemax
       do itbbpp=1,ntbbpp
          icalc=parababar%image*ntbbpp+itbbpp
          if (lmaster) babarcur=>babarloc(itbbpp)
          atconfb=>config_atom_b(itbbpp)
          cellb=>config_cell_b(itbbpp)
          boxb=>config_box_b(itbbpp)
          pscbb=>config_psc_b(itbbpp)
#ifdef PARA
          if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
             call maj_atomes_frt_ftm(atconfb,cellb,boxb,pscbb)
          end if
#endif
          call initbabarloop1 (atconfb,cellb,boxb,pscbb,unitwb(itbbpp))
       end do
       iteration=0
       do while (iteration.lt.itloopmax)
          if (lwrt) write(uwrt,*)'ieration',iteration
          iteration=iteration+4
          replica_energie(:)=0.
          do itbbpp=1,ntbbpp
             icalc=parababar%image*ntbbpp+itbbpp
             if (lmaster)then
                babarcur=>babarloc(itbbpp)
             end if
             atconfb=>config_atom_b(itbbpp)
             cellb=>config_cell_b(itbbpp)
             boxb=>config_box_b(itbbpp)
             pscbb=>config_psc_b(itbbpp)
             if (lmaster) then 
                write(unitwb(itbbpp),*)'BABAR RUN',babarcur%temp,itbbpp,unitwb(itbbpp)
                lwrtb=.true.
                text=babarcur%temp
             end if
             call comm_space%bcast(0,text)
             call babarloop1 (atconfb,cellb,boxb,pscbb,unitwb(itbbpp),itbbpp)
             if (lmaster)then
                babarcur%energie=boxb%Htot
                replica_energie(babarcur%indice)=babarcur%energie
             end if
             call analyseT (atconfb,cellb,boxb,pscbb,lwrtanaR=lmaster,uwrtanaR=unitwb(itbbpp))        
          end do
          call mpi_world%barrier
          if (lmaster) then
             call comm_babar%sum(replica_energie)
             if (lbigmaster) then
                babartot(:)%energie=replica_energie(:)
                CALL rex_id_exchange
             end if
             call babardistrib
          endif
       end do


       do itbbpp=1,ntbbpp
          if (lmaster)then
             babarcur=>babarloc(itbbpp)
             tempwrt=int(babarcur%temp)
          end if
          call comm_space%bcast(0,tempwrt)
          atconfb=>config_atom_b(itbbpp)
          cellb=>config_cell_b(itbbpp)
          boxb=>config_box_b(itbbpp)
          pscbb=>config_psc_b(itbbpp)

          nameo=fnam(1:lenfnam)
          write(extension,'(i6.6)') tempwrt
          fnamecout=trim(nameo)//'.'//trim(extension)//'K.cout'
          itapp=tempwrt; write (unitwb(itbbpp),*)'tempwrt',itapp,ivisu
          if ((lspacendm).and.(nprocspace.gt.1)) then
             call rasmolT(atconfb,boxb,itapp,latcomp=.false.)
             call sauvegardeT(atconfb,cellb,boxb,formatsauv,fnamecout,latcomp=.false.)
          else
             call rasmolT(atconfb,boxb,itapp,latcomp=.true.)
             call sauvegardeT(atconfb,cellb,boxb,formatsauv,fnamecout,latcomp=.true.)
          end if
       end do
       do itemp=1,ntempbabar
          write(uwrt,*)'temp accept', id_temp(itemp),exchange_accepted(itemp)
       end do

    end if

  end subroutine babar


  subroutine initbabarloop1(atpr,celndm,boxndm,psc,unitwl)
    use Parrinello_Rahman, only:initlpr,tinitbox
    type(para_space_config)::psc
    type(box_config_lpr)::boxndm
    class(atom_config_d)::atpr
    type(cell_config):: celndm
    integer,intent(in)::unitwl
    if (lmaster) then
       tinitbox=babarcur%temp
    end if
    tinitbox=0
    if(lmaster)write(unitwl,*)'TBOX',tinitbox
    call comm_space%bcast(0,tinitbox)
    call initlpr(atpr,celndm,boxndm,psc,lwrtprR=lmaster,uwrtprR=unitwl)

  end subroutine initbabarloop1

  subroutine babarloop1(atpr,celndm,boxndm,psc,unitwl,itloc)
    use Parrinello_Rahman,only:pr1
    type(para_space_config)::psc
    type(box_config_lpr)::boxndm
    class(atom_config_d)::atpr
    type(cell_config):: celndm
    integer,intent(in)::itloc
    integer::unitwl
    integer::ite4
    integer::ustk
    logical::lstk
    integer,save::icall
    ite4=0
    icall=icall+1
    do while (ite4.lT.4)
       ite4=ite4+1
!       write(6,*)'ITER ICAL',ite4,iteration,icall
       call pr1(atpr,celndm,boxndm,psc)
       timel=timel+tstep
       boxndm%Htot= kine+potist+boxndm%EcellPR
       if (itloc==ntbbpp) timel=timel+tstep
    end do

    return
  end subroutine babarloop1


  subroutine rex_id_rank
    integer :: icalc , err,id

    id2calc(:)=0
    do icalc=1,ntempbabar
       id=babartot(icalc)%indice
       id2calc(id)=icalc
    end do
    if (any(id2calc==0)) then
       write(6,*)'id2calc=0',id2calc
       call arret_ndm
    end if

  end subroutine rex_id_rank

  subroutine rex_id_exchange
    integer ::peer_calc , peer_id,lead_id,lead_calc
    integer ::last_id,idl,icalc,nts2,itemp,id,exchang_accept_tot
    logical,allocatable::deja(:)
    real(double)::peer_energy,lead_energy,lead_beta,peer_beta,u_rand,metropolis
    exchang_accept_tot=0
    call rex_id_rank
    allocate(deja(ntempbabar))
    deja=.false.
    nts2=int((ntempbabar+0.1)/2)
    do idl=1,nts2
       
       if (mod(exchange_attempts, 2) == 0) then
          lead_id=2*idl+1
          peer_id=lead_id-1
       else
          lead_id=2*idl
          peer_id=lead_id-1
       end if
       if ((peer_id.gt.ntempbabar) .or.(lead_id.gt.ntempbabar)) cycle
       last_id = lead_id
!       write(6,*)'lead_peer_id', lead_id,peer_id
       lead_calc=id2calc(lead_id)
       peer_calc=id2calc(peer_id)
!       write(6,*)'lead_peer_CACL', lead_calc,peer_calc
       lead_energy=replica_energie(lead_id)
       lead_beta=id_beta(lead_id)
       peer_energy=replica_energie(peer_id)
       peer_beta=id_beta(peer_id)
       call random_number(u_rand)
       metropolis = exp((lead_beta - peer_beta) * (lead_energy - peer_energy))
!       write(uwrt,*)u_rand,metropolis
       if (u_rand < metropolis) then
          babartot(lead_calc)%indice=peer_id
          babartot(peer_calc)%indice=lead_id
          exchange_accepted(peer_id) = exchange_accepted(peer_id) + 1
          exchange_accepted(lead_id) = exchange_accepted(lead_id) + 1
          exchang_accept_tot=exchang_accept_tot+2
          id2calc(peer_id)=lead_calc
          id2calc(lead_id)=peer_calc
!          write(uwrt,*)'exchange accepted'
       else
!          write(uwrt,*)'exchange rejected'
       end if
       exchange_attempted(peer_id) = exchange_attempted(peer_id) + 1
       exchange_attempted(lead_id) = exchange_attempted(lead_id) + 1

    end do
    exchange_attempts = exchange_attempts + 1
    do icalc=1,ntempbabar
       id=babartot(icalc)%indice
       babartot(icalc)%temp=id_temp(id)
       babartot(icalc)%beta=id_beta(id)
    end do
    write(uwrt,*)'NTEMP BABAR, ACCEPTED EXCHANGES ',ntempbabar,exchang_accept_tot
 

  end subroutine rex_id_exchange




end module babar_mod
  
