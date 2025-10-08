module babar_mod
  USE T_kind_param_m, ONLY:  double
  USE atomconfig,only:atom_config,atom_config_d, switch_atom
  USE cellconfig, only:cell_config, caltabtC
  USE arret_ndm_mod,only:arret_ndm
  USE boxconfig,only:box_config,periodbox,box_config_lpr,updatebox
  use paraconfig,only:para_config,commconstr,initparapuresp

#ifdef PARA
  USE mod_para,only:maj_atomes_frt_ftm
  use Tpara,only:grp_world,nprocs,myidsp,MPI_COMM_space,nprocspace,ierr,mpi_comm_world,&
       &NDM_MPI_REAL_DOUBLE,para_space_config,status,comm_space,mpi_world
  USE init_vois_mod,only: init_voisinage
#else
  use Tpara,only:myidsp,nprocspace,para_space_config
#endif
  USE gen_com_m,only: rang,fnam,lenfnam,lmultin,imm_glob,dmtype,itloopmax,itmax,timemax,timeloopmax,iteration,lwrtb,unitwb,&
       &fnam,lenfnam,lmasterb,lspacendm,latcomp
   use read_val,only:ltabvois
 USE init_simple_mod,only:init_simple
     USE init_pot_mod,only:init_pot
  implicit none
  
  type babar_config
     real(double)::temp,energie
     integer::indice
  end type babar_config


  
  type(babar_config), allocatable,target:: babartot(:),babarloc(:)
  type(babar_config), pointer:: babarcur,babar1,babar2
  
  type(para_config)::parababar,parabbsp ! division de parapath en 2*espace
  type(para_space_config)::pscbabar

  integer::ntempbabar ! nombre de température
  integer::nbabarprocs! nombre de process gérant les températures (par défaut= ntempbabar)
  real(double)::bbtempmin,bbtempmax



  integer::ntbbpp,itbbpp! nombre de température par process
  logical::lbigmaster !(true= master du calcul complet)
  
  type(atom_config_d),pointer::atconfb,atconfb1,atconfb2 !type derive atom_config du systeme a n atomes
  type(cell_config),pointer:: cellb !type derive cell_config du systeme a n atomes
  type(box_config_lpr),pointer::boxb
  type(atom_config_d),allocatable,target::config_atom_b(:) !type derive atom_config du systeme a n atomes
  type(cell_config),allocatable,target:: config_cell_b(:) !type derive cell_config du systeme a n atomes
  type(box_config_lpr),allocatable,target:: config_box_b(:) !type derive cell_config du systeme a n atomes

  type(para_space_config),pointer::pscbb
  type(para_space_config),allocatable,target::config_psc_b(:)
  
contains
  subroutine init_mpi_babar

    integer::itbbtot,itbbpp
    character*80::namef,nameo
    character*6::extension
    
    allocate(babartot(ntempbabar))
    do itbbtot=1,ntempbabar
       babartot(itbbtot)%indice=itbbtot
       babartot(itbbtot)%temp=bbtempmin+(itbbtot-1)*(bbtempmax-bbtempmin)/(ntempbabar-1)
    end do
    if (mod(nprocs,nbabarprocs).ne.0) then
       if (rang==0)  write(6,*)'nprocs/nbabarprocs <>0 STOP'
       call MPI_FINALIZE(ierr)
       call arret_ndm
    else

    end if
       
    if (mod(ntempbabar,nbabarprocs).ne.0) then
       if (rang==0)      write(6,*)'ntempbabar/nbabarprocs <>0 STOP'
       call MPI_FINALIZE(ierr)
       call arret_ndm
    else
       ntbbpp=ntempbabar/nbabarprocs
       allocate( babarloc(ntbbpp))
       allocate(config_atom_b(ntbbpp))
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

    
    if (parababar%mpi_orig%rank==0) lbigmaster=.true.
    do itbbpp=1,ntbbpp
       babarcur=> babarloc(itbbpp)
       babarcur%indice=ntbbpp*parababar%image+itbbpp
       babarcur%temp=bbtempmin+(babarcur%indice-1)*(bbtempmax-bbtempmin)/(ntempbabar-1)
       unitwb=1000+babarcur%temp
       nameo=fnam(1:lenfnam)
       write(extension,'(i6.6)') int(babarcur%temp)
       namef=trim(nameo)//'.'//trim(extension)//'K.out'
       open(unit=unitwb, file=namef, status='unknown')
    end do

  end subroutine init_mpi_babar

  subroutine init_babar(rv)
    real(double),intent(in)::rv
    integer::itbbpp,itemp
    character*80::name1
    character*84::filename
    character :: extension*4
    real(double)::tinitb
    name1=fnam(1:lenfnam)
    call init_pot
    do itbbpp=1,ntbbpp
       babarcur=>babarloc(itbbpp)
       itemp=babarcur%indice
       atconfb=>config_atom_b(itbbpp)
       cellb=>config_cell_b(itbbpp)
       boxb=>config_box_b(itbbpp)
       pscbb=>config_psc_b(itbbpp)
       call atconfb%init(0,imm_glob,ltabvois,rvois=rv)
       if (lmultin.eqv..true.) then
          write(extension,'(i4.4)') itemp
          filename=trim(name1)//trim(extension)
       else
          filename=trim(name1)
       end if

       tinitb=babarcur%temp
       unitwb=1000+babarcur%temp
       if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
          latcomp=.false.
       else
          latcomp=.true.
       end if
       call init_simple(atconfb,cellb,boxb,filename=trim(filename),psc=pscbb,linitpot=.false.,tinitr=tinitb)
    end do
    
  end subroutine init_babar
    
  subroutine babar
    use dmloop_pilot_mod,only:dmloop_pilot
    integer::itbbpp,iter
    dmtype=4
    itloopmax=itmax
    timeloopmax=timemax

    do itbbpp=1,ntbbpp
    iteration=0
       babarcur=>babarloc(itbbpp)
       atconfb=>config_atom_b(itbbpp)
       cellb=>config_cell_b(itbbpp)
       boxb=>config_box_b(itbbpp)
       pscbb=>config_psc_b(itbbpp)
       if (lmasterb.eqv..true.)write(6,*)'BABAR RUN',babarcur%temp
       if (lmasterb.eqv..true.)lwrtb=.true.
       unitwb=1000+babarcur%temp
#ifdef PARA
       if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
          call maj_atomes_frt_ftm(atconfb,cellb,boxb,pscbb)
       end if
#endif
       
       call dmloop_pilot(atconfb,cellb,boxb,pscbb,linit=.true.)
    end do
    
  end subroutine babar
  
end module babar_mod
  
