module ndm2art2ndm

  USE arret_ndm_mod,only:arret_ndm
  USE atomconfig,only : atom_config
  USE cellconfig, only:cell_config
  USE boxconfig,only:box_config_lpr,box_config
  USE gen_com_m, ONLY: rang,latcomp,angst,inv_angst,lperiod,iteration,sig,potist,erg2ev
  use Tpara,only:para_space_config,nprocs,myidsp,nprocspace
#ifdef PARA
  use Tpara,only:para_space_config,nprocs,myidsp,nprocspace,comm_space ,mpi_comm_space,&
       &mpi_world,mpi_comm_world
#endif
  
  use endrunT_mod,only:endrunT
  use montecarlo_mod,only:lparapath,nparapath
  use defs,only:NATOMS,VECSIZE,typat,force,posref,boxref,use_local_forces,FCOUNTER,mincounter,&
       &new_event,restart,boundary,pos,box,typat,cell,invcell,iproc,nproc,t1,constr,COUNTER,unit6P,&
       &atdisp
  use var_pot,only:rumax
  use parautils,only:initloc,depeche_mode
    use paraconfig,only:para_config,commconstr,initparapuresp
  use update_invcell_mod,only: update_invcell

  implicit none
  type(para_config),target::parapath ! division de tous les procs en nparapath chemins calculés simultanément
  type(box_config),pointer::boxart
  type(atom_config),pointer::atcfart
  type(cell_config),pointer:: celart
  type(para_space_config),target::pscart
  class(atom_config),pointer::atcfartloc
  type(cell_config),pointer::celartloc
  type(cell_config),target:: cellcible ! ne sert qu'à faire pointer cellnebloc sur quelquechose
  type(atom_config),target::atcible
  logical,target:: lchg,lcalcvois

  logical :: lmaster,lbigmaster
  integer::ierr
contains

  subroutine init_mpi_art

#ifdef PARA
    if (rang==0)write(unit6P,*)'INITMPI_ART********************************************'
       write(unit6P,*)'JPNPNPP',nprocs,nparapath
    if (lparapath) then
       if (mod(nprocs,nparapath).ne.0) then
          write(unit6P,*)'nprocs/nparapath <>0 STOP'
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


    myidsp=parapath%mpi_image%rank  
    call MPI_COMM_free(mpi_comm_space,ierr)
    MPI_COMM_space=parapath%mpi_image%comm
    call comm_space%init(MPI_COMM_SPACE)
    nprocspace=parapath%mpi_image%nproc
    lmaster=parapath%lmaster ! ==myidsp=0 ==master of space calculations
    if (parapath%mpi_orig%rank==0) lbigmaster=.true. !== global master

#else

    parapath%mpi_orig%nproc=1
    parapath%mpi_orig%rank=0
    parapath%mpi_image%nproc=1
    parapath%lmaster=.true.

    myidsp=0
    nprocspace=1
    lbigmaster=.true.
    lmaster=.true.

#endif
    return
  end subroutine init_mpi_art

  subroutine ndm2art

!    type(para_space_config)::psc
!!$    class(box_config)::boxcf
!!$    class(atom_config)::atcf
!!$    type(cell_config):: celcf

    character(len=20) :: dummy, fname
    logical ::flag

    integer::ierror

    if (.not. restart) then
       inquire( file = COUNTER, exist = flag )
       if ( flag .and. iproc == 0 ) then
          open(unit=FCOUNTER,file=COUNTER,status='old',action='read',iostat=ierror)
          read(FCOUNTER,'(A12,I6)') dummy, mincounter
          close(FCOUNTER)
       else
          mincounter = 1000
       end if
    end if
    
    use_local_forces = .false.
    constr=0
    typat(1:NATOMS)=atcfart%ityp(1:NATOMS)
    pos(1:NATOMS)=angst*atcfart%xp(1,1:NATOMS)
    pos(1+NATOMS:2*NATOMS)=angst*atcfart%xp(2,1:NATOMS)
    pos(1+2*NATOMS:3*NATOMS)=angst*atcfart%xp(3,1:NATOMS)
    boundary='T'
    cell(:,:)=boxart%at(:,:)*angst
    box=0 ; boxref=0 ! initilisation à 0 pour provoquer un plantage
    allocate(atdisp(natoms))
    atdisp=0
    call update_invcell( )


  end subroutine ndm2art
  subroutine art2ndm(atcf,boxcf,celcf,linit)

!    type(para_space_config)::psc
    class(box_config)::boxcf
    class(atom_config)::atcf
    type(cell_config):: celcf
    logical,optional ::linit
    logical ::lini=.false.
    character(len=20) :: dummy, fname
    logical ::flag

    integer::ierror,ipbc(3)
    ipbc=1
    if(present(linit))lini=linit

    if (lini) then
       call atcf%init(NATOMS,lreallocate=.true.)
    else
       if (atcf%im.ne.NATOMS) then
          write(unit6P,*)'ART2NDM incosistency between NATOMS and atcf%im stop'
          call arret_ndm
       end if
    end if
    atcf%ityp(1:NATOMS)=typat(1:NATOMS)
    atcf%xp(1,1:NATOMS)=pos(1:NATOMS)/angst
    atcf%xp(2,1:NATOMS)=   pos(1+NATOMS:2*NATOMS)/angst
    atcf%xp(3,1:NATOMS)=    pos(1+2*NATOMS:3*NATOMS)/angst
    boxcf%at(:,:)=cell(:,:)/angst
    call boxcf%init(boxcf%at,ipbc)
    call celcf%init(boxcf,celart%nox,celart%noy,celart%noz)



  end subroutine art2ndm

  subroutine calcforce_ndm(nat, posa,  forca, energy)
    use defs, only :  use_local_forces, local_ref_energy, global_ref_energy
    implicit none

    !Arguments
    integer,      intent(in)                            :: nat
    real(kind=8), intent(in),  dimension(3*nat)         :: posa
    real(kind=8), intent(out), dimension(3*nat)         :: forca
    real(kind=8), intent(out)                           :: energy

    integer,save::iteration

!    atcfart%ityp(1:NATOMS)=typ_a(1:NATOMS)
    atcfart%xp(1,1:NATOMS)=posa(1:NATOMS)*inv_angst
    atcfart%xp(2,1:NATOMS)=posa(1+NATOMS:2*NATOMS)*inv_angst
    atcfart%xp(3,1:NATOMS)=posa(1+2*NATOMS:3*NATOMS)*inv_angst
!    write(unit6P,*)'JPp1',rang,iproc
    call depeche_mode(parapath,.false.)

    energy=potist*erg2ev
!    write(unit6P,*)'JPenergy',energy
    forca(1:NATOMS)=(erg2ev/angst)*atcfart%fp(1,1:NATOMS)
    forca(1+NATOMS:2*NATOMS)=(erg2ev/angst)*atcfart%fp(2,1:NATOMS)
    forca(1+2*NATOMS:3*NATOMS)=(erg2ev/angst)*atcfart%fp(3,1:NATOMS)


    iteration=iteration+1
    return
  end subroutine calcforce_ndm

  subroutine set_pointers_art
    use parautils,only: psc_p,sig_p,potist_p,atcomp_p,cellcomp_p,box_p,div_p,atloc_p,celloc_p,lcv_p,&
         lperiod_p,lchg_p,lperiod_p,it_p
    character,target::carac(3)
    carac='xft'
    psc_p=>pscart
    sig_p=>sig
    potist_p=>potist
    atcomp_p=>atcfart
    cellcomp_p=>celart
    box_p=>boxart
    div_p=>parapath
    atloc_p=>atcfartloc
    celloc_p=>celartloc
    lcv_p=>lcalcvois
    it_p=>iteration
    lperiod_p=>lperiod
    lchg_p=>lchg
  end subroutine set_pointers_art

  subroutine init_mpi_art2(atcf,celcf,boxcf,psc,parapath)
    type(para_config),intent(in)::parapath
    type(atom_config),intent(in),target::atcf
    type(cell_config),target::celcf
    type(box_config),target::boxcf
    type(para_space_config)::psc
    integer :: ierr

    ! Call CPU_TIME to use for analysis of efficienty
    call CPU_TIME( t1 )

    ! These lines are used in a MPI environment. To access them, define MPI in the environment
    iproc=parapath%mpi_image%rank
    nproc=parapath%mpi_image%nproc
    lcalcvois=.false.
    atcfart=>atcf
    boxart=>boxcf
    celart=>celcf

    atcfartloc=> atcible
    celartloc=>cellcible
    pscart= psc
    call initloc(atcfart,celart,atcfartloc,celartloc,boxart,parapath,rumax,lperiod,ldistrib=.true.,psc=pscart,&
         &lcalcvois=lcalcvois,lboxchange=.false.)

    call set_pointers_art
!!$
!!$#ifdef MPI_VERSION_ART
!!$    iproc=0
!!$    nproc=1
!!$    ! This subroutine is in def.f90 and creates all mpi instances needed.
!!$    call mpi_group_creation()
!!$    write(unit6P,*) "Mpi_group_creation done check if the nproc and iproc are correct"
!!$    write(unit6P,*) iproc,nproc
!!$    call MPI_Barrier( MPI_COMM_WORLD, ierr )
!!$#else
!!$    iproc=0
!!$    nproc=1
!!$#endif


  end subroutine init_mpi_art2

end module ndm2art2ndm
