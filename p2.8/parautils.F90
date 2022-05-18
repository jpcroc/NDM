  module parautils
   USE arret_ndm_mod,only:arret_ndm
  use paraconfig,only:para_config
#ifdef PARA
  USE mod_para,only:maj_atomes_frt_ftm
    use Tpara,only:mpi_comm_world
#endif
  use Tpara,only:para_space_config
  use T_kind_param_m, ONLY:  double
  USE decoupage_mod,only: decoupage
  use gen_com_m,only:lspacendm,itetabvois,rang
  use atomconfig,only: atom_config,atom_config_d,atom_config_e
  USE boxconfig,only:box_config,periodbox,updatebox
  USE cellconfig,only:cell_config,caltabtC
  use calfo_mod,only:calfo
  use constrconf_mod,only:repartition
   USE caltabi_mod,only: caltabi
  implicit none

    type(para_space_config),pointer::psc_p
    real(double),pointer::sig_p(:,:),potist_p
    class(atom_config),pointer::atcomp_p
    type(cell_config),pointer::cellcomp_p
    type(box_config),pointer::box_p
    type(para_config),pointer::div_p
    class(atom_config),pointer::atloc_p
    type(cell_config),pointer::celloc_p
    integer,pointer::it_p
    logical,pointer::lcv_p
    logical,pointer::lchg_p,lperiod_p


    integer,parameter::STOP_TAG=0
    integer,parameter::FORCE_TAG=1
    integer,parameter::WORKER_TAG=-1

  contains
  subroutine initloc(atcomp,cellcomp,atloc,celloc,box,div,rum,lperiod,ldistrib,psc,lcalcvois)
    USE setcell,only:setcellconf
    class(atom_config),intent(in),target::atcomp
    type(cell_config),intent(in),target::cellcomp
    type(box_config)::box
    type(para_space_config)::psc
    class(atom_config),pointer::atloc
    type(cell_config),pointer::celloc
    type(para_config),intent(in)::div
    real(double),intent(in)::rum
    logical,intent(in),optional :: lcalcvois
    logical,intent(in)::lperiod
    logical,optional,intent(in)::ldistrib
    logical::ldistr
    integer::ierr,iun
    logical::lcalcv
    ldistr=.false.
    if (present(lcalcvois)) then
       lcalcv=lcalcvois
    else
       lcalcv=atcomp%ltabvois
    end if

    
    if (present(ldistrib))ldistr=ldistrib
!!$    call atcomp%print
!!$    write(6,*)'TOTO'
!!$    call atloc%print
    atloc%im_glob=atcomp%im_glob
!    if ((div%mpi_image%nproc.gt.1).and.(lspaceNDM.eqv..true.)) then
    if (div%mpi_image%nproc.gt.1) then
       if (ldistr) then
          call atcomp%send2all(0,div%mpi_image)
       endif
       if (lspaceNDM.eqv..true.) then
          call cellcomp%copy(celloc,box)
          call decoupage(div%mpi_image%nproc,0,celloc,atloc,lverbose=.false.,psc=psc)
          call repartition(atcomp,atloc,box,celloc) ! mettre les éléments de la répartition dans un type
          call setcellconf(celloc,atloc,box,rum,lverbose=.false.)
       else
          atloc=>atcomp
          celloc=>cellcomp
       end if
    else
       atloc=>atcomp
       celloc=>cellcomp
    end if
    if (lspacendm.and.div%mpi_image%nproc.gt.1) then
       call caltabtC(celloc,atloc,lperiod,box,psc=psc)
    else
       call caltabtC(celloc,atloc,lperiod,box)
    end if
    if ((lcalcv).and.(atloc%ltabvois)) call caltabi(atloc,celloc,box)
#ifdef PARA
    if ((div%mpi_image%nproc.gt.1).and.(lspaceNDM.eqv..true.)) then
       call maj_atomes_frt_ftm(atloc,celloc,box,psc)
    end if
#endif

  end subroutine initloc

!******************************************
  subroutine pointer_caltabt_calfo(sig,potist,atcomp,cellcomp,box,atloc,celloc,div,&
       &lperiod,lupdate,psc,lcalcvois)
    ! driver routine for caltabt calfo, period, caltabi, maj_atomes_frt_ftm
    ! lupdate= true ==> positions have changed update is needed:
    !   step1 :if atcomp is a complete set on master node update includes sending back to slaves (either all the configuration or only the (possibly ex-)local atoms)
    !    atloc, celloc are then the atoms and cells for each proc, either distributed or shared
    !   step2 : updates of atloc list of updates : PBC (period), cells (caltabtc), neighbours(if lcalcvois=true),
    !    passing of atoms from cpu to neihbours (maj_atomes_frt_ftm)
    !lupdate=false : no need to update atloc. That can occur when :
    !   -atloc has just been initiated (by initloc (e.g. neb.F90, montecarlo)
    !   -has not changed since last exit from this routine (strange but possible, no exemple yet at time of writing)
    !   -the change is dealt with outside of this routine : regular MD calls to driver_caltabt_para in dmloop, dyn_vverlet
    ! then
    !  call to calfo
    !  sending of the forces to the master if the atomic configurations are gathered on the master
    

#ifdef PARA
    use mpi
#endif
    type(para_space_config)::psc    
    real(double)::sig(3,3),potist
    class(atom_config),target::atcomp
    type(cell_config),target::cellcomp
    type(box_config)::box
    type(para_config)::div
    class(atom_config),pointer::atloc
    type(cell_config),pointer::celloc
    logical::lperiod
    logical::lupdate
    class(atom_config),pointer::atcalc
    type(cell_config),pointer::cellcalc
!    character(len=*),optional,intent(in)::caracT
    logical,optional::lcalcvois
    logical::lcalcv
!    character(len=26)::caracm2l,caracvm
    integer::ierr,i,ierror
    real(double)::atl(3,3)
    integer::iun
    integer,save::ncall=0
    lcalcv=.false.
    ncall=ncall+1
    if (present(lcalcvois))lcalcv=lcalcvois

#ifdef PARA
       if (div%mpi_image%nproc.gt.1)then
          if (lspaceNDM.eqv..true.) then
             if (lupdate) then
                call atcomp%master2loc(atloc,div)
             end if
             atcalc=>atloc
             cellcalc=>celloc
          else
             if (lupdate) then
                call atcomp%send2all(0,div%mpi_image)
             end if
             atcalc=>atcomp
             cellcalc=>cellcomp
          end if
       else
          atcalc=>atcomp
          cellcalc=>cellcomp                 
       end if
#else
       atcalc=>atcomp
       cellcalc=>cellcomp
#endif
       if (lupdate)   call driver_caltabt_para(atcalc,cellcalc,box,psc,lperiod,lcalcv)
       
    CALL CalFo(sig,potist,atcalc,cellcalc,box,t_sigma=.true.,psc=psc)
#ifdef PARA
    if ((div%mpi_image%nproc.gt.1).and.(lspaceNDM.eqv..true.)) then
       call atloc%vers_master(atcomp,div)
    end if
#endif
    return
  end subroutine pointer_caltabt_calfo


  
  subroutine initcomp(atcomp,cellcomp,atlocin,cellocin,box,div,lperiod,caracT,lorder)
    
    class(atom_config),intent(in)::atlocin
    type(cell_config),intent(in)::cellocin
    class(atom_config)::atcomp
    type(cell_config)::cellcomp
    type(box_config)::box
    type(para_config),intent(in)::div
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    logical,optional::lorder
    class(atom_config),pointer::atcdes
    type(atom_config),target::atb
    type(atom_config_d),target::atd
    type(atom_config_e),target::ate
    logical::lord
    integer::ierr,iun,i,j
    logical::lperiod
    lord=.false.
    if (present(lorder))lord=lorder
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugas'
    else
       carac=caracT//'np'       
    end if

       if ((div%mpi_image%nproc.gt.1).and.(lspaceNDM.eqv..true.)) then
          !    if (div%mpi_image%nproc.gt.1) then
          if (lord) then
             call cellcomp%init(box,cellocin%nox,cellocin%noy,cellocin%noz,cellocin%natperc,cellocin%ltpcel)
             select type (atcomp)
             type is (atom_config)
                atb=atcomp
                atcdes=>atb
             type is (atom_config_d)
                atd=atcomp
                atcdes=>atd
             type is (atom_config_e)
                ate=atcomp
                atcdes=> ate
             end select
             call atlocin%vers_master(atcdes,div,carac)

             if (div%mpi_image%rank==0) then
                
                do i=1,atcdes%im
                   j=atcdes%num_at_glob(i)
                   call atcdes%copy_atom(i,atcomp,j)
                end do
                
                
                call caltabtC(cellcomp,atcomp,lperiod,box)
             end if
             
          else
             call cellcomp%init(box,cellocin%nox,cellocin%noy,cellocin%noz,cellocin%natperc,cellocin%ltpcel)
             call atlocin%vers_master(atcomp,div,carac)
             if (div%mpi_image%rank==0) then
                call caltabtC(cellcomp,atcomp,lperiod,box)
             end if
             
          end if
    else
       call atlocin%copy_config(atcomp,lrescl=.false.)
       cellcomp=cellocin
       call caltabtC(cellcomp,atcomp,lperiod,box)
       if (atlocin%ltabvois) then
          call caltabi(atcomp,cellcomp,box)
       end if
    
    end if
  end subroutine initcomp
  
  

  subroutine tolstoi (tag,div)
    !https://www.youtube.com/watch?v=IsvfofcIE1Q

!    character,intent(in)::carac
    integer,intent(in)::tag
    type(para_config),intent(in)::div
    integer::newtag

    do
       call div%mpi_image%barrier !Servants wait for the master
       if (div%lmaster) newtag=tag  ! master is there from depeche_mode
       call div%mpi_image%bcast(0,newtag) ! all have newtag
       select case(newtag)
       case( STOP_TAG)
          return !back to "main" call  
       case(FORCE_TAG)
          if (div%lmaster) then
             return !master is back to pointer_caltabt_calfo_driver
          else
             call depeche_mode (div)! Servants enter pointer_caltabt_calfo_driver ! 
             cycle ! upon return cycle to wait next call
          end if
       end select
       write(6,*)'you shoulndt be here', div%mpi_orig%rank
       call arret_ndm
    end do
       
  end subroutine tolstoi

  subroutine depeche_mode(div,lchgboxT)
!    use mpi
    use gen_com_m,only:rang
    type(para_config)::div
!    character (len=*)::carac
    logical, optional,intent(in)::lchgboxT
    logical::lchgbox
    real(double)::atl(3,3),ex
    integer::ierror,i1,i2
    lchgbox=.false.
    ex=rang
    if (present(lchgboxT))lchgbox=lchgboxT
    
    
    if (div%lmaster) then  ! Go in tolstoi get the servants
       call tolstoi(FORCE_TAG,div)
    end if
    call div%mpi_image%bcast(0,lchgbox)
    !    call div%mpi_image%bcast(0,box_p%at) CA MARCHE PAS AVEC LE POINTEUR !
    if (lchgbox)then
       atl=box_p%at(1:3,1:3)
       call div%mpi_image%bcast(0,atl)
       call updatebox(box_p,atl)
    end if
    lcv_p=.false.
    it_p=it_p+1
    if ((atloc_p%ltabvois).and.(mod(it_p,itetabvois)==0)) lcv_p=.true.
    call pointer_caltabt_calfo(sig_p,potist_p,atcomp_p,cellcomp_p,box_p,atloc_p,celloc_p,div_p,lperiod_p&
         &,lchg_p,psc_p,lcalcvois=lcv_p)

    return ! master returns to "main", servants return to tolstoi to wait for next call
  end subroutine depeche_mode
    

subroutine driver_caltabt_para(atcf,celcf,boxcf,psc,lperiod,lcalcvois)

  use Tpara,only:nprocspace
  use gen_com_m,only:iteration,itetabvois,itesigma
    class(atom_config),intent(inout),target::atcf
    type(cell_config),intent(inout),target::celcf
    type(box_config),intent(inout)::boxcf
    type(para_space_config)::psc!    real(double),intent(in)::potistcf,sigcf(3,3)
    logical,intent(in)::lperiod
    logical,intent(in),optional::lcalcvois
    logical::lcalcv=.false.
    logical  ::test_sigma
    integer::i
    !conditions periodiques
    if (present(lcalcvois))lcalcv=lcalcvois
    call periodbox (boxcf,atcf)
    ! repartition des atomes dans la nouvelle boite
    if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       call caltabtC(celcf,atcf,lperiod,boxcf,psc=psc)
    else
       call caltabtC(celcf,atcf,lperiod,boxcf)
    end if
    if (atcf%ltabvois.and.(lcalcv)) then
       call caltabi(atcf,celcf,boxcf)
    end if
#ifdef PARA
if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs
       call maj_atomes_frt_ftm(atcf,celcf,boxcf,psc)
    end if
#endif
    
    return
  end subroutine driver_caltabt_para
    
end module parautils



