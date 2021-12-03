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
  use gen_com_m,only:lspacendm
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
    logical,pointer::ltabvois_p
    integer,pointer::itetabvois_p,it_p
    logical,pointer::lperiod_p
    logical,pointer::lchg_p


    integer,parameter::STOP_TAG=0
    integer,parameter::FORCE_TAG=1
    integer,parameter::WORKER_TAG=-1

  contains
  subroutine initloc(atcomp,cellcomp,atloc,celloc,box,div,rum,lperiod,ldistrib,psc)
    USE setcell,only:setcellconf
    class(atom_config_d),intent(in),target::atcomp
    type(cell_config),intent(in),target::cellcomp
    type(box_config)::box
    type(para_space_config)::psc
    class(atom_config),pointer::atloc
    type(cell_config),pointer::celloc
    type(para_config),intent(in)::div
    real(double),intent(in)::rum
    logical::lperiod
    logical,optional,intent(in)::ldistrib
    logical::ldistr
    integer::ierr,iun
    ldistr=.false.
    if (present(ldistrib))ldistr=ldistrib
    
!    if ((div%mpi_image%nproc.gt.1).and.(lspaceNDM.eqv..true.)) then
    if (div%mpi_image%nproc.gt.1) then
       if (ldistr) then
          call atcomp%send2all(0,div%mpi_image)
       endif
       if (lspaceNDM.eqv..true.) then
          call cellcomp%copy_cell(celloc,box)
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

    call caltabtC(celloc,atloc,lperiod,box)

  end subroutine initloc

!*****
  subroutine pointer_caltabt_calfo(sig,potist,atcomp,cellcomp,box,atloc,celloc,div,&
       &lperiod,ltabvois,it,itetabvois,lchg,psc,caracT)
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
    logical,optional,intent(in)::ltabvois
    integer,optional,intent(in)::itetabvois,it
    logical::lperiod
    logical,optional::lchg
    character(len=*),optional,intent(in)::caracT
    character(len=26)::caracm2l,caracvm
    integer::ierr,i,ierror
    logical::lchange=.true.
    real(double)::atl(3,3)
    integer::iun

    if (.not.present(caracT)) then
       caracm2l='xfniewdlpvrugas'
       caracvm=caracm2l
    else
       caracm2l=caracT//'npft'
       caracvm=caracT//'npxt'       
    end if

    if(present(lchg))lchange=lchg
#ifdef PARA
    if (lchange) then
       if (div%mpi_image%nproc.gt.1)then

          if (lspaceNDM.eqv..true.) then
             call atcomp%master2loc(atloc,div)
   
          else
             call atcomp%send2all(0,div%mpi_image)
             atloc=>atcomp
             celloc=>cellcomp
             
       
          end if
       else
          atloc=>atcomp
          celloc=>cellcomp                 
       end if
    end if
#else
    atloc=>atcomp
    celloc=>cellcomp

#endif
     call periodbox (box,atloc)

    call caltabtC(celloc,atloc,lperiod,box)
    if (present(ltabvois)) then
       if (ltabvois.and.((it==1).or.(mod(it,itetabvois)==0)))&
            &call caltabi(atloc,celloc,box)
    end if

#ifdef PARA
!    if (lchange) then  ! même sans changement il faut mettre à jour pour initialisze les tableaux NDM like de mod_para pour maj_tab_density

    if ((div%mpi_image%nproc.gt.1).and.(lspaceNDM.eqv..true.)) then
       call maj_atomes_frt_ftm(atloc,celloc,box,psc)
    end if

!    end if
#endif
!!$

    CALL CalFo(sig,potist,atloc,celloc,box,t_sigma=.true.,psc=psc)
#ifdef PARA

    if ((div%mpi_image%nproc.gt.1).and.(lspaceNDM.eqv..true.)) then
          call atloc%vers_master(atcomp,div)
          
       end if
             
#else
!    atcomp=atloc
!    cellcomp=celloc
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
  
  

  subroutine tolstoi (tag,div,carac)
    !https://www.youtube.com/watch?v=IsvfofcIE1Q

    character,intent(in)::carac
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
             call depeche_mode (div,carac)! Servants enter pointer_caltabt_calfo_driver ! 
             cycle ! upon return cycle to wait next call
          end if
       end select
       write(6,*)'you shoulndt be here', div%mpi_orig%rank
       call arret_ndm
    end do
       
  end subroutine tolstoi

  subroutine depeche_mode(div,carac,lchgboxT)
!    use mpi
    use gen_com_m,only:rang
    type(para_config)::div
    character (len=*)::carac
    logical, optional,intent(in)::lchgboxT
    logical::lchgbox
    real(double)::atl(3,3),ex
    integer::ierror,i1,i2
    integer,save::nc=0
    lchgbox=.false.
    ex=rang
    nc=nc+1
    if (present(lchgboxT))lchgbox=lchgboxT
    
    
    if (div%lmaster) then  ! Go in tolstoi get the servants
       call tolstoi(FORCE_TAG,div,carac)
    end if
    call div%mpi_image%bcast(0,lchgbox)
    !    call div%mpi_image%bcast(0,box_p%at) CA MARCHE PAS AVEC LE POINTEUR !
    if (lchgbox)then
       atl=box_p%at(1:3,1:3)
       call div%mpi_image%bcast(0,atl)
       call updatebox(box_p,atl)
    end if
    call pointer_caltabt_calfo(sig_p,potist_p,atcomp_p,cellcomp_p,box_p,atloc_p,celloc_p,div_p,lperiod_p,&
         &ltabvois_p,it_p,itetabvois_p,lchg_p,psc_p,carac)
    it_p=it_p+1
    return ! master returns to "main", servants return to tolstoi to wait for next call
  end subroutine depeche_mode
    

subroutine driver_caltabt_DM(sigcf,potistcf,atcf,celcf,boxcf,psc,lperiod)

  use Tpara,only:nprocspace
  use gen_com_m,only:it,itetabvois,itesigma
    class(atom_config_d),intent(inout),target::atcf
    type(cell_config),intent(inout),target::celcf
    type(box_config),intent(inout)::boxcf
    type(para_space_config)::psc
    real(double),intent(out)::potistcf,sigcf(3,3)
    logical,intent(in)::lperiod

    logical  ::test_sigma
    
    !conditions periodiques

    call periodbox (boxcf,atcf)
    ! repartition des atomes dans la nouvelle boite
!!$    if (.not.lprahman) then
!!$       if (itab/=0) then
!!$          if (mod(it,itab)==0) then
             call caltabtC(celcf,atcf,lperiod,boxcf)
!!$          endif
!!$       endif
!!$    end if
    if (atcf%ltabvois.and.mod(it,itetabvois)==0) then
       call caltabi(atcf%atom_config,celcf,boxcf)
    end if


#ifdef PARA
if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs
       call maj_atomes_frt_ftm(atcf,celcf,boxcf,psc)
    end if
#endif


!!$    ! Force calculation
!!$!    jq=0.0
!!$    if (itesigma>0) test_sigma=(mod(it,itesigma)==0)
!!$    CALL CalFo(sigcf,potistcf,atcf,celcf,boxcf,t_sigma=test_sigma,psc=psc)
    return
  end subroutine driver_caltabt_DM
    
end module parautils



