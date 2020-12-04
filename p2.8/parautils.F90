  module parautils

  use paraconfig,only:para_config
#ifdef PARA
  use Tpara,only:NDM_MPI_REAL_DOUBLE
  USE mod_para,only:maj_atomes_frt_ftm,mpi_comm_world
#endif
  use T_kind_param_m, ONLY:  double
  USE decoupage_mod,only: decoupage

  use atomconfig,only: atom_config,atom_config_d
  USE boxconfig,only:box_config
  USE cellconfig,only:cell_config,caltabtC
  use calfo_mod,only:calfo
  use constrconf_mod,only:repartition
   USE caltabi_mod,only: caltabi
  implicit none
contains
  subroutine initloc(atcomp,cellcomp,atloc,celloc,box,div,rum,lperiod)
    USE setcell,only:setcellconf
    type(atom_config_d),intent(in),target::atcomp
    type(cell_config),intent(in),target::cellcomp
    type(box_config)::box
    type(atom_config_d),pointer::atloc
    type(cell_config),pointer::celloc
    type(para_config),intent(in)::div
    real(double),intent(in)::rum
    logical::lperiod
    integer::ierr,iun
    if (div%npim.gt.1) then
       call cellcomp%copy_cell(celloc)
       call decoupage(div%npim,0,celloc,atloc)
       call repartition(atcomp,atloc,box,celloc,div=div) ! mettre les éléments de la répartition dans un type
    call setcellconf(celloc,atloc,box,atcomp%im,rum)
    else
       atloc=>atcomp
       celloc=>cellcomp
    end if
    call caltabtC(celloc,atloc,lperiod,box)
     end subroutine initloc
  
  subroutine pointer_caltabt_calfo(sig,potist,atcomp,cellcomp,box,atloc,celloc,div,lperiod,ltabvois,it,itetabvois,lchg,ii)

    
    real(double)::sig(3,3),potist
    type(atom_config_d),target::atcomp
    type(cell_config),target::cellcomp
    type(box_config)::box
    type(para_config)::div
    type(atom_config_d),pointer::atloc
    type(cell_config),pointer::celloc
    logical,intent(in)::lperiod
    logical,optional,intent(in)::ltabvois
    integer,optional,intent(in)::itetabvois,it,ii
    
    logical,optional::lchg
    integer::ierr,i
    logical::lchange=.true.
    integer::iun
    
    if(present(lchg))lchange=lchg
#ifdef PARA
    if (lchange) then
       if (div%npim.gt.1) then
          call atcomp%master2loc(atloc,div)
       else
          atloc=>atcomp
          celloc=>cellcomp                 
       end if
    end if

    
#else
    atloc=>atcomp
    celloc=>cellcomp

#endif


    call caltabtC(celloc,atloc,lperiod,box)
    if (present(ltabvois)) then
       if (ltabvois.and.((it==1).or.(mod(it,itetabvois)==0)))&
            &call caltabi(atloc,celloc,box)
    end if

#ifdef PARA
!    if (lchange) then  ! même sans changement il faut mettre à jour pour initialisze les tableaux NDM like de mod_para pour maj_tab_density
    if (div%npim.gt.1) then
       call maj_atomes_frt_ftm(atloc,celloc)
    end if

!    end if
#endif
       
       CALL CalFo(sig,potist,atloc,celloc,box,t_sigma=.true.)
    
#ifdef PARA
    if (div%npim.gt.1) then
       call atloc%vers_master(atcomp,div)
    else
!       atcomp=atloc    
!       cellcomp=celloc
    end if
    
#else
    atcomp=atloc
    cellcomp=celloc
#endif
    return
  end subroutine pointer_caltabt_calfo



    
end module parautils



