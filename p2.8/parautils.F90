  module parautils

  use paraconfig,only:para_config
#ifdef PARA
  USE mod_para,only:maj_atomes_frt_ftm
#endif
  use Tpara,only:para_space_config,mpi_comm_world
  use T_kind_param_m, ONLY:  double
  USE decoupage_mod,only: decoupage
  use gen_com_m,only:lspacendm
  use atomconfig,only: atom_config,atom_config_d
  USE boxconfig,only:box_config,periodbox
  USE cellconfig,only:cell_config,caltabtC
  use calfo_mod,only:calfo
  use constrconf_mod,only:repartition
   USE caltabi_mod,only: caltabi
  implicit none
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
    
    if ((div%mpi_image%nproc.gt.1).and.(lspaceNDM.eqv..true.)) then
       if (ldistr) then
          call atcomp%send2all(0,div%mpi_image%comm)
       endif
       call cellcomp%copy_cell(celloc)
       call decoupage(div%mpi_image%nproc,0,celloc,atloc,lverbose=.false.,psc=psc)
       call repartition(atcomp,atloc,box,celloc) ! mettre les éléments de la répartition dans un type
       call setcellconf(celloc,atloc,box,atcomp%im,rum,lverbose=.false.)
    else
       atloc=>atcomp
       celloc=>cellcomp
    end if

    call caltabtC(celloc,atloc,lperiod,box)

  end subroutine initloc
  
  subroutine initcomp(atcomp,cellcomp,atlocin,cellocin,box,div,lperiod)
    
    class(atom_config),intent(in)::atlocin
    type(cell_config),intent(in)::cellocin
    class(atom_config)::atcomp
    type(cell_config)::cellcomp
    type(box_config)::box
    type(para_config),intent(in)::div
    integer::ierr,iun
    logical::lperiod


       if ((div%mpi_image%nproc.gt.1).and.(lspaceNDM.eqv..true.)) then
!    if (div%mpi_image%nproc.gt.1) then
       call cellcomp%init(cellocin%nox,cellocin%noy,cellocin%noz,cellocin%natperc,cellocin%ltpcel)
       call atlocin%vers_master(atcomp,div)
       if (div%mpi_image%rank==0) then
          call caltabtC(cellcomp,atcomp,lperiod,box)
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
  
  subroutine pointer_caltabt_calfo(sig,potist,atcomp,cellcomp,box,atloc,celloc,div,lperiod,ltabvois,it,itetabvois,lchg,psc)

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
    integer::ierr,i
    logical::lchange=.true.
    integer::iun
    if(present(lchg))lchange=lchg
#ifdef PARA
    if (lchange) then
       if (div%mpi_image%nproc.gt.1)then

          if (lspaceNDM.eqv..true.) then
             call atcomp%master2loc(atloc,div)
          else
             call atcomp%send2all(0,div%mpi_image%comm)
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
    if (lperiod)   call periodbox (box,atloc)

    call caltabtC(celloc,atloc,lperiod,box)
    if (present(ltabvois)) then
       if (ltabvois.and.((it==1).or.(mod(it,itetabvois)==0)))&
            &call caltabi(atloc,celloc,box)
    end if

#ifdef PARA
!    if (lchange) then  ! même sans changement il faut mettre à jour pour initialisze les tableaux NDM like de mod_para pour maj_tab_density

    if ((div%mpi_image%nproc.gt.1).and.(lspaceNDM.eqv..true.)) then
       call maj_atomes_frt_ftm(atloc,celloc,psc)
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



    
end module parautils



