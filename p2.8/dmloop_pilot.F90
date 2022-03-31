module dmloop_pilot_mod
  USE arret_ndm_mod,only:arret_ndm
  USE atomconfig,only : atom_config_d, atom_config_e
  USE cellconfig, only:cell_config
  USE boxconfig,only:box_config
  USE gen_com_m, ONLY: dmtype,lcdp,rang,latcomp

  use Tpara,only:para_space_config
  use endrunT_mod,only:endrunT

  USE dmloop_vverlet_mod,only: dmloop_vverlet
  USE dmloop_mod,only: dmloop
  USE dmloop_lpr_mod,only: dmloop_lpr
 USE arret_ndm_mod,only: arret_ndm
  implicit none 
contains
  ! boucle de DM pour velocity Verlet
  ! ************************************************

  subroutine dmloop_pilot(atdml,celndm,boxndm,psc,linit)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    implicit none
    type(para_space_config)::psc
    type(box_config)::boxndm
    class(atom_config_d)::atdml
    type(cell_config):: celndm
    logical,optional::linit
    logical::lini=.false.
    if (present(linit))lini=linit

    select case (dmtype) 
    case(4,10)
       call dmloop_vverlet (atdml,celndm,boxndm,psc)
    case(8)
       call dmloop_lpr (atdml,celndm,boxndm,psc,linit=lini)
    case (1)
       call dmloop (atdml,celndm,boxndm,psc)
    case (21)
       call dmloop(atdml,celndm,boxndm,psc)
    case(22)
       call dmloop_vverlet (atdml,celndm,boxndm,psc)
    end select

    if (lcdp) then
       if (rang==0) write (6, *) '*******return CDP**** '
       return
    else
       if (rang==0) write (6, *) '*******Derniere iteration **** '
       call endrunT(atdml,celndm,boxndm,latcomp)
       write (6, *) 'predeal '
       !       call DeallocateAll
       call arret_ndm
    end if

    return
  end subroutine dmloop_pilot
end module dmloop_pilot_mod
