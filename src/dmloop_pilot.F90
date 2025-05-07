module dmloop_pilot_mod
  USE arret_ndm_mod,only:arret_ndm
  use arps_mod,only:dmloop_arps
  USE atomconfig,only : atom_config_d, atom_config_e,atom_config_arps
  USE cellconfig, only:cell_config,cell_config_arps
  USE boxconfig,only:box_config,box_config_lpr
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
    class(box_config)::boxndm
    class(atom_config_d)::atdml
    class(cell_config):: celndm
    logical,optional::linit
    logical::lini=.false.
    if (present(linit))lini=linit
!    call celndm%print(unit=100)
!    stop
    select case (dmtype)
    case(41,42)
       select type(atdml)
       type is (atom_config_arps) !special case ARPS
          select type(celndm)
          type is(cell_config_arps)
             call dmloop_arps(atdml,celndm,boxndm,psc)
          end select
       end select

    case(4)
       call dmloop_vverlet (atdml,celndm,boxndm,psc)
    case(8,22,24,88)
       select type (boxndm)
       type is (box_config_lpr)
          call dmloop_lpr (atdml,celndm,boxndm,psc,linit=lini)
       type is (box_config)
          write(6,*)'WTF dmloop_pilot call lpr'
          call arret_ndm
       end select
                
    case (1,21,23)
       call dmloop (atdml,celndm,boxndm,psc)
    case default
       write(6,*)'DMTYPE ?? dmloop_pilot'
       call arret_ndm
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
