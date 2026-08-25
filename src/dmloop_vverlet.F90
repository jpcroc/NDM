module dmloop_vverlet_mod
  USE calfo_mod,only: calfo
  USE analyseT_mod,only: analyseT
  USE controleT_mod,only: controleT
  USE dyn_vverlet_mod,only: dyn_vverlet
  USE atomconfig,only : atom_config_d, atom_config_e
  USE cellconfig, only:cell_config
  USE boxconfig,only:box_config
  use var_pot,only:ntyp,cm
  USE gen_com_m, only:uwrt,lwrt, itesauvforce,itesauvposition,ev2erg,rang,iteration,l2t,lTberendsen,potist,sig,sigtot,&
       &tstep,itesauv,itesigma,lsigat,ltpcel,lspaceNDM,itloopmax,sigkine,timeloopmax,timel,lpcube,lcdp

  USE eloss, ONLY : calceloss,ibrake !, tcelec,ecelec,ibrake,elstopforce,elosselectot,elosselectot1,elosselec1,ngrdel,elosselec
  USE elec_cell, ONLY :i2t       
  USE calfoberend_mod,only:calfoberend
  use Tpara,only:para_space_config
  use endrunT_mod,only:endrunT
  use sigkinetot_mod,only:sigkinetot

  implicit none 
contains
  ! boucle de DM pour velocity Verlet
  ! ************************************************

  subroutine dmloop_vverlet(atdml,celndm,boxndm,psc)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

#ifdef PARA
    USE Tpara,only:COMM_space,nprocspace
#else
    USE Tpara,only:nprocspace
#endif
    implicit none
    type(para_space_config)::psc
    type(box_config)::boxndm
    class(atom_config_d)::atdml
    type(cell_config):: celndm
    integer::ic
    integer::ilocal
    real(double) ::pint
    !    real(double) :: temptyp(ntyp)

    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    logical :: test_sigma,lreturn
    if (rang==0) write (uwrt, *) '***** PREMIERE ITERATION  VVERLET****',itloopmax,timeloopmax
    ! Appel de la routine generale des forces
    test_sigma=((mod(iteration,itesigma)==0).or.(iteration==0))

    CALL CalFo(sig,potist,atdml,celndm,boxndm,t_sigma=test_sigma,psc=psc)
    if (l2t)then
       if (i2t==1)  call calceloss(celndm,atdml)
    else
       if(ibrake.gt.0) call calceloss(celndm,atdml)
    end if
    if (lTberendsen) call calfoberend(atdml)
    if (itloopmax==0) then
       call sigkinetot(atdml,boxndm,sig,sigkine,sigtot,celndm)
       call analyseT (atdml,celndm,boxndm,psc)
       if (lcdp) return
       call endrunT(atdml,celndm,boxndm,.true.)
    end if
    call analyseT (atdml,celndm,boxndm,psc)
     do while ((iteration.le.itloopmax).and.(timel.lt.timeloopmax))
       iteration = iteration+1
       test_sigma=(mod(iteration,itesigma)==0)
       call dyn_vverlet(atdml,celndm,boxndm,psc)
       ! les positions et les vitesses sont synchrones en ce point ; les atomes sont bien r�partis en cellules

       if (test_sigma) then
          call sigkinetot(atdml,boxndm,sig,sigkine,sigtot,celndm)
       end if
       call analyseT (atdml,celndm,boxndm,psc)
       call controleT(atdml,celndm,boxndm,psc,lreturn)
       if (lreturn) return

    end do

    return
  end subroutine dmloop_vverlet
end module dmloop_vverlet_mod
