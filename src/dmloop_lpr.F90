module dmloop_lpr_mod
  USE analyseT_mod,only: analyseT
  USE controleT_mod,only: controleT
  USE gen_com_m, ONLY: itesauvforce, itesauvposition,itesauv,ltnose,lperiod,lspacendm,itloopmax,pi,l2t,&
       &ltberendsen,potist,iteration,tstep,sig,rang,timel,timeloopmax,rang
  USE calfo_mod,only: calfo

  USE atomconfig,only : atom_config_d
  USE cellconfig, only:cell_config,caltabtc
  USE boxconfig,only:box_config,box_config_lpr
  use Tpara,only:para_space_config

  USE eloss, ONLY : calceloss,ibrake !, tcelec,ecelec,ibrake,elstopforce,elosselectot,elosselectot1,elosselec1,ngrdel,elosselec
  USE elec_cell, ONLY :i2t
  USE calfoberend_mod,only:calfoberend
  implicit none
contains
  ! boucle de DM pour velocity Verlet
  ! ************************************************

  subroutine dmloop_lpr(atpr,celndm,boxndm,psc,linit)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE Parrinello_Rahman,only:pr1,initlpr
    USE Parrinello_Rahman_Nose,only:prnose,initlprnose

#ifdef PARA

  USE mod_para,only:maj_atomes_frt_ftm

#else
  use Tpara,only:nprocspace
#endif
    implicit none

    type(para_space_config)::psc
    type(box_config_lpr)::boxndm
    class(atom_config_d)::atpr
    type(cell_config):: celndm
    logical,optional::linit
    logical::lini=.false.
    logical:: lreturn
!    integer::icall
!    icall=0
    if (present(linit))lini=linit
    

    if (rang==0) write (6, *) '***** PREMIERE ITERATION LPR  ****', itloopmax,timeloopmax

    if(lini) then 
       ! Initialization -------------------------------------------------------
       IF (lTNose) THEN ! Parrinello-Rahman with Nose thermostat
          call initlprNose(atpr,celndm,boxndm%box_config)
       ELSE ! Parinello-Rahman with Nose-Hoover thermostat or constant energy
          call initlpr(atpr,celndm,boxndm,psc)
       END IF
    end if
    ! MD loop -------------------------------------------------------------

    do while ((iteration.le.itloopmax).and.(timel.lt.timeloopmax))
    iteration = iteration+1
    IF (lTNose) THEN ! Parrinello-Rahman with Nose thermostat
       CALL CalFo(sig,potist,atpr,celndm,boxndm%box_config,t_sigma=.true.,psc=psc)
       !  CALL CalFo(sig,potist,atpr,celndm)
       if (l2t)then
          if (i2t==1)  call calceloss (celndm,atpr)
       else
          if(ibrake.gt.0) call calceloss(celndm,atpr)
       end if
       if (lTberendsen) call calfoberend(atpr)
       !  write(6,*)'dml potist ',potist,atpr%potist
       call prNose(atpr,celndm,boxndm%box_config,psc)
    ELSE ! Parinello-Rahman with Nose-Hoover thermostat or constant energy
!       if (icall==0)call analyseT(atpr,celndm,boxndm%box_config,psc)
!       icall=icall+1
       call pr1(atpr,celndm,boxndm,psc)
       timel=timel+tstep
    END IF
    call analyseT(atpr,celndm,boxndm%box_config,psc)
     call controleT(atpr,celndm,boxndm%box_config,psc,lreturn)
     if (lreturn) return

  end do

    return
  end subroutine dmloop_lpr

end module dmloop_lpr_mod
