module dmloop_lpr_mod
  USE analyseT_mod,only: analyseT
  USE controleT_mod,only: controleT
  USE sauveposition_mod,only: sauveposition
  USE sauveforce_mod,only: sauveforce
  USE gen_com_m, ONLY: itesauvforce, itesauvposition,itesauv,ltnose,lperiod,lspacendm
  USE calfo_mod,only: calfo

  USE atomconfig,only : atom_config_d
  USE cellconfig, only:cell_config,caltabtc
  USE boxconfig,only:box_config
  use Tpara,only:para_space_config


  USE eloss, ONLY : calceloss,ibrake !, tcelec,ecelec,ibrake,elstopforce,elosselectot,elosselectot1,elosselec1,ngrdel,elosselec
  USE elec_cell, ONLY :i2t
USE calfoberend_mod,only:calfoberend

#ifdef PARA
  USE recips_mod,only: recips
#endif
  implicit none
contains
  ! boucle de DM pour velocity Verlet
  ! ************************************************

  subroutine dmloop_lpr(atpr,celndm,boxndm,psc)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE Parrinello_Rahman
    USE Parrinello_Rahman_Nose

#ifdef PARA
  USE Tpara,only:nprocspace
  USE mod_para,only:maj_atomes_frt_ftm

#else
  use Tpara,only:nprocspace
#endif
    implicit none

    type(para_space_config)::psc
    type(box_config)::boxndm
    class(atom_config_d)::atpr
    type(cell_config):: celndm
#ifdef PARA
    real(double)::wbox_tot
    real(double) sigkine_tot(3,3)
    integer :: nb1, nb2, nb3, i1, l,noxn,noyn,nozn
    real(double) :: zlx, zly, zlz, ux, uy, uz,  pi2, fact, fact1&
         , fact2, hk2, ex, ex1, ex2
    !real(double), external :: calcvol

#endif

    if (rang==0) write (6, *) '***** PREMIERE ITERATION  ****'


    ! Initialization -------------------------------------------------------
    IF (lTNose) THEN ! Parrinello-Rahman with Nose thermostat
       call initlprNose(atpr,celndm,boxndm)
    ELSE ! Parinello-Rahman with Nose-Hoover thermostat or constant energy
       call initlpr(atpr,celndm,boxndm,psc)
    END IF

!    call analyseT(atpr,celndm,boxndm)

    ! MD loop -------------------------------------------------------------
1   continue
    it = it+1
    IF (lTNose) THEN ! Parrinello-Rahman with Nose thermostat
  CALL CalFo(sig,potist,atpr,celndm,boxndm,t_sigma=.true.,psc=psc)

!  CALL CalFo(sig,potist,atpr,celndm)
      if (l2t)then
       if (i2t==1)  call calceloss (atpr%im,atpr%fp,atpr%vp,atpr%ityp,atpr%ielat,atpr%num_at_glob)
    else
       if(ibrake.gt.0) call calceloss (atpr%im,atpr%fp,atpr%vp,atpr%ityp,atpr%ielat,atpr%num_at_glob)
    end if
    if (lTberendsen) call calfoberend(atpr%im,atpr%imm,atpr%xp,atpr%vp,atpr%fp,atpr%ityp)
!  write(6,*)'dml potist ',potist,atpr%potist



       call prNose(atpr,celndm,boxndm)

#ifdef PARA
if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then

       boxndm%zl(1) = Sqrt( Sum(boxndm%at(1:3,1)**2 ) )
       boxndm%zl(2) = Sqrt( Sum(boxndm%at(1:3,2)**2 ) )
       boxndm%zl(3) = Sqrt( Sum(boxndm%at(1:3,3)**2 ) )
       boxndm%volu=calcvol(boxndm%at(1:3,1),boxndm%at(1:3,2),boxndm%at(1:3,3))
       boxndm%zls2(1:3) = 0.5d0*boxndm%zl(1:3)
             call caltabtC(celndm,atpr,lperiod,boxndm)

       ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs

       call maj_atomes_frt_ftm(atpr,celndm,psc)


       if (iewald>0) then

          ! --- Tableaux des troisiemes termes de la sommation d'Ewald ---
          auxe = 23.06134575D-20                  ! en erg.cm (charge electron^2/4*pi*permitivite vide)
          pi2 = pi*pi
          boxndm%volu=calcvol(boxndm%at(1:3,1),boxndm%at(1:3,2),boxndm%at(1:3,3))
          fact = pi2/alpha**2
          fact1 = auxe/2./pi/boxndm%volu
          fact2 = auxe*2./boxndm%volu
          do nb1 = -ncoucx, ncoucx
             do nb2 = -ncoucy, ncoucy
                do nb3 = -ncoucz, ncoucz
                   if (nb1==0.and.nb2==0.and.nb3==0) cycle
                   hk2 = nb1*nb1/boxndm%zl(1)**2+nb2*nb2/boxndm%zl(2)**2+nb3*nb3/boxndm%zl(3)**2
                   ex = exp((-hk2*fact))/hk2
                   ex1 = ex*fact1
                   ex2 = ex*fact2
                   tabv3(nb1,nb2,nb3) = ex1
                   tabf3(:,nb1,nb2,nb3) = ex2*q(:)
                end do
             end do
          end do

       endif
    else
           CALL ScaleBox(atpr,celndm,boxndm)

    end if
#else

       CALL ScaleBox(atpr,celndm,boxndm)

#endif
       timel=timel+fNose*tstep
    ELSE ! Parinello-Rahman with Nose-Hoover thermostat or constant energy
       call pr(atpr,celndm,boxndm,psc)
       timel=timel+tstep
    END IF

    call analyseT(atpr,celndm,boxndm)


!!$    if (itesauv.GT.0) then
!!$       if (mod(it,itesauv)==0) call sauvegarde
!!$    endif
!!$    if (itesauvposition.GT.0) then
!!$       if (mod(it,itesauvposition)==0) call sauveposition (it)
!!$    endif
!!$    if (itesauvforce.GT.0) then
!!$       if (mod(it,itesauvforce)==0) call sauveforce (it)
!!$    endif


     call controleT(atpr,celndm,boxndm)

    go to 1

    return
  end subroutine dmloop_lpr

end module dmloop_lpr_mod
