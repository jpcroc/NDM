module dyn_vverlet_mod
  USE calfo_mod,only: calfo
  USE calfoberend_mod,only: calfoberend 
  use var_pot,only:ntyp
  USE gen_com_m, ONLY:ilangevin,itab,dmtype,fnemd,lcalcjq,lnemd,lperiod,lprahman,&
       &l2T,llangevin,itesigma,it,itetabvois,ltberendsen,potist,sig,timel,tstep,&
       lspaceNDM
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e!,ndm2config, config2ndm
  USE cellconfig, only:cell_config,caltabtC
  USE boxconfig,only:box_config,periodbox
  use var_pot,only : cm
  USE eloss, only:ibrake, calceloss
  use Tpara,only:para_space_config
  USE parautils,only:driver_caltabt_DM
  implicit none
contains
  ! *************************************************************
  subroutine dyn_vverlet(atdml,celndm,boxndm,psc)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE jqmod
    USE elec_cell,ONLY: dynelec,i2t
#ifdef PARA
    USE mod_para,only:nprocspace,maj_atomes_frt_ftm
#else
    USE Tpara,only:nprocspace
#endif
    USE caltabi_mod,only:caltabi
    USE elec_cell, ONLY:TTlangevin
    !    USE Parrinello_Rahman
    use calfoberend_mod,only:dynlangevin
    use period_mod,only:period
    implicit none
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    type(para_space_config)::psc
    type(box_config)::boxndm
    class(atom_config_d)::atdml
    type(cell_config):: celndm
    integer :: i, ia,ic,il
    real(double), dimension(ntyp) :: aux
    real(double), save :: tmoyinst, imesureT
    !-----------------------------------------------
    !real(double), external :: tempinst
    real(double)::eatommoy


    logical::test_sigma=.false.

    timel = timel+tstep
    
    aux(:ntyp) = tstep/cm(:ntyp)/2.d0


    ! step 1 First half-step velocities update, v(t) -> v(t+dt/2)
    


    if (dmtype==22) then
       do i = 1, atdml%im
          do ic = 1, 3
             if (atdml%vp(ic,i)*atdml%fp(ic,i)<0) then
                atdml%vp(ic,i)=0.
             end if
          end do
       end do
    end if

    select type (atdml)
    class is (atom_config_e)
       if (lLangevin) then
          il=2*(ilangevin-1)+1
          call dynlangevin(atdml,il)
       elseif (l2T) then
          il=2*(ilangevin-1)+1
          call TTlangevin(atdml,il,psc,celndm)
       else
          DO i=1, atdml%im
             atdml%vp(1:3,i) = atdml%vp(1:3,i) + aux(atdml%iTyp(i))*atdml%fp(1:3,i)
          END DO
       end if
    class default
       DO i=1, atdml%im
          atdml%vp(1:3,i) = atdml%vp(1:3,i) + aux(atdml%iTyp(i))*atdml%fp(1:3,i)
       END DO
    end select


    !step 2  Coordinate update, x(t)-> x(t+dt)

    DO i=1, atdml%im
       atdml%xpp(1:3,i)=atdml%xp(1:3,i)
       atdml%xp(1:3,i) = atdml%xp(1:3,i) + tstep*atdml%vp(1:3,i)
    END DO


!!$    !conditions periodiques
!!$    if (lperiod)  call periodbox (boxndm,atdml)
!!$
!!$    ! repartition des atomes dans la nouvelle boite
!!$    if (.not.lprahman) then
!!$       if (itab/=0) then
!!$          if (mod(it,itab)==0) then
!!$             call caltabtC(celndm,atdml,lperiod,boxndm)
!!$          endif
!!$       endif
!!$    end if
!!$    if (atdml%ltabvois.and.mod(it,itetabvois)==0) then
!!$       call caltabi(atdml%atom_config,celndm,boxndm)
!!$    end if
!!$
!!$
!!$#ifdef PARA
!!$if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
!!$       ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs
!!$       call maj_atomes_frt_ftm(atdml,celndm,psc)
!!$    end if
!!$#endif
!!$
!!$
!!$    ! Force calculation
!!$
    call  driver_caltabt_DM(sig,potist,atdml,celndm,boxndm,psc,lperiod)

  ! a été déplacé après calfo . Etait situé juste avant calfo :
    if (l2T) then
       call dynelec(celndm)
    end if
  ! Force calculation
    jq=0.0
    if (itesigma>0) test_sigma=(mod(it,itesigma)==0)
    CALL CalFo(sig,potist,atdml,celndm,boxndm,t_sigma=test_sigma,psc=psc)

    
    if (l2t)then
       if (i2t==1)  call calceloss(celndm,atdml)
    else
       if(ibrake.gt.0) call calceloss(celndm,atdml)
    end if
    if (lTberendsen) call calfoberend(atdml%im,atdml%imm,atdml%xp,atdml%vp,atdml%fp,atdml%ityp)
    !  write(6,*)'dml potist ',potist,atdml%potist


    !    if (lnemd) then
    !       eatommoy=0.
    !       do i=1,imd
    !          eatommoy=eatommoy+eatom(i)/float(imd)
    !       end do
    !       do i=1,imd
    !          atdml%fp(1,i)=atdml%fp(1,i)+(eatom(i)-eatommoy)*Fnemd
    !       end do
    !    end if




    ! Second half-step velocities update, v(t+1/2dt) -> v(t+dt)
    select type (atdml)
    class is (atom_config_e)
       if (lLangevin) then
          il=2*(ilangevin-1)+2
          call dynlangevin(atdml,il)
       elseif (l2T) then
          il=2*(ilangevin-1)+3
          call TTlangevin(atdml,il,psc,celndm)
       else
          DO i=1, atdml%im
             atdml%vp(1:3,i) = atdml%vp(1:3,i) + aux(atdml%iTyp(i))*atdml%fp(1:3,i)
          END DO
       end if
    class default
       DO i=1, atdml%im
          atdml%vp(1:3,i) = atdml%vp(1:3,i) + aux(atdml%iTyp(i))*atdml%fp(1:3,i)
       END DO
    end select



    return
  end subroutine dyn_vverlet
end module dyn_vverlet_mod
