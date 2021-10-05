module dmloop_mod
  USE calfo_mod,only: calfo
  USE dyn_mod,only: dyn
  USE analyseT_mod,only: analyseT
  USE controleT_mod,only: controleT
  USE trempe_mod,only: trempe
  USE gen_com_m, ONLY:itesauvforce,itesauvposition,lfire
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE cellconfig, only:cell_config,caltabtC
  USE boxconfig,only:box_config
  USE eloss, ONLY : calceloss,ibrake 
  USE elec_cell, ONLY :i2t       
  USE calfoberend_mod,only:calfoberend
  USE caltabi_mod,only: caltabi
  USE parautils,only:driver_caltabt_DM

  USE gen_com_m,only: dmtype,it,itesauv, potist,rang,sig,l2t,sigkine,sigtot,itesigma,ltberendsen,itab, &
       & itetabvois,lperiod,lspaceNDM,itmax
  use var_pot, only: cm
   use Tpara,only:nprocspace,para_space_config,comm_space

  implicit none
contains
  ! ************************************************
  !           Sous-programme dmloop.f
  ! ************************************************

  subroutine dmloop (atdml,celndm,boxndm,psc)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE FireModule

    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    type(para_space_config)::psc
    integer :: i, iti,ilocal
    REAL(double) :: fire_dt, fire_alph
    INTEGER :: fire_nstep
    !-----------------------------------------------
    !
    class(atom_config_d)::atdml
    type(cell_config):: celndm
    type(box_config)::boxndm
    logical:: test_sigma=.false.


    if (rang==0) write (6, *) '***** PREMIERE ITERATION  VERLET STD ***'

    ! Initialization
    IF ((dmtype.EQ.21).AND.lFire) THEN
       CALL init_trempe_fire(fire_dt, fire_nstep, fire_alph)
    END IF

    !      write(6,*)'im',im
    do while (it.le.itmax)

       it = it+1


       ! appel de la routine generale des forces
       if (itesigma>0)      test_sigma=(mod(it,itesigma)==0)
       if (test_sigma) then
          sig(:,:)=0.d0 ; if (celndm%ltpcel.EQV..true.) celndm%sigc=0
          select type(atdml)
          type is (atom_config_e)
             if(atdml%lSigat) atdml%sigat(:,:,:)=0. ;
          end select
       end if
       CALL CalFo(sig,potist,atdml,celndm,boxndm,t_sigma=test_sigma,psc=psc)
!!$    if (l2t)then
!!$       if (i2t==1)  call calceloss (atdml%im,atdml%fp,atdml%vp,atdml%ityp,atdml%ielat,atdml%num_at_glob)
!!$    else
       if(ibrake.gt.0) call calceloss (celndm,atdml)
!!$    end if

       if (test_sigma) then                   
          sigkine=0.
          do ilocal = 1, atdml%im
             sigkine(1:3,1) = sigkine(1:3,1) + &
                  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(1,ilocal)
             sigkine(1:3,2) = sigkine(1:3,2) + &
                  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(2,ilocal)
             sigkine(1:3,3) = sigkine(1:3,3) + &
                  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(3,ilocal)
             select type(atdml)
             type is (atom_config_e)
                if(atdml%lSigat) atdml%sigat(:,:,:)=0. ;

                if (atdml%lsigat) then 
                   atdml%sigat(1:3,1,ilocal) = atdml%sigat(1:3,1,ilocal) +  &
                        &cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(1,ilocal)
                   atdml%sigat(1:3,2,ilocal) = atdml%sigat(1:3,2,ilocal) + &
                        &cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(2,ilocal)
                   atdml%sigat(1:3,3,ilocal) = atdml%sigat(1:3,3,ilocal) +  &
                        &cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(3,ilocal)
                end if
             end select
          end do

          sigkine(1:3,1:3) = sigkine(1:3,1:3)/boxndm%volu

#ifdef PARA
          if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
             call comm_space%sum(sigkine)
          end if
#endif

          sigtot = sigkine+sig
       end if
       if (lTberendsen) call calfoberend(atdml)

       select case (dmtype)

       case(1)
          call dyn  (atdml)



       case (21) 
          IF (lFire) THEN
             call trempe_fire (atdml,fire_dt, fire_nstep, fire_alph)

          ELSE
             call trempe (atdml)
          END IF

       case default
          write (6, *) 'ne sait pas quoi faire stop'
          stop
       end select

       call  driver_caltabt_DM(sig,potist,atdml,celndm,boxndm,psc,lperiod)

       call analyseT (atdml,celndm,boxndm)    
       call controleT(atdml,celndm,boxndm)

    end do


    return

  end subroutine dmloop
end module dmloop_mod
