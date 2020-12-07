module dmloop_mod
  USE calfo_mod,only: calfo
  USE dyn_mod,only: dyn
  USE analyseT_mod,only: analyseT
  USE controleT_mod,only: controleT
  USE trempe_mod,only: trempe
  USE correl_mod,only: correlvp
  USE gen_com_m, ONLY:itesauvforce,itesauvposition,lcorrelvp,lfire
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE cellconfig, only:cell_config,caltabtC
  USE boxconfig,only:box_config,periodbox
  USE eloss, ONLY : calceloss,ibrake 
  USE elec_cell, ONLY :i2t       
  USE calfoberend_mod,only:calfoberend
  USE caltabi_mod,only: caltabi

  USE gen_com_m,only: dmtype,it,itesauv, potist,rang,sig,l2t,sigkine,sigtot,itesigma,ltberendsen,itab, &
       & itetabvois,lperiod 
  use var_pot, only: cm! iewald,l3c,npotmax,potiseam,lpotentiel,cm,ipotentiel,potisglue,potisrep,potiseam
#ifdef PARA
  use mpi
  use mod_para,only:NDM_MPI_REAL_DOUBLE,MPI_COMM_space,ierr,nprocspace
#else
  
#endif
  
  implicit none
contains
  ! ************************************************
  !           Sous-programme dmloop.f
  !          Version MPI du 21 fevrier 2001
  ! ************************************************

  subroutine dmloop (atdml,celndm,boxndm)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    !  USE tab_imm_m,only:xp,xpp,vp,fp,iwmax,ityp,ielat,num_at_glob,ax
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
    integer :: i, iti,ilocal
    REAL(double) :: fire_dt, fire_alph
    INTEGER :: fire_nstep
    !-----------------------------------------------
    !
    class(atom_config_d)::atdml
    type(cell_config):: celndm
    type(box_config)::boxndm
#ifdef PARA
    real(double), dimension(3,3) :: sig_tot,sigkine_tot

#endif
    logical:: test_sigma=.false.


    ! MPI
    if (rang==0) write (6, *) '***** PREMIERE ITERATION  VERLET STD ***'

    ! Initialization
    IF ((dmtype.EQ.2).AND.lFire) THEN
       CALL init_trempe_fire(fire_dt, fire_nstep, fire_alph)
    END IF

    !      write(6,*)'im',im
1   continue
    it = it+1

    write(6,*)
    write(6,*)'***** ITERATION  ****', it

    ! appel de la routine generale des forces
    if (itesigma>0)      test_sigma=(mod(it,itesigma)==0)
    if (test_sigma) then
       sig(:,:)=0.d0 ; if (celndm%ltpcel.EQV..true.) celndm%sigc=0
       select type(atdml)
       type is (atom_config_e)
          if(atdml%lSigat) atdml%sigat(:,:,:)=0. ;
       end select
    end if
    CALL CalFo(sig,potist,atdml,celndm,boxndm,t_sigma=test_sigma)
    if (l2t)then
       if (i2t==1)  call calceloss (atdml%im,atdml%fp,atdml%vp,atdml%ityp,atdml%ielat,atdml%num_at_glob)
    else
       if(ibrake.gt.0) call calceloss (atdml%im,atdml%fp,atdml%vp,atdml%ityp,atdml%ielat,atdml%num_at_glob)
    end if

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
    if (nprocspace.gt.1) then

       !  call MPI_ALLREDUCE(sig,sig_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
       !  sig=sig_tot
       call MPI_ALLREDUCE(sigkine,sigkine_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
       sigkine=sigkine_tot
    end if
#endif

       sigtot = sigkine+sig
    end if
    if (lTberendsen) call calfoberend(atdml%im,atdml%imm,atdml%xp,atdml%vp,atdml%fp,atdml%ityp)

    select case (dmtype)

    case(1)
       call dyn  (atdml)

      
!       if (lcorrelvp) call correlvp(atdml%xp,atdml%xpp,atdml%vp,atdml%ax,atdml%fp,atdml%ax, atdml%ityp)

    case (2) 
       IF (lFire) THEN
          call trempe_fire (atdml,fire_dt, fire_nstep, fire_alph)

       ELSE
          call trempe (atdml)
       END IF

    case default
       write (6, *) 'ne sait pas quoi faire stop'
       stop
    end select
    if (lperiod)       call periodbox(boxndm,atdml)

    call analyseT (atdml,celndm,boxndm)


       if (itab/=0) then
          if (mod(it,itab)==0) then
             call caltabtC(celndm,atdml,lperiod,boxndm)
          endif
       endif

    if (atdml%ltabvois.and.mod(it,itetabvois)==0) then
       call caltabi(atdml%atom_config,celndm,boxndm)
    end if

    
    call controleT(atdml,celndm,boxndm)


    !     write(6,*)' controle ->'


    go to 1


    return

  end subroutine dmloop
end module dmloop_mod
