module dmloop_vverlet_mod
  USE calfo_mod,only: calfo
  USE analyseT_mod,only: analyseT
  USE controleT_mod,only: controleT
  USE dyn_vverlet_mod,only: dyn_vverlet
  USE correl_mod,only: correlvp
  USE atomconfig,only : atom_config_d, atom_config_e
  USE cellconfig, only:cell_config
  USE boxconfig,only:box_config
  use var_pot,only:ntyp
  USE gen_com_m, ONLY: itesauvforce,itesauvposition,lcorrelvp,ecyl,ev2erg,lgc,rang,rayonc,&
       &tstep,vdc,pc,vdc,itdes,itesauv,itesigma,ldesinteg,lsigat,ltpcel,lsuivinonpbc,lspaceNDM

  USE eloss, ONLY : calceloss,ibrake !, tcelec,ecelec,ibrake,elstopforce,elosselectot,elosselectot1,elosselec1,ngrdel,elosselec
  USE elec_cell, ONLY :i2t       
USE calfoberend_mod,only:calfoberend


  implicit none 
contains
  ! boucle de DM pour velocity Verlet
  ! ************************************************

  subroutine dmloop_vverlet(atdml,celndm,boxndm)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE Parrinello_Rahman
    USE suivinonpbc


#ifdef PARA
    use mpi
    USE Tpara,only:MPI_COMM_space,NDM_MPI_REAL_DOUBLE
    USE mod_para,only:temps_debpara,temps_para
#else
  USE Tpara,only:nprocspace
#endif
    implicit none
    type(box_config)::boxndm
    class(atom_config_d)::atdml
    type(cell_config):: celndm
    character :: extension*2
    integer::lenfn2,i
    integer::ilocal
    real(double) sigkine_tot(3,3)
    !    real(double) :: temptyp(ntyp)
#ifdef PARA
    ! declarations supplementaires pour MPI
    real(double), dimension(3,3,celndm%noxyz) :: sigc_tot
    integer::ierr2
#endif

    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    ! MPI
    logical :: test_sigma
    if (rang==0) write (6, *) '***** PREMIERE ITERATION  VVERLET****'
#ifdef PARA
    temps_para=0.
#endif

    if (lsuivinonpbc) call init_suivinonpbc()
    ! Appel de la routine generale des forces
    test_sigma=(mod(it,itesigma)==0)

    CALL CalFo(sig,potist,atdml,celndm,boxndm,t_sigma=test_sigma)
    if (l2t)then
       if (i2t==1)  call calceloss (atdml%im,atdml%fp,atdml%vp,atdml%ityp,atdml%ielat,atdml%num_at_glob)
    else
       if(ibrake.gt.0) call calceloss (atdml%im,atdml%fp,atdml%vp,atdml%ityp,atdml%ielat,atdml%num_at_glob)
    end if
    if (lTberendsen) call calfoberend(atdml%im,atdml%imm,atdml%xp,atdml%vp,atdml%fp,atdml%ityp)




    !  call analyse
    !    call calctemp (temptyp) 
1   continue
    it = it+1


    if (ldesinteg)itdes=itdes+1
    call dyn_vverlet(atdml,celndm,boxndm)
    ! les positions et les vitesses sont synchrones en ce point ; les atomes sont bien r�partis en cellules

    if (test_sigma) then

       sigkine=0.
       do ilocal = 1, atdml%im
          sigkine(1:3,1) = sigkine(1:3,1) + &
               cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(1,ilocal)
          sigkine(1:3,2) = sigkine(1:3,2) + &
               cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(2,ilocal)
          sigkine(1:3,3) = sigkine(1:3,3) + &
               cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(3,ilocal)
          select type (atdml)
          type is (atom_config_e)

             if (atdml%lsigat) then 
                atdml%sigat(1:3,1,ilocal) = atdml%sigat(1:3,1,ilocal) +  &
                     &cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(1,ilocal)
                atdml%sigat(1:3,2,ilocal) = atdml%sigat(1:3,2,ilocal) + &
                     &cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(2,ilocal)
                atdml%sigat(1:3,3,ilocal) = atdml%sigat(1:3,3,ilocal) +  &
                     &cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(3,ilocal)
             end if
          end select
          if ((mod(it,itesigma)==0).and.(lTPcel.EQV..true.)) then
             celndm%sigc(1:3,1,atdml%ielat(ilocal)) = celndm%sigc(1:3,1,atdml%ielat(ilocal)) + &
                  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(1,ilocal)*celndm%noxyz/boxndm%volu
             celndm%sigc(1:3,2,atdml%ielat(ilocal)) = celndm%sigc(1:3,2,atdml%ielat(ilocal)) + &
                  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(2,ilocal)*celndm%noxyz/boxndm%volu
             celndm%sigc(1:3,3,atdml%ielat(ilocal)) = celndm%sigc(1:3,3,atdml%ielat(ilocal)) + &
                  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(3,ilocal)*celndm%noxyz/boxndm%volu
          end if
       end do
       sigkine(1:3,1:3) = sigkine(1:3,1:3)/boxndm%volu

#ifdef PARA

       !  call MPI_ALLREDUCE(sig,sig_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
       !  sig=sig_tot
       
if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       call MPI_ALLREDUCE(sigkine,sigkine_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr2)
       sigkine=sigkine_tot
       if (allocated(celndm%sigc)) then
          call MPI_ALLREDUCE(celndm%sigc,sigc_tot, 9*celndm%noxyz,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr2)
          celndm%sigc=sigc_tot
       end if
    end if
#endif
       sigtot = sigkine+sig
    end if
    call analyseT (atdml,celndm,boxndm)
    call controleT(atdml,celndm,boxndm)

    go to 1

    return
  end subroutine dmloop_vverlet
end module dmloop_vverlet_mod
