module dmloop_vverlet_mod
  USE calfo_mod,only: calfo
  USE analyseT_mod,only: analyseT
  USE controleT_mod,only: controleT
  USE dyn_vverlet_mod,only: dyn_vverlet
  USE atomconfig,only : atom_config_d, atom_config_e
  USE cellconfig, only:cell_config
  USE boxconfig,only:box_config
  use var_pot,only:ntyp
  USE gen_com_m, ONLY: itesauvforce,itesauvposition,ecyl,ev2erg,lgc,rang,rayonc,&
       &tstep,vdc,pc,vdc,itdes,itesauv,itesigma,ldesinteg,lsigat,ltpcel,lspaceNDM,itmax

  USE eloss, ONLY : calceloss,ibrake !, tcelec,ecelec,ibrake,elstopforce,elosselectot,elosselectot1,elosselec1,ngrdel,elosselec
  USE elec_cell, ONLY :i2t       
USE calfoberend_mod,only:calfoberend
use Tpara,only:para_space_config
use endrunT_mod,only:endrunT


  implicit none 
contains
  ! boucle de DM pour velocity Verlet
  ! ************************************************

  subroutine dmloop_vverlet(atdml,celndm,boxndm,psc)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE Parrinello_Rahman


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
    character :: extension*2
    integer::lenfn2,i
    integer::ilocal
    real(double) sigkine_tot(3,3)
    !    real(double) :: temptyp(ntyp)

    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    logical :: test_sigma
    if (rang==0) write (6, *) '***** PREMIERE ITERATION  VVERLET****'

    ! Appel de la routine generale des forces
    test_sigma=(mod(it,itesigma)==0)

    CALL CalFo(sig,potist,atdml,celndm,boxndm,t_sigma=test_sigma,psc=psc)
    if (l2t)then
       if (i2t==1)  call calceloss(celndm,atdml)
    else
       if(ibrake.gt.0) call calceloss(celndm,atdml)
    end if
    if (lTberendsen) call calfoberend(atdml)
    if (itmax==0) then
       call analyseT (atdml,celndm,boxndm)
       call endrunT(atdml,celndm,boxndm,.true.)
    end if


    !  call analyse
    !    call calctemp (temptyp) 
1   continue
    it = it+1


    if (ldesinteg)itdes=itdes+1
    call dyn_vverlet(atdml,celndm,boxndm,psc)
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

       
if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
   call comm_space%sum(sigkine)
       if (allocated(celndm%sigc)) then
          call comm_space%sum(celndm%sigc)
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
