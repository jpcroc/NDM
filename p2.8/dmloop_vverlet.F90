module dmloop_vverlet_mod
  USE calfo_mod,only: calfo
  USE analyse_mod,only: analyse
  USE controle_mod,only: controle
  USE dyn_vverlet_mod,only: dyn_vverlet
!  USE calctemp_mod,only: calctemp
  USE sauvegarde_mod,only: sauvegarde
  USE sauveposition_mod,only: sauveposition
  USE sauveforce_mod,only: sauveforce
  USE correl_mod,only: correlvp
  USE atomconfig,only : atom_config_d,ndm2config, config2ndm
  USE cellconfig, only:cell_config,ndm2cellconfig,cellconfig2ndm
  USE boxconfig,only:box_config,boxconfig2ndm,ndm2boxconfig
  use var_pot,only:ntyp
  USE gen_com_m, ONLY: itesauvforce,itesauvposition,lcorrelvp,at,ecyl,ev2erg,im,lgc,rang,rayonc,&
       &tstep,vdc,pc,vdc,itdes,itesauv,itesigma,ldesinteg,lsigat,ltpcel,sigat,sigc&
       &,noxyz,lsuivinonpbc

  USE eloss, ONLY : calceloss,ibrake !, tcelec,ecelec,ibrake,elstopforce,elosselectot,elosselectot1,elosselec1,ngrdel,elosselec
  USE elec_cell, ONLY :i2t       
USE calfoberend_mod,only:calfoberend


  implicit none 
contains
  ! boucle de DM pour velocity Verlet
  ! ************************************************

  subroutine dmloop_vverlet
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE Parrinello_Rahman
    USE tab_imm_m,only:xp,xpp,vp,fp,iwmax,ityp,ielat,num_at_glob,ax
    USE suivinonpbc


#ifdef PARA
  use mpi
  USE mod_para,only:MPI_COMM_space,ierr,NDM_MPI_REAL_DOUBLE,status,nprocs,temps_debpara,temps_para

#endif
    implicit none
    character :: extension*2
    integer::lenfn2,i
    integer::ilocal
    real(double) sigkine_tot(3,3)
!    real(double) :: temptyp(ntyp)
#ifdef PARA
    ! declarations supplementaires pour MPI
    real(double), dimension(3,3,noxyz) :: sigc_tot

#endif

    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    ! MPI
    type(box_config)::boxndm
    type(atom_config_d)::atdml
    type(cell_config):: celndm
    logical :: test_sigma
    if (rang==0) write (6, *) '***** PREMIERE ITERATION  VVERLET****'
#ifdef PARA
    temps_para=0.
#endif

    if (lsuivinonpbc) call init_suivinonpbc()
    ! Appel de la routine generale des forces
    !    call calfo
  call ndm2boxconfig(at,bg,zl,zls2,nzl,volu,normat,boxndm)
      call ndm2cellconfig(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
  call ndm2config(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,&
       &iwmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
  test_sigma=(mod(it,itesigma)==0)
  CALL CalFo(sig,potist,atdml,celndm,boxndm,t_sigma=test_sigma)
  if (l2t)then
     if (i2t==1)  call calceloss (atdml%im,atdml%fp,atdml%vp,atdml%ityp,atdml%ielat,atdml%num_at_glob)
    else
       if(ibrake.gt.0) call calceloss (atdml%im,atdml%fp,atdml%vp,atdml%ityp,atdml%ielat,atdml%num_at_glob)
    end if
    if (lTberendsen) call calfoberend(atdml%im,atdml%imm,atdml%xp,atdml%vp,atdml%fp,atdml%ityp)
    call config2ndm(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,xpp=xpp)
    call cellconfig2ndm(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize) !inutile (calfo ne change pas celndm) mais laissé par sécurite




    !  call analyse
!    call calctemp (temptyp) 
1   continue
    it = it+1


    if (ldesinteg)itdes=itdes+1


    call dyn_vverlet
    ! les positions et les vitesses sont synchrones en ce point ; les atomes sont bien r�partis en cellules




    sigkine=0.
    do ilocal = 1, im
       sigkine(1:3,1) = sigkine(1:3,1) + &
            cm(ityp(ilocal))*vp(1:3,ilocal)*vp(1,ilocal)
       sigkine(1:3,2) = sigkine(1:3,2) + &
            cm(ityp(ilocal))*vp(1:3,ilocal)*vp(2,ilocal)
       sigkine(1:3,3) = sigkine(1:3,3) + &
            cm(ityp(ilocal))*vp(1:3,ilocal)*vp(3,ilocal)
       if (lsigat) then 
          sigat(1:3,1,ilocal) = sigat(1:3,1,ilocal) +  cm(ityp(ilocal))*vp(1:3,ilocal)*vp(1,ilocal)
          sigat(1:3,2,ilocal) = sigat(1:3,2,ilocal) +  cm(ityp(ilocal))*vp(1:3,ilocal)*vp(2,ilocal)
          sigat(1:3,3,ilocal) = sigat(1:3,3,ilocal) +  cm(ityp(ilocal))*vp(1:3,ilocal)*vp(3,ilocal)
       end if
       if ((mod(it,itesigma)==0).and.(lTPcel.EQV..true.)) then
          sigc(1:3,1,ielat(ilocal)) = sigc(1:3,1,ielat(ilocal)) + &
               cm(ityp(ilocal))*vp(1:3,ilocal)*vp(1,ilocal)*noxyz/volu
          sigc(1:3,2,ielat(ilocal)) = sigc(1:3,2,ielat(ilocal)) + &
               cm(ityp(ilocal))*vp(1:3,ilocal)*vp(2,ilocal)*noxyz/volu
          sigc(1:3,3,ielat(ilocal)) = sigc(1:3,3,ielat(ilocal)) + &
               cm(ityp(ilocal))*vp(1:3,ilocal)*vp(3,ilocal)*noxyz/volu
       end if
    end do
    sigkine(1:3,1:3) = sigkine(1:3,1:3)/volu
#ifdef PARA

    !  call MPI_ALLREDUCE(sig,sig_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
    !  sig=sig_tot
    call MPI_ALLREDUCE(sigkine,sigkine_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
    sigkine=sigkine_tot
    if (allocated(sigc)) then
       call MPI_ALLREDUCE(sigc,      sigc_tot,      9*noxyz,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
       sigc=sigc_tot
    end if
#endif
    sigtot = sigkine+sig

    call analyse 
    if (lcorrelvp) call correlvp(xp,xpp,vp,ax,fp,ityp)
    ! MPI
    !     write(6,*)'analyse -> sauvegarde'XS
    if (itesauv.GT.0) then
       if (mod(it,itesauv)==0) call sauvegarde
    endif

    !     write(6,*)'analyse -> sauveposition'
    if (itesauvposition.GT.0) then
       if (mod(it,itesauvposition)==0) then
          call sauveposition (it)
          if (lsuivinonpbc) then
             call reset_suivinonpbc
             call sauvepositionnonpbc (it)
          end if
       end if
    endif
    if (itesauvforce.GT.0) then
       if (mod(it,itesauvforce)==0) call sauveforce (it)
    end if
    !     write(6,*)'sauvposition -> control'

    call controle
    !     write(6,*)' controle ->'
    go to 1

    return
  end subroutine dmloop_vverlet
end module dmloop_vverlet_mod
