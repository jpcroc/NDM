module prog_mod
  USE init_mod,only: init
  USE calfo_mod,only: calfo
  USE analyse_mod,only: analyse
  USE controle_mod,only: controle
  USE endrunT_mod,only: endrunT
  USE neb_mod,only: neb
  USE dmloop_lpr_mod,only: dmloop_lpr
  USE loopforcetest_mod,only: loopforcetest
  USE gcII_mod,only: gcII
  USE dmloop_vverlet_mod,only: dmloop_vverlet
  USE dmloop_mod,only: dmloop
  USE analyseT_mod,only: analyseT
  USE controleT_mod,only: controleT
 
  USE montecarlo_mod, only: montecarlo
  USE boxconfig,only:box_config,boxconfig2ndm,ndm2boxconfig
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e,ndm2config, config2ndm
  USE cellconfig, only:cell_config,ndm2cellconfig,cellconfig2ndm
#if defined ML || defined PARAML    
  USE ml_main_mod,only: ml_main
#endif
      USE cellconfig, only:cell_config,ndm2cellconfig,cellconfig2ndm
  implicit none
contains
  subroutine prog
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:dmtype,im,imm,indi,ltabvois,parallele,potist,rang,sig,nvois ,&
         &nox,noy,noz,noxyz,natperc,nato,ncel,atincel,deltadist,celsize,&
             &at,bg,zl,zls2,nzl,volu,normat,lax,lprteat,lsigat

    USE tab_imm_m
    
!    USE montecarlo_mod, ONLY: config_atom_n, cells_n

#ifdef PARA
    use mpi
    USE mod_para,only:MPI_COMM_space,TEMPS_INIT_DEB,TEMPS_INIT,MYID,TEMPS_DEB,TEMPS_DMLOOP_DEB,maj_atomes_frt_ftm
#endif

    implicit none
    character :: extension*2
    integer::lenfn2,i,ko
    class(atom_config_d),pointer::atdml
!    type(atom_config),target:: atdm
    type(atom_config_d),target:: atdmd
    type(atom_config_e),target:: atdme
    type(cell_config)::celndm
    type(box_config)::boxndm


    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e rs
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------

    ! Allocation des tableaux dimensionnes sur le nombre d'atomes
    call alloc_all_tab_imm(imm)
    if ((lax).or.(lsigat).or.(lprteat).or.(llangevin))then
       atdml=>atdme
    else
       atdml=>atdmd
    end if
    im=0 ; nvois=0
    call atdml%init(im,imm,ltabvois,nvois,lsigat,lprteat,llangevin,lax)
    
#ifdef PARA
    temps_init_deb = MPI_Wtime()
#endif
    ! Initilisation
    call init(atdml,boxndm,celndm)
#ifdef PARA
    temps_init=MPI_Wtime()-temps_init_deb
#endif

#ifdef DECOUP
    ! Dans ce cas, pas la peine d'aller plus loin on peut terminer le programme
    return
#endif

#ifdef PARA
    ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs
    if (rang==0) write(6,*)'PARA-T avant MAJ'


    call maj_atomes_frt_ftm
    if (rang==0) write(6,*)'PARA-T apres MAJ'




    ! Affichage du temps d'initialisation
#ifdef PARA
    if (myid==0) then
       print *, 'Temps d''initialisation : ', MPI_Wtime() - temps_deb
    endif
#endif

#endif

#ifdef PARA
    temps_dmloop_deb = MPI_Wtime()
#endif
    ! Actuellement uniquement le cas dmloop_vverlet est traite en parallele


    select case (dmtype) 
    case(5)
       if (.not.parallele)    call loopforcetest (xp, xpp, vp, ax, fp, ielat, iwmax, ityp,num_at_glob)
    case(4,10)
       call dmloop_vverlet (atdml,celndm,boxndm)
    case(8)
       call dmloop_lpr 
    case (1)
       if (.not.parallele)  call dmloop (atdml,celndm,boxndm)
    case (2)
       if (.not.parallele)  then
          call dmloop(atdml,celndm,boxndm)
       else
#ifdef PARA
          if(rang==0) write (6,*)'DMTYPE 2 +PARA=DMLOOP_VVERLET_+OPTION'
          call dmloop_vverlet ! (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
#endif
       endif
    case (3,30)
       !    stop

       call gcII ! (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
    case (9)
       if (.not.parallele)   call neb  ! (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
    case(11)
       if (rang==0) write (6, *) '***** PREMIERE ET UNIQUE ITERATION  ****'
!       call ndm2boxconfig(at,bg,zl,zls2,nzl,volu,normat,boxndm)
!       call ndm2cellconfig(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
!       call ndm2config(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,&
!            &iwmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
       CALL CalFo(sig,potist,atdml,celndm,boxndm) !(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
!       call config2ndm(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,xpp=xpp)
!    call cellconfig2ndm(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize) !inutile (calfo ne change pas celndm) mais laissé par sécurite
       call analyseT(atdml,celndm,boxndm)
       call controleT(atdml,celndm,boxndm)
       call endrunT(atdml,celndm,boxndm)

#ifdef ART    
    case (12) 
       call art90
#endif

#ifdef SUNDAE    
    case (16) 
       call sundae
#endif

#ifdef MAB    
    case (17) 
       call mab
#endif

#if defined PHONDY || defined PARAPH    
    case (7) 
       call phondy
#endif

#if defined ML || defined PARAML    
    case (18) 
       call ml
#endif


    case (15)
!       call ndm2cellconfig(cells_n,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
!       call ndm2config(config_atom_n,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,i&
!         &wmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
       call montecarlo(atdml,celndm,boxndm)

    end select


  end subroutine prog
end module
