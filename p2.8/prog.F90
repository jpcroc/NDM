module prog_mod
        USE init_mod
        USE calfo_mod
        USE analyse_mod
        USE controle_mod
        USE endrun_mod
        USE neb_mod
        USE dmloop_lpr_mod
        USE loopforcetest_mod
        USE gcII_mod
        USE dmloop_vverlet_mod
        USE dmloop_mod
        USE atomconfig
#if defined ML || defined PARAML    
        USE ml_main_mod
#endif 
        implicit none
        contains
subroutine prog
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:
  USE tab_imm_m

#ifdef PARA
  USE mod_para
#endif

  implicit none
             character :: extension*2
    integer::lenfn2,i,ko
    type(atom_config_d)::atdml
    integer, allocatable ::iwmaxCF(:),indiCF(:)


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

#ifdef PARA
  temps_init_deb = MPI_Wtime()
#endif
  ! Initilisation
  call init
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
       call dmloop_vverlet 
  case(8)
       call dmloop_lpr 
  case (1)
       if (.not.parallele)  call dmloop 
  case (2)
       if (.not.parallele)  then
          call dmloop
       else
#ifdef PARA
          if(rang==0) write (6,*)'DMTYPE+PARA=DMLOOP_VVERLET_+OPTION'
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
       call ndm2config(atdml,im,imm,xp,fp,vp,xpp,ityp,ielat,num_at_glob,ltabvois,iwmax,indi)
    CALL CalFo(atdml) !(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
!    call config2ndm(atdml,im,imm,potist,sig,xp,fp,vp,xpp,ityp,ielat,ltabvois,iwmaxCF,indiCF)
    iwmax=iwmaxCF
    indi=indiCF
!     call calfo()
     call analyse()
     call controle()
     call endrun()

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



  end select


end subroutine prog
end module
