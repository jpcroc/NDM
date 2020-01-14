subroutine prog
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m

#if(PARA)
  use mod_mpi
#endif

  implicit none
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

#if(PARA)
  temps_init_deb = MPI_Wtime()
#endif
  ! Initilisation
  call init
#if(PARA)
  temps_init=MPI_Wtime()-temps_init_deb
#endif

#if(DECOUP)
  ! Dans ce cas, pas la peine d'aller plus loin on peut terminer le programme
  return
#endif
#if(PARA)
  ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs
!	if (rang==0) write(6,*)'PARA-T avant MAJ'
  call maj_atomes_frt_ftm


  ! Affichage du temps d'initialisation
#if(PARA)
  if (myid==0) then
     print *, 'Temps d''initialisation : ', MPI_Wtime() - temps_deb
  endif
#endif

#endif

#if(PARA)
  temps_dmloop_deb = MPI_Wtime()
#endif
  ! Actuellement uniquement le cas dmloop_vverlet est traite en parallele


  select case (dmtype) 
  case(5)
   if (.not.parallele)    call loopforcetest (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
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
      if(rang==0) write (6,*)'DMTYPE+PARA=DMLOOP_VVERLET_+OPTION'
      call dmloop_vverlet (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
   endif
  case (3,30) 
    if (.not.parallele)   call gcII(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
  case (9)
     if (.not.parallele)   call neb(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
  case(11)
     if (rang==0) write (6, *) '***** PREMIERE ET UNIQUE ITERATION  ****'
     CALL calfo()
     CALL analyse()
     CALL controle()
     CALL endrun()

#if(ART)    
     case (12) 
          call art90
#endif

#if(SUNDAE)    
     case (16) 
          call sundae
#endif

#if(MAB)    
     case (17) 
          call mab
#endif

#if(PHONDY | PARAPH)    
     case (7) 
          call phondy
#endif

#if(ML || PARAML)    
     case (18) 
          call ml
#endif



  end select


end subroutine prog
