! ****************************************
! Programme de dynamique moleculaire DM
!   Programme parallele/sequentiel
! ****************************************

program ndm
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
#if(PARA)
  use mod_mpi
#endif

#if(MAB)
  use mod_mpi_mab
#endif



#if(PHONDY && PARAPH)
 use mod_mpi_phondy
#endif
#if(ML && PARAML)
 use mod_mpi_ml
#endif


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
  character :: a1*20
  !-----------------------------------------------
  !
  !Initialisation MPI
#if(PARA)
  call init_MPI()

 PRINT *, 'Process ', myid, ' of ', nprocs, ' is alive'
  rang = myid
  parallele = .true.
#else
  rang = 0
  parallele = .false.
#endif

#if(ML || PARAML)
rangml=0
#endif
#if(ML && PARAML)
  call init_mpi_ml()
  rang=rangml
#endif


#if(PHONDY || PARAPH)
rangph=0
#endif
#if(PHONDY && PARAPH)
  call init_mpi_phondy()
  rang=rangph
#endif



  if (rang==0) write(6,*)'*** NDM628 ***'
#if(ART)
  if (rang==0) write(6,*)'*** NDM628+ ART ***'
#endif

#if(PHONDY)
  if (rang==0) write(6,*)'*** NDMP628 +  PHONDY ***'
#endif

#if(MAB)
  rangph=0
#if(LAMMPS_VERSION)
  call init_mpi_mab()
#endif
  if (rang==0) write(6,*)'*** NDMP628 +   MAB ***'
#endif

#if(ML)
  if (rang==0) write(6,*)'*** NDMP628 +   ML ***'
#endif



     
  open(29, file='name.in', status='unknown')
  read (29, *) a1
  fnam = a1
  lenfnam = index(fnam,' ')-1
  !     write(6,*) 'main -> readdm'
  call readdm

  call prog

#if (PARAPH || MAB)
continue
#else
  call arret_ndm
#endif

end program ndm



! ca c'est du programme, papa 
