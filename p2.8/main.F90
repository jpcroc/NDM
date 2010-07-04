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

!  PRINT *, 'Process ', myid, ' of ', nprocs, ' is alive'
  rang = myid
  parallele = .true.
#else
  rang = 0
  parallele = .false.
#endif

#if(ART)
  if (rang==0) write(6,*)'*** NDMP-2.4 + ART ***'
#else
  if (rang==0) write(6,*)'*** NDMP-2.4 ***'
#endif

  !     read(5,*)a1
  ! modif pour compaq


     
  open(29, file='name.in', status='unknown')
  read (29, *) a1
  fnam = a1
  lenfnam = index(fnam,' ')-1
  !     write(6,*) 'main -> readdm'
  call readdm

  call prog

  call arret_ndm

end program ndm



! ca c'est du programme, papa !
