! ****************************************
! Programme de dynamique moleculaire DM
!   Programme parallele/sequentiel
! ****************************************

program ndm
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY: 
  USE prog_mod
  USE readdm_mod
  USE arret_ndm_mod
#ifdef PARA
  USE mod_para
  USE init_mpi_mod
#endif

#ifdef MAB
  USE mod_mpi_mab
#endif

#if defined PHONDY && defined PARAPH
 USE mod_mpi_phondy
#endif
#if defined ML && defined PARAML
 USE mod_mpi_ml
 USE init_mpi_ml_mod
 USE gen_init_mpi_mod

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
#ifdef PARA
  call init_MPI()

 PRINT *, 'Process ', myid, ' of ', nprocs, ' is alive'
  rang = myid
  parallele = .true.
#else
  rang = 0
  parallele = .false.
#endif



#if defined PARAML || defined PARAPH || defined MAB
  call gen_init_mpi
#endif

#if defined ML || defined PARAML
  rangml=0
#endif

#if defined ML && defined PARAML
  call init_mpi_ml()
  rang=rangml
#endif


#if defined PHONDY || defined PARAPH
  rangph=0
#endif

#if defined PHONDY && defined PARAPH
  call init_mpi_phondy()
  rang=rangph
#endif



  if (rang==0) write(6,*)'*** NDM859 ***'
#ifdef ART
  if (rang==0) write(6,*)'*** NDM859+ ART ***'
#endif

#ifdef PHONDY
  if (rang==0) write(6,*)'*** NDMP859 +  PHONDY ***'
#endif


!if MAB .....
#ifdef MAB
  rangmab=0
  if (rang==0) write(6,*)'*** NDMP859 +   MAB ***'
#if defined ML && defined PARAML
  if (rang==0) write(6,*)'*** NDMP859 +   MAB + ML + PARAML ***'
  call init_mpi_mab()
  rang=rangmab
#endif
#ifdef LAMMPS_VERSION
  call init_mpi_mab()
#endif
#endif
!endif MAB ......

#ifdef ML
  if (rang==0) write(6,*)'*** NDMP859 +   ML ***'
#endif




  open(29, file='name.in', status='unknown')
  read (29, *) a1
  fnam = a1
  lenfnam = index(fnam,' ')-1
  !     write(6,*) 'main -> readdm'
  call readdm

  call prog

#if defined PARAPH || defined MAB
continue
#else
  call arret_ndm
#endif

end program ndm



! ca c'est du programme, papa
