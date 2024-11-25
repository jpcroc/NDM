! ****************************************
! Programme de dynamique moleculaire DM
!   Programme parallele/sequentiel
! ****************************************

program ndm
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY: fnam,lenfnam,rang,low_limit

  USE prog_mod,only: prog
  USE readdm_mod,only: readdm
  USE arret_ndm_mod,only: arret_ndm
  USE init_mpi_mod,only: init_mpi
#ifdef PARA
  USE Tpara,only:myidsp,nprocs,mpi_comm_world
  USE neb_module,only:init_mpi_neb
#else
  USE Tpara,only:myidsp,nprocs,nprocspace
#endif

#ifdef ML
  use mld_interface_mod, only: mld_copy_fnam, mld_init_mpi
#endif

#ifdef MAB
  USE mod_mpi_mab
#endif

#if defined PHONDY && defined PARAPH
 USE mod_mpi_phondy
#endif

 
  implicit none
  character :: a1*20
  integer::ierr
  !
  !Initialisation MPI
  call init_MPI()
#ifdef PARA


  write(6,*) 'Process ', rang, ' of ', nprocs, ' is alive'
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  myidsp=rang
#else
  rang = 0
  myidsp=rang
  nprocs=1
  nprocspace=nprocs
  
#endif

#ifdef ML

  call mld_init_mpi()
#endif




  open(29, file='name.in', status='unknown')
  read (29, *) a1
  fnam = a1
  lenfnam = index(fnam,' ')-1
#ifdef ML
  call mld_copy_fnam(fnam)
#endif
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
