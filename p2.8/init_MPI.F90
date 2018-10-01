subroutine init_mpi()

  use mod_mpi
  implicit none

  ! Routine d'initialisation de MPI pour le code NDM


  !--------------------------------------------------
  !Variables de la routine

  !--------------------------------------------------
  !Variables locales

  !--------------------------------------------------
  !Corps de la routine

  call MPI_INIT(ierr)
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

  temps_deb = MPI_Wtime()

end subroutine init_mpi



