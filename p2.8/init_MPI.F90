#ifdef PARA
module init_mpi_mod
        implicit none
        contains
subroutine init_mpi()
  use mpi
  use mod_para,only:ierr,myid,nprocs,temps_deb,MPI_COMM_space
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
  MPI_COMM_space=MPI_COMM_WORLD
  temps_deb = MPI_Wtime()

end subroutine init_mpi
end module

#endif

