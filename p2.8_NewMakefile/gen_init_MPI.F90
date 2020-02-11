#if defined PARA || defined PARAML || defined PARAPH 
module gen_mpi
     integer:: code_mpi, nb_procs_mpi, rang_mpi
end module gen_mpi

subroutine gen_init_mpi()

  use mpi
  use gen_mpi, only: code_mpi, nb_procs_mpi, rang_mpi
  implicit none

  ! Routine d'initialisation de MPI pour le code NDM


  !--------------------------------------------------
  !Variables de la routine

  !--------------------------------------------------
  !Variables locales

  !--------------------------------------------------
  !Corps de la routine

  call MPI_INIT(code_mpi)
  call MPI_COMM_RANK( MPI_COMM_WORLD, rang_mpi, code_mpi )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nb_procs_mpi, code_mpi )



end subroutine gen_init_mpi
#endif
