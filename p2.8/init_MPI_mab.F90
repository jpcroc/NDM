module mod_mpi_mab 
 use mpi
     integer, dimension(MPI_STATUS_SIZE) :: statut
     integer :: nb_procsph,codeph
end module mod_mpi_mab




subroutine init_mpi_mab()

  use mpi
  use mod_mpi_mab
  use gen_com_m , ONLY: rangph
  implicit none

  ! Routine d'initialisation de MPI pour le code NDM


  !--------------------------------------------------
  !Variables de la routine

  !--------------------------------------------------
  !Variables locales

  !--------------------------------------------------
  !Corps de la routine

call MPI_INIT (codeph)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procsph, codeph)
call MPI_COMM_RANK(MPI_COMM_WORLD,rangph,codeph)


end subroutine init_mpi_mab
