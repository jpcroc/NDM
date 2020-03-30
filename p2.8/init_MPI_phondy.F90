#ifdef PHONDY
module mod_para_phondy 
 use mpi
     integer, dimension(MPI_STATUS_SIZE) :: statut
     integer :: nb_procsph,codeph
end module mod_para_phondy




subroutine init_mpi_phondy()

  use mpi
  use mod_para_phondy
  use gen_mpi
  use gen_com_m , ONLY: rangph
  implicit none

  ! Routine d'initialisation de MPI pour le code NDM


  !--------------------------------------------------
  !Variables de la routine

  !--------------------------------------------------
  !Variables locales

  !--------------------------------------------------
  !Corps de la routine

codeph=code_mpi
!call MPI_INIT (codeph)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procsph, codeph)
call MPI_COMM_RANK(MPI_COMM_WORLD,rangph,codeph)


end subroutine init_mpi_phondy
#endif
