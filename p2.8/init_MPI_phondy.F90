#ifdef PHONDY
module mod_para_phondy 
 USE mpi
     integer, dimension(MPI_STATUS_SIZE) :: statut
     integer :: nb_procsph,codeph
end module mod_para_phondy




subroutine init_mpi_phondy()

  USE mpi
  USE mod_para,only:MPI_COMM_space,_phondy
  USE gen_mpi
  USE gen_com_m, ONLY: , ONLY: rangph
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
call MPI_COMM_SIZE(MPI_COMM_space,nb_procsph, codeph)
call MPI_COMM_RANK(MPI_COMM_space,rangph,codeph)


end subroutine init_mpi_phondy
#endif
