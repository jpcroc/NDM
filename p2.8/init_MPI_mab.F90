#ifdef MAB
module mod_para_mab 
 USE mpi
     integer, dimension(MPI_STATUS_SIZE) :: statut
     integer :: nb_procsmab,codemab
end module mod_para_mab




subroutine init_mpi_mab()

  USE mpi
  USE mod_para,only:_mab
  USE gen_mpi
  USE gen_com_m, ONLY: , ONLY: rangmab
  implicit none

  ! Routine d'initialisation de MPI pour le code NDM


  !--------------------------------------------------
  !Variables de la routine

  !--------------------------------------------------
  !Variables locales

  !--------------------------------------------------
  !Corps de la routine
codemab=code_mpi
!call MPI_INIT (codemab)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procsmab, codemab)
call MPI_COMM_RANK(MPI_COMM_WORLD,rangmab,codemab)


end subroutine init_mpi_mab
#endif
