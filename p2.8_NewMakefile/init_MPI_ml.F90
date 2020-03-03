#ifdef ML

module time_measure
real(kind(0.d0)) :: temps_energy, temps_force, temps_descripteurs, temps_neigh, temps_stress
end module time_measure


module mod_mpi_ml 
 use mpi
     integer, dimension(MPI_STATUS_SIZE) :: statut
     integer :: nb_procsml,codeml
end module mod_mpi_ml



module init_mpi_ml_mod 
  implicit none
  contains  

subroutine init_mpi_ml()

  use mpi
  use mod_mpi_ml
  use gen_mpi
  use gen_com_m , ONLY: rangml
  implicit none

  ! Routine d'initialisation de MPI pour le code NDM+ML


  !--------------------------------------------------
  !Variables de la routine

  !--------------------------------------------------
  !Variables locales

  !--------------------------------------------------
  !Corps de la routine

codeml=code_mpi
!call MPI_INIT (codeml)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procsml, codeml)
call MPI_COMM_RANK(MPI_COMM_WORLD,rangml,codeml)


end subroutine init_mpi_ml
end module


#endif
