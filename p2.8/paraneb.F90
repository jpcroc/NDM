
module paraneb_mod
#ifdef PARANEB
  USE T_kind_param_m, ONLY:  double
  !USE mpi
  !  USE T_kind_param_mpi_m
  implicit none

  include 'mpif.h'
  integer :: NDM_MPI_REAL_DOUBLE = MPI_REAL8

  ! Module de declaration des variables MPI pour le code NDM
!  Integer :: COMM_NEB = MPI_COMM_WORLD
  !Entiers :
  real(double)::temps_deb
  integer :: myid 			! numero de process
  integer :: np 			! nombre de process
  integer :: ierr  ! erreur MPI
  integer::nprocs
  integer,dimension(MPI_STATUS_SIZE):: status  ! statut de la communication


contains

  subroutine init_mpi_neb()


    ! Routine d'initialisation de MPI pour la NEB


    !--------------------------------------------------
    !Variables de la routine

    !--------------------------------------------------
    !Variables locales

    !--------------------------------------------------
    !Corps de la routine

    call MPI_INIT(ierr)
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )
    write(6,*)'INPNEB', myid,nprocs

    temps_deb = MPI_Wtime()
  end subroutine init_mpi_neb
#endif
end module paraneb_mod

