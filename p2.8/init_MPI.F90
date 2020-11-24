module init_mpi_mod


#ifdef PARA
  use Tpara, ONLY:   NDM_MPI_REAL_DOUBLE
#endif
  use T_kind_param_m, ONLY:  double

!  use atomconfig,only: atom_config

  implicit none
contains
#ifdef PARA
  subroutine init_mpi()
    use mpi
    use mod_para,only:ierr,nprocs,temps_deb,MPI_COMM_space,rang,grp_world,nprocspace

    implicit none

    ! Routine d'initialisation de MPI pour le code NDM


    !--------------------------------------------------
    !Variables de la routine

    !--------------------------------------------------
    !Variables locales

    !--------------------------------------------------
    !Corps de la routine

    call MPI_INIT(ierr)
    call MPI_COMM_RANK( MPI_COMM_WORLD, rang, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )
    call MPI_COMM_GROUP( MPI_COMM_WORLD, grp_world, ierr )
    !  MPI_COMM_space=MPI_COMM_WORLD
    temps_deb = MPI_Wtime()
    nprocspace=nprocs

  end subroutine init_mpi

#endif


    
end module init_mpi_mod



