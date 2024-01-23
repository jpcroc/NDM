module init_mpi_mod


  use T_kind_param_m, ONLY:  double

!  use atomconfig,only: atom_config

  implicit none
contains

  subroutine init_mpi()
    use gen_com_m,only:rang

    use Tpara,only:comm_space
#ifdef PARA
    use mpi
    use Tpara,only:ierr,nprocs,MPI_COMM_space,grp_world,nprocspace,comm_space,NDM_MPI_REAL_DOUBLE,mpi_world


    implicit none

    ! Routine d'initialisation de MPI pour le code NDM


    call MPI_INIT(ierr)
    call MPI_COMM_RANK( MPI_COMM_WORLD, rang, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )
    call MPI_COMM_GROUP( MPI_COMM_WORLD, grp_world, ierr )
    call mpi_world%init(MPI_COMM_WORLD)

    call MPI_COMM_DUP(MPI_COMM_WORLD,MPI_COMM_SPACE,ierr)
!    MPI_COMM_space=MPI_COMM_WORLD

    call comm_space%init(MPI_COMM_SPACE)
    nprocspace=comm_space%nproc
#else
    call comm_space%init(1)
#endif
  end subroutine init_mpi


   
end module init_mpi_mod



