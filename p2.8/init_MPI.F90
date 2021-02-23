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
    use gen_com_m,only:rang
    use Tpara,only:ierr,nprocs,MPI_COMM_space,grp_world,nprocspace,myidsp
    use mod_para,only:temps_deb

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
    call MPI_COMM_DUP(MPI_COMM_WORLD,MPI_COMM_SPACE,ierr)
!    MPI_COMM_space=MPI_COMM_WORLD
    temps_deb = MPI_Wtime()
    nprocspace=nprocs

  end subroutine init_mpi

#endif
   
end module init_mpi_mod



