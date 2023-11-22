module mpi_initialization
  use defs
#ifdef PARA
  use mpi
#endif
  implicit none
contains
  subroutine mpi_group_creation
    use defs
    integer :: error,i,istat
#ifdef PARA
    call mpi_init(error)
    call mpi_comm_rank(MPI_COMM_WORLD, iproc,error)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, error)
    call mpi_comm_group(MPI_COMM_WORLD, orig_group, error)
    call mpi_comm_dup(MPI_COMM_WORLD, orig_comm, error)
#endif

    lammps_ranksize = nproc     ! SLAVE per group should be declare in KMC.sh
    allocate(all_rank(0:nproc-1),STAT=istat)
    do i = 0 , nproc-1
       all_rank(i) = i
    enddo

  end subroutine mpi_group_creation


end module mpi_initialization
