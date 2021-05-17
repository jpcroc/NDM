
module paraconfig


#ifdef PARA
  use Tpara,only:NDM_MPI_REAL_DOUBLE,mpi_communicator

#endif
  use T_kind_param_m, ONLY:  double

use mpi
  implicit none
#ifdef PARA
!  include 'mpif.h'
#endif
  contains
    subroutine commconstr
      integer::ierr
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    return

  end subroutine commconstr



    
end module paraconfig



