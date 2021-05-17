
module paraconfig

use mpi
  implicit none
#ifdef PARA
!  include 'mpif.h'

  contains
  subroutine commconstr
    call MPI_BARRIER(MPI_COMM_WORLD)
    
    return

  end subroutine commconstr



    
end module paraconfig



