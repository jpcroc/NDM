
module paraconfig

use mpi
 ! implicit none

!  include 'mpif.h'

  contains
  subroutine commconstr
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    return

  end subroutine commconstr



    
end module paraconfig



