module Tpara
#ifdef PARA
  use mpi
  integer, parameter :: NDM_MPI_REAL_DOUBLE = MPI_REAL8
#endif
  
end module Tpara
