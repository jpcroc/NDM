program mpitest
  
  use mpi
  implicit none
  integer, parameter :: double = 8
 integer, parameter :: NDM_MPI_REAL_DOUBLE = MPI_REAL8
 real*8 :: d
 real (double)::d2
  integer::i,ierr,rang
  integer, dimension(MPI_STATUS_SIZE) :: status
  call MPI_INIT(ierr)
  call MPI_COMM_RANK( MPI_COMM_WORLD, rang, ierr )  
  i=rang
  d=dfloat(rang)
  d2=rang
  if (rang==0) then
     call MPI_SEND(d, 1, MPI_REAL8, 2,104,MPI_COMM_WORLD,ierr)
  elseif (rang==2) then
     call MPI_RECV(d, 1, MPI_REAL8, 0,104,MPI_COMM_WORLD,status, ierr)
  end if
  if (rang==4) then
     call MPI_SEND(i, 1, MPI_INTEGER, 5,104,MPI_COMM_WORLD,ierr)
  elseif (rang==5) then
     call MPI_RECV(i, 1, MPI_INTEGER, 4,104,MPI_COMM_WORLD,status, ierr)
  end if
  if (rang==1) then
     call MPI_SEND(d2, 1, NDM_MPI_REAL_DOUBLE , 3,105,MPI_COMM_WORLD,ierr)
  elseif (rang==3) then
     call MPI_RECV(d2, 1,  NDM_MPI_REAL_DOUBLE, 1,105,MPI_COMM_WORLD,status, ierr)
  end if
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  write(6,*)rang,d,d2,i
  call mpi_finalize (ierr)
  stop
end program mpitest
