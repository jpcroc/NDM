module atomconfig
  use mpi
  implicit none


  integer,dimension(MPI_STATUS_SIZE):: status  ! statut de la communication
  integer::ierr
  integer, parameter :: NDM_MPI_REAL_DOUBLE = MPI_REAL8
  integer::rang

contains
  !initialisations
  
  subroutine init_atom_config
    real::x
    integer::i
    logical:: lt
    call MPI_INIT(ierr)
  call MPI_COMM_RANK( MPI_COMM_WORLD, rang, ierr )  
  
    if (rang==0) then 
       call MPI_SEND(x,1, NDM_MPI_REAL_DOUBLE,1 ,10002,MPI_COMM_WORLD,status,ierr)
       call MPI_SEND(i,1, MPI_INTEGER, 1, 10005, MPI_COMM_WORLD, status, ierr)
       call MPI_SEND(lt,1, MPI_LOGICAL, 1, 10006, MPI_COMM_WORLD, status, ierr)
    elseif (rang==1) then
       call MPI_RECV(x,1, NDM_MPI_REAL_DOUBLE,1 ,10002,MPI_COMM_WORLD,ierr)
       call MPI_RECV(i,1, MPI_INTEGER, 1, 10005, MPI_COMM_WORLD, ierr)
       call MPI_RECV(lt,1, MPI_LOGICAL, 1, 10006, MPI_COMM_WORLD,  ierr)
    end if
     end subroutine init_atom_config
    
end module atomconfig


