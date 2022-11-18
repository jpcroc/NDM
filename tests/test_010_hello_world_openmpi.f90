! https://www.open-mpi.org/faq/?category=mpi-apps!
! mpifort is a new name for the Fortran wrapper compiler that debuted in Open MPI v1.7.
! Many people base their "wrapper compilers suck!" mentality on bad behavior from poorly-implemented
! wrapper compilers in the mid-1990's.
! Things are much better these days; wrapper compilers can handle almost any situation,
! and are far more reliable than you attempting to hard-code
! the Open MPI-specific compiler and linker flags manually.

! mpifort -v || module load mpi    # useful centos
! afile=test_010_hello_world_openmpi ; mpifort -o ${afile}.exe ${afile}.f90 && mpirun --np 4 ${afile}.exe

! mode ifort
! afile=test_010_hello_world_openmpi
! rm *.exe
! I_MPI_FC=ifort mpifc -o ${afile}.exe ${afile}.f90 && mpirun --np 4 ${afile}.exe
! mpif90 -o ${afile}.exe ${afile}.f90 && mpirun --np 4 ${afile}.exe
! just for test...
! ldd ${afile}.exe

program main
  use mpi
  implicit none

  integer rank, siz, ierror, tag, status(mpi_status_size)

  call mpi_init(ierror)
  write (*, "(a)") "hello world !"
  write (*, "(a, i10)") "mpi_comm_world ", mpi_comm_world
  call mpi_comm_size(mpi_comm_world, siz, ierror)
  call mpi_comm_rank(mpi_comm_world, rank, ierror)
  write (*, "(a, i2, a, i2)") "hello world ! from rank ", rank, ' of ', siz
  call mpi_finalize(ierror)
end program
