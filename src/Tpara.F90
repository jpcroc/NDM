module Tpara
  
  use T_kind_param_m
#ifdef PARA
  use mpi

  integer, parameter :: NDM_MPI_REAL_DOUBLE = MPI_REAL8
  integer, parameter :: NDM_MPI_COMPLEX_DOUBLE = MPI_COMPLEX16
  integer::MPI_COMM_space
  integer:: grp_world
  integer,dimension(MPI_STATUS_SIZE):: status

#else
  integer:: status 
#endif
  integer,target :: nprocs 			! numero de process mis là pour être utilisé en sequentiesl
  integer :: myidsp,nprocspace 			! numero de process mis là pour être utilisé en sequentiesl

  integer::ierr
  type para_space_config
     integer, allocatable :: res_cpu(:,:)   	!stocke le nombre de cellules de chaques decoupages pour le meilleur decoupage
     integer, allocatable :: proc_voisin(:)        ! liste des processeurs voisins du processeur courant
     integer, allocatable :: cell_frontiere(:,:)   ! (i,j) jeme cellule frontiere associee au ieme processeur voisin
     integer, allocatable :: nbr_cell_frontiere(:) ! nbre de cellules frontieres associees au ieme processeur voisin
     integer :: nbr_cell_ftm                       ! nbr de cellules fantomes du processeur courant
     integer, allocatable :: cell_ftm(:)           ! liste des cellules fantomes du processeur courant
     integer :: nbr_proc_voisin            ! nbre de processeurs voisins du processeur courant
     integer :: cell_debx, cell_deby, cell_debz     !numero de la premiere cellule locale suivant x, y et z
     integer :: cell_finx, cell_finy, cell_finz     !numero de la derniere cellule locale  suivant x, y et z
     integer :: nb_cell_x, nb_cell_y, nb_cell_z     !nb de cel locales suivant x y z
   contains
     procedure::print=>printpsc
  end type para_space_config


  type mpi_communicator
     integer    :: comm       ! MPI communicator
     integer    :: nproc      ! number of procs in the communicator comm
     integer    :: rank       ! index           in the communicator comm
     integer    :: group       ! group          of the communicator comm
   contains
     procedure :: init => mpic_init
     procedure :: init0 => mpic_init0
     procedure :: probe => mpic_probe
     procedure :: barrier => mpic_barrier
     procedure :: print
     ! sum
     generic :: send  => mpic_send_dp,mpic_send_cdp,mpic_send_i,mpic_send_l
     generic :: recv  => mpic_recv_dp,mpic_recv_cdp,mpic_recv_i,mpic_recv_l
     generic :: sum  => mpic_sum_dp
     generic :: sum  => mpic_sum_cdp
     generic :: sum  => mpic_sum_i
     procedure :: mpic_sum_dp
     procedure :: mpic_sum_cdp
     procedure :: mpic_sum_i
     procedure:: mpic_send_dp,mpic_send_cdp,mpic_send_i,mpic_send_l
     procedure:: mpic_recv_dp,mpic_recv_cdp,mpic_recv_i,mpic_recv_l
     ! min ALLREDUCE !!
     generic :: min  => mpic_min_dp
     generic :: min  => mpic_min_i
     procedure :: mpic_min_dp
     procedure :: mpic_min_i
     ! max ALLREDUCE !!
     generic :: max  => mpic_max_dp
     generic :: max  => mpic_max_i
     generic::maxloc=>mpic_maxloc_dp
     procedure::mpic_maxloc_dp
     procedure :: mpic_max_dp
     procedure :: mpic_max_i
     ! and
     procedure :: and => mpic_and_l
     ! broadcast
     generic :: bcast  => mpic_bcast_dp
     generic :: bcast  => mpic_bcast_l
     generic :: bcast  => mpic_bcast_i
     generic :: bcast  => mpic_bcast_cdp
!!$    generic :: average  => mpic_av_dp
!!$    generic :: average  => mpic_av_i
!!$    generic :: average  => mpic_av_cdp
!!$    procedure :: mpic_av_dp,mpic_av_cdp,mpic_av_i
     procedure :: mpic_bcast_dp,mpic_bcast_l
     procedure:: mpic_bcast_i
     procedure :: mpic_bcast_cdp
     generic :: build=>build_i,build_dp,build_cdp,build_l
     procedure :: build_i,build_dp,build_cdp,build_l

  end type mpi_communicator

  type(mpi_communicator),target::comm_space
  type(mpi_communicator)::mpi_world

contains

  subroutine endmpi
#ifdef PARA    
    call MPI_finalize(ierr)
    stop
#endif
    return
  end subroutine endmpi

  subroutine print (mpic,unit)
    class(mpi_communicator),intent(in) :: mpic
    integer,intent(in)::unit
    write(unit,*)'MPICOMM rank comm group nproc',mpic%rank,mpic%comm,mpic%group,mpic%nproc
    flush(unit)
  end subroutine print
  !=========================================================================
  subroutine mpic_init(mpic,comm_in)
    implicit none

    class(mpi_communicator),intent(inout) :: mpic
    integer,intent(in)                    :: comm_in
    !=====
    integer :: ierror
    !=====

#if defined(PARA)
    mpic%comm = comm_in
    call MPI_COMM_SIZE(mpic%comm,mpic%nproc,ierror)
    call MPI_COMM_RANK(mpic%comm,mpic%rank,ierror)
    call MPI_COMM_GROUP(mpic%comm, mpic%group, ierr )
#else
    mpic%comm  = 1
    mpic%nproc = 1
    mpic%rank  = 0
#endif

  end subroutine mpic_init

  !=========================================================================
  subroutine mpic_init0(mpic)
    implicit none

    class(mpi_communicator),intent(inout) :: mpic
    !=====
    integer :: ierror
    !=====

#if defined(PARA)
#else
    mpic%comm  = 1
    mpic%nproc = 1
    mpic%rank  = 0
#endif

  end subroutine mpic_init0


  !=========================================================================
  subroutine mpic_probe(mpic,tag,sourceout,sourcein)
    implicit none

    class(mpi_communicator),intent(inout) :: mpic
    integer,intent(in)                    :: tag
    integer,optional,intent(in)::sourcein
    integer,optional,intent(out)::sourceout
    integer :: ierror
#if defined(PARA)
    integer,dimension(MPI_STATUS_SIZE):: statut
    ierror=0
    if (present(sourcein)) then
       call MPI_PROBE(sourcein, tag, mpic%comm,statut,ierror)
       sourceout=statut(MPI_SOURCE)
       if (sourceout.ne.sourcein) then
          write(6,*) 'error in mpic_probe sourceout<> sourcein',sourceout,sourcein
       endif
    else
       call MPI_PROBE(MPI_ANY_SOURCE, tag, mpic%comm,statut,ierror)
       sourceout=statut(MPI_SOURCE)
    end if

#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in mpic_probe'
    endif

  end subroutine mpic_probe


  !=========================================================================
  subroutine mpic_barrier(mpic)
    implicit none

    class(mpi_communicator),intent(in) :: mpic
    !=====
    integer :: ierror
    !=====

#if defined(PARA)
    call MPI_BARRIER(mpic%comm,ierror)
#endif

  end subroutine mpic_barrier


  !=========================================================================
  subroutine mpic_sum_dp(mpic,array,torank)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    real(double),intent(inout) :: array(..)
    integer,intent(in),optional::torank
    !=====
    integer :: nsize,trk
    integer :: ierror=0
    !=====
    if(present(torank))trk=torank
    if( mpic%nproc == 1 ) return

    nsize = SIZE(array)

#if defined(PARA)
    if (present(torank))then
       if (mpic%rank==torank) then
          call MPI_REDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_PRECISION, MPI_SUM, trk,mpic%comm, ierror)
       else
          call MPI_REDUCE(array, array, nsize, MPI_DOUBLE_PRECISION, MPI_SUM, trk,mpic%comm, ierror)
       end if
    else
       call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_PRECISION, MPI_SUM, mpic%comm, ierror)
    end if
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_ALLREDUCE'
    endif

  end subroutine mpic_sum_dp


  !=========================================================================
  subroutine mpic_sum_cdp(mpic,array,torank)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    complex(double),intent(inout) :: array(..)

    integer,intent(in),optional::torank
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====

    if( mpic%nproc == 1 ) return

    nsize = SIZE(array)

#if defined(PARA)
    if (present(torank))then
       if (mpic%rank==torank) then
          call MPI_REDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_COMPLEX, MPI_SUM, torank,mpic%comm, ierror)
       else
          call MPI_REDUCE( array, array, nsize, MPI_DOUBLE_COMPLEX, MPI_SUM, torank,mpic%comm, ierror)
       end if
    else
       call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_COMPLEX, MPI_SUM, mpic%comm, ierror)
    end if
#endif

    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_ALLREDUCE'
    endif

  end subroutine mpic_sum_cdp


  !=========================================================================
  subroutine mpic_sum_i(mpic,array,torank)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    integer,intent(inout) :: array(..)

    integer,intent(in),optional::torank
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====

    if( mpic%nproc == 1 ) return

    nsize = SIZE(array)

#if defined(PARA)
    if (present(torank))then
       if (mpic%rank==torank) then
          call MPI_REDUCE( MPI_IN_PLACE, array, nsize, MPI_INTEGER, MPI_SUM, torank,mpic%comm, ierror)
       else
          call MPI_REDUCE( array, array, nsize, MPI_INTEGER, MPI_SUM, torank,mpic%comm, ierror)
       end if
    else
       call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_INTEGER, MPI_SUM, mpic%comm, ierror)
    end if
#endif

    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_ALLREDUCE'
    endif

  end subroutine mpic_sum_i

  !=========================================================================
  subroutine mpic_maxloc_dp(mpic,array,torank)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    real(double),intent(inout) :: array(2)
    real(double)::array_glob(2)
    integer,intent(in),optional::torank

    !=====
    integer :: nsize
    integer :: ierror=0
    !=====

    if( mpic%nproc == 1 ) return

    !  nsize = SIZE(array)

#if defined(PARA)

    if (present(torank))then
       call MPI_REDUCE(array,array_glob,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,torank,MPIc%COMM,ierr) 
    else
       call MPI_ALLREDUCE(array,array_glob,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPIc%COMM,ierr)
    end if


    !  call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_PRECISION, MPI_MAX, mpic%comm, ierror)
    array=array_glob
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_ALLREDUCE'
    endif

  end subroutine mpic_maxloc_dp


  !=========================================================================

  !=========================================================================
  subroutine mpic_max_dp(mpic,array,torank)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    real(double),intent(inout) :: array(..)
    integer,intent(in),optional::torank
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====
    if( mpic%nproc == 1 ) return
    nsize = SIZE(array)
#if defined(PARA)
    if (present(torank))then
       if (mpic%rank==torank) then
          call MPI_REDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_PRECISION, MPI_MAX, torank,mpic%comm, ierror)
       else
          call MPI_REDUCE( array, array, nsize, MPI_DOUBLE_PRECISION, MPI_MAX, torank,mpic%comm, ierror)
       end if
    else
       call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_PRECISION, MPI_MAX, mpic%comm, ierror)
    end if
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_ALLREDUCE'
    endif

  end subroutine mpic_max_dp


  !=========================================================================
  subroutine mpic_max_i(mpic,array,torank)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    integer,intent(inout) :: array(..)
    integer,intent(in),optional::torank
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====

    if( mpic%nproc == 1 ) return
    nsize = SIZE(array)
#if defined(PARA)
    if (present(torank))then
       if (mpic%rank==torank) then
          call MPI_REDUCE( MPI_IN_PLACE, array, nsize, MPI_INTEGER, MPI_MAX, torank,mpic%comm, ierror)
       else
          call MPI_REDUCE( array, array, nsize, MPI_INTEGER, MPI_MAX, torank,mpic%comm, ierror)
       end if
    else
       call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_INTEGER, MPI_MAX, mpic%comm, ierror)
    end if
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_ALLREDUCE'
    endif
  end subroutine mpic_max_i


  !=========================================================================
 subroutine mpic_min_dp(mpic,array,torank)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    real(double),intent(inout) :: array(..)
    integer,intent(in),optional::torank
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====
    if( mpic%nproc == 1 ) return
    nsize = SIZE(array)
#if defined(PARA)
    if (present(torank))then
       if (mpic%rank==torank) then
          call MPI_REDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_PRECISION, MPI_MIN,TORANK, mpic%comm, ierror)
       else
          call MPI_REDUCE( array, array, nsize, MPI_DOUBLE_PRECISION, MPI_MIN,TORANK, mpic%comm, ierror)
       end if
      
    else
       call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_PRECISION, MPI_MIN, mpic%comm, ierror)
    end if
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_ALLREDUCE'
    endif

  end subroutine mpic_min_dp


  !=========================================================================
  subroutine mpic_min_i(mpic,array,torank)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    integer,intent(inout) :: array(..)
    integer,intent(in),optional::torank
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====
    if( mpic%nproc == 1 ) return
    nsize = SIZE(array)
#if defined(PARA)
    if (present(torank))then
       if (mpic%rank==torank) then
          call MPI_REDUCE( MPI_IN_PLACE, array, nsize, MPI_INTEGER, MPI_MIN, torank,mpic%comm, ierror)
       else
          call MPI_REDUCE(array, array, nsize, MPI_INTEGER, MPI_MIN, torank,mpic%comm, ierror)
       end if
    else
       call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_INTEGER, MPI_MIN, mpic%comm, ierror)
    end if
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_ALLREDUCE'
    endif

  end subroutine mpic_min_i


  !=========================================================================
  subroutine mpic_and_l(mpic,array,torank)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    logical,intent(inout) :: array(..)
    integer,intent(in),optional::torank
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====
    if( mpic%nproc == 1 ) return
    nsize = SIZE(array)
#if defined(PARA)
    if (present(torank))then
       if (mpic%rank==torank) then
          call MPI_REDUCE( MPI_IN_PLACE, array, nsize, MPI_LOGICAL, MPI_LAND, torank,mpic%comm, ierror)
       else
          call MPI_REDUCE(array, array, nsize, MPI_INTEGER, MPI_MIN, torank,mpic%comm, ierror)
       end if
       

    else
       call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_LOGICAL, MPI_LAND, mpic%comm, ierror)
    end if
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_ALLREDUCE'
    endif

  end subroutine mpic_and_l


  !=========================================================================
  subroutine mpic_bcast_dp(mpic,rank,array)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    integer,intent(in)     :: rank
    real(double),intent(inout) :: array(..)
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====
    !  write(6,*)'INBCAST', mpic%nproc
    if( mpic%nproc == 1 ) return

    nsize = SIZE(array)
#if defined(PARA)
    call MPI_BCAST(array,nsize,MPI_DOUBLE_PRECISION,rank,mpic%comm,ierror)
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_BCAST'
    endif

  end subroutine mpic_bcast_dp


  subroutine mpic_bcast_i(mpic,rank,array)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    integer,intent(in)     :: rank
    integer,intent(inout) :: array(..)
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====

    if( mpic%nproc == 1 ) return

    nsize = SIZE(array)

#if defined(PARA)
    call MPI_BCAST(array,nsize,MPI_INTEGER,rank,mpic%comm,ierror)
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_BCAST I'
    endif

  end subroutine mpic_bcast_i

  subroutine mpic_bcast_l(mpic,rank,array)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    integer,intent(in)     :: rank
    logical,intent(inout) :: array(..)
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====

    if( mpic%nproc == 1 ) return

    nsize = SIZE(array)

#if defined(PARA)
    call MPI_BCAST(array,nsize,MPI_LOGICAL,rank,mpic%comm,ierror)
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_BCAST L'
    endif

  end subroutine mpic_bcast_l


  !=========================================================================
  subroutine mpic_bcast_cdp(mpic,rank,array)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    integer,intent(in)     :: rank
    complex(ext_complex),intent(inout) :: array(..)
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====

    if( mpic%nproc == 1 ) return

    nsize = SIZE(array)

#if defined(PARA)
    call MPI_BCAST(array,nsize,MPI_DOUBLE_COMPLEX,rank,mpic%comm,ierror)
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_BCAST'
    endif

  end subroutine mpic_bcast_cdp
!!$

!!$!=========================================================================
!!$subroutine mpic_av_i(mpic,array,nval,torank)
!!$  implicit none
!!$  class(mpi_communicator),intent(in) :: mpic
!!$  integer,intent(in)     :: nval
!!$  integer,intent(in),optional     :: torank
!!$  integer,intent(inout) :: array(..)
!!$  !=====
!!$  integer :: nsize
!!$  integer :: ierror=0
!!$  integer::nvaltot
!!$  !=====
!!$
!!$  if( mpic%nproc == 1 ) return
!!$  nsize = SIZE(array)
!!$  array=array*nval
!!$#if defined(PARA)
!!$  if (present(torank)) then
!!$     call mpic_sum_i(mpic,nval,torank)
!!$     call mpic_sum_i(mpic,array,torank)   
!!$     if (mpic%rank==torank) then
!!$        array=array/nval
!!$     end if
!!$  else
!!$     call mpic_sum_i(mpic,nval)
!!$     call mpic_sum_i(mpic,array)   
!!$     array=array/nval
!!$  end if
!!$  
!!$#endif
!!$end subroutine mpic_av_i
!!$
!!$subroutine mpic_av_dp(mpic,array,nval,torank)
!!$  implicit none
!!$  class(mpi_communicator),intent(in) :: mpic
!!$  integer,intent(in)     :: nval
!!$  integer,intent(in),optional     :: torank
!!$  real(double),intent(inout) :: array(..)
!!$  !=====
!!$  integer :: nsize
!!$  integer :: ierror=0
!!$  integer::nvaltot
!!$  !=====
!!$
!!$  if( mpic%nproc == 1 ) return
!!$  nsize = SIZE(array)
!!$  array=array*nval
!!$#if defined(PARA)
!!$  if (present(torank)) then
!!$     call mpic_sum_i(mpic,nval,torank)
!!$     call mpic_sum_dp(mpic,array,torank)   
!!$     if (mpic%rank==torank) then
!!$        array=array/nval
!!$     end if
!!$  else
!!$     call mpic_sum_i(mpic,nval)
!!$     call mpic_sum_dp(mpic,array)   
!!$     array=array/nval
!!$  end if
!!$  
!!$#endif
!!$end subroutine mpic_av_dp
!!$
!!$subroutine mpic_av_cdp(mpic,array,nval,torank)
!!$  implicit none
!!$  class(mpi_communicator),intent(in) :: mpic
!!$  integer,intent(in)     :: nval
!!$  integer,intent(in),optional     :: torank
!!$  complex(double),intent(inout) :: array(..)
!!$  !=====
!!$  integer :: nsize
!!$  integer :: ierror=0
!!$  integer::nvaltot
!!$  !=====
!!$
!!$  if( mpic%nproc == 1 ) return
!!$  nsize = SIZE(array)
!!$  array=array*nval
!!$#if defined(PARA)
!!$  if (present(torank)) then
!!$     call mpic_sum_i(mpic,nval,torank)
!!$     call mpic_sum_cdp(mpic,array,torank)   
!!$     if (mpic%rank==torank) then
!!$        array=array/nval
!!$     end if
!!$  else
!!$     call mpic_sum_i(mpic,nval)
!!$     call mpic_sum_cdp(mpic,array)   
!!$     array=array/nval
!!$  end if
!!$  
!!$#endif
!!$end subroutine mpic_av_cdp

!!$

  !=========================================================================
  subroutine mpic_send_dp(mpic,array,rgcib,tag)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    real(double),intent(in) :: array(..)
    integer,intent(in)::rgcib
    integer,optional,intent(in)::tag
    !=====
    integer :: nsize,tagv
    integer :: ierror=0
    !=====

    if( mpic%nproc == 1 ) return

    nsize = SIZE(array)

#if defined(PARA)
    if (present(tag)) then
       call MPI_SEND( array, nsize,  NDM_MPI_REAL_DOUBLE, rgcib, tag,mpic%comm, ierror)
    else
       call MPI_SEND( array, nsize,  NDM_MPI_REAL_DOUBLE, rgcib, MPI_ANY_TAG,mpic%comm, ierror)
    end if
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_SEND_DP'
    endif

  end subroutine mpic_send_dp

  subroutine mpic_send_l(mpic,array,rgcib,tag)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    logical,intent(in) :: array(..)
    integer,intent(in)::rgcib
    integer,optional,intent(in)::tag
    !=====
    integer :: nsize,tagv
    integer :: ierror=0
    !=====

    if( mpic%nproc == 1 ) return

    nsize = SIZE(array)

#if defined(PARA)
    if (present(tag)) then
       call MPI_SEND( array, nsize, MPI_LOGICAL, rgcib, tag,mpic%comm, ierror)
    else
       call MPI_SEND( array, nsize,  MPI_LOGICAL, rgcib, MPI_ANY_TAG,mpic%comm, ierror)
    end if
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_SEND_DP'
    endif

  end subroutine mpic_send_l


  !=========================================================================
  subroutine mpic_send_cdp(mpic,array,rgcib,tag)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    complex(ext_complex),intent(in) :: array(..)
    integer,intent(in)::rgcib
    integer,optional,intent(in)::tag
    !=====
    integer :: nsize,tagv
    integer :: ierror=0
    !=====

    if( mpic%nproc == 1 ) return

    nsize = SIZE(array)

#if defined(PARA)
    if (present(tag)) then
       call MPI_SEND( array, nsize, MPI_DOUBLE_COMPLEX, rgcib, tag,mpic%comm, ierror)
    else
       call MPI_SEND( array, nsize,MPI_DOUBLE_COMPLEX, rgcib, MPI_ANY_TAG,mpic%comm, ierror)
    end if
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_SEND_CDP'
    endif

  end subroutine mpic_send_cdp


  !=========================================================================
  subroutine mpic_send_i(mpic,array,rgcib,tag)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    integer,intent(inout) :: array(..)
    integer,intent(in)::rgcib
    integer,optional,intent(in)::tag
    !=====
    integer :: nsize,tagv
    integer :: ierror=0
    !=====

    if( mpic%nproc == 1 ) return
    nsize = SIZE(array)
#if defined(PARA)
    if (present(tag)) then
       call MPI_SEND( array, nsize,  MPI_INTEGER, rgcib, tag,mpic%comm, ierror)
    else
       call MPI_SEND( array, nsize,  MPI_INTEGER, rgcib, MPI_ANY_TAG,mpic%comm, ierror)
    end if
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_SEND_I'
    endif

  end subroutine mpic_send_i


  !=========================================================================
  subroutine mpic_recv_dp(mpic,array,rgem,tag)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    real(double),intent(inout) :: array(..)
    integer,intent(in)::rgem
    integer,optional,intent(in)::tag
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====

    if( mpic%nproc == 1 ) return

    nsize = SIZE(array)
#if defined(PARA)
    if (present(tag)) then
       call MPI_RECV( array, nsize, NDM_MPI_REAL_DOUBLE, rgem, tag,mpic%comm, status,ierror)
    else
       call MPI_RECV( array, nsize, NDM_MPI_REAL_DOUBLE, rgem, MPI_ANY_TAG,mpic%comm,status, ierror)
    end if
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_RECV_DP'
    endif

  end subroutine mpic_recv_dp

  subroutine mpic_recv_l(mpic,array,rgem,tag)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    logical,intent(inout) :: array(..)
    integer,intent(in)::rgem
    integer,optional,intent(in)::tag
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====

    if( mpic%nproc == 1 ) return

    nsize = SIZE(array)
#if defined(PARA)
    if (present(tag)) then
       call MPI_RECV( array, nsize, MPI_LOGICAL, rgem, tag,mpic%comm, status,ierror)
    else
       call MPI_RECV( array, nsize, MPI_LOGICAL, rgem, MPI_ANY_TAG,mpic%comm,status, ierror)
    end if
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_RECV_DP'
    endif

  end subroutine mpic_recv_l


  !=========================================================================
  subroutine mpic_recv_cdp(mpic,array,rgem,tag)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    complex(ext_complex),intent(inout) :: array(..)
    integer,intent(in)::rgem
    integer,optional,intent(in)::tag
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====

    if( mpic%nproc == 1 ) return

    nsize = SIZE(array)
#if defined(PARA)
    if (present(tag)) then
       call MPI_RECV( array, nsize, MPI_DOUBLE_COMPLEX,rgem, tag,mpic%comm, status,ierror)
    else
       call MPI_RECV( array, nsize, MPI_DOUBLE_COMPLEX, rgem, MPI_ANY_TAG,mpic%comm,status, ierror)
    end if
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_RECV_CDP'
    endif

  end subroutine mpic_recv_cdp

  !=========================================================================
  subroutine mpic_recv_i(mpic,array,rgem,tag)
    implicit none
    class(mpi_communicator),intent(in) :: mpic
    integer,intent(inout) :: array(..)
    integer,intent(in)::rgem
    integer,optional,intent(in)::tag
    !=====
    integer :: nsize
    integer :: ierror=0
    !=====
    if( mpic%nproc == 1 ) return
    nsize = SIZE(array)
#if defined(PARA)
    if (present(tag)) then
       call MPI_RECV( array, nsize,MPI_INTEGER, rgem, tag,mpic%comm, status,ierror)
    else
       call MPI_RECV( array, nsize,MPI_INTEGER, rgem, MPI_ANY_TAG,mpic%comm,status, ierror)
    end if
#endif
    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_RECV_I'
    endif

  end subroutine mpic_recv_i


  subroutine printpsc(psc,unit)
    class(para_space_config),intent(in)::psc
    integer, optional::unit
    integer::unitw
    if (present(unit) )then
       unitw=unit
    else
       unitw=6
    end if
    if (allocated(psc%res_cpu))then
       write(unitw,*)'res_cpu',psc%res_cpu
    else
       write(unitw,*)'res_cpu NOT ALLOCATED'
    end if
!!$    if (allocated(psc%proc_cell))then
!!$       write(unitw,*)'proc_cell',psc%proc_cell
!!$    else
!!$       write(unitw,*)'proc_cell NOT ALLOCATED'
!!$    end if
    if (allocated(psc%proc_voisin))then
       write(unitw,*)'proc_voisin',psc%proc_voisin
    else
       write(unitw,*)'proc_voisin NOT ALLOCATED'
    end if
    if (allocated(psc%cell_frontiere))then
       write(unitw,*)'cell_frontiere',psc%cell_frontiere
    else
       write(unitw,*)'cell_frontiere NOT ALLOCATED'
    end if
    if (allocated(psc%nbr_cell_frontiere))then
       write(unitw,*)'nbr_cell_frontiere',psc%nbr_cell_frontiere
    else
       write(unitw,*)'nbr_cell_frontiere NOT ALLOCATED'
    end if
    write(unitw,*)'nbr_cell_ftm',psc%nbr_cell_ftm
    if(allocated(psc%cell_ftm))then
       write(unitw,*)'cell_ftm',psc%cell_ftm
    else
       write(unitw,*)'cell_ftm NOT ALLOCATED'
    end if
    write(unitw,*)'nbr_proc_voisin',psc%nbr_proc_voisin
    write(unitw,*)'cell_deb', psc%cell_debx, psc%cell_deby, psc%cell_debz     !numero de la premiere cellule locale suivant x, y et z
    write(unitw,*)'cell_fin', psc%cell_finx, psc%cell_finy, psc%cell_finz     !numero de la premiere cellule locale suivant x, y et z
    write(unitw,*)'nbr_cell', psc%nb_cell_x, psc%nb_cell_y, psc%nb_cell_z     !nb de cel locales suivant x y z
  end subroutine printpsc

  subroutine build_i(mpic,val,arrayval,torank)
    implicit none
    class(mpi_communicator) :: mpic
    integer,optional,intent(in)::torank

    integer,allocatable :: arrayval(:)
    integer,intent(in)::val
    integer::valp,proc_source
    integer::iproc,ierror,sourceout,statut

    arrayval(:)=0
#ifdef PARA    
    if (present (torank)) then
       if (mpic%rank==torank) then
          do iproc=0,mpic%rank-1
             if (iproc==torank)then 
                arrayval(torank)=val
             else
!                call MPI_PROBE(MPI_ANY_SOURCE, 114, mpic%comm,statut,ierror)
                call mpic%probe(114,sourceout=proc_source)
!                sourceout=statut(MPI_SOURCE)
                !                call MPI_SEND( val, 1,  MPI_INTEGER, torank, 114,mpic%comm, ierror)
                call MPI_RECV( valp, 1, MPI_INTEGER, sourceout, 114,mpic%comm, status,ierror)
!                call mpic%recv(valp,proc_source,999)
                arrayval(proc_source)=valp
             end if
          end do
       else
          call MPI_SEND( val, 1,  MPI_INTEGER, torank, 114,mpic%comm, ierror)
          !call mpic%send(val,torank,999)
       end if
    else
       arrayval(mpic%rank)=val
       call mpic%sum(arrayval)
    end if
#else
    arrayval(0)=val
#endif
    return

  end subroutine build_i
  subroutine build_dp(mpic,val,arrayval,torank)
    implicit none
    class(mpi_communicator) :: mpic
    integer,optional,intent(in)::torank

    real(double),allocatable :: arrayval(:)
    real(double),intent(in)::val
    real(double)::valp

    integer::iproc,proc_source,ierror,sourceout,statut
    arrayval(:)=0
#ifdef PARA

    if (present (torank)) then
       if (mpic%rank==torank) then
          do iproc=0,mpic%rank-1
             if (iproc==torank)then 
                arrayval(torank)=val
             else
                call mpic%probe(999,sourceout=proc_source)
                call MPI_RECV( valp, 1, NDM_MPI_REAL_DOUBLE, sourceout, 115,mpic%comm, status,ierror)
!                call mpic%recv(valp,proc_source,999)
                arrayval(proc_source)=valp
             end if
          end do
       else
          call MPI_SEND( val, 1,  NDM_MPI_REAL_DOUBLE, torank, 115,mpic%comm, ierror)
       end if
    else
       arrayval(mpic%rank)=val
       call mpic%sum(arrayval)
    end if
#else
    arrayval(0)=val
#endif
    return


  end subroutine build_dp
  subroutine build_cdp(mpic,val,arrayval,torank)
    implicit none
    class(mpi_communicator) :: mpic
    integer,optional,intent(in)::torank

    complex(double),allocatable :: arrayval(:)
    complex(double),intent(in)::val
    complex(double)::valp
    integer::proc_source,ierror,sourceout,statut
    integer::iproc
    arrayval(:)=0
#ifdef PARA

    if (present (torank)) then
       if (mpic%rank==torank) then
          do iproc=0,mpic%rank-1
             if (iproc==torank)then 
                arrayval(torank)=val
             else
                call mpic%probe(999,sourceout=proc_source)
                call MPI_RECV( valp, 1, MPI_DOUBLE_COMPLEX, sourceout, 116,mpic%comm, status,ierror)
!                call MPI_RECV( valp, 1,MPI_DOUBLE_COMPLEX, proc_source, MPI_ANY_TAG,mpic%comm,status, ierror)
                !call mpic%recv(valp,proc_source,999)
                arrayval(proc_source)=valp
             end if
          end do
       else
          !call mpic%send(val,torank,999)
          call MPI_SEND( val, 1,  MPI_DOUBLE_COMPLEX, torank, 116,mpic%comm, ierror)
!          call MPI_SEND( val, 1,  MPI_DOUBLE_COMPLEX, torank, MPI_ANY_TAG,mpic%comm, ierror)
       end if
    else
       arrayval(mpic%rank)=val
       call mpic%sum(arrayval)
    end if
#else
    arrayval(0)=val
#endif
    return


  end subroutine build_cdp
  subroutine build_l(mpic,val,arrayval,torank)
    implicit none
    class(mpi_communicator) :: mpic
    integer,optional,intent(in)::torank

    logical,allocatable :: arrayval(:)
    logical,intent(in)::val
    logical::valp

    integer::iproc,proc_source,ierror,sourceout,statut

#ifdef PARA    
    if (present (torank)) then
       if (mpic%rank==torank) then
          do iproc=0,mpic%rank-1
             if (iproc==torank)then 
                arrayval(torank)=val
             else
                call mpic%probe(999,sourceout=proc_source)
                !                call mpic%recv(valp,proc_source,999)
                call MPI_RECV( valp, 1, MPI_LOGICAL, sourceout, 118,mpic%comm, status,ierror)
                arrayval(proc_source)=valp
             end if
          end do
       else
          call mpic%send(val,torank,999)
       end if
    else
       arrayval(mpic%rank)=val
!       call mpic%bcast(mpic%rank,arrayval(mpic%rank))
       call MPI_SEND( val, 1,  MPI_LOGICAL, torank, 118,mpic%comm, ierror)
    end if
#else
    arrayval(0)=val
#endif
    return


  end subroutine build_l
end module Tpara
