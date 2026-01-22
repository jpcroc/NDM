module Tpara_io
  
  use T_kind_param_m
  use Tpara, only: mpi_communicator, endmpi
#ifdef PARA
  use mpi
#endif

implicit none

#ifdef PARA
  ! Deux possibilités principales : "native" et "external32"
  character(len=*), parameter :: NDM_MPI_DATA_REPRESENTATIONS = "external32"
#endif

  interface file_write
     module procedure file_write_i
     module procedure file_write_dp
     module procedure file_write_char
  end interface

  interface file_write_all
     module procedure file_write_all_i
     module procedure file_write_all_dp
     module procedure file_write_all_char
  end interface

  interface file_read
     module procedure file_read_i
     module procedure file_read_dp
     module procedure file_read_char
  end interface

  interface file_read_all
     module procedure file_read_all_i
     module procedure file_read_all_dp
     module procedure file_read_all_char
  end interface

contains



  !=========================================================================
  integer function type_size(type)

    integer, intent(in) :: type
    !=====
    integer :: ierror=0
    !=====
    ! TO DO: ajouter type_size sans mpi

#if defined(PARA)
    call MPI_Type_size(type, type_size, ierror)
    call error_check(ierror)
#endif
    return
  end function type_size


  !=========================================================================
  subroutine mpic_file_open(mpic,file,fh)
    class(mpi_communicator),intent(in) :: mpic
    character, intent(in) :: file*80
    integer, intent(inout) :: fh
    !=====
    integer :: ierror=0
    !=====

    ! TO DO: ajouter mode : MPI_MODE_CREATE + MPI_MODE_RDWR / status = 'unknown'

#if defined(PARA)
    call MPI_File_open(mpic%comm, file, MPI_MODE_CREATE + MPI_MODE_RDWR, MPI_INFO_NULL, fh, ierror)
    call error_check(ierror)
#else
    open(unit=fh, file=file, form='unformatted', status='unknown')
#endif
  end subroutine mpic_file_open


!=========================================================================
  subroutine file_close(fh)
    integer, intent(inout) :: fh
    !=====
    integer :: ierror=0
    !=====

#if defined(PARA)
    call MPI_File_close(fh, ierror)
    call error_check(ierror)
#else
    close(unit=fh)
#endif
  end subroutine file_close


  !=========================================================================
#if defined(PARA)
  subroutine file_set_view(fh, offset, etype, filetype)
    integer, intent(in) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(in):: offset
    integer, intent(in) :: etype, filetype
    !=====
    integer :: ierror=0
    !=====

    call MPI_File_set_view(fh, offset, etype, filetype, NDM_MPI_DATA_REPRESENTATIONS, MPI_INFO_NULL, ierror)
    call error_check(ierror)

  end subroutine file_set_view
#endif


  !=========================================================================
#if defined(PARA)
  subroutine file_write_i(fh, array)
    integer, intent(in) :: fh
    integer,intent(in) :: array(..)
    !=====
    integer :: ierror=0
    !=====

    call MPI_File_write(fh, array, size(array), MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
    call error_check(ierror)

  end subroutine file_write_i
#endif


  !=========================================================================
#if defined(PARA)
  subroutine file_write_dp(fh, array)
    integer, intent(in) :: fh
    real(double), intent(in) :: array(..)
    !=====
    integer :: ierror=0
    !=====

    call MPI_File_write(fh, array, size(array), MPI_REAL8, MPI_STATUS_IGNORE, ierror)
    call error_check(ierror)

  end subroutine file_write_dp
#endif


  !=========================================================================
#if defined(PARA)
  subroutine file_write_char(fh, array)
    integer, intent(in) :: fh
    character(*),intent(in) :: array(:)
    !=====
    integer :: ierror=0
    integer :: nsize,longueur,nsizetot
    !=====

    nsize = SIZE(array)
    longueur=len(array)
    nsizetot=nsize*longueur
    call MPI_File_write(fh, array, nsizetot, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
    call error_check(ierror)

  end subroutine file_write_char
#endif


  !=========================================================================
#if defined(PARA)
  subroutine file_write_all_i(fh, array)
    integer, intent(in) :: fh
    integer,intent(in) :: array(..)
    !=====
    integer :: ierror=0
    !=====

    call MPI_File_write_all(fh, array, size(array), MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
    call error_check(ierror)

  end subroutine file_write_all_i
#endif


  !=========================================================================
#if defined(PARA)
  subroutine file_write_all_dp(fh, array)
    integer, intent(in) :: fh
    real(double), intent(in) :: array(..)
    !=====
    integer :: ierror=0
    !=====

    call MPI_File_write_all(fh, array, size(array), MPI_REAL8, MPI_STATUS_IGNORE, ierror)
    call error_check(ierror)

  end subroutine file_write_all_dp
#endif


  !=========================================================================
#if defined(PARA)
  subroutine file_write_all_char(fh, array)
    integer, intent(in) :: fh
    character(*),intent(in) :: array(:)
    !=====
    integer :: ierror=0
    integer :: nsize,longueur,nsizetot
    !=====

    nsize = SIZE(array)
    longueur=len(array)
    nsizetot=nsize*longueur
    call MPI_File_write_all(fh, array, nsizetot, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
    call error_check(ierror)

  end subroutine file_write_all_char
#endif

  !=========================================================================
#if defined(PARA)
  subroutine file_read_i(fh, array)
    integer, intent(in) :: fh
    integer,intent(in) :: array(..)
    !=====
    integer :: ierror=0
    !=====

    call MPI_File_read(fh, array, size(array), MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
    call error_check(ierror)

  end subroutine file_read_i
#endif


  !=========================================================================
#if defined(PARA)
  subroutine file_read_dp(fh, array)
    integer, intent(in) :: fh
    real(double), intent(in) :: array(..)
    !=====
    integer :: ierror=0
    !=====

    call MPI_File_read(fh, array, size(array), MPI_REAL8, MPI_STATUS_IGNORE, ierror)
    call error_check(ierror)

  end subroutine file_read_dp
#endif


  !=========================================================================
#if defined(PARA)
  subroutine file_read_char(fh, array)
    integer, intent(in) :: fh
    character(*),intent(in) :: array(:)
    !=====
    integer :: ierror=0
    integer :: nsize,longueur,nsizetot
    !=====

    nsize = SIZE(array)
    longueur=len(array)
    nsizetot=nsize*longueur
    call MPI_File_read(fh, array, nsizetot, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
    call error_check(ierror)

  end subroutine file_read_char
#endif


  !=========================================================================
#if defined(PARA)
  subroutine file_read_all_i(fh, array)
    integer, intent(in) :: fh
    integer,intent(in) :: array(..)
    !=====
    integer :: ierror=0
    !=====

    call MPI_File_read_all(fh, array, size(array), MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
    call error_check(ierror)

  end subroutine file_read_all_i
#endif


  !=========================================================================
#if defined(PARA)
  subroutine file_read_all_dp(fh, array)
    integer, intent(in) :: fh
    real(double), intent(in) :: array(..)
    !=====
    integer :: ierror=0
    !=====

    call MPI_File_read_all(fh, array, size(array), MPI_REAL8, MPI_STATUS_IGNORE, ierror)
    call error_check(ierror)

  end subroutine file_read_all_dp
#endif


  !=========================================================================
#if defined(PARA)
  subroutine file_read_all_char(fh, array)
    integer, intent(in) :: fh
    character(*),intent(in) :: array(:)
    !=====
    integer :: ierror=0
    integer :: nsize,longueur,nsizetot
    !=====

    nsize = SIZE(array)
    longueur=len(array)
    nsizetot=nsize*longueur
    call MPI_File_read_all(fh, array, nsizetot, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
    call error_check(ierror)

  end subroutine file_read_all_char
#endif


  !=========================================================================
  subroutine error_check(errorcode)

    integer, intent(in) :: errorcode
    !=====
    integer :: resultlen
    integer :: ierror=0
    !=====

#if defined(PARA)
    character(len=MPI_MAX_ERROR_STRING) :: error_string

    if (errorcode/=0) then
      call MPI_Error_string(errorcode, error_string, resultlen, ierror)
      write(6,*) error_string, errorcode
      call endmpi
    end if
#endif

    if( ierror /= 0 ) then
       write(6,*) 'error in MPI_Error_string'
    endif

  end subroutine

end module Tpara_io