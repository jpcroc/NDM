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