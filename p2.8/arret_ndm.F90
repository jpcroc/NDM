module arret_ndm_mod
#ifdef PARA
    USE mpi
!    use TPara,only: MPI_COMM_space,status,ierr,myidsp,NDM_MPI_REAl_DOUBLE
    USE mod_para,only:MPI_COMM_space,status,ierr,myidsp,NDM_MPI_REAl_DOUBLE
    use gen_com_m ,only:rang
   
#endif
        implicit none
        contains
subroutine arret_ndm()

  USE T_kind_param_m

  !USE mpi
  implicit none
#ifdef PARA  
  include "mpif.h"
#endif

  ! Routine d'arret du code NDM

  !--------------------------------------------------
  !Variables de la routine

  !--------------------------------------------------
  !Variables locales

  !--------------------------------------------------
  !Corps de la routine

#ifdef PARA
  call MPI_FINALIZE(ierr)
#endif

  stop

end subroutine arret_ndm

end module
