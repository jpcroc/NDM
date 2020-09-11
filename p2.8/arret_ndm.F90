module arret_ndm_mod
#ifdef PARA
    USE mpi
    USE mod_para,only:MPI_COMM_space,status,ierr,nprocs,myid,NDM_MPI_REAl_DOUBLE,temps_initspeed,temps_para,temps_dmloop,temps_init,temps_input,&
         &temps_config,temps_deb
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
  real(double) :: temps_exe

  ! Routine d'arret du code NDM

  !--------------------------------------------------
  !Variables de la routine

  !--------------------------------------------------
  !Variables locales

  !--------------------------------------------------
  !Corps de la routine

#ifdef PARA
!  temps_dmloop=MPI_Wtime() - temps_dmloop_deb

  temps_exe = MPI_Wtime() - temps_deb
  if (myid==0) then
     print *, 'Temps d''execution : ', temps_exe
     print *, 'Temps d''init      : ', temps_init
     print *, 'Temps d''input     : ', temps_input
     print *, 'Temps de config   : ', temps_config
     print *, 'Temps d''initspeed : ', temps_initspeed
     print *, 'Temps para estime : ', temps_para
     print *, 'Temps dmloop : ', temps_dmloop
  endif
  call MPI_FINALIZE(ierr)
#endif

  stop

end subroutine arret_ndm

end module
