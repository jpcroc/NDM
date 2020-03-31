subroutine arret_ndm()

  USE T_kind_param_m
#if(PARA)
  use mod_mpi
#endif
  implicit none
  real(double) :: temps_exe

  ! Routine d'arret du code NDM

  !--------------------------------------------------
  !Variables de la routine

  !--------------------------------------------------
  !Variables locales

  !--------------------------------------------------
  !Corps de la routine

#if(PARA)
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
