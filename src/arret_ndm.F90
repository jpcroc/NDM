module arret_ndm_mod
#ifdef PARA
    USE Tpara,only:endmpi
!    USE mod_para,only:MPI_COMM_space,status,ierr,myidsp,NDM_MPI_REAl_DOUBLE
    use gen_com_m ,only:uwrt,lwrt,rang
   
#endif
        implicit none
        contains
subroutine arret_ndm(lforcestop)

  USE T_kind_param_m
  implicit none
  logical,optional::lforcestop
  logical::lfstp

  lfstp=.false.
  if (present(lforcestop)) lfstp=lforcestop

  ! Routine d'arret du code NDM


#ifdef PARA
  call endMPI(lfstp)
#endif
  stop
end subroutine arret_ndm

end module arret_ndm_mod
