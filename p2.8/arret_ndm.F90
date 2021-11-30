module arret_ndm_mod
#ifdef PARA
    USE Tpara,only:endmpi
!    USE mod_para,only:MPI_COMM_space,status,ierr,myidsp,NDM_MPI_REAl_DOUBLE
    use gen_com_m ,only:rang
   
#endif
        implicit none
        contains
subroutine arret_ndm()

  USE T_kind_param_m

  
  implicit none


  ! Routine d'arret du code NDM


#ifdef PARA
  call endMPI
#endif
  stop
end subroutine arret_ndm

end module
