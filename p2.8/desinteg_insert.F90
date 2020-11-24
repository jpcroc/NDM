module desinteg_insert_mod
  USE caltabt_mod,only: caltabt
  USE gen_com_m, ONLY:bk,deltaespr,deltaf,erg2ev,it,itdes,itmax,nstepdes,pm1des,rang,&
       &tempdes,typspr,vpchdeb,vpchdn,vpchup,xpchdeb,xpchdn,xpchup,xpspr,xpspr0,im,&
       &imdesup,imdesdeb,imdesdn,itichdeb,itichdn,itichup,num_at_globdesdeb,num_at_globdesup,num_at_globdesdn,im_glob

!  USE atomconfig
  implicit none
contains

  subroutine desinteg_insert
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE tab_imm_m,only:xp,vp,ielat,ityp,num_at_glob
#ifdef PARA
    USE mpi
    USE mod_para,only:MPI_COMM_space,status,ierr,nprocs,myid,NDM_MPI_REAl_DOUBLE
#endif

    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    !------------------,-----------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------

    real(double),allocatable,save:: Wchemin(:)
    real(double),save :: Wch0
    real(double)::wch1,testval,u1,taup,taum,mug
    real(double), parameter::bkev=8.617385d-5 
    real(double),save::taupm=0,taumm=0

    integer:: i,k,l,m,n,ic,iaccept
    integer,save:: nchemin=0,nchacc,nchup=0,nchdn=0


    logical:: laccept


    return
  end subroutine desinteg_insert
end module desinteg_insert_mod
