module mab_in_ndm_module
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use var_pot,  ONLY: cm
      use jqmod
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      implicit none

       
      integer                                       :: it_phondy
      integer :: nmat
      real(double),dimension(:),allocatable,save    :: tmass_mab
      real(double) ::  convert_phondy
      real(double),save :: epot0
      real(double),dimension(:),allocatable, save:: w
      real(double)   :: avogadro, electron,two_pi,unit_nu
 contains
 

subroutine allocate_mab()

   implicit none
   
   nmat=3*im   
   allocate( tmass_mab(nmat))
!  pi=4*datan(1.D0)
   two_pi=2.0d0*4.d0*datan(1.d0)
   avogadro=6.0221367
   electron=1.60217733
   !hplanck=6.62618
   unit_nu=dsqrt(avogadro*electron/1.D3)
   unit_nu=unit_nu/two_pi         




return
end   subroutine allocate_mab


 subroutine  init_mab_in_ndm_module               &
             (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
!
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: ielat(imm)
    integer  :: iwmax(imm)
    integer  :: ityp(imm)
    real(double)  :: xp(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: ax(3,imm)
    real(double)  :: fp(3,imm)
 !
    integer   :: ic_local,iatom
     

   it_phondy=0
   
   
   
   do ic_local=1,3*im
      iatom = MOD(ic_local,im)
      if (iatom==0) iatom=im       
      tmass_mab(ic_local) = (cm(ityp(iatom))/umass)
   end do
   

 return
!
end subroutine init_mab_in_ndm_module 
 

end module mab_in_ndm_module
