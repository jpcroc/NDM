module art_in_ndm_module
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use jqmod
      use var_pot
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      implicit none
!      integer,       dimension(:),allocatable,save  :: ielat_c   ,ielat_i
!      integer,       dimension(:),allocatable,save  :: iwmax_c   ,iwmax_i
!      integer,       dimension(:),allocatable,save  :: ityp_c    ,ityp_i
!      real(double),dimension(:,:),allocatable,save  :: xp_c      ,xp_i
!      real(double),dimension(:,:),allocatable,save  :: xpp_c     ,xpp_i
!      real(double),dimension(:,:),allocatable,save  :: vp_c      ,vp_i
!      real(double),dimension(:,:),allocatable,save  :: ax_c      ,ax_i
!      real(double),dimension(:,:),allocatable,save  :: fp_c      ,fp_i


      integer                                       :: it_art
      real(double),dimension(:),allocatable,save    :: tmass_art
      real(double) ::  tstep_art,usdh_art, convert_art


contains
 

subroutine allocate_art()

   implicit none
   
!     allocate( ielat_i(imm),   &
!               iwmax_i(imm),   &
!               ityp_i(imm),    & 
!               ielat_c(imm),   &
!               iwmax_c(imm),   &
!               ityp_c(imm),    &
!		xp_i(3,imm),    &
!                xpp_i(3,imm),   &
!                vp_i(3,imm),    &
!                ax_i(3,imm),    &
!                fp_i(3,imm),    &
!                xp_c(3,imm),    &
!                xpp_c(3,imm),   &
!                vp_c(3,imm),    &
!                ax_c(3,imm),    &
!                fp_c(3,imm)    
   allocate( tmass_art(3*im))



  return

end   subroutine allocate_art


 subroutine  init_art_in_ndm_module               &
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
     
!       ielat_i   (:)  = ielat    (:)     
!       iwmax_i   (:)  = iwmax    (:)  
!       ityp_i    (:)  = ityp     (:)  
!       xp_i      (:,:)= xp       (:,:)
!       xpp_i     (:,:)= xpp      (:,:)
!       vp_i      (:,:)= vp       (:,:)
!       ax_i      (:,:)= ax       (:,:)
!       fp_i      (:,:)= fp       (:,:)
 !
 !
!       ielat_c   (:)  = ielat    (:)     
!       iwmax_c   (:)  = iwmax    (:)  
!       ityp_c    (:)  = ityp     (:)  
!       xp_c      (:,:)= xp       (:,:)
!       xpp_c     (:,:)= xpp      (:,:)
!       vp_c      (:,:)= vp       (:,:)
!       ax_c      (:,:)= ax       (:,:)
!       fp_c      (:,:)= fp       (:,:)
! 



   it_art=0
   
   tstep_art=tstep / utemps
   
   convert_art=9.6485d0/10000.d0
   
   do ic_local=1,3*im
      iatom = MOD(ic_local,im)
      if (iatom==0) iatom=im       
      tmass_art(ic_local) = tstep_art**2*convert_art/(cm(ityp(iatom))/umass)
   end do
   
   usdh_art= 1.d0/(2.d0*tstep_art)
   
   


 return
!
end subroutine init_art_in_ndm_module 
 
 


end module art_in_ndm_module
