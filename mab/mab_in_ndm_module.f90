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

      real(double), dimension(:,:), allocatable :: sig_i,rga_i,xp0,m_i
      real(double)  :: dtlang,temperature,Ecinetique,m_tot,a0bcc
      

      real(double), parameter :: KtoERG=1.3791946308724831d-16 
      integer :: nlangevin 
      integer                                       :: it_mab,it_langevin
      real(double),save :: epot0
      real(double),dimension(:),allocatable, save:: w
      real(double)   :: pinumber,dcsi,normxlac,deltasph,radiussph
      real(double),dimension(3) :: xbar,xbarini,xlaci,xlacf,rfilac 

 contains
 

subroutine allocate_mab()

   implicit none
   
   allocate (sig_i(3,imm),rga_i(3,imm),m_i(3,imm),xp0(3,imm)) 
  

return
end   subroutine allocate_mab


 subroutine  init_mab_in_ndm_module()
    use tab_imm_m
    implicit none
    integer :: ic

     do ic=1,3
      m_i(ic,1:im) = cm(ityp(1:im))
     end do
     m_tot=SUM(m_i(1,1:im))
    
     do ic=1,3
       xbarini(ic)  = sum(xp(ic,1:im)*m_i(ic,1:im))/m_tot
     enddo
 
     pinumber=4.d0*datan(1.D0)
    
    xp0(:,:)=xp(:,:)
    
    xlaci(1:3)=(/a0bcc,a0bcc,a0bcc/)/angst
    xlacf(1:3)=(/a0bcc/2.d0,a0bcc/2.d0,a0bcc/2.d0 /)/angst
    rfilac(1:3)=xlacf(1:3)-xlaci(1:3)
    normxlac=sqrt(SUM((rfilac(:)**2)))

    


 return
!
end subroutine init_mab_in_ndm_module 

end module mab_in_ndm_module
