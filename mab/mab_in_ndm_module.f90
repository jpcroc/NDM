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

      real(double), dimension(:,:), allocatable :: sig_i,sig_ll,rga_i,xp0,m_i
      real(double)  :: dtlang,temperature,Ecinetique,m_tot,a0bcc,omega_abf
      

      real(double), parameter :: KtoERG=1.3791946308724831d-16 
      integer :: nlangevin,abf_type,sim_mode,langevin_type 
      integer                                       :: it_mab
      real(double),save :: epot0,gamma
      real(double),dimension(:),allocatable, save:: w
      real(double)   :: pinumber,dcsi,normxlac,deltasph,radiussph,rtestlac
      real(double),dimension(3) :: xbar,xbarini,xlaci,xlacf,rfilac 
      real(double)  :: deltar1,deltar2,delta_z
      real(double)  :: xi_min,xi_max,a_Fermi,Fermi_percent
      integer       :: nhisto,nhisto1,nhisto2,icsi,nwrite_histo
      integer,dimension(:), allocatable :: histo,histo1,histo2,histo_temp,histo_temp1
      real(double),dimension(:), allocatable :: mean_force,mean_force1
      real(double),dimension(:),allocatable::x_mol,cumul_force_denom1,cumul_force1
      real(double),dimension(:),allocatable::Free_energy
      real(double)::sigma_eta
      integer::ecart_eta
      real(double)::eta_mab
       
      
      logical :: block,test_end

 contains
 

subroutine allocate_mab()

   implicit none
   
   allocate (sig_i(3,imm),sig_ll(3,imm),rga_i(3,imm),m_i(3,imm),xp0(3,imm)) 
  

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
   
    normxlac=sqrt(SUM((xlacf(1:3)-xlaci(1:3))**2))
    rfilac(1:3)=(xlacf(1:3)-xlaci(1:3))/normxlac
  
        
    a_Fermi=a_Fermi/normxlac 
    xi_min=-normxlac*Fermi_percent
    xi_max=(1.d0+Fermi_percent)*normxlac

    delta_z=normxlac/dble(nhisto)

    !sigma_eta=sqrt(eta_mab)
    sigma_eta=eta_mab*delta_z ! choisir largeur de gaussienne
    
    ecart_eta=nint(3.d0*sigma_eta/delta_z)
    write(*,*),'ecart_eta,delta_z,sigma_eta',ecart_eta,delta_z,sigma_eta
    nhisto1=deltar1/dble(delta_z)
    nhisto2=deltar2/dble(delta_z)
    allocate(histo(nhisto),histo1(-nhisto1:nhisto+nhisto1),histo2(-nhisto2:nhisto+nhisto2),& 
             histo_temp(nhisto),histo_temp1(-nhisto1:nhisto+nhisto1))
    allocate(mean_force(nhisto),mean_force1(-nhisto1:nhisto+nhisto1))
    allocate(cumul_force1(-nhisto1:nhisto+nhisto1),cumul_force_denom1(-nhisto1:nhisto+nhisto1))
    allocate(x_mol(-nhisto1:nhisto+nhisto1))
    allocate(Free_energy(-nhisto1:nhisto+nhisto1))
    forall(ic=-nhisto1:nhisto+nhisto1) x_mol(ic)=ic*delta_z 
    cumul_force1(:)=0.d0 
    cumul_force_denom1(:)=0.d0
    histo(1:nhisto)=0
    histo_temp(1:nhisto)=0
    histo1(-nhisto1:nhisto+nhisto1)=0
    histo2(-nhisto2:nhisto+nhisto2)=0
   
    
    

 return
!
end subroutine init_mab_in_ndm_module 

end module mab_in_ndm_module
