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
      real(double)  :: dtlang,temperature,Ecinetique,m_tot,a0bcc,omega_abf,maxforce
      

      real(double), parameter :: KtoERG=1.3791946308724831d-16 
      real(double), parameter :: THZtoK=47.99407332506204d0
      real(double), parameter :: unit_omega_to_erg=1.d+24*(2.d0*pi)**2!*umass, umass is not included because the 
                                                                      !cm are already in  multiplied by 
                                                                      ! umass (in g) in the main NDM program. 

      integer :: nlangevin,abf_type,sim_mode,langevin_type,n_equilibre,abf_mode
      integer                                       :: it_mab
      real(double),save :: epot0,gamma
      real(double),dimension(:),allocatable, save:: w
      real(double)   :: pinumber,dcsi,normxlac,deltasph,radiussph,rtestlac
      real(double),dimension(3) :: xbar,xbarini,xlaci,xlacf,rfilac 
      real(double)  :: deltar1,deltar2,delta_z
      real(double)  :: xi_min,xi_max
      integer       :: nhisto,nhisto1,nhisto2,icsi,nwrite_histo
      real(double),dimension(:), allocatable :: histo,histo1,histo2,histo_temp,histo_temp1,histo_xi
      real(double),dimension(:),allocatable::histo_zeta
      real(double),dimension(:), allocatable :: mean_force,mean_force1,mean_force2,mean_force_ABFee
      real(double),dimension(:),allocatable::x_mol,cumul_force_denom1,cumul_force1
      real(double),dimension(:),allocatable::Free_energy
      real(double),dimension(:),allocatable::A_dev_ee,A_ee,P_ee,P_ee_num,P_ee_denom,A_bar_ee
      real(double),dimension(:),allocatable::exp_A_bar
      real(double),dimension(:),allocatable::A_theo,error_A,error_A_bar

      real(double) :: omega_einstein,ene_einstein,ene0,einstein_free_3N, einstein_correction, pbc_correction
      real(double), dimension(:,:), allocatable :: omega_veinstein,fpeinstein
      integer :: it_en
      real(double)::sigma_eta,sigma_carre,eta_ABFee,sum_error_A,sum_error_A_bar
      integer::ecart_eta,nom_deconvo
      real(double)::eta_mab
      integer::compute_mode,error_step
      
      logical :: block,test_end,histo_equi


 contains
 

subroutine allocate_mab()

   implicit none
   
   allocate (sig_i(3,imm),sig_ll(3,imm),rga_i(3,imm),m_i(3,imm),xp0(3,imm)) 
   allocate (omega_veinstein(3,imm),fpeinstein(3,imm))  

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
    
     xp0(:,:) = xp(:,:)
     do ic=1,3
       xbarini(ic)  = sum(xp(ic,1:im)*m_i(ic,1:im))/m_tot
     enddo
 
     pinumber=4.d0*datan(1.D0)

    it_en=-1
    
! set-up the reaction coordinate case
   if (abf_mode==1) then
     xlaci(1:3)=(/a0bcc,a0bcc,a0bcc/)/angst
     xlacf(1:3)=(/a0bcc/2.d0,a0bcc/2.d0,a0bcc/2.d0 /)/angst
     normxlac=sqrt(SUM((xlacf(1:3)-xlaci(1:3))**2))
     rfilac(1:3)=(xlacf(1:3)-xlaci(1:3))/normxlac
     delta_z=normxlac/dble(nhisto)
     deltar1=dsqrt(3.d0)*a0bcc*deltar1/(angst*2.d0)
     deltar2=dsqrt(3.d0)*a0bcc*deltar2/(angst*2.d0)
     rtestlac=dsqrt(3.d0)*a0bcc*rtestlac/(angst*2.d0)
   end if 
!set-up the alchemical case ...

   if (abf_mode==2) then
    xi_min=0.d0
    xi_max=1.d0
    normxlac=xi_max-xi_min
    delta_z=normxlac/dble(nhisto)
   end if 

    !sigma_eta=sqrt(eta_mab)
     sigma_eta=eta_mab*delta_z ! choisir largeur de gaussienne
     sigma_carre=sigma_eta**2
     ecart_eta=nint(3.d0*sigma_eta/delta_z)
      write(*,*),'ecart_eta,delta_z,sigma_eta',ecart_eta,delta_z,sigma_eta
     nhisto1=deltar1/dble(delta_z)
     nhisto2=deltar2/dble(delta_z)
 
    if (abf_mode==2) then
        omega_veinstein(:,:)=omega_einstein
        !instead that I will a file with all the einstein  frequencies 
    end if 


   
    allocate(histo(nhisto),histo1(-nhisto1:nhisto+nhisto1),histo2(-nhisto2:nhisto+nhisto2),& 
             histo_temp(nhisto),histo_temp1(-nhisto1:nhisto+nhisto1))
    allocate(histo_xi(-nhisto1:nhisto+nhisto1),histo_zeta(-nhisto1:nhisto+nhisto1))
    allocate(mean_force(nhisto),mean_force1(-nhisto1:nhisto+nhisto1),mean_force2(-nhisto2:nhisto+nhisto2))
    allocate(mean_force_ABFee(-nhisto2:nhisto+nhisto2))
    allocate(cumul_force1(-nhisto1:nhisto+nhisto1),cumul_force_denom1(-nhisto1:nhisto+nhisto1))
    allocate(x_mol(-nhisto2:nhisto+nhisto2))
    allocate(Free_energy(-nhisto2:nhisto+nhisto2))
    allocate(A_theo(-nhisto1:nhisto+nhisto1),error_A(-nhisto1:nhisto+nhisto1))
    allocate(error_A_bar(-nhisto1:nhisto+nhisto1))
    forall(ic=-nhisto2:nhisto+nhisto2) x_mol(ic)=ic*delta_z 
    cumul_force1(:)=0.d0 
    cumul_force_denom1(:)=0.d0
    histo(1:nhisto)=0
    histo_temp(1:nhisto)=0
    histo1(-nhisto1:nhisto+nhisto1)=0
    histo2(-nhisto2:nhisto+nhisto2)=0
    histo_xi(-nhisto1:nhisto+nhisto1)=0
    histo_zeta(-nhisto1:nhisto+nhisto1)=0.d0
    Free_energy(:)=0
    mean_force(:)=0
    mean_force1(:)=0
    mean_force2(:)=0
    mean_force_ABFee(:)=0
    allocate(A_ee(-nhisto2:nhisto+nhisto2),A_dev_ee(-nhisto2:nhisto+nhisto2),&
            P_ee(-nhisto2:nhisto+nhisto2),P_ee_num(-nhisto2:nhisto+nhisto2),&
            P_ee_denom(-nhisto2:nhisto+nhisto2),A_bar_ee(-nhisto2:nhisto+nhisto2),&
            exp_A_bar(-nhisto2:nhisto+nhisto2))
  A_ee(:)=0.d0
  A_theo(:)=0.d0
  error_A(:)=0.d0
  error_A_bar(:)=0.d0
  A_dev_ee(:)=0.d0
  P_ee(:)=0.d0
  P_ee_num(:)=0.d0
  P_ee_denom(:)=0.d0    
  A_bar_ee(:)=0.d0
  exp_A_bar(:)=1.d0 

 return
!
end subroutine init_mab_in_ndm_module 

end module mab_in_ndm_module
