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
      integer :: rangmab=0
      real(double), dimension(:,:), allocatable :: sig_i,sig_ll,rga_i,xp0,m_i
      real(double)  :: dtlang,dtlang_ini, temperature,Ecinetique,m_tot,a0bcc,omega_abf,maxforce,lang_factor
      

      real(double), parameter :: KtoERG=1.3791946308724831d-16 
      real(double), parameter :: THZtoK=47.99407332506204d0
      real(double), parameter :: unit_omega_to_erg=1.d+24*(2.d0*pi)**2!*umass, umass is not included because the 
                                                                      !cm are already in  multiplied by 
                                                                      ! umass (in g) in the main NDM program. 

      real(double) :: crit=0.01/1d+8
      integer :: nlangevin,abf_type,sim_mode,langevin_type,n_equilibre,abf_mode,mode_zeta_potential
      integer :: it_mab,it_stop,itest_stop,it_calc_brute
      integer :: ntestvacancyjump
      real(double),save :: epot0,gamma
      real(double),dimension(:),allocatable, save:: w
      real(double)   :: pinumber,dcsi,normxlac,deltasph,radiussph,rtestlac
      real(double),dimension(3) :: xbar,xbarini,xlaci,xlacf,rfilac 
      real(double)  :: deltar1,deltar2,delta_z,rests
      real(double)  :: xi_min,xi_max,alpha_zeta
      integer       :: nhisto,nhisto1,nhisto2,icsi,nwrite_histo
      real(double),dimension(:), allocatable :: histo,histo1,histo2
      real(double),dimension(:), allocatable :: histo_xi, histo_xi1, histo_xi2
      real(double),dimension(:),allocatable::histo_zeta
      real(double),dimension(:), allocatable :: mean_force,mean_force1,mean_force2,mean_force_ABFee
      real(double),dimension(:),allocatable::x_mol,cumul_force_denom1,cumul_force1
      real(double),dimension(:),allocatable::Free_energy
      real(double),dimension(:), allocatable :: unit_histo2
      real(double),dimension(:),allocatable::A_dev_ee,A_ee,P_ee,P_ee_num,P_ee_denom,A_bar_ee
      real(double),dimension(:),allocatable::exp_A_bar
      real(double),dimension(:),allocatable::A_theo,error_A,error_A_bar
      real(double) :: limit1m, limit1p, limit2m, limit2p,limit1,limit2
      real(double) :: Free_energy_brute
      real(double) :: omega_einstein,ene_einstein,ene0,einstein_free_3N, einstein_correction, pbc_correction
      real(double), dimension(:,:), allocatable :: omega_veinstein,fpeinstein
      real(double) , dimension(:,:,:,:) , allocatable :: matfor
      integer :: it_en, neq_lang
      real(double)::sigma_eta,sigma_carre,eta_ABFee,sum_error_A,sum_error_A_bar
      integer::ecart_eta,nom_deconvo
      real(double)::eta_mab,ha_mix,temperature_zeta_min,temperature_zeta_max,equit,units_phondy
      integer::compute_mode,error_step
      integer:: nsite_block,atom_to_jump,itype_reaction,itype_einstein
      integer, dimension(:), allocatable :: isite_block
      logical :: block,block_file, test_end,histo_equi
      integer,parameter :: abf_mode_reaction=1,  &
                           abf_mode_alchemical=2,&
                           abf_mode_temperature=22
      !neb part:
      integer :: nimage_neb, nimage_lambda

 contains
 

subroutine allocate_mab()

   implicit none
   
   allocate (sig_i(3,imm),sig_ll(3,imm),rga_i(3,imm),m_i(3,imm),xp0(3,imm)) 
   allocate (fpeinstein(3,imm))  
  return
end   subroutine allocate_mab


 subroutine  init_mab_in_ndm_module()
    use tab_imm_m
    implicit none
    integer :: ic
     
     dtlang_ini=dtlang
     m_tot=SUM(m_i(1,1:im))
    
     xp0(:,:) = xp(:,:)
     do ic=1,3
       xbarini(ic)  = sum(xp(ic,1:im)*m_i(ic,1:im))/m_tot
     enddo
 
     pinumber=4.d0*datan(1.D0)

    it_en=-1
    it_stop=0
    itest_stop=0
 
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
     nhisto1=deltar1/dble(delta_z)
     nhisto2=deltar2/dble(delta_z)
   end if 
!set-up the alchemical case ...

   if (abf_mode==2) then
    xi_min=0.d0
    xi_max=1.d0
    normxlac=xi_max-xi_min
    delta_z=normxlac/dble(nhisto)
   end if 



   if (abf_mode==22) then
                                            ! 500 is the temperature in order to made the simulation
    xi_max=temperature/temperature_zeta_min ! this is xi_max in order to have 100K at the reference temperature 500 K 
    xi_min=temperature/temperature_zeta_max ! this is xi_min in order to have melting temperature;  500 K is the reference temperature
    normxlac=xi_max-xi_min
    delta_z=normxlac/dble(nhisto)
   end if 
   !
   if ((abf_mode==2).or.(abf_mode==22)) then
    deltar1=deltar1*normxlac  
    deltar2=deltar2*normxlac  
    !
    nhisto1=deltar1/dble(delta_z)
    nhisto2=deltar2/dble(delta_z)
    !
    limit1m=xi_min-deltar1
    limit1p=xi_max+deltar1
    limit2m=xi_min-deltar2
    limit2p=xi_max+deltar2
    !
     if (mode_zeta_potential==0) then
      limit1=limit1m
      limit2=limit1p
     else if (mode_zeta_potential==1) then 
      limit1=limit2m
      limit2=limit2p
     end if 
   end if  
   !
    if (abf_mode==22) equit=0.d0
     equit=0.d0

    write(*,*) 'Equit correction ', equit*erg2ev
    !sigma_eta=sqrt(eta_mab)
     sigma_eta=eta_mab*delta_z ! choisir largeur de gaussienne
     sigma_carre=sigma_eta**2
     ecart_eta=nint(3.d0*sigma_eta/delta_z)
      write(*,*),'ecart_eta,delta_z,sigma_eta',ecart_eta,delta_z,sigma_eta
    write(*,'("The distribution over the histogram............:")')

    write(*,'("...-nhisto2=",i6,"...-nhisto1=",i6,"..0..........nhisto=",i6,"......nhisto+nhisto1=",i6, &
          "....nhisto+nhisto2=",i6,"...")') -nhisto2, -nhisto1,nhisto,nhisto+nhisto1,nhisto+nhisto2
     write(*,'("...-nhisto2=",f6.2,"...-nhisto1=",f6.2,"..",f6.2,".....nhisto=",f6.2,"......nhisto+nhisto1=",f6.2, &
          "....nhisto+nhisto2=",f6.2,"...")') -deltar2, -deltar1,xi_min,xi_max,xi_max+deltar1,xi_max+deltar2




  

    


    allocate(histo(0:nhisto),histo1(-nhisto1:nhisto+nhisto1),histo2(-nhisto2:nhisto+nhisto2))
    allocate(histo_xi(0:nhisto),histo_xi1(-nhisto1:nhisto+nhisto1),histo_xi2(-nhisto2:nhisto+nhisto2))
    allocate(histo_zeta(-nhisto1:nhisto+nhisto1))
    allocate(mean_force(0:nhisto),mean_force1(-nhisto1:nhisto+nhisto1),mean_force2(-nhisto2:nhisto+nhisto2))
    allocate(mean_force_ABFee(-nhisto2:nhisto+nhisto2))
    allocate(cumul_force1(-nhisto1:nhisto+nhisto1),cumul_force_denom1(-nhisto1:nhisto+nhisto1))
    allocate(x_mol(-nhisto2:nhisto+nhisto2))
    allocate(Free_energy(-nhisto2:nhisto+nhisto2))
    allocate(A_theo(-nhisto1:nhisto+nhisto1),error_A(-nhisto1:nhisto+nhisto1))
    allocate(error_A_bar(-nhisto1:nhisto+nhisto1))
    forall(ic=-nhisto2:nhisto+nhisto2) x_mol(ic)=xi_min+ic*delta_z 
    allocate (unit_histo2(-nhisto2:nhisto+nhisto2))
    
    if (abf_mode==1) unit_histo2(:)=x_mol(:)/A2cm
    if ((abf_mode==2).or.(abf_mode==22)) then 
      do ic=-nhisto2, nhisto+nhisto2
       unit_histo2(ic)=xi_min+delta_z*dble(ic)
     end do
    end if 
    cumul_force1(:)=0.d0 
    cumul_force_denom1(:)=0.d0

    histo(0:nhisto)=0
    histo1(-nhisto1:nhisto+nhisto1)=1
    histo2(-nhisto2:nhisto+nhisto2)=1

    histo_xi(0:nhisto)=0
    histo_xi1(-nhisto1:nhisto+nhisto1)=0
    histo_xi2(-nhisto2:nhisto+nhisto2)=0
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
