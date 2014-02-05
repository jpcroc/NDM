!-------------------------------------------------------------------
!--------------------Calculate errors for ABF at the convergence
!-----------------------------------------------------------------


subroutine error_ABF_direct()!!---------Calculate error for ABF

 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: erg2eV,A2cm
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY: nhisto,nhisto1,nhisto2,       & 
                              histo,histo1,histo2,delta_z,  &
                              histo_xi,histo_equi,abf_type, &
                              A_ee,A_bar_ee,histo_zeta,x_mol,&
                              mean_force1,mean_force2,abf_type,&
                              Free_energy,temperature,A_theo,error_A,&
                              sum_error_A



implicit none
integer::i_loop,i_iter
real(double)::Free_temp(-nhisto2:nhisto+nhisto2)
real(double)::renorm_f
Free_temp(:)=0
Free_energy(:)=0
error_A(:)=0

  Free_energy(-nhisto1)=0.d0
  do i_loop=-nhisto1+1,nhisto+nhisto1
    Free_energy(i_loop)=Free_energy(i_loop-1)+0.5d0*delta_z*(mean_force1(i_loop-1)+mean_force1(i_loop))
  enddo


 forall(i_loop=-nhisto1:nhisto+nhisto1) Free_temp(i_loop)=exp(-Free_energy(i_loop)/temperature)
 
   renorm_f=temperature*log(sum(Free_temp)*delta_z)
 
 forall(i_loop=-nhisto1:nhisto+nhisto1) Free_energy(i_loop)=Free_energy(i_loop)+renorm_f


 forall(i_iter=-nhisto1+1:nhisto+nhisto1) error_A(i_iter)=abs(Free_energy(i_iter)-A_theo(i_iter))

error_A(:)=error_A(:)/(2*nhisto1+nhisto)

sum_error_A=sum(error_A(:))

end subroutine error_ABF_direct


subroutine error_ABFee_proba_corrected()
!!---------Calculate error for ABFee with A_bar corrected by histograms on \xi, and A_tilde corrected by \zeta
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: erg2eV,A2cm
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY: nhisto,nhisto1,  & 
                              histo1,A_bar_ee,         &
                              temperature,A_theo,      &
                              error_A_bar, sum_error_A_bar,delta_z

implicit none
integer::i_loop,i_iter
real(double)::sum_histo1,temp_A_bar
real(double),dimension(:),allocatable::A_bar_corrige
real(double),dimension(:),allocatable::histo1_temp
real(double)::Free_temp(-nhisto1:nhisto+nhisto1)
allocate(A_bar_corrige(-nhisto1:nhisto+nhisto1))
allocate(histo1_temp(-nhisto1:nhisto+nhisto1))
Free_temp(:)=0
A_bar_corrige(:)=0
error_A_bar(:)=0
histo1_temp(:)=0

!----normailiser histogramme
sum_histo1=sum(histo1)
histo1_temp(:)=histo1(:)/sum_histo1

do i_iter=-nhisto1,nhisto+nhisto1

if (histo1(i_iter)==0) then

A_bar_corrige(i_iter)=A_bar_ee(i_iter)*erg2eV

else

A_bar_corrige(i_iter)=A_bar_ee(i_iter)*erg2eV-log(histo1_temp(i_iter))*temperature*erg2eV

end if

end do

forall(i_loop=-nhisto1:nhisto+nhisto1) Free_temp(i_loop)=exp(-A_bar_corrige(i_loop)/temperature)


temp_A_bar=temperature*log(sum(Free_temp)*delta_z)


A_bar_corrige(:)=A_bar_corrige(:)+temp_A_bar


forall(i_iter=-nhisto1+1:nhisto+nhisto1) error_A_bar(i_iter)=abs(A_bar_corrige(i_iter)-A_theo(i_iter))


error_A_bar(:)=error_A_bar(:)/(2*nhisto1+nhisto)


sum_error_A_bar=sum(error_A_bar(:))

end subroutine error_ABFee_proba_corrected




