!--------------------------------------------------------------------------------------------------
!------Subroutines for histograms and reconstruction of free energy from meanforce observed. 
!--------------------------------------------------------------------------------------------------

subroutine fill_histo()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: 
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY: dcsi,icsi, delta_z,nhisto,nhisto1,nhisto2,& 
                              histo,histo1,histo2, &
                              histo_temp,histo_temp1
implicit none
histo_temp(:)=0
histo_temp1(:)=0
 !The border zone. Over this border the mean force is zero
 if ((icsi > -nhisto2).and.(icsi<nhisto+nhisto2)) then
   histo2(icsi)=histo2(icsi) + 1
  ! The intermediate region 
   if ((icsi > -nhisto1).and.(icsi<nhisto+nhisto1)) then
    histo1(icsi)=histo1(icsi) + 1
    histo_temp1(icsi)=histo_temp1(icsi) + 1
    ! The inner regin  
     if ((icsi > 0).and.(icsi<=nhisto)) then
       histo(icsi)=histo(icsi) + 1
       histo_temp(icsi)=histo_temp(icsi) + 1
     end if
   end if  

 end if

return
end subroutine fill_histo


subroutine Free_energy_ABF()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: erg2eV,A2cm
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY: delta_z,nhisto,nhisto1,nhisto2,& 
                              histo,histo1,histo2,Free_energy,&
                              mean_force1,abf_type,x_mol,temperature,&
                              A_ee,A_dev_ee
implicit none
integer::i_loop
real(double)::Free_temp(-nhisto2:nhisto+nhisto2)
real(double)::renorm_f,minfreeval
Free_temp(:)=0

if (abf_type .NE. 5) then ! pour ABFee, on va calculer autrement l'énergie libre 
 Free_energy(-nhisto1)=0.d0
 do i_loop=-nhisto1+1,nhisto+nhisto1
  Free_energy(i_loop)=Free_energy(i_loop-1)+0.5d0*delta_z*(mean_force1(i_loop-1)+mean_force1(i_loop))
 enddo


 forall(i_loop=-nhisto1:nhisto+nhisto1) Free_temp(i_loop)=exp(-free_energy(i_loop)/temperature)
 renorm_f=temperature*log(sum(Free_temp)*delta_z)
 forall(i_loop=-nhisto1:nhisto+nhisto1) Free_energy(i_loop)=Free_energy(i_loop)+renorm_f

 minfreeval=minval(Free_energy(:))

endif


 select case (abf_type)
 case(1)
  write(*,*),'Free energy computation....Langevin Dynamics'
  open(unit=992,file='Free_energy_Langevin',status='unknown')
   do i_loop=-nhisto1+1,nhisto+nhisto1
    write(992,*), x_mol(i_loop)/A2cm,(Free_energy(i_loop)-minfreeval)*erg2eV
   enddo
  close(992)

 case(2)
  write(*,*),'Free energy computation....ABF BIN'
  open(unit=993,file='Free_energy_ABFBIN',status='unknown')
   do i_loop=-nhisto1+1,nhisto+nhisto1
    write(993,*), x_mol(i_loop)/A2cm,(Free_energy(i_loop)-minfreeval)*erg2eV
   enddo
  close(993)

 case(3)
  write(*,*),'Free energy computation....ABF BIN OMEGA'
  open(unit=994,file='Free_energy_ABFBIN_OMEGA',status='unknown')
   do i_loop=-nhisto1+1,nhisto+nhisto1
    write(994,*), x_mol(i_loop)/A2cm,(Free_energy(i_loop)-minfreeval)*erg2eV
   enddo
  close(994)

 case(4)
  write(*,*),'Free energy computation....ABF Gaussian'
  open(unit=995,file='Free_energy_ABFGaussian',status='unknown')
   do i_loop=-nhisto1+1,nhisto+nhisto1
    write(995,*), x_mol(i_loop)/A2cm,(Free_energy(i_loop)-minfreeval)*erg2eV
   enddo
 close(995)

 case(5)
  write(*,*),'Free energy computation....ABF ee'
  !--------Pour ABFee, on a besoin des A_ee pour calculer l'énergie libre.
 forall(i_loop=-nhisto2:nhisto+nhisto2) Free_temp(i_loop)=exp(-A_ee(i_loop)/temperature)
 
 renorm_f=temperature*log(sum(Free_temp)*delta_z)

 forall(i_loop=-nhisto2:nhisto+nhisto2) Free_energy(i_loop)=A_ee(i_loop)+renorm_f

 minfreeval=minval(Free_energy(:))
 
open(unit=996,file='Free_energy_ABFee',status='unknown')
   do i_loop=-nhisto2+1,nhisto+nhisto2
    write(996,*), x_mol(i_loop)/A2cm,(Free_energy(i_loop)-minfreeval)*erg2eV
   enddo
  close(996)

end select

end subroutine Free_energy_ABF


