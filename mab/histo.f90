!--------------------------------------------------------------------------------------------------
!------Subroutines for histograms. 
!--------------------------------------------------------------------------------------------------

subroutine fill_histo()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: 
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY: dcsi,icsi, delta_z,nhisto,nhisto1,nhisto2,& 
                              histo,histo1,histo2
implicit none
! ....-nhisto2......-nhisto1....0.........nhisto....nhisto1......nhisto2......

 !The border zone. Over this border the mean force is zero
 if ((icsi >= -nhisto2).and.(icsi <= nhisto+nhisto2)) then
   histo2(icsi)=histo2(icsi) + 1
  ! The intermediate region 
   if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then
    histo1(icsi)=histo1(icsi) + 1
    ! The inner regin  
     if ((icsi > 0).and.(icsi <= nhisto)) then
       histo(icsi)=histo(icsi) + 1
     end if
   end if  

 end if
return
end subroutine fill_histo

subroutine fill_histo_equilibre()! we fill histogram only if it is after n_equilibre steps ( the biais is considered stable) 
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: 
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY: icsi,nhisto,nhisto1,nhisto2,& 
                              histo_equi,histo_xi
                   
implicit none
  
! The intermediate region 
   if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then
    histo_xi(icsi)=histo_xi(icsi) + 1
   end if  
return
end subroutine fill_histo_equilibre


subroutine fill_final_histo()!To create the final histogram for the simulation
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: erg2eV,A2cm
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY: nhisto,nhisto1,nhisto2,& 
                              histo,histo1,histo2,delta_z,&
                              histo_xi,histo_equi,abf_type,abf_mode, histo_zeta,&
                              x_mol

implicit none
integer::i_iter
real(double)::sum_histo1,sum_histo,sum_histo_xi,sum_histo_zeta

if ((abf_type .ne. 6) .and. (abf_type .ne. 7)) then ! when the biais is updated 

 open(unit=989,file='histogram1',status='unknown')
 open(unit=990,file='histogram',status='unknown')

  sum_histo1=sum(histo1)
 sum_histo=sum(histo)

   do i_iter=-nhisto1,nhisto+nhisto1
    if (abf_mode==1) write(989,*),x_mol(i_iter)/A2cm, histo1(i_iter)/sum_histo1
    if (abf_mode==2) write(989,*),dble(i_iter)*delta_z, histo1(i_iter)/sum_histo1
   enddo
   
   do i_iter=1,nhisto
    if (abf_mode==1) write(990,*),x_mol(i_iter)/A2cm, histo(i_iter)/sum_histo
    if (abf_mode==2) write(990,*),dble(i_iter)*delta_z, histo(i_iter)/sum_histo
   enddo

  close(989)
  close(990)

 if (histo_equi == .True.) then
 open(unit=970,file='histogram_xi',status='unknown')  ! this file is available only if we only fill histogram after n_equilibre. 
 sum_histo_xi=sum(histo_xi)

    do i_iter=-nhisto1,nhisto+nhisto1
     if (abf_mode==1) write(970,*),x_mol(i_iter)/A2cm, histo_xi(i_iter)/sum_histo_xi
     if (abf_mode==2) write(970,*),dble(i_iter)*delta_z, histo_xi(i_iter)/sum_histo_xi
    enddo

  close(970)
 endif


endif


if (abf_type == 5 .or. abf_type == 8 ) then ! if it is ABFee, we need to fill histogram on zeta, for the case when the biais is updated or fixed

open(unit=969,file='histogram_zeta',status='unknown')

sum_histo_zeta=sum(histo_zeta)
  
 do i_iter=-nhisto1,nhisto+nhisto1
    if (abf_mode==1) write(969,*),x_mol(i_iter)/A2cm, histo_zeta(i_iter)/sum_histo_zeta
    if (abf_mode==2) write(969,*),dble(i_iter)*delta_z, histo_zeta(i_iter)/sum_histo_zeta
 enddo


 close(969)

endif


if (abf_type== 6) then
open(unit=1105,file='histogram_zeta_const_biais',status='unknown')
open(unit=1106,file='histogram_xi_const_biais',status='unknown')
  sum_histo_xi=sum(histo1)
 sum_histo_zeta=sum(histo_zeta)
  
 do i_iter=-nhisto1,nhisto+nhisto1
    if (abf_mode==1) write(1105,*),x_mol(i_iter)/A2cm, histo_zeta(i_iter)/sum_histo_zeta
    if (abf_mode==2) write(1105,*),dble(i_iter)*delta_z, histo_zeta(i_iter)/sum_histo_zeta
    if (abf_mode==1) write(1106,*),x_mol(i_iter)/A2cm, histo1(i_iter)/sum_histo_xi
    if (abf_mode==2) write(1106,*),dble(i_iter)*delta_z, histo1(i_iter)/sum_histo_xi
 enddo



 close(1105)
 close(1106)
endif


if (abf_type == 7) then

 open(unit=1107,file='histogram_ABFbin_const_biais',status='unknown')
 
 sum_histo_xi=sum(histo1)
  
 do i_iter=-nhisto1,nhisto+nhisto1
   if (abf_mode==1)  write(1107,*),x_mol(i_iter)/A2cm, histo1(i_iter)/sum_histo_xi
   if (abf_mode==2)  write(1107,*),dble(i_iter)*delta_z, histo1(i_iter)/sum_histo_xi
 enddo



 close(1107)
 

endif

end subroutine fill_final_histo



