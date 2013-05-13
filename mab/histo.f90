subroutine fill_histo()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: 
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY: dcsi,icsi, delta_z,nhisto,nhisto1,nhisto2,& 
                              histo,histo1,histo2, &
                              histo_temp,histo_temp1
 
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


