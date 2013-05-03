subroutine calfo_ABF_BIN()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:dcsi,icsi,rfilac,histo,histo_temp,     &
                             mean_force,cumul_force,nhisto,nhisto1, &
                             mean_force1,histo_temp1,histo1
  
 real(double), dimension(3,imm) :: fpabf
 real(double) :: force

 fpabf(:,:) = zero
! Computing the forces from the ABF bins ...
force = - DOT_PRODUCT(fp(:,7),rfilac(:))

if ((icsi>-nhisto1).and.(icsi<nhisto+nhisto1)) then 
cumul_force =  cumul_force + histo_temp1(icsi)*force
mean_force1 (icsi) = cumul_force/histo1(icsi)
fpabf(1:3,7)=rfilac(1:3)*mean_force1(icsi)

end if

! Updating the forces ...
!write (41,*) fp(1,7),fpabf(1,7)
!write (42,*) fp(2,7),fpabf(2,7)
!write (43,*) fp(3,7),fpabf(3,7)

  fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)
  

return
end subroutine calfo_ABF_BIN
