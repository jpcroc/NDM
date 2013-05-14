subroutine calfo_mab()
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
  USE mab_in_ndm_module, only: it_mab,abf_type,block
  implicit none
  integer :: it_langevin
 
  it_langevin=it_mab
! one force calculation ....

! NDM part ...

         if (itab/=0) then
          if (mod(it_langevin,itab)==0) then
           call caltabt
          endif
         endif
         if (ltabvois.and.mod(it_langevin,itetabvois)==0) call caltabi
        call calfo
        if (block) call calfoblock()

!ABF part ...

        call fill_histo()
        
        select case (abf_type)
           case (1)
                  continue
           case (2)  
                  call calfo_ABF_BIN
           case (3) 
                  call calfo_ABF_BIN_OMEGA
          ! case (4)  
          !        call calfo_ABF_GAUSSIAN ! ABF Gaussian
          ! case (5) 
          !        call calfo_ABF_EE
        end select 

return
end subroutine calfo_mab



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

subroutine calfo_ABF_BIN_OMEGA()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:dcsi,icsi,rfilac,histo,histo_temp,     &
                             mean_force,cumul_force,nhisto,nhisto1, &
                             mean_force1,histo_temp1,histo1,omega_abf
  
 real(double), dimension(3,imm) :: fpabf
 real(double) :: force,omega


 
 fpabf(:,:) = zero
 ! Computing the forces from the ABF bins ...
 force = - DOT_PRODUCT(fp(:,7),rfilac(:))

 if ((icsi>-nhisto1).and.(icsi<nhisto+nhisto1)) then 
  cumul_force =  cumul_force + histo_temp1(icsi)*force
  mean_force1 (icsi) = cumul_force/(1.d0/omega_abf+histo1(icsi))
  fpabf(1:3,7)=rfilac(1:3)*mean_force1(icsi)
 end if

! Updating the forces ...
!write (41,*) fp(1,7),fpabf(1,7)
!write (42,*) fp(2,7),fpabf(2,7)
!write (43,*) fp(3,7),fpabf(3,7)

  fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)
  
return
end subroutine calfo_ABF_BIN_OMEGA
