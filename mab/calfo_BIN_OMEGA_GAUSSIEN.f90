!This file contains the:
! - calfo_ABF_BIN
! - calfo_ABF_BIN_OMEGA
! - calfo_ABF_BIN_GAUSSIEN 




subroutine calfo_ABF_BIN() ! this concerns only the force applied on atoms 
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev,potist
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:dcsi,icsi,rfilac,histo,     &
                             mean_force,cumul_force1,nhisto,nhisto1, &
                             mean_force1,histo1,ene_einstein,abf_mode,fpeinstein,ene0, &
                             ha_mix,equit,atom_to_jump
 implicit none

 real(double), dimension(3,imm) :: fpabf
 real(double) :: force

 fpabf(:,:) = zero
! Computing the forces from the ABF bins in the reaction coordinate case ...
if (abf_mode==1) then
 force = - DOT_PRODUCT(fp(:,atom_to_jump),rfilac(:))   ! dU(q)/dq
 !
 if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then 
  cumul_force1(icsi) =  cumul_force1(icsi) + force
  mean_force1 (icsi) = cumul_force1(icsi)/histo1(icsi)
  fpabf(1:3,atom_to_jump)=rfilac(1:3)*mean_force1(icsi)
 end if
 !
 fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)
end if 

if (abf_mode==2) then              
  force = potist - ene_einstein - ene0 !  - d U(zeta,q)/d zeta  
  fpabf(:,:) = (1.d0-dcsi)*fpeinstein(:,:) + dcsi*fp(:,:) ! -d U(zeta,q)/d q
  fp(:,:)=fpabf(:,:)
  ! write(*,*) 'fpabf', fp(:,1)
  !
  if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then 
  cumul_force1(icsi) =  cumul_force1(icsi) + force
  mean_force1 (icsi) = cumul_force1(icsi)/histo1(icsi)
  end if
  !
end if 


if (abf_mode==22) then              
  force = potist + ha_mix*ene_einstein - ene0 - equit !  d U(zeta,q)/d zeta  
  fpabf(:,:) = dcsi*ha_mix*fpeinstein(:,:) + dcsi*fp(:,:) ! -d U(zeta,q)/d q
  fp(:,:)=fpabf(:,:)
  !
  if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then 
   cumul_force1(icsi) =  cumul_force1(icsi) + force
   mean_force1 (icsi) = cumul_force1(icsi)/histo1(icsi)
  end if
  !
end if 



 !debug write(*,'("deb",i7,3D15.4,2D15.1)') icsi,dcsi, potist*erg2ev, ene_einstein*erg2ev,fpeinstein(1,1),fp(1,1)
! Updating the forces ...
!write (41,*) fp(1,7),fpabf(1,7)
!write (42,*) fp(2,7),fpabf(2,7)
!write (43,*) fp(3,7),fpabf(3,7)

  

return
end subroutine calfo_ABF_BIN


subroutine calfo_ABF_BIN_OMEGA()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev,potist
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:dcsi,icsi,rfilac,histo,     &
                             mean_force,cumul_force1,nhisto,nhisto1, &
                             mean_force1,histo1,omega_abf,&
                             abf_mode,ene_einstein,ene0,fpeinstein,ha_mix,equit
 implicit none 
 real(double), dimension(3,imm) :: fpabf
 real(double) :: force

 
 fpabf(:,:) = zero
 ! Computing the forces from the ABF bins ...
if (abf_mode==1) then
 force = - DOT_PRODUCT(fp(:,7),rfilac(:))
 if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then 
  cumul_force1(icsi) =  cumul_force1(icsi) + force
  mean_force1(icsi) = cumul_force1(icsi)/(1.d0/omega_abf+histo1(icsi))
  fpabf(1:3,7)=rfilac(1:3)*mean_force1(icsi)
 end if
  fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)
end if

if (abf_mode==2) then
  !the mean force is - force = - [ - d U(zeta,q)/d zeta ] 
  force = potist - ene_einstein - ene0 !  d U(zeta,q)/d zeta  
  !  - d_q U(csi,q) 
  fpabf(:,:) = (1.d0-dcsi)*fpeinstein(:,:) + dcsi*fp(:,:)
  fp(:,:)=fpabf(:,:)
 if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then 
  cumul_force1(icsi) =  cumul_force1(icsi) + force
  mean_force1(icsi) = cumul_force1(icsi)/(1.d0/omega_abf+histo1(icsi))
 end if
end if 


if (abf_mode==22) then              
  force = potist + ha_mix*ene_einstein - ene0 - equit !  d U(zeta,q)/d zeta  
  fpabf(:,:) = dcsi*ha_mix*fpeinstein(:,:) + dcsi*fp(:,:) ! -d U(zeta,q)/d q
  fp(:,:)=fpabf(:,:)
  !
  if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then 
   cumul_force1(icsi) =  cumul_force1(icsi) + force
   mean_force1 (icsi) = cumul_force1(icsi)/histo1(icsi)
  end if
  !
end if 


 
! Updating the forces ...
!write (41,*) fp(1,7),fpabf(1,7)
!write (42,*) fp(2,7),fpabf(2,7)
!write (43,*) fp(3,7),fpabf(3,7)

  
return
end subroutine calfo_ABF_BIN_OMEGA


subroutine calfo_ABF_BIN_const_biais()! This subroutine fixes the biais of ABF bin, then launch the dynamics and save histogram files
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:dcsi,icsi,rfilac,histo,     &
                             mean_force,cumul_force1,nhisto,nhisto1, &
                             mean_force1,histo1
 implicit none

 real(double), dimension(3,imm) :: fpabf
 real(double) :: force

 fpabf(:,:) = zero

! Computing the forces from the ABF bins ...
force = - DOT_PRODUCT(fp(:,7),rfilac(:))

if ((icsi >= -nhisto1).and.(icsi <= nhisto+nhisto1)) then 
fpabf(1:3,7)=rfilac(1:3)*mean_force1(icsi)
end if

! Updating the forces ...
!write (41,*) fp(1,7),fpabf(1,7)
!write (42,*) fp(2,7),fpabf(2,7)
!write (43,*) fp(3,7),fpabf(3,7)

  fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)
  

return
end subroutine calfo_ABF_BIN_const_biais



subroutine calfo_ABF_GAUSSIEN()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY:dcsi,icsi,rfilac,histo,     &
                             mean_force,cumul_force1,nhisto,nhisto1, &
                             mean_force1,histo1,omega_abf,&
                             xi_min,xi_max,ecart_eta,eta_mab,&
                             cumul_force_denom1,x_mol,sigma_carre
 implicit none
 real(double), dimension(3,imm) :: fpabf
 real(double) :: force
 integer:: indice_gaussian

 
 fpabf(:,:) = zero
 ! Computing the forces from the ABF bins ...
 force = - DOT_PRODUCT(fp(:,7),rfilac(:))
do indice_gaussian=icsi-ecart_eta,icsi+ecart_eta
 if ((indice_gaussian >=-nhisto1).and.(indice_gaussian <= nhisto+nhisto1)) then 
  cumul_force1(indice_gaussian)=cumul_force1(indice_gaussian) +force*dexp(-(dcsi-x_mol(indice_gaussian))**2/2.d0*sigma_carre)

  cumul_force_denom1(indice_gaussian)=cumul_force_denom1(indice_gaussian)+dexp(-(dcsi-x_mol(indice_gaussian))**2/2.d0*sigma_carre)
  mean_force1 (indice_gaussian) = cumul_force1(indice_gaussian)/cumul_force_denom1(indice_gaussian)

  if(indice_gaussian == icsi) then
  fpabf(1:3,7)=rfilac(1:3)*mean_force1(icsi)
  endif
 end if
enddo
! Updating the forces ...
!write (41,*) fp(1,7),fpabf(1,7)
!write (42,*) fp(2,7),fpabf(2,7)
!write (43,*) fp(3,7),fpabf(3,7)

  fp(1:3,:)=fp(1:3,:)+fpabf(1:3,:)
  
return
end subroutine calfo_ABF_GAUSSIEN


