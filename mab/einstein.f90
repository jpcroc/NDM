subroutine init_einstein_solid ()
USE  T_kind_param_m, ONLY:  double
USE gen_com_m
USE tab_imm_m
USE var_pot
USE  mab_in_ndm_module, ONLY: temperature,omega_einstein,KtoERG,THZtoK
implicit none 

integer :: tmp_dmtype


tinit=temperature/KtoERG
tempdeplainit=tinit
debyetemp=omega_einstein*THZtoK
tmp_dmtype=dmtype
dmtype=1 
ldeplainit=.true.

call initspeed ()

dmtype=tmp_dmtype

return

end subroutine init_einstein_solid




subroutine calfo_einstein_solid ()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev
 USE tab_imm_m
 USE var_pot, ONLY : cm
 USE mab_in_ndm_module, ONLY:xbar,xbarini,xp0,maxforce,unit_omega_to_erg,ene_einstein, &
                             omega_veinstein,fpeinstein
 implicit none                             
 
 integer :: ic

    do ic=1,3
       xbar(ic)=sum(xp(ic,1:im))/dble(im) ! barycentre sur les particules
    enddo


 fpeinstein (:,:) = zero
! Computing the forces on the protevtives spheres...
 ene_einstein=0.d0
 do ic=1,im
  fpeinstein(:,ic)=-omega_veinstein(:,ic)**2*unit_omega_to_erg*cm(ityp(ic))*(xp(:,ic)-xp0(:,ic)-xbar(:)+xbarini(:))
  ene_einstein=ene_einstein-SUM(fpeinstein(:,ic)*(xp(:,ic)-xp0(:,ic)-xbar(:)+xbarini(:)))
 end do
 ene_einstein=0.5d0*ene_einstein


end subroutine calfo_einstein_solid
