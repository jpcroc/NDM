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
                             omega_veinstein,fpeinstein,it_mab,n_equilibre
 implicit none                             
 
 integer :: ic
 real(double) :: ddepla_temp(3)

    do ic=1,3
       xbar(ic)=sum(xp(ic,1:im))/dble(im) ! barycentre sur les particules
    enddo


 fpeinstein (:,:) = zero
! Computing the forces on the protevtives spheres...
 ene_einstein=0.d0
 do ic=1,im
  ddepla_temp(:)=xp(:,ic)-xp0(:,ic)-xbar(:)+xbarini(:)
  fpeinstein(:,ic)=-omega_veinstein(:,ic)**2*unit_omega_to_erg*cm(ityp(ic))*ddepla_temp(:)
  ene_einstein=ene_einstein-DOT_PRODUCT(fpeinstein(:,ic),ddepla_temp(:))
 end do
 ene_einstein=0.5d0*ene_einstein


end subroutine calfo_einstein_solid

subroutine free_and_correction_einstein()
USE T_kind_param_m, ONLY :double
use gen_com_m, ONLY: im, imm, erg2ev,pi,hbar,volu
use mab_in_ndm_module, ONLY : m_i,omega_einstein,     &
                              einstein_free_3N, einstein_correction, &
                              pbc_correction,temperature,omega_veinstein,rests
implicit none
real(double) :: oerg, hval,mtot,m2tot
integer :: ia,jx

! oerg=omega_einstein*THZtoK*KtoERG
 hval=2.d0*hbar*pi
 oerg=omega_einstein*hbar*1.d+12
 mtot=SUM(m_i(1,1:im))

 m2tot=0.d0
 do ia=1,im
  m2tot=m2tot+m_i(1,ia)/(mtot*(omega_veinstein(1,ia)*2.d0*pi*1.d+12)**2)
 end do

 einstein_free_3N=0.d0
 do ia=1,im
   do jx=1,3
     einstein_free_3N= einstein_free_3N -temperature*erg2ev*log(temperature/(omega_veinstein(jx,ia)*hbar*2.d0*pi*1.d+12))
   end do 
 end do
  ! This is not general formula
  ! We should replace m_i(1,1) with the appropiate factor. 
!  einstein_correction = - temperature*erg2ev*1.5d0*   &
!   log(m_i(1,1)*oerg**2*hval**2*mtot/(4.d0*temperature**2*pi**2*m2tot) ) 
   einstein_correction = - 3.d0/2.d0*temperature*erg2ev*log(hbar**2/(temperature**2*m2tot))  
 
  rests= -temperature*erg2ev*(im-1)*log(volu)

  if (volu==0.d0) then
  write(*,*) 'MALHEUUUUUUR volume NUL! '
  end if 
!Frenkel PBC
!  pbc_correction=temperature*erg2ev*log(dble(im)/volu*sqrt((hval**2/(2.d0*pi*mtot*temperature))**3) )
!Almarza Correction
pbc_correction=temperature*erg2ev*log(dble(im)/volu*sqrt((2.d0*pi*temperature/(m_i(1,1)*omega_einstein**2*1.d+24)**3)))


end subroutine free_and_correction_einstein

