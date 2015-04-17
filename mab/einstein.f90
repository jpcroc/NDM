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
!NOTworkingNOW... ldeplainit=.true.

call initspeed ()

dmtype=tmp_dmtype

return

end subroutine init_einstein_solid




subroutine calfo_einstein_solid ()
 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: zero,im,imm,low_limit,angst,ev2erg,erg2ev, hbar,pi
 USE tab_imm_m
 USE var_pot, ONLY : cm
 USE mab_in_ndm_module, ONLY:xbar,xbarini,xp0,maxforce,unit_omega_to_erg,ene_einstein, &
                             omega_veinstein,fpeinstein,it_mab,n_equilibre,            &
                             itype_einstein,matfor,m_i,temperature
 implicit none                             
 
 integer :: ic,ia,jc,ja
 real(double) :: ddepla_temp(3),temp
 real(double) :: hval,mtot,m2tot,corr 

    do ic=1,3
       xbar(ic)=sum(xp(ic,1:im))/dble(im) ! barycentre sur les particules
    enddo


 fpeinstein (:,:) = zero
! Computing the forces on the protectives spheres...
 ene_einstein=0.d0
if (itype_einstein==0) then
 do ic=1,im
  ddepla_temp(:)=xp(:,ic)-xp0(:,ic)-xbar(:)+xbarini(:)
  fpeinstein(:,ic)=-omega_veinstein(:,ic)**2*unit_omega_to_erg*cm(ityp(ic))*ddepla_temp(:)
  ene_einstein=ene_einstein-DOT_PRODUCT(fpeinstein(:,ic),ddepla_temp(:))
 end do
 ene_einstein=0.5d0*ene_einstein
end if   ! Einstein model ....

if (itype_einstein==1) then 
 ene_einstein=0
 do ic=1,im
   do  ia=1,3
     temp=0.d0
     do jc=1,im
       do ja=1,3
        ! KEEP the sum of the symmetric terms to remove the numerical noise ...
        temp = temp + 0.5d0*(matfor(ic,jc,ia,ja) + matfor(jc,ic,ja,ia))*(xp(ja,jc)-xp0(ja,jc)-xbar(ja)+xbarini(ja))
        ene_einstein = ene_einstein + 0.5d0*matfor(ic,jc,ia,ja)              & 
                          *(xp(ja,jc)-xp0(ja,jc)-xbar(ja)+xbarini(ja))       &
                          *(xp(ia,ic)-xp0(ia,ic)-xbar(ia)+xbarini(ia))

       end do
      end do
     fpeinstein(ia,ic)= - temp
   end do 
 end do
end if   ! HA approximation ....  




end subroutine calfo_einstein_solid

subroutine free_and_correction_einstein()
USE T_kind_param_m, ONLY :double
use gen_com_m, ONLY: im, imm, erg2ev,pi,hbar,volu
use mab_in_ndm_module, ONLY : m_i,omega_einstein,     &
                              einstein_free_3N, einstein_correction, &
                              pbc_correction,temperature,omega_veinstein,rests, itype_einstein
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
! Einstein case .... 
 if (itype_einstein==0) then

  do ia=1,im
    do jx=1,3
       einstein_free_3N= einstein_free_3N -temperature*erg2ev*log(temperature/(omega_veinstein(jx,ia)*hbar*2.d0*pi*1.d+12))
   end do 
  end do
   einstein_correction = - 3.d0/2.d0*temperature*erg2ev*log(hbar**2/(temperature**2*m2tot))  
  rests= -temperature*erg2ev*(im-1)*log(volu)

  if (volu==0.d0) then
  write(*,*) 'MALHEUUUUUUR volume NUL! '
  end if 
!Frenkel PBC
!  pbc_correction=temperature*erg2ev*log(dble(im)/volu*sqrt((hval**2/(2.d0*pi*mtot*temperature))**3) )
!Almarza Correction
pbc_correction=temperature*erg2ev*log(dble(im)/volu*sqrt((2.d0*pi*temperature/(m_i(1,1)*omega_einstein**2*1.d+24)**3)))

end if 

! HA case ....
if (itype_einstein == 1) then
  do ia=1,im
    do jx=1,3
     if (dabs(omega_veinstein(jx,ia)).gt.1d-2) then
       einstein_free_3N= einstein_free_3N -temperature*erg2ev*log(temperature/(omega_veinstein(jx,ia)*hbar*2.d0*pi*1.d+12))
     end if
   end do 
  end do
 einstein_correction=0.d0
 pbc_correction= 0.d0
 rests = 0.d0
end if 



end subroutine free_and_correction_einstein

