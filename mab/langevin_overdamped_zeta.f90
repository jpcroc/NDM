subroutine langevin_overdamped_zeta()
   !On input ... the parameters 
   !dt -        integration step size (units, internal units ndm ?)
   !temperature (units ? )
   !rga         (?)
   !schema Euler

   !-----------------------------------------------
   !   M o d u l e s
   !-----------------------------------------------
   USE T_kind_param_m, ONLY:  double
   use gen_com_m
   use tab_imm_m
   use var_pot
   USE mab_in_ndm_module, only: dtlang,abf_mode,block,gamma,dcsi, &
                                ene_einstein,ene0, temperature,mean_force1,icsi, &
                                limit1,limit2,langevin_type,m_i,it_mab,delta_z,nhisto,nhisto1, &
                                lang_factor,it_mab,it_stop,limit1m, limit1p,ha_mix,abf_type,equit, &
                                xi_max, xi_min,dtlang_ini


   implicit none
   real(double)             :: noise,tmp_dcsi,tmp_force
   real(double) :: dcsi_ini, force_zeta, tmp_factor,the_moise
   real(double), save :: average_temp=0.d0
   integer, save :: it_temp = 0
   integer :: it_count
   real(8)  :: crit_zeta,dtlang_zeta
    

    !ddd ffactor=5.d+8
    if (it_mab==1) then
     dcsi=(xi_max-xi_min)/2.d0
!     dcsi=1.0
    end if 
    if ((dcsi>=limit1m).and.(dcsi<=limit1p)) then
      tmp_force=mean_force1(icsi)
     else 
      if (dcsi<limit1m) tmp_force=mean_force1(-nhisto1)
      if (dcsi>limit1p) tmp_force=mean_force1(nhisto+nhisto1)
    end if
    dcsi_ini=dcsi
    tmp_dcsi=0.d0
    noise=0
    it_count=0
    it_temp=0

    tmp_factor=lang_factor*dtlang_ini*angst*angst/(gamma*m_i(1,1))
  if (abf_mode==2) then
    crit_zeta = 0.05
    dtlang_zeta=1.0*crit_zeta**2*gamma*m_i(1,1)/(6.d0*temperature)
!     dtlang_zeta=crit_zeta**2*6.0*gamma*m_i(1,1)*temperature/(ev2erg)**2
!    write(*,*) dtlang_zeta/(gamma*m_i(1,1)),  sqrt(2.d0*temperature*dtlang_zeta/(gamma*m_i(1,1)))
     
  end if 

   10 continue
      call genere_bruit_one_value(noise)
      call zeta_potential (dcsi,force_zeta)

!   write(*,'("lang ", i6, 4e18.4)') it_mab,  &
!            (potist-ene0)*dtlang_zeta/(gamma*m_i(1,1)), &
!            (potist-ene0-ene_einstein-tmp_force)*dtlang_zeta/(gamma*m_i(1,1)), &
!                                                         noise*sqrt(2.d0*temperature*dtlang_zeta/(gamma*m_i(1,1))) , &
!            (potist-ene0)*dtlang_zeta/(gamma*m_i(1,1))/sqrt(2.d0*temperature*dtlang_zeta/(gamma*m_i(1,1)))




  if (abf_mode==2) then
     tmp_dcsi = -(potist -ene0 -ene_einstein-tmp_force  )*dtlang_zeta/(gamma*m_i(1,1))      &
                  + noise*sqrt(2.d0*temperature*dtlang_zeta/(gamma*m_i(1,1)))              &
                  + force_zeta*tmp_factor*10.d0  
                  ! *100 for lang with mass under
                  ! * 10 for   over
  end if 


  if (abf_mode==22) then
     crit_zeta=0.01
      dtlang_zeta=crit_zeta**2*(gamma*m_i(1,1))/(6.d0*(temperature/dcsi))
  end if 

  if (abf_mode==22) then
!oldw     tmp_dcsi = -(potist+ha_mix*ene_einstein -ene0 - equit- tmp_force )*tmp_factor  &
!oldw                + noise*sqrt(2.d0*temperature*tmp_factor)  &
!oldw                + force_zeta* tmp_factor*10.d0
     tmp_dcsi = -(potist+ha_mix*ene_einstein -ene0 - equit- tmp_force )*dtlang_zeta/(gamma*m_i(1,1))  &
                 + noise*sqrt(2.d0*temperature*dtlang_zeta/(gamma*m_i(1,1)))  &
                 + force_zeta* tmp_factor*10.d0
   end if 
    the_moise=noise*sqrt(2.d0*lang_factor*temperature*dtlang_zeta*angst*angst/(gamma*m_i(1,1)))
! the best choise is factor = 5*d+8. actually for this value the 
!dtlang_zeta * angst * ffactor / gamma (for ffactor=5.d+8 gamma=1.d+14 and dtlag=2*1d-15) is nothing else than 25* dtlan*dtlang_zeta * angst * angst  !!!!! 
     dcsi=tmp_dcsi  + dcsi_ini
     if ((it_mab< 2000).and.(dabs(tmp_dcsi)>dabs(xi_max-xi_min)/20.d0)) then
      if (dcsi<xi_min) dcsi=xi_min
      if (dcsi>xi_max) dcsi=xi_max
      tmp_dcsi=0
    end if 
 
!debug    if (it_mab==1) write(747,*) " it_mab, dcsi_ini, nint(dcsi_ini/delta_z), tmp_dcsi, &
!debug       nint(tmp_dcsi/delta_z),  potist-ene0, ene_einstein, force_zeta, tmp_force  "


!debug   if (it_count==0)    write(747,'("  ",i6, E16.2, i6, E16.2, i6, 4D23.7)') it_mab, dcsi_ini, nint(dcsi_ini/delta_z), tmp_dcsi, &
!debug        nint(tmp_dcsi/delta_z),  potist-ene0, (potist-ene0)*dtlang_zeta/(gamma*m_i(1,1)) , force_zeta, tmp_force

!debug    if (it_count==0)    write(744,'("  ",i6,  D23.7, i6, D23.7, f12.6)') it_mab, noise, nint(noise*sqrt(2.d0*temperature/(gamma*m_i(1,1)))/delta_z), noise*sqrt(2.d0*temperature/(gamma*m_i(1,1))), delta_z 

  !debug   write(*,*) 'ddd', nint(potist*tmp_factor/delta_z) , nint(ene0*tmp_factor/delta_z),ha_mix*ene_einstein,equit,nint(tmp_force*tmp_factor/delta_z)


    if (it_count/=0)    write(767,'("  ",i6, 2E16.2, i8, 5D23.7)') it_mab, dcsi,tmp_dcsi, &
                        nint(dcsi/delta_z), -(potist-ene_einstein -ene0 - tmp_force )*tmp_factor , &
                        the_moise, force_zeta*tmp_factor, tmp_force

       average_temp=average_temp+ dabs(tmp_dcsi)

       if ((dcsi<limit1).or.(dcsi>limit2)) then
            it_temp=it_temp+1
            it_count = it_count+1
            if (it_temp==10000) then
              write(*,*) 'WLANGEVIN <langevin_overdamped_csi>: Too much rejection, the program will stop'
!              it_stop=1
              dcsi=dcsi_ini
              go to 11
            end if 
            go to 10
       end if 
11 continue
  
!debug    if (it_count/=0)  write(757,'(2i7,2D15.7)') it_mab, it_temp, tmp_dcsi, &
!debug                      average_temp/dble(it_mab + it_temp)
    return
 end subroutine langevin_overdamped_zeta



subroutine zeta_potential(x,force_zeta)

USE T_kind_param_m, ONLY:  double
USE mab_in_ndm_module, ONLY: limit1m,limit1p,limit2m,limit2p,alpha_zeta,mode_zeta_potential

implicit none
real(double), intent (out) :: force_zeta
real(double),intent(in) :: x

if (mode_zeta_potential==0) then
    force_zeta=0.d0
end if 

!it assumes that the potential is 1/2*alpha_zeta*(x-l1m)**2 for x < l1m and 1/2*alpha_zeta*(x-l1p)**2 for x > l1p

if (mode_zeta_potential==1) then  
 if (x<=limit2m) then
     !force_zeta=-alpha_zeta*(limit2m-limit1m)
     force_zeta=-alpha_zeta*(x-limit1m)**3
  else if ((x<limit1m).and.(x>limit2m)) then
    force_zeta=-alpha_zeta*(x-limit1m)
  else if ((x>=limit1m).and.(x<=limit1p)) then
    force_zeta=0.d0 
  else if ((x>limit1p).and.(x<limit2p)) then
    force_zeta=-alpha_zeta*(x-limit1p)
  else if (x>=limit2p) then
    !force_zeta=-alpha_zeta*(limit2p-limit1p)
    force_zeta=-alpha_zeta*(x-limit1p)**3
 end if 
end if 
!write(*,'("csi potential", 5d11.2,2D23.7)')  limit2m,limit1m,limit1p,limit2p, x, force_zeta, alpha_zeta
return
end subroutine zeta_potential
