subroutine langevin_overdamped_csi()
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
    USE mab_in_ndm_module, only: dtlang,abf_mode,block,gamma,dcsi, &
                                 ene_einstein,ene0, temperature,mean_force1,icsi, &
                                 limit1,limit2,langevin_type,m_i,it_mab,delta_z,nhisto,nhisto1, &
                                 lang_factor,it_mab,it_stop,limit1m, limit1p


    implicit none
    real(double)             :: noise,tmp_dcsi,tmp_force
    real(double) :: dcsi_ini, force_csi
    real(double), save :: average_temp=0.d0
    integer, save :: it_temp = 0
    integer :: it_count 
     

     !ddd ffactor=5.d+8
     if ((dcsi>=limit1m).and.(dcsi<=limit1p)) then
       tmp_force=mean_force1(icsi)
      else 
       if (dcsi<limit1m) tmp_force=mean_force1(-nhisto1)
       if (dcsi>limit1p) tmp_force=mean_force1(nhisto+nhisto1)
       !tmp_force=0.d0
     end if
     dcsi_ini=dcsi
     tmp_dcsi=0.d0
     noise=0
     it_count=0
10 continue
      call genere_bruit_one_value(noise)
      call csi_potential (dcsi,force_csi)

   !  if (langevin_type==2) then
   !   tmp_dcsi = -(potist-ene_einstein -ene0 - tmp_force )*dtlang/(gamma*umass*10.0) + noise*sqrt(2.d0*temperature*dtlang/(gamma*umass*10.0))      
    ! end if 
    ! if (langevin_type==1) then  !overdamped
     !ddd tmp_dcsi = -(potist-ene_einstein -ene0 - tmp_force )*dtlang*angst*ffactor/(gamma*m_i(1,1)) + noise*sqrt(2.d0*temperature*dtlang*angst*ffactor/(gamma*m_i(1,1)))      
     tmp_dcsi = -(potist-ene_einstein -ene0 - tmp_force )*lang_factor*dtlang*angst*angst/(gamma*m_i(1,1)) + noise*sqrt(2.d0*lang_factor*temperature*dtlang*angst*angst/(gamma*m_i(1,1))) -force_csi* lang_factor*dtlang*angst*angst/(gamma*m_i(1,1))     
   !  end if

! the best choise is factor = 5*d+8. actually for this value the 
!dtlang * angst * ffactor / gamma (for ffactor=5.d+8 gamma=1.d+14 and dtlag=2*1d-15) is nothing else than 25* dtlan*dtlang * angst * angst  !!!!! 
 
     dcsi=tmp_dcsi  + dcsi_ini
    if (it_count==0)    write(747,'("  ",i6, 2E16.2, i8, 4D23.7)') it_mab, dcsi,tmp_dcsi, nint(dcsi/delta_z), potist-ene0, ene_einstein, force_csi, tmp_force
    if (it_count/=0)    write(767,'("  ",i6, 2E16.2, i8, 5D23.7)') it_mab, dcsi,tmp_dcsi, nint(dcsi/delta_z), potist-ene0, ene_einstein, force_csi, tmp_force,noise

       average_temp=average_temp+ dabs(tmp_dcsi)

       if ((dcsi<limit1).or.(dcsi>limit2)) then
            it_temp=it_temp+1
            it_count = it_count+1
            if (it_temp==100000) then
              write(*,*) 'WLANGEVIN <lagevin_overdamped_csi>: Too much rejection, the program will stop'
              it_stop=1
              dcsi=dcsi_ini
              go to 11
            end if 
            go to 10
       end if 
11 continue
  
    if (it_count/=0)  write(757,'(2i7,2D15.7)') it_mab, it_temp, tmp_dcsi, average_temp/dble(it_mab + it_temp)
    return
 end subroutine langevin_overdamped_csi



subroutine csi_potential(x,force_csi)

USE T_kind_param_m, ONLY:  double
USE mab_in_ndm_module, ONLY: limit1m,limit1p,limit2m,limit2p,alpha_csi,mode_csi_potential

implicit none
real(double), intent (out) :: force_csi
real(double),intent(in) :: x

if (mode_csi_potential==0) then
    force_csi=0.d0
end if 

if (mode_csi_potential==1) then  
 if (x<=limit2m) then
     force_csi=alpha_csi*(limit2m-limit1m)
  else if ((x<limit1m).and.(x>limit2m)) then
    force_csi=alpha_csi*(x-limit1m)
  else if ((x>=limit1m).and.(x<=limit1p)) then
    force_csi=0.d0 
  else if ((x>limit1p).and.(x<limit2p)) then
    force_csi=alpha_csi*(x-limit1p)
  else if (x>=limit2p) then
    force_csi=alpha_csi*(limit2p-limit1p)
 end if 
end if 
!force_csi=-force_csi
!write(*,'("csi potential", 5d11.2,2D23.7)')  limit2m,limit1m,limit1p,limit2p, x, force_csi, alpha_csi
return
end subroutine csi_potential
