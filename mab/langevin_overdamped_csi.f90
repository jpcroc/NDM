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
                                 deltar1,deltar2,langevin_type,m_i,it_mab,delta_z


    implicit none
    real(double)             :: noise,tmp_dcsi,tmp_force
    real(double) :: dcsi_ini,limit1,limit2, ffactor
    integer :: it_test

     limit1=0.d0-deltar1
     limit2=1.d0+deltar1
     ffactor=1.d+8
     if ((dcsi<=limit1).and.(dcsi>=limit2)) then
       tmp_force=mean_force1(icsi)
      else 
       tmp_force=0.d0
     end if
     dcsi_ini=dcsi
     tmp_dcsi=0.d0
     noise=0
     it_test=0
10 continue
    
 !    it_test=it_test+1
 !    if (it_test==1) then
      call genere_bruit_one_value(noise)
 !    end if 

 !    if (it_test > 1) noise=-noise

     if (langevin_type==2) then
      tmp_dcsi = -(potist-ene_einstein -ene0 - tmp_force )*dtlang/(gamma*umass*10.0) + noise*sqrt(2.d0*temperature*dtlang/(gamma*umass*10.0))      
     end if 
     if (langevin_type==1) then  !overdamped
      !tmp_dcsi = -(potist-ene_einstein -ene0 - tmp_force )*dtlang/(gamma*umass) + noise*sqrt(2.d0*temperature*dtlang/(gamma*umass))      
      tmp_dcsi = -(potist-ene_einstein -ene0 - tmp_force )*dtlang*angst*ffactor/(gamma*m_i(1,1)) + noise*sqrt(2.d0*temperature*dtlang*angst*ffactor/(gamma*m_i(1,1)))      
     end if 
! g*cm2/s2*s/g

! dt/(g*m) gamma=cm-2

![gamma]=s
      !dcsi=tmp_dcsi*angst  + dcsi_ini
      dcsi=tmp_dcsi  + dcsi_ini
      ! write(*,*) dcsi,deltar1
      !write(*,'("lang  ",2E16.2, 5D23.7)') dcsi,tmp_dcsi, potist, ene_einstein,noise,temperature,dtlang 
      write(747,'("  ",i6, 2E16.2, i6, 2D23.7)') it_mab, dcsi,tmp_dcsi, nint(dcsi/delta_z), potist, ene_einstein
      if ((dcsi<limit1).or.(dcsi>limit2)) go to 10
   
    return
 end subroutine langevin_overdamped_csi




