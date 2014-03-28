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
                                 deltar1,deltar2,langevin_type,m_i,it_mab,delta_z,nhisto,nhisto1, &
                                 lang_factor,it_mab


    implicit none
    real(double)             :: noise,tmp_dcsi,tmp_force
    real(double) :: dcsi_ini,limit1,limit2, ffactor
    real(double), save :: average_temp=0.d0
    integer, save :: it_temp = 0
    
     

     limit1=0.d0-deltar1
     limit2=1.d0+deltar1
     !ddd ffactor=5.d+8
     if ((dcsi>=limit1).and.(dcsi<=limit2)) then
       tmp_force=mean_force1(icsi)
      else 
       !debug if (dcsi<=limit1) tmp_force=mean_force1(-nhisto1)
       !debug if (dcsi>=limit2) tmp_force=mean_force1(nhisto+nhisto1)
       tmp_force=0.d0
     end if
     dcsi_ini=dcsi
     tmp_dcsi=0.d0
     noise=0
10 continue
      call genere_bruit_one_value(noise)


     if (langevin_type==2) then
      tmp_dcsi = -(potist-ene_einstein -ene0 - tmp_force )*dtlang/(gamma*umass*10.0) + noise*sqrt(2.d0*temperature*dtlang/(gamma*umass*10.0))      
     end if 
     if (langevin_type==1) then  !overdamped
     !ddd tmp_dcsi = -(potist-ene_einstein -ene0 - tmp_force )*dtlang*angst*ffactor/(gamma*m_i(1,1)) + noise*sqrt(2.d0*temperature*dtlang*angst*ffactor/(gamma*m_i(1,1)))      
     tmp_dcsi = -(potist-ene_einstein -ene0 - tmp_force )*lang_factor*dtlang*angst*angst/(gamma*m_i(1,1)) + noise*sqrt(2.d0*lang_factor*temperature*dtlang*angst*angst/(gamma*m_i(1,1)))      
     end if

! the best choise is factor = 5*d+8. actually for this value the 
!dtlang * angst * ffactor / gamma (for ffactor=5.d+8 gamma=1.d+14 and dtlag=2*1d-15) is nothing else than 25* dtlan*dtlang * angst * angst  !!!!! 
 
     dcsi=tmp_dcsi  + dcsi_ini
       !write(747,'("  ",i6, 2E16.2, i6, 3D23.7)') it_mab, dcsi,tmp_dcsi, nint(dcsi/delta_z), potist-ene0, ene_einstein, tmp_force

       average_temp=average_temp+ dabs(tmp_dcsi)
       if ((dcsi<limit1).or.(dcsi>limit2)) then
            it_temp=it_temp+1
            if (it_temp==100000) then
              write(*,*) '<lagevin_overdamped_csi>: Too much rejection, the program will stop'
              stop
            end if 
            go to 10
       end if 
    if (mod(it_mab,100)==0)  write(757,'(2i7,2D15.7)') it_mab, it_temp, tmp_dcsi, average_temp/dble(it_mab + it_temp)
    return
 end subroutine langevin_overdamped_csi




