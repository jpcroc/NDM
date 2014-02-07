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
                                 deltar1,deltar2


    implicit none
    real(double)             :: noise,tmp_dcsi,tmp_force
    real(double) :: dcsi_ini,limit1,limit2

     limit1=0.d0-deltar1
     limit2=1.d0+deltar1

     if ((dcsi<=limit1).and.(dcsi>=limit2)) then
       tmp_force=mean_force1(icsi)
      else 
       tmp_force=0.d0
     end if
     dcsi_ini=dcsi
     tmp_dcsi=0.d0
     noise=0
10 continue
     call genere_bruit_one_value(noise)
      tmp_dcsi = -(potist-ene_einstein -ene0 - tmp_force )*dtlang/(gamma*umass*10.0) + noise*sqrt(2.d0*temperature*dtlang/(gamma*umass*10.0))      
       dcsi=tmp_dcsi*angst  + dcsi_ini
      if ((dcsi<limit1).or.(dcsi>limit2)) go to 10
     !debug write(*,'("lang  ",2E16.2, 5D23.7)') dcsi,tmp_dcsi, potist, ene_einstein,noise,temperature,dtlang 
   
    return
 end subroutine langevin_overdamped_csi




