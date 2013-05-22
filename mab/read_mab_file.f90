!  
subroutine read_mab_file()
 USE T_kind_param_m, ONLY:  double
 use gen_com_m, ONLY: lenfnam,fnam,angst
 !use tab_imm_m
 USE mab_in_ndm_module, ONLY: dtlang,nlangevin,temperature,KtoERG,a0bcc,deltasph,  &
                              radiussph,nhisto,deltar1,deltar2,abf_type,block,     &
                              sim_mode,rtestlac,langevin_type,gamma,omega_abf,     & 
                              nwrite_histo,Fermi_percent,a_Fermi,sigma_eta,ecart_eta, &
                              eta_mab

 namelist /input_mab/ dtlang,nlangevin,temperature,a0bcc,deltasph,radiussph,       &
                      nhisto,deltar1,deltar2,block,abf_type,sim_mode,rtestlac,     &
                      langevin_type,gamma,omega_abf,nwrite_histo,Fermi_percent,    &
                      a_Fermi,eta_mab

 character(len=128) :: fnamtin
 integer :: lumab


 rtestlac=0.1d0
 nwrite_histo=1000000
 omega_abf=1.d0
 
 fnamtin = fnam(1:lenfnam)//'.mab'
 write(*,*) 'file name', fnamtin
 lumab = 778
open(unit=lumab, file=fnamtin, status='unknown')
read (lumab, nml=input_mab)

if (block) write(6,*) 'WARNING: Some spheres are in protective domains!'
      if (langevin_type==2) then
        select case (abf_type)
           case (1)
                  write(6,*) ' Dumped Langevin dynamics'
           case (2)  
                  write(6,*) ' Dumped Langevin + ABF BIN dynamics'
           case (3) 
                  write(6,*) ' Dumped Langevin + ABF BIN dynamics with Omega'
           case (4) 
                  write(6,*) ' Dumped Langevin + ABF GAUSSIAN dynamics'
           case (5) 
                  write(6,*) ' Dumped Langevin + ABF EE dynamics'
        end select 
      end if

      if (langevin_type==1) then
        select case (abf_type)
           case (1)
                  write(6,*) ' Overdumped Langevin dynamics'
           case (2)  
                  write(6,*) ' Overdumped Langevin + ABF BIN dynamics'
           case (3) 
                  write(6,*) ' Overdumped Langevin + ABF BIN dynamics with Omega'
           case (4) 
                  write(6,*) ' Dumped Langevin + ABF GAUSSIAN dynamics'
           case (5) 
                  write(6,*) ' Dumped Langevin + ABF EE dynamics'
        end select 
      end if

     select case (sim_mode)
      case (1) 
           write(6,*) 'The simulatuion check the first passage time and '
           write(6,*) 'will stop once the vacacy reach the final postion'
      case (2) 
           write(6,*) ' The simulation stops after nlangevin steps =',nlangevin
     end select 




write(*,'("a0 of the cubic unit cell....................:",D15.4)') a0bcc 
write(*,'("Langevin time step in s......................:",D15.4)') dtlang 
write(*,'("Total number of steps .......................:",I9)')  nlangevin
write(*,'("Langevin temperature in K....................:",F8.1)') temperature
write(*,'("Langevin dumping coefficient ................:",D15.4)') gamma


if (abf_type==3) then
 write(*,'("Omega ABF BIN............... ................:",D15.4)') omega_abf
end if


write(*,'("Radius of the blocking spheres (1nn units) ..:",D15.4)') radiussph
write(*,'("Width of the FD function in A................:",D15.4)') deltasph
write(*,'("Number of the bins of histo..................:",i7)') nhisto
write(*,'("The frequency of writing histo...............:",i7)') nwrite_histo
write(*,'("The first shell of the histo (1nn units).....:",D15.4)') deltar1
write(*,'("The second shell of the histo (1nn units)....:",D15.4)') deltar2
write(*,'("The cutoff radius for ending sim (1nn unit)..:",D15.4)') rtestlac

 temperature=temperature*KtoERG
 deltasph=deltasph/angst
 radiussph=dsqrt(3.d0)*a0bcc*radiussph/(angst*2.d0)

    if (deltar2 < deltar1) then
     write(6,*) 'Increase deltar2. deltar2 should be greater than deltar1'
     write(6,*) 'deltar2,deltar1',deltar2,deltar1
     stop
    end if 

 deltar1=dsqrt(3.d0)*a0bcc*deltar1/(angst*2.d0)
 deltar2=dsqrt(3.d0)*a0bcc*deltar2/(angst*2.d0)
 rtestlac=dsqrt(3.d0)*a0bcc*rtestlac/(angst*2.d0)
 


 



close (lumab)

end subroutine read_mab_file
