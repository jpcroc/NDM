!  
subroutine read_mab_file()
 USE T_kind_param_m, ONLY:  double
 use gen_com_m, ONLY: lenfnam,fnam,angst
 !use tab_imm_m
 USE mab_in_ndm_module, ONLY: dtlang,nlangevin,temperature,KtoERG,a0bcc,deltasph,  &
                              radiussph,nhisto,deltar1,deltar2,abf_type,block,     &
                              sim_mode,rtestlac,langevin_type,gamma,omega_abf,     &
                              omega_einstein,abf_mode,            & 
                              nwrite_histo,sigma_eta,ecart_eta, &
                              eta_mab,eta_ABFee,histo_equi,n_equilibre,           & 
                              maxforce,compute_mode,error_step,nom_deconvo

 namelist /input_mab/ dtlang,nlangevin,temperature,a0bcc,deltasph,radiussph,       &
                      nhisto,deltar1,deltar2,block,abf_type,sim_mode,rtestlac,     &
                      langevin_type,gamma,omega_abf,omega_einstein, nwrite_histo,  &
                      eta_mab,eta_ABFee,histo_equi,n_equilibre,            &
                      maxforce,compute_mode,abf_mode, error_step,nom_deconvo

 character(len=128) :: fnamtin
 integer :: lumab


 rtestlac=0.1d0
 nwrite_histo=1000000
 omega_abf=1.d0
 maxforce=2.d0   ! in order to enhance the max force on the protective domains 
 nom_deconvo=10
 abf_mode = 1
 abf_type = 1
 omega_einstein=5.d0 ! einstein frequecy in THz

 fnamtin = fnam(1:lenfnam)//'.mab'
 write(*,*) 'file name', fnamtin
 lumab = 778
open(unit=lumab, file=fnamtin, status='unknown')
read (lumab, nml=input_mab)

if (block) write(6,*) 'WARNING: Some spheres are in protective domains!'
      if (langevin_type==2) then
        if (abf_mode==2) then
           write(6,*) 'Alchemical transition not yet implemented with the underdamped Langevin'
           write(6,*) 'put langevin_type = 1 and restart'
           write(6,*) 'stop in <read_mab_file>'
           stop 
        end if 
        select case (abf_type)
           case (1)
                  write(6,*) ' Dumped Langevin dynamics'
           case (2)  
                  if (abf_mode==1) write(6,*) ' Dumped Langevin + ABF BIN dynamics'
                  if (abf_mode==2) write(6,*) ' Dumped Langevin + ABF BIN dynamics + External parameter'
           case (3) 
                  write(6,*) ' Dumped Langevin + ABF BIN dynamics with Omega'
           case (4) 
                  write(6,*) ' Dumped Langevin + ABF GAUSSIAN dynamics'
           case (5) 
                  write(6,*) ' Dumped Langevin + ABF EE dynamics'
           case (6) 
                  write(6,*) ' Dumped Langevin + ABF EE dynamics + constant biais'
           case (7) 
                  write(6,*) ' Dumped Langevin + ABF BIN dynamics + constant biais'
           case (8) 
                  write(6,*) ' Dumped Langevin + ABF EE dynamics + iterative'
        end select 
      end if

      if (langevin_type==1) then
        select case (abf_type)
           case (1)
                  write(6,*) ' Overdumped Langevin dynamics'
           case (2)
                  if (abf_mode==1) write(6,*) ' Dumped Langevin + ABF BIN dynamics'
                  if (abf_mode==2) write(6,*) ' Dumped Langevin + ABF BIN dynamics + External parameter'
           case (3) 
                  write(6,*) ' Overdumped Langevin + ABF BIN dynamics with Omega'
           case (4) 
                  write(6,*) ' Overdumped Langevin + ABF GAUSSIAN dynamics'
           case (5) 
                  write(6,*) ' OverDumped Langevin + ABF EE dynamics'
           case (6) 
                  write(6,*) ' OverDumped Langevin + ABF EE dynamics + constant biais'
           case (7) 
                  write(6,*) ' OverDumped Langevin + ABF BIN dynamics + constant biais'
           case (8) 
                  write(6,*) ' OverDumped Langevin + ABF EE dynamics + iterative'
        end select 
      end if

     if ((abf_mode==2) .and. ((abf_type==1).or.(abf_type==3).or.(abf_type==4).or. &
                             (abf_type==5).or.(abf_type==6).or.(abf_type==7).or. &
                             (abf_type==9) ) ) then
         write(6,*) 'The is no ABF implementation for this mode'
         write(6,*) 'abf_type .....',abf_type
         write(6,*) 'abf_mode .....',abf_mode
         write(6,*) '<stop in read_mab_file>'
         stop
    end if 

    if (abf_mode==2) then
     write(*,'("Einstein frequency (omega_einstein)..........:",D15.4)') omega_einstein
        !instead that I will a file with all the einstein  frequencies 
    end if 

     select case (sim_mode)
      case (1) 
           write(6,*) 'The simulatuion check the first passage time and '
           write(6,*) 'will stop once the vacacy reach the final postion'
      case (2) 
       write(6,'("The simulation stops after nlangevin steps....:",i9)')  nlangevin
     end select 




write(*,'("a0 of the cubic unit cell....................:",D15.4)') a0bcc 
write(*,'("Langevin time step in s......................:",D15.4)') dtlang 
write(*,'("Total number of steps .......................:",I9)')  nlangevin
write(*,'("Langevin temperature in K....................:",F8.1)') temperature
write(*,'("Langevin dumping coefficient ................:",D15.4)') gamma


if (abf_type==3) then
 write(*,'("Omega ABF BIN............... ................:",D15.4)') omega_abf
end if

if (abf_type==4) then
 write(*,'("eta_mab the width of the Gaussian in bins.....:",D15.4)') eta_mab 
end if

if (block) then
  write(*,'("Radius of the blocking spheres (1nn units) ..:",D15.4)') radiussph
  write(*,'("Width of the FD function in A................:",D15.4)') deltasph
end if

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
