!  
subroutine read_mab_file()
 USE T_kind_param_m, ONLY:  double
 use gen_com_m, ONLY: lenfnam,fnam,angst,ev2erg
 !use tab_imm_m
 USE mab_in_ndm_module, ONLY: dtlang,nlangevin,temperature,KtoERG,a0bcc,deltasph,  &
                              radiussph,nhisto,deltar1,deltar2,abf_type,block,     &
                              sim_mode,rtestlac,langevin_type,gamma,omega_abf,     &
                              omega_einstein,abf_mode,            & 
                              nwrite_histo,sigma_eta,ecart_eta, &
                              eta_mab,eta_ABFee,histo_equi,n_equilibre,           & 
                              maxforce,compute_mode,error_step,nom_deconvo,lang_factor, &
                              mode_zeta_potential, alpha_zeta,ntestvacancyjump,ha_mix,  &
                              temperature_zeta_min,temperature_zeta_max

 namelist /input_mab/ dtlang,nlangevin,temperature,a0bcc,deltasph,radiussph,       &
                      nhisto,deltar1,deltar2,block,abf_type,sim_mode,rtestlac,     &
                      langevin_type,gamma,omega_abf,omega_einstein, nwrite_histo,  &
                      eta_mab,eta_ABFee,histo_equi,n_equilibre,            &
                      maxforce,compute_mode,abf_mode, error_step,nom_deconvo,lang_factor, &
                      mode_zeta_potential, alpha_zeta,ntestvacancyjump,ha_mix,temperature_zeta_min,temperature_zeta_max

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
 lang_factor=1.d0
 mode_zeta_potential=0
 alpha_zeta=2.d0
 ntestvacancyjump=200
 ha_mix=0.d0
 temperature_zeta_min=150.d0
 temperature_zeta_max=1800.d0

 fnamtin = fnam(1:lenfnam)//'.mab'
 write(*,*) 'file name', fnamtin
 lumab = 778
open(unit=lumab, file=fnamtin, status='unknown')
read (lumab, nml=input_mab)

if (block) write(6,*) 'WARNING: Some spheres are in protective domains!'
      if (langevin_type==2) then
        !if (abf_mode==2) then
        !   write(6,*) 'Alchemical transition not yet implemented with the underdamped Langevin'
        !   write(6,*) 'put langevin_type = 1 and restart'
        !   write(6,*) 'stop in <read_mab_file>'
        !   stop 
        !end if 
        select case (abf_type)
           case (1)
                  write(6,*) ' Dumped Langevin dynamics'
                  if (abf_mode==2) write(6,*) ' Dumped Langevin + ABF BIN dynamics + External parameter'
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
                  if (abf_mode==1) write(6,*) ' Overdumped Langevin + ABF BIN dynamics with Omega'
                  if (abf_mode==2) write(6,*) ' Overdumped Langevin + ABF BIN dynamics with Omega + External parameter'
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

     if  ((abf_mode==2) .and. ((abf_type==1).or.(abf_type==4).or. &
                             (abf_type==6).or.(abf_type==7).or.  &
                             (abf_type==9))) then
         write(6,*) 'The is no ABF implementation for this mode'
         write(6,*) 'abf_type .....',abf_type
         write(6,*) 'abf_mode .....',abf_mode
         if (abf_type==1) write(6,*) 'MESSAGE: You cannot use Langevin dynamics having abf_mode in the input file.' 
         if (abf_type==1) write(6,*) 'MESSAGE: Put abf_mode=1 and restart the calcultations' 
         write(6,*) '<stop in read_mab_file>'
         stop
    end if 


     if ((abf_mode==22) .and. ((abf_type==1).or.(abf_type==4).or. &
                             (abf_type==6).or.(abf_type==7).or.  &
                             (abf_type==9))) then
         write(6,*) 'The is no ABF implementation for this mode'
         write(6,*) 'abf_type .....',abf_type
         write(6,*) 'abf_mode .....',abf_mode
         if (abf_type==1) write(6,*) 'MESSAGE: You cannot use Langevin dynamics having abf_mode in the input file.' 
         if (abf_type==1) write(6,*) 'MESSAGE: Put abf_mode=1 and restart the calcultations' 
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
write(*,'("Factor for Langevin coefficient ................:",D15.4)') lang_factor


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
if ((abf_mode==2).or.(abf_mode==22)) then
 if (mode_zeta_potential==1) then
  write(*,'("===============THERE IS AN EXTRA POTENTIAL FOR ZETA===========")') 
  write(*,'("The prefactor of zeta potential...............:",D15.4)') alpha_zeta
  alpha_zeta=alpha_zeta*ev2erg
 end if 
end if

if (abf_mode==22) then
write(*,'("ha_mix, U(\zeta,q)=\zeta*[U(q)+ha_mix*U_HA(q))]....:",D15.4)') ha_mix
write(*,'("Temperature min \zeta..............................:",D15.4)') temperature_zeta_min
write(*,'("Temperature max \zeta..............................:",D15.4)') temperature_zeta_max
end if 
 
write(*,'("The cutoff radius for ending sim (1nn unit)..:",D15.4)') rtestlac

 temperature=temperature*KtoERG
 temperature_zeta_min=temperature_zeta_min*KtoERG
 temperature_zeta_max=temperature_zeta_max*KtoERG
 deltasph=deltasph/angst
 radiussph=dsqrt(3.d0)*a0bcc*radiussph/(angst*2.d0)

    if (deltar2 < deltar1) then
     write(6,*) 'Increase deltar2. deltar2 should be greater than deltar1'
     write(6,*) 'deltar2,deltar1',deltar2,deltar1
     stop
    end if 



close (lumab)

end subroutine read_mab_file
