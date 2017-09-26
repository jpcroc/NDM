subroutine mab
!-----------------------------------------------

      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use var_pot
      use tab_imm_m
#if(LAMMPS_VERSION) 
! use mkl_service
 use mpi
 use mod_mpi_mab 
 use LAMMPS
 use vars_lammps
#endif



!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      !use DEFS
      use mab_in_ndm_module
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      implicit none
 integer :: iloc
  
! 
! Copyleft NDM Dec-2014
  
  integer:: i_iter
  real(double)::temp_read,tmp2,corr3N,corr3Nm3
  logical :: dir_e
  integer ::  stop_sim, abf_mode_r, abf_type_r, it_mab_r, nhisto_r, nhisto1_r, nhisto2_r
  real(double) :: delta_z_r 
  character(len=1) :: ac_read  


  write(6,*)
  write(6,*)
  write(6,'("MAB: ************ DEBUT DE MAB ****************")')
  write(6,*)
#if(LAMMPS_VERSION)
  if (rangph==0) write(6,*)'---------------LAMMPS FORCES-----------------'
#endif 
 write(6,*)
  if (rangmab==0) call  print_mab(rangmab)

!  call random_seed()

  if (rangmab==0) write(6,'("MAB:           ......ALLOCATE....           ")') 
  call allocate_mab()

  if (rangmab==0) write(6,'("MAB:           ......READING.....           ")') 
  call read_mab_file()

  if (rangmab==0) write(6,'("MAB:           .......INIT.......           ")') 

  call read_gin_file()
  call init_mab_in_ndm_module()

#if(LAMMPS_VERSION)
   if (rangph==0) write(6,*) 'init communicators for lammps...'
   call define_communicators()
!v1   call MPI_BARRIER(orig_group,codeph)
   call init_potential_lammps()
!v1   call MPI_BARRIER(orig_group,codeph)
   if (rangph==0) write(6,*) '...init LAMMPS forces'
   call test_minimum()
#endif




  if (itype_reaction==1) then
   call init_neb_reaction()
   call read_neb_from_ndm()
   call interpolate_the_neb_images()
   call test_energy_along_reaction()
!   call build_defect_force_along_reaction()
   stop
  end if 
  if (rangmab==0) write(6,'("MAB:           ......PREPARE.....          ")') 

neq_lang=0
if (abf_mode==abf_mode_reaction) neq_lang=0
if (abf_mode==abf_mode_alchemical) neq_lang=4000

call prepare_langevin()

if (rangmab==0) then
     write(6,'("MAB: --------Few consitencies check--for--langevin steps------")')
     write(6,'("MAB: crit_lanngevin_dist (in ang) ........................:",E25.12E3)') crit_langevin_dist*angst

   if (abf_mode==abf_mode_reaction) then
     write(6,'("MAB: (nhisto2-nhisto1)*delta_z (in ang)...................:",E25.12E3)') dble(nhisto2-nhisto1)*delta_z*angst
     write(6,'("MAB: (nhisto1        )*delta_z (in ang)...................:",E25.12E3)') dble(nhisto1)*delta_z*angst
   end if

   if (crit_langevin_dist / (dble(nhisto2-nhisto1)*delta_z) > 1.d0 ) then 
     write(6,'("MAB:  WARNING !!!! Langevin  vs histo settings  !!! WARNING ")')
     write(6,'("MAB: the Langevins step is too large compared with (nhisto2-nhisto1)*delta_z, possible solutions ( or / and ):")')
     write(6,'("MAB:  - decrease the Langevins step")')
     write(6,'("MAB:  - increase the difference nhisto2 - nhisto1 ")')
   end if 

   if (crit_langevin_dist / (dble(nhisto1)*delta_z) > 1.d0 ) then 
     write(6,'("MAB:  WARNING !!!! Langevin vs nhisto1  !!! WARNING ")') 
     write(6,'("MAB: the Langevins step is too large compared with nhisto1*delta_z, possible solutions ( or / and ):")')
     write(6,'("MAB:  - decrease the Langevins step")')
     write(6,'("MAB:  - increase the difference nhisto2 - nhisto1 ")')
   end if 


   if (crit_langevin_dist / (dble(nhisto2)*delta_z) > 1.d0 ) then 
     write(6,'("MAB: the Langevins step is huge. The system will leave the +/- nhisto2 histograms extensions  ")')
     write(6,'("MAB:  - decrease the Langevins step")')
     write(6,'("MAB:  - increase the difference nhisto2  ")') 
     ! stop '< mab.f90 >'
   end if  
end if 


if (block) then
   write(6,'("MAB: width of the FD block  (in ang).....................:",E25.12E3)') deltasph*angst

   if ((crit_langevin_dist*angst / deltasph ) > 1.d0)  then 
     write(6,'("MAB:  WARNING !!!! Langevin  vs FD width  !!! WARNING ")') 
     write(6,'("MAB: the Langevins step is huge or the width of the FD function for protectives domains are too narrow  ")')
     write(6,'("MAB:  - decrease the Langevins step .....................:", E25.12E3)') crit_langevin_dist
     write(6,'("MAB:  - increase the width of the FD function............:", E25.12E3)')  deltasph
     write(6,'("MAB:  ratio Lang / FD width ", E25.12E3)')  crit_langevin_dist/deltasph
     !stop '< mab.f90 >'
   end if 
   write(6,'("MAB: --------------END  Consistencies ---------------")')
end if !rangmab 



!old if (abf_type==1)  call test_minimum_abf()
if (abf_type==1)  call test_minimum()
if ((abf_mode==abf_mode_alchemical).or.(abf_mode==abf_mode_temperature)) then
 call test_minimum () 
 !old call test_minimum_abf () 
 !call init_einstein_solid  ()
  do it_en=1,neq_lang
   if (langevin_type==1) call langevin_overdamped ()
   if (langevin_type==2) call langevin()
   it=it_en
   call analyse 
   call controle
  end do

 it_en=-1
end if 



  write(6,*)'......LANGEVIN.....' 
select case(compute_mode)

case(1)!--------one simulation
  write(6,*)'------COMPUTE MODE IS ONE SIMULATION!!!!!!--------------'
 
!if (abf_type == 8) then ! This calculate iterally ABFee ( process to calculate \bar A2 given \bar A1
if (abf_restart) then
if (abf_mode==abf_mode_reaction) then
  inquire( file="meanforce_input", exist=dir_e )
  if ( dir_e ) then
   open(unit=133,file='meanforce_input',action='read')! this require meanforce data of ABFee!!!!!
    do i_iter=-nhisto2,nhisto+nhisto2
     read(133,*)temp_read,mean_force_ABF_restart(i_iter)
       if (temp_read-dble(i_iter)/=zero) then
        if (rangmab==0) write(6,'("MAB: Problems in reading meanforce_input file ... stop in mab.f90",i6, 2E25.12E3)') & 
                      i_iter, temp_read, dble(i_iter)
        stop
       end if 
    enddo
    close(133)
  else
    write(*,*)'Input file meanforce_input does not exsit!! Verify your input files!!'
    stop
  end if
end if 


if ((abf_mode==abf_mode_reaction).or.(abf_mode==abf_mode_alchemical)) then
  inquire( file="A_ee_input", exist=dir_e )
  if ( dir_e ) then
    open(unit=133,file='A_ee_input',action='read')
    
    read(133,*) ac_read, abf_mode_r, abf_type_r, it_mab_r, nhisto_r, nhisto1_r, nhisto2_r, delta_z_r   
    stop_sim=1 
    if ( nhisto   /= nhisto_r) stop_sim=0
    if ( nhisto1  /= nhisto1_r) stop_sim=0
    if ( nhisto2  /= nhisto2_r) stop_sim=0
    if ( abf_mode /= abf_mode_r) stop_sim=0
    if ( dabs(delta_z - delta_z_r)*angst>1.d-8 ) stop_sim=0
    if ( abf_type/=abf_type_r) then
      if (rangmab==0) write(6,'("MAB: WARNING !!!! the abf_type from old simulation and the new simulation is not the same ")')
    end if 
    if (stop_sim==0) then
       if (rangmab==0) then
           write(6,'("MAB: the A_ee_restart file is not compatible with one of the following things. Check ...")')
           write(6,'("MAB: Check that                   abf_mode  ,  nhisto  , nhisto1  , nhisto2  , delta_z")')
           write(6,'("MAB: present simulation         ",2i3,3i5,E25.12E3)') abf_mode  ,  abf_type, &
                nhisto  , nhisto1  , nhisto2  , delta_z   
           write(6,'("MAB: old A_ee_restart simulation",2i3,3i5,E25.12E3)') abf_mode_r,  abf_type_r, &
                nhisto_r, nhisto1_r, nhisto2_r, delta_z_r   
       end if 
       stop '< mab.f90 >'
    end if 

    do i_iter=-nhisto2,nhisto+nhisto2
      read(133,*) temp_read, A_ee_restart(i_iter), P_ee_num_restart(i_iter), P_ee_denom_restart(i_iter)
       if (dabs((temp_read-dble(i_iter)*delta_z)) >  1.d-10 ) then
        if (rangmab==0) write(6,'("MAB: Problems in reading meanforce_input file ... stop in mab.f90")') &
         i_iter, temp_read, dble(i_iter)
        stop
       end if 


    enddo
    close(133)
    A_ee_restart(:)=A_ee_restart(:)/erg2ev
  else 
    write(*,*)'Input file A_ee_input does not exsit!! Verify your input files!!'
    stop
  end if 
end if 

end if 


! Check whether the file exists or not.
!endif
 
do it_mab=neq_lang+1,neq_lang+nlangevin

 !write(*,'(i6,5d20.8)') it_mab, xp(1:3,1)

   select case (langevin_type)
    case (1)
      call langevin_overdamped ()
    case (2)
      call langevin()
    end select 
   if ((abf_mode==abf_mode_alchemical).or.(abf_mode==abf_mode_temperature)) then
         if (.NOT.(abf_type==5)) call langevin_overdamped_zeta()
   end if

   if (.NOT.(abf_type==1))  call reaction()
   ! if (mod(it_mab,40)==0) then 
    write(36,'(i6,3d15.7)') it_mab,dcsi,xbar(1)-xbarini(1),xp(1,7)*1.d+08
   !debug write(35,*) it_mab,it_en,(2.d0*Ecinetique)/(KtoERG*3.d0*dble(im))
   ! end if
    it=it_mab
    if (mod(it_mab,ntestvacancyjump)==0)  call test_displacement ()
    call analyse 
    call controle
     select case (sim_mode)
      case (1) 
           call test_vacancy_position
           
           if (test_end) then
             write(6,*) 'First passage time (step)....:',it_mab 
             write(6,*) 'First passage time (ps)......:',it_mab*dtlang*1.d12 
             stop 
           end if
      case (2) 
           continue
     end select 
  if (abf_type==1) then
    call brute_force_free_energy (itest_stop)
    if (itest_stop==1) then
     write(6,*) 'An unwanted hop happens the simulation will stop'
     write(6,*) ' it_mab   ', it_mab
    end if  
  end if 

! Writing intermediar steps ....
  if (mod(it_mab,nwrite_histo)==0) then
     call on_run_writting_histo_freetemp(it_mab)
  end if 
! End writting intermediare steps ....

 !write(*,*) 'main', it_mab, A_ee(1)
    if (itesauvposition.GT.0) then
        if (mod(it_mab,itesauvposition)==0) then
         call sauveposition (it_mab)
	 if (lsuivinonpbc) then
          write(*,*) 'lsuivinonpbc true not yet implemented with mab'
          stop
	  !call reset_suivinonpbc
	  !call sauvepositionnonpbc (it_mab)
	 end if	 
        end if
     endif



  if (it_stop==1) exit
  if (itest_stop==1) exit

 end do   !end it_mab, langevin
 if (abf_type /= 1) then
  call fill_final_histo()! Attention!! pour ABFee il faut l'histogramme pour calculer l'énergie libre
  call Free_energy_ABF()! Calculate energy landscape for ABF
  call create_files()! Create files needed

 end if 

 if (abf_type==1) then
   if (itest_stop==1)   write(6,*) '----------WLANGEVIN NOT CONVERGED-----------'
      call correct_free_energy_brute(corr3N,corr3Nm3)
      write(*,*) 'outsub', corr3N, corr3Nm3,corr3N*erg2ev
      write(6,*) '----------FREE ENERGY FINAL RESULTS---------' 
      write(6,'("DeltaFreeMD     (eV) ................:  ", F15.7)') Free_energy_brute*erg2ev
      write(6,*) ' '
      write(6,'("DeltaFree3N     (eV) ................:  ", F15.7)') (Free_energy_brute+corr3N)*erg2ev
      write(6,'("FreeTOT3N       (eV) ................:  ", F15.7)') ene0*erg2ev+(Free_energy_brute+corr3N)*erg2ev
      write(6,*) ' '
      write(6,'("DeltaFree(3N-3) (eV) ................:  ", F15.7)') (Free_energy_brute+corr3Nm3)*erg2ev
      write(6,'("FreeTOT(3N-3)   (eV) ................:  ", F15.7)') ene0*erg2ev+(Free_energy_brute+corr3Nm3)*erg2ev
      write(6,'("Average_over N steps ..................:  ", i7)') it_calc_brute
  end if 

 if (abf_mode==2) then
!debug  Write (*,*) 'temp', temperature, temperature/KtoERG, omega_einstein
!debug  Why should be - 
!debug  because it should be: F = F^HA + [ A(1) - A(0) ]
   tmp2= -(Free_energy(0)-Free_energy(nhisto))*erg2ev
   call free_and_correction_einstein()
   if (it_stop==1) write(6,*) '----------WLANGEVIN NOT CONVERGED-----------'
    write(6,*) '----------FREE ENERGY FINAL RESULTS---------' 
    write(6,'("F(Einstein) 3N                   (eV) .......:  ", F15.7)') einstein_free_3N+einstein_correction
    write(6,'("Einstein PBC config  correction  (eV) .......:  ", F15.7)') einstein_correction
    write(6,'("Einstein PBC kinetic correction  (eV) .......:  ", F15.7)') pbc_correction 
    write(6,'("F(Einstein)                      (eV) .......:  ", F15.7)') einstein_free_3N+einstein_correction+pbc_correction
    write(6,'("F(Einstein) - F(Full)            (eV) .......:  ", F15.7)') tmp2
    write(6,'("F(Full3N-3)                      (eV) .......:  ", F15.7)') einstein_free_3N+einstein_correction+tmp2 
    write(6,'("F(Full3N-6)                      (eV) .......:  ", F15.7)') einstein_free_3N+einstein_correction+ &
                                                                           pbc_correction+tmp2 
  endif  !abf_mode==2

  call sauveposition (it_mab)

  write(6,*)
  write(6,*)
  write(6,*)'************  FIN  DE MAB ****************'
  write(6,*)
  write(6,*)
 
case (2)!------error analysis for ABF bin et ABF ee ( abf_type == 2 or abf_type == 5)

write(*,*) '*******************************************************************'
write(*,*) '************Error Analysis!!!!!!!!*********************************'
write(*,*) '*******************************************************************'


call init_random_seed()

write(6,*) '------COMPUTE MODE: Compute statistical error!!!!!!--------------'
!-------It needs entry file: Free_energy_theo.

if ( (abf_type .ne. 2) .and. (abf_type .ne. 5)) then
write(*,*) 'Error !!!! Verify your abf_type!!!! This subroutine calculate only when abf_type == 2 or 5'
stop
endif



inquire( file="Free_energy_theo", exist=dir_e )
if ( dir_e ) then

open(unit=1099,file='Free_energy_theo',action='read')

do i_iter=-nhisto1+1,nhisto+nhisto1
  read(1099,*)temp_read, A_theo(i_iter)
enddo

else
  write(*,*) 'Input file does not exsit!! Verify your input files!!'
  stop

end if



!xp=xp0


  A_ee(:)=0.d0
  error_A(:)=0.d0
  A_dev_ee(:)=0.d0
  P_ee(:)=0.d0
  P_ee_num(:)=0.d0
  P_ee_denom(:)=0.d0    
  A_bar_ee(:)=0.d0
  exp_A_bar(:)=1.d0 
 cumul_force1(:)=0.d0 
    cumul_force_denom1(:)=0.d0
    histo(1:nhisto)=0
    histo1(-nhisto1:nhisto+nhisto1)=0
    histo2(-nhisto2:nhisto+nhisto2)=0
    histo_xi(0:nhisto)=0
    histo_xi1(-nhisto1:nhisto+nhisto1)=0
    histo_zeta(-nhisto1:nhisto+nhisto1)=0.d0
    Free_energy(:)=0
    mean_force(:)=0
    mean_force1(:)=0
    mean_force_ABFee_dyn(:)=0

if (abf_type == 2) then

open(unit=1102,file='error_ABFdirect',status='unknown')

endif

if (abf_type == 5) then
 
open(unit=1103,file='error_Abar_corrected',status='unknown')

endif

  do it_mab=1, nlangevin! BEGIN OF ONE SIMULATION

     select case (langevin_type)
      case (1)
        call langevin_overdamped ()
      case (2)
        call langevin()
    end select 
    
    call reaction()

     it=it_mab
     call analyse 
     call controle
     select case (sim_mode)
     case (1) 
           call test_vacancy_position
           
           if (test_end) then
             write(6,*) 'First passage time (step)....:',it_mab 
             write(6,*) 'First passage time (ps)......:',it_mab*dtlang*1.d12 
             stop 
           end if
    case (2) 
           continue
    end select 


    if(mod(it_mab,error_step)==0) then
        call error_ABF_direct()!calculate directly error using meanforce
        call error_ABFee_proba_corrected()!calculate indirectly error with proba corrected

        temp_read=sqrt(dble(it_mab))

!-----ERROR FILE: 1. steps of simulation, 2. absolute error 3. absolute error*sqrt(steps of simulation
    
        if(abf_type == 2) then
           write(1102,*)it_mab,sum_error_A,temp_read*sum_error_A
        endif

        if (abf_type == 5) then! for ABFee only
          write(1103,*)it_mab,sum_error_A_bar,temp_read*sum_error_A_bar
        endif

    endif
    

  enddo!---END OF ONE SIMULATION


 if (abf_type == 2) then
 close(1102)
 endif

 if (abf_type == 5) then
 close(1103)
 endif

call fill_final_histo()! Attention!! pour ABFee il faut l'histogramme pour calculer l'énergie libre

call Free_energy_ABF()! Calculate energy landscape for ABF
 
call create_files()! Create files needed

  write(6,*)
  write(6,*)
  write(6,*)'************  FIN  DE MAB ****************'
  write(6,*)
  write(6,*)

case (3)! Deconvolution and dynamics with constant biais
!-----this is used when we want to generate other files or launch simulation with data provided/calculated.

write(*,*) '*******************************************************************'
write(*,*) '*********Deconvolution and dynamics with constant biais************'
write(*,*) '*******************************************************************'

if ((abf_type .ne. 7) .and. (abf_type .ne. 6) .and. (abf_type .ne. 5)) then

write (*,*) 'Error !!! This process of Deconvolution or constant biais works only when abf_type = 5, 6, 7 !'
stop

endif

!-------------Dynamics with constant biais for ABFbin
! you need the previous ABFee free energy (A_ee) and the mean force 
if (abf_type == 7) then

! Check whether the file exists or not.
inquire( file="meanforce", exist=dir_e )
if ( dir_e ) then

open(unit=991,file='meanforce',action='read')

do i_iter=-nhisto1,nhisto+nhisto1

read(991,*)temp_read, mean_force1(i_iter)

enddo

else
  write(*,*) 'Input file does not exsit!! Verify your input files!!'
  stop

end if

inquire( file="Free_energy_mollifiee_ABFee", exist=dir_e )
if ( dir_e ) then

open(unit=996,file='Free_energy_mollifiee_ABFee',action='read')

do i_iter=-nhisto1,nhisto+nhisto1

read(996,*) temp_read, A_ee(i_iter)

enddo


else
  write(*,*) 'Input file does not exsit!! Verify your input files!!'
  stop

end if


A_ee(:)=A_ee(:)/erg2eV

do it_mab=1,nlangevin

   select case (langevin_type)
    case (1)
      call langevin_overdamped ()
    case (2)
      call langevin()
    end select 
    
   call reaction()
    it=it_mab
    call analyse 
    call controle
     select case (sim_mode)
      case (1) 
           call test_vacancy_position
           
           if (test_end) then
             write(6,*) 'First passage time (step)....:',it_mab 
             write(6,*) 'First passage time (ps)......:',it_mab*dtlang*1.d12 
             stop 
           end if
      case (2) 
           continue
     end select 

  if (mod(it_mab,nwrite_histo)==0) then
   
  open(unit=1104,file='histogram1_const_biais',status='unknown')
   do i_iter=-nhisto1,nhisto+nhisto1
    write(1104,*)i_iter, histo1(i_iter)
   enddo
  close(1104)

  end if 



end do


  call fill_final_histo()

endif


!---------for ABFee dynamics

if (abf_type == 6 .or. abf_type == 5) then

call deconvolution_ABFee()! Deconvolution process

if (abf_type == 6) then ! This is the process of ABFee with constant biais.



! Check whether the file exists or not.
inquire( file="meanforce_input", exist=dir_e )
if ( dir_e ) then

open(unit=13,file='meanforce_input',action='read')! this require meanforce data of ABFee!!!!!
do i_iter=-nhisto2,nhisto+nhisto2
read(13,*) temp_read,mean_force_ABF_restart(i_iter)
enddo

else
  write(*,*) 'Input file does not exsit!! Verify your input files!!'
  stop

end if


do it_mab=1,nlangevin

   select case (langevin_type)
    case (1)
      call langevin_overdamped ()
    case (2)
      call langevin()
    end select 
    
   call reaction()
    it=it_mab
    call analyse 
    call controle
     select case (sim_mode)
      case (1) 
           call test_vacancy_position
           
           if (test_end) then
             write(6,*) 'First passage time (step)....:',it_mab 
             write(6,*) 'First passage time (ps)......:',it_mab*dtlang*1.d12 
             stop 
           end if
      case (2) 
           continue
     end select 

  if (mod(it_mab,nwrite_histo)==0) then
   
  open(unit=1104,file='histogram1_const_biais',status='unknown')
   do i_iter=-nhisto1,nhisto+nhisto1
    write(1104,*)i_iter, histo1(i_iter)
   enddo
  close(1104)

  end if 



end do


  call fill_final_histo()

endif
 
endif


 

  write(6,*)
  write(6,*)
  write(6,*)'************  FIN  DE MAB ****************'
  write(6,*)
  write(6,*)

end select

  return

  
  end subroutine mab
