subroutine mab
!-----------------------------------------------

      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use var_pot
      use tab_imm_m
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

  
! 
! Copyright LL Cao and all NDM band, April- 2013
  
  integer:: i_iter
  real(double)::temp_read
  logical :: dir_e

  write(6,*)
  write(6,*)
  write(6,*)'************ DEBUT DE MAB ****************'
  write(6,*)
  write(6,*)

  write(6,*)'......READING.....' 
  call read_mab_file()

  write(6,*)'......ALLOCATE....' 
  call allocate_mab()
  write(6,*)'.......INIT.......' 
  !debug .... write(*,*) 'temperature ',temperature/KtoERG 
  call init_mab_in_ndm_module()
  write(6,*)'......PREPARE.....' 

 call prepare_langevin()

  write(6,*)'......LANGEVIN.....' 
 

select case(compute_mode)

case(1)!--------one simulation
  write(6,*)'------COMPUTE MODE IS ONE SIMULATION!!!!!!--------------'
 
if (abf_type == 8) then ! This calculate iterally ABFee ( process to calculate \bar A2 given \bar A1

! Check whether the file exists or not.
  inquire( file="meanforce_input", exist=dir_e )
  if ( dir_e ) then
   open(unit=13,file='meanforce_input',action='read')! this require meanforce data of ABFee!!!!!
    do i_iter=-nhisto2,nhisto+nhisto2
     read(13,*),temp_read,mean_force_ABFee(i_iter)
    enddo
  else
    write(*,*), 'Input file does not exsit!! Verify your input files!!'
   stop
end if


endif
 
do it_mab=1,nlangevin

   select case (langevin_type)
    case (1)
      call langevin_overdamped ()
    case (2)
      call langevin()
    end select 
    
   call reaction()
   !debug  if (mod(it_mab,40)==0) then 
   !debug   write(36,*) it_mab,dcsi,xbar(1)-xbarini(1),xp(1,7)
   !debug  write(35,*) it_mab,(2.d0*Ecinetique)/(KtoERG*3.d0*dble(im))
   !debub  end if
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
   open(unit=989,file='histogram1',status='unknown')
   open(unit=990,file='histogram',status='unknown')
   do i_iter=-nhisto1,nhisto+nhisto1
    write(989,*),i_iter, histo1(i_iter)
   enddo
    do i_iter=1,nhisto
    write(990,*),i_iter, histo(i_iter)
   enddo
 
  close(989)
  close(990)

  end if 

end do
  !call force_constant(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)


  call fill_final_histo()! Attention!! pour ABFee il faut l'histogramme pour calculer l'énergie libre

  call Free_energy_ABF()! Calculate energy landscape for ABF
 
  call create_files()! Create files needed




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
  read(1099,*),temp_read, A_theo(i_iter)
enddo

else
  write(*,*), 'Input file does not exsit!! Verify your input files!!'
  stop

end if



xp=xp0


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
    histo_temp(1:nhisto)=0
    histo1(-nhisto1:nhisto+nhisto1)=0
    histo2(-nhisto2:nhisto+nhisto2)=0
    histo_xi(-nhisto1:nhisto+nhisto1)=0
    histo_zeta(-nhisto1:nhisto+nhisto1)=0.d0
    Free_energy(:)=0
    mean_force(:)=0
    mean_force1(:)=0
    mean_force2(:)=0

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
           write(1102,*),it_mab,sum_error_A,temp_read*sum_error_A
        endif

        if (abf_type == 5) then! for ABFee only
          write(1103,*),it_mab,sum_error_A_bar,temp_read*sum_error_A_bar
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
if (abf_type == 7) then

! Check whether the file exists or not.
inquire( file="meanforce", exist=dir_e )
if ( dir_e ) then

open(unit=991,file='meanforce',action='read')

do i_iter=-nhisto1,nhisto+nhisto1

read(991,*),temp_read, mean_force1(i_iter)

enddo

else
  write(*,*), 'Input file does not exsit!! Verify your input files!!'
  stop

end if

inquire( file="Free_energy_mollifiée_ABFee", exist=dir_e )
if ( dir_e ) then

open(unit=996,file='Free_energy_mollifiée_ABFee',action='read')

do i_iter=-nhisto1,nhisto+nhisto1

read(996,*),temp_read, A_ee(i_iter)

enddo


else
  write(*,*), 'Input file does not exsit!! Verify your input files!!'
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
    write(1104,*),i_iter, histo1(i_iter)
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
read(13,*),temp_read,mean_force_ABFee(i_iter)
enddo

else
  write(*,*), 'Input file does not exsit!! Verify your input files!!'
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
    write(1104,*),i_iter, histo1(i_iter)
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
