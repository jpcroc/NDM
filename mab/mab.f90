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
  do it_mab=1,nlangevin
    call langevin()
    call reaction()
    if (mod(it_mab,40)==0) then 
     !write(*,*) it_mab
     write(36,*) it_mab,dcsi,xbar(1)-xbarini(1),xp(1,7)
     write(35,*) it_mab,(2.d0*Ecinetique)/(KtoERG*3.d0*dble(im))
    end if
    it=it_mab
    call analyse 
    call controle
     select case (sim_mode)
      case (1) 
           call test_vacancy_position
           
           if (test_end) then
             write(6,*) 'First passage time (step)....:',it_mab 
             write(6,*) 'First passage time (ps)....:',it_mab*dtlang*1.d12 
             stop
           end if
      case (2) 
           continue
     end select 
  end do
  !call force_constant(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)


  write(6,*)
  write(6,*)
  write(6,*)'************  FIN  DE MAB ****************'
  write(6,*)
  write(6,*)
 

  return

  
  end subroutine mab
