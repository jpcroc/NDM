!  
subroutine read_mab_file()
 USE T_kind_param_m, ONLY:  double
 use gen_com_m, ONLY: lenfnam,fnam,angst
 !use tab_imm_m
 USE mab_in_ndm_module, ONLY: dtlang,nlangevin,temperature,KtoERG,a0bcc,deltasph,&
                              radiussph

 namelist /input_mab/ dtlang,nlangevin,temperature,a0bcc,deltasph,radiussph

 character(len=128) :: fnamtin
 integer :: lumab



 fnamtin = fnam(1:lenfnam)//'.mab'
 write(*,*) 'file name', fnamtin
 lumab = 778
open(unit=lumab, file=fnamtin, status='unknown')
read (lumab, nml=input_mab)
write(*,'("a0 of the cubic unit cell....................:",D15.4)') a0bcc 
write(*,'("Langevin time step in s......................:",D15.4)') dtlang 
write(*,'("Total number of steps .......................:",I9)')  nlangevin
write(*,'("Langevin temperature in K....................:",F8.1)') temperature
write(*,'("Radius of the blocking spheres (1nn units) ..:",D15.4)') radiussph
write(*,'("Width of the FD function in A................:",D15.4)') deltasph

 temperature=temperature*KtoERG
 deltasph=deltasph/angst
 radiussph=dsqrt(3.d0)*a0bcc*radiussph/(angst*2.d0)



close (lumab)

end subroutine read_mab_file
