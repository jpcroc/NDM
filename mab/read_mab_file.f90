!  
subroutine read_mab_file()
 USE T_kind_param_m, ONLY:  double
 use gen_com_m, ONLY: lenfnam,fnam
 !use tab_imm_m
 USE mab_in_ndm_module, ONLY: dtlang,nlangevin,temperature,KtoERG

 namelist /input_mab/ dtlang,nlangevin,temperature

 character(len=128) :: fnamtin
 integer :: lumab



 fnamtin = fnam(1:lenfnam)//'.mab'
 write(*,*) 'file name', fnamtin
 lumab = 778
open(unit=lumab, file=fnamtin, status='unknown')
read (lumab, nml=input_mab)
write(*,'("The Langevin time step in s....:",D15.4)') dtlang 
write(*,'("The total number of steps .....:",I9)')  nlangevin
write(*,'("The Langevin temperature in K..:",F8.1)') temperature

 temperature=temperature*KtoERG


close (lumab)

end subroutine read_mab_file
