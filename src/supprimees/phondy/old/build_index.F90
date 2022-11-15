! ************************************************
!           Sous-programme force_constant
! ************************************************
subroutine force_constant 
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use tab_imm_m, ONLY : im
  use phondy_in_ndm_module 
#if(PARAPH) 
  use mpi
#endif 

  implicit none
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, j
  real(double)::time
  integer :: count1,count2,count_rate,count_max 
  integer :: i_start,i_final
#if(PARAPH)
integer :: nb_procsph,codeph
integer :: nratio
#endif

i_start=1
i_final=im

#if(PARAPH)
call MPI_INIT (codeph)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procsph, codeph)
call MPI_COMM_RANK(MPI_COMM_WORLD,rangph,codeph)
 
nratio=int(im/nb_procsph)
  
i_start=rangph*nratio+1

if (rang /= (nb_procsph-1) ) i_final=(rangph+1)*nratio
if (rang == (nb_procsph-1) ) i_final= im 
#endif 

  
  call system_clock (count1,count_rate,count_max)

  matfor(:,:)=0.d0

  call partial_hessian (i_start,i_final,HessianOrder)


#if(PARAPH)
call MPI_FINALIZE(codeph) 
#endif

  call system_clock (count2,count_rate,count_max)
  time=real((count2-count1))/real(count_rate)
  if (rangph==0) write(*,"(' PHONDY: MATFOR  was filled in.......:  ',f16.8,' s')") time
  if (rangph==0) write(*,*) 'PHONDY: MATFOR', matfor (1,2) ,matfor(2,1)
  if (rangph==0) write(*,*) 'PHONDY: MATFOR', matfor (1,5) ,matfor(5,1)


  do i=1,3*im
    do j=i,3*im
    matfor(i,j)=0.5d0*(matfor(i,j)+matfor(j,i))
    matfor(j,i)=matfor(i,j)
    end do
  end do
  call system_clock (count1,count_rate,count_max)
  time=real((count1-count2))/real(count_rate)
  if (rangph==0) write(*,"(' PHONDY: MATFOR  was symmetrized in...:  ',f16.8,' s')") time

  if (rangph==0) write(*,*) 'PHONDY: ...the force constants were  filled'



end subroutine force_constant


