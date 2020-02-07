module send_data_mod
        implicit none
        contains
subroutine send_data(xdata,imm)
  USE T_kind_param_m, ONLY:  double
  implicit none
  integer :: imm 
  real(double),dimension(3,imm) :: xdata
 

  

 return
 end subroutine send_data


 end module
