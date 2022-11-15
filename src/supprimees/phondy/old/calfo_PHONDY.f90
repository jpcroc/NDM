subroutine calfo_phondy()
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
  USE phondy_in_ndm_module, only: it_phondy
  implicit none
 
! one force calculation ....

! NDM part ...

         if (itab/=0) then
          if (mod(it_phondy,itab)==0) then
           call caltabt
          endif
         endif
         if (ltabvois.and.mod(it_phondy,itetabvois)==0) call caltabi
        call calfo

return
end subroutine calfo_phondy

