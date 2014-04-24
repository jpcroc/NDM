subroutine calfo_mab()
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
  USE mab_in_ndm_module, only: it_mab,abf_type,abf_mode, block,histo_equi,n_equilibre,ene0
  implicit none
  integer :: it_langevin
 
  it_langevin=it_mab
! one force calculation ....

! NDM part ...

          if (itab/=0) then
           if (mod(it_langevin,itab)==0) then
           call caltabt
           endif
          endif
          if (ltabvois.and.mod(it_langevin,itetabvois)==0)  call caltabi
        call calfo
        if (block) call calfoblock()

!ABF part ...

      if ((abf_mode==2).or.(abf_mode==22)) call calfo_einstein_solid ()
        call fill_histo()

!---if we only want to fill histogram after n_equilibre steps-------------
        if ((histo_equi == .true.) .And. (it_mab >= n_equilibre)) then
        call fill_histo_equilibre()
        endif
!up to here we have fp(:,:) - the forces on atomic configurations
!--------------------------------------------------------------------------

      
    select case (abf_type)
           case (1)
                  continue
           case (2)  
                  call calfo_ABF_BIN
           case (3) 
                  call calfo_ABF_BIN_OMEGA
          case (4)  
                  call calfo_ABF_Gaussien ! ABF Gaussian
          case (5) 
                  call calfo_ABFee
          case (6)
                  call calfo_ABFee_const_biais
          case (7)
                  call calfo_ABF_BIN_const_biais
          case (8)
                  call calfo_ABFee_iter
   end select 
!here we have fp(:,)+fpabf(:,:), the forces from ABF part. 
return
end subroutine calfo_mab



