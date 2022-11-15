
subroutine calfo_phondy (ia_s, ia_f)
!  use gen_com_m, ONLY : erg2ev,  A2cm, at, potist, rangph
  USE T_kind_param_m, ONLY:  double
  use tab_imm_m
  use gen_com_m
  use tab_imm_m
  use contrainte
  use jqmod

  use phondy_in_ndm_module, ONLY: nat, posa,  forca, boxl, it_phondy, passage, passage_inv

  implicit none
  integer :: ia_s, ia_f
  real(8) :: energy
  real(8) :: posa_out(3*nat)
  integer :: i
  real(kind=8)  :: ttemp(3)




  call calfo_phondy_ndm (ia_s, ia_f)



return
end subroutine calfo_phondy

subroutine calfo_phondy_ndm(ia_s,ia_f)
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
  use contrainte
  use jqmod
  USE phondy_in_ndm_module, only: it_phondy
  implicit none
  integer :: ia_s,ia_f
! one force calculation ....

! NDM part ...

 if(lPrtSigat) sigat(:,:,:)=0. ;
  if (lsigtyp) then
     sigtyp=0. ; sigtyptyp=0.
  end if

  fp(:,:) = zero
  jq=0.0
 if (associated(eatom)) eatom(:)=0

   select case (ipotentiel)

      case (10,11)
         if (itab/=0) then
          if (mod(it_phondy,itab)==0) then
           call caltabt
          endif
         endif
         if (ltabvois.and.mod(it_phondy,itetabvois)==0) call caltabi
         !call calfo()
         call calfoeamtabvois_ph(ia_s,ia_f, xp,  vp,  fp, ielat, iwmax, ityp)
      case default
         if (itab/=0) then
          if (mod(it_phondy,itab)==0) then
           call caltabt
          endif
         endif
         if (ltabvois.and.mod(it_phondy,itetabvois)==0) call caltabi
        !ph if (rangph==0) write(*,*) 'partial_hessian.f90 tesf1.........:', rangph, ia_s
        call calfo
        !ph if (rangph==0) write(*,*) 'partial_hessian.f90 tesf2.........:', rangph, ia_s
   end select


return
end subroutine calfo_phondy_ndm
