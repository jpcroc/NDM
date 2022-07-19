
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

#ifdef LAMMPS_VERSION
interface

  subroutine calcforce_lammps(nat,posa,boxl,forca, posa_out, energy)
   integer,      intent(in)                            :: nat
   real(kind=8), intent(in),  dimension(3*nat)         :: posa
   real(kind=8), dimension(3), intent(inout)           :: boxl
   real(kind=8), intent(out), dimension(3*nat)         :: forca,posa_out
   real(kind=8), intent(out)                           :: energy
  end subroutine

end interface

#endif


#ifdef LAMMPS_VERSION

  do i=1,nat
    ttemp(:)=matmul(passage, xp(:,i))
    posa(i)=ttemp(1)/A2cm
    posa(nat+i)=ttemp(2)/A2cm
    posa(2*nat+i)=ttemp(3)/A2cm
  end do


  boxl(1)=at(1,1)/A2cm
  boxl(2)=at(2,2)/A2cm
  boxl(3)=at(3,3)/A2cm

  call calcforce_lammps(nat,posa,boxl,forca,posa_out, energy)

  it_phondy=it_phondy+1

  do i=1,nat
    fp(1,i)=forca(i) !/ (A2cm*erg2ev)
    fp(2,i)=forca(nat+i) !/ (A2cm*erg2ev)
    fp(3,i)=forca(2*nat+i) !/ (A2cm*erg2ev)
    ttemp(:)=matmul(passage_inv, fp(:,i))
    fp(:,i) = ttemp(:) / (A2cm*erg2ev)
  end do

  potist=energy/erg2ev

#else

  call calfo_phondy_ndm (ia_s, ia_f)


#endif

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
