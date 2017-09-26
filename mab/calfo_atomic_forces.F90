
subroutine calfo_atomic_forces (it_langevin)
!  use gen_com_m, ONLY : erg2ev,  A2cm, at, potist, rangph
  USE T_kind_param_m, ONLY:  double
  use tab_imm_m
  use gen_com_m
  use var_pot
  use mab_in_ndm_module, ONLY:nat, posa,  forca, boxl, it_mab

  implicit none
  integer :: it_langevin
  real(8) :: energy
  real(8) :: posa_out(3*nat)
  integer :: i 


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

!  call calfo_phondy_lammps (ia_s,ia_f,nat,box_at, xp, fp,energy)
  !call caltabt
  !call caltabi
  
!  posa(1      :  nat)=xp(1,1:nat)/A2cm
!  posa(nat+1  :2*nat)=xp(2,1:nat)/A2cm
!  posa(2*nat+1:3*nat)=xp(3,1:nat)/A2cm
  do i=1,nat
  posa(i)=xp(1,i)/A2cm
  posa(nat+i)=xp(2,i)/A2cm
  posa(2*nat+i)=xp(3,i)/A2cm
  end do


  boxl(1)=at(1,1)/A2cm
  boxl(2)=at(2,2)/A2cm
  boxl(3)=at(3,3)/A2cm

  call calcforce_lammps(nat,posa,boxl,forca,posa_out, energy) 

  !fp(1,1:nat)=forca(1:nat) / (A2cm*erg2ev)
  !fp(2,1:nat)=forca(nat+1:2*nat) / (A2cm*erg2ev)
  !fp(3,1:nat)=forca(2*nat+1:3*nat) / (A2cm*erg2ev)
  do i=1,nat
  fp(1,i)=forca(i) / (A2cm*erg2ev)
  fp(2,i)=forca(nat+i) / (A2cm*erg2ev)
  fp(3,i)=forca(2*nat+i) / (A2cm*erg2ev)
  end do

  !write(*,*) it_phondy, rangph, dsqrt(SUM((posa(:)-posa_out(:))**2))

  potist=energy/erg2ev 
   !if (rangph==0) then 
   !write(*,'("1p ", 3G15.7)')   posa(2), posa(nat+2), posa(2*nat+2)
   !write(*,'("2p ", 3G15.7)')  xp(1:3,2)/A2cm
   !write(*,'("1f ", 3G15.7)')   forca(2), forca(nat+2), forca(2*nat+2)
   !write(*,'("2f ", 3G15.7)')  fp(1:3,2)
   !write(*,'(i6,2D20.10)')  it_phondy, potist,energy
   !end if 
   !if (it_phondy > 1000) then
   !write(*,*) 'MAX number of iteration', it_phondy
   !stop
   !end if 

#else

  call calfo_atomic_forces_ndm(it_langevin)

#endif  

!debug   it_mab=it_mab+1
 !if (rangph==0) then
 !write(*,'(i6,2D20.10)')  it_phondy, potist,potist*erg2ev
 !write(*,'("2f ", 3G15.7)')  fp(1:3,2)
 !end if 
 !stop

return
end subroutine calfo_atomic_forces
  





subroutine calfo_atomic_forces_ndm(it_force)
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
  use var_pot
  USE mab_in_ndm_module, only: block
  implicit none
  integer, intent(in) :: it_force

! NDM part ...
          if (itab/=0) then
           if (mod(it_force,itab)==0) then
           call caltabt
           endif
          endif
          if (ltabvois.and.mod(it_force,itetabvois)==0)  call caltabi
        call calfo
!deubg         write(*,*) 'potist', potist
        if (block) call calfoblock()



!here we have only the atomic forces fp(:,) 
return
end subroutine calfo_atomic_forces_ndm



