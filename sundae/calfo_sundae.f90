subroutine calcforce_lanc(N,pos,nforce,ene_out, it_art) !calcul des forces POUR LANCZOS 
! cette subroutine prends N, les positions, et doit rendre position, forces et l'energie de la configuration
      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      !use art_in_ndm_module
      use tab_imm_m

  implicit none
  integer, intent(in):: N,it_art
  real(double), intent(out) :: ene_out
  real(double), dimension(3*N), target, intent(in) :: pos
  real(double), dimension(:), pointer :: x , y , z
  real(double), intent(out)  :: nforce(3*N)
  integer  :: i,ic
  real(double), parameter:: cmTOang=1.d8


    x => pos(1:N)
    y => pos(N+1:2*N)
    z => pos(2*N+1:3*N)
 
do i=1,N
  xp(1,i)=x(i)/cmTOang
  xp(2,i)=y(i)/cmTOang
  xp(3,i)=z(i)/cmTOang
enddo

!write(*,*) 'force:itart', it_art

        if (itab/=0) then
            if (mod(it_art,itab)==0) then
               call caltabt
            endif
         endif
         if (ltabvois.and.mod(it_art,itetabvois)==0) call caltabi 
         call calfo 

ene_out=potist

do i=1,N
  nforce(i)=fp(1,i)*erg2ev/cmTOang
  nforce(i+N)=fp(2,i)*erg2ev/cmTOang
  nforce(i+2*N)=fp(3,i)*erg2ev/cmTOang
enddo

end subroutine calcforce_lanc

subroutine calfo_teledyn(it_counter)
 USE T_kind_param_m, ONLY:  double 
 use gen_com_m
!use art_in_ndm_module
 use tab_imm_m
 integer , intent(in) :: it_counter
 
        !write(*,*) ltabvois,it_counter ,itab,itetabvois
        if (itab/=0) then
            if (mod(it_counter,itab)==0) then
               call caltabt
            endif
         endif
         if (ltabvois.and.mod(it_counter,itetabvois)==0) call caltabi 
         call calfo 

return
end subroutine

