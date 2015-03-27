subroutine calcforce_lanc(N,pos,nforce,ene_out, it_art) !calcul des forces POUR LANCZOS 
! cette subroutine prends N, les positions, et doit rendre position, forces et l'energie de la configuration
	USE T_kind_param_m, ONLY:  double
	use gen_com_m
	!use art_in_ndm_module
	use tab_imm_m
        use var_pot

	implicit none
	integer, intent(in):: N,it_art
	real(double), intent(out) :: ene_out
	real(double), dimension(3*N),  intent(in) :: pos
	real(double), intent(out)  :: nforce(3*N)
	integer  :: i
	real(double), parameter:: cmTOang=1.d8


	xp(1,1:N) = pos(1:N)
	xp(2,1:N) = pos(N+1:2*N)
	xp(3,1:N) = pos(2*N+1:3*N)

	xp(1:3,1:N)=xp(1:3,1:N)/cmTOang

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
        use var_pot
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

