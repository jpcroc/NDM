! This file contains initialisation's subroutines for sundae:
! - allocate_tele_vac
! - init_tele_vac


subroutine allocate_tele_vac 
	use sundae_module
	use gen_com_m
	implicit none

	allocate ( m_i(3,imm),gau(6,imm),rga_i(3,imm) )
	allocate (ipovois(15), d2vois(15),xpvois(3,15),xpvoisini(3,15),xtransla(3))

end subroutine allocate_tele_vac
  


subroutine init_tele_vac ()
	use sundae_module
	use gen_com_m
	!-----------------------------------------------
	!   M o d u l e s
	!-----------------------------------------------
	implicit none
	!-----------------------------------------------
	!   D u m m y   A r g u m e n t s
	!-----------------------------------------------

	

	! initialisation des parametres du langevin
	gamma_sundae=one/(tstep*1d2)
	m_i(1:3,1:im) = cm(1)

	rga_i(:,:) = exp(-gamma_sundae*tstep/two)

	write(6,*) ' rga_i  ', rga_i(:,1:1)


	return

end subroutine init_tele_vac
