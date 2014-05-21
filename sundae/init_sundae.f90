! This file contains initialisation's subroutines for sundae:
! - allocate_tele_vac
! - init_tele_vac


subroutine allocate_tele_vac 
	use sundae_module
	use gen_com_m
	implicit none

	allocate ( m_i(3,imm),gau(6,imm),sig_i(3,imm),rga_i(3,imm),fadd(3,imm),       &
		       xbarini(3))
	allocate ( tmass_tele_vac(imm) )
	allocate (ipovois(15), d2vois(15),xpvois(3,15),xpvoisini(3,15),xtransla(3))

end subroutine allocate_tele_vac
  


subroutine init_tele_vac (ityp)
	use sundae_module
	use gen_com_m
	!-----------------------------------------------
	!   M o d u l e s
	!-----------------------------------------------
	implicit none
	!-----------------------------------------------
	!   D u m m y   A r g u m e n t s
	!-----------------------------------------------
	integer       :: ityp(imm)

	!-----------------------------------------------
	integer   :: ic_local,iatom
	!-----------------------------------------------

	it_tele_vac=0
	tstep_tele_vac=tstep / utemps
	convert_tele_vac=9.6485d0/10000.d0
	im1=im+1
	do ic_local=1,im
		iatom = ic_local
		if (iatom==0) iatom=im     
		tmass_tele_vac(ic_local) = tstep_tele_vac**2*convert_tele_vac/(cm(ityp(iatom))/umass)
	enddo

	usdh_tele_vac= 1.d0/(2.d0*tstep_tele_vac)

	sigma=1.d-2 
	write(6,*) ' kappa ',kappa 
	kappas2 = kappa/two

	! initialisation des parametres du langevin
	gamma_sundae=one/(tstep*1d2)
	m_i(1:3,1:im) = cm(1)

	rga_i(:,:) = exp(-gamma_sundae*tstep/two)

	write(6,*) ' rga_i  ', rga_i(:,1:1)

	!   parametres additionnels

	m_tot = sum(m_i(1,1:im))

	fnampout = fnam(1:lenfnam)//'.pout'
	fnamwout = fnam(1:lenfnam)//'.wout'
	fnamfout = fnam(1:lenfnam)//'.fout'
	fnamsout = fnam(1:lenfnam)//'.sout'
	fnamvout = fnam(1:lenfnam)//'.vout'

	write(6,*) ' fnampout',fnampout


	write(6,*) ' xbarini =',xbarini, im


	return

end subroutine init_tele_vac
