subroutine sundae
!-----------------------------------------------
      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use tab_imm_m
      use sundae_module 

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      implicit none
! This is the main program for SUNDAE nouveau version 2012
! 
! Copyleft M. Athenes & M.-C. Marinica 
!

	real(double) :: xalea
	!character (len=80) :: fnamtin

	write(6,*)'************ DEBUT DE tele_vacancy ****************'
	write(6,*)
	write(6,*)
	write(6,*)
	write(6,*)
	write(6,*) 'im imm',im,imm
	lta=.true.


	fnamtin = fnam(1:lenfnam)//'.tin'
	! variables de dynamique
	write(*,*) 'file name', fnamtin

	itab=10
	text_teledyn=0.0
	!open(unit=lutin, file=fnamtin, status='unknown', err=567)
	! read (lutin, nml=input_teledyn)

	call random_number(xalea)
	write(6,*) ' xalea ', xalea


	call allocate_tele_vac

	call index_premier_voisin()
	call distance_premier_voisin()

	call init_tele_vac  (ityp)

	call distance_premier_voisin () 

	call LyapLanczos_vac 


	stop


end subroutine  sundae
