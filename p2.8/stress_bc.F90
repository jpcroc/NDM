module stress_bc_mod
        implicit none 
        contains
subroutine stress_bc
	!-----------------------------------------------
 	!   M o d u l e s
 	!-----------------------------------------------
	USE T_kind_param_m, ONLY: double
	use gen_com_m
	use tab_imm_m
	use eam
	use posana
        USE cfg_module
	
	implicit none
	!-----------------------------------------------
 	!   G l o b a l   P a r a m e t e r s
 	!-----------------------------------------------
 	!-----------------------------------------------
 	!   L o c a l   P a r a m e t e r s
 	!-----------------------------------------------
	!-----------------------------------------------
	!   L o c a l   V a r i a b l e s
	!-----------------------------------------------
	integer      :: i, j1, j2
	logical      :: lwrite, lprint
	real(double) :: mo_strain_inf, mo_strain_sup, fact2Ddyn, ef_strainXZ

	call it_countdown

	if (ibound==3) then
		if (it /= 1) call recalc_stress  ! on appelle la routine pour corriger la contrainte appliquée
	end if					 ! (but: avoir une vitesse de déformation constante)

	lwrite = .TRUE.

	fact2Ddyn = 1.
	IF (ldyn2D.EQV..true.) THEN
		fact2Ddyn = 0.
	END IF

	IF (lwrite.EQV..TRUE.) THEN  ! si on veut : genere un fichier de sortie (it, deformation)

	  mo_strain_inf = 0.
	  mo_strain_sup = 0.
	  DO i=1,i_surfMIN
		j1 = b2sINF(i)
		fp(3,j1) = forceatinf + fp(3,j1)
		fp(2,j1) = fact2Ddyn*fp(2,j1)
		j2 = b2sSUP(i)
		fp(3,j2) = forceatsup + fp(3,j2)
	        fp(2,j2) = fact2Ddyn*fp(2,j2)
		mo_strain_inf = mo_strain_inf + xp(3,j1)-ax(3,j1)
		mo_strain_sup = mo_strain_sup + xp(3,j2)-ax(3,j2)
	  END DO
	  IF (i_surfSUP == i_surfMAX .AND. i_surfINF /= i_surfMAX) THEN
		Do i=i_surfMIN+1,i_surfMAX
		   j2 = b2sSUP(i)
		   fp(3,j2) = forceatsup + fp(3,j2)
		   fp(2,j2) = fact2Ddyn*fp(2,j2)
		   mo_strain_sup = mo_strain_sup + xp(3,j2)-ax(3,j2)
		End do
	  ELSE IF (i_surfINF == i_surfMAX .AND. i_surfSUP /= i_surfMAX) THEN
		Do i=i_surfMIN+1,i_surfMAX
		   j1 = b2sINF(i)
		   fp(3,j1) = forceatinf + fp(3,j1)
		   fp(2,j1) = fact2Ddyn*fp(2,j1)
		   mo_strain_inf = mo_strain_inf + xp(3,j1)-ax(3,j1)
		End do
	  END IF
	  mo_strain_inf = mo_strain_inf / i_surfINF
	  mo_strain_sup = mo_strain_sup / i_surfSUP
	  ef_strain     = (mo_strain_sup - mo_strain_inf)/gap

 	  If (it==1) Open (unit=121, file='stress_it.dat', status='unknown',action='write')
	  Write (121,'(i7,2x,e24.16)') it,ef_strain

	ELSE IF (lwrite.EQV..false.) THEN  ! ici, on applique nos CL sans generer de fichier de sortie (itération, déformation)

	  DO i=1,i_surfMIN
		j1 = b2sINF(i)
		fp(3,j1) = forceatinf + fp(3,j1)
		fp(2,j1) = fact2Ddyn*fp(2,j1)
		j2 = b2sSUP(i)
		fp(3,j2) = forceatsup + fp(3,j2)
	        fp(2,j2) = fact2Ddyn*fp(2,j2)
	  END DO
	  IF (i_surfSUP == i_surfMAX .AND. i_surfINF /= i_surfMAX) THEN
		Do i=i_surfMIN+1,i_surfMAX
		   j2 = b2sSUP(i)
		   fp(3,j2) = forceatsup + fp(3,j2)
		   fp(2,j2) = fact2Ddyn*fp(2,j2)
		End do
	  ELSE IF (i_surfINF == i_surfMAX .AND. i_surfSUP /= i_surfMAX) THEN
		Do i=i_surfMIN+1,i_surfMAX
		   j1 = b2sINF(i)
		   fp(3,j1) = forceatinf + fp(3,j1)
		   fp(2,j1) = fact2Ddyn*fp(2,j1)
		End do
	  END IF
	
	END IF

	! ------------------------------------------------------------------------------------------

  ! -------- POUR EDITER UN FICHIER DE SORTIE (strainXZ en fonction de l'it)

!	  mo_strain_inf = 0.
!	  mo_strain_sup = 0.
!	  ef_strainXZ = 0.
!	  DO i=1,i_surfMIN
!		j1 = b2sINF(i)
!		j2 = b2sSUP(i)
!		mo_strain_inf = mo_strain_inf + xp(1,j1)-ax(1,j1)
!		mo_strain_sup = mo_strain_sup + xp(1,j2)-ax(1,j2)
!	  END DO
!	  IF (i_surfSUP == i_surfMAX .AND. i_surfINF /= i_surfMAX) THEN
!		Do i=i_surfMIN+1,i_surfMAX
!		   j2 = b2sSUP(i)
!		   mo_strain_sup = mo_strain_sup + xp(1,j2)-ax(1,j2)
!		End do
!	  ELSE IF (i_surfINF == i_surfMAX .AND. i_surfSUP /= i_surfMAX) THEN
!		Do i=i_surfMIN+1,i_surfMAX
!		   j1 = b2sINF(i)
!		   mo_strain_inf = mo_strain_inf + xp(1,j1)-ax(1,j1)
!		End do
!	  END IF
!	  mo_strain_inf = mo_strain_inf / i_surfINF
!	  mo_strain_sup = mo_strain_sup / i_surfSUP
!	  ef_strainXZ     = (mo_strain_sup - mo_strain_inf)/at(1,1)
!
!	  If (it==1) Open (unit=721, file='strainXZ.dat', status='unknown',action='write')
!	  Write (721,'(i7,2x,e24.16)') it,ef_strainXZ

end subroutine stress_bc





! ---------------------------------------------------------------------
! SUBROUTINE : on fonctionne en CL controlees en contrainte, on
! veut imposer une vitesse de deformation, calcul de la nouvelle vitesse de
! deformation à chaque pas de temps
! ---------------------------------------------------------------------

subroutine recalc_stress
	!-----------------------------------------------
 	!   M o d u l e s
 	!-----------------------------------------------
	USE T_kind_param_m, ONLY: double
	use gen_com_m
	use tab_imm_m
	use eam
	use posana
	
	implicit none
	!-----------------------------------------------
 	!   G l o b a l   P a r a m e t e r s
 	!-----------------------------------------------
 	!-----------------------------------------------
 	!   L o c a l   P a r a m e t e r s
 	!-----------------------------------------------
	!-----------------------------------------------
	!   L o c a l   V a r i a b l e s
	!-----------------------------------------------
	!real(double), external :: calcvol
	integer 	:: i,j1,j2
	real (double)	:: dmoy, dmoy1, dmoy2, gamma_pt
	real (double)	:: sigma_new, inst1, inst2
	logical		:: lecrire

  ! calcul de gamma_pt sur cette iteration
	dmoy1  = 0.
	dmoy2  = 0.
	dmoy   = 0.
	inst1 = 0.
	inst2 = 0.

	Do i=1,i_surfMIN
		j1 = b2sINF(i)
		dmoy1  = dmoy1 - (xp(3,j1) - ax(3,j1))
		inst1  = inst1 - (xp(3,j1) - xpp(3,j1))
		j2 = b2sSUP(i)
		dmoy2  = dmoy2 + (xp(3,j2) - ax(3,j2))
		inst2  = inst2 + (xp(3,j2) - xpp(3,j2))
	End do
	IF (i_surfSUP == i_surfMAX .AND. i_surfINF /= i_surfMAX) THEN
		Do i=i_surfMIN+1,i_surfMAX
		   j2 = b2sSUP(i)
		   dmoy2  = dmoy2 + (xp(3,j2) - ax(3,j2))
		   inst2  = inst2 + (xp(3,j2) - xpp(3,j2))
		End do
	Else if (i_surfINF == i_surfMAX .AND. i_surfSUP /= i_surfMAX) THEN
		Do i=i_surfMIN+1,i_surfMAX
		   j1 = b2sINF(i)
		   dmoy1  = dmoy1 - (xp(3,j1) - ax(3,j1))
		   inst1  = inst1 - (xp(3,j1) - xpp(3,j1))
		End do
	END IF

	dmoy = (dmoy1/i_surfINF + dmoy2/i_surfSUP) / 2

	gamma_pt = 2.*(inst1/i_surfINF + inst2/i_surfSUP)/(Lz_cm*tstep) ! strainrate instantanné

   ! calcul de l'erreur sur gamma_pt - on corrige sigma en consequence
	!currentstress = 0.2e-23
	currentstress = currentstress - fdbkcoef*(gamma_pt - user_strainrate)
	!currentstress en Gpa

   ! on recalcule maitenant les forces à appliquer sur les atomes	
	sigma_new = currentstress*1e+9      !(on a maintenant sigma_new en Pa)

	forceatinf = - sigma_new*layer_surf*1e-04*(y_max - y_min)*1e2*1e3 / (imm * thickness)
	forceatsup =   sigma_new*layer_surf*1e-04*(y_max - y_min)*1e2*1e3 / (imm * thickness)
		! (force en cm.g.s-2)

   ! si on veut : ecriture d'un fichier de sortie pour controle le comportement du materiau	
!	lecrire = .false.
!	if (lecrire == .true.)  then
!	   if (it==2) then
!		open(unit=138,file='hooke_stress_bc.dat',status='unknown',action='write')
		open(unit=131,file='stress_evol.dat',status='unknown',action='write')
!		open(unit=132,file='strain_evol.dat',status='unknown',action='write')
!		open(unit=133,file='strrate_evol.dat',status='unknown',action='write')
!	   end if
!	   Write (138, '(e24.16,3x,e24.16)') 2.*dmoy/Lz_cm, - currentstress
!	   Write (131, '(i5,3x,e24.16)') it, - currentstress
!	   Write (132, '(i5,3x,e24.16)') it, 2.*dmoy/Lz_cm
!	   Write (133, '(i5,3x,e24.16)') it, gamma_pt
!        end if

end subroutine recalc_stress





subroutine it_countdown

	!-----------------------------------------------
 	!   M o d u l e s
 	!-----------------------------------------------
	USE T_kind_param_m, ONLY: double
	use gen_com_m
	use tab_imm_m
	use eam
	use posana
	
	implicit none
	!-----------------------------------------------
 	!-----------------------------------------------

	! CETTE ROUTINE SERT A EDITER LE FICHIER last_iteration.dat, necessaire pour 
	! fonctionner en conditions aux limites controlees en contrainte, avec increment de
	! contrainte

	Open(unit=543, file='last_iteration.dat',status='replace', action='write')
	Write(543,*) it
	Close(543)

end subroutine it_countdown
end module
