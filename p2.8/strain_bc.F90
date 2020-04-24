module strain_bc_mod
   USE gen_com_m, ONLY:ef_strain,gap,i_surfinf,i_surfmax,i_surfmin,i_surfsup,&
        &layer_surf,tstep,user_strainrate,b2sinf,b2ssup,b2ssup,b2sinf,b2sinf,b2ssup,b2ssup,b2sinf,b2sinf,&
        &currentstress,fdbkcoef,forceatinf,forceatsup,i_surfinf,i_surfmax,i_surfmin,i_surfsup,layer_surf,&
        lz_cm,thickness,tstep,user_strainrate,y_max,y_min,ef_strain,forceatinf,forceatsup,gap,i_surfinf,&
        i_surfmax,i_surfmin,i_surfsup,ldyn2d,b2sinf,b2ssup,b2ssup,b2sinf,b2sinf,b2ssup,b2ssup,b2sinf,&
        &b2sinf,b2ssup,b2ssup,b2sinf

  use cryst_to_cart_mod,only:cryst_to_cart


        implicit none
        contains
subroutine strain_bc
	!-----------------------------------------------
 	!   M o d u l e s
 	!-----------------------------------------------
	USE T_kind_param_m, ONLY: double

	USE tab_imm_m
        USE eam
  	USE posana
	
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
	real(double) :: moy_inf, moy_sup, th_strain, ef_strrate
	logical      :: lwrite
	real(double) :: mo_strain_inf, mo_strr_inf
	!-----------------------------------------------

	lwrite = .FALSE.		! pour controle (écriture des parametres de simulation (espilon, epsilon point, ..)

	moy_inf = 0.
	moy_sup = 0.
  ! calcul de la force sur le CDM de chaque couche de surface
	DO i=1,i_surfMIN
		j1 = b2sINF(i)
		moy_inf = moy_inf + fp(3,j1)
		j2 = b2sSUP(i)
		moy_sup = moy_sup + fp(3,j2)
	END DO
	IF (i_surfSUP == i_surfMAX .AND. i_surfINF /= i_surfMAX) THEN
		Do i=i_surfMIN+1,i_surfMAX
		   j2 = b2sSUP(i)
		   moy_sup = moy_sup + fp(3,j2)
		End do
	ELSE IF (i_surfINF == i_surfMAX .AND. i_surfSUP /= i_surfMAX) THEN
		Do i=i_surfMIN+1,i_surfMAX
		   j1 = b2sINF(i)
		   moy_inf = moy_inf + fp(3,j1)
		End do
	END IF

	moy_inf = moy_inf/i_surfINF
	moy_sup = moy_sup/i_surfSUP

  ! on impose V_cdm_z = 0, à chaque itération
	DO i=1,i_surfMIN
		j1 = b2sINF(i)
		fp(3,j1) = fp(3,j1) - moy_inf
		j2 = b2sSUP(i)
		fp(3,j2) = fp(3,j2) - moy_sup
	END DO
	IF (i_surfSUP == i_surfMAX .AND. i_surfINF /= i_surfMAX) THEN
		Do i=i_surfMIN+1,i_surfMAX
		   j2 = b2sSUP(i)
		   fp(3,j2) = fp(3,j2) - moy_sup
		End do
	ELSE IF (i_surfINF == i_surfMAX .AND. i_surfSUP /= i_surfMAX) THEN
		Do i=i_surfMIN+1,i_surfMAX
		   j1 = b2sINF(i)
		   fp(3,j1) = fp(3,j1) - moy_inf
		End do
	END IF

	If (lwrite.EQV..TRUE.) Then	! pour obtenir les lois de comportement
					! ( sigma = f(espilon) par exemple )
	   If (it==1) Open (unit=119, file='strain_it.dat', status='unknown',action='write') ! destiné à etre tracé
	   mo_strain_inf = 0.
	   mo_strr_inf   = 0.
	   Do i=1,i_surfINF
		j1 = b2sINF(i)
		mo_strain_inf = mo_strain_inf + xp(3,j1)-ax(3,j1)
		mo_strr_inf   = mo_strr_inf   + xp(3,j1)-xpp(3,j1)
	   End do
	   mo_strain_inf = - 2.*mo_strain_inf  /(i_surfINF*gap)
	   mo_strr_inf   = - 2.*mo_strr_inf    /(i_surfINF*gap*tstep)
	   ef_strain     = mo_strain_inf
	   ef_strrate    = mo_strr_inf
	   th_strain     = USEr_strainrate * tstep * it
	   Write(119, '(e24.16,2x,e24.16)')  ef_strain, (-(moy_sup*1E-5)*i_surfSUP/layer_surf*1e-04)*1E-9
	End if
end subroutine strain_bc







subroutine surf_calc
	! pour recalculer l'aire des surfaces sur lesquelles on applique nos CL, à 
	! chaque itération
	!-----------------------------------------------
 	!   M o d u l e s
 	!-----------------------------------------------
	USE T_kind_param_m, ONLY: double
	USE gen_com_m, ONLY:
	USE tab_imm_m
        USE eam
  	USE posana
	
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
	real (double) :: x_min, x_max, Lx
	integer :: i
	!-----------------------------------------------

	! calcul de la surface Y=cte
	call cryst_to_cart (imm, xp, bg, -1) !cart vers cryst

  	x_min = xp(1,1)
  	x_max = xp(1,1)
 	Do i=2,imm
   	    If     (xp(1,i) < x_min) Then
		x_min = xp(1,i) 
	    Elseif (xp(1,i) > x_max) Then
		x_max = xp(1,i)
	    End if
 	End do

  	call cryst_to_cart (imm, xp, at, 1)  !cryst vers cart

  	Lx = (x_max-x_min) * sqrt(at(1,1)**2+at(2,1)**2)   

end subroutine surf_calc

end module
