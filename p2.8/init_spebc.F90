
subroutine init_spebc
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m
  use gen_com_m
  use tab_imm_m
  use eam
  use posana

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  logical :: lwriter
  integer :: i, stock_read
  character (LEN=7) :: f_name
  character (LEN=12):: f_name2
  real(double), dimension(imm) :: b2s_real
  real(double) :: theta, ray_cut
  real(double) :: x_max, x_min, x_2nd_max
  !-----------------------------------------------

  flag_fin = .false.

 ! ------------- Cas CL spécifiques : detection des surfaces du cristal ---------!*!
 ! ------------------------------------------------------------------------------!*!
  lwriter = .true.
  if (ibound==1 .OR. ibound==2 .OR. ibound==3) then				!*!

	! calculation of the gap (needed parameter for v=cte boundary cond.)
	! 	(gap: distance selon Y entre l'interface " couche de surface supérieure / bulk"
	! 	et "couche de surface intérieure / bulk"
  	y_max = xp(2,1)
  	y_min = xp(2,1)
  	Do i = 2,imm
  	   If     (xp(2,i) < y_min) Then
  		y_min = xp(2,i) 
  	   Elseif (xp(2,i) > y_max) Then
  		y_max = xp(2,i)
  	   End if
  	End do

	y_2nd_max = (y_min + y_max)/2
	Do i = 2,imm
  	   If     (xp(2,i) > y_2nd_max .AND. xp(2,i) < y_max) Then
  		y_2nd_max = xp(2,i)
  	   End if
  	End do

  	!thickness = rcut			! our way to define the surface
	!on définit l'épaisseur des couches de surface comme étant le rayon de coupure
	!potentiel d'interaction qu'on utilise
	Open(unit=322, file='eamtab.potin', status='old', action='read')
	Read (322,*) stock_read
	Read (322,*) ray_cut
	thickness = ray_cut * 1e-8
	Close(322)
	
  	gap = (2.*y_max-y_2nd_max-y_min) - 2.*thickness

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
	  x_2nd_max = (x_max+x_min)/2.
	  Do i=2,imm
	   	If (xp(1,i) > x_2nd_max .AND. xp(1,i)<x_max) x_2nd_max = xp(1,i)
	  End do
	call cryst_to_cart (imm, xp, at, 1)  !cryst vers cart

	Lx_cm = (2.*x_max-x_2nd_max-x_min) * sqrt(at(1,1)**2+at(2,1)**2)
	Lz_cm = at(3,3)
	layer_surf = Lx_cm * Lz_cm    ! en cm2



	! calcul des paramètres en prenant en compte l'orientation du cristal
	! (indépendant de l'angle d'une éventuelle rotation d'axe Oz)
  ! PAS UTILISé : IL FAUT UTILISER UN FICHIER D'ENTREE .gin PRENANT DEJA EN COMPTE
  !               LA ROTATION
!	theta = 0.
!	theta = at(1,1)*at(1,2)+at(2,1)*at(2,2)+at(3,1)*at(3,2)
!	theta = theta / sqrt(at(1,1)**2+at(2,1)**2+at(3,1)**2)
!	theta = theta / sqrt(at(1,2)**2+at(2,2)**2+at(3,2)**2)
!	theta = acos (theta)

!	thick_cryst = thickness / (sin(theta))
!	thick_cryst = thick_cryst / sqrt(at(1,2)**2+at(2,2)**2+at(3,2)**2)




	thick_cryst = thickness / sqrt(at(1,2)**2+at(2,2)**2+at(3,2)**2)

	call surface_detect  ! appel routine de détection de la surface

  end if

 ! après avoir détecté la surface, on va appeler les routines d'initilisation
 ! de notre étude qui sont propres à chaque type de conditions aux limites

! ----------- CL strain controlées - deformation imposée   --------------------- !*!
! ----------- Initialisation des parametres et des vitesses   ------------------ !*!




! ----------- Rq : fonctionnement non validé avec le code de David -------------
! -----------      problème identifié : plan de glissement (211) et      -------
! -----------      contrainte de Peierls trop faible
  If (ibound==1) Then 
	Call inisurfspeed
	f_name2 = "in_speed.cfg"
        Call writing (117, f_name2, vp(1,:))
  End if


! ----------- CL stress controlées - contrainte imposée constante   ------------ !*!
! ----------- Initialisation des parametres        ----------------------------- !*!

! ----------- Rq : fonctionnement en partie validé avec le code de David -------
! -----------      à utiliser avec des surfaces libres selon X       -----------
  if (ibound==2) call initstress
	


! ----------- CL stress controlées - deformation imposée   --------------------- !*!
! ----------- Initialisation des parametres     -------------------------------- !*!

! ----------- Rq : fonctionnement non validé
  if (ibound==3)  then
	currentstress = user_stress_yz
	call initstress
  end if

end subroutine init_spebc








subroutine surface_detect

	! *****************************************************
	! THIS PROGRAM DETECT THE ATOMS THAT BELONG TO THE
	! CRYSTAL SURFACE (case ibound = 1 or 2)
	! *****************************************************
	!-----------------------------------------------
 	!   M o d u l e s
 	!-----------------------------------------------
	USE T_kind_param_m
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
	real(double) :: distance, theta, y_max_cry, y_min_cry		! extrema 
	integer :: i
	integer, dimension(imm) :: temp_tabSUP, temp_tabINF
	real, dimension(imm) :: plottab !degub!
	!-----------------------------------------------

 	call cryst_to_cart (imm, xp, bg, -1) !cart vers cryst

    ! calculation of the local and absolute extrema (cryst. coord.)
	y_max_cry = xp(2,1)
	y_min_cry = xp(2,1)
	Do i = 2,imm !natperc
		If     (xp(2,i) < y_min_cry) Then
			y_min_cry = xp(2,i) 
		Elseif (xp(2,i) > y_max_cry) Then
			y_max_cry = xp(2,i)
		End if
	End do

	distance = 0.	
				! new !
	if (y_max_cry-y_min_cry <= 2 * thick_cryst) write(6,*) '########### Cellule_trop_fine'
	i_surfINF = 0
	i_surfSUP = 0
	Do i = 1,imm
	   plottab(i)=0. !debug!
	   distance = xp(2,i)-y_min_cry		
	   If (distance < thick_cryst) then
		i_surfINF = i_surfINF + 1
		temp_tabINF(i_surfINF) = i
		plottab(i)=1. !debug!
	   End if	   
   	   distance = y_max_cry - xp(2,i)	   
   	   If (distance < thick_cryst) Then
		i_surfSUP = i_surfSUP + 1
		temp_tabSUP(i_surfSUP) = i
		plottab(i)=2. !debug!
	   Endif   
	   distance = 0.
	End do			

	allocate (b2sSUP(i_surfSUP))
	allocate (b2sINF(i_surfINF))
	i_surfMAX = max(i_surfSUP,i_surfINF)
	i_surfMIN = min(i_surfSUP,i_surfINF)

	Do i=1,i_surfMIN
		b2sSUP(i) = temp_tabSUP(i)
		b2sINF(i) = temp_tabINF(i)
	End do
	IF (i_surfSUP == i_surfMAX .AND. i_surfINF /= i_surfMAX) THEN
		Do i=i_surfMIN+1,i_surfMAX
		   b2sSUP(i) = temp_tabSUP(i)
		End do
	Else if (i_surfINF == i_surfMAX .AND. i_surfSUP /= i_surfMAX) THEN
		Do i=i_surfMIN+1,i_surfMAX
		   b2sINF(i) = temp_tabINF(i)
		End do
	END IF			! new !
				! new !	
				! new !

	! ---------------------------------------------
	! ---------------------------------------------

	
!	Open(623, file='surfaces.cfg', status='replace', action='write') !debug!
!
!	Write(623,'(a,i)')'Number of particles =', imm
!	Write(623,'(a)')'A = 1.000 Angstrom (basic length-scale)'
!	Write(623,'(a)')'# Unit cell vector #1'
!	Write(623,'(a,f,a)')'H0(1,1) = ', at(1,1)*1e+8, '  A'
!	Write(623,'(a,f,a)')'H0(1,2) = ', at(2,1)*1e+8, '  A'
!	Write(623,'(a,f,a)')'H0(1,3) = ', at(3,1)*1e+8, '  A'
!	Write(623,'(a)')'# Unit cell vector #2'
!	Write(623,'(a,f,a)')'H0(2,1) = ', at(1,2)*1e+8, '  A'
!	Write(623,'(a,f,a)')'H0(2,2) = ', at(2,2)*1e+8, '  A'
!	Write(623,'(a,f,a)')'H0(2,3) = ', at(3,2)*1e+8, '  A'
!	Write(623,'(a)')'# Unit cell vector #3'
!	Write(623,'(a,f,a)')'H0(3,1) = ', at(1,3)*1e+8, '  A'
!	Write(623,'(a,f,a)')'H0(3,2) = ', at(2,3)*1e+8, '  A'
!	Write(623,'(a,f,a)')'H0(3,3) = ', at(3,3)*1e+8, '  A'
!	Write(623,'(a)')'.NO_VELOCITY.'
!	Write(623,'(a)')'entry_count = 4'
!	Write(623,'(a)')'auxiliary[0] = surface_atoms'
!	Write(623,'(a)')'0.000'
!	Write(623,'(a)')'Fe'
!
!	DO i=1,imm
!		Write(623,'(f,f,f,f)') xp(1,i),xp(2,i),xp(3,i),plottab(i) 
!	END DO
!
!        CLOSE(623)     !debug!

	! ---------------------------------------------
	! ---------------------------------------------

	call cryst_to_cart (imm, xp, at, 1)  !cryst vers cart	

end subroutine surface_detect






subroutine inisurfspeed
  ! initialisation dans le cas ibound == 1 (CL controlees en deformation)

  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m
  use gen_com_m
  use tab_imm_m
  use eam
  use posana

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  real 	  :: v_sup, v_inf, v_mil
  integer :: i, j1, j2
  integer, dimension(imm) :: b2bulk
  !-----------------------------------------------

! calculation of surface atom speed (needed param. for v=cte boundary cond.)
  speed_user = gap*user_strainrate/(2.)

! surface layer atom speed initialization
  v_sup = 0.
  v_inf = 0.
  v_mil = 0.

  b2bulk(:) = 1

  Do i=1,i_surfMIN
	j1 = b2sINF(i)
	v_inf = v_inf + vp(3,j1)
	b2bulk(j1) = 0
	j2 = b2sSUP(i)
	v_sup = v_sup + vp(3,j2)
	b2bulk(j2) = 0
  End do
  IF (i_surfSUP == i_surfMAX .AND. i_surfINF /= i_surfMAX) THEN
	Do i=i_surfMIN+1,i_surfMAX
	   j2 = b2sSUP(i)
	   v_sup = v_sup + vp(3,j2)
	   b2bulk(j2) = 0
	End do
  Else if (i_surfINF == i_surfMAX .AND. i_surfSUP /= i_surfMAX) THEN
	Do i=i_surfMIN+1,i_surfMAX
	   j1 = b2sINF(i)
	   v_inf = v_inf + vp(3,j1)
	   b2bulk(j1) = 0
	End do
  END IF
								
  Do i=1,imm
	if (b2bulk(i) == 1) v_mil = v_mil + vp(3,i)
  End do

  v_sup = v_sup / i_surfSUP
  v_inf = v_inf / i_surfINF
  v_mil = v_mil / (imm - i_surfSUP - i_surfINF)

  ! initialisation des vitesses et du tableau xpp ("position à l'itération précédente)
  ! Rmq : on veut v_cdm = 0 pour chaque surface, ainsi que pour le bulk

  Do i=1,i_surfMIN
	j1 = b2sINF(i)
	vp(3,j1) = vp(3,j1) - v_inf - speed_user
        xpp(3,j1)= xp(3,j1) - vp(3,j1)*tstep
	j2 = b2sSUP(i)
        vp(3,j2) = vp(3,j2) - v_sup + speed_user
	xpp(3,j2)= xp(3,j2) - vp(3,j2)*tstep
  End do
  IF (i_surfSUP == i_surfMAX .AND. i_surfINF /= i_surfMAX) THEN
	Do i=i_surfMIN+1,i_surfMAX
	   j2 = b2sSUP(i)
           vp(3,j2) = vp(3,j2) - v_sup + speed_user
	   xpp(3,j2)= xp(3,j2) - vp(3,j2)*tstep
	End do
  Else if (i_surfINF == i_surfMAX .AND. i_surfSUP /= i_surfMAX) THEN
	Do i=i_surfMIN+1,i_surfMAX
	   j1 = b2sINF(i)
	   vp(3,j1) = vp(3,j1) - v_inf - speed_user
           xpp(3,j1)= xp(3,j1) - vp(3,j1)*tstep
	End do
  END IF
								
  Do i=1,imm
	if (b2bulk(i) == 1) then
	     vp(3,i) = vp(3,i) - v_mil
	     xpp(3,i)= xp(3,i) - vp(3,i)*tstep
	end if
  End do

end subroutine inisurfspeed






subroutine initstress
   !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m
  use gen_com_m
  use tab_imm_m

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  real (double) :: sigma_user
  integer :: i, j1, j2
  !-----------------------------------------------

  ! l'utilisateur rentre simga (en GPa)
  sigma_user = user_stress_yz*1e+9      !(on a maintenant sigma_user en Pa)

  forceatinf = - (sigma_user*layer_surf*1e-04 / i_surfINF) 
  forceatinf = forceatinf *1e2*1e3
  forceatsup =   (sigma_user*layer_surf*1e-04 / i_surfSUP)
  forceatsup = forceatsup *1e2*1e3 ! calcul des forces atomiques

	! (force en cm.g.s-2)

  IF (ldyn2D==.true.) THEN   	! initialisation vp(2,i)=0. pour les atomes en surface
			 	! (cas DYN 2D)
	  Do i=1,i_surfMIN
		j1 = b2sINF(i)
		vp(2,j1) = 0.
		xpp(2,j1)=xp(2,j1)
		j2 = b2sSUP(i)
		vp(2,j2) = 0.
		xpp(2,j2)=xp(2,j2)
	  End do
	  IF (i_surfSUP == i_surfMAX .AND. i_surfINF /= i_surfMAX) THEN
		Do i=i_surfMIN+1,i_surfMAX
		   j2 = b2sSUP(i)
		   vp(2,j2) = 0.
		   xpp(2,j2)=xp(2,j2)
		End do
	  Else if (i_surfINF == i_surfMAX .AND. i_surfSUP /= i_surfMAX) THEN
		Do i=i_surfMIN+1,i_surfMAX
		   j1 = b2sINF(i)
		   vp(2,j1) = 0.
		   xpp(2,j1)=xp(2,j1)
		End do
	  END IF

  END IF

end subroutine initstress



