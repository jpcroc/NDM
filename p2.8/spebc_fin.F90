module spebc_fin_mod
  USE gen_com_m, ONLY:b2ssup,y_max,USEr_stress_yz,Lz_cm,thick_cryst,itespebcout,inXMdis1,&
  &tstep,inXMdis2,i_surfSUP,i_surfINF,thick_cryst,b2sSUP,Lz_cm,USEr_stress_yz,y_max,y_2nd_max,y_min,&
       &y_max,y_2nd_max,y_min,USEr_stress_yz,tstep, ef_strain,erg2eV,forceatinf,forceatsup,gap, i_surfINF,&
       &inXMdis1,inXMdis2,b2sINF,tstep,ef_strain,erg2ev,gap,y_max,y_min, USEr_stress_yz,&
       &y_2nd_max, layer_surf, Lx_cm


  implicit none
contains
subroutine spebc_fin (flagfinloc) !(energietotale, flagfinloc)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double

  USE tab_imm_m,only:
#ifdef PARA
  USE mod_para,only:
#endif
  USE fcc_module
  USE cfg_module
  USE posana
  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  real (double) :: moyenne_ener, moyenne_ener2, XP1MAX, XP1MIN, XP2MIN, XP2MAX !energietotale
  integer :: i, compteur, compteur2, j1, numcfg
  integer, dimension(imm) :: indiciz, indiciz2
  logical :: flagfinloc
  character (LEN=19) :: cfgname
  logical :: llaurent, ldislospeed
  real :: XMdislo1, XMdislo2
  integer, dimension(imm) :: tab_bulk
  real(double), dimension(3) :: xprovi
  !-----------------------------------------------

  CALL cryst_to_cart (imm, xp, bg, -1) !cart vers cryst
  CALL cryst_to_cart (imm, xpp, bg, -1)

  IF (flagfinloc.EQV..true.) THEN

 ! GENERE LE FICHIER DES PARAMETRES DE SORTIE ----------------------------------
    Open (unit=388, file='param.out', status='replace', action='write')
      Write(388, '(a,i10)') ' -- derniere iteration ---- ', it
      Write(388, '(a,e24.16)') ' -- Lx -- (Ang) ----------- ', Lx_cm*1e+08
      Write(388, '(a,e24.16)') ' -- Ly -- (Ang) ----------- ', (2*y_max-y_2nd_max-y_min)*1e+08
      Write(388, '(a,e24.16)') ' -- Lz -- (Ang) ----------- ', Lz_cm*1e+08
      Write(388, '(a,e24.16)') ' -- Surface -- (Ang2) ----- ', layer_surf*1e+16
      Write(388, '(a,e24.16)') ' -- Gap - (Ang) ----------- ', gap*1e+08
      Write(388, '(a,e24.16)') ' -- deformation ----------- ', ef_strain
      Write(388, '(a,e24.16)') ' -- SigmaYZ -- (GPa) ------ ', USEr_stress_yz
      Write(388, '(a,i10)') ' -- Nb at dans S+ --------- ', i_surfSUP
      Write(388, '(a,i10)') ' -- Nb at dans S- --------- ', i_surfINF
      Write(388, '(a,e24.16)') ' -- f at dans S+ (eV/Ang) - ', forceatsup*1e-05*1e-10/(1.60217646e-19)
      Write(388, '(a,e24.16)') ' -- f at dans S- (eV/Ang) - ', forceatinf*1e-05*1e-10/(1.60217646e-19)
     ! Write(388, '(a,G10.3)' ) ' -- Etot -- (eV) ---------- ', energietotale
    Close(388)

 ! GENERER UN FICHIER DE SORTIE DANS LE FORMAT UTILISE PAR LE LOGICIEL DE LAURENT/DAVID RODNEY ------------
    llaurent=.FALSE.                          
    IF (llaurent.EQV..TRUE.) THEN			! ATTENTION : non fini, non testé !

      OPEN(unit=366, file='laurent.xyz', status='replace', action='write')
        Write(366,*) imm, USEr_stress_yz*1000
      	Write(366,'(3e20.12)') at(1,1)*1e+8,at(2,2)*1e+8,at(3,3)*1e+8
	tab_bulk(:) = 0
        DO i=1,i_surfINF
           j1 =  b2sINF(i)
	   tab_bulk(j1) = 1
           Write(366,'(a2,4e18.10)') 'Ni',(xp(1,j1)-1/2)*at(1,1)*1e+8 , & 
           (xp(2,j1)-1/2)*at(2,2)*1e+8 , (xp(3,j1)-1/2)*at(3,3)*1e+8 , eatom(j1)*erg2eV	 
        END DO
        DO i=1,i_surfSUP
           j1 =  b2sSUP(i)
	   tab_bulk(j1) = 1
           Write(366,'(a2,4e18.10)') 'Ni',(xp(1,j1)-1/2)*at(1,1)*1e+8 , & 
           (xp(2,j1)-1/2)*at(2,2)*1e+8 , (xp(3,j1)-1/2)*at(3,3)*1e+8 , eatom(j1)*erg2eV	 
        END DO 
        DO i=1,imm
	   IF (tab_bulk(i) == 0) Write(366,'(a2,4e18.10)') 'Si',(xp(1,i)-1/2)*at(1,1)*1e+8 ,&
           (xp(2,i)-1/2)*at(2,2)*1e+8 , (xp(3,i)-1/2)*at(3,3)*1e+8 , eatom(i)*erg2eV
        END DO
      Close(366)

   END IF ! (llaurent)

  END IF !! (flagfinloc==.true.)




 ! GENERE LE FICHIER CFG NE REPRESENTANT QUE LES ATOMES DE LA DISLO -------------
	

	IF (flagfinloc.EQV..true.) THEN
		OPEN(377, file='fin.dislo.cfg', status='replace', action='write')
	ELSE IF (flagfinloc.EQV..false.) THEN
		cfgname = 'dislo.000000000.cfg'
		numcfg  = INT(Real(it)/Real(itespebcout))
		Write(cfgname(7:15), '(i9.9)') numcfg
		Open(377, file=cfgname, status='replace', action='write')
	END IF


	moyenne_ener = 0.
	XP1MAX = xp(1,1)
	XP1MIN = xp(1,1)
	XP2MAX = xp(2,1)
	XP2MIN = xp(2,1)
	DO i=1,imm
		moyenne_ener = moyenne_ener + eatom(i)
		IF (xp(1,i) > XP1MAX) XP1MAX=xp(1,i)
		IF (xp(1,i) < XP1MIN) XP1MIN=xp(1,i)
		IF (xp(2,i) > XP2MAX) XP2MAX=xp(2,i)
		IF (xp(2,i) < XP2MIN) XP2MIN=xp(2,i)
	END DO
	moyenne_ener = moyenne_ener/(real(imm))
	compteur = 0
	moyenne_ener2 = 0.
	compteur2 = 0
	DO i=1,imm
   ! detection des surface grace a leur energie
		IF ( ABS((eatom(i)-moyenne_ener)/moyenne_ener) > 0.02 ) THEN
			compteur = compteur + 1
			indiciz(compteur) = i
   ! detection de la dislo
		! ELSE IF (XP1MAX-xp(1,i) > 0.15 .AND. xp(1,i)-XP1MIN > 0.15) THEN
		ELSE IF (XP1MAX-xp(1,i) > thick_cryst*1.5 .AND. xp(1,i)-XP1MIN > thick_cryst*1.5) THEN
			moyenne_ener2 = moyenne_ener2 + eatom(i)
			compteur2 = compteur2 + 1
		END IF
	END DO
	indiciz(compteur+1) = 0
	moyenne_ener2 = moyenne_ener2 / compteur2
	compteur2 = 0
	DO i=1,imm
		!IF (XP1MAX-xp(1,i) > 0.15 .AND. xp(1,i)-XP1MIN > 0.15) THEN
		IF (XP1MAX-xp(1,i) > thick_cryst*1.5 .AND. xp(1,i)-XP1MIN > thick_cryst*1.5) THEN
		IF (ABS((eatom(i)-moyenne_ener2)/moyenne_ener2) > 0.02 ) THEN
			compteur2 = compteur2 + 1
			indiciz2(compteur2) = i
		END IF
		END IF
	END DO
	indiciz2(compteur2+1) = 0

	Write(377,'(a,i10)')'Number of particles =', compteur+compteur2
	Write(377,'(a)')'A = 1.000 Angstrom (basic length-scale)'
	Write(377,'(a)')'# Unit cell vector #1'
	Write(377,'(a,f25.12,a)')'H0(1,1) = ', at(1,1)*1e+8, '  A'
	Write(377,'(a,f25.12,a)')'H0(1,2) = ', at(2,1)*1e+8, '  A'
	Write(377,'(a,f25.12,a)')'H0(1,3) = ', at(3,1)*1e+8, '  A'
	Write(377,'(a)')'# Unit cell vector #2'
	Write(377,'(a,f25.12,a)')'H0(2,1) = ', at(1,2)*1e+8, '  A'
	Write(377,'(a,f25.12,a)')'H0(2,2) = ', at(2,2)*1e+8, '  A'
	Write(377,'(a,f25.12,a)')'H0(2,3) = ', at(3,2)*1e+8, '  A'
	Write(377,'(a)')'# Unit cell vector #3'
	Write(377,'(a,f25.12,a)')'H0(3,1) = ', at(1,3)*1e+8, '  A'
	Write(377,'(a,f25.12,a)')'H0(3,2) = ', at(2,3)*1e+8, '  A'
	Write(377,'(a,f25.12,a)')'H0(3,3) = ', at(3,3)*1e+8, '  A'
	Write(377,'(a)')'.NO_VELOCITY.'
	Write(377,'(a)')'entry_count = 4'
	Write(377,'(a)')'auxiliary[0] = Energie_pat'
	Write(377,'(a)')'0.000'
	Write(377,'(a)')'Fe'

	DO i=1,compteur
		Write(377,'(4f25.12)') xp(1,indiciz(i)),xp(2,indiciz(i)),xp(3,indiciz(i)), 1.0 !eatom(indiciz(i))*erg2eV
	END DO	
	DO i=1,compteur2
		Write(377,'(4f25.12)') xp(1,indiciz2(i)),xp(2,indiciz2(i)),xp(3,indiciz2(i)), 2.0 !eatom(indiciz(i))*erg2eV
	END DO
      CLOSE(377)



  CALL cryst_to_cart (imm, xp, at, 1)  !cryst vers cart
  CALL cryst_to_cart (imm, xpp, at, 1)  !cryst vers cart



  ! POUR MESURER LA VITESSE DE PROPAGATION DE LA DISLO ----------------------

	ldislospeed=.false.

	IF (ldislospeed.EQV..true.) THEN
	   IF (it==1) Open(unit=446, file='PositionDislo.dat', status='replace', action='write')
	   XMdislo1 = 0.
	   XMdislo2 = 0.
	  DO i=1,compteur2
		XMdislo1 = XMdislo1 + xp(1,indiciz2(i))
		XMdislo1 = XMdislo1 + xp(1,indiciz2(i))
	  END DO
           XMdislo1 = XMdislo1 / Real(compteur2)
	   XMdislo2 = XMdislo2 / Real(compteur2)
	   IF(it==1) THEN
		inXMdis1 = XMdislo1
		inXMdis2 = XMdislo2
	   END IF
	   Write(446,'(2f25.12)') it*tstep*1e+15,  &
           SQRT((XMdislo1-inXMdis1)*(XMdislo1-inXMdis1)+(XMdislo2-inXMdis2)*(XMdislo2-inXMdis2))*1e+8/0.8164970
	END IF



! TEST 2 ecriture OUT.GIN ----------------------------------------------------

	OPEN(399, file='out.gin', status='replace', action='write')
	      Write(399,'(3(i0,1x))') 1, 1, 1
	      Write(399,'(3(g24.16,1x))') at(1:3,1)*1e+8
	      Write(399,'(3(g24.16,1x))') at(1:3,2)*1e+8
	      Write(399,'(3(g24.16,1x))') at(1:3,3)*1e+8
	      Write(399,'(i0)') imm
	      DO i=1, imm
		    WRITE(399,'(3g24.16,2x,i0)') xp(1,i)/at(1,1), xp(2,i)/at(2,2),xp(3,i)/at(3,3), ityp(i)
	      END DO

    	CLOSE(399)

end subroutine spebc_fin








end module
