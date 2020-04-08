MODULE FireModule

  ! Quench using FIRE algorithm
  ! Ref.: Bitzek, E., Koskinen, P., Gähler, F., Moseler, M., and Gumbsch, P.
  !       "Structural Relaxation Made Simple"
  !       Phys. Rev. Lett. 97, 170201 (2006).

  USE T_kind_param_m, ONLY:  double
  USE period_mod
  USE gen_com_m, ONLY:imm,im,lperiod,tstep,usdh,tstep
  ! --- Paramètres de l'algorithme fire -----------------------
  real(double), parameter, private :: finc=1.1
  real(double), parameter, private :: fdec=0.5
  real(double), parameter, private :: alph_start=0.1
  real(double), parameter, private :: f_alph=0.99
  real(double), parameter, private :: tstep_MM=10
  integer, parameter, private:: nStepMin=5

CONTAINS

SUBROUTINE init_trempe_fire(dt, nstep, alph)



  implicit none
  REAL(double), intent(out) :: dt
  INTEGER, intent(out) :: nstep
  REAL(double), intent(out) :: alph
  

  alph = alph_start
  nstep = 0
  dt = tstep

END SUBROUTINE init_trempe_fire

! **************************************************************
subroutine trempe_fire(xp, xpp, vp, ax, fp, ielat, iwmax, ityp, &
        dt, nstep, alph)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE gen_com_m, ONLY:
  USE var_pot, ONLY:ntyp,cm
  implicit none

  integer  :: ielat(imm)
  integer  :: iwmax(imm)
  integer  :: ityp(imm)
  real(double)  :: xp(3,imm)
  real(double)  :: xpp(3,imm)
  real(double)  :: vp(3,imm)
  real(double)  :: ax(3,imm)
  real(double)  :: fp(3,imm)
  REAL(double), intent(inout) :: dt
  INTEGER, intent(inout) :: nstep
  REAL(double), intent(inout) :: alph


  real(double):: norme_de_fp, norme_de_vp, pscal
  integer::i

  real(double), dimension(ntyp) :: aux
  real(double), dimension(1:3) :: xprov

  ! 1/ Intégration de l'équation de mouvement
  aux(:ntyp) = dt**2/(2.d0*cm(:ntyp))
  usdh = 1.d0/(2.d0*dt)
  DO i=1, im
     xprov(:) = xp(:,i) + vp(:,i)*dt + fp(:,i)*aux(iTyp(i))   
     vp(:,i) = (xprov(:) - xpp(:,i))*usdh
     xpp(:,i) = xp(:,i)
     xp(:,i) = xprov(:)
  END DO

  
  IF (lperiod) call period


  ! 2/ Renormalisation des vitesses par l'algorithme fire

  ! Puissance dissipée
  pScal = Sum( vp(:,1:im)*fp(:,1:im) )

  ! Modification du vecteur vitesse
  if (pScal.gt.0) then

          ! Norme du vecteur force
          norme_de_fp = Sqrt( Sum( fp(:,1:im)**2 ) )

          ! Norme du vecteur vitesse
          norme_de_vp = Sqrt( Sum( vp(:,1:im)**2 ) )

          ! Nouveau vecteur vitesse
          vp(:,1:im) = (1.d0-alph)*vp(:,1:im) + alph*norme_de_vp/norme_de_fp*fp(:,1:im)

          nStep = nStep + 1
          if (nStep.gt.nStepMin) then
                  dt=min(dt*finc,tstep_MM*tstep)
                  alph=alph*f_alph
          end if
  else
          vp(:,:)=0.
          dt=dt*fdec
          alph=alph_start
          nstep=0
  end if
  !write(6,'(A,4g14.5)')'FIRE p, dt, alpha, v: ', pscal,dt,alph,norme_de_vp
  return
end subroutine trempe_fire

END MODULE FireModule

