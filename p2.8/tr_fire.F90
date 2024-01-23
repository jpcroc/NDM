MODULE FireModule
!PRL 97, 170201 (2006)
  ! Quench using FIRE algorithm
  ! Ref.: Bitzek, E., Koskinen, P., Gähler, F., Moseler, M., and Gumbsch, P.
  !       "Structural Relaxation Made Simple"
  !       Phys. Rev. Lett. 97, 170201 (2006).
  use gen_com_m,only:iteration
  USE T_kind_param_m, ONLY:  double
  use atomconfig,only:atom_config_d
  ! --- Paramètres de l'algorithme fire -----------------------
  real(double), parameter :: finc=1.1
  real(double), parameter :: fdec=0.5
  real(double), parameter :: alph_start=0.1
  real(double), parameter :: f_alph=0.99
  real(double), parameter :: tstep_MM=10
  integer, parameter:: nStepMin=5

  real(double)::tstep0
  
CONTAINS

SUBROUTINE init_trempe_fire(dt, nstep, alph)



  implicit none
  REAL(double), intent(out) :: dt
  INTEGER, intent(out) :: nstep
  REAL(double), intent(out) :: alph
  
  tstep0=dt
  alph = alph_start
  nstep = 0

  

END SUBROUTINE init_trempe_fire

! **************************************************************
subroutine trempe_fire(atdml, dt, nstep, alph)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------

  USE var_pot, ONLY:ntyp,cm
  implicit none
  class(atom_config_d)::atdml
    integer::im
  REAL(double), intent(inout) :: dt
  INTEGER, intent(inout) :: nstep
  REAL(double), intent(inout) :: alph


  real(double):: norme_de_fp, norme_de_vp, pscal,usdh
  integer::i

  real(double), dimension(ntyp) :: aux
  real(double), dimension(1:3) :: xprov
  im=atdml%im
  ! 1/ Intégration de l'équation de mouvement
  aux(:ntyp) = dt**2/(2.d0*cm(:ntyp))
  usdh = 1.d0/(2.d0*dt)
  DO i=1, im
     xprov(:) = atdml%xp(:,i) + atdml%vp(:,i)*dt + atdml%fp(:,i)*aux(atdml%iTyp(i))   
     atdml%vp(:,i) = (xprov(:) - atdml%xpp(:,i))*usdh
     atdml%xpp(:,i) = atdml%xp(:,i)
     atdml%xp(:,i) = xprov(:)
  END DO

  ! 2/ Renormalisation des vitesses par l'algorithme fire
  ! Puissance dissipée
  pScal = Sum( atdml%vp(:,1:im)*atdml%fp(:,1:im) )
!  write(6,*)'PSCAL',iteration,nstep,pscal,dt
  ! Modification du vecteur vitesse
  if (pScal.gt.0) then
          ! Norme du vecteur force
          norme_de_fp = Sqrt( Sum( atdml%fp(:,1:im)**2 ) )
          ! Norme du vecteur vitesse
          norme_de_vp = Sqrt( Sum( atdml%vp(:,1:im)**2 ) )
          ! Nouveau vecteur vitesse
          atdml%vp(:,1:im) = (1.d0-alph)*atdml%vp(:,1:im) + alph*norme_de_vp/norme_de_fp*atdml%fp(:,1:im)
          nStep = nStep + 1
          if (nStep.gt.nStepMin) then
             dt=min(dt*finc,tstep_MM*tstep0)
             alph=alph*f_alph
          end if
  else
          atdml%vp(:,:)=0.
          dt=dt*fdec
          alph=alph_start
          nstep=0
  end if
  !write(6,'(A,4g14.5)')'FIRE p, dt, alpha, v: ', pscal,dt,alph,norme_de_vp
  return
end subroutine trempe_fire

END MODULE FireModule

