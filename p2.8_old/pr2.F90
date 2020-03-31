module Parrinello_Rahman

  ! Algorithme de Parrinello-Rahman [1] combiné ou non à un thermostat de Nosé-Hoover [4]
  ! basé sur les article de Ray et Rahman [2,3]
  ! Si le thermostat de Nosé-Hoover n'est pas branché (lTHoover=.false.), la
  ! viscosité zHoover restera nulle. 
  !
  ! Quelques remarques
  ! * La vitesse des particules ne tient pas compte de la dérivée du tenseur h:
  ! * On impose une tension thermodynamique plutôt qu'une contrainte constante:
  !   eq. (2.22) de [2]
  !   Du coup, la matrice h0 définissant l'état de référence n'apparaît plus
  !   nulle part dans l'algorithme. Cette matrice est nécessaire seulement si à
  !   un instant donné on souhaite calculer la déformation epsilon de la boîte
  !   et également son énergie potentielle Ucell, mais ces 2 quantités ne sont
  !   pas nécessaires à l'algorithme.
  ! * Le thermostat est imposé sur l'énergie cinétique des atomes et de la
  !   cellule.
  ! * L'algorithme utilisé pour intégrer les équations de mouvement est
  !   "velocity Verlet". Les positions et les vitesses sont donc correctes
  !   jusqu'à respectivement l'ordre 3 et 2 inclus.
  ! * Pour que l'algorithme soit correcte, les forces ne doivent pas faire
  !   intervenir la vitesse des particules. Des thermostats du type lTcon=.true.
  !   ou lTBerendsen=.true. ne sont donc pas autorisés.
  !  
  !
  ! [1] Parrinello, M. & Rahman, 
  !     A. Polymorphic Transitions in Single rcystals: A New Molecular Dynamics Method
  !     J. Appl. Phys., 1981, 52, 7182-7190
  ! [2] Ray, J.R. & Rahman, A. 
  !     Statistical Ensembles and Molecular Dynamics Studies of Anisotropic Solids
  !     J. Chem. Phys., 198bulk.crcin4, 80, 4423-4428
  ! [3] Ray, J.R. & Rahman, A. 
  !     Statistical Ensembles and Molecular Dynamics Studies of Anisotropic Solids. II 
  !     J. Chem. Phys., 1985, 82, 4243-4247
  ! [4] Nosé, S. 
  !     A Molecular Dynamics Method for Simulations in the Canonical Ensemble
  !     Mol. Phys., 1984, 52, 255-268

  USE T_kind_param_m
  use gen_com_m    
  use var_pot
  use mat_util
#if(PARA)
  use mod_mpi
#endif
  implicit none
  ! Vecteurs de la boîte et leurs dérivées
  real(double), dimension(3,3), save , private :: h, hDot
  real(double), dimension(3,3), save , private :: trh, invh, invtrh, Gmat, invGmat, Gdot
  real(double), save, private :: invVolu

  ! Coordonnées réduites des atomes et leurs dérivées
  real(double), pointer, save, private :: sp(:,:), sdot(:,:), sdot_new(:,:)

  ! Variable associée au thermostat de Nosé-Hoover
  !  (zHoover est défini dans gen_com_m.F90)
  REAL(double), dimension(:), allocatable, save, private :: wHoover   ! Poids associé au thermostat de Hoover
  REAL(double), dimension(:), allocatable, save, private :: zOld, zNew, zDot !  Viscosité et dérivée
  REAL(double), dimension(:), allocatable, save, private :: KHoover
  !REAL(double), dimension(:), allocatable, save, private :: UHoover, UHoover_new, UHoover_old Énergies

  ! Nombre de degrés de liberté
  REAL(double), save, private :: gNose

  ! Variables uniquement nécessaires au calcul de l'énergie potentielle de la
  ! boîte
  real(double), dimension(3,3), save , private ::trh0,invh0,invtrh0,epsi, tension
  real(double), save , private ::volu0, invVolu0


contains

  subroutine initlpr (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

    implicit none
    ! Variables utiles
    integer, intent(in)  :: ityp(imm)
    real(double), intent(inout)  :: xp(3,imm), xpp(3,imm), vp(3,imm), fp(3,imm)
    ! Variables inutiles
    integer  :: ielat(imm), iwmax(imm)
    real(double) :: ax(3,imm)

    INTEGER :: ia, i, j
    real(double), external :: calcvol
#if(PARA)
    real(double)::wbox_tot
    real(double) sigkine_tot(3,3)
#endif


    IF(RANG==0) WRITE(6,*)
    IF(RANG==0) WRITE(6,*) 'Algorithme de Parrinello-Rahman (V2)'
    IF(RANG==0) WRITE(6,'(a)') '  -> la vitesse de la boîte ne prend pas en compte la dérivée du tenseur h'
    IF(RANG==0) WRITE(6,*)

    IF (lUcell) THEN
       IF(RANG==0) WRITE(6,'(a)') "Repère de référence pour Parrinello-Rahman  (A):"
       IF(RANG==0) WRITE(6,'(a,3(f0.5,1x))') ' h0(1:3,1) = ', 1e8*h0(1:3,1)
       IF(RANG==0) WRITE(6,'(a,3(f0.5,1x))') ' h0(1:3,2) = ', 1e8*h0(1:3,2)
       IF(RANG==0) WRITE(6,'(a,3(f0.5,1x))') ' h0(1:3,3) = ', 1e8*h0(1:3,3)
       IF(RANG==0) WRITE(6,*)
    ELSE
       h0 = at
    END IF
    IF(RANG==0) WRITE(6,'(a)') "Repère actuel  (A):"
    IF(RANG==0) WRITE(6,'(a,3(f0.5,1x))') ' h (1:3,1) = ', 1e8*at(1:3,1)
    IF(RANG==0) WRITE(6,'(a,3(f0.5,1x))') ' h (1:3,2) = ', 1e8*at(1:3,2)
    IF(RANG==0) WRITE(6,'(a,3(f0.5,1x))') ' h (1:3,3) = ', 1e8*at(1:3,3)
    IF (wbox==0.0) THEN
       wbox = sum(0.5*cm(ityp(:im)))       ! La moitié de la masse totale des atomes
#if(PARA)
       call MPI_ALLREDUCE(wbox,wbox_tot,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
       wbox=wbox_tot
#endif


    END IF
    IF(RANG==0) WRITE(6,'(a,g20.12)')'Masse de la boîte pour Parrinello-Rahman: wbox=',wbox

    ! État de référence défini par la matrice h0
    !   Cet état de référence doit correspondre à un tenseur de contrainte nul.
    !   Il n'est utile que pour calculer la déformation et l'énergie potentielle
    !   de la boîte.
    volu0 = calcvol(h0(1:3,1),h0(1:3,2),h0(1:3,3))
    invVolu0 = 1.d0/volu0
    trh0=Transpose(h0)
    CALL MatInv(h0,invh0)
    invtrh0=Transpose(invh0)

    ! Vecteurs de la boîte et grandeurs associées à l'instant initial
    h(:,:)=at(:,:)
    trh=Transpose(h)
    Gmat = MatMul(trh,h)
    CALL MatInv(Gmat,invGmat)
    call MatInv(h,invh)
    invtrh = Transpose(invh)
    volu = calcvol(h(1:3,1),h(1:3,2),h(1:3,3))
    invVolu = 1.d0/volu

    ! Initialisation de la vitesse de la boîte
    IF(RANG==0) WRITE(6,'(a,f0.3,a)') 'Initialisation de la vitesse de la boîte pour la température ', 0.d0, ' K'
    hdot(:,:) = 0.d0
    ! Time derivative of the metric tensor Gmat
    DO i=1, 3
       DO j=1, 3
          hdot(1:3,1:3)=hdot(1:3,1:3)*ihbox0(1:3,1:3)
          Gdot(i,j) = Sum( h(1:3,i)*hdot(1:3,j) + hdot(1:3,i)*h(1:3,j) )
       END DO
    END DO

    ! Kinetic energy of the cell (Eq. 2.14 of Ref. [2])

    Kcell = 0.5d0*wbox*Sum( hDot(1:3,1:3)**2 )


    if(rang==0) write(6,'(I7,D10.3,A,D21.12,A,a,f0.3,a)') 0,0.d0,'*Kcell = ',Kcell*unitE,cunitE, &
         '  (', 2.d0*Kcell/(9.d0*bk), ' K)'

    ! Nombre de thermostats de Hoover
    IF (nHoover.LT.0) nHoover=0 
    IF (.NOT.lTHoover) nHoover=0
    ALLOCATE(zHoover(1:nHoover+1))


    IF (lTHoover) THEN

       ! Allocation des tableaux nécessaires au thermostat de Hoover 
       IF (nHoover.LE.0) THEN
          WRITE(0,'(a)') 'Si le thermostat de Nosé-Hoover est utilisé (lTHoover=.true.), '
          WRITE(0,'(a)') 'le nombre de thermostats doit être supérieur à 0'
          WRITE(0,'(a,i0)') '  nHoover = ', nHoover
          STOP '< init_lpr >'
       END IF
       ALLOCATE(zOld(1:nHoover), zNew(1:nHoover), zDot(1:nHoover), wHoover(1:nHoover), &
            KHoover(1:nHoover))
       !ALLOCATE(UHoover(1:nHoover), UHoover_new(1:nHoover), UHoover_old(1:nHoover))

       ! Nombre de degrés de liberté pour le thermostat de Nosé-Hoover
       !crc       gNose=dble(3*imana)
       gNose=dble(3*im_glob)

       ! Masse de chaque thermostat
       IF (wNose.EQ.0) THEN
          ! On veut qu'une variation de la température de 10K corresponde à
          ! une variation de fNose de 1% 
          wNose = gNose*bk*10.d0*tstep**2/1.d-2**2
       END IF
       wHoover(1)=wNose
       DO i=2, nHoover
          wHoover(i)=wNose/gNose
       END DO

       DO i=1, nHoover
          KHoover(i) = 0.5d0*bk*Text
          zHoover(i) = Sqrt( 2.d0*KHoover(i)/wHoover(i) )
          zDot(i) = 0.d0
          zOld(i) = zHoover(i) - tstep*zdot(i)
          !UHoover(i) = 0.d0
          !UHoover_old(i) = UHoover(i) - tstep*bk*Text*zHoover(i)   
       END DO
       !UHoover(1) = gNose*UHoover(1)
       KNose = Sum( KHoover(1:nHoover) )
       !UNose = Sum( UHoover(1:nHoover) )
       zHoover(nHoover+1)=0.d0

       IF(RANG==0) WRITE(6,*)
       IF(RANG==0) WRITE(6,'(a)') 'Thermostat de Nosé-Hoover (V2)'
       IF(RANG==0) WRITE(6,'(a)') "  -> l'énergie cinétique des atomes et de la boîte est thermalisée"
       IF(RANG==0) WRITE(6,'(a,g20.12)')'Masse de la boîte pour thermostat de Nosé-Hoover: wHoover=',wNose
       IF(RANG==0) WRITE(6,'(a,g20.12)')'Nombre de degrés de liberté: gNose=',gNose              
       IF(RANG==0) WRITE(6,'(a,i0)')    'Nombre de thermostats: nHoover=', nHoover
       IF(RANG==0) WRITE(6,'(a,f0.3,a)') 'Initialisation du thermostat de Nosé-Hoover pour la température ', &
            2.d0*KNose/(bk*dble(nHoover)), ' K'
       IF(RANG==0) WRITE(6,'(I7,D10.3,A,D21.12,A,a,f0.3,a)') it,timel,'*KNose = ',KNose*unitE,cunitE, &
            '  (', 2.d0*KNose/(bk*nHoover), ' K)'
       IF(RANG==0) WRITE(6,*)
    ELSE
       gNose=0.d0; zHoover(1)=0.d0 
    END IF

    ALLOCATE(sp(1:3,1:imm), sdot(1:3,1:imm), sdot_new(1:3,1:imm))
    ! Coordonnées réduites des atomes et leurs dérivées à l'instant initial
    sp(:,1:im) = MatMul(invh(:,:), xp(:,1:im) )
    sdot(:,1:im) = MatMul(invh(:,:), vp(:,1:im) )

    ! Forces à l'instant initial
    CALL CalFo (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

    !  Contrainte thermique à l'instant initial
    sigkine(:,:)=0.d0
    do ia = 1, im
       do j = 1,3
          sigkine(1:3,j) = sigkine(1:3,j) + cm(ityp(ia))*vp(1:3,ia)*vp(j,ia)
       enddo
    enddo
    sigkine(1:3,1:3) = invVolu*sigkine(1:3,1:3)

#if(PARA)
    call MPI_ALLREDUCE(sigkine,sigkine_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
    sigkine=sigkine_tot

#endif

    ! Énergie cinétique des atomes à l'instant initial
    kine = 0.5d0*volu*( sigKine(1,1) + sigKine(2,2) + sigKine(3,3) )

    ! Contrainte totale à l'instant initial
    sigtot = 0.5d0*(sigkine + Transpose(sigkine) + sig + Transpose(sig) )

    RETURN

  end subroutine initlpr

  !-----------------------------------------------

  subroutine pr (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

    implicit none
    ! Variables utiles
   
    integer, intent(in)  :: ityp(imm)
    real(double), intent(inout)  :: xp(3,imm), xpp(3,imm), vp(3,imm), fp(3,imm)
    ! Variables inutiles
    integer  :: ielat(imm), iwmax(imm)
    real(double) :: ax(3,imm)

    real(double),dimension(3,3)::mf,mfi, grsig, hdot_new, hdot_last,forcebox
    REAL(double) :: diff, tdiff
    integer:: i,j,ia, iter,ic
    real(double) , external ::  calcvol 
    ! Parameter for Parrinello-Rahman self consistency loop
    REAL(double), parameter :: tol=1.0d-12        ! Tolerance for h convergency
    INTEGER, parameter :: max_Iter=100            ! Maximal number of iterations in self-consistency loop
#if(PARA)
    real(double)::wbox_tot
    real(double) sigkine_tot(3,3)
  integer :: nb1, nb2, nb3, i1, l,noxn,noyn,nozn
  real(double) :: zlx, zly, zlz, ux, uy, uz,  pi2, fact, fact1&
       , fact2, hk2, ex, ex1, ex2


#endif

  if (lprtrp) then
     do i = 1, im
        do ic = 1, 3
           if (vp(ic,i)*fp(ic,i)<0) then
              vp(ic,i)=0.
           end if
        end do
     end do
!     do i = 1, im
!        do ic = 1, 3
!           if (vp(ic,i)*fp(ic,i)<0) then
!              vp(ic,i)=0.
!           end if
!        end do
!     end do
     forcebox(:,:)=MatMul( sigtot(:,:) - sigext(:,:), invtrh(:,:) )
     do i = 1, 3
        do ic = 1, 3
           if (hdot(ic,i)*forcebox(ic,i)<0) then
              hdot(ic,i)=0.
           end if
        end do
     end do
  end if

    ! Coordonnées réduites des atomes (au cas où elles ont été modifiées à l'extérieur)
    sp(:,1:im) = MatMul(invh(:,:), xp(:,1:im) )
    ! De même pour les vitesses au cas où, par exemple, on utilise le thermostat 

    sdot(:,1:im) = MatMul(invh(:,:), vp(:,1:im) )

    ! Dérivée des coordonnées réduites des atomes à l'instant t+dt/2
    mf(:,:) = -0.5d0*tstep*MatMul(invGmat,Gdot)
    DO i=1, 3
       mf(i,i) = 1.d0 - 0.5d0*tstep*zHoover(1) + mf(i,i)
    END DO
    DO ia=1, im
       sdot(:,ia) = MatMul( mf(:,:), sdot(:,ia) ) &
            + tstep/(2.d0*cm(ityp(ia))) * MatMul( invh(:,:), fp(:,ia) ) 
    END DO

    ! Dérivée du tenseur h à l'instant t+dt/2
    IF (lpcon2.EQV..true.) THEN    ! On ajoute une force de friction
       hdot(:,:) = ( 1.d0 - 0.5d0*tstep*zHoover(1) - 0.5d0/tbox )* hdot(:,:)*ihbox0(:,:) &
            + tstep/(2.d0*wBox)*volu*MatMul( sigtot(:,:) - sigext(:,:), invtrh(:,:) )*ihbox0(:,:)
    ELSE        ! Équation sans force de friction supplémentaire
       hdot(:,:) = ( 1.d0 - 0.5d0*tstep*zHoover(1) )* hdot(:,:)*ihbox0(:,:) &
            + tstep/(2.d0*wBox)*volu*MatMul( sigtot(:,:) - sigext(:,:), invtrh(:,:) )*ihbox0(:,:)
    END IF

    ! Coordonnées réduites des atomes à l'instant t+dt
    sp(:,1:im) = sp(:,1:im) + sdot(:,1:im)*tstep

    ! Tenseur h à l'instant t+dt
    h(:,:) = h(:,:) + hdot(:,:)*tstep*ihbox0(:,:)

    ! Coordonnées réelles à l'instant t+dt
    xpp(:,1:im) = xp(:,1:im)
    xp(:,1:im) = MatMul( h, sp(:,1:im) )
    at(:,:) = h(:,:)                            ! Vecteur de périodicité
    call recips (h(:,1),h(:,2),h(:,3), bg(:,1),bg(:,2),bg(:,3)) ! Vecteurs réciproques
    volu = calcvol(h(1:3,1),h(1:3,2),h(1:3,3))  ! Volume
    invVolu = 1.d0/volu
    trh=Transpose(h)                            ! Matrices associées à h
    Gmat = MatMul(trh,h)
    CALL MatInv(Gmat,invGmat)
    call MatInv(h,invh)
    invtrh = Transpose(invh)


#if(PARA)
     zl(1) = Sqrt( Sum(at(1:3,1)**2 ) )
     zl(2) = Sqrt( Sum(at(1:3,2)**2 ) )
     zl(3) = Sqrt( Sum(at(1:3,3)**2 ) )
     volu=calcvol(at(1:3,1),at(1:3,2),at(1:3,3))
     zls2(1:3) = 0.5d0*zl(1:3)


     call caltabt
!  temps_debpara=MPI_Wtime()
	   ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs
  call maj_atomes_frt_ftm
!  temps_para=temps_para+MPI_Wtime()-temps_debpara


  if (iewald>0) then

     ! --- Tableaux des troisiemes termes de la sommation d'Ewald ---
     auxe = 23.06134575D-20                  ! en erg.cm (charge electron^2/4*pi*permitivite vide)
     pi2 = pi*pi
     volu=calcvol(at(1:3,1),at(1:3,2),at(1:3,3))
     fact = pi2/alpha**2
     fact1 = auxe/2./pi/volu
     fact2 = auxe*2./volu
     do nb1 = -ncoucx, ncoucx
        do nb2 = -ncoucy, ncoucy
           do nb3 = -ncoucz, ncoucz
              if (nb1==0.and.nb2==0.and.nb3==0) cycle
              hk2 = nb1*nb1/zl(1)**2+nb2*nb2/zl(2)**2+nb3*nb3/zl(3)**2
              ex = exp((-hk2*fact))/hk2
              ex1 = ex*fact1
              ex2 = ex*fact2
              tabv3(nb1,nb2,nb3) = ex1
              tabf3(:,nb1,nb2,nb3) = ex2*q(:)
           end do
        end do
     end do

  endif



#else
    ! On recalcule et réalloue les cellules, puis on applique les conditions aux
  ! limites périodiques sur les positions des atomes

    CALL ScaleBox(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
#endif


    ! Calcul des forces et des contraintes à l'instant t+dt
    CALL CalFo (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

    ! Calcul de la viscosité à l'instant ...
    DO i=1, nHoover
       zNew(i) = zOld(i) + 2.d0*zDot(i)*tstep    ! ... t+dt
       zOld(i) = zHoover(i)                      ! ... t
       zHoover(i) = zNew(i)                      ! ... t+dt
    END DO

    ! Estimation de la contrainte totale à l'instant t+dt
    ! (la contrainte cinétique est calculée à l'instant t 
    ! et la contrainte potentielle à l'instant t+dt)
    sigtot = 0.5d0*(sigkine + Transpose(sigkine) + sig + Transpose(sig) )

    ! Estimation de la dérivée du tenseur h à l'instant t+dt
    IF (lpcon2.EQV..true.) THEN    ! On ajoute une force de friction
       hdot_new(:,:) = 1.d0/( 1.d0 + 0.5d0*tstep*zHoover(1) + 0.5d0/tbox )*( hdot(:,:)*ihbox0(:,:) &
            + tstep/(2.d0*wBox)*volu*MatMul( sigtot(:,:) - sigext(:,:), invtrh(:,:) ) )*ihbox0(:,:)
    ELSE        ! Équation sans force de friction supplémentaire
       hdot_new(:,:) = 1.d0/( 1.d0 + 0.5d0*tstep*zHoover(1) )*( hdot(:,:)*ihbox0(:,:) &
            + tstep/(2.d0*wBox)*volu*MatMul( sigtot(:,:) - sigext(:,:), invtrh(:,:) ) )*ihbox0(:,:)
    END IF

    ! Estimation de la dérivée du tenseur Gmat à l'instant t+dt
    DO i=1, 3
       DO j=1, 3
          Gdot(i,j) = Sum( h(1:3,i)*hdot_new(1:3,j) + hdot_new(1:3,i)*h(1:3,j) )
       END DO
    END DO

    ! Cycle autocohérent (hdot -> sdot -> sigkine -> hdot)
    DO iter=1, Max_Iter

       ! Valeurs de la dernière itération du cycle d'autocohérence
       hdot_last(:,:) = hdot_new(:,:)*ihbox0(:,:)

       ! Dérivée des coordonnées réduites des atomes à l'instant t+dt
       mf(:,:) = 0.5d0*tstep*MatMul(invGmat,Gdot)
       DO i=1, 3
          mf(i,i) = 1.d0 + 0.5d0*tstep*zHoover(1) + mf(i,i)
       END DO
       CALL MatInv(mf, mfi)
       DO ia=1, im
          sdot_new(:,ia) = MatMul( mfi(:,:), sdot(:,ia)  &
               + tstep/(2.d0*cm(ityp(ia))) * MatMul( invh(:,:), fp(:,ia) ) )
       END DO

       ! Vitesse des atomes à l'instant t+dt
       vp(:,1:im) = MatMul( h(:,:), sdot_new(:,1:im) )

       !  Contrainte thermique à l'instant t+dt
       sigkine(:,:)=0.d0
       do ia = 1, im
          do j = 1,3
             sigkine(1:3,j) = sigkine(1:3,j) + cm(ityp(ia))*vp(1:3,ia)*vp(j,ia)
          enddo
       enddo
       sigkine(1:3,1:3) = invVolu*sigkine(1:3,1:3)
#if(PARA)
       call MPI_ALLREDUCE(sigkine,sigkine_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
       sigkine=sigkine_tot
#endif
       ! Contrainte totale à l'instant t+dt
       sigtot = 0.5d0*(sigkine + Transpose(sigkine) + sig + Transpose(sig) )

       ! Dérivée du tenseur h à l'instant t+dt
       IF (lpcon2.EQV..true.) THEN    ! On ajoute une force de friction
          hdot_new(:,:) = 1.d0/( 1.d0 + 0.5d0*tstep*zHoover(1) + 0.5d0/tbox )*( hdot(:,:)*ihbox0(:,:) &
               + tstep/(2.d0*wBox)*volu*MatMul( sigtot(:,:) - sigext(:,:), invtrh(:,:) ) )*ihbox0(:,:)
       ELSE        ! Équation sans force de friction supplémentaire
          hdot_new(:,:) = 1.d0/( 1.d0 + 0.5d0*tstep*zHoover(1) )*( hdot(:,:) &
               + tstep/(2.d0*wBox)*volu*MatMul( sigtot(:,:) - sigext(:,:), invtrh(:,:) ) )*ihbox0(:,:)
       END IF

       ! Dérivée du tenseur Gmat à l'instant t+dt
       DO i=1, 3
          DO j=1, 3
             Gdot(i,j) = Sum( h(1:3,i)*hdot_new(1:3,j) + hdot_new(1:3,i)*h(1:3,j) )
          END DO
       END DO

       ! Vérifie l'autocohérence de h
       diff = Sum( abs( hdot_new(1:3,1:3) - hdot_last(1:3,1:3) ) )
       tdiff = Sum( abs( hdot_last(1:3,1:3) ) )
       if (tdiff .eq. 0.0d0) then
          if (diff .LE. tol) exit
       else
          if (diff/tdiff .LE. tol) exit
       endif

    END DO

    ! Valeurs convergées de hdot et sdot à l'instant t+dt
    hdot(:,:) = hdot_new(:,:)*ihbox0(:,:)
    sdot(:,:) = sdot_new(:,:)

    IF (iter.GE.Max_Iter) THEN
       WRITE(0,'(a,i0,a)') 'Maximal number of iterations (', Max_Iter, &
            ') in Parrinello-Rahman / Nosé-Hoover self consistency loop has been reached'
       WRITE(0,*)
       WRITE(0,'(a)') 'Last hdot proposed:'
       WRITE(0,'(a,3(f0.5,1x))') ' hdot(1,1:3) = ', 1e8*hdot_last(1,1:3)
       WRITE(0,'(a,3(f0.5,1x))') ' hdot(2,1:3) = ', 1e8*hdot_last(2,1:3)
       WRITE(0,'(a,3(f0.5,1x))') ' hdot(3,1:3) = ', 1e8*hdot_last(3,1:3)
       WRITE(0,*)
       WRITE(0,'(a)') 'New hdot proposed:'
       WRITE(0,'(a,3(f0.5,1x))') ' hdot(1,1:3) = ', 1e8*hdot_new(1,1:3)
       WRITE(0,'(a,3(f0.5,1x))') ' hdot(2,1:3) = ', 1e8*hdot_new(2,1:3)
       WRITE(0,'(a,3(f0.5,1x))') ' hdot(3,1:3) = ', 1e8*hdot_new(3,1:3)
       STOP '< PRNose >'
    END IF

    ! Énergie cinétique des atomes à l'instant t+dt
    kine = 0.5d0*volu*( sigKine(1,1) + sigKine(2,2) + sigKine(3,3) )

    ! Déformation (Eq. 2.16, Ref.2)
    epsi=0.5d0*MatMul( MatMul( invtrh0, Gmat ), invh0 )
    DO i=1, 3
       epsi(i,i) = epsi(i,i) - 1.d0
    END DO

    ! Tension thermodynamique (Eq. 2.22 et 2.26, Ref.2)
    grsig = volu * MatMul(invh, MatMul( sigext, invtrh) )
    tension = invVolu0*MatMul( MatMul( h0, grsig), trh0 )

    ! Énergie potentielle de la cellule (Eq. 2.25, Ref.2)
    Ucell = volu0*Sum( tension(1:3,1:3) * epsi(1:3,1:3) ) 

    ! Énergie cinétique de la cellule (Eq. 2.14, Ref.2)
    Kcell = 0.5d0*wbox*Sum( hDot(1:3,1:3)**2 )
    EcellPR = Kcell + Ucell

    IF (lTHoover) THEN
       ! Dérivée de la viscosité et énergie cinétique du thermostat
       zDot(1) = (2.d0*(kine + Kcell) - gNose*bk*Text)/wHoover(1) &
            - zHoover(2)*zHoover(1)
       KHoover(1) = 0.5d0*wHoover(1)*zHoover(1)**2         ! t+dt
       DO i=2, nHoover
          zDot(i) = ( 2.d0*KHoover(i-1) - bk*Text )/wHoover(i) &
               - zHoover(i+1)*zHoover(i)
          KHoover(i) = 0.5d0*wHoover(i)*zHoover(i)**2
       END DO
       KNose = Sum( KHoover(1:nHoover) )
       ! Énergie potentielle du thermostat 
       !DO i=1, nHoover
       !   UHoover_new(i) = UHoover_old(i) + 2.d0*tstep*bk*Text*zHoover(i)   ! t+dt
       !   UHoover_old(i) = UHoover(i)              ! t
       !   UHoover(i) = UHoover_new(i)              ! t+dt
       !END DO
       !UHoover(1) = gNose*UHoover(1)               ! t+dt
       !UNose = Sum( UHoover(1:nHoover) )           ! t+dt
       !ENose = KNose + UNose                       ! t+dt
    END IF

  end subroutine pr

end module Parrinello_Rahman





!******************
