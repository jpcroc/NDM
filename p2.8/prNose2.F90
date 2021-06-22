module Parrinello_Rahman_Nose

  ! Algorithme de Parrinello-Rahman [1] combiné à un thermostat de Nosé [4]
  ! basé sur les article de Ray et Rahman [2,3]
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
  !
  ! En résumé, l'algorithme est basé sur les équations (3.2), (3.3) et (3.4) de [3]
  !
  ! [1] Parrinello, M. & Rahman,
  !     A. Polymorphic Transitions in Single rcystals: A New Molecular Dynamics Method
  !     J. Appl. Phys., 1981, 52, 7182-7190
  ! [2] Ray, J.R. & Rahman, A.
  !     Statistical Ensembles and Molecular Dynamics Studies of Anisotropic Solids
  !     J. Chem. Phys., 1984, 80, 4423-4428
  ! [3] Ray, J.R. & Rahman, A.
  !     Statistical Ensembles and Molecular Dynamics Studies of Anisotropic Solids. II
  !     J. Chem. Phys., 1985, 82, 4243-4247
  ! [4] Nosé, S.
  !     A Molecular Dynamics Method for Simulations in the Canonical Ensemble
  !     Mol. Phys., 1984, 52, 255-268


  USE T_kind_param_m
  USE gen_com_m, ONLY:   ecellpr,enose,fnose,kcell,kine,knose,lpcon2,sigext,sigtot,tbox,text,&
       &tstep,ucell,unose,wboxf,wnose,enose,erg2ev,fnose,im_glob,it,kcell,knose,leev,&
       &lucell,rang,timel,tstep,unose,wnose,sigkine,rang,sig,bk,lspaceNDM,h0
  USE var_pot, ONLY:cm
  USE tempinstT_mod,only: tempinstT
  USE Mat_utils_mod,only:  matinv
  USE recips_mod,only: recips,calcvol
  USE boxconfig,only:box_config
  use atomconfig,only:atom_config_d
  use cellconfig,only:cell_config
#ifdef PARA
  USE Tpara,only:COMM_space,nprocspace
#else
  use Tpara,only:nprocspace

#endif

  implicit none

  real(double), dimension(3,3), save , private ::h,trh,invh,invtrh,Gmat,invGmat,Area,hnew,hlast,hold,invhold
  real(double), dimension(3,3), save , private ::hpoint, h2point, whpointpoint, Gpoint
  real(double), allocatable, save, private :: sp(:,:),sold(:,:),snew(:,:),sdot(:,:)

  ! Variables uniquement nécessaires au calcul de l'énergie potentielle de la
  ! boîte
  real(double), dimension(3,3), save , private ::trh0,invh0,invtrh0,epsi, tension
  real(double), save , private ::volu0, invVolu0

  ! Variable associée au thermostat de Nosé
  !  (fNose est défini dans gen_com_m.F90)
  REAL(double), save, private :: fnew, flast, fold, fpoint, f2point

  ! Nombre de degrés de liberté
  REAL(double), save, private :: gNose
  real(double)::wbox

contains

  subroutine initlprNose(atpr,celndm,boxndm)
    type(box_config)::boxndm
    class(atom_config_d)::atpr
    type(cell_config):: celndm

    real(double), dimension(1:3,1:3) :: maux2
    real(double), external :: detmat
    real(double)::temp0, unitE
    character*5 :: cunitE


    ! Calcul de la température initiale
    temp0=tempinstT(atpr)

    if (rang==0) WRITE(6,*)
    if (rang==0) WRITE(6,'(a)') 'Algorithme de Parrinello-Rahman couplé au thermostat de Nosé (V2)'
    if (rang==0) WRITE(6,'(a)') '  -> la vitesse de la boîte ne prend pas en compte la dérivée du tenseur h'
    if (rang==0) WRITE(6,'(a)') "  -> l'énergie cinétique des atomes et de la boîte est thermalisée"
    if (rang==0) WRITE(6,*)
    IF (lUcell) THEN
       if (rang==0) WRITE(6,'(a)') "Repère de référence pour Parrinello-Rahman  (A):"
       if (rang==0) WRITE(6,'(a,3(f0.5,1x))') ' h0(1:3,1) = ', 1e8*h0(1:3,1)
       if (rang==0) WRITE(6,'(a,3(f0.5,1x))') ' h0(1:3,2) = ', 1e8*h0(1:3,2)
       if (rang==0) WRITE(6,'(a,3(f0.5,1x))') ' h0(1:3,3) = ', 1e8*h0(1:3,3)
       if (rang==0) WRITE(6,*)
    ELSE
       h0 = boxndm%at
    END IF
    if (rang==0) WRITE(6,'(a)') "Repère actuel  (A):"
    if (rang==0) WRITE(6,'(a,3(f0.5,1x))') ' h (1:3,1) = ', 1e8*boxndm%at(1:3,1)
    if (rang==0) WRITE(6,'(a,3(f0.5,1x))') ' h (1:3,2) = ', 1e8*boxndm%at(1:3,2)
    if (rang==0) WRITE(6,'(a,3(f0.5,1x))') ' h (1:3,3) = ', 1e8*boxndm%at(1:3,3)
       wbox =wboxf*sum(0.5*cm(atpr%ityp(:atpr%im)))       ! La moitié de la masse totale des atomes
#ifdef PARA
if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
  call comm_space%sum(wbox)
end if
#endif
    if (rang==0) WRITE(6,'(a,g20.12)')'Masse de la boîte pour Parrinello-Rahman: wbox=',wbox

    ! Nombre de degrés de liberté pour le thermostat de Nosé
    gNose = dble(3*im_glob+1)

    IF (wNose.EQ.0) THEN
       ! On suppose que la fréquence de vibration typique du solide est
       !  1 THz = 1e-12 s¯¹
       !  => ça ne marche pas
       ! On veut qu'une variation de la température de 10K corresponde à
       ! une variation de f de 1% avec f~1
       wNose = gNose*bk*10.d0*tstep**2/1.d-2**2
    END IF
    if (rang==0) WRITE(6,'(a,g20.12)')'Masse de la boîte pour thermostat de Nosé: wNose=',wNose
    if (rang==0) WRITE(6,'(a,g20.12)')'Nombre de degrés de liberté: gNose=',gNose


    ! État de référence défini par la matrice h0
    !   Cet état de référence doit correspondre à un tenseur de contrainte nul.
    !   Il n'est utile que pour calculer la déformation et l'énergie potentielle
    !   de la boîte.
    volu0 = calcvol(h0(1:3,1),h0(1:3,2),h0(1:3,3))
    invVolu0 = 1.d0/volu0
    trh0=Transpose(h0)
    CALL MatInv(h0,invh0)
    invtrh0=Transpose(invh0)

    ! Vecteurs de la boîte et matrice inverse
    h = boxndm%at

    ! Coordonnées réduites des atomes
    allocate(sp(3,atpr%imm),sdot(3,atpr%imm),sold(3,atpr%imm),snew(3,atpr%imm))

    ! Initialisation de la vitesse de la boîte
    !if (rang==0) WRITE(6,'(a,f0.3,a)') 'Initialisation de la vitesse de la boîte pour la température ', temp0, ' K'
    !do i=1,3
    !do j=1,3
    !call random_number(z1)
    !call random_number(z2)
    !if(z1.eq.0.d0) z1=0.000000001d0
    !hpoint(i,j)=sqrt(2*bk*temp0/wbox)*sqrt(-log(z1))*(1.-2*z2)
    !end do
    !end do
    if (rang==0) WRITE(6,'(a,f0.3,a)') 'Initialisation de la vitesse de la boîte pour la température ', 0.d0, ' K'
    hpoint(:,:) = 0.d0
    hold(:,:) = h(:,:) - tstep*hpoint(:,:)
    ! Kinetic energy of the cell (Eq. 2.14 of Ref. [2])
    maux2 = MatMul( Transpose(hpoint), hpoint )
    Kcell = 0.5d0*wbox*( maux2(1,1) + maux2(2,2) + maux2(3,3) )
    if(lEev) then
       unitE=erg2eV
       cunitE='  eV'
    else
       unitE=1.0
       cunitE=' erg'
    end if
    write(6,'(I7,D10.3,A,D21.12,A,a,f0.3,a)') 0,0.d0,'*Kcell = ',Kcell*unitE,cunitE, &
         '  (', 2.d0*Kcell/(9.d0*bk), ' K)'

    ! Kinetic and potential energies of Nosé thermostat (Eq. 3.1 Ref. [3])
    !KNose = 0.5d0*bk*temp0
    KNose = 0.d0
    fpoint = Sqrt(2.d0*KNose/wNose)
    UNose = 0.d0
    ENose = KNose + UNose
    fNose=1.d0
    fold = fNose - fpoint*tstep
    if (rang==0) WRITE(6,'(3(a,g22.12))') 'fNose = ', fNose, '  fold = ', fold, '  fpoint = ', fpoint
    if (rang==0) WRITE(6,'(a,f0.3,a)') 'Initialisation du thermostat de Nosé pour la température ', &
         0.d0, ' K'
    !if (rang==0) WRITE(6,'(a,f0.3,a)') 'Initialisation du thermostat de Nosé pour la température ', &
    !temp0, ' K'
    if (rang==0) WRITE(6,'(I7,D10.3,A,D21.12,A,a,f0.3,a)') it,timel,'*KNose = ',KNose*unitE,cunitE, &
         '  (', 2.d0*KNose/bk, ' K)'
    if (rang==0) WRITE(6,'(3(a,g22.12))') 'fNose = ', fNose, ' -  fold = ', fold, &
         ' -  fpoint * tstep = ', fpoint*tstep
    if (rang==0) WRITE(6,*)

    return
  end subroutine initlprNose

  !-----------------------------------------------

  subroutine prNose(atpr,celndm,boxndm)

    implicit none
    type(box_config)::boxndm
    class(atom_config_d)::atpr
    type(cell_config):: celndm


    real(double),dimension(3,3)::maux1,maux2,mf,mfi, grsig
    REAL(double), dimension(1:3,1:3) ::  Gpoint
    real(double):: diff,tdiff, invVolu, fNose2, f2point
    integer:: i,j,ia, iter
    !real(double) , external ::  calcvol

    ! Parameter for Parrinello-Rahman self consistency loop
    REAL(double), parameter :: tol=1.0d-12        ! Tolerance for h convergency
    INTEGER, parameter :: max_Iter=100            ! Maximal number of iterations in self-consistency loop


    ! Paramètres du thermostat
    fNose2=fNose*fNose

    ! Vecteurs de la boîte
    h(:,:)=boxndm%at(:,:)
    trh=Transpose(h)
    ! Métrique de la boîte
    Gmat = MatMul(trh,h)
    call MatInv(Gmat,invGmat)
    call MatInv(h,invh)
    call MatInv(hold,invhold)
    invtrh = Transpose(invh)
    boxndm%volu = calcvol(h(1:3,1),h(1:3,2),h(1:3,3))
    invVolu = 1.d0/boxndm%volu
    Area(:,:)=boxndm%volu*invtrh(:,:)    ! Correspond à sigma dans l'article de Parrinello Rahman

    ! Déduit de la tension thermodynamique correspondant à la contrainte imposée
    ! la matrice grsig
    ! éq. 2.22 et 2.26 dans l'article de Ray et Rahman
    grsig = boxndm%volu * MatMul(invh, MatMul( sigext, invtrh) )

    ! Strain tensor (Eq. 2.16)
    epsi=0.5d0*MatMul( MatMul( invtrh0, Gmat ), invh0 )
    DO i=1, 3
       epsi(i,i) = epsi(i,i) - 1.d0
    END DO

    ! Thermodynamic tension (Eq. 2.22)
    tension = invVolu0*MatMul( MatMul( h0, grsig), trh0 )

    ! Potential energy of the cell (Eq. 2.25)
    maux1 = MatMul( tension, epsi )
    Ucell = volu0*( maux1(1,1) + maux1(2,2) + maux1(3,3) )

    ! Coordonnées réduites des atomes
    sp(1:3,1:atpr%imm) = MatMul(invh(1:3,1:3), atpr%xp(1:3,1:atpr%imm) )
    sold(1:3,1:atpr%imm) = MatMul( invhold(1:3,1:3), atpr%xpp(1:3,1:atpr%imm) )
    ! snew est le propagé de s avec seulement fp unc==uncorrected
    do ia = 1,atpr%im
       snew(1:3,ia) = -sold(1:3,ia) + 2.d0*sp(1:3,ia) &
            + tstep**2/(fNose2*cm(atpr%ityp(ia)))*MatMul( invH(1:3,1:3), atpr%fp(1:3,ia))
    enddo

    ! Compute initial guess for Parrinello-Rahman and Nosé
    !   variables at next time step ...........................................
    ! ce hnew est le premier h(in) de la boucle autocohérente (à noter pas de force sur h)
    hnew = 2.d0*h - hold
    fnew = 2.d0*fNose - fold
    iter = 0

    ! Start selfconsistency loop to calculate P-R and Nosé variables ..........
10  continue
    iter = iter + 1

!!$    if (rang==0) WRITE(6,*)                                                          ! DEBUG
!!$    if (rang==0) WRITE(6,'(a)') 'Parrinello-Rahman / Nosé self consistent loop'      ! DEBUG
!!$    if (rang==0) WRITE(6,'(a,i0)') '  iter = ', iter                                 ! DEBUG
!!$    if (rang==0) WRITE(6,'(a,3(f0.5,1x))') ' h(1:3,1) = ', 1e8*hnew(1:3,1)           ! DEBUG
!!$    if (rang==0) WRITE(6,'(a,3(f0.5,1x))') ' h(1:3,2) = ', 1e8*hnew(1:3,2)           ! DEBUG
!!$    if (rang==0) WRITE(6,'(a,3(f0.5,1x))') ' h(1:3,3) = ', 1e8*hnew(1:3,3)           ! DEBUG
!!$    if (rang==0) WRITE(6,*)                                                          ! DEBUG
!!$    if (rang==0) WRITE(6,'(4(a,g22.12))') 'fnew = ', fnew, '  fNose = ', fNose, &
!!$        '  fold = ', fold, '  fpoint = ', fpoint                        ! DEBUG

    IF (iter.GT.Max_Iter) THEN
       WRITE(0,'(a,i0,a)') 'Maximal number of iterations (', Max_Iter, &
            ') in Parrinello-Rahman / Nosé self consistency loop has been reached'
       WRITE(0,*)
       WRITE(0,'(a)') 'Last h proposed:'
       WRITE(0,'(a,3(f0.5,1x))') ' h(1:3,1) = ', 1e8*hlast(1:3,1)
       WRITE(0,'(a,3(f0.5,1x))') ' h(1:3,2) = ', 1e8*hlast(1:3,2)
       WRITE(0,'(a,3(f0.5,1x))') ' h(1:3,3) = ', 1e8*hlast(1:3,3)
       WRITE(0,*)
       WRITE(0,'(a)') 'New h proposed:'
       WRITE(0,'(a,3(f0.5,1x))') ' h(1:3,1) = ', 1e8*hnew(1:3,1)
       WRITE(0,'(a,3(f0.5,1x))') ' h(1:3,2) = ', 1e8*hnew(1:3,2)
       WRITE(0,'(a,3(f0.5,1x))') ' h(1:3,3) = ', 1e8*hnew(1:3,3)
       WRITE(0,*)
       WRITE(0,'(a,f0.5)') 'Last fNose proposed: ', flast
       WRITE(0,'(a,f0.5)') 'New fNose proposed:  ', fnew
       STOP '< PRNose >'
    END IF

    ! Valeurs de la dernière itération du cycle d'autocohérence
    hlast = hnew
    flast = fnew

    ! hpoint guess for time derivatives at current time, and related stuff
    hpoint = (hnew - hold)/(2.d0*tstep)
    ! Kinetic energy of the cell (Eq. 2.14, Ref.[2])
    maux2 = MatMul( Transpose(hpoint), hpoint )
    Kcell = 0.5d0*wbox*fNose2*( maux2(1,1) + maux2(2,2) + maux2(3,3) )
    ! Time derivative of the metric tensor Gmat
    DO i=1, 3
       DO j=1, 3
          Gpoint(i,j) = Sum( h(1:3,i)*hpoint(1:3,j) + hpoint(1:3,i)*h(1:3,j) )
       END DO
    END DO
    fpoint = (fnew - fold)/(2.d0*tstep)

    ! Résolution de l'équation 3.2 de la Ref. [3]
    mf = 0.5d0*tstep*MatMul(invGmat,Gpoint)
    do i = 1,3
       mf(i,i) = mf(i,i) + 1.d0 + tstep*fpoint/fNose
    enddo
    call matinv(mf,mfi)
    mfi = 0.5d0/tstep*mfi
    sdot(1:3,1:atpr%im) = MatMul(mfi(1:3,1:3), snew(1:3,1:atpr%im) - sold(1:3,1:atpr%im) )

    ! avec ce sdot on peut calculer la vitesse des particules
    atpr%vp(1:3,1:atpr%imm) = fNose*MatMul(h(1:3,1:3),sdot(1:3,1:atpr%imm))
    !  ... la contrainte thermique associée
    sigkine(:,:)=0.d0
    do ia = 1, atpr%im
       do j = 1,3
          sigkine(1:3,j) = sigkine(1:3,j) + cm(atpr%ityp(ia))*atpr%vp(1:3,ia)*atpr%vp(j,ia)
       enddo
    enddo
    sigkine(1:3,1:3) = invVolu*sigkine(1:3,1:3)
#ifdef PARA
if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
   call comm_space%sum(sigkine)
 end if

#endif

    ! ... et l'énergie cinétique
    kine = 0.5d0*boxndm%volu*( sigKine(1,1) + sigKine(2,2) + sigKine(3,3) )

    ! Contrainte totale
    sigtot = 0.5d0*(sigkine + Transpose(sigkine) + sig + Transpose(sig) )

    ! Résolution de l'équation (3.3) de la Ref. [3]
    if(lpcon2.EQV..true.)then
       whpointpoint(:,:) = MatMul(sigtot,Area)- MatMul(h, grsig) &
            - 2.d0*wbox*fNose*fpoint*hpoint(:,:) &
            - wbox/(tstep*tbox)*hpoint(:,:)
    else
       whpointpoint(:,:) = MatMul(sigtot,Area) - MatMul(h, grsig)&
            - 2.d0*wbox*fNose*fpoint*hpoint(:,:)
    end if
    hnew(:,:) = 2.d0*h(:,:) - hold(:,:) + whpointpoint(:,:)*tstep**2/(fNose2*wbox)

    ! Résolution de l'équation (3.4) de la Réf. [3]
    f2point = (2.d0*(kine+Kcell) - gNose*bk*Text)/(fNose*wNose)
    fnew = 2.d0*fNose - fold + f2point*tstep**2

    ! Vérifie l'autocohérence
    ! --- de h
    diff = Sum( abs( hnew(1:3,1:3) - hlast(1:3,1:3) ) )
    tdiff = Sum( abs( hlast(1:3,1:3) ) )
    if (tdiff .eq. 0.0d0) then
       if (diff .gt. tol) goto 10
    else
       if (diff/tdiff .gt. tol) goto 10
    endif
    ! --- de fNose
    diff = abs( fnew - flast )
    tdiff = abs( flast )
    if (tdiff .eq. 0.0d0) then
       if (diff .gt. tol) goto 10
    else
       if (diff/tdiff .gt. tol) goto 10
    endif

    ! Ici hnew et fnew sont convergés
!!$  if (rang==0) WRITE(6,'(a,i0)') 'PR: iter = ', iter
    snew(1:3,1:atpr%imm) = sold(1:3,1:atpr%imm) + 2.d0*tstep*sdot(1:3,1:atpr%imm)

    ! Save current atomic positions as old ones,
    !   and next positions as current ones
    !on connait hdot et sdot
    sold(1:3,1:atpr%imm) = sp(1:3,1:atpr%imm)
    sp(1:3,1:atpr%imm) = snew(1:3,1:atpr%imm)

    hold = h
    h = hnew
    fold = fNose
    fNose = fnew

    ! Transform back to absolute coordinates
    boxndm%at(:,:) =  h(:,:)
    call recips (boxndm%at(1:3,1), boxndm%at(1:3,2), boxndm%at(1:3,3), boxndm%bg(1:3,1), boxndm%bg(1:3,2), boxndm%bg(1:3,3))
    boxndm%volu = calcvol(boxndm%at(1:3,1),boxndm%at(1:3,2),boxndm%at(1:3,3))

    atpr%xp(1:3,1:atpr%imm) = MatMul(h(1:3,1:3), sp(1:3,1:atpr%imm) )
    atpr%xpp(1:3,1:atpr%imm) = MatMul(hold(1:3,1:3), sold(1:3,1:atpr%imm) )

    ! Total energy of the cell
    EcellPR = Kcell + Ucell

    ! Kinetic and potential energies of Nosé thermostat (Eq. 3.1 Ref. [3])
    KNose = 0.5d0*wNose*fpoint**2
    UNose = gNose*bk*Text*log(fNose)
    ENose = KNose + UNose


  end subroutine prNose

end module ! Parrinello_Rahman_Nose





!******************
