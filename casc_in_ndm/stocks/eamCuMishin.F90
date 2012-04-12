module eam
  USE T_kind_param_m

  !=============================================================================!
  !  Eam potential parameters corresponding to Cu potential                     !
  !  Y. Mishin, M.J. Mehl, D.A. Papaconstantopoulos, A.F. Voter, and J.D. Kress !
  !  "Structural Stability and Lattice Defects in Copper:                       !
  !    Ab Initio, Tight-Binding and Embedded-Atom Calculations",                !
  !    Phys. Rev. B 63 (2001) 224106                                            !
  !=============================================================================!

  implicit none


  type :: EamT
     real(double),dimension (:),pointer :: feam,dfeam,xg
     real(double)::deltaEAM
  end type EamT

  type :: RepT
     real(double),dimension (:),pointer :: potr,dpotr,xr
     real(double)::deltaREP
  end type Rept

  type :: DensityT
     real(double),dimension (:),pointer::rho,drho,xd
     real(double)::deltaRHO
  end type DensityT

  type(DensityT),dimension(:), pointer :: rhotyp
  type(EamT),dimension(:), pointer :: embtyp
  type(repT),dimension(:), pointer :: reppair


  integer :: nptmax ! nombre de points dans la grille lue

  public ::  extrapolateRho, extrapolateRep, extrapolateEam,inputeam


  PRIVATE :: RepMishin, RhoMishin, EmbMishin

  REAL(double), parameter, private :: RcutMishin=5.50679d0        ! A
  REAL(double), parameter, private :: h=0.50037d0           ! A
  REAL(double), parameter, private :: E1=2.01458d2          ! eV
  REAL(double), parameter, private :: E2=6.59288d-3         ! eV
  REAL(double), parameter, private :: r01=0.83591d0         ! A
  REAL(double), parameter, private :: r02=4.46867d0         ! A
  REAL(double), parameter, private :: alpha1=2.97758d0      ! A¯¹
  REAL(double), parameter, private :: alpha2=1.54927d0      ! A¯¹
  REAL(double), parameter, private :: delta=0.86225d-2      ! A
  REAL(double), parameter, private :: Rs1=2.24d0            ! A
  REAL(double), parameter, private :: Rs2=1.8d0             ! A
  REAL(double), parameter, private :: Rs3=1.2d0             ! A
  REAL(double), parameter, private :: S1=4.d0               ! eV/A^4
  REAL(double), parameter, private :: S2=40.d0              ! eV/A^4
  REAL(double), parameter, private :: S3=1.15d3             ! eV/A^4
  REAL(double), parameter, private :: a=3.80362d0
  REAL(double), parameter, private :: r03=-2.19885d0        ! A
  REAL(double), parameter, private :: r04=-2.61984d2        ! A
  REAL(double), parameter, private :: beta1=0.17394d0       ! A¯²
  REAL(double), parameter, private :: beta2=5.35661d2       ! A¯²
  REAL(double), parameter, private :: F0=-2.28235d0         ! eV
  REAL(double), parameter, private :: F2=1.35535d0          ! eV
  REAL(double), parameter, private :: q1=-1.27775d0         ! eV
  REAL(double), parameter, private :: q2=-0.86074d0         ! eV
  REAL(double), parameter, private :: q3=1.78804d0          ! eV
  REAL(double), parameter, private :: q4=2.97571d0          ! eV
  REAL(double), parameter, private :: bigQ1=0.4d0
  REAL(double), parameter, private :: bigQ2=0.3d0
  REAL(double), parameter, private :: rhomax=2.5d0       ! Densité maximale pour la fonction glue

contains


  !---------------------------------------------------------------------------
  subroutine inputeam(ntyp,npair,ntrip,cm,catom,ty,umass,rue,rumax,iewald,l3c,rang,r3cm,roff1,roff2)

    !

    integer, intent(out) :: ntyp                           !nb de type
    integer, intent(out) :: npair                  ! = ntyp*(ntyp+1)/2
    integer, intent(out) :: ntrip                  ! = ntyp*ntyp *(ntyp+1)/2
    real(double), dimension(:), pointer :: cm, catom,roff1,roff2
    character , dimension(:), pointer  :: ty*3
    real(double), intent(in) ::umass 
    real(double), intent(out) :: rue,rumax,r3cm
    integer, intent(out) :: iewald
    logical, intent(out) :: l3c
    integer , intent(in) ::rang



    !local variables
    integer:: i,l,k,n,npt


    iewald=0; l3c=.false.; r3cm=0.

    ntyp=1        ! nombre de types
    npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
    call  alloc_typ

    rue=RcutMishin ! rayon de coupure en Ang
    rue=rue*1.0d-8
    rumax=max(rue,rumax)

    cm(1)=63.54d0 ! masse
    cm(:ntyp) = cm(:ntyp)*umass
    catom(1)=29   ! numéro atomique
    ty(1)='Cu'    ! symbole

    roff1(1)=-0.1 ; roff2(1)=-0.2   ! raccordement Ziegler (en Ang)
    roff1=roff1*1.0d-8
    roff2=roff2*1.0d-8
    
    nptmax=500000 ! nombre de points dans les grilles
    write(6,*)'nptmax',nptmax
    allocate(rhotyp(ntyp)) 
    allocate(embtyp(ntyp)) 
    allocate(reppair(npair)) 

    !lecture de Glue
    npt=nptmax
    embtyp(1)%deltaEAM=rhomax/dble(npt-1)      ! Incrément de densité entre chaque point
    allocate(embtyp(1)%feam(nptmax)) 
    allocate(embtyp(1)%dfeam(nptmax)) 
    allocate(embtyp(1)%xg(nptmax)) 
    do i=1,nptmax
       if(i.le.npt) then
          ! Densité où on calcule la fonction glue
          embtyp(1)%xg(i) = dble(i-1)*embtyp(1)%deltaEAM
          ! Fonction glue
          embtyp(1)%feam(i) = EmbMishin(embtyp(1)%xg(i))
       else
          embtyp(1)%xg(i)=(i-npt)*embtyp(1)%deltaEAM+ embtyp(1)%xg(i)
          embtyp(1)%feam(i)=embtyp(1)%feam(npt)
       end if

    end do
    do i=1,nptmax-1
       embtyp(1)%dfeam(i)=embtyp(1)%feam(i+1)-embtyp(1)%feam(i)
    end do
    embtyp(1)%dfeam(nptmax)=0.

    !lecture de dens
    npt=nptmax
    rhotyp(1)%deltaRHO = RcutMishin/dble(npt-1)        ! Incrément de distance entre chaque point
    allocate(rhotyp(1)%rho(nptmax)) 
    allocate(rhotyp(1)%drho(nptmax)) 
    allocate(rhotyp(1)%xd(nptmax)) 
    do i=1,nptmax
       if(i.le.npt) then
          ! Distance où on calcule la fonction densité
          rhotyp(1)%xd(i) = dble(i-1)*rhotyp(1)%deltaRHO
          ! Fonction densité
          rhotyp(1)%rho(i) = RhoMishin(rhotyp(1)%xd(i))
       else
          rhotyp(1)%xd(i)=(i-npt)*rhotyp(1)%deltaRHO+ rhotyp(1)%xd(npt)
          rhotyp(1)%rho(i)=rhotyp(1)%rho(npt)
          rhotyp(1)%drho(i)=0.
       end if
    end do
    do i=1,nptmax-1
       rhotyp(1)%drho(i)=rhotyp(1)%rho(i+1)-rhotyp(1)%rho(i)
    end do
    rhotyp(1)%drho(nptmax)=0.

    !lecture de la fonction de paire
    npt=nptmax
    reppair(1)%deltaREP = RcutMishin/dble(npt-1)        ! Incrément de distance entre chaque point
    allocate(reppair(1)%potr(nptmax)) 
    allocate(reppair(1)%dpotr(nptmax)) 
    allocate(reppair(1)%xr(nptmax)) 
    do i=1,nptmax
       if(i.le.npt) then
          ! Distance où on calcule l'énergie de paire
          reppair(1)%xr(i) = dble(i-1)*reppair(1)%deltaREP
          ! Énergie de paire
          reppair(1)%potr(i) = RepMishin(reppair(1)%xr(i))
       else
          reppair(1)%xr(i)=(i-npt)*reppair(1)%deltaREP+ reppair(1)%xr(i)
          reppair(1)%potr(i)=reppair(1)%potr(npt)
       end if
    end do
    do i=1,nptmax-1
       reppair(1)%dpotr(i)=reppair(1)%potr(i+1)-reppair(1)%potr(i)
    end do
    reppair(1)%dpotr(nptmax)=0.

    return

  end subroutine inputeam

  !-------------------------------------------

  subroutine extrapolateRho(density, r2, rho)
    ! calculate electronic density at distance sqrt(r)

    implicit none

    type(DensityT), intent(in) :: density 
    real(kind(0.d0)), intent(in) :: r2
    real(kind(0.d0)), intent(out) :: rho
    !local 
    integer :: iti,n,k
    real(double) :: xmax,xmin,ktor,r,drk

!!$    xmax=density%xd(nptmax)
!!$    r=sqrt(r2)/A2cm
!!$
!!$    if(r.gt.xmax) then
!!$       Rho=0.0
!!$      else
!!$       k=Int(r/density%deltaRHO)+1
!!$       drk=r/density%deltaRHO+1-k
!!$       Rho=density%rho(k)+drk*density%drho(k)
!!$    end if

    r=sqrt(r2)*1.d8
    rho = RhoMishin(r)
    RETURN


  end subroutine extrapolateRho


  !---------------------------------------------------------------------------


  subroutine extrapolateEam(eam, rho, embF, err)
    ! calculate eam function for electronic density rho
    !   err= 0 if everyting ok
    !       -1 if density too large for extrapolation
    USE gen_com_m
    implicit none

    type(EamT), intent(in) :: eam

    real(double), intent(in) :: rho
    real(double), intent(out), optional :: embF
    integer, intent(out), optional :: err

    !local 
    integer :: iti,n,k
    real(double) :: xmax,xmin,ktor,drk

!!$    xmax=eam%xg(nptmax)
!!$    
!!$    if (rho.lt.0.) then 
!!$       embf=0
!!$       return
!!$    end if
!!$    if(rho.gt.xmax) then
!!$       embF=eam%feam(nptmax)*ev2erg
!!$
!!$    else
!!$       k=Int(rho/eam%deltaEAM)+1
!!$       drk=rho/eam%deltaEAM +1 -k
!!$       embf=  ev2erg*(eam%feam(k)+drk*eam%dfeam(k))
!!$       
!!$    end if


    embF = ev2erg*EmbMishin(rho)
    RETURN

  end subroutine extrapolateEam

  !----------------------------------------------

  subroutine extrapolateRep(rep, r2, Erep)
    ! calculate repulsive potential at distance sqrt(r2)
    USE gen_com_m
    implicit none

    type(RepT), intent(in) :: rep
    real(kind(0.d0)), intent(in) :: r2
    real(kind(0.d0)), intent(out), optional :: Erep


    !-----------------------------------

    integer :: k
    real(double) :: xmax,r,drk

!!$    xmax=rep%xr(nptmax)
!!$    r=sqrt(r2)/A2cm
!!$    if(r.gt.xmax) then
!!$       Erep=0.0
!!$       return
!!$    else
!!$       k=Int(r/rep%deltaREP)+1
!!$       drk=r/rep%deltaREP+1-k
!!$       Erep=ev2erg*(rep%potr(k)+drk*rep%dpotr(k))
!!$    end if

    r=sqrt(r2)*1.d8
    Erep=ev2erg*RepMishin(r)
    RETURN
    !-----------------------------------
  end subroutine extrapolateRep



  !===========================================================================

  FUNCTION RepMishin(r) RESULT(Erep)
    ! Calculate repulsive part of the potential
    IMPLICIT NONE
    REAL(kind(0.d0)), intent(in) :: r
    REAL(kind(0.d0)) :: Erep

    REAL(kind(0.d0)) :: x4

    IF (r.LT.RcutMishin) THEN
       x4=((r-RcutMishin)/h)**4
       Erep = ( E1*( exp(-2.d0*alpha1*(r-r01)) - 2.d0*exp(-alpha1*(r-r01)) ) &
            + E2*( exp(-2.d0*alpha2*(r-r02)) - 2.d0*exp(-alpha2*(r-r02)) ) &
            + delta )*x4/(1.d0+x4)
       IF (r.LE.Rs1) THEN
          Erep = Erep - S1*(Rs1-r)**4
       END IF
       IF (r.LE.Rs2) THEN
          Erep = Erep - S2*(Rs2-r)**4
       END IF
       IF (r.LE.Rs3) THEN
          Erep = Erep - S3*(Rs3-r)**4
       END IF
    ELSE
       Erep = 0.d0
    END IF

  END FUNCTION RepMishin

  !===========================================================================

  FUNCTION RhoMishin(r) RESULT(rho)
    ! Calculate electronic density
    IMPLICIT NONE
    REAL(kind(0.d0)), intent(in) :: r
    REAL(kind(0.d0)) :: rho

    REAL(kind(0.d0)) :: x4

    IF (r.LT.RcutMishin) THEN
       x4=((r-RcutMishin)/h)**4
       rho = ( a*exp(-beta1*(r-r03)**2) + exp(-beta2*(r-r04)) )*x4/(1.d0+x4)
    ELSE
       rho = 0.d0
    END IF

  END FUNCTION RhoMishin

  !===========================================================================

  FUNCTION EmbMishin(rho) RESULT(Eemb)

    IMPLICIT NONE
    REAL(kind(0.d0)), intent(in) :: rho
    REAL(kind(0.d0)) :: Eemb

    IF (rho.LE.1.d0) THEN
       Eemb = F0 + 0.5d0*F2*(rho-1.d0)**2 &
            + q1*(rho-1.d0)**3 &
            + q2*(rho-1.d0)**4 &
            + q3*(rho-1.d0)**5 &
            + q4*(rho-1.d0)**6 
    ELSE
       Eemb = ( F0 + 0.5d0*F2*(rho-1.d0)**2 &
            + q1*(rho-1.d0)**3 + bigQ1*(rho-1.d0)**4 ) &
            / ( 1.d0 + bigQ2*(rho-1.d0)**3 )
    END IF

  END FUNCTION EmbMishin


end module eam
     
