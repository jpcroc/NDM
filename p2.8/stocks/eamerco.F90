module eam
  USE T_kind_param_m

  implicit none

  !Eamtype, Reptype et DensityType definissent les éléments dont sont censés dépendre les fonctions eam, répulsion et densité. Dans un cas d'alliage on peut avoir la meme forme analytique mais des valeurs différentes des coefficients pour les différents types. La liste des coefficients est definie dans les types. la valeur des coefficients pour ces dfférents types sont définis dans les routines generRho, generEAm, generRep. Ils sont ensuite utilisé dans les routines extrapolate

  type :: EamT
     integer :: toto
  end type EamT

  type :: Rept
     integer :: toto
  end type Rept


  type :: DensityT
     integer::toto
  end type DensityT

  type(DensityT),dimension(:), pointer :: rhotyp
  type(EamT),dimension(:), pointer :: embtyp
  type(repT),dimension(:), pointer :: reppair


  public ::  extrapolateRho, extrapolateRep, extrapolateEam,inputeam

  real(double) :: ev2erg=1.602d-12, & !conversion eV ->erg
       evA2dyn=1.602d-4 ,& ! conversion ev/A -> dyn
       A2cm =1.0d-8     !conversion A->cm
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
    integer:: i,l
    integer :: lupotin=95
    character ::  fnampotin*80


    fnampotin = 'alerco.potin'

    lupotin = 95
    open(unit=lupotin, file=fnampotin)

    iewald=0; l3c=.false.; r3cm=0.

    read(lupotin,*)ntyp
    npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
    call  alloc_typ
    allocate(rhotyp(ntyp)) 
    allocate(embtyp(ntyp)) 
    allocate(reppair(npair)) 

    read(lupotin,*) rue
    rue=rue*1.0d-8
    rumax=max(rue,rumax)
    if(rang==0)  write(6,*) 'Types d_atomes :'
    do i = 1, ntyp
       read (lupotin,*) cm(i),catom(i),ty(i)
       if (rang/=0) cycle
       write (6, '(I4,2F9.3,A5)') i, cm(i),catom(i),ty(i)
    end do
    do i = 1, npair
       read (lupotin,*) roff1(i),roff2(i)
       if (rang/=0) cycle
       write (6, '(2F9.3)') roff1(i),roff2(i)
    end do
    roff1=roff1*1.0d-8
    roff2=roff2*1.0d-8

    cm(:ntyp) = cm(:ntyp)*umass


    close(lupotin)

    return

  end subroutine inputeam

  !-------------------------------------------

  subroutine extrapolateRho(density, r2, rho, drho, ddrho)
    ! calculate electronic density at distance sqrt(r)
    ! or its first and second derivatives

    implicit none

    type(DensityT), intent(in) :: density
    real(kind(0.d0)), intent(in) :: r2
    real(kind(0.d0)), intent(out), optional :: rho, drho, ddrho

    !  Ercolessi potential
    REAL(kind(0.d0)) :: r, func, dfunc, d2func
    r=sqrt(r2)/A2cm

    CALL rh(r,func,dfunc,d2func)
    IF (present(rho)) rho=func
    IF (present(drho)) drho=dfunc/A2cm
    IF (present(ddrho)) ddrho=d2func/(A2cm**2)

    RETURN


  end subroutine extrapolateRho


  !---------------------------------------------------------------------------


  subroutine extrapolateEam(eam, rho, embF, dembF, ddembF, err)
    ! calculate eam function for electronic density rho
    !   or its first and second derivatives
    !   err= 0 if everyting ok
    !       -1 if density too large for extrapolation
    implicit none

    type(EamT), intent(in) :: eam

    real(double), intent(in) :: rho
    real(double), intent(out), optional :: embF, dembF, ddembF
    integer, intent(out), optional :: err



    ! Debugging with Ercolessi potential
    REAL(kind(0.d0)) :: func, dfunc, d2func
    CALL uu(rho, func, dfunc, d2func)
    IF (present(embF)) embF=func*ev2erg
    IF (present(dembF)) dembF=dfunc*ev2erg
    IF (present(ddembF)) ddembF=d2func*ev2erg
    IF (present(err)) err=0
    RETURN

  end subroutine extrapolateEam

  !----------------------------------------------

  subroutine extrapolateRep(rep, r2, Erep, dErep, ddErep)
    ! calculate repulsive potential at distance sqrt(r2)
    ! or its first and second derivatives

    implicit none

    type(RepT), intent(in) :: rep
    real(kind(0.d0)), intent(in) :: r2
    real(kind(0.d0)), intent(out), optional :: Erep, dErep, ddErep


    !-----------------------------------
    ! Debugging with Ercolessi potential
    REAL(kind(0.d0)) :: r, func, dfunc, d2func
    r=sqrt(r2)/A2cm
    CALL v2(r, func, dfunc, d2func)
    IF (present(Erep)) Erep=func*ev2erg
    IF (present(dErep)) dErep=dfunc*ev2erg/A2cm
    IF (present(ddErep)) ddErep=d2func*ev2erg/(A2cm**2)
    RETURN
    !-----------------------------------
  end subroutine extrapolateRep
end module eam
