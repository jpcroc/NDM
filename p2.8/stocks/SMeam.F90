module eam
  USE T_kind_param_m

  implicit none

  !Eamtype, Reptype et DensityType definissent les éléments dont sont censés dépendre les fonctions eam, répulsion et densité. Dans un cas d'alliage on peut avoir la meme forme analytique mais des valeurs différentes des coefficients pour les différents types. La liste des coefficients est definie dans les types. la valeur des coefficients pour ces dfférents types sont définis dans les routines generRho, generEAm, generRep. Ils sont ensuite utilisé dans les routines extrapolate

  type :: EamT
     real(double)::ksi
  end type EamT

  type :: RepT
     real(double)::A
     real(double)::r0
     real(double)::p
  end type RepT


  type :: DensityT
     real(double)::r0
     real(double)::q
  end type DensityT

  type(DensityT),dimension(:), pointer :: rhotyp  !(ntyp)
  type(EamT),dimension(:), pointer :: embtyp      !(ntyp)
  type(repT),dimension(:), pointer :: reppair     !npair


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
    integer:: i
    integer :: lupotin=95
    character ::  fnampotin*80


    fnampotin = 'SM.potin'

    lupotin = 95
    open(unit=lupotin, file=fnampotin, status='old')

    iewald=0; l3c=.false.; r3cm=0.

    read(lupotin,*)ntyp
    npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
    call  alloc_typ
    allocate(rhotyp(ntyp)) 
    allocate(embtyp(ntyp)) 
    allocate(reppair(npair)) 



    read(lupotin,*) rue
    rue=rue*1.0d-8
    rumax=max(rumax,rue)
    write(6,*) 'Types d_atomes :'
    do i = 1, ntyp
       write(6,*) i
       read (lupotin,*) cm(i),catom(i),ty(i)
       if (rang/=0) cycle
       write (6,*)'type i, cm(i),catom(i),ty(i)'
       write (6, '(I4,2F9.3,A5)') i, cm(i),catom(i),ty(i)
    end do
    if (ntyp.eq.1 )then
       write(6,*) i
       read (lupotin,*) rhotyp(1)%r0,reppair(1)%A,embtyp(1)%ksi,reppair(1)%p,rhotyp(1)%q,roff1(1),roff2(1)
       reppair(1)%r0= rhotyp(1)%r0

       !     if (rang/=0) cycle
       write (6,*)'type i, r0 A ksi P q roff1 roff2'
       write (6, '(I4,7F9.3)') i,rhotyp(1)%r0,reppair(1)%A,embtyp(1)%ksi,reppair(1)%p,rhotyp(1)%q,roff1(1),roff2(1)
       roff1=roff1*1.0d-8 ;roff2=roff2*1.0d-8 
    else
       write(6,*) 'pas programmé ntypw<>1 SM stop'
       stop
    end if
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

    REAL(kind(0.d0)) :: r
    r=sqrt(r2)/A2cm


    IF (present(rho)) rho=exp(2.0*density%q*(1.0-r/density%r0))
    IF (present(drho)) drho=(-2.0*density%q/density%r0*exp(2.0*density%q*(1.0-r/density%r0))) /A2cm
    IF (present(ddrho))then
       write(6,*) 'pas programmé!'
       stop
    end IF

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

    IF (present(embF)) embF=-1.0*eam%ksi*dsqrt(rho)*ev2erg
    IF (present(dembF)) dembF=(-0.5*eam%ksi/dsqrt(rho))*ev2erg
    IF (present(ddembF)) ddembF=(0.25*eam%ksi/rho**1.5)*ev2erg
    IF (present(err)) err=0
    RETURN

  end subroutine extrapolateEam

  !----------------------------------------------

  subroutine extrapolateRep(rep, r2, Erep, dErep, ddErep)
    ! calculate repulsive potential at distance sqrt(r2)
    ! or its first and second derivatives

    implicit none

    type(Rept), intent(in) :: rep
    real(double), intent(in) :: r2
    real(double), intent(out), optional :: Erep, dErep, ddErep

    REAL(double) :: r
    r=sqrt(r2)/A2cm

    IF (present(Erep))  Erep=2.0*rep%A*exp((1.0-r/rep%r0)*rep%p)*ev2erg
    IF (present(dErep)) dErep=-2.0*(rep%A*rep%p/rep%r0)*exp((1.0-r/rep%r0)*rep%p)*ev2erg/A2cm
    IF (present(ddErep)) then
       write(6,*) 'pas programmé!'
       stop
    end IF
    RETURN
    !-----------------------------------
  end subroutine extrapolateRep
end module eam
