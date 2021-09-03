module eamerco
  USE T_kind_param_m
  USE gen_com_m, ONLY: ev2erg,A2cm
  USE alloc_typ_mod,only: alloc_typ
  use var_pot,only:rhomax,rhomin
  implicit none


  public ::  extrapolateRhoerco, extrapolateReperco, extrapolateEamerco,inputeamerco

       
contains


  !---------------------------------------------------------------------------
  subroutine inputeamerco(ntyp,npair,ntrip,cm,catom,ty,umass,rue,rumax,&
                iewald,l3c,rang,r3cm,roff1,roff2,&
                &typ_and_pot,npotmax,&
                ipotentiel,typ_pot_pair,lue_typ,lue_paire,lu_roff_pair,&
                npotentiel,ipo)



    !

    integer, intent(out) :: ntyp                           !nb de type
    integer, intent(out) :: npair                  ! = ntyp*(ntyp+1)/2
    integer, intent(out) :: ntrip                  ! = ntyp*ntyp *(ntyp+1)/2
    real(double), dimension(:), allocatable :: cm, catom,roff1,roff2
    character , dimension(:), allocatable  :: ty*3
    real(double), intent(in) ::umass 
    real(double), intent(out) :: rue,rumax,r3cm
    integer, intent(out) :: iewald
    logical, intent(out) :: l3c
    integer , intent(in) ::rang

    integer::ipotentiel,npotmax,npotentiel
    logical, allocatable :: typ_and_pot(:,:) ! typ_and_pot(iti,ipot)=.true. si le type iti interagit (en autres) par le potentiel ipot
    integer,allocatable:: typ_pot_pair(:)
    integer, dimension(:,:), allocatable  :: ipo			! indice des paires d'atomes
    logical,allocatable::lue_typ(:),lue_paire(:),lu_roff_pair(:)


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
!    allocate(rhotyp(ntyp)) 
!    allocate(embtyp(ntyp)) 
!    allocate(reppair(npair)) 


    allocate (typ_and_pot(1,npotmax))

    typ_and_pot(1,ipotentiel)=.true.
    typ_pot_pair(1)=ipotentiel


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
    rhomin(:)=1d30;rhomax(:)=0


    close(lupotin)

    return

  end subroutine inputeamerco

  !-------------------------------------------

  subroutine extrapolateRhoerco(r2, rho)
    ! calculate electronic density at distance sqrt(r)
    ! or its first and second derivatives

    implicit none

    real(kind(0.d0)), intent(in) :: r2
    real(kind(0.d0)), intent(out) :: rho

    !  Ercolessi potential
    REAL(kind(0.d0)) :: r, func, dfunc, d2func
    r=sqrt(r2)/A2cm

    CALL rh(r,func,dfunc,d2func)
     rho=func
!    IF (present(drho)) drho=dfunc/A2cm
!    IF (present(ddrho)) ddrho=d2func/(A2cm**2)

    RETURN


  end subroutine extrapolateRhoerco


  !---------------------------------------------------------------------------


  subroutine extrapolateEamerco(rho, embF)
    ! calculate eam function for electronic density rho
    !   or its first and second derivatives
    !   err= 0 if everyting ok
    !       -1 if density too large for extrapolation
    implicit none


    real(double), intent(in) :: rho
    real(double), intent(out) :: embF



    ! Debugging with Ercolessi potential
    REAL(kind(0.d0)) :: func, dfunc, d2func
    CALL uu(rho, func, dfunc, d2func)
     embF=func*ev2erg
!    IF (present(dembF)) dembF=dfunc*ev2erg
!    IF (present(ddembF)) ddembF=d2func*ev2erg
!    IF (present(err)) err=0
    RETURN

  end subroutine extrapolateEamerco

  !----------------------------------------------

  subroutine extrapolateReperco(r2, Erep)
    ! calculate repulsive potential at distance sqrt(r2)
    ! or its first and second derivatives

    implicit none

    real(kind(0.d0)), intent(in) :: r2
    real(kind(0.d0)), intent(out) :: Erep


    !-----------------------------------
    ! Debugging with Ercolessi potential
    REAL(kind(0.d0)) :: r, func, dfunc, d2func
    r=sqrt(r2)/A2cm
    CALL v2(r, func, dfunc, d2func)
     Erep=func*ev2erg
!    IF (present(dErep)) dErep=dfunc*ev2erg/A2cm
!    IF (present(ddErep)) ddErep=d2func*ev2erg/(A2cm**2)
    RETURN
    !-----------------------------------
  end subroutine extrapolateReperco
end module eamerco
