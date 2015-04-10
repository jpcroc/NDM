module SMjuli
  USE T_kind_param_m
  USE gen_com_m, ONLY: ev2erg,A2cm
  implicit none


  type :: EamTjl  ! liste les param�tres des fonctions glue
     real(double)::ksi
  end type EamTjl

  type :: RepTjl ! liste les param�tres des fonctions r�pulsions
     real(double)::A
     real(double)::B
     real(double)::K
     real(double)::rc
  end type RepTjl


  type :: DensityTjl ! liste les param�tres des fonctions densit�s
     real(double)::C
     real(double)::D
     real(double)::K
     real(double)::rc
  end type DensityTjl

  type(DensityTjl),dimension(:), pointer :: rhotypjl  !
  type(EamTjl),dimension(:), pointer :: embtypjl      !
  type(repTjl),dimension(:), pointer :: reppairjl     !

  real(double):: alphaPbeta, beta

  public ::  extrapolateRhojl, extrapolateRepjl, extrapolateEamjl,inputeamjl

       
contains


  !---------------------------------------------------------------------------
  subroutine inputeamjl(ntyp,npair,ntrip,cm,catom,ty,umass,rue,&
rumax,iewald,l3c,rang,r3cm,roff1,roff2,typ_and_pot,npotmax,&
ipotentiel,typ_pot_pair)

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
    integer , intent(in) ::rang,ipotentiel,npotmax
    logical, pointer :: typ_and_pot(:,:) ! typ_and_pot(iti,ipot)=.true. si le type iti interagit (en autres) par le potentiel ipot
    integer, pointer:: typ_pot_pair(:) ! donne le type d'interaction de la paire
    !local variables
    integer:: i
    integer :: lupotin=95
    character ::  fnampotin*80


    !  fnampotin = 'SM.potin'

    !  lupotin = 95
    !  open(unit=lupotin, file=fnampotin, status='old')

    iewald=0; l3c=.false.; r3cm=0.

    !  read(lupotin,*)ntyp
    ntyp=2
    npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
    call  alloc_typ
    allocate(rhotypjl(npair)) 
    allocate(embtypjl(ntyp)) 
    allocate(reppairjl(npair)) 
    allocate (typ_and_pot(1:ntyp,npotmax))
    typ_and_pot=.false. ; typ_and_pot(:,ipotentiel)=.true.
    allocate (typ_pot_pair(1:npair))
    typ_pot_pair(1)=ipotentiel
    typ_pot_pair(2)=ipotentiel

    embtypjl%ksi=1.0

    reppairjl(1)%A=2.9296875
    reppairjl(1)%B=2.58787395638939
    reppairjl(1)%K=0.1
    reppairjl(1)%rc=7.0

    rhotypjl(1)%C=1.32308467741935
    rhotypjl(1)%D=4.3672464262563
    rhotypjl(1)%K=reppairjl(1)%K
    rhotypjl(1)%rc=reppairjl(1)%rc

    reppairjl(2)%A=3.24589393669854
    reppairjl(2)%B=2.05679804919117
    reppairjl(2)%K=0.1
    reppairjl(2)%rc=3.5

    !rhotypjl(2)%C=0.82303818052368
    rhotypjl(2)%C=2*0.82303818052368
    rhotypjl(2)%D=4.15482225815134
    !rhotypjl(2)%K=reppairjl(2)%K
    rhotypjl(2)%K=2*reppairjl(2)%K
    rhotypjl(2)%rc=reppairjl(2)%rc

    rhotypjl(3)%C=0.
    rhotypjl(3)%D=0.0
    rhotypjl(3)%K=0.0
    rhotypjl(3)%rc=0.0

    reppairjl(3)%A=0.0
    reppairjl(3)%B=0.0
    reppairjl(3)%K=0.0
    reppairjl(3)%rc=1.5



    beta=14.59345494373451
    !beta=0.0
    alphaPbeta=1.80853303846249**beta


    !alphaPbeta=1.0


    !  read(lupotin,*) rue
    rue=7.0*1.0d-8
    rumax=max(rue,rumax)
    write(6,*) 'Types d_atomes :'

    cm(1)=91.22 ; catom(1)=40. ; ty(1)='Zr'
    cm(2)=12.01115 ; catom(2)=6. ; ty(2)='C '
    roff1(1)=1.1 ; roff2(1)=1.4
    roff1(2)=1.0 ; roff2(2)=1.2
!    roff1(3)=0.95 ; roff2(3)=1.12
    roff1=roff1*1.0d-8 ; roff2=roff2*1.0d-8

    do i = 1, ntyp
       write(6,*) i
       !     read (lupotin,*) cm(i),catom(i),ty(i)
       if (rang/=0) cycle
       write (6,*)'type i, cm(i),catom(i),ty(i)'
       write (6, '(I4,2F9.3,A5)') i, cm(i),catom(i),ty(i)
    end do
    write(6,*) '****potentiel de Ju Li pour ZrC ****'
    write(6,*) '****Zr=1 C =2 ****'

    cm(:ntyp) = cm(:ntyp)*umass
    close(lupotin)
    return

  end subroutine inputeamjl

  !-------------------------------------------

  subroutine extrapolateRhojl(density, r2, rho, drho, ddrho)
    ! calculate electronic density at distance sqrt(r2)
    ! or its first and second derivatives

    implicit none

    type(DensityTjl), intent(in) :: density
    real(kind(0.d0)), intent(in) :: r2
    real(kind(0.d0)), intent(out), optional :: rho, drho, ddrho

    REAL(kind(0.d0)) :: r,aux1
    r=sqrt(r2)/A2cm
    r=max(r,1.d0)


    IF (present(rho))then 
       if (r.ge.density%rc) then
          rho=0.
       else
          
          aux1=density%K/(r-density%rc)
          if(aux1.le.50.0)then
             rho=exp(density%C*(density%D-r)+density%K/(r-density%rc))
          else
             rho=0.0
          end if
       end if
    end IF
    IF (present(drho)) then
       write(6,*) 'pas programm�!'
       stop
    end IF
    IF (present(ddrho))then
       write(6,*) 'pas programm�!'
       stop
    end IF

    RETURN
  end subroutine extrapolateRhojl
  !---------------------------------------------------------------------------
  subroutine extrapolateEamjl(eam, rho, embF, dembF, ddembF, err)
    ! calculate eam function for electronic density rho
    !   or its first and second derivatives
    !   err= 0 if everyting ok
    !       -1 if density too large for extrapolation
    implicit none

    type(EamTjl), intent(in) :: eam

    real(double), intent(in) :: rho
    real(double), intent(out), optional :: embF, dembF, ddembF
    integer, intent(out), optional :: err

    IF (present(embF)) embF=-1.0*eam%ksi*dsqrt(rho)*ev2erg
    IF (present(dembF)) dembF=(-0.5*eam%ksi/dsqrt(rho))*ev2erg
    IF (present(ddembF)) ddembF=(0.25*eam%ksi/rho**1.5)*ev2erg
    IF (present(err)) err=0
    RETURN

  end subroutine extrapolateEamjl

  !----------------------------------------------

  subroutine extrapolateRepjl(rep, r2, Erep, dErep, ddErep)
    ! calculate repulsive potential at distance sqrt(r2)
    ! or its first and second derivatives

    implicit none

    type(Reptjl), intent(in) :: rep
    real(double), intent(in) :: r2
    real(double), intent(out), optional :: Erep, dErep, ddErep

    REAL(double) :: r,aux1,aux2
    r=sqrt(r2)/A2cm

    IF (present(Erep)) then
          if (r.ge.rep%rc) then
             Erep=0.
          else
             
             aux1=abs(rep%K/(r-rep%rc))
             if(aux1.le.50.0)then
!                write(6,*)r-rep%rc,aux1
                Erep=exp(rep%A*(rep%B-r)+rep%K/(r-rep%rc))*ev2erg
             else
                Erep=0.0
             end if
          end if

       end IF
    IF (present(dErep))  then
       write(6,*) 'pas programm�!'
       stop
    end IF
    IF (present(ddErep)) then
       write(6,*) 'pas programm�!'
       stop
    end IF
    RETURN
    !-----------------------------------
  end subroutine extrapolateRepjl
end module SMjuli
