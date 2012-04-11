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


  real(double),pointer, dimension(:,:):: xg,xr,xd ! tablezau construit à partir de la grille lue
  real(double),pointer, dimension(:,:,:):: eg,vr,dd ! tablezau construit à partir de la grille lue
  real(double),pointer, dimension(:):: xspb,yspb,bspb,cspb,dspb
  integer :: nptmax ! nombre de points dans la grille lue

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
    integer:: i,l,k,iti,n,npt
    integer :: lupotin=95
    character ::  fnampotin*80
    real(double) :: xmin,xmax,xdum
    real(double):: deltaEAM, deltaRHO,deltaREP

    fnampotin = 'eamtab.potin'

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
    rumax=max(rue,rumax)
    write(6,*) 'Types d_atomes :'
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

    read(lupotin,*)nptmax
    write(6,*)'nptmax',nptmax
    allocate(xg(nptmax,ntyp) )
    allocate(eg(4,ntyp,nptmax))
    allocate(xr(nptmax,ntyp))
    allocate(Vr(4,npair,nptmax))
    allocate(xd(nptmax,npair))
    allocate(dd(4,ntyp,nptmax))

    allocate(xspb(nptmax))  
    allocate(yspb(nptmax))  
    allocate(bspb(nptmax))  
    allocate(cspb(nptmax))  
    allocate(dspb(nptmax))
    xg=0.0 ; eg=0.
    xr=0.0 ; Vr=0.
    xd=0.0 ; dd=0.

    do iti=1,ntyp
       !lecture de Glue
       read(lupotin,*)n
       write(6,*)'eam'
       if(n.ne.iti)then
          write(6,*)' ordre de lecture de EAM stop'
          stop
       end if
       embtyp(iti)%toto=iti
       read(lupotin,*)npt,deltaEAM
       write(6,*)'npt',npt    
       if(npt.gt.nptmax)then
          write(6,*) 'nb de points de grille  EAM stop'
          stop
       end if
       !     read(lupotin,*)xmin
       !     write(6,*)xmin
       !     read(lupotin,*)xmax
       !     write(6,*)xmax
       xmin=0.d0
       do i=1,nptmax
          xg(i,iti)=xmin+float(i-1)*deltaEAM
          if(i.le.npt) then
             read(lupotin,*)xdum,eg(1,iti,i)
             write(6,*)xdum,  xg(i,iti)
          else
             eg(1,iti,i)=eg(1,iti,npt)
          end if
       end do
       xspb(:)=xg(:,iti)
       yspb(:)=eg(1,iti,:)
       call cspline(nptmax,xspb,yspb,bspb,cspb,dspb)
       eg(2,iti,:)=bspb(:) ; eg(3,iti,:)=cspb(:) ; eg(4,iti,:)=dspb(:) 

       !lecture de dens
       read(lupotin,*)n
       write(6,*)'dens'
       if(n.ne.iti)then
          write(6,*)' ordre de lecture de EAM densstop'
          stop
       end if
       rhotyp(iti)%toto=iti
       read(lupotin,*)npt,deltaRHO
       if(npt.ne.nptmax)then
          write(6,*) 'nb de points de grille  EAM stop'
          stop
       end if
       !     read(lupotin,*)xmin

       !     read(lupotin,*)xmax
       xmin=0.
       do i=1,npt
          xd(i,iti)=xmin+float(i-1)*deltaRHO

          !        if(i.le.npt) then
          read(lupotin,*)xdum, dd(1,iti,i)
          write(6,*)xdum,   dd(1,iti,i)
          !        end if
       end do

       xspb(:)=xd(:,iti)
       yspb(:)=dd(1,iti,:)
       !     do i=1,npt
       !        write(6,*)xd(i,iti),xspb(i),dd(1,iti,i),yspb(i)
       !     end do

       call cspline(nptmax,xspb,yspb,bspb,cspb,dspb)
       dd(2,iti,:)=bspb(:) ; dd(3,iti,:)=cspb(:) ; dd(4,iti,:)=dspb(:) 

    end do

    do iti=1,npair
       read(lupotin,*)n
       write(6,*)'rep'
       if(n.ne.iti)then
          write(6,*)' ordre de lecture de EAM rep stop'
          stop
       end if
       reppair(iti)%toto=iti
       read(lupotin,*)npt,deltaREP
       if(npt.gt.nptmax)then
          write(6,*) 'nb de points de grille  EAM stop'
          stop
       end if
       !     read(lupotin,*)xmin
       !     read(lupotin,*)xmax
       xmin=0.
       do i=1,nptmax
          xr(i,iti)=xmin+float(i-1)*deltaREP
          if(i.le.npt) then
             read(lupotin,*)xdum,vr(1,iti,i)
             !           write(6,*) i,xr(i,iti),vr(1,iti,i)
          end if
       end do


       xspb(:)=xr(:,iti)
       yspb(:)=vr(1,iti,:)
       call cspline(nptmax,xspb,yspb,bspb,cspb,dspb)
       vr(2,iti,:)=bspb(:) ; vr(3,iti,:)=cspb(:) ; vr(4,iti,:)=dspb(:) 
       !     do i=1,nptmax
       !        write(6,*) i,xr(i,iti),vr(1,iti,i),vr(2,iti,i)
       !     end do
       !     call cspline(nptmax,xr(:,iti),vr(1,iti,:),vr(2,iti,:),vr(3,iti,:),vr(4,iti,:))
    end do



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
    !local 
    integer :: iti,n,k
    real(double) :: xmax,xmin,ktor,r,drk

    iti=density%toto
    xmax=xd(nptmax,iti) ; xmin=xd(1,iti)
    ktor=(xmax-xmin)/nptmax


    r=sqrt(r2)/A2cm

    if(r.gt.xmax) then
       rho=0.
       return
    else
       k=Int(r/ktor)+1
       if (k.eq.0) then
          rho=0.0
          return
       else
          drk=r-(k-1)*ktor
          rho=dd(1,iti,k)+dd(2,iti,k)*drk+dd(3,iti,k)*drk**2+dd(4,iti,k)*drk**3 
          !          write(6,*) 'r k rho dd',r, k,rho, dd(1,iti,k)!,dd(2,iti,k),dd(3,iti,k),dd(4,iti,k)

          write(6,*) 'extra rh ',r,k,rho,dd(1,iti,k)
       end if
    end if
    !    CALL rh(r,func,dfunc,d2func)
    !    IF (present(rho)) rho=func
    !    IF (present(drho)) drho=dfunc/A2cm
    !    IF (present(ddrho)) ddrho=d2func/(A2cm**2)


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

    !local 
    integer :: iti,n,k
    real(double) :: xmax,xmin,ktor,drk

    iti=eam%toto
    xmax=xg(nptmax,iti) ; xmin=xg(1,iti)
    ktor=(xmax-xmin)/nptmax

    if(rho.gt.xmax) then
       embF=eg(1,iti,nptmax)
       return
    else
       k=Int(rho/ktor)+1
       drk=rho-(k-1)*ktor
       embf=(eg(1,iti,k)+eg(2,iti,k)*drk+eg(3,iti,k)*drk**2+eg(4,iti,k)*drk**3 )*ev2erg

       write(6,*) 'extra emb ',rho,k,embf/ev2erg,eg(1,iti,k)
    end if



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

    integer :: iti,n,k
    real(double) :: xmax,xmin,ktor,r,drk

    iti=rep%toto
    xmax=xr(nptmax,iti) ; xmin=xr(1,iti)
    ktor=(xmax-xmin)/nptmax

    r=sqrt(r2)/A2cm

    if(r.gt.xmax) then
       Erep=0.0
       return
    else
       k=Int(r/ktor)+1

       drk=r-(k-1)*ktor
       Erep=ev2erg*(vr(1,iti,k)+vr(2,iti,k)*drk+vr(3,iti,k)*drk**2+vr(4,iti,k)*drk**3) 
       !       write(6,*) 'extra rrep ',r,k,Erep/ev2erg,vr(1,iti,k)
    end if




    RETURN
    !-----------------------------------
  end subroutine extrapolateRep



end module eam
     
