module eam
  USE T_kind_param_m

  implicit none

  !Eamtype, Reptype et DensityType definissent les éléments dont sont censés dépendre 
  !les fonctions eam, répulsion et densité. Dans un cas d'alliage on peut avoir 
  !la meme forme analytique mais des valeurs différentes des coefficients pour les 
  !différents types. La liste des coefficients est definie dans les types.
  ! la valeur des coefficients pour ces dfférents types sont définis dans les
  ! routines generRho, generEAm, generRep. Ils sont ensuite utilisé dans les routines extrapolate

  type :: EamT
     real(double),dimension (:),pointer :: feam,dfeam,xg
     real(double)::deltaEAM
     !     integer :: typ
  end type EamT

  type :: RepT
     real(double),dimension (:),pointer :: potr,dpotr,xr
     real(double)::deltaREP
     !     integer :: pair
  end type Rept

  type :: DensityT
     real(double),dimension (:),pointer::rho,drho,xd
     real(double)::deltaRHO
     !     integer :: typ
  end type DensityT

  type(DensityT),dimension(:), pointer :: rhotyp
  type(EamT),dimension(:), pointer :: embtyp
  type(repT),dimension(:), pointer :: reppair


  !  real(double),pointer, dimension(:,:):: xg,xr,xd ! tablezau construit à partir de la grille lue
  !  real(double),pointer, dimension(:,:,:):: eg,vr,dd ! tablezau construit à partir de la grille lue
  !  real(double),pointer, dimension(:):: xspb,yspb,bspb,cspb,dspb
  integer :: nptmax ! nombre de points dans la grille lue

  public ::  extrapolateRho, extrapolateRep, extrapolateEam,inputeam

  real(double) :: ev2erg=1.602d-12, & !conversion eV ->erg
       evA2dyn=1.602d-4 ,& ! conversion ev/A -> dyn
       A2cm =1.0d-8     !conversion A->cm

  !  real(double) :: deltaEAM,deltaREP,deltaRHO
contains


  !---------------------------------------------------------------------------
  subroutine inputeam(ntyp,npair,ntrip,cm,catom,ty,umass,&
       rue,rumax,iewald,l3c,rang,r3cm,roff1,roff2,typ_and_pot,&
       npotmax,ipotentiel,typ_pot_pair,lue_typ,lue_paire,&
       lu_roff_pair,npotentiel,ipo)

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
    integer , intent(in) ::rang,ipotentiel,npotmax,npotentiel
    logical, pointer :: typ_and_pot(:,:) ! typ_and_pot(iti,ipot)=.true. si le type iti interagit (en autres) par le potentiel ipot
    integer,pointer:: typ_pot_pair(:)
    integer, dimension(:,:), pointer  :: ipo			! indice des paires d'atomes
    logical,pointer::lue_typ(:),lue_paire(:),lu_roff_pair(:)

    !local variables
    integer:: i,l,k,iti,n,npt,ipr
    integer :: lupotin=95
    character ::  fnampotin*80
    real(double) :: xmin,xmax,xdum,ruelu,cmr,catomr
    integer,pointer :: typtyp(:),ind_pair(:)
    integer::itir,npair_r,ipair,ntypr,j,itj
    character :: tyr*3

    !    real(double):: deltaEAM, deltaRHO,deltaREP

    fnampotin = 'eamtab.potin'

    lupotin = 95
    open(unit=lupotin, file=fnampotin, status='old')

    iewald=0; l3c=.false.; r3cm=0.

    if (npotentiel.gt.1) then
       read(lupotin,*)ntypr
       allocate (typtyp(ntypr))
       npair_r=  ntypr*(ntypr+1)/2 
       allocate (ind_pair(npair_r))
       write(6,*)'ntypr pour ce pot',ntypr
       read(lupotin,*) rue
       rue=rue*1.0d-8
       if (rang==0)    write(6,*) 'Types d_atomes pour ce potentiel:'
       do i = 1, ntypr
          read (lupotin,*) cmr,catomr,tyr,iti
          typtyp(i)=iti
          if (lue_typ(iti).eqv..true.) then 
             if (rang==0)write(6,*) 'type',iti,'deja lu ; verification de la cohérence'
             !if (cmr*umass.ne.cm(iti))then 
             if (abs(cmr*umass-cm(iti)) > 100.d0*spacing(cm(iti))) then
                if (rang==0)write(6,*) 'pb avec cm'
                stop
             end if
             if (tyr.ne.ty(iti))then
                if (rang==0)write(6,*) 'pb avec ty'
                stop
             end if
          else
             cm(iti)=cmr*umass;ty(iti)=tyr; catom(iti)=catomr; lue_typ(iti)=.true.
          endif

          if (rang/=0) cycle
          write(6,*)'type        cm      catom    ty'
          write (6, '(I4,E12.3,F9.3,A5)') iti,cm(iti),catom(iti),ty(iti)
          typ_and_pot(iti,ipotentiel)=.true.
       end do
       !lecture des roff des paires EAM
       ipair=0
       do i=1,ntypr
          iti=typtyp(i)
          do j=i,ntypr       
             itj=typtyp(j)
             ipr=ipo(iti,itj)
             ipair=ipair+1
             ind_pair(ipair)=ipr
             if (rang==0) write(6,*)'paire l active ipotentiel: ',ipr, ipotentiel
             if(lue_paire(ipr).eqv..true.) then
                write(6,*) rang,'paire l lue deux fois ', ipr,iti,itj
                stop
             end if
             read (lupotin,*) roff1(ipr),roff2(ipr)
             if (rang==0) write(6,*)'roff1 et 2 pour cette paire',ipair,ipr,roff1(ipr),roff2(ipr)
             lu_roff_pair(ipr)=.true.;typ_pot_pair(ipr)=ipotentiel
          end do
       end do

    else
       read(lupotin,*)ntyp
       ntypr=ntyp
       allocate (typtyp(ntypr))     
       npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
       npair_r=npair
       allocate (ind_pair(npair_r))
       call  alloc_typ
       read(lupotin,*) rue
       rue=rue*1.0d-8
       if (rang==0)    write(6,*) 'Types d_atomes :'
       do i = 1, ntyp
          typtyp(i)=i
          read (lupotin,*) cm(i),catom(i),ty(i)
          if (rang/=0) cycle
          write (6, '(I4,2F9.3,A5)') i, cm(i),catom(i),ty(i)
       end do
       do i = 1, npair
          ind_pair(i)=i
          read (lupotin,*) roff1(i),roff2(i)
          if (rang/=0) cycle
          write (6, '(2F9.3)') roff1(i),roff2(i)
       end do

       roff1=roff1*1.0d-8
       roff2=roff2*1.0d-8
       lu_roff_pair(1:npair)=.true. ;typ_pot_pair(:)=ipotentiel
       cm(:ntyp) = cm(:ntyp)*umass
       allocate (typ_and_pot(ntyp,npotmax))
       typ_and_pot(:,:)=.false.
       typ_and_pot(1:ntyp,ipotentiel)=.true.

    end if


    read(lupotin,*)nptmax
    if (rang==0) write(6,*)'nptmax',nptmax
    allocate(rhotyp(ntyp)) 
    allocate(embtyp(ntyp)) 
    allocate(reppair(npair)) 


    !  rhotyp(:,:)%rho=0.;  rhotyp(:,:)%drho=0.;  rhotyp(:,:)%xr=0.
    !  embtyp(:,:)%eam=0.;  embtyp(:,:)%dfeam=0.;  embtyp(:,:)%xg=0.
    !  reppair(:,:)%potr=0.;  reppair(:,:)%dpotr=0.;  reppair(:,:)%xr=0.


    do itir=1,ntypr
       iti=typtyp(itir)
       !lecture de Glue
       read(lupotin,*)n
       if (rang==0) write(6,*)'eam',n,iti
       if(n.ne.itir)then
          write(6,*) rang,' ordre de lecture de EAM stop'
          call arret_ndm
       end if
       read(lupotin,*)npt,embtyp(iti)%deltaEAM
!       if (rang==0) write(6,*)'npt',npt    
       if(npt.gt.nptmax)then
          write(6,*) rang,'nb de points de grille  EAM stop'
          call arret_ndm
       end if
       allocate(embtyp(iti)%feam(nptmax)) 
       allocate(embtyp(iti)%dfeam(nptmax)) 
       allocate(embtyp(iti)%xg(nptmax)) 
       do i=1,nptmax
          if(i.le.npt) then
             read(lupotin,*)embtyp(iti)%xg(i),embtyp(iti)%feam(i),xdum
             !           write(6,*)i,embtyp(iti)%xg(i),embtyp(iti)%feam(i),xdum
          else
             embtyp(iti)%xg(i)=(i-npt)*embtyp(iti)%deltaEAM+ embtyp(iti)%xg(i)
             embtyp(iti)%feam(i)=embtyp(iti)%feam(npt)
          end if

       end do
       do i=1,nptmax-1
          embtyp(iti)%dfeam(i)=embtyp(iti)%feam(i+1)-embtyp(iti)%feam(i)
          !        write(6,*)'eam',iti,i,embtyp(iti)%xg(i),embtyp(iti)%feam(i),embtyp(iti)%dfeam(i)
       end do
       embtyp(iti)%dfeam(nptmax)=0.



       !lecture de dens
       read(lupotin,*)n
       if (rang==0) write(6,*)'dens',n,iti
       if(n.ne.itir)then
          write(6,*) rang,' ordre de lecture de EAM densstop'
          call arret_ndm
       end if
       !     rhotyp(iti)%toto=iti
       read(lupotin,*)npt,rhotyp(iti)%deltaRHO
       if (rang==0) write(6,*)npt,rhotyp(iti)%deltaRHO
       if(npt.ne.nptmax)then
          write(6,*) rang,'nb de points de grille  EAM stop'
          call arret_ndm
       end if
       allocate(rhotyp(iti)%rho(nptmax)) 
       allocate(rhotyp(iti)%drho(nptmax)) 
       allocate(rhotyp(iti)%xd(nptmax)) 

       do i=1,nptmax
          if(i.le.npt) then
             read(lupotin,*)rhotyp(iti)%xd(i),rhotyp(iti)%rho(i),xdum
             !           write(6,*)i,rhotyp(iti)%xd(i),rhotyp(iti)%rho(i),xdum
          else
             rhotyp(iti)%xd(i)=(i-npt)*rhotyp(iti)%deltaRHO+ rhotyp(iti)%xd(npt)
             rhotyp(iti)%rho(i)=rhotyp(iti)%rho(npt)
             rhotyp(iti)%drho(i)=0.
          end if

       end do
       do i=1,nptmax-1
          rhotyp(iti)%drho(i)=rhotyp(iti)%rho(i+1)-rhotyp(iti)%rho(i)
          !        write(6,*)'rho',iti,i,rhotyp(iti)%xd(i),rhotyp(iti)%rho(i),rhotyp(iti)%drho(i)
       end do
       rhotyp(iti)%drho(nptmax)=0.
    end do

    do ipair=1,npair_r
       read(lupotin,*)n
       if(n.ne.ipair)then
          write(6,*) rang, ' ordre de lecture de EAM rep stop'
          call arret_ndm
       end if
       ipr=ind_pair(ipair)
       if (rang==0) write(6,*)'paire eam ; paire complete',n,ipr
       if (rang==0) write(6,*)'rep'

       read(lupotin,*)npt,reppair(ipr)%deltaREP
       if(npt.gt.nptmax)then
          write(6,*) rang,'nb de points de grille  EAM stop'
          call arret_ndm
       end if
       allocate(reppair(ipr)%potr(nptmax)) 
       allocate(reppair(ipr)%dpotr(nptmax)) 
       allocate(reppair(ipr)%xr(nptmax)) 

       do i=1,nptmax
          if(i.le.npt) then
             read(lupotin,*)reppair(ipr)%xr(i),reppair(ipr)%potr(i),xdum
             !          write(6,*)i,reppair(ipr)%xr(i),reppair(ipr)%potr(i),xdum
          else
             reppair(ipr)%xr(i)=(i-npt)*reppair(ipr)%deltaREP+ reppair(ipr)%xr(i)
             reppair(ipr)%potr(i)=reppair(ipr)%potr(npt)
          end if
       end do
       do i=1,nptmax-1
          reppair(ipr)%dpotr(i)=reppair(ipr)%potr(i+1)-reppair(ipr)%potr(i)
          !        write(6,*)'rep',ipr,i,reppair(ipr)%xr(i),reppair(ipr)%potr(i),reppair(ipr)%dpotr(i)
       end do
       reppair(ipr)%dpotr(nptmax)=0.
    end do

    close(lupotin)
!if(associated* (typ_and_pot).eqv..false.), i.e. if npotentiel==1 
    if(associated (typ_and_pot).eqv..false.) then
       allocate (typ_and_pot(ntyp,npotmax))
       typ_and_pot(:,:)=.false.
       typ_and_pot(1:ntyp,ipotentiel)=.true.
    end if


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

    xmax=density%xd(nptmax)
    r=sqrt(r2)/A2cm

    if(r.gt.xmax) then
       Rho=0.0
    else
       k=Int(r/density%deltaRHO)+1
       drk=r/density%deltaRHO+1-k
       Rho=density%rho(k)+drk*density%drho(k)
    end if
    !    write(6,*)'rho',r,k,drk,rho

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

    xmax=eam%xg(nptmax)

    if (rho.lt.0.) then 
       embf=0
       return
    end if
    if(rho.gt.xmax) then
       embF=eam%feam(nptmax)*ev2erg

    else
       k=Int(rho/eam%deltaEAM)+1
       drk=rho/eam%deltaEAM +1 -k
       embf=  ev2erg*(eam%feam(k)+drk*eam%dfeam(k))

    end if
    !    write(6,*)'eam',rho,k,drk,embf/ev2erg


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

    integer :: k
    real(double) :: xmax,r,drk

    xmax=rep%xr(nptmax)

    r=sqrt(r2)/A2cm

    if(r.gt.xmax) then
       Erep=0.0
       return
    else
       k=Int(r/rep%deltaREP)+1
       drk=r/rep%deltaREP+1-k
       Erep=ev2erg*(rep%potr(k)+drk*rep%dpotr(k))
    end if
    !    write(6,*)'rep',r,k,drk,Erep/ev2erg



    RETURN
    !-----------------------------------
  end subroutine extrapolateRep



end module eam
     
