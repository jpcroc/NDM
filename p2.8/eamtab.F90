module eam
  USE T_kind_param_m
  USE gen_com_m, ONLY: A2cm
  USE var_pot, ONLY:rhomin,rhomax,lforcetabulate
  USE spline_mod
  USE alloc_typ_mod
  USE arret_ndm_mod

  implicit none

  !Eamtype, Reptype et DensityType definissent les éléments dont sont censés dépendre 
  !les fonctions eam, répulsion et densité. Dans un cas d'alliage on peut avoir 
  !la meme forme analytique mais des valeurs différentes des coefficients pour les 
  !différents types. La liste des coefficients est definie dans les types.
  ! la valeur des coefficients pour ces dfférents types sont définis dans les
  ! routines generRho, generEAm, generRep. Ils sont ensuite utilisé dans les routines extrapolate

  type :: EamT
     real(double),dimension (:),pointer :: feam,xg
     real(double)::deltaEAM
     !     integer :: typ
  end type EamT

  type :: EamTsp
     real(double),dimension (:),pointer :: beam,ceam,deam
  end type EamTsp

  type :: RepT
     real(double),dimension (:),pointer :: potr,xr
     real(double)::deltaREP
  end type Rept

  type :: RepTsp
     real(double),dimension (:),pointer :: bpotr,cpotr,dpotr
  end type RepTsp

  type :: DensityT
     real(double),dimension (:),pointer::rho,xd
     real(double)::deltaRHO
  end type DensityT

  type :: DensityTsp
     real(double),dimension (:),pointer::brho,crho,drho
  end type DensityTsp

  type(DensityT),dimension(:), pointer :: rhotyp
  type(EamT),dimension(:), pointer :: embtyp
  type(repT),dimension(:), pointer :: reppair

  type(DensityTsp),dimension(:), pointer :: SPrhotyp
  type(EamTsp),dimension(:), pointer :: SPembtyp
  type(repTsp),dimension(:), pointer :: SPreppair

  type(DensityT),dimension(:), pointer :: rhotyp_d
  type(EamT),dimension(:), pointer :: embtyp_d
  type(repT),dimension(:), pointer :: reppair_d

  type(DensityTsp),dimension(:), pointer :: SPrhotyp_d
  type(EamTsp),dimension(:), pointer :: SPembtyp_d
  type(repTsp),dimension(:), pointer :: SPreppair_d
  !  real(double),pointer, dimension(:,:):: xg,xr,xd ! tablezau construit à partir de la grille lue
  !  real(double),pointer, dimension(:,:,:):: eg,vr,dd ! tablezau construit à partir de la grille lue
  !  real(double),pointer, dimension(:):: xspb,yspb,bspb,cspb,dspb
  integer :: nptmax ! nombre de points dans la grille lue

  public ::  extrapolateRho, extrapolateRep, extrapolateEam,inputeam


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
    integer:: i,iti,n,npt,ipr
    integer :: lupotin=95
    character ::  fnampotin*80
    real(double) :: xdum,cmr,catomr,drk,erep
    integer,pointer :: typtyp(:),ind_pair(:)
    integer::itir,npair_r,ipair,ntypr,j,itj,k
    character :: tyr*3

    !    real(double):: deltaEAM, deltaRHO,deltaREP

    fnampotin = 'eamtab.potin'
    rhomin=1d30;rhomax=0
    lupotin = 95
    open(unit=lupotin, file=fnampotin, status='old')

    iewald=0; l3c=.false.; r3cm=0.


!PAIR PART 
    if (npotentiel.gt.1) then
       read(lupotin,*)ntypr
       allocate (typtyp(ntypr))
       npair_r=  ntypr*(ntypr+1)/2 
       allocate (ind_pair(npair_r))
       write(6,*)'ntypr for this pot',ntypr
       read(lupotin,*) rue
       rue=rue*A2cm
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

          typ_and_pot(iti,ipotentiel)=.true.
          if (rang/=0) cycle
          write(6,*)'type        cm      catom    ty'
          write (6, '(I4,E12.3,F9.3,A5)') iti,cm(iti),catom(iti),ty(iti)
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
       rue=rue*A2cm
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
          write (6, '(A,2F9.3)') 'ROFF1_2', roff1(i),roff2(i)
       end do

       roff1=roff1*A2cm
       roff2=roff2*A2cm
       lu_roff_pair(1:npair)=.true. ;typ_pot_pair(:)=ipotentiel
       cm(:ntyp) = cm(:ntyp)*umass
       allocate (typ_and_pot(ntyp,npotmax))
       typ_and_pot(:,:)=.false.
       typ_and_pot(1:ntyp,ipotentiel)=.true.

    end if


    read(lupotin,*)nptmax
    if (rang==0) write(6,*)'nptmax in the max number of points on grid  ',nptmax
    allocate(rhotyp(ntyp)) 
    allocate(embtyp(ntyp)) 
    allocate(reppair(npair)) 
    allocate(SPrhotyp(ntyp)) 
    allocate(SPembtyp(ntyp)) 
    allocate(SPreppair(npair)) 

   if (lforcetabulate) then
    allocate(rhotyp_d(ntyp)) 
    allocate(embtyp_d(ntyp)) 
    allocate(reppair_d(npair)) 
    allocate(SPrhotyp_d(ntyp)) 
    allocate(SPembtyp_d(ntyp)) 
    allocate(SPreppair_d(npair))
  end if 

!EMBD EAM PART
    do itir=1,ntypr
       iti=typtyp(itir)
       !lecture de Glue
       read(lupotin,*)n
       if (rang==0) write(6,*)'EAM',n,iti
       if(n.ne.itir)then
          write(6,*) rang,' ordre de lecture de EAM stop'
          call arret_ndm
       end if
       read(lupotin,*)npt,embtyp(iti)%deltaEAM
       if (lforcetabulate) embtyp_d(iti)%deltaEAM=embtyp(iti)%deltaEAM
       if (rang==0) write(6,*)'EAM number of points in potin and the step',npt,embtyp(iti)%deltaEAM     
       if(npt.gt.nptmax)then
          write(6,*) rang,'nb de points de grille  EAM stop'
          call arret_ndm
       end if
     
       allocate(embtyp(iti)%xg(nptmax)) 
       allocate(embtyp(iti)%feam(nptmax)) 
       allocate(SPembtyp(iti)%beam(nptmax)) 
       allocate(SPembtyp(iti)%ceam(nptmax)) 
       allocate(SPembtyp(iti)%deam(nptmax)) 

      if (lforcetabulate) then
       allocate(embtyp_d(iti)%xg(nptmax)) 
       allocate(embtyp_d(iti)%feam(nptmax)) 
       allocate(SPembtyp_d(iti)%beam(nptmax)) 
       allocate(SPembtyp_d(iti)%ceam(nptmax)) 
       allocate(SPembtyp_d(iti)%deam(nptmax))
      end if
! write(6,*)'NPT NPTMAX',npt,nptmax
       do i=1,nptmax
          if(i.le.npt) then
             !newCOS
            if (lforcetabulate) then
             read(lupotin,*)embtyp(iti)%xg(i),embtyp(iti)%feam(i),embtyp_d(iti)%feam(i)
            else 
             read(lupotin,*)embtyp(iti)%xg(i),embtyp(iti)%feam(i)
            end if
          else
             embtyp(iti)%xg(i)=(i-npt)*embtyp(iti)%deltaEAM+ embtyp(iti)%xg(i)
             embtyp(iti)%feam(i)=embtyp(iti)%feam(npt)
             !newCOS
             if (lforcetabulate) then
                embtyp_d(iti)%feam(i)=embtyp_d(iti)%feam(npt)
             end if
          end if

       end do
       rhomin=min(rhomin,embtyp(iti)%xg(1))
       rhomax=max(rhomax,embtyp(iti)%xg(npt))
!       write(6,*)'rhomin rhomax',rhomin,rhomax

       call cspline (nptmax,embtyp(iti)%xg,embtyp(iti)%feam,SPembtyp(iti)%beam,SPembtyp(iti)%ceam,SPembtyp(iti)%deam)
       if (lforcetabulate) then
        embtyp_d(iti)%xg=embtyp(iti)%xg
        call cspline (nptmax,embtyp_d(iti)%xg,embtyp_d(iti)%feam,SPembtyp_d(iti)%beam,SPembtyp_d(iti)%ceam,SPembtyp_d(iti)%deam)
       end if




!DENS PART
       !lecture de dens
       read(lupotin,*)n
       if (rang==0) write(6,*)'dens',n,iti
       if(n.ne.itir)then
          write(6,*) rang,' ordre de lecture de EAM densstop'
          call arret_ndm
       end if
       !     rhotyp(iti)%toto=iti
       read(lupotin,*)npt,rhotyp(iti)%deltaRHO
       if (lforcetabulate) rhotyp_d(iti)%deltaRHO=rhotyp(iti)%deltaRHO

       if (rang==0) write(6,*)'RHO potin points and the step: ', npt,rhotyp(iti)%deltaRHO
       if(npt.ne.nptmax)then
          write(6,*) rang,'nb de points de grille  EAM stop'
          call arret_ndm
       end if
     
       allocate(rhotyp(iti)%xd(nptmax)) 
       allocate(rhotyp(iti)%rho(nptmax)) 
       allocate(SPrhotyp(iti)%brho(nptmax)) 
       allocate(SPrhotyp(iti)%crho(nptmax)) 
       allocate(SPrhotyp(iti)%drho(nptmax)) 
      
      if (lforcetabulate) then
       allocate(rhotyp_d(iti)%xd(nptmax)) 
       allocate(rhotyp_d(iti)%rho(nptmax)) 
       allocate(SPrhotyp_d(iti)%brho(nptmax)) 
       allocate(SPrhotyp_d(iti)%crho(nptmax)) 
       allocate(SPrhotyp_d(iti)%drho(nptmax)) 
      end if
     do i=1,nptmax
          if(i.le.npt) then
            if (lforcetabulate) then
             read(lupotin,*)rhotyp(iti)%xd(i),rhotyp(iti)%rho(i),rhotyp_d(iti)%rho(i)
            else 
             read(lupotin,*)rhotyp(iti)%xd(i),rhotyp(iti)%rho(i)
            end if
           else
             rhotyp(iti)%xd(i)=(i-npt)*rhotyp(iti)%deltaRHO+ rhotyp(iti)%xd(npt)
             rhotyp(iti)%rho(i)=rhotyp(iti)%rho(npt)
            if (lforcetabulate) then
              rhotyp_d(iti)%rho(i)=rhotyp_d(iti)%rho(npt)
            end if 
          end if

       end do
        call cspline (nptmax,rhotyp(iti)%xd,rhotyp(iti)%rho,  SPrhotyp(iti)%brho,  SPrhotyp(iti)%crho,  SPrhotyp(iti)%drho  )
       if (lforcetabulate) then
        rhotyp_d(iti)%xd=rhotyp(iti)%xd
        call cspline (nptmax,rhotyp_d(iti)%xd,rhotyp_d(iti)%rho,SPrhotyp_d(iti)%brho,SPrhotyp_d(iti)%crho,SPrhotyp_d(iti)%drho)
       end if
    end do

!PAIR PART 

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
       if (lforcetabulate) reppair_d(ipr)%deltaREP=reppair(ipr)%deltaREP

       if(npt.gt.nptmax)then
          write(6,*) rang,'nb de points de grille  EAM stop'
          call arret_ndm
       end if
       allocate(reppair(ipr)%xr(nptmax)) 
       allocate(reppair(ipr)%potr(nptmax)) 
       allocate(SPreppair(ipr)%bpotr(nptmax)) 
       allocate(SPreppair(ipr)%cpotr(nptmax)) 
       allocate(SPreppair(ipr)%dpotr(nptmax)) 
       
  
      if (lforcetabulate) then
       allocate(reppair_d(ipr)%xr(nptmax)) 
       allocate(reppair_d(ipr)%potr(nptmax)) 
       allocate(SPreppair_d(ipr)%bpotr(nptmax)) 
       allocate(SPreppair_d(ipr)%cpotr(nptmax)) 
       allocate(SPreppair_d(ipr)%dpotr(nptmax)) 
      end if

       do i=1,nptmax
          if(i.le.npt) then
           if (lforcetabulate) then        
             read(lupotin,*)reppair(ipr)%xr(i),reppair(ipr)%potr(i),reppair_d(ipr)%potr(i)
           else
             read(lupotin,*)reppair(ipr)%xr(i),reppair(ipr)%potr(i)
           end if
          else
             reppair(ipr)%xr(i)=(i-npt)*reppair(ipr)%deltaREP+ reppair(ipr)%xr(i)
             reppair(ipr)%potr(i)=reppair(ipr)%potr(npt)
             if (lforcetabulate) then
              reppair_d(ipr)%potr(i)=reppair_d(ipr)%potr(npt)
             end if

          end if
       end do
       call cspline (nptmax,reppair(ipr)%xr,reppair(ipr)%potr,SPreppair(ipr)%bpotr,& 
             SPreppair(ipr)%cpotr,SPreppair(ipr)%dpotr)
       if (lforcetabulate) then
        reppair_d(ipr)%xr=reppair(ipr)%xr
        call cspline (nptmax,reppair_d(ipr)%xr,reppair_d(ipr)%potr,SPreppair_d(ipr)%bpotr,&
            SPreppair_d(ipr)%cpotr,SPreppair_d(ipr)%dpotr)
       end if

    end do

    close(lupotin)
!if(associated* (typ_and_pot).eqv..false.), i.e. if npotentiel==1 
    if(associated (typ_and_pot).eqv..false.) then
       allocate (typ_and_pot(ntyp,npotmax))
       typ_and_pot(:,:)=.false.
       typ_and_pot(1:ntyp,ipotentiel)=.true.
    end if
    rumax=max(rue,rumax)


    return

  end subroutine inputeam

  !-------------------------------------------

  subroutine extrapolateRho(density,SPdensity, r2, rho)
    ! calculate electronic density at distance sqrt(r)
    ! or its first and second derivatives

    USE gen_com_m, ONLY:  ev2erg
    implicit none

    type(DensityT), intent(in) :: density 
    type(DensityTsp), intent(in) :: SPdensity 
    real(kind(0.d0)), intent(in) :: r2
    real(kind(0.d0)), intent(out), optional :: rho
    !local
    integer:: kr 
    real(double) :: xmax,r,drk
    
    kr=0
    xmax=density%xd(nptmax)
    r=sqrt(r2)/A2cm

    if(r.gt.xmax) then
       Rho=0.0
    else
       kr=Int(r/density%deltaRHO)+1
       drk=r+(1-kr)*density%deltaRHO
       rho = density%rho(kr)+drk*(SPdensity%brho(kr)+drk*(SPdensity%crho(kr)+drk*SPdensity%drho(kr)))

    end if
    RETURN


  end subroutine extrapolateRho


  !---------------------------------------------------------------------------


  subroutine extrapolateEam(eam, SPeam,rho, embF)

    USE gen_com_m, ONLY:  ev2erg
    implicit none

    type(EamT), intent(in) :: eam
    type(EamTsp), intent(in) :: SPeam

    real(double), intent(in) :: rho
    real(double), intent(out), optional :: embF

    !local 
    integer :: k
    real(double) :: xmax,drk

    xmax=eam%xg(nptmax)
    k=0
    if (rho.lt.0.) then 
       embf=0
       return
    end if
    if(rho.gt.xmax) then
       embF=eam%feam(nptmax)*ev2erg

    else
!       write(6,*)
       k=Int((rho-eam%deltaEAM/1d10)/eam%deltaEAM)+1
!       write(6,*) k,rho,rho/eam%deltaEAM
       drk=rho +(1 -k)*eam%deltaEAM
       Embf = ev2erg*(eam%feam(k)+drk*(SPeam%beam(k)+drk*(SPeam%ceam(k)+drk*SPeam%deam(k))))
    end if


    RETURN

  end subroutine extrapolateEam

  !----------------------------------------------

  subroutine extrapolateRep(rep, SPrep,r2,Erep)
    ! calculate repulsive potential at distance sqrt(r2)
    ! or its first and second derivatives

    USE gen_com_m, ONLY:  ev2erg
    implicit none

    type(RepT), intent(in) :: rep
    type(RepTsp), intent(in) :: SPrep
    real(kind(0.d0)), intent(in) :: r2
    real(kind(0.d0)), intent(out), optional :: Erep


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
       drk=r+(1-k)*rep%deltaREP
       Erep = ev2erg*(rep%potr(k) +drk*(SPrep%bpotr(k) +drk*(SPrep%cpotr(k) +drk*SPrep%dpotr(k))))

    end if



    RETURN
    !-----------------------------------
  end subroutine extrapolateRep



end module eam
     
