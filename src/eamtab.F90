module eam
  USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m
  USE gen_com_m, ONLY: A2cm,rang,lopt,ev2erg,low_limit
  USE var_pot, ONLY:rhomin,rhomax,lforcetabulate,q,alpha,precisew,ncouc3,ncoucx,ncoucy,ncoucz,&
       &  kpmex, kpmey, kpmez,ipotrep,rcwolf,ngrid
  USE spline_mod,only: cspline
  USE alloc_typ_mod,only: alloc_typ
  USE arret_ndm_mod,only: arret_ndm

  USE calerf_mod
  use input_pair_mod,only:checklu,checkewald
  implicit none

  !Eamtype, Reptype et DensityType definissent les éléments dont sont censés dépendre 
  !les fonctions eam, répulsion et densité. Dans un cas d'alliage on peut avoir 
  !la meme forme analytique mais des valeurs différentes des coefficients pour les 
  !différents types. La liste des coefficients est definie dans les types.
  ! la valeur des coefficients pour ces dfférents types sont définis dans les
  ! routines generRho, generEAm, generRep. Ils sont ensuite utilisé dans les routines extrapolate
  type crg_T
     real(double),allocatable,dimension(:):: D,gam,A,rho,C,r0,G,n,densmax
   contains
     procedure,pass::alloc=> alloc_crg
  end type crg_T
  
  type :: EamT
     real(double),dimension (:),allocatable :: feam,xg
     real(double)::deltaEAM
     !     integer :: typ
  end type EamT

  type :: EamTsp
     real(double),dimension (:),allocatable :: beam,ceam,deam
  end type EamTsp

  type :: RepT
     real(double),dimension (:),allocatable :: potr,xr
     real(double)::deltaREP
  end type Rept

  type :: RepTsp
     real(double),dimension (:),allocatable :: bpotr,cpotr,dpotr
  end type RepTsp

  type :: DensityT
     real(double),dimension (:),allocatable::rho,xd
     real(double)::deltaRHO
  end type DensityT

  type :: DensityTsp
     real(double),dimension (:),allocatable::brho,crho,drho
  end type DensityTsp

  type(DensityT),dimension(:), allocatable :: rhotyp
  type(EamT),dimension(:), allocatable :: embtyp
  type(repT),dimension(:), allocatable :: reppair

  type(DensityTsp),dimension(:), allocatable :: SPrhotyp
  type(EamTsp),dimension(:), allocatable :: SPembtyp
  type(repTsp),dimension(:), allocatable :: SPreppair

  type(DensityT),dimension(:), allocatable :: rhotyp_d
  type(EamT),dimension(:), allocatable :: embtyp_d
  type(repT),dimension(:), allocatable :: reppair_d

  type(DensityTsp),dimension(:), allocatable :: SPrhotyp_d
  type(EamTsp),dimension(:), allocatable :: SPembtyp_d
  type(repTsp),dimension(:), allocatable :: SPreppair_d
  !  real(double),allocatable, dimension(:,:):: xg,xr,xd ! tablezau construit à partir de la grille lue
  !  real(double),allocatable, dimension(:,:,:):: eg,vr,dd ! tablezau construit à partir de la grille lue
  !  real(double),allocatable, dimension(:):: xspb,yspb,bspb,cspb,dspb
  integer :: nptmax ! nombre de points dans la grille lue

  public ::  extrapolateRho, extrapolateRep, extrapolateEam,inputeam

  type(crg_T)::crg

  !  real(double) :: deltaEAM,deltaREP,deltaRHO
contains

  subroutine alloc_crg(crg,ntyp)
    class(crg_T)::crg
    integer,intent(in)::ntyp

    integer::npair
    npair=ntyp*(ntyp+1)/2
    allocate(crg%G(ntyp))
    allocate(crg%densmax(ntyp))
    allocate(crg%n(ntyp))
    
    allocate(crg%D(npair))
    allocate(crg%gam(npair))
    allocate(crg%A(npair))
    allocate(crg%rho(npair))
    allocate(crg%C(npair))
    allocate(crg%r0(npair))
  end subroutine alloc_crg
  !---------------------------------------------------------------------------
  subroutine inputeam(ntyp,npair,ntrip,cm,catom,ty,umass,&
       rue,rumax,iewald,l3c,r3cm,roff1,roff2,typ_and_pot,&
       npotmax,ipotentiel,typ_pot_pair,lue_typ,lue_paire,&
       lu_roff_pair,npotentiel,ipo)

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
    integer , intent(in) ::ipotentiel,npotmax,npotentiel
    logical, allocatable :: typ_and_pot(:,:) ! typ_and_pot(iti,ipot)=.true. si le type iti interagit (en autres) par le potentiel ipot
    integer,allocatable:: typ_pot_pair(:)
    integer, dimension(:,:), allocatable  :: ipo			! indice des paires d'atomes
    logical,allocatable::lue_typ(:),lue_paire(:),lu_roff_pair(:)
    real(double)::rf1,rf2
    integer::l,lect_paire,nb_paire_a_lire,it1,it2

    !local variables
    integer:: i,iti,n,npt,ipr
    integer :: lupotin=95
    character ::  fnampotin*80
    real(double) :: cmr,catomr,precis,qr
    integer,allocatable :: typtyp(:)
    integer::itir,ipair,ntypr
    logical::lok
    character :: tyr*3
    namelist /ewald/ rue, alpha, precis, ncouc3, ncoucx, ncoucy, ncoucz,&
         kpmex, kpmey, kpmez, lopt,iewald,ipotrep
    !EWALD
    !  rumax=0.0
    r3cm=0.0
    rue = 0.0
    alpha = 0.0
    precis=0.0
    precisew = 0.0
    ncouc3 = 0
    ncoucx = 0
    ncoucy = 0
    ncoucz = 0
    kpmex = 0
    kpmey = 0
    kpmez = 0
    lopt=.FALSE.
    if (ipotentiel==10) then
       ipotrep=1
    else
       ipotrep=0
    end if

    !    real(double):: deltaEAM, deltaRHO,deltaREP
    select case (ipotentiel)
    case(10)
       fnampotin = 'eamtab.potin'
    case(16)
       fnampotin = 'CRG.potin'
    end select

    lupotin = 95
    open(unit=lupotin, file=fnampotin, status='old')

    iewald=0; l3c=.false.; r3cm=0.


    !PAIR PART 
    if (npotentiel.gt.1) then
       read(lupotin,*)ntypr 
       allocate (typtyp(ntypr))
       if (rang==0) write(6,*)'ntypr for this pot',ntypr
    else
       read(lupotin,*)ntyp
       ntypr=ntyp
       allocate (typtyp(ntypr))     
       npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
       nb_paire_a_lire=npair
       call  alloc_typ
       typ_pot_pair(1:npair)=ipotentiel
       if (rang==0) write(6,*)'NTYP NPAIR ',ntyp,npair
    end if
    
    rhomin(:)=1d30;rhomax(:)=0
    select case(ipotentiel)
    case(16)
       call crg%alloc(ntyp)

       read (lupotin, nml=ewald)            ! lecture de la namelist ewald
       call checkewald(iewald)
       if (iewald==0) then
          if (rang==0)then
             write (6, *) '-*-*-*-* PAS DE SOMMATION D-EWALD *-*-*-*-'
             write (6, *) '-*-*-*-* IPOTENTIEL=16. STOP ! ipotentiel-> 10 ! *-*-*-*-'
          end if
          call arret_ndm
       elseif (iewald==1) then
          if (rang==0) then
             write (6, *) '-*-*-*-*-* SOMMATION D-EWALD CLASSIQUE *-*-*-*-*-'
             !                if (npotentiel.gt.1) write(6,*)'FONCTIONNEMENT NON GARANTI!!!'
          end if
       elseif (iewald==2) then
          if (rang==0) write (6, *) '-*-*-*-*-* SOMMATION D-EWALD METHODE PME *-*-*-*-*-'
          call arret_ndm
       elseif (iewald==3) then
          if (rang==0) write (6, *) '-*-*-*-*-* SOMMATION DE WOLF *-*-*-*-*-'
          if ((rue==0).or.(alpha==0)) then
             if (rang==0) write (6, *) 'RUE and alpha must be set'
             call arret_ndm
          else
             if (rang==0)write(6,*)'RcWolf=',rcwolf,' alpha= ',alpha
          end if
          if (RcWolf==0) RcWolf=Rue
       else
          write (6, *) rang, 'Valeur de iewald erronee : iewald=',iewald
          call arret_ndm
       endif
       
       precisew=precis
       if (rue==0) then
          if (rang==0) write (6, *) '-*-*-*-*-* RUE must be NON ZERO *-*-*-*-*-'
          call arret_ndm
       end if
       
    case(10)
       read(lupotin,*) rue
       call checkewald(iewald)
    end select
    
    rue=rue*A2cm
    if (rang==0)    write(6,*) 'Types d_atomes pour ce potentiel:'
    if (iewald==0) then
       if (rang==0) write (6, *) 'CM, masse,type, NUMERO DU TYPE D ATOME'
    else
       if (rang==0) write (6, *) 'CHARGE,CM, masse,type, NUMERO DU TYPE D ATOME'
    end if
    if (npotentiel.gt.1) then
       do i = 1, ntypr
          select case(ipotentiel)
          case(10)
             if (iewald.gt.0) then
                read (lupotin,*) qr,cmr,catomr,tyr,iti
                call checklu(iti,tyr,cmr,catomr,qr)
             else
                read (lupotin,*) cmr,catomr,tyr,iti
                call checklu(iti,tyr,cmr,catomr)
                cm(iti)=cmr*umass;ty(iti)=tyr; catom(iti)=catomr; lue_typ(iti)=.true.
             end if
          case(16)
             read (lupotin,*) qr,cmr,catomr,tyr,iti
             call checklu(iti,tyr,cmr,catomr,qr)
             ty(iti)=tyr; catom(iti)=catomr; lue_typ(iti)=.true.;q(iti)=qr;cm(iti)=cmr*umass
          end select
          typtyp(i)=iti
          typ_and_pot(iti,ipotentiel)=.true.
          if (rang/=0) cycle
          if (rang==0)write(6,*)'type        cm      catom    ty'
          if (rang==0) write (6, '(I4,E12.3,F9.3,A5)') iti,cm(iti),catom(iti),ty(iti)
          if((iewald.ne.0) .and.(rang==0)) write(6,*)'charge = ',q(iti)
       end do
    else
       do i = 1, ntyp
          typtyp(i)=i
          select case(ipotentiel)
          case(10)
             read (lupotin,*) cm(i),catom(i),ty(i)
          case(16)
             read (lupotin,*) q(i),cm(i),catom(i),ty(i)
          end select
          if (rang/=0) cycle
          write (6, '(I4,2F9.3,A5)') i, cm(i),catom(i),ty(i)
          if((ipotentiel==16) .and.(rang==0)) write(6,*)'charge = ',q(i)
       end do
       cm(:ntyp) = cm(:ntyp)*umass
       allocate (typ_and_pot(ntyp,npotmax))
       typ_and_pot(:,:)=.false.
       typ_and_pot(1:ntyp,ipotentiel)=.true.
    end if
    if (npotentiel.gt.1) then   
       read(lupotin,*)nb_paire_a_lire
    else
       nb_paire_a_lire=npair
    end if
    if (rang==0) write(6,*)'nb of readed pairs',  nb_paire_a_lire
    if (ipotrep==1) then
       do lect_paire=1,nb_paire_a_lire
          if (npotentiel.gt.1) then
             read (lupotin,*) rf1,rf2,l
             roff1(l)=rf1;roff2(l)=rf2
             roff1(l)=roff1(l)*A2cm
             roff2(l)=roff2(l)*A2cm
          else
             read (lupotin,*) roff1(lect_paire),roff2(lect_paire)
             l=lect_paire
             roff1(l)=roff1(l)*A2cm
             roff2(l)=roff2(l)*A2cm
          end if
       end do
    end if
    read(lupotin,*)nptmax

    if (rang==0) write(6,*)'nptmax in the max number of points on grid  ',nptmax
!    if(nptmax.gt.ngrid) ngrid=nptmax
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

    do iti=1,ntyp
       allocate(embtyp(iti)%xg(0:nptmax)) 
       allocate(embtyp(iti)%feam(0:nptmax)) 
       allocate(SPembtyp(iti)%beam(0:nptmax)) 
       allocate(SPembtyp(iti)%ceam(0:nptmax)) 
       allocate(SPembtyp(iti)%deam(0:nptmax)) 

       if (lforcetabulate) then
          allocate(embtyp_d(iti)%xg(0:nptmax)) 
          allocate(embtyp_d(iti)%feam(0:nptmax)) 
          allocate(SPembtyp_d(iti)%beam(0:nptmax)) 
          allocate(SPembtyp_d(iti)%ceam(0:nptmax)) 
          allocate(SPembtyp_d(iti)%deam(0:nptmax))
       end if

       allocate(rhotyp(iti)%xd(0:nptmax)) 
       allocate(rhotyp(iti)%rho(0:nptmax)) 
       allocate(SPrhotyp(iti)%brho(0:nptmax)) 
       allocate(SPrhotyp(iti)%crho(0:nptmax)) 
       allocate(SPrhotyp(iti)%drho(0:nptmax)) 

       if (lforcetabulate) then
          allocate(rhotyp_d(iti)%xd(0:nptmax)) 
          allocate(rhotyp_d(iti)%rho(0:nptmax)) 
          allocate(SPrhotyp_d(iti)%brho(0:nptmax)) 
          allocate(SPrhotyp_d(iti)%crho(0:nptmax)) 
          allocate(SPrhotyp_d(iti)%drho(0:nptmax)) 
       end if
    end do
    do ipr=1,npair
       allocate(reppair(ipr)%xr(0:nptmax)) 
       allocate(reppair(ipr)%potr(0:nptmax)) 
       allocate(SPreppair(ipr)%bpotr(0:nptmax)) 
       allocate(SPreppair(ipr)%cpotr(0:nptmax)) 
       allocate(SPreppair(ipr)%dpotr(0:nptmax)) 


       if (lforcetabulate) then
          allocate(reppair_d(ipr)%xr(0:nptmax)) 
          allocate(reppair_d(ipr)%potr(0:nptmax)) 
          allocate(SPreppair_d(ipr)%bpotr(0:nptmax)) 
          allocate(SPreppair_d(ipr)%cpotr(0:nptmax)) 
          allocate(SPreppair_d(ipr)%dpotr(0:nptmax)) 
       end if
    end do

    select case(ipotentiel)
    case(10)

       !EMBD EAM PART
       do itir=1,ntypr
          iti=typtyp(itir)
          !lecture de Glue
          read(lupotin,*)n
          if (rang==0) write(6,*)'EAM',n,iti
          if(n.ne.iti)then
             write(6,*) rang,' incohérence eamtab iti'
             call arret_ndm
          end if

          ! write(6,*)'NPT NPTMAX',npt,nptmax
          read(lupotin,*)npt,embtyp(iti)%deltaEAM
          if (lforcetabulate) embtyp_d(iti)%deltaEAM=embtyp(iti)%deltaEAM
          if (rang==0) write(6,*)'EAM number of points in potin and the step',npt,embtyp(iti)%deltaEAM     
          if(npt.gt.nptmax)then
             write(6,*) rang,'nb de points de grille  EAM stop'
             call arret_ndm
          end if
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
          embtyp(iti)%xg(0)=0
          embtyp(iti)%feam(0)=0
          if (lforcetabulate) embtyp_d(iti)%feam(0)=0
          call checkround(embtyp(iti)%xg(1),embtyp(iti)%deltaEAM,lok)
            if (lok.eqv..false.) then
               write(6,*) 'error in glue grid  EAM stop', embtyp(iti)%xg(1), embtyp(iti)%xg(2),embtyp(iti)%deltaEAM
               write(6,*)'will not work for negative densities'
               call arret_ndm
            end if
          rhomin(iti)=min(rhomin(iti),embtyp(iti)%xg(1))
          rhomax(iti)=max(rhomax(iti),embtyp(iti)%xg(npt))
          !       write(6,*)'rhomin rhomax',rhomin,rhomax

          call cspline (nptmax,embtyp(iti)%xg,embtyp(iti)%feam,SPembtyp(iti)%beam,SPembtyp(iti)%ceam,SPembtyp(iti)%deam)
          if (lforcetabulate) then
             embtyp_d(iti)%xg=embtyp(iti)%xg
             call cspline (nptmax,embtyp_d(iti)%xg,embtyp_d(iti)%feam,SPembtyp_d(iti)%beam,&
                  &SPembtyp_d(iti)%ceam,SPembtyp_d(iti)%deam)
          end if
          !DENS PART
          !lecture de dens
          read(lupotin,*)n
          if (rang==0) write(6,*)'dens',n,iti
          if(n.ne.iti)then
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
          rhotyp(iti)%xd(0)= 0
          rhotyp(iti)%rho(0)= 0
          if(lforcetabulate) rhotyp_d(iti)%rho(0)=0
          call checkround(rhotyp(iti)%xd(1),rhotyp(iti)%deltaRHO,lok)
          if (lok.eqv..false.) then
             write(6,*) 'error in dens grid  EAM stop', rhotyp(iti)%xd(1), rhotyp(iti)%xd(2),rhotyp(iti)%deltaRHO
             write(6,*)'will not work for negative densities'
             call arret_ndm
          end if
          
          call cspline (nptmax,rhotyp(iti)%xd,rhotyp(iti)%rho,  SPrhotyp(iti)%brho,  SPrhotyp(iti)%crho,  SPrhotyp(iti)%drho  )
          if (lforcetabulate) then
             rhotyp_d(iti)%xd=rhotyp(iti)%xd
             call cspline (nptmax,rhotyp_d(iti)%xd,rhotyp_d(iti)%rho,SPrhotyp_d(iti)%brho,SPrhotyp_d(iti)%crho,SPrhotyp_d(iti)%drho)
          end if
       end do

       !PAIR PART 

       do ipair=1,nb_paire_a_lire
          if (npotentiel.GT.1) then
             read(lupotin,*)it1,it2
             if (rang==0)write(6,*)'pair it1 it2',it1,it2
             l=ipo(it1,it2)
             if (typ_pot_pair(l).ne.0)then
                write(6,*)'pair deja lue typ_pot_pair,ipotentiel',l,typ_pot_pair(l),ipotentiel
             end if
             typ_pot_pair(l)=ipotentiel
             ipr=l
          else
             read(lupotin,*)n
             if(n.ne.ipair)then
                write(6,*) rang, ' ordre de lecture de EAM rep stop'
                call arret_ndm
             end if
             ipr=n
          end if


          if (rang==0) write(6,*)'paire eam ; paire complete',n,ipr
          if (rang==0) write(6,*)'rep'

          read(lupotin,*)npt,reppair(ipr)%deltaREP
!          reppair(ipr)%deltaREP=reppair(ipr)%deltaREP*1.00000000000000000000001 !+low_limit ! *1.00000000000000000000000000000010.9999999
          if (lforcetabulate) reppair_d(ipr)%deltaREP=reppair(ipr)%deltaREP

          if(npt.gt.nptmax)then
             write(6,*) rang,'nb de points de grille  EAM stop'
             call arret_ndm
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
          reppair(ipr)%xr(0)= 0  
          reppair(ipr)%potr(0)= 0
          if (lforcetabulate) reppair_d(ipr)%potr(0)=0
          call checkround(reppair(iti)%xr(1),reppair(iti)%deltaREP,lok)
          if (lok.eqv..false.) then
             write(6,*) 'error in glue grid  EAM stop', reppair(iti)%xr(1),reppair(iti)%xr(2),reppair(iti)%deltaREP

             call arret_ndm
          end if

          call cspline (nptmax,reppair(ipr)%xr,reppair(ipr)%potr,SPreppair(ipr)%bpotr,& 
               SPreppair(ipr)%cpotr,SPreppair(ipr)%dpotr)
          if (lforcetabulate) then
             reppair_d(ipr)%xr=reppair(ipr)%xr
             call cspline (nptmax,reppair_d(ipr)%xr,reppair_d(ipr)%potr,SPreppair_d(ipr)%bpotr,&
                  SPreppair_d(ipr)%cpotr,SPreppair_d(ipr)%dpotr)
          end if

       end do
    case(16)
       rhomin(:)=0
       do itir=1,ntypr
          iti=typtyp(itir)
          read(lupotin,*)n
          if (rang==0) write(6,*)'EAM',n,iti
          if(n.ne.iti)then
             write(6,*) rang,' incohérence eamtab iti'
             call arret_ndm
          end if
          read(lupotin,*)crg%G(iti),crg%n(iti),crg%densmax(iti)
          rhomax(iti)=crg%densmax(iti)
          rhomin(iti)=0
       end do
       do ipair=1,nb_paire_a_lire
          if (npotentiel.GT.1) then
             read(lupotin,*)it1,it2
             l=ipo(it1,it2)
             if (typ_pot_pair(l).ne.0)then
                write(6,*)'pair deja lue typ_pot_pair,ipotentiel',l,typ_pot_pair(l),ipotentiel
             end if
             typ_pot_pair(l)=ipotentiel
             ipr=l
          else
             read(lupotin,*)it1,it2
             n=ipo(it1,it2)
             !             read(lupotin,*)n
             if (rang==0) write(6,*)'paire ',n,ipair
             if(n.ne.ipair)then
                write(6,*) rang, ' ordre de lecture de EAM rep stop'
                call arret_ndm
             end if
             ipr=ipair
             typ_pot_pair(ipr)=ipotentiel
          end if

          read(lupotin,*)crg%D(ipr),crg%gam(ipr),crg%A(ipr),crg%rho(ipr),crg%C(ipr),crg%R0(ipr)
       end do
!       crg%D(:)=crg%D(:)*ev2erg
!       crg%A(:)=crg%A(:)*ev2erg
!       crg%C(:)=crg%C(:)*ev2erg*1d48
       
    end select
    close(lupotin)
    !if(allocated* (typ_and_pot).eqv..false.), i.e. if npotentiel==1 
    if(allocated (typ_and_pot).eqv..false.) then
       allocate (typ_and_pot(ntyp,npotmax))
       typ_and_pot(:,:)=.false.
       typ_and_pot(1:ntyp,ipotentiel)=.true.
    end if
    rumax=max(rue,rumax)

    alpha=alpha*1d8
    return

  end subroutine inputeam


  subroutine extrapolateRepCRG(rcm,l,Erep)
    integer,intent(in)::l
    real(double),intent(in)::rcm

    real(double),intent(out)::Erep

    real(double)::r,fhi_M,fhi_B,damp
    real(double):: rc,wc,rd

    rc=0.7; wc=10.0
    r=rcm/A2cm

    fhi_M=crg%D(l)*(exp(-2*crg%gam(l)*(r-crg%r0(l)))-2*exp(-1*crg%gam(l)*(r-crg%r0(l))))

    rd=(r-rc)*wc
    call calerf(rd,damp,0)
    damp=0.5*(1+erf(rd))
!    damp=1
    if (crg%rho(l).ne.0) then
       fhi_B=crg%A(l)*exp(-1*r/crg%rho(l))
    else
       fhi_B=0
    end if
    fhi_B=fhi_B-damp*crg%C(l)/r**6

    Erep=ev2erg*(fhi_B+fhi_M)

  end subroutine extrapolateRepCRG

  subroutine extrapolateRhoCRG(rcm,iti,rho)
    integer,intent(in)::iti
    real(double),intent(in)::rcm

    real(double),intent(out)::rho

    real(double)::r,damp
    real(double):: rc,wc,rd

    rc=1.5; wc=20.0
    r=rcm/A2cm

    rd=(r-rc)*wc
    call calerf(rd,damp,0)
    damp=0.5*(1+erf(rd))
    if (crg%n(iti).ne.0) then
       rho=crg%n(iti)*(damp/r**8+(1-damp)/rc**8)
    else
       rho=0
    end if
!    write(6,*)r,rho,damp
  end subroutine extrapolateRhoCRG

  
  !-------------------------------------------

  subroutine extrapolateRho(density,SPdensity, r2, rho)
    ! calculate electronic density at distance sqrt(r)
    ! or its first and second derivatives

    USE gen_com_m, ONLY:  ev2erg
    implicit none

    type(DensityT), intent(in) :: density 
    type(DensityTsp), intent(in) :: SPdensity 
    real(double), intent(in) :: r2
    real(double), intent(out), optional :: rho
    !local
    integer:: kr 
    real(double) :: xmax,r,drk
    
    kr=0
    xmax=density%xd(nptmax)
    r=sqrt(r2)/A2cm

    if(r.gt.xmax) then
       Rho=0.0
    else
       kr=Int(r/density%deltaRHO)
       drk=r-kr*density%deltaRHO
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
       k=Int((rho-eam%deltaEAM/1d10)/eam%deltaEAM)
!       write(6,*) k,rho,rho/eam%deltaEAM
       drk=rho -k*eam%deltaEAM
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
    real(double), intent(in) :: r2
    real(double), intent(out), optional :: Erep


    !-----------------------------------

    integer :: k
    real(double) :: xmax,drk,r

    xmax=rep%xr(nptmax)

    r=r2/A2cm

    if(r.gt.xmax) then
       Erep=0.0
       return
    else

       k=Int(r/rep%deltaREP)
       drk=r-k*rep%deltaREP
       Erep = ev2erg*(rep%potr(k) +drk*(SPrep%bpotr(k) +drk*(SPrep%cpotr(k) +drk*SPrep%dpotr(k))))

    end if



    RETURN
    !-----------------------------------
  end subroutine extrapolateRep

  subroutine checkround(x1,x2,lok)
    real(double),intent(in)::x1,x2
    logical,intent(out)::lok
    real(double)::x1r,x2r
    integer::n
    x1r=nint(x1*1d6)/1d6
    x2r=nint(x2*1d6)/1d6
    if (x1r==x2r) then
       lok=.true.
    else
       lok=.false.
    end if
  end subroutine checkround
end module eam
     
