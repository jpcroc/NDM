module posana
  USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m
  use setnoxsimple_mod,only:setnoxsimple
  USE sic
  USE var_pot, ONLY:ty,ntyp,rumax
  !  USE period_mod,only: period
  USE recips_mod,only: recips
  use cryst_to_cart_mod,only:cryst_to_cart
  USE gen_com_m, ONLY: fnam,rang,lperiod,pi,dmtype,lenfnam,iteration,timel,ivisu,unit6P
  USE atomconfig,only:atom_config
  use cellconfig,only: cell_config,caltabtC
  use boxconfig,only:box_config,periodbox
  use rasmolT_mod,only:rasmolT
  USE constrconf_mod,only:gin2ndm,read_cin
  use vect_dist_mod,only:vect_dist
  use neb_module, only : npath
  use newunit_mod,only:newunit
  use read_val,only:ipbc
  implicit none

  CHARACTER(len=89) :: fnamcr,namecr,fnamperfdef
  logical :: lcomp, & ! comparaison ou non avec un cristal de dÃ©part
       ldecal, & ! decalage en tre boite cr et boite ana
       ldesord, &      ! vielle variable historique
       lnbvois, &      ! analyse des nombres de voisins
       lc15, &      ! lacunes et int ensembles
       lvac, &           ! analyse en lacune
       ldetdec,&          !determination du decalage
       lrescale ,&          ! rescale des posistions de dÃ©part sur la boite d'arrivÃ©e
       lpstruct ,&          
       lallint ,&          ! dans ws : si true dumbbal=2 ints ; si false dumbbal =1 int
       lpdep ,&          
       ldeptest ,&          
       lws ,&          ! Wigner-Seitz pour INT et VAC
       lperfdef ,&          ! defectr eference frame for perfect config : defectif !
       lpdef, &
       ldefcat, &        ! defauts sur les cations seulement
       ldepla, &        ! nouvelle analyse des déplacements à partir de atana0
       lsubc, &        ! analyse en sous cascade BLOB
       distordflag    ! analyse des angles dans le cristal si flag==.true.
  integer::         idistord,iprtnvi        ! analyse des diff angulaires
  integer::immcr
  !   integer :: imcr ! nb d'atomes dans le cristal de reference
  integer::imdum,ivisuana
  !    real(double)::rdum
  integer::ivuana
  real(double)::pstmax(20)
  !    real(double),allocatable:: xpcr(:,:)
  !    integer, allocatable :: itypcr(:)
  integer:: ndvblob,ndvmin ! voisins blob pour SC
  real(double)::rdv ! distance entre defauts pour SC
  integer, allocatable:: indws(:),indatsit(:,:),natsit(:),indint(:),indvac(:),indas(:) ! indices des atomes deplaces et plottes
  ! **************************************************************
  type(atom_config)::atcr,atrefdep,atperfdef
  type(cell_config)::celcr,celrefdep,celperfdef
  type(box_config)::boxcr,boxrefdep,boxperfdef
  integer::igencr
  integer,allocatable::na(:)
  real(double) :: plmin(3),plmax(3) ! bords de la portion afichÃ©e de la boite
  real(double) :: tvac ,tint,deltx,delty,deltz ! distance pour les lacunes et les int
  real(double) :: tdep ! seuil de deplacement
  real(double),allocatable::rc(:)
contains

  subroutine initanapos(atana0,celana0,boxana0)
    USE T_kind_param_m

    ! **************************************************************
    implicit none

    class(atom_config),intent(in)::atana0
    class(cell_config),intent(in)::celana0
    class(box_config),intent(in)::boxana0
    integer :: i, ic,nbvoisparf(20,20)
    real(double)  :: decal(3)


    integer :: idecal    ! alignement des posistions sur l'atome idecal
    integer, save:: icall=0
    integer::unitana
    real(double)::rclu(20)

    real (double) :: plmin1,plmin2,plmin3, plmax1,plmax2,plmax3 ! bord de plot lu dans la namelist
    namelist /analyse/ldecal,ldesord,idistord,idecal,tdep,lvac,tvac,tint,lcomp,plmin1,igencr,fnamcr, &
         plmin2,plmin3, plmax1,plmax2,plmax3,ldetdec,lrescale,lpstruct,lpdef,lpdep,ldeptest,namecr, &
         ldefcat,rclu, lnbvois,nbvoisparf,pstmax,iprtnvi,lc15,lws,deltx,delty,deltz,lallint,ndvblob,ndvmin,rdv,&
         &lsubc,ivisuana,ldepla,immcr,lperfdef

    allocate(rc(ntyp))
    allocate (na(ntyp))

    !
    !   set default values for variables in namelist
    !
    !-----------------------------------------------
    namecr='ZZ'
    lperfdef=.false.
    immcr=-1
    ivisuana=-1
    igencr=-1
    rclu(:)=2.8
    pstmax(:)=0.
    lcomp=.false.
    !   iprtnvi=.false.
    iprtnvi=0      
    lws=.false.
    ldetdec=.false.
    idecal= -1
    ldesord=.false.
    lnbvois=.false.
    idistord=0
    ldeptest=.true.
    lpdep=.false.
    lpstruct=.true.
    lpdef=.true.
    lvac=.false.
    lc15=.false.
    tdep=2.0
    tvac=1.1
    tint=1.1
    plmin1=0;plmin2=0.;plmin3=0.
    plmax1=0.;plmax2=0.;plmax3=0.
    lrescale=.true.
    ldefcat=.false.
    deltx=0. ; delty=0.0; deltz=0.0
    lallint=.false.
    nbvoisparf(:,:)=0
    lsubc=.false.
    ndvblob=0
    ndvmin=0
    ldepla=.false.
    if (rang==0)write(6,*)'*** analyse du crystal'
    call newunit(unitana)
    open(unitana, file='analyse.in')
    read(unitana,nml=analyse)

    if (ldepla) then
       if (rang==0) write(6,*)'detection des deplacements par rapport a atana0'
       !       if (.not.present(atana0)) then
       !          if (rang==0) write(6,*)'atana0 pas definie STOP'
       !          call arret_ndm
       !       end if
       atrefdep=atana0
       celrefdep=celana0
       boxrefdep=boxana0
    end if

    if(.not.(lnbvois).and.(ldesord))lnbvois=.true.
    if(lws) then
       lcomp=.true.
       !       lvac=.true.

    end if
    if (lcomp) then
       if (immcr==-1) then
          if (rang==0) write(6,*)'IMMcr must ne specified'
          call arret_ndm
       end if
       if(ivisuana==-1) ivisuana=ivisu
    end if
    if (ivisuana==4) ivisuana=40
    if (ivisuana==6) ivisuana=60

    if ((idistord.gt.0).and.(.not.lperiod)) then
       write(6,*) 'distord seulement avec lperiod =true.'
       call arret_ndm
    end if
    if ((idistord.gt.0).and.(.not.lperiod)) then
       write(6,*) 'distord ne fonctionne pas '
       call arret_ndm
    end if
    if ((idistord.ge.3).and.(maxval(pstmax)==0))then
       write(6,*) 'isdistrod=3 preciser pstmax'
       call arret_ndm
    end if
    if(namecr=='ZZ') then
       namecr=fnam(1:lenfnam)
    end if
    ! write(6,*)'distordflag',distordflag
    rc(1:ntyp)=rclu(1:ntyp)*1.0d-8

    ! pourquoi plotpart ici ?
    plmin(1)=plmin1; plmax(1)=plmax1
    plmin(2)=plmin2; plmax(2)=plmax2
    plmin(3)=plmin3; plmax(3)=plmax3
    plmin=plmin*1.0d-8 ; plmax=plmax*1.0d-8
!!$    if (any(plmin.ne.0.).or.any(plmax.ne.0.)) then
!!$       call plotpart(atcf,plmin)
!!$    end if

    if (lcomp) then
       select case (igencr)
       case(1)
          fnamcr=namecr(1:len(namecr))//'crcin'
          call read_cin(boxcr,1,atcr,immcr,fnamcr)
          !             atcf%im_glob=atcr%im
          call setnoxsimple (atcr,boxcr,celcr,rumax)
       case(0)
          fnamcr=trim(namecr)//'.crgin'
          if (rang==0)write(6,*)'fnamcr ',len(fnamcr),fnamcr
          call gin2ndm(atcr,celcr,boxcr,fnamcr,rumax,lrepartition=.false.,immread=immcr,lconstrsimple=.true.)
          call celana0%copy(celcr,boxana0)
          !            call setnoxsimple (atcr,boxcr,celcr,rumax)
          fnamperfdef=trim(fnam)//'.perfdef.gin'
          if (rang==0)write(6,*)'fnamperfdef ',len(fnamperfdef),fnamperfdef
          call gin2ndm(atperfdef,celperfdef,boxperfdef,fnamperfdef,rumax,lrepartition=.false.,immread=immcr,lconstrsimple=.true.)
          !             call setnoxsimple (atperfdef,boxperfdef,celperfdef,rumax)
          call celana0%copy(celperfdef,boxana0)
          call caltabtC(celperfdef,atperfdef,lperiod,boxcr,lchktrav=.false.)
       case default
          write(6,*)'set igencr to 1 or 0 for .crcin or .crgin file respectively'
          call arret_ndm
       end select
       call caltabtc(celcr,atcr,lperiod,boxcr,lchktrav=.false.)
       !       write(6,*)' celana0 celcr celperfdef ',celana0%nox,celcr%nox,celperfdef%nox
       if (lws) then
          if (rang==0)write(6,*)'analyse de Wigner-Seitz'
       else
          if(ldeptest) then
             if (rang==0)write(6,'(A,F6.1)') 'seuil deplacement pour detection de defauts ',tdep
             tvac=tdep
          end if
          if(lpdep) then
             if (rang==0)write(6,'(A,F6.1)') 'ecriture des deplacés ',tdep
          end if
          tint=tvac
          if (rang==0)write(6,'(A,F6.1)') 'seuil lacune ',tvac 
          if (rang==0)write(6,'(A,F6.1)') 'seuil interstitiel ',tint 
          tdep=tdep*1.0d-8
          tvac=(tvac*1.0d-8)
          tint=(tint*1.0d-8)
       end if
       


!!$    if (any(boxcf%at.ne.boxcr%at) )then
!!$       write(6,*)'boxf <> boxcr'
!!$       write(6,*)'atcf',boxcf%at
!!$       write(6,*)'atcr',boxcr%at
!!$!       call arret_ndm
!!$    end if

       if (ldecal) then
          if (idecal.gt.0) then
             if (idecal.gt.atana0%im) then
                write(6,*)'idecal >atcf%im ; stop'
                call arret_ndm
             end if
             do ic=1,3
                decal(ic)=atcr%xp(ic,idecal)-atana0%xp(ic,idecal)
             end do
             do i=1,atcr%im
                do ic=1,3
                   atcr%xp(ic,i)=atcr%xp(ic,i)-decal(ic)
                end do
             end do
          else
             do i=1,atcr%im
                atcr%xp(1,i)=atcr%xp(1,i)-deltx*1d-8
                atcr%xp(2,i)=atcr%xp(2,i)-delty*1d-8
                atcr%xp(3,i)=atcr%xp(3,i)-deltz*1d-8
             end do
          end if

       end if
       !       if (lperiod) then
       !          call periodbox (boxcf,atcf)
       call periodbox(boxcr,atcr)
       !       end if
       call caltabtc(celcr,atcr,lperiod,boxcr,lchktrav=.false.)
    else
       if (ldepla) then
          if (rang==0)write(6,'(A,F6.1)') 'test depla ',tdep
          tdep=tdep*1d-8
       end if
    end if
    close (unitana)
  end subroutine initanapos

  subroutine anapos(atana,celana,boxana,itapp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m

    ! **************************************************************
    implicit none

    class(atom_config),intent(in)::atana
    class(cell_config),intent(in)::celana
    class(box_config),intent(in)::boxana
    integer,optional::itapp
    integer::itap
    integer :: i,nbvoisparf(20,20)




    integer, save:: icall=0
    type(atom_config)::atcf
    type(cell_config)::celcf
    type(box_config)::boxcf
    real(double)::rclu(20)
    real(double),dimension(:,:),allocatable::vectdep
    real(double),dimension(:),allocatable::valdepla
    integer,dimension(:),allocatable::inddep

    integer::ndep
    real (double) :: plmin1,plmin2,plmin3, plmax1,plmax2,plmax3 ! bord de plot lu dans la namelist
!!$    namelist /analyse/ldecal,ldesord,idistord,idecal,tdep,lvac,tvac,tint,lcomp,plmin1,igencr,fnamcr, &
!!$         plmin2,plmin3, plmax1,plmax2,plmax3,ldetdec,lrescale,lpstruct,lpdef,lpdep,ldeptest,namecr, &
!!$         ldefcat,rclu, lnbvois,nbvoisparf,pstmax,iprtnvi,lc15,lws,deltx,delty,deltz,lallint,ndvblob,ndvmin,rdv,&
!!$         &lsubc,ivisuana

    atcf=atana
    boxcf=boxana
    celcf=celana
    allocate (na(ntyp))

    do i=1,atcf%im
       na(atcf%ityp(i))= na(atcf%ityp(i))+1
    end do
    if (present(itapp)) then
       itap=itap
    else
       itap=0
    end if
    !
    !   set default values for variables in namelist
    !
    !-----------------------------------------------
    if (lcomp)then
       if(ldetdec) then
          if (atcf%im.ne.atcr%im) then
             write(6,*)'detdec impossible'
             call arret_ndm
          end if
          do i=1,atcf%im
             write(490,'(I8,3E16.5)')i,atcr%xp(1,i)-atcf%xp(1,i),atcr%xp(2,i)-atcf%xp(2,i),atcr%xp(3,i)-atcf%xp(3,i)
             write(487,'(I8,3E16.5)')i,atcf%xp(1,i),atcr%xp(1,i),atcr%xp(1,i)-atcf%xp(1,i)
             write(488,'(I8,3E16.5)')i,atcf%xp(2,i),atcr%xp(2,i),atcr%xp(2,i)-atcf%xp(2,i)
             write(489,'(I8,3E16.5)')i,atcf%xp(3,i),atcr%xp(3,i),atcr%xp(3,i)-atcf%xp(3,i)
             deltx=deltx+(atcr%xp(1,i)-atcf%xp(1,i))/atcf%im
             delty=delty+(atcr%xp(2,i)-atcf%xp(2,i))/atcf%im
             deltz=deltz+(atcr%xp(3,i)-atcf%xp(3,i))/atcf%im
          end do
          write(6,*)deltx*1d8,delty*1d8,deltz*1d8
          call arret_ndm
       end if

       !comparaison avec cristal
       call ws (atcf,celcf,boxcf,atcr,celcr,boxcr,lws)  ! WS and detcec !
!!$       if (lws) then
!!$
!!$       else
!!$          call depcr (atcf,celcf,boxcf,atcr,celcr,boxcr)  
!!$       end if
    end if
    if (ldepla) then
       call depladet(atcf,celcf,boxcf,atrefdep,celrefdep,boxrefdep,tdep,ndep,inddep,vectdep,valdepla)
    end if
    !analyse des voisins (sans comparaison avec le cristal de référence)
    if (lnbvois)  call nbvois(atcf,boxcf,celcf,icall,nbvoisparf,itapp)

    !    close(175)     

    !  if (dmtype==6)   call arret_ndm
    !    write(6,*)'fin compcr'

    deallocate(na)

    return
  end subroutine anapos
  !**********************************************************

  subroutine depladet(atc,celc,boxc,atr,celr,boxr,tdep,ndep,inddep,vectdepla,valdepla)
    !cette routine ne fonctionne :
    ! 1/ qu'en séquentiel
    ! 2/ que pour des configurations atcf atrefdep directement comparables (mêmes atomes même ordre)
    ! il faudra programmer une détection des déplacements à la caldepla avec ax pour les cascades en para
    ! cette routine remplace depcr qui supposait que la même configuration servait pour les défauts et les dépalcements
    USE T_kind_param_m

    !------------pnp-----------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------

    class(atom_config)::atc,atr
    type(cell_config)::celc,celr
    type(box_config)::boxc,boxr
    real(double),intent(in)::tdep
    real(double),dimension(:,:),allocatable::vectdepla
    real(double),dimension(:),allocatable::valdepla
    integer, dimension(:), allocatable:: inddep(:)
    integer::ndep
    !Local variables
    integer :: i,idp,im


    real(double) :: tdep2
    real(double) :: a1,a2,a3,c1,c2,c3,r2
    real(double), dimension(1,3) :: cv

    real(double) :: r2min

    real(double) :: c3p,c2p,c1p,c1abs,c2abs,c3abs,r,XJI(3)
    integer:: koo,i2,i1,ncelvois,ko1,id


    if(atc%im.ne.atr%im) then
       write(6,*)'not depladet works for confs with equal number of atoms',atc%im,atr%im
       call arret_ndm
    end if
    im=min(atc%im,atr%im)

    ncelvois = min(celc%noxyz,27)-1
    !    atc%lgul=.false.
    !    atr%lgul=.false.

    tdep2=tdep*tdep
    ndep=0
    plmin=1000.0 ; plmax=-1000.0




    !    allocate(inddep(im))
    !    atc%lgul=.false. ! lgul = true pour les déplacés
    !    atr%lgul=.false.
    call cryst_to_cart (atc%im, atc%xp, boxc%bg, -1)    !cart vers cryst
    call cryst_to_cart (atr%im, atr%xp, boxr%bg, -1)    !cart vers cryst

    do i=1,im
       !     if(ityp(i)==2) cycle
       c1 = atc%xp(1,i)-atr%xp(1,i)
       c2 = atc%xp(2,i)-atr%xp(2,i)
       c3 = atc%xp(3,i)-atr%xp(3,i)
       if (c1>0.5) c1 = c1-1.
       if (c1<(-0.5)) c1 = c1+1.
       if (c2>0.5) c2 = c2-1.
       if (c2<(-0.5)) c2 = c2+1.
       if (c3>0.5) c3 = c3-1.
       if (c3<(-0.5)) c3 = c3+1.
       cv(1,1) = c1
       cv(1,2) = c2
       cv(1,3) = c3
       call cryst_to_cart (1, cv, boxc%at, 1) !cryst vers cart sur cv
       r2 = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)
       if(r2.gt.tdep2) then
          ndep=ndep+1
          !          inddep(ndep)=i
          !          atc%lgul(i)=.true.
          !          atr%lgul(i)=.true.
       end if

    end do  !boucle i

    if (allocated(vectdepla))deallocate(vectdepla)
    if (allocated(inddep))deallocate(inddep)
    if (allocated(valdepla))deallocate(valdepla)
    allocate(inddep(ndep))
    allocate(valdepla(ndep))
    allocate(vectdepla(3,ndep))
    id=0

    do i=1,im
       c1 = atc%xp(1,i)-atr%xp(1,i)
       c2 = atc%xp(2,i)-atr%xp(2,i)
       c3 = atc%xp(3,i)-atr%xp(3,i)
       if (c1>0.5) c1 = c1-1.
       if (c1<(-0.5)) c1 = c1+1.
       if (c2>0.5) c2 = c2-1.
       if (c2<(-0.5)) c2 = c2+1.
       if (c3>0.5) c3 = c3-1.
       if (c3<(-0.5)) c3 = c3+1.
       cv(1,1) = c1
       cv(1,2) = c2
       cv(1,3) = c3
       call cryst_to_cart (1, cv, boxc%at, 1) !cryst vers cart sur cv
       r2 = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)
       if(r2.gt.tdep2) then
          id=id+1
          inddep(id)=i
          vectdepla(:,id)=cv(1,:)
          valdepla(id)=sqrt(r2)
       end if

    end do  !boucle i

    call cryst_to_cart (atc%im, atc%xp, boxc%at, 1)    !cart vers cryst
    call cryst_to_cart (atr%im, atr%xp, boxr%at, 1)    !cart vers cryst

  end subroutine depladet


  subroutine nbvois(atcf,boxcf,celcf,icall,nbvoisparf,itapp)

    USE T_kind_param_m
    implicit none
    !variables transmises
    type(atom_config)::atcf
    type(box_config)::boxcf
    type(cell_config)::celcf
    integer , intent(in)  :: nbvoisparf(20,20)
    integer :: icall,itapp

    type(box_config)::boxplt

    !Variables locales
    integer :: koo,iti,i,ko1,i1,i2,lenfn2,lenfn,j
    integer :: natvi(0:20),natvityp(0:20,20),iatvi,latvi
    real(double):: c1,c2,c3,r2,xp1,xp2,xp3
    character*9 :: extension
    character(len=2) :: extension2
    integer, save:: lurasmol
    integer :: maxvois ,nana,itj,i5,iwr
    integer,allocatable, save :: nvi(:),nvityp(:,:),ivois(:,:)

    real(double) :: dx(3,20)
    real(double)::cv(1,3)
    logical::lvoisOK,permut
    real(double)::pst,psta,pstp,mtheta(200),costheta,moytet(200,20),thetajik,tetreg,temptri
    integer,allocatable::iatd(:)
    integer::natd,k11,j11,nteta,ntetmax,natdtyp(20)
    real(double)::rdist,XJI(3),rcm
    !    type(atom_config)::atplt
    class(atom_config),allocatable::atplt,atpltpp
    real(double),allocatable:: rccar(:)

    integer::naux
    real(double),allocatable::vaux(:,:)
    character (len=80),allocatable::charaux(:)
    logical::linter
    CHARACTER(len=89) :: namedes,namedis
    natd=0
    if (idistord.ge.3)allocate(iatd(atcf%im))
    lvoisOK=.true.
    atcf%lgul=.false.

    allocate (rccar(ntyp))
    rccar(1:ntyp)=rc(1:ntyp)**2

    rcm=maxval(rc(:))
    write(6,*)'rcm ', rcm

    call setnoxsimple(atcf,boxcf,celcf,rcm)
    allocate(nvi(atcf%im))
    allocate(nvityp(atcf%im,ntyp))
    natvi(:)=0
    !  itapp=it

    !calcul en deux temps
    !calcul du nombre de voisins par atome

    call caltabtC(celcf,atcf,lperiod,boxcf,lchktrav=.false.)



    natvityp(:,:)=0
    ntetmax=-1
    LOOPAT1I:    do i = 1,atcf%im
       if(ldefcat.and.atcf%ityp(i)==2) cycle 
       pst=0;psta=0;pstp=0.
       nvi(i)=0
       nvityp(i,:)=0
       koo = atcf%ielat(i)
       iti=atcf%ityp(i)
       do i1 = 0, celcf%ncelvois(koo)
          ko1 = celcf%ncel(koo,i1)
          !       write(6,*)'koo,ko1',koo,ko1
          do i2 = 1, celcf%nato(ko1)
             j = celcf%atincel(i2,ko1)
             !         write(6,*)'j',j
             !         write(6,*)i,xp(:,i)
             !         write(6,*)j,xp(:,j)
             if (i==j) cycle

             call vect_dist(atcf,celcf,boxcf,i,j,XJI,i1,lperiod,rcm,rdist,linter)

             !        write(6,*)c1,c2,c3,r2,rc(iti)

             if (linter) then 
                if (rdist<rc(iti))then
                   nvi(i)=nvi(i)+1
                   nvityp(i,atcf%ityp(j))=nvityp(i,atcf%ityp(j))+1
                   dx(:,nvi(i))=XJI(:)
                end if
             end if

          end do
       end do
       natvi(nvi(i))=natvi(nvi(i))+1
       natvityp(nvi(i),iti)=natvityp(nvi(i),iti)+1

       if (idistord.gt.0) then
          do itj=1,ntyp
             if (nvityp(i,itj).ne.nbvoisparf(iti,itj))then
                lvoisok=.false.
                !                write(6,*)'pour i de type mauvais nb de voisins de type j'
                !                write(6,*)i,iti, itj, nvityp(i,itj),nbvoisparf(iti,itj)
             end if
          end do


          if (lvoisok.EQV..false.) then
             lvoisok=.true.
             cycle ! disordered atoms cannot be distorted
          end if

          nteta=0
          loopvoisi1 :do j11=1,nvi(i)-1
             loopvoisi2 :   do k11= j11+1, nvi(i)
                nteta=nteta+1
                if (ntetmax.lt.nteta) ntetmax=nteta

                costheta = (dx(1,j11)*dx(1,k11)+dx(2,j11)*dx(2,k11)+dx(3,j11)*dx(3,k11))/ &
                     (sqrt(dx(1,j11)**2+dx(2,j11)**2+dx(3,j11)**2)*sqrt(dx(1,k11)**2+ &
                     dx(2,k11)**2+dx(3,k11)**2))
                !		write(6,*)i,j11,k11,indice_mtheta1 , costheta OK

                if (costheta>-1.0000001 .and. costheta<-0.9999999) then
                   thetajik=pi
                else
                   thetajik = dacos(costheta)
                end if

                !		write(6,*) 'thetajik', thetajik OK

                mtheta(nteta)=thetajik

             end do loopvoisi2
          end do loopvoisi1
          !************************classement par ordre croissant des angles********************
          permut =.true.

          do while (permut.EQV..true.)
             permut=.false.
             do i5=1, nteta-1
                if (mtheta(i5)>mtheta(i5+1)) then
                   temptri=mtheta(i5)
                   mtheta(i5)=mtheta(i5+1)
                   mtheta(i5+1)=temptri
                   permut=.true.
                end if
             end do
          end do



          select case (idistord)
          case(1)    ! cristal = 0K
             do i5=1,nteta
                moytet(i5,iti)=moytet(i5,iti)+mtheta(i5)/float(na(iti))
             end do


          case(2,3)

             iwr=350+iti
             open(unit=iwr)
             do i5=1,nteta
                read(iwr,*)moytet(i5,iti)
                tetreg=pi/nteta
                psta=psta+(mtheta(i5)-moytet(i5,iti))**2
                pstp=pstp+(moytet(i5,iti)-tetreg)**2
                write(666,*)i,mtheta(i5),moytet(i5,iti),tetreg
             end do
             close(iwr)
             pst=sqrt(psta)/sqrt(pstp)
             if( idistord==2) then
                write(380,*)i,iti,pst
             else
                if (pst.gt.pstmax(iti)) then
                   natd=natd+1
                   iatd(natd)=i
                   atcf%lgul(i)=.true.
                   natdtyp(iti)=natdtyp(iti)+1
                end if
             end if

          case default
          end select


       end if
    end do LOOPAT1I


    select case (idistord)
    case(1)    ! cristal = 0K
       do iti=1,ntyp
          if (na(iti).ne.0) then
             iwr=350+iti
             do i5=1,ntetmax
                write(iwr,*)moytet(i5,iti)
             end do
          end if
       end do
    case(3)

       !    §REMPLACER PAR UN APPEL A RASMOLT
       !       open(92,file="distordimag.mol")
       if (natd .NE. 0) then
          call atcf%fab(atplt,lback=.false.)
          atcf%lgul=.false.
          write(6,*)natd,' atomes distordus'
          do iti=1,ntyp
             if (natdtyp(iti).ne.0)write(6,*)natdtyp(iti),' atomes distordus de type ',iti
          end do
          if (itapp.gt.0) then
             lenfn = 9
             write(extension,'(i9.9)') itapp
             namedis='distort.'//extension
          else
             namedis='distort.'
          end if
          if (any(plmin.ne.0.).or.any(plmax.ne.0.)) then
             call atcf%deftype(atpltpp)
             call plotpart(atplt,plmin,plmax,atpltpp,boxplt)
             call rasmolT(atpltpp,boxplt,namefr=namedis,latcomp=.true.,ivisumol=ivisuana)
          else
             boxplt=boxcf
             call rasmolT(atplt,boxplt,namefr=namedis,latcomp=.true.,ivisumol=ivisuana)
          end if

       else
          write(6,*) 'Pas d atomes distordus'
       end if

    end select
    !ouverture des fichiers NVI

    atcf%lgul=.false.
    do iatvi=1,20
       if(natvi(iatvi).ne.0) then
          write(6,*)
          write(6,*)'Nb d_at. avec',iatvi,'vois. =',natvi(iatvi)!
          do iti=1,20
             if (natvityp(iatvi,iti).ne.0) write(6,*)'Nb d_at. de type',iti,' avec',iatvi,'vois. =',natvityp((iatvi),iti)
          end do
       end if
    end do
    select case (iprtnvi)
    case(1)
       call atcf%deftype(atplt)
       atcf%lgul=.false.
       do iatvi=1,20
          if(natvi(iatvi).ne.0) then
             !          write(6,*)it,'Nb d_at. avec',iatvi,'vois. =',natvi(iatvi)!
             lenfn2 = 2
             write(extension2,'(i2.2)') iatvi
             if (itapp.gt.0) then
                lenfn = 9
                write(extension,'(i9.9)') itapp
                namedes='disorder.'//extension//'.NVI.'//extension2(1:lenfn2)
             else
                namedes='disorder.'//'.NVI.'//extension2(1:lenfn2)
             end if
             do i=1,atcf%im
                if (nvi(i)==iatvi) atcf%lgul(i)=.true.
             end do
             call atcf%fab(atplt,lback=.false.)
             atcf%lgul=.false.

             if (any(plmin.ne.0.).or.any(plmax.ne.0.)) then
                call atcf%deftype(atpltpp)
                call plotpart(atplt,plmin,plmax,atpltpp,boxplt)
                call rasmolT(atpltpp,boxplt,namefr=namedes,latcomp=.true.,ivisumol=ivisuana)
             else
                boxplt=boxcf
                call rasmolT(atplt,boxplt,namefr=namedes,latcomp=.true.,ivisumol=ivisuana)
             end if
             atcf%lgul=.false.
          end if
       end do



    case(2)
       atcf%lgul=.false.
       if (itapp.gt.0) then
          lenfn = 9
          write(extension,'(i9.9)') itapp
          namedes='disorder.'//extension
       else
          namedes='disorder.'
       end if
       naux=1 ; allocate (vaux(naux,atcf%im)) ; allocate(charaux(naux))
       vaux(1,:)=float(nvi(:))
       charaux(1)='neighbours '

       if (any(plmin.ne.0.).or.any(plmax.ne.0.)) then
          call atcf%deftype(atplt)
          call plotpart(atcf,plmin,plmax,atplt,boxplt)
          call rasmolT(atplt,boxplt,namefr=namedes,latcomp=.true.,ivisumol=ivisuana,naux=naux,charaux=charaux,vaux=vaux)
       else
          boxplt=boxcf
          call rasmolT(atcf,boxplt,namefr=namedes,latcomp=.true.,ivisumol=ivisuana,naux=naux,charaux=charaux,vaux=vaux)
       end if

    case default
    end select


    !recalcul puis ecriture des positions
    !  do i = 1, im
    !     latvi=840+nvi(i)
    !     write(latvi, 135) ty(ityp(i)),xp(1,i)*1.0d8,xp(2,i)*1.0d8,xp(3,i)*1.0d8
    !  end do



    ! write(6,*) 'coucou',ldistord,distordflag, icall
  end subroutine nbvois

  !  ------------------------------------------------------------
  !**********************************************************
  !**********************************************************
!!$  subroutine depcr(atc,celc,boxc,atr,celr,boxr)
!!$    USE T_kind_param_m
!!$
!!$    USE tabcr
!!$    implicit none
!!$    !------------pnp-----------------------------------
!!$    !   D u m m y   A r g u m e n t s
!!$    !-----------------------------------------------
!!$
!!$    class(atom_config)::atc,atr
!!$    type(cell_config)::celc,celr
!!$    type(box_config)::boxc,boxr
!!$    
!!$    !Local variables
!!$    integer :: i,j,k,ndep,ic,nplt,idp,iplt
!!$    integer, dimension(:), allocatable:: inddep(:),indplt(:) ! indices des atomes deplaces et plottes
!!$
!!$    real(double) :: tdep2
!!$    real(double) :: a1,a2,a3,c1,c2,c3,r2
!!$    real(double), dimension(1,3) :: cv
!!$
!!$    integer :: nvac,nint,nremp,nanti,ivac,iint,iremp,ias
!!$    integer, dimension(:), allocatable :: indremp
!!$!    logical :: vac
!!$
!!$    real(double) :: r2min
!!$
!!$    real(double) :: c3p,c2p,c1p,c1abs,c2abs,c3abs,r,XJI(3)
!!$    integer:: koo,i2,i1,ncelvois,ko1,immin,immax,iat,jtrf,nas
!!$
!!$    immin=min(atc%im,atr%im)
!!$    immax=max(atc%im,atr%im)
!!$
!!$    write(6,*)'comparison with reference structure', atc%im,atr%im,immin,immax
!!$    ncelvois = min(celc%noxyz,27)-1
!!$
!!$    allocate(indplt(atc%imm))
!!$    allocate(natsit(immax))
!!$    allocate(indatsit(immax,4))
!!$    natsit=0
!!$
!!$    tdep2=tdep*tdep
!!$    ndep=0
!!$    plmin=1000.0 ; plmax=-1000.0
!!$    ndep=0;nvac=0;nanti=0;nint=0;nas=0; 
!!$    nremp=0 ; nplt=0
!!$
!!$
!!$    !write(6,*)'ecr pos'
!!$    !  do i=1,im
!!$    !  write(852,*)i
!!$    !  write(852,*)xp(1,i),xp(2,i),xp(3,i)
!!$    !  write(852,*)xpcr(1,i),xpcr(2,i),xpcr(3,i)
!!$    !    enddo
!!$
!!$    allocate(inddep(immax))
!!$    atc%lgul=.false. ! lgul = true pour les déplacés
!!$    atr%lgul=.false.
!!$    if (ldeptest.EQV..true.) then
!!$       write(6,*)' displacement detection assumes that the atoms are identically sorted in currect and reference state'
!!$       call cryst_to_cart (atc%im, atc%xp, boxc%bg, -1)    !cart vers cryst
!!$       call cryst_to_cart (atr%im, atr%xp, boxr%bg, -1)    !cart vers cryst
!!$
!!$       do i=1,immin
!!$          !     if(ityp(i)==2) cycle
!!$          c1 = atc%xp(1,i)-atr%xp(1,i)
!!$          c2 = atc%xp(2,i)-atr%xp(2,i)
!!$          c3 = atc%xp(3,i)-atr%xp(3,i)
!!$          if (c1>0.5) c1 = c1-1.
!!$          if (c1<(-0.5)) c1 = c1+1.
!!$          if (c2>0.5) c2 = c2-1.
!!$          if (c2<(-0.5)) c2 = c2+1.
!!$          if (c3>0.5) c3 = c3-1.
!!$          if (c3<(-0.5)) c3 = c3+1.
!!$          cv(1,1) = c1
!!$          cv(1,2) = c2
!!$          cv(1,3) = c3
!!$          call cryst_to_cart (1, cv, boxc%at, 1) !cryst vers cart sur cv
!!$          r2 = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)
!!$          if(r2.gt.tdep2) then
!!$             ndep=ndep+1
!!$             inddep(ndep)=i
!!$             atc%lgul(i)=.true.
!!$             atr%lgul(i)=.true.
!!$          else
!!$             indatsit(i,1)=i
!!$             natsit(i)=1
!!$          end if
!!$
!!$       end do  !boucle i
!!$       write(6,*)
!!$       write(6,*)'number of true displaced atoms nombres d atomes deplaces de plus deby more than ',tdep*1.0d8,' = ',ndep
!!$
!!$       if (immin.lt.immax)then
!!$          do i=immin,immax
!!$             ndep=ndep+1
!!$             inddep(ndep)=i
!!$             if (atc%im.lt.immax) then
!!$                atr%lgul(immin+1:immax)=.true.
!!$             else 
!!$                atr%lgul(immin+1:immax)=.true.
!!$             end if
!!$          end do
!!$       end if
!!$       call cryst_to_cart (atc%im, atc%xp, boxc%at, 1)    !cryst vers cart
!!$       call cryst_to_cart (atr%im, atr%xp, boxr%at, 1)    !cryst vers cart
!!$       write(6,*)
!!$       write(6,*)'nombres d atomes deplaces de plus de ',tdep*1.0d8,' = ',ndep
!!$       !       end if
!!$
!!$
!!$    else
!!$       do i=1,immax
!!$          inddep(i)=i
!!$       end do
!!$       ndep=immax
!!$       natsit(:)=0
!!$    end if
!!$
!!$    !determination des atomes plottes
!!$
!!$    if(lvac) then      ! calcul des lacunes et interstitiels
!!$       !ATTENTION CE CALCUL EST LIMITE AU ATOMES DEPLACES si lpdep=.true. CHOISIR TDEP EN CONSEQUENCE
!!$       allocate(indint(immax))
!!$       allocate(indremp(immax))
!!$       allocate(indvac(immax))
!!$       allocate(indas(immax))
!!$       call cryst_to_cart (atc%im, atc%xp, boxc%bg, -1)    !cart vers cryst
!!$       call cryst_to_cart (atr%im, atr%xp, boxr%bg, -1)    !cart vers cryst
!!$
!!$       iloop0:do idp=1,ndep
!!$          r2min =10.0
!!$          i=inddep(idp)
!!$          if (i.gt.atc%im)cycle
!!$          if(ldefcat.and.atc%ityp(i)==2) cycle iloop0
!!$          koo = atc%ielat(i)                          ! Numero de la cellule
!!$!          write(6,*)'i idp ',i,idp
!!$          ! pour chaque cel. voisine
!!$          do i1 = 0, ncelvois
!!$             ko1=celr%ncel(koo,i1)
!!$             !              write(6,*)i,idp,koo,i1,ko1,natocr(ko1)
!!$             do i2 = 1, celr%nato(ko1) !atomes dans la cel dans la conf. init.
!!$                j = celr%atincel(i2,ko1)
!!$
!!$                c1 = atc%xp(1,i)-atr%xp(1,j)
!!$                c2 = atc%xp(2,i)-atr%xp(2,j)
!!$                c3 = atc%xp(3,i)-atr%xp(3,j)
!!$                if (c1>0.5) c1 = c1-1.
!!$                if (c1<(-0.5)) c1 = c1+1.
!!$                if (c2>0.5) c2 = c2-1.
!!$                if (c2<(-0.5)) c2 = c2+1.
!!$                if (c3>0.5) c3 = c3-1.
!!$                if (c3<(-0.5)) c3 = c3+1.
!!$                cv(1,1) = c1
!!$                cv(1,2) = c2
!!$                cv(1,3) = c3
!!$                call cryst_to_cart (1, cv,boxc%at, 1) !cryst vers cart sur cv
!!$                r = sqrt(cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3))
!!$!                write(6,*)'j', j,r
!!$                if(r.lt.tvac) then ! i est sur le site d'un atome du crystal de depart
!!$                   natsit(j)=natsit(j)+1
!!$                   indatsit(j,natsit(j))=i
!!$                   cycle iloop0 
!!$                end if
!!$             end do
!!$          end do
!!$       end do iloop0
!!$
!!$       nvac=0 ;nint=0;nremp=0;nas=0
!!$       iloop1: do j=1,atr%im
!!$          select case (natsit(j))
!!$          case(0) 
!!$             nvac=nvac+1
!!$             indvac(nvac)=j
!!$          case(1)
!!$             if (atr%ityp(j)==atc%ityp(indatsit(j,1)) )then
!!$                if (j.ne.i) then
!!$                   nremp=nremp+1
!!$                   indremp(nremp)=i
!!$                end if
!!$             else
!!$                nas=nas+1
!!$                indas(nas)=i
!!$             end if
!!$          case default
!!$             l3:do iat=1,natsit(j)
!!$                if (atc%ityp(indatsit(j,iat))==atr%ityp(j)) then
!!$                   jtrf=indatsit(j,iat)
!!$                   indatsit(j,iat)=indatsit(j,1)
!!$                   indatsit(j,1)=jtrf
!!$                   exit l3
!!$                end if
!!$             end do l3
!!$             if (atr%ityp(j)==atc%ityp(indatsit(j,1))) then 
!!$                indint(nint+1:nint+natsit(j)-1)=indatsit(j,2:natsit(j))
!!$                nint=nint+natsit(j)-1
!!$             else
!!$                nas=nas+1
!!$                indas(nas)=indatsit(j,1)
!!$             end if
!!$          end select
!!$       end do iloop1
!!$       loopint: do i=1,atc%im
!!$          do j=1,atr%im
!!$             if (natsit(j).ne.0) then
!!$                if (any(indatsit(j,1:natsit(j))==i))then
!!$                   exit loopint
!!$                end if
!!$             end if
!!$          end do
!!$          nint=nint+1
!!$          indint(nint)=i
!!$       end do loopint
!!$
!!$       call cryst_to_cart (atc%im, atc%xp, boxc%at, 1)     !cryst vers cart
!!$       call cryst_to_cart (atr%im, atr%xp, boxr%at, 1)     !cryst vers cart
!!$       write(6,*)'IT = ',it,' nombres de lacunes ',nvac
!!$       write(6,*)'IT = ',it,'nombres d_interstitiels ',nint
!!$       write(6,*)'IT = ',it,'nombres d_antisites ',nanti
!!$    end if
!!$    if(lpdep)     write(6,*)'IT = ',it,'nombres de remplacements ',nremp
!!$    if (lpdef) then
!!$       call  plt_extr(ndep,inddep,'displaced',atc,boxc)
!!$       call  plt_extr(nvac,indvac,'vacancies',atc,boxc)
!!$       call  plt_extr(nint,indint,'interstitials',atc,boxc)
!!$       call  plt_extr(nas,indas,'antisites',atc,boxc)
!!$       call  plt_extr(nremp,indremp,'replacements',atc,boxc)
!!$    end if
!!$    
!!$    if(lvac) then
!!$       deallocate (indvac) ; deallocate (indint) ; deallocate (indas) ;deallocate(indremp)
!!$       deallocate (indatsit) ; deallocate(indatsit)
!!$    end if
!!$    
!!$  end subroutine depcr




  subroutine ws(atc,celc,boxc,atr,celr,boxr,lws)
    USE T_kind_param_m

    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    class(atom_config)::atc,atr
    type(cell_config)::celc,celr
    type(box_config)::boxc,boxr
    logical,intent(in)::lws


    !Local variables
    integer :: i,j,nplt,idp,iplt,immax,immin,koo


    real(double) :: a1,a2,a3,c1,c2,c3,r2
    real(double), dimension(1,3) :: cv
    integer, dimension(:), allocatable:: inddep(:)
    integer :: nvac,nint,nremp,nanti,ivac,iint,iremp,ias,nas
    logical :: vacfl
    integer::iat,jtrf
    real(double) :: r2min,r2ali

    real(double) :: c3p,c2p,c1p,c1abs,c2abs,c3abs,r,tdep2
    integer:: i1,ko1,indws0,i2,ndep
    !  integer,allocatable :: lastcr (:,:),natocr(:),ielatcr(:)
    integer, dimension(:), allocatable:: indplt(:),indremp(:) ! indices des atomes deplaces et plottes
    immin=min(atc%im,atr%im)
    immax=max(atc%im,atr%im)
    if (lws) then
       if (rang==0)write(6,*)'comparison with reference structure WIGNER SEITZ', atc%im,atr%im,immin,immax
    else
       write(6,*)'comparison with reference structure cut-off radius', atc%im,atr%im,immin,immax,tvac
       if (ldeptest) then
          if (rang==0)write(6,*)'comparison based on displacements', tdep
       end if
    end if

    allocate(indplt(atc%imm))
    allocate(natsit(immax))
    allocate(indatsit(immax,4))
    natsit=0
    indatsit=0
    if (lvac) then 
       allocate(indint(atc%imm))
       allocate(indvac(atc%imm))
       allocate(indas(atc%imm))
       allocate(indremp(atc%imm))
       allocate(indws(atc%imm))
    end if

    nremp=0;nplt=0

    nvac=0;nanti=0;nint=0;nas=0; nremp=0

    if (lws) then
       call wsb(atc,celc,boxc,atr,celr,boxr,natsit,indatsit,indws)


    else ! test depcr

       tdep2=tdep*tdep
       allocate(inddep(immax))
       atc%lgul=.false. ! lgul = true pour les déplacés
       atr%lgul=.false.
       if (ldeptest.EQV..true.) then
          call cryst_to_cart (atc%im, atc%xp, boxc%bg, -1)    !cart vers cryst
          call cryst_to_cart (atr%im, atr%xp, boxr%bg, -1)    !cart vers cryst

          write(6,*)' displacement detection assumes that the atoms are identically sorted in currect and reference state'

          do i=1,immin
             !     if(ityp(i)==2) cycle
             c1 = atc%xp(1,i)-atr%xp(1,i)
             c2 = atc%xp(2,i)-atr%xp(2,i)
             c3 = atc%xp(3,i)-atr%xp(3,i)
             if (boxc%ipbc(1)==1) then
                if (c1>0.5) c1 = c1-1.
                if (c1<(-0.5)) c1 = c1+1.
             end if
             if (boxc%ipbc(2)==1) then
                if (c2>0.5) c2 = c2-1.
                if (c2<(-0.5)) c2 = c2+1.
             end if
             if (boxc%ipbc(3)==1) then
                if (c3>0.5) c3 = c3-1.
                if (c3<(-0.5)) c3 = c3+1.
             end if
             cv(1,1) = c1
             cv(1,2) = c2
             cv(1,3) = c3
             call cryst_to_cart (1, cv, boxc%at, 1) !cryst vers cart sur cv
             r2 = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)
             if(r2.gt.tdep2) then
                ndep=ndep+1
                inddep(ndep)=i
                atc%lgul(i)=.true.
                atr%lgul(i)=.true.
             else
                indatsit(i,1)=i
                natsit(i)=1
             end if

          end do  !boucle i
          write(6,*)
          write(6,*)'number of true displaced atoms nombres d atomes deplaces de plus deby more than ',tdep*1.0d8,' = ',ndep

          if (immin.lt.immax)then
             do i=immin,immax
                ndep=ndep+1
                inddep(ndep)=i
                if (atc%im.lt.immax) then
                   atr%lgul(immin+1:immax)=.true.
                else 
                   atr%lgul(immin+1:immax)=.true.
                end if
             end do
          end if
          write(6,*)
          write(6,*)'nombres d atomes deplaces de plus de ',tdep*1.0d8,' = ',ndep
          !       end if

          call cryst_to_cart (atc%im, atc%xp, boxc%at, 1)    !cryst vers cart
          call cryst_to_cart (atr%im, atr%xp, boxr%at, 1)    !cryst vers cart
       else
          do i=1,immax
             inddep(i)=i
          end do
          ndep=immax
          natsit(:)=0
       end if

       !determination des atomes plottes

       if(lvac) then      ! calcul des lacunes et interstitiels
          !ATTENTION CE CALCUL EST LIMITE AU ATOMES DEPLACES si lpdep=.true. CHOISIR TDEP EN CONSEQUENCE
          call cryst_to_cart (atc%im, atc%xp, boxc%bg, -1)    !cart vers cryst
          call cryst_to_cart (atr%im, atr%xp, boxr%bg, -1)    !cart vers cryst

          iloop00:do idp=1,ndep
             r2min =10.0
             i=inddep(idp)
             if (i.gt.atc%im)cycle
             if(ldefcat.and.atc%ityp(i)==2) cycle iloop00
             koo = atc%ielat(i)                          ! Numero de la cellule
             !          write(6,*)'i idp ',i,idp
             ! pour chaque cel. voisine
             do i1 = 0, celr%ncelvois(koo)
                ko1=celr%ncel(koo,i1)
                !              write(6,*)i,idp,koo,i1,ko1,natocr(ko1)
                do i2 = 1, celr%nato(ko1) !atomes dans la cel dans la conf. init.
                   j = celr%atincel(i2,ko1)

                   c1 = atc%xp(1,i)-atr%xp(1,j)
                   c2 = atc%xp(2,i)-atr%xp(2,j)
                   c3 = atc%xp(3,i)-atr%xp(3,j)
                   if (boxc%ipbc(1)==1) then
                      if (c1>0.5) c1 = c1-1.
                      if (c1<(-0.5)) c1 = c1+1.
                   end if
                   if (boxc%ipbc(2)==1) then
                      if (c2>0.5) c2 = c2-1.
                      if (c2<(-0.5)) c2 = c2+1.
                   end if
                   if (boxc%ipbc(3)==1) then                      
                      if (c3>0.5) c3 = c3-1.
                      if (c3<(-0.5)) c3 = c3+1.
                   end if
                   cv(1,1) = c1
                   cv(1,2) = c2
                   cv(1,3) = c3
                   call cryst_to_cart (1, cv,boxc%at, 1) !cryst vers cart sur cv
                   r = sqrt(cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3))
                   !                write(6,*)'j', j,r
                   if(r.lt.tvac) then ! i est sur le site d'un atome du crystal de depart
                      natsit(j)=natsit(j)+1
                      indatsit(j,natsit(j))=i
                      cycle iloop00 
                   end if
                end do
             end do
          end do iloop00
          call cryst_to_cart (atc%im, atc%xp, boxc%at, 1)    !cryst vers cart
          call cryst_to_cart (atr%im, atr%xp, boxr%at, 1)    !cryst vers cart

       end if
    end if
    if (lvac) then 

       nvac=0 ;nint=0;nremp=0;nas=0
       iloop1: do j=1,atr%im
          select case (natsit(j))
          case(0) 
             nvac=nvac+1
             indvac(nvac)=j
          case(1)

             if (atr%ityp(j)==atc%ityp(indatsit(j,1)) )then
                i=indatsit(j,1)
                if (j.ne.i) then
                   nremp=nremp+1
                   indremp(nremp)=i
                end if
             else
                nas=nas+1
                indas(nas)=i
             end if
          case default
             l31:do iat=1,natsit(j)
                if (atc%ityp(indatsit(j,iat))==atr%ityp(j)) then
                   jtrf=indatsit(j,iat)
                   indatsit(j,iat)=indatsit(j,1)
                   indatsit(j,1)=jtrf
                   exit l31
                end if
             end do l31
             if (atr%ityp(j)==atc%ityp(indatsit(j,1))) then 
                indint(nint+1:nint+natsit(j)-1)=indatsit(j,2:natsit(j))
                nint=nint+natsit(j)-1
             else
                nas=nas+1
                indas(nas)=indatsit(j,1)
             end if
          end select
       end do iloop1
       loopint: do i=1,atc%im
          do j=1,atr%im
             if (natsit(j).ne.0) then
                if (any(indatsit(j,1:natsit(j))==i) )then
                   exit loopint
                end if
             end if
          end do
          nint=nint+1
          indint(nint)=i
       end do loopint

       write(6,*)'IT = ',iteration,' nombres de lacunes ',nvac
       write(6,*)'IT = ',iteration,'nombres d_interstitiels ',nint
       write(6,*)'IT = ',iteration,'nombres d_antisites ',nanti
       write(6,*)'IT = ',iteration,'nombres d_remplacments ',nremp

       if(lpdep)     write(6,*)'IT = ',iteration,'nombres de remplacements ',nremp
       if (lpdef) then
          call  plt_extr(nvac,indvac,'vacancies',atc,boxc)
          call  plt_extr(nint,indint,'interstitials',atc,boxc)
          call  plt_extr(nas,indas,'antisites',atc,boxc)
          call  plt_extr(nremp,indremp,'replacements',atc,boxc)
       end if
    end if

    if(lvac) then
       deallocate (indvac) ; deallocate (indint) ; deallocate (indas) ;deallocate(indremp)
       deallocate (indatsit) 
    end if
  end subroutine ws


!!$  subroutine subc(nvac,indvac,nint,indint)
!!$    USE T_kind_param_m
!!$
!!$    type:: deftype
!!$       real(double)::xd(3)
!!$       integer::typ ! 1=vac ; 2=int ; 3=AS
!!$       integer::attyp 
!!$       integer:: sc ! sous cascade
!!$       integer::nvd !nombre de d�fauts voisins
!!$       integer::nvdvois !nombre max de d�fauts parmi les atomes voisins
!!$       integer,dimension(:),allocatable::indvd !indice des d�fauts voisins
!!$    end type deftype
!!$
!!$    integer ::nvac,nint,ntypdefsc(3),nattypdefsc(ntyp),ivac,iint,nsc,id12,sc1,sc2,id3sc,id3,&
!!$         & isc,jrsc,jsc,krsc,ksc,irsc
!!$    integer, dimension(:), allocatable:: indvac,indint
!!$    integer,allocatable::ndefsc (:),inddefsc(:,:),indrg(:),rgsc(:),ndefvois(:)
!!$    type(deftype), allocatable :: deft(:), deft2(:)
!!$    real(double)::c1,c2,c3,cv(1,3),dist
!!$    integer::ndeft,ndeft2,id1,id2,id
!!$    character*2::ch2
!!$    integer, dimension (1000):: nscIn,nscVn
!!$!    write(6,*)'indvac',indvac(1:nvac)
!!$!    write(6,*)'indint',indint(1:nint)
!!$
!!$    ndeft=nvac+nint
!!$    allocate (deft(ndeft)) 
!!$    allocate (ndefvois(ndeft))
!!$    do id=1,ndefT
!!$       allocate(defT(id)%indvd(ndeft))
!!$    end do
!!$    id=0
!!$    do ivac=1,nvac
!!$       id=id+1  
!!$
!!$       do ic=1,3
!!$          defT(id)%xd(ic)=xp(ic,indvac(ivac))
!!$!          write(6,*)xp(ic,indvac(ivac))
!!$       end do
!!$       deft(id)%typ=1
!!$       deft(id)%nvd=0
!!$
!!$    end do
!!$
!!$    do iint=1,nint
!!$       id=id+1
!!$       do ic=1,3
!!$          defT(id)%xd(ic)=xp(ic,indint(iint))
!!$       end do
!!$       deft(id)%typ=2
!!$       deft(id)%nvd=0
!!$       deft(id)%attyp=ityp(indint(iint))
!!$    end do
!!$    write(6,*)'nb de defauts pour SC=',ndeft,nvac
!!$!    do id=1,ndeft
!!$!       write(6,*)id,deft(id)%xd
!!$!    end do
!!$    ndefvois=0
!!$    write(6,*)'TATA'
!!$    do id1=1,ndeft
!!$       do id2=id1+1,ndeft
!!$          cv(1,1) = deft(id1)%xd(1)-deft(id2)%xd(1)
!!$          cv(1,2) = deft(id1)%xd(2)-deft(id2)%xd(2)
!!$          cv(1,3) = deft(id1)%xd(3)-deft(id2)%xd(3)
!!$          call cryst_to_cart (1, cv, bg, -1) !cryst vers cart sur cv
!!$          if (cv(1,1)>0.5) cv(1,1) = cv(1,1)-1.
!!$          if (cv(1,1)<(-0.5)) cv(1,1) = cv(1,1)+1.
!!$          if (cv(1,2)>0.5) cv(1,2) = cv(1,2)-1.
!!$          if (cv(1,2)<(-0.5)) cv(1,2) = cv(1,2)+1.
!!$          if (cv(1,3)>0.5) cv(1,3) = cv(1,3)-1.
!!$          if (cv(1,3)<(-0.5)) cv(1,3) = cv(1,3)+1.
!!$          call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
!!$          dist= sqrt(cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3))
!!$!          write(6,*)'dist',dist,rdv
!!$          if (dist.le.(rdv*1d-8)) then
!!$             deft(id1)%nvd=deft(id1)%nvd+1
!!$             deft(id1)%indvd(deft(id1)%nvd)=id2
!!$             deft(id2)%nvd=deft(id2)%nvd+1
!!$             deft(id2)%indvd(deft(id2)%nvd)=id1
!!$          end if
!!$       end do
!!$    end do
!!$    write(6,*)'TOTO'
!!$
!!$    do id1=1,ndeft
!!$       deft(id1)%nvdvois=deft(id1)%nvd
!!$       do id2=1,deft(id1)%nvd
!!$          deft(id1)%nvdvois=max(deft(id1)%nvdvois,deft(deft(id1)%indvd(id2))%nvd)
!!$       end do
!!$    end do
!!$
!!$    ndeft2=0
!!$    do id1=1,ndeft
!!$!       write(6,*)deft(id1)%nvd,deft(id1)%nvdvois
!!$       if(deft(id1)%nvdvois.ge.ndvblob)then
!!$          ndeft2=ndeft2+1
!!$       end if
!!$    end do
!!$    write(6,*)'nb de defauts dans les SC ',ndeft2
!!$    allocate (deft2(ndeft2))
!!$    do id2=1,ndefT2
!!$       allocate(defT2(id2)%indvd(ndeft2))
!!$    end do
!!$    id2=0
!!$
!!$    do id1=1,ndeft
!!$       if(deft(id1)%nvdvois.ge.ndvblob)then
!!$          id2=id2+1
!!$          deft2(id2)=deft(id1) 
!!$          deft(id1)%sc=id2
!!$          deft2(id2)%indvd(:)=0
!!$          deft2(id2)%nvd=0
!!$       end if
!!$    end do
!!$
!!$    do id1=1,ndeft2
!!$       do id2=id1+1,ndeft2
!!$          cv(1,1) = deft2(id1)%xd(1)-deft2(id2)%xd(1)
!!$          cv(1,2) = deft2(id1)%xd(2)-deft2(id2)%xd(2)
!!$          cv(1,3) = deft2(id1)%xd(3)-deft2(id2)%xd(3)
!!$          call cryst_to_cart (1, cv, bg, -1) !cryst vers cart sur cv
!!$          if (cv(1,1)>0.5) cv(1,1) = cv(1,1)-1.
!!$          if (cv(1,1)<(-0.5)) cv(1,1) = cv(1,1)+1.
!!$          if (cv(1,2)>0.5) cv(1,2) = cv(1,2)-1.
!!$          if (cv(1,2)<(-0.5)) cv(1,2) = cv(1,2)+1.
!!$          if (cv(1,3)>0.5) cv(1,3) = cv(1,3)-1.
!!$          if (cv(1,3)<(-0.5)) cv(1,3) = cv(1,3)+1.
!!$          call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
!!$          dist= cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)
!!$          if (dist.le.(rdv*1d-8)**2) then
!!$             deft2(id1)%nvd=deft2(id1)%nvd+1
!!$             deft2(id1)%indvd(deft2(id1)%nvd)=id2
!!$             deft2(id2)%nvd=deft2(id2)%nvd+1
!!$             deft2(id2)%indvd(deft2(id2)%nvd)=id1
!!$          end if
!!$       end do
!!$    end do
!!$    do id1=1,ndeft2
!!$       deft2(id1)%sc=id1
!!$    end do
!!$
!!$    nsc=ndeft2
!!$    allocate(ndefsc(nsc))
!!$    ndefsc(:)=1
!!$    allocate(inddefsc(ndeft2,nsc))
!!$    inddefsc(:,:)=0
!!$    do id=1,nsc
!!$       inddefsc(1,id)=id
!!$    end do
!!$
!!$    do id1=1,ndeft2
!!$       do id12=1,deft2(id1)%nvd
!!$          id2=deft2(id1)%indvd(id12)
!!$          if(deft2(id1)%sc.ne.deft2(id2)%sc)then
!!$             if (deft2(id1)%sc.lt.deft2(id2)%sc) then
!!$                sc1=deft2(id1)%sc
!!$                sc2=deft2(id2)%sc
!!$             else
!!$                sc1=deft2(id2)%sc
!!$                sc2=deft2(id1)%sc
!!$             end if
!!$             !                          write(6,*)'ndefc'
!!$             !                          write(6,*)ndefsc(sc1)
!!$             !                          write(6,*)ndefsc(sc2)
!!$             do id3sc=1,ndefsc(sc2)
!!$                id3=inddefsc(id3sc,sc2)
!!$                ndefsc(sc1)=ndefsc(sc1)+1
!!$                inddefsc(ndefsc(sc1),sc1)=id3
!!$                deft2(id3)%sc=sc1
!!$             end do
!!$             do isc=sc2+1,nsc
!!$                ndefsc(isc-1)=ndefsc(isc)
!!$                do id3sc=1,ndefsc(isc-1)
!!$                   id3=inddefsc(id3sc,isc)
!!$                   inddefsc(id3sc,isc-1)=inddefsc(id3sc,isc)
!!$                   deft2(id3)%sc=isc-1
!!$                end do
!!$             end do
!!$             nsc=nsc-1
!!$          end if
!!$       end do
!!$    end do
!!$
!!$
!!$    write(6,'(A,G14.5,I4,A,I6,A)')'POUR RDV/ndvblob =',RDV,ndvblob,' il y a ', nsc,' sous cascades'
!!$
!!$    allocate(indrg(nsc))
!!$    allocate(rgsc(nsc))
!!$    do isc=1,nsc
!!$       rgsc(isc)=isc
!!$       indrg(isc)=isc
!!$    end do
!!$    do isc=1,nsc  !on ordonne les cascades
!!$       !       write(6,*)'ISC',isc
!!$       do jrsc=1,isc-1  ! cascades ordonn�es
!!$          jsc=indrg(jrsc)  ! indices de la jrsc �me cascade
!!$          if (ndefsc(jsc).ge.ndefsc(isc))cycle
!!$          do krsc=isc-1,jrsc,-1  ! krsc cascades suivantes 
!!$             ksc=indrg(krsc) 
!!$             rgsc(ksc)=krsc+1
!!$             indrg(krsc+1)=ksc
!!$          end do
!!$          rgsc(isc)=jrsc
!!$          indrg(jrsc)=isc
!!$          exit
!!$       end do
!!$       !       do jrsc=1,isc
!!$       !          write(6,*)'ndef scR',jrsc,indrg(jrsc),ndefsc(indrg(jrsc))
!!$       !       end do
!!$    end do
!!$    nclustI=0; nclustV=0
!!$    do irsc=1,nsc
!!$       isc=indrg(irsc)
!!$       ntypdefsc(1:3)=0
!!$       nattypdefsc(:)=0
!!$       do id=1,ndefsc(isc)
!!$          id2=inddefsc(id,isc)
!!$          ntypdefsc(deft2(id2)%typ)=ntypdefsc(deft2(id2)%typ)+1
!!$          nattypdefsc(deft2(id2)%attyp)=nattypdefsc(deft2(id2)%attyp)+1
!!$          !             write(6,*)
!!$       end do
!!$       write(6,'(A,I7,I7,A,8I7)')'sous cascade ',isc,ndefsc(isc),' defauts',ntypdefsc(1:3),nattypdefsc(1:ntyp)
!!$       if (ntypdefsc(2)==ndefsc(isc))          nscIn(ndefsc(isc))=nscIn(ndefsc(isc))+1
!!$       if (ntypdefsc(1)==ndefsc(isc))          nscVn(ndefsc(isc))=nscVn(ndefsc(isc))+1
!!$    end do
!!$    if (ndeft.ne.ndeft2) write(6,*)'defauts isoles ', ndeft-ndeft2
!!$    write(6,*)
!!$       do i=1,maxval(ndefsc)
!!$       write(6,*)'nclustI ', i,' = ',nscIn(i)
!!$       enddo
!!$       write(6,*)
!!$       do i=1,maxval(ndefsc)
!!$       write(6,*)'nclustV ', i,' = ',nscVn(i)
!!$       enddo
!!$       write(6,*)
!!$
!!$
!!$
!!$    open(file='def_et_SC.mol', unit=182)
!!$    write(182,*)ndeft,'IT = ',it,' DEF ET SC'
!!$    at=at*1.d8
!!$    write (182,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),&
!!$         at(1,3),at(2,3),at(3,3)
!!$    at=at/1.d8
!!$119 format(a2,1x,3(G17.8,1x),1x,i5)       
!!$    do id=1,ndeft
!!$       if (deft(id)%typ==1)then
!!$          ch2='V '
!!$       else
!!$          ch2=ty(deft(id)%attyp)
!!$       end if
!!$
!!$       if(deft(id)%sc==0) then
!!$          write(182, 119) ch2 ,deft(id)%xd(1)*1d8,deft(id)%xd(2)*1d8,deft(id)%xd(3)*1d8,deft(id)%sc
!!$       else
!!$          write(182, 119) ch2 ,deft(id)%xd(1)*1d8,deft(id)%xd(2)*1d8,deft(id)%xd(3)*1d8,deft2(deft(id)%sc)%sc
!!$       end if
!!$    end do
!!$
!!$
!!$    deallocate (deft) ; 
!!$    deallocate (deft2);
!!$    deallocate (ndefvois);
!!$    deallocate(ndefsc);
!!$    deallocate(inddefsc);
!!$    deallocate(indrg);
!!$    deallocate(rgsc);
!!$  end subroutine subc

  subroutine plt_extr(nprt,indprt,nameprt,atprt,boxprt,itapp,plmin,plmax)
    integer,intent(in)::nprt,indprt(:)
    integer,optional::itapp
    character(len=*)::nameprt
    class(atom_config)::atprt
    type(box_config)::boxprt
    real(double),optional :: plmin(3,3),plmax(3,3)
    integer::iprt,iat
    class(atom_config),allocatable::atplt,atpltpp
    type(box_config)::boxplt
    character*9 :: extension
    character(len=89)::nameplt

    if (nprt.gt.0) then
       !         write(6,*)'nameprt',nameprt
       call atprt%deftype(atplt)
       atprt%lgul=.false.
       do iat=1,nprt
          atprt%lgul(indprt(iat))=.true.
       end do
       call atprt%fab(atplt,lback=.false.)
       atprt%lgul=.false.
       if (present(itapp)) then
          write(extension,'(i9.9)') itapp
          nameplt=trim(nameprt)//extension
       else
          nameplt=trim(nameprt)
       end if
       !          write(6,*)'nameplt ',nameprt,' ',nameplt
       if (present(plmin)) then
          call atprt%deftype(atpltpp)
          call plotpart(atplt,plmin,plmax,atpltpp,boxplt)
          call rasmolT(atpltpp,boxplt,namefr=nameplt,latcomp=.true.,ivisumol=ivisuana)
       else
          boxplt=boxprt
          call rasmolT(atplt,boxplt,namefr=nameplt,latcomp=.true.,ivisumol=ivisuana)
       end if
    end if
  end subroutine plt_extr


  !**************** PLOT PART****

  subroutine plotpart(atcf,plmin,plmax,atplt,boxplt)
    USE T_kind_param_m

    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    class(atom_config)::atcf
    !    integer :: ityp(imm)
    class(atom_config)::atplt
    class(box_config)::boxplt

    real(double) ::  plmin(3),plmax(3),at_plt(3,3)
    !local variables 
    integer :: i,nplt,iplt,ic
    real(double),dimension(:,:),allocatable::vectdepla
    !    integer, allocatable :: indplt(:) ! indices des atomes  plottes
    !    allocate(indplt(imm))
    !    plmin=plmin*1.0d-8 ; plmax=plmax*1.0d-8
    write(6,'(A,3F7.2,A,3F7.2)')'portion affichee entre ', plmin*1.0d8 ,' et ',plmax*1.0d8  
    !    write(6,*)'portion affichee entre ', plmin*1.0d8 ,' et ',plmax*1.0d8  
    nplt=0
    atcf%lgul=.false.
    do i=1,atcf%im     
       if((atcf%xp(1,i).gt.plmin(1).and.atcf%xp(1,i).lt.plmax(1)).and.&
            & (atcf%xp(2,i).gt.plmin(2).and.atcf%xp(2,i).lt.plmax(2)).and. &
            & (atcf%xp(3,i).gt.plmin(3).and.atcf%xp(3,i).lt.plmax(3))) then
          atcf%lgul(i)=.true.
       end if
    end do
    call atcf%fab(atplt,lback=.false.)
    atcf%lgul=.false.
    do i=1,atplt%im
       atplt%xp(1:3,i)=atplt%xp(1:3,i)-plmin(1:3)
    end do
    at_plt(:,:)=0
    do ic=1,3
       at_plt(ic,ic)=plmax(ic)-plmin(ic)
    end do
    call boxplt%init(at_plt,ipbc)
    !    namepltpart='partial'
    !    itapp=it
    !    call rasmolT(atplt,boxplt,namefr=namepltpart,latcomp=.true.,ivisumol=ivisuana)
    !    call rasmolT(atplt,boxplt,itap,namepltpart,latcomp=.true.,ivisumol=ivisuana)
  end subroutine plotpart


  subroutine anaposart(atcf,celcf,boxcf,lperfw,namemol,lnewref)

    class(atom_config)::atcf
    class(cell_config)::celcf
    class(box_config)::boxcf

    logical::lperfw
    character(len=60) :: namemol
    character(len=60) :: namemolperf

    logical,optional::lnewref

    type(atom_config)::atperf
    type(cell_config)::celperf
    integer, allocatable:: indws(:),indatsit(:,:),natsit(:)
    real(double),dimension(:,:),allocatable::vectdep
    real(double),dimension(:),allocatable::valdepla
    integer,dimension(:),allocatable::inddep

    integer::i,j,ioc,ndep,id
    logical::lnrf
    logical,save::lpsp=.true.
    lnrf=.false.
    if(present(lnewref))lnrf=lnewref

    if (lnrf) then
       atrefdep=atcf
       celrefdep=celcf
       boxrefdep=boxcf
    end if


    allocate(indws(atcf%im))
    allocate(natsit(atcr%im))
    allocate(indatsit(atcr%im,4))
    natsit=0
    indatsit=0


    if (.Not.lnrf) then
       if (lcomp) then
          call wsb(atcf,celcf,boxcf,atcr,celcr,boxcr,natsit,indatsit,indws)
          write(unit6P,*)
          do j=1,atcr%im
             if ((atcr%ityp(j).ne.0).and.(natsit(j)==0)) then ! WS vacancy
                write(unit6P,'(A,I7,3G15.7,I2)')'vacancy in site ',j,atcr%xp(1,j),atcr%xp(2,j),atcr%xp(3,j),atcr%ityp(j)
             end if
          end do
          do j=1,atcr%im
             if (natsit(j).gt.1) then ! WS multiple occupation
                do ioc=1,natsit(j)
                   i=indatsit(j,ioc)
                   write(unit6P,'(A,I7,3G15.7,A,I2,A,I7,A,I2)')'multiple occupancy in site ',j,  1d8*atcr%xp(1,j)&
                        &,1d8*atcr%xp(2,j),1d8*atcr%xp(3,j),&
                        &'of ref type ',atcr%ityp(j),'with atom ',i, ' of type ',atcf%ityp(i)
                end do
             end if
          end do
          do i=1,atcf%im
             j=indws(i)
             if (atcr%ityp(j)==0) then
                write(unit6P,'(A,I7,3G15.7,A,I7,A,I2,3G15.7)')'interstitial in site ',j,&
                     &1d8*atcr%xp(1,j),1d8*atcr%xp(2,j),1d8*atcr%xp(3,j)&
                     &, 'of index',i,' and type ', atcf%ityp(i), 1d8*atcr%xp(1,j)/27.273,&
                     &1d8*atcr%xp(2,j)/27.273,1d8*atcr%xp(3,j)/27.273
             end if
          end do
       end if


       write(unit6P,*)'DEP0',tdep
       call depladet(atcf,celcf,boxcf,atrefdep,celrefdep,boxrefdep,tdep,ndep,inddep,vectdep,valdepla)
       write(unit6P,*)ndep, ' displaced atoms'
       do id=1,ndep
          write(unit6P,'(A,I7,A,A,4G15.7)')'DISPLACED ATOM',inddep(id),'  ',ty(atcf%ityp(inddep(id))),valdepla(id)*1d8,&
               &vectdep(1,id)*1d8,vectdep(2,id)*1d8,vectdep(3,id)*1d8
       end do
    end if

!    write(unit6P,*)
    if (lnrf) then
       atrefdep=atcf
       celrefdep=celcf
       boxrefdep=boxcf

!!$       if (lcomp.all(natsit.le.1)) then
!!$          write(unit6P,*)'new perfect structure is possible'
!!$          lpsp=.true.
!!$          !          atperfdef=atcf
!!$          !          celperfdef=celcf
!!$          !          boxperfdef=boxcf
!!$          !          atperfdef%xp(:,1:atcf%im)=atcr%xp(:,1:atcf%im)
!!$       else
!!$          lpsp=.false.
!!$       end if
       return
    end if


    if (lperfdef) then 
       if (all(natsit.le.1)) then
          write(unit6P,*)'output of perfect structure is possible'
          call atperf%init(atcf%im)
          do i=1,atcf%im
             atperf%ityp(i)=atcf%ityp(i)
             atperf%xp(:,i)=atcr%xp(:,indws(i))
          end do
          namemolperf=trim(namemol)//'perf'
          call rasmolT(atperf,boxcr,namefr=namemolperf,latcomp=.true.,ivisumol=5)


          if (lpsp) then          
             call celcf%copy(celperf,boxcf)

             call caltabtc(celperf,atperf,lperiod,boxcf,lchktrav=.false.)
             call depladet(atperf,celperf,boxcf,atperfdef,celperfdef,boxperfdef,tdep,ndep,inddep,vectdep,valdepla)
             write(unit6P,*)ndep, ' displaced atoms perfect structure'
             valdepla=valdepla*(1d8)*5./27.273
             vectdep=vectdep*(1d8)*5./27.273
             do id=1,ndep
                i=inddep(id)
                write(unit6P,'(A,I7,A,A,4G15.7,A,3G15.7,A,3G15.7)')'DISPLACED ATOM',inddep(id),'  '&
                     &,ty(atcf%ityp(inddep(id))), valdepla(id),&
                     &vectdep(1,id),vectdep(2,id),vectdep(3,id),' INIT ',&
                     &atperfdef%xp(1,i)*(1d8)/27.273,atperfdef%xp(2,i)*(1d8)/27.273,atperfdef%xp(3,i)*(1d8)/27.273,&
                     &' FINAL ',atperf%xp(1,i)*(1d8)/27.273,atperf%xp(2,i)*(1d8)/27.273,atperf%xp(3,i)*(1d8)/27.273
             end do
             write(unit6P,*)
          end if
       end if
    end if


    return

  end subroutine anaposart

  subroutine wsb(atc,celc,boxc,atr,celr,boxr,natsit,indatsit,indws)

    class(atom_config)::atc,atr
    type(cell_config)::celc,celr
    type(box_config)::boxc,boxr

    integer,dimension(:),intent(out)::natsit,indws
    integer,dimension(:,:),intent(out)::indatsit

    real(double), dimension(1,3) :: cv
    integer::i,koo,j,i2,i1,ko1
    real(double)::r2min,c1,c2,c3,r
    natsit=0;indws=0;indatsit=0


    call cryst_to_cart (atc%im, atc%xp, boxc%bg, -1)    !cart vers cryst
    call cryst_to_cart (atr%im, atr%xp, boxr%bg, -1)    !cart vers cryst

    iloop0:do i=1,atc%im
       r2min =100.0
       koo = atc%ielat(i)                          ! Numero de la cellule
       ! pour chaque cel. voisine
       do i1 = 0, celr%ncelvois(koo)
          ko1=celr%ncel(koo,i1)

          !              write(6,*)i,idp,koo,i1,ko1,natocr(ko1)
          do i2 = 1, celr%nato(ko1) !atomes dans la cel dans la conf. init.
             j = celr%atincel(i2,ko1)

             c1 = atc%xp(1,i)-atr%xp(1,j)
             c2 = atc%xp(2,i)-atr%xp(2,j)
             c3 = atc%xp(3,i)-atr%xp(3,j)
             if (boxc%ipbc(1)==1) then
                if (c1>0.5) c1 = c1-1.
                if (c1<(-0.5)) c1 = c1+1.
             end if
             if (boxc%ipbc(2)==1) then
                if (c2>0.5) c2 = c2-1.
                if (c2<(-0.5)) c2 = c2+1.
             end if
             if (boxc%ipbc(3)==1) then
                if (c3>0.5) c3 = c3-1.
                if (c3<(-0.5)) c3 = c3+1.
             end if
             cv(1,1) = c1
             cv(1,2) = c2
             cv(1,3) = c3
             call cryst_to_cart (1, cv, boxc%at, 1) !cryst vers cart sur cv
             r = sqrt(cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3))

             if(r.lt.r2min) then ! j est pour l'instant le site le plus proche de i
                !                write(6,*)r,i,j
                r2min=r
                indws(i)=j
             end if
          end do
       end do


       j=indws(i)
       natsit(j)=natsit(j)+1
       indatsit(j,natsit(j))=i
    end do iloop0
    call cryst_to_cart (atc%im, atc%xp, boxc%at, 1)    !cryst vers cart
    call cryst_to_cart (atr%im, atr%xp, boxr%at, 1)    !cryst vers cart


  end subroutine wsb

end module posana
