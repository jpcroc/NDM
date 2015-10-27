module posana
    USE T_kind_param_m
  use sic
  use var_pot
  logical :: lsic
    logical :: lcomp, & ! comparaison ou non avec un cristal de dÃ©part
         ldecal, & ! decalage en tre boite cr et boite ana
         ldesord, &      ! vielle variable historique
         lnbvois, &      ! analyse des nombres de voisins
         lc15, &      ! lacunes et int ensembles
         lvac, &           ! analyse en lacune
         ldetdec,&          !determination du decalage
         lrescale ,&          ! rescale des posistions de dÃ©part sur la boite d'arrivÃ©e
         lpstruct ,&          
         lpdep ,&          
         ldeptest ,&          
         lws ,&          ! Wigner-Seitz pour INT et VAC
         lpdef, &
         ldefcat, &        ! defauts sur les cations seulement
         distordflag    ! analyse des angles dans le cristal si flag==.true.
    integer::         idistord,iprtnvi        ! analyse des diff angulaires
    integer :: imcr ! nb d'atomes dans le cristal de reference
    real(double)::pstmax(20)
    real(double),pointer:: xpcr(:,:)
    integer, pointer :: itypcr(:)

  ! **************************************************************
contains

  subroutine anapos(itapp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m
    use gen_com_m
    use tab_imm_m
#if(PARA)
    use mod_mpi
#endif

    ! **************************************************************


    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: i, ic,j,k,l,m,n,itapp,nbvoisparf(20,20)
    real(double)  :: decal(3)
    real(double) :: tdep ! seuil de deplacement
    real(double) :: plmin(3),plmax(3) ! bords de la portion afichÃ©e de la boite
    integer :: idecal    ! alignement des posistions sur l'atome idecal
    integer, save:: icall=0

    real(double) :: tvac ,tint,deltx,delty,deltz ! distance pour les lacunes et les int
    real (double) :: plmin1,plmin2,plmin3, plmax1,plmax2,plmax3 ! bord de plot lu dans la namelist
    namelist /analyse/ldecal,ldesord,idistord,idecal,tdep,lvac,tvac,tint,lcomp,plmin1, &
         plmin2,plmin3, plmax1,plmax2,plmax3,ldetdec,lrescale,lpstruct,lpdef,lpdep,ldeptest, &
         ldefcat,rclu, lnbvois,lsic,nbvoisparf,pstmax,iprtnvi,lc15,lws,deltx,delty,deltz




    !
    !   set default values for variables in namelist
    !

    !-----------------------------------------------
    rclu(:)=2.475
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
    lsic=.false.
    lc15=.false.
    tdep=2.0
    tvac=1.1
    tint=1.1
    plmin1=0;plmin2=0.;plmin3=0.
    plmax1=0.;plmax2=0.;plmax3=0.
    lrescale=.true.
    ldefcat=.false.
    deltx=0. ; delty=0.0; deltz=0.0
    nbvoisparf(:,:)=0
    write(6,*)'*** analyse du crystal'
    open(175, file='analyse.in')
    read(175,nml=analyse)
    if(.not.(lnbvois).and.(ldesord))lnbvois=.true.
    if(lws) then
       lcomp=.true.
       lvac=.true.
    end if
    if(lsic)lnbvois=.true.
    ! rc(ntyp) est par defaut 2.0 Ang OU vaut rclu dans readdm OU est reprecise ici


    !  if(ldetdec.and.ldecal) then
    !     write(6,*)'ldetdec ET ldecal STOP'
    !     stop
    !  end if

    if ((idistord.gt.0).and.(.not.lperiod)) then
       write(6,*) 'distord seulement avec lperiod =true.'
       stop
    end if
    if ((idistord.ge.3).and.(maxval(pstmax)==0))then
       write(6,*) 'isdistrod=3 preciser pstmax'
       stop
    end if
       

    ! write(6,*)'distordflag',distordflag
    where(rclu(:).gt.0)
       rc(:)=rclu(:)*1.0d-8
    end where

    icall=icall+1
    if(icall==1)then
       open(file='defauts', unit=71)
       open(file='vac', unit=73)
       open(file='int', unit=74)
       open(file='as', unit=75)
       open(file='remp', unit=77)
       open(file='depla', unit=76)
       open(file='structure', unit=72)
       if (lc15)       open(file='vac_et_int', unit=172)
    end if

    if (lcomp) then
       allocate (xpcr(3,imm)) ; allocate (itypcr(imm))
       write(6,*)
       write(6,*)'COMPARAISON crystal it = ' ,itapp
       write(6,*)

       if (lws) then
          write(6,*)'analyse de Wigner-Seitz'
       else
          if(ldeptest) write(6,'(A,F6.1)') 'seuil deplacement pour detection de defauts ',tdep 
          if(lpdep) write(6,'(A,F6.1)') 'ecriture des deplacés ',tdep 
          write(6,'(A,F6.1)') 'seuil lacune ',tvac 
          write(6,'(A,F6.1)') 'seuil interstitiel ',tint 
          tdep=tdep*1.0d-8
          tvac=(tvac*1.0d-8)
          tint=(tint*1.0d-8)
       end if

    end if

    plmin(1)=plmin1; plmax(1)=plmax1
    plmin(2)=plmin2; plmax(2)=plmax2
    plmin(3)=plmin3; plmax(3)=plmax3

    if (any(plmin.ne.0.).or.any(plmax.ne.0.)) then
       call plotpart(xp,plmin,plmax,ityp)
    end if
    if (lcomp) then 
       if (lperiod) call period
       call configcr (xpcr,ityp,lrescale,itypcr)

       !     write(6,*)'apres configcr'


       if (ldecal) then
          if (idecal.gt.0) then
             do ic=1,3
                decal(ic)=xpcr(ic,idecal)-xp(ic,idecal)
             end do
             do i=1,im
                do ic=1,3
                   xpcr(ic,i)=xpcr(ic,i)-decal(ic)
                end do
             end do
          else
             do i=1,im
                xpcr(1,i)=xpcr(1,i)-deltx
                xpcr(2,i)=xpcr(2,i)-delty
                xpcr(3,i)=xpcr(3,i)-deltz
             end do
          end if
          if (lperiod)         call period
          
       end if
       if(ldetdec) then
        
          do i=1,im
             write(490,'(I8,3E16.5)')i,xpcr(1,i)-xp(1,i),xpcr(2,i)-xp(2,i),xpcr(3,i)-xp(3,i)
             write(487,'(I8,3E16.5)')i,xp(1,i),xpcr(1,i),xpcr(1,i)-xp(1,i)
             write(488,'(I8,3E16.5)')i,xp(2,i),xpcr(2,i),xpcr(2,i)-xp(2,i)
             write(489,'(I8,3E16.5)')i,xp(3,i),xpcr(3,i),xpcr(3,i)-xp(3,i)
             deltx=deltx+(xpcr(1,i)-xp(1,i))/im
             delty=delty+(xpcr(2,i)-xp(2,i))/im
             deltz=deltz+(xpcr(3,i)-xp(3,i))/im
          end do
          write(6,*)deltx,delty,deltz
          stop
       end if



       if (lws) then
          call ws
       else
          call depcr (tdep,plmin,plmax,tvac,tint,lvac,lpstruct,lpdef,lpdep, ldeptest)  
       end if
       deallocate (xpcr) ; deallocate (itypcr)
    end if

    if (lnbvois)  call nbvois(xp,ityp,ielat,icall,nbvoisparf,itapp)

    close(175)     

    !  if (dmtype==6)   stop
    !    write(6,*)'fin compcr'

    return
  end subroutine anapos
  !**********************************************************

  subroutine nbvois(xp, ityp,ielat,icall,nbvoisparf,itapp)

    USE T_kind_param_m
    use gen_com_m
    implicit none
    !variables transmises
    integer , intent(in)  :: ielat(imm),nbvoisparf(20,20)
    integer , intent(in) :: ityp(imm)
    real(double) , intent(in) :: xp(3,imm)
    integer :: icall,itapp

    !Variables locales
    integer :: koo,iti,i,ko1,i1,i2,lenfn2,lenfn,j
    integer :: natvi(20),natvityp(20,20),iatvi,latvi,ncelvois
    real(double):: c1,c2,c3,r2,xp1,xp2,xp3
    character*9 :: extension
    character(len=2) :: extension2
    integer, save:: lurasmol

    integer :: maxvois ,nana,itj,i5,iwr
    integer,pointer, save :: nvi(:),nvityp(:,:),ivois(:,:)
    real(double),pointer:: rccar(:)
    real(double) :: xpnp(3,imm),dx(3,20)
    real(double)::cv(1,3)
    logical::lvoisOK,permut
    real(double)::pst,psta,pstp,mtheta(200),costheta,moytet(200,20),thetajik,tetreg,temptri
    integer,pointer::iatd(:)
    integer::natd,k11,j11,nteta,ntetmax,natdtyp(20)

    natd=0
    if (idistord.ge.3)allocate(iatd(im))
    lvoisOK=.true.
    if (lsic) then
       maxvois=8
    end if


    allocate (rccar(ntyp))
    rccar(:)=rc(:)**2
    allocate(nvi(im))
    allocate(nvityp(im,ntyp))
    natvi(:)=0
    !  itapp=it

    !calcul en deux temps
    !calcul du nombre de voisins par atome
    call caltabt 


    if (lperiod) then
       xpnp(:,:)=xp(:,:)
    else 
       call notperiod(xp,xpnp)
    end if

    natvityp(:,:)=0
    ntetmax=-1
    do i = 1, imd
       if(ldefcat.and.ityp(i)==2) cycle 
       pst=0;psta=0;pstp=0.
       nvi(i)=0
       nvityp(i,:)=0
       koo = ielat(i)
       iti=ityp(i)
       ncelvois = min(noxyz,27)-1
       do i1 = 0, ncelvois
          ko1 = ncel(koo,i1)
          !       write(6,*)'koo,ko1',koo,ko1
          do i2 = 1, nato(ko1)
             j = last(i2,ko1)
             !         write(6,*)'j',j
             !         write(6,*)i,xp(:,i)
             !         write(6,*)j,xp(:,j)
             if (i==j) cycle
             c1 = xpnp(1,i)-xpnp(1,j)
             c2 = xpnp(2,i)-xpnp(2,j)
             c3 = xpnp(3,i)-xpnp(3,j)

             if (ldecal_bc .eqv. .true.) then
                if (ibound.gt.1) then		
                   IF (c1>0.5) THEN
                      c3 = c3 - decal_bc
                   ELSE IF (c1 < -0.5) THEN
                      c3 = C3 + decal_bc
                   END IF
                end if
             end if

             c1 = c1+sum(at(1,:)*deltadist(:,i1,koo))
             c2 = c2+sum(at(2,:)*deltadist(:,i1,koo))
             c3 = c3+sum(at(3,:)*deltadist(:,i1,koo))
             if(noxyz==1) then
                cv(1,1) = c1
                cv(1,2) = c2
                cv(1,3) = c3
                call cryst_to_cart (1, cv, bg, -1) !cryst vers cart sur cv
                WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                   cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
                END WHERE
                call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
                c1=cv(1,1)
                c2=cv(1,2)
                c3=cv(1,3)
             end if
             
             r2 = c1*c1+c2*c2+c3*c3
             !        write(6,*)c1,c2,c3,r2,rc(ityp(i))

             if (r2<rccar(ityp(i)))then
                nvi(i)=nvi(i)+1
                nvityp(i,ityp(j))=nvityp(i,ityp(j))+1
                dx(1,nvi(i))=c1
                dx(2,nvi(i))=c2
                dx(3,nvi(i))=c3
             end if


          end do
       end do
       !         write(6,*)i,nvi(i)

       natvi(nvi(i))=natvi(nvi(i))+1
       natvityp(nvi(i),ityp(i))=natvityp(nvi(i),ityp(i))+1



       if (idistord.gt.0) then
          do itj=1,ntyp
             if (nvityp(i,itj).ne.nbvoisparf(ityp(i),itj))then
                lvoisok=.false.
!                write(6,*)'pour i de type mauvais nb de voisins de type j'
!                write(6,*)i,ityp(i), itj, nvityp(i,itj),nbvoisparf(ityp(i),itj)
             end if
          end do


          if (lvoisok.EQV..false.) then
             lvoisok=.true.;cycle
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
                moytet(i5,ityp(i))=moytet(i5,ityp(i))+mtheta(i5)/float(na(ityp(i)))
             end do


          case(2,3)

             iwr=350+ityp(i)
             open(unit=iwr)
             do i5=1,nteta
                read(iwr,*)moytet(i5,ityp(i))
                tetreg=pi/nteta
                psta=psta+(mtheta(i5)-moytet(i5,ityp(i)))**2
                pstp=pstp+(moytet(i5,ityp(i))-tetreg)**2
             end do
             close(iwr)
             pst=sqrt(psta)/sqrt(pstp)
             if( idistord==2) then
                write(380,*)i,ityp(i),pst
             else
                if (pst.gt.pstmax(ityp(i))) then
                   natd=natd+1
                   iatd(natd)=i
                   natdtyp(ityp(i))=natdtyp(ityp(i))+1
                end if
             end if

          case default
          end select


       end if
    end do


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
       open(92,file="distordimag.mol")
       if (natd .NE. 0) then
          write(6,*)natd,' atomes distordus'
          do iti=1,ntyp
             if (natdtyp(iti).ne.0)write(6,*)natdtyp(iti),' atomes distordus de type ',iti
          end do
          write(92,*) natd
          at=at*1.d8
          write (92,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),&
               at(1,3),at(2,3),at(3,3)
          at=at/1.d8
          do i5=1,natd
             write(92,*) ty(ityp(iatd(i5))), xp(:,iatd(i5))*1.0d8
          end do
       else
          write(6,*) 'Pas d atomes distordus'
       end if
       
    end select
    !ouverture des fichiers utiles



    lenfn = 9
    write(extension,'(i9.9)') itapp

134 format(i6)
135 format(A,3f10.4,I7)
136 format(A,3f10.4,D14.5,I7)

200 format(i2)
300 format(i3)
400 format(i4)
500 format(i5)
600 format(i6)
700 format(i7)
800 format(i8)
900 format(i9)
101 format(a1)
201 format(a2)
301 format(a3)
401 format(a4)
501 format(a5)
601 format(a6)
701 format(a7)
801 format(a8)
901 format(a9)



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
       do iatvi=1,20
          if(natvi(iatvi).ne.0) then
             !          write(6,*)it,'Nb d_at. avec',iatvi,'vois. =',natvi(iatvi)!

             lenfn2 = 2
             write(extension2,'(i2.2)') iatvi
             
             latvi=840+iatvi
             !	write(6,*)fnam,extension2
             open(latvi, file=fnam(1:lenfnam)//'.'//extension(1:lenfn)//'.NVI.'//extension2(1:lenfn2)//'.mol', form='formatted', &
                  status='unknown')
             write (latvi, '(I7,A,I7,A,F12.6)') natvi(iatvi), ' IT =', itapp, ' Time = ', timel
             at=at*1.d8
             write (latvi,'(9F10.4)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
             at=at/1.d8 
             do i=1,im
                if (nvi(i)==iatvi)then
                   xp1 = xp(1,i)*1D+08
                   xp2 = xp(2,i)*1D+08
                   xp3 = xp(3,i)*1D+08
                   write(latvi, 138) ty(ityp(i)),xp1, xp2, xp3,nvi(i)
                endif
             end do
          end if

       end do

    case(2)
421    continue
138    format(A,3f10.4,I7,I1)
       lurasmol=840
       open(lurasmol, file=fnam(1:lenfnam)//'.'//extension(1:lenfn)//'.VOIS.mol', form='formatted', &
            status='unknown')

       write (lurasmol, '(I9,A,I7,A,F12.6)') im_glob, ' IT =', itapp, ' Time = ', timel
       at=at*1.d8
       write (lurasmol,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
       at=at/1.d8 

       do i=1,im
          xp1 = xp(1,i)*1D+08
          xp2 = xp(2,i)*1D+08
          xp3 = xp(3,i)*1D+08
          write(lurasmol, 138) ty(ityp(i)),xp1, xp2, xp3,nvi(i)
       end do
    case default
    end select
    if (lsic) then
       if(dmtype==9) then 
          nana=npath
       else
          nana=0
       end if
       call anasic(im,xp,ityp,nvi,ivois,at,bg,itapp,maxvois,eatom,nana)

    end if


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
  subroutine depcr(tdep,plmin,plmax,tvac,tint,lvac,lpstruct,lpdef,lpdep,ldeptest)
    USE T_kind_param_m
    use gen_com_m
    use tab_imm_m
    use tabcr
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double)  :: tdep, plmin(3),plmax(3),tvac,tint
    logical :: lvac,lpstruct,lpdep,lpdef,ldeptest

    !Local variables
    integer :: i,j,k,ndep,ic,nplt,idp,iplt
    integer, dimension(:), pointer:: inddep(:),indplt(:) ! indices des atomes deplaces et plottes

    real(double) :: tdep2
    real(double) :: a1,a2,a3,c1,c2,c3,r2
    real(double), dimension(1,3) :: cv

    integer :: nvac,nint,nremp,nanti,ivac,iint,iremp,ias
    integer, dimension(:), pointer :: indvac,indint,indremp,indas
    logical :: vacfl

    real(double) :: r2min

    real(double) :: c3p,c2p,c1p,c1abs,c2abs,c3abs,r
    integer:: koo,i2,i1,ncelvois,ko1,immin,immax
    !  integer,pointer :: lastcr (:,:),natocr(:),ielatcr(:)

    immin=min(im,imcr)
    immax=max(im,imcr)
    write(6,*)'im imcr ', im,imcr,immin,immax
    ncelvois = min(noxyz,27)-1

    allocate(indplt(imm))

    allocate (lastcr(natperc,0:noxyz))
    allocate(ielatcr(imm))
    allocate(natocr(0:noxyz))

    natocr(:noxyz) = 0
    lastcr(natperc,:noxyz) = 0

    call caltabtcr (natperc,nox,noy,noz,xpcr,imcr,imm,bg,at)
    tdep2=tdep*tdep
    ndep=0
    plmin=1000.0 ; plmax=-1000.0
    ndep=0;nvac=0;nanti=0;nint=0;nas=0; 
    nremp=0 ; nplt=0


    !write(6,*)'ecr pos'
    !  do i=1,im
    !  write(852,*)i
    !  write(852,*)xp(1,i),xp(2,i),xp(3,i)
    !  write(852,*)xpcr(1,i),xpcr(2,i),xpcr(3,i)
    !    enddo

    allocate(inddep(imm))
    if (ldeptest.EQV..true.) then

       call cryst_to_cart (imm, xp, bg, -1)    !cart vers cryst
       call cryst_to_cart (imm, xpcr, bg, -1)    !cart vers cryst

       do i=1,immin
          !     if(ityp(i)==2) cycle
          c1 = xp(1,i)-xpcr(1,i)
          c2 = xp(2,i)-xpcr(2,i)
          c3 = xp(3,i)-xpcr(3,i)
          if (c1>0.5) c1 = c1-1.
          if (c1<(-0.5)) c1 = c1+1.
          if (c2>0.5) c2 = c2-1.
          if (c2<(-0.5)) c2 = c2+1.
          if (c3>0.5) c3 = c3-1.
          if (c3<(-0.5)) c3 = c3+1.
          cv(1,1) = c1
          cv(1,2) = c2
          cv(1,3) = c3
          call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
          r2 = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)
          if(r2.gt.tdep2) then
             ndep=ndep+1
             inddep(ndep)=i
          end if

       end do  !boucle i
       if (immin.lt.immax)then
          do i=immin,immax
             ndep=ndep+1
             inddep(ndep)=i
          end do
       end if

       call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
       call cryst_to_cart (imm, xpcr, at, 1)     !cryst vers cart

       !determination de plmin plmax
       do idp=1,ndep

          do ic=1,3
             if(xp(ic,inddep(idp)).gt.plmax(ic))plmax(ic)=xp(ic,inddep(idp))
             if(xp(ic,inddep(idp)).lt.plmin(ic))plmin(ic)=xp(ic,inddep(idp))
          end do
       end do

       plmin=plmin-1.0d-8 ; plmax=plmax+1.0d-8
       !       if (ndep.gt.0) then
       write(6,*)
       write(6,*)'nombres d atomes deplaces de plus de ',tdep*1.0d8,' = ',ndep
       !       end if


    else
       do i=1,immax
          inddep(i)=i
       end do
       ndep=immax
    end if

    !determination des atomes plottes



    if(lvac) then      ! calcul des lacunes et interstitiels
       !ATTENTION CE CALCUL EST LIMITE AU ATOMES DEPLACES si lpdep=.true. CHOISIR TDEP EN CONSEQUENCE


       allocate(indint(imm))
       allocate(indremp(imm))
       allocate(indvac(imm))
       allocate(indas(imm))






       call cryst_to_cart (imm, xp, bg, -1)    !cart vers cryst
       call cryst_to_cart (imm, xpcr, bg, -1)    !cart vers cryst
       iloop0:do idp=1,ndep
          r2min =10.0
          i=inddep(idp)
          if (i.gt.im)cycle
          if(ldefcat.and.ityp(i)==2) cycle iloop0
          koo = ielat(i)                          ! Numero de la cellule
!          write(6,*)'i idp ',i,idp
          ! pour chaque cel. voisine
          do i1 = 0, ncelvois
             ko1=ncel(koo,i1)
             !              write(6,*)i,idp,koo,i1,ko1,natocr(ko1)
             do i2 = 1, natocr(ko1) !atomes dans la cel dans la conf. init.
                j = lastcr(i2,ko1)

                c1 = xp(1,i)-xpcr(1,j)
                c2 = xp(2,i)-xpcr(2,j)
                c3 = xp(3,i)-xpcr(3,j)
                if (c1>0.5) c1 = c1-1.
                if (c1<(-0.5)) c1 = c1+1.
                if (c2>0.5) c2 = c2-1.
                if (c2<(-0.5)) c2 = c2+1.
                if (c3>0.5) c3 = c3-1.
                if (c3<(-0.5)) c3 = c3+1.
                cv(1,1) = c1
                cv(1,2) = c2
                cv(1,3) = c3
                call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
                r = sqrt(cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3))
!                write(6,*)'j', j,r
                if(r.lt.tint) then ! i est sur le site d'un atome du crystal de depart
                   !                 write(852,*)'remp',i,j,r
                   if(itypcr(j)==ityp(i)) then ! simple remplacement
                      if(lpdep.and.(i.ne.j))then
                         nremp=nremp+1
                         indremp(nremp)=i
                      end if
                   else                      ! antisite
                      nanti=nanti+1
                      indas(nanti)=i
                   end if
                   cycle iloop0 
                end if
             end do
          end do
          !        write(852,*)'INT',i
          nint=nint+1 
          indint(nint)=i

       end do iloop0



       iloop20:  do idp=1,ndep
          r2min =10.0
          i=inddep(idp)
          if (i.gt.imcr)cycle
          if(ldefcat.and.ityp(i)==2) cycle iloop20
          koo = ielatcr(i)                          ! Numero de la cellule
          !         write(6,*)idp,i,koo
          ! pour chaque cel. voisine
          do i1 = 0, ncelvois
             ko1=ncel(koo,i1)
             do i2 = 1, nato(ko1) !atomes dans la cel dans la conf. init.
                j = last(i2,ko1)
                c1 = xpcr(1,i)-xp(1,j)
                c2 = xpcr(2,i)-xp(2,j)
                c3 = xpcr(3,i)-xp(3,j)
                if (c1>0.5) c1 = c1-1.
                if (c1<(-0.5)) c1 = c1+1.
                if (c2>0.5) c2 = c2-1.
                if (c2<(-0.5)) c2 = c2+1.
                if (c3>0.5) c3 = c3-1.
                if (c3<(-0.5)) c3 = c3+1.
                cv(1,1) = c1
                cv(1,2) = c2
                cv(1,3) = c3
                call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
                r = sqrt(cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3))

                if(r.lt.r2min) r2min=r
                if(r2min.lt.tvac) then
                   !                write(852,*)'NOVAC',i,j,r,tvac 
                   cycle iloop20 ! i n'est pas un site de lacune
                end if
             end do
          end do
          !      write(852,*)'VAC',i
          nvac=nvac+1
          indvac(nvac)=i
       end do iloop20

       call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
       call cryst_to_cart (imm, xpcr, at, 1)     !cryst vers cart

       !determination de plmin plmax
       do ivac=1,nvac
          !     write(6,*)plmin(1),plmax(1),xp(1,indvac(ivac))
          do ic=1,3
             if(xp(ic,indvac(ivac)).gt.plmax(ic))plmax(ic)=xp(ic,indvac(ivac))
             if(xp(ic,indvac(ivac)).lt.plmin(ic))plmin(ic)=xp(ic,indvac(ivac))
          end do
       end do
       !          write(6,*)'portion affichee entre ', plmin*1.0d8 ,' et ',plmax*1.0d8
       do iint=1,nint
          do ic=1,3
             if(xp(ic,indint(iint)).gt.plmax(ic))plmax(ic)=xp(ic,indint(iint))
             if(xp(ic,indint(iint)).lt.plmin(ic))plmin(ic)=xp(ic,indint(iint))
          end do
       end do

       write(6,*)'IT = ',it,' nombres de lacunes ',nvac
       write(6,*)'IT = ',it,'nombres d_interstitiels ',nint
       write(6,*)'IT = ',it,'nombres d_antisites ',nanti
       if(lpdep)     write(6,*)'IT = ',it,'nombres de remplacements ',nremp
       if (all(plmin.le.plmax))then

          write(6,'(A,3F7.2,A,3F7.2)')'portion affichee entre ', plmin*1.0d8 ,' et ',plmax*1.0d8
       end if

    end if

    do i=1,im

       if((xp(1,i).gt.plmin(1).and.xp(1,i).lt.plmax(1)).and.&
            & (xp(2,i).gt.plmin(2).and.xp(2,i).lt.plmax(2)).and. &
            & (xp(3,i).gt.plmin(3).and.xp(3,i).lt.plmax(3))) then
          nplt=nplt+1
          indplt(nplt)=i
       end if
    end do
    !----------------------------------------------
    if (lpdep) then

       write(71,*)ndep+2,'IT = ',it,' atomes deplaces de + de ', tdep*1.0d8
       write(71,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3),1d8*at(2,3),&
            1d8*at(3,3)     

       write(71,113)'H', 0. ,0. ,0. 
       write(71,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)

       do idp=1,ndep
          write(71, 113) ty(ityp(inddep(idp))), xp(1,inddep(idp))*1D+8, xp(2,&
               inddep(idp))*1D+8, xp(3,inddep(idp))*1D+8, inddep(idp)
       end do
       write(76,*)ndep+2,'IT = ',it,' atomes deplaces de + de ', tdep*1.0d8
       write(76,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3),1d8*at(2,3),&
            1d8*at(3,3)     

       write(76,113)'H', 0. ,0. ,0. 
       write(76,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)

       do idp=1,ndep
          write(76, 113) ty(ityp(inddep(idp))), xp(1,inddep(idp))*1D+8, xp(2,&
               inddep(idp))*1D+8, xp(3,inddep(idp))*1D+8, inddep(idp)
       end do
    end if
    if (lpdef) then
       !     if (nvac.ne.0) then
       write(71,*)nvac+2,'IT = ',it,' lacunes'
       write(71,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3),1d8*at(2,3),&
            1d8*at(3,3)     

       write(71,113)'H', 0. ,0. ,0. 
       write(71,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)
       do ivac=1,nvac
          write(71, 113) ty(ityp(indvac(ivac))), xpcr(1,indvac(ivac))*1D+8, xpcr(2,&
               indvac(ivac))*1D+8, xpcr(3,indvac(ivac))*1D+8, indvac(ivac)
       end do
       write(73,*)nvac+2,'IT = ',it,' lacunes'
       write(73,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3),1d8*at(2,3),&
            1d8*at(3,3)     

       write(73,113)'H', 0. ,0. ,0. 
       write(73,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)
       do ivac=1,nvac
          write(73, 113) ty(ityp(indvac(ivac))), xpcr(1,indvac(ivac))*1D+8, xpcr(2,&
               indvac(ivac))*1D+8, xpcr(3,indvac(ivac))*1D+8, indvac(ivac)
       end do
       !     end if

       !     if (nint.ne.0) then
       write(71,*)nint+2,'IT = ',it,' interstitiels'
       write(71,'(9F11.5)'),1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3),1d8*at(2,3),&
            1d8*at(3,3)     

       write(71,113)'H', 0. ,0. ,0. 
       write(71,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)
       do iint=1,nint
          write(71, 113) ty(ityp(indint(iint))), xp(1,indint(iint))*1D+8, xp(2,&
               indint(iint))*1D+8, xp(3,indint(iint))*1D+8, indint(iint)
       end do
       write(74,*)nint+2,'IT = ',it,' interstitiels'
       write(74,'(9F11.5)'),1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3),1d8*at(2,3),&
            1d8*at(3,3)     

       write(74,113)'H', 0. ,0. ,0. 
       write(74,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)
       do iint=1,nint
          write(74, 113) ty(ityp(indint(iint))), xp(1,indint(iint))*1D+8, xp(2,&
               indint(iint))*1D+8, xp(3,indint(iint))*1D+8, indint(iint)
       end do

       if (lc15) then

          write(172,*)nint+nvac+2,'IT = ',it,' lacunes et int'
          write(172,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
               1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3),1d8*at(2,3),&
               1d8*at(3,3)     

          write(172,113)'H', 0. ,0. ,0. 
          write(172,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)
          do ivac=1,nvac
             write(172, 113) 'V ', xpcr(1,indvac(ivac))*1D+8, xpcr(2,&
                  indvac(ivac))*1D+8, xpcr(3,indvac(ivac))*1D+8, indvac(ivac)
          end do
          do iint=1,nint
             write(172, 113) 'I ', xp(1,indint(iint))*1D+8, xp(2,&
                  indint(iint))*1D+8, xp(3,indint(iint))*1D+8, indint(iint)
          end do
       end if



       !     end if
       !     if (nanti.ne.0) then
       write(71,*)nanti+2,'IT = ',it,' antisites'
       write(71,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3), & 
            & 1d8*at(2,3),1d8*at(3,3)             

       write(71,113)'H', 0. ,0. ,0. 
       write(71,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)

       do ias=1,nanti
          write(71, 113) ty(ityp(indas(ias))), xp(1,indas(ias))*1D+8, xp(2,&
               indas(ias))*1D+8, xp(3,indas(ias))*1D+8, indas(ias)
       end do
       write(75,*)nanti+2,'IT = ',it,' antisites'
       write(75,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3), & 
            & 1d8*at(2,3),1d8*at(3,3)             

       write(75,113)'H', 0. ,0. ,0. 
       write(75,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)

       do ias=1,nanti
          write(75, 113) ty(ityp(indas(ias))), xp(1,indas(ias))*1D+8, xp(2,&
               indas(ias))*1D+8, xp(3,indas(ias))*1D+8, indas(ias)
       end do
       !     end if

       if(lpdep) then
          write(71,*)nremp+2,'IT = ',it,' remplacements'
          write(71,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),&
               1d8*at(3,1),1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3), &
               1d8*at(2,3),1d8*at(3,3)             

          write(71,113)'H', 0. ,0. ,0. 
          write(71,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)

          do iremp=1,nremp
             write(71, 113) ty(ityp(indremp(iremp))), xp(1,indremp(iremp))*1D+8, xp(2,&
                  indremp(iremp))*1D+8, xp(3,indremp(iremp))*1D+8, indremp(iremp)
          end do
          write(77,*)nremp+2,'IT = ',it,' remplacements'
          write(77,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),&
               1d8*at(3,1),1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3), &
               1d8*at(2,3),1d8*at(3,3)             

          write(77,113)'H', 0. ,0. ,0. 
          write(77,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)

          do iremp=1,nremp
             write(77, 113) ty(ityp(indremp(iremp))), xp(1,indremp(iremp))*1D+8, xp(2,&
                  indremp(iremp))*1D+8, xp(3,indremp(iremp))*1D+8, indremp(iremp)
          end do
       end if


    end if
    if(lpstruct) then
       write(72,*)nplt+2,'IT = ',it,' structure finale'
       write(72,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3),1d8*at(2,3),&
            1d8*at(3,3)             

       write(72,113)'H', 0. ,0. ,0. 
       write(72,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)

       do iplt=1,nplt
          write(72, 113) ty(ityp(indplt(iplt))), xp(1,indplt(iplt))*1D+8, xp(2,&
               indplt(iplt))*1D+8, xp(3,indplt(iplt))*1D+8, indplt(iplt)
       end do
    end if

    deallocate (inddep) 
    deallocate (indplt)
    if(lvac) then
       deallocate (indvac) ; deallocate (indint) ; deallocate (indas) ;deallocate (indremp)
    end if

    !  if (nvac==0) then
    !     write(6,*)'plus de lacunes : STOP'
    !     write(6,*)'it = ', it, ' timel= ', timel
    !     stop
    !  end if

    close(2000); close(2001)

113 format(a2,1x,3(f10.4,1x),1x,1x,i6)

    deallocate (lastcr)
    deallocate(ielatcr)
    deallocate(natocr)



  end subroutine depcr

  !**************** PLOT PART****
  subroutine plotpart(xp,plmin,plmax,ityp)
    USE T_kind_param_m
    use gen_com_m
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double) :: xp(3,imm)
    integer :: ityp(imm)
    real(double) ::  plmin(3),plmax(3)
    !local variables 
    integer :: i,nplt,iplt
    integer, allocatable :: indplt(:) ! indices des atomes  plottes

    allocate(indplt(imm))
    plmin=plmin*1.0d-8 ; plmax=plmax*1.0d-8
    write(6,'(A,3F7.2,A,3F7.2)')'portion affichee entre ', plmin*1.0d8 ,' et ',plmax*1.0d8  
    write(6,*)'portion affichee entre ', plmin*1.0d8 ,' et ',plmax*1.0d8  
    nplt=0
    do i=1,im     
       if((xp(1,i).gt.plmin(1).and.xp(1,i).lt.plmax(1)).and.&
            & (xp(2,i).gt.plmin(2).and.xp(2,i).lt.plmax(2)).and. &
            & (xp(3,i).gt.plmin(3).and.xp(3,i).lt.plmax(3))) then
          nplt=nplt+1
          indplt(nplt)=i
       end if
    end do

    open(file='structurepart', unit=72)

    write(72,*)nplt
    write(72,*)'IT = ',it,' structure partielle'
    do iplt=1,nplt
       write(72, 113) ty(ityp(indplt(iplt))), xp(1,indplt(iplt))*1D+8, xp(2,&
            indplt(iplt))*1D+8, xp(3,indplt(iplt))*1D+8, indplt(iplt)
    end do

    close(72)
113 format(a2,1x,3(f10.4,1x),1x,1x,i6)
  end subroutine plotpart



  subroutine ws
    USE T_kind_param_m
    use gen_com_m
    use tab_imm_m
    use tabcr
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------


    !Local variables
    integer :: i,j,k,ndep,ic,nplt,idp,iplt
    integer, dimension(:), pointer:: indws(:),indatsit(:,:),natsit(:),indint(:),indvac(:),indas(:) ! indices des atomes deplaces et plottes

    real(double) :: a1,a2,a3,c1,c2,c3,r2
    real(double), dimension(1,3) :: cv

    integer :: nvac,nint,nremp,nanti,ivac,iint,iremp,ias
    logical :: vacfl

    real(double) :: r2min

    real(double) :: c3p,c2p,c1p,c1abs,c2abs,c3abs,r
    integer:: koo,i2,i1,ncelvois,ko1,immin,immax
    !  integer,pointer :: lastcr (:,:),natocr(:),ielatcr(:)

    immin=min(im,imcr)
    immax=max(im,imcr)
    write(6,*)'im imcr ', im,imcr,immin,immax
    ncelvois = min(noxyz,27)-1

    allocate(indint(imm))
    allocate(indvac(imm))
    allocate(indas(imm))

    allocate(indws(imm))

    allocate (lastcr(natperc,0:noxyz))
    allocate(ielatcr(imm))
    allocate(natocr(0:noxyz))

    natocr(:noxyz) = 0
    lastcr(natperc,:noxyz) = 0

    call caltabtcr (natperc,nox,noy,noz,xpcr,imcr,imm,bg,at)
    ndep=0;nvac=0;nanti=0;nint=0;nas=0; 


    allocate(indws(imm)) ! inddep note le site j le plus proche de i

    allocate(natsit(imm))
    allocate(indatsit(imm,10))
    natsit=0 ; indws=0 ; indatsit=0

    call cryst_to_cart (imm, xp, bg, -1)    !cart vers cryst
    call cryst_to_cart (imm, xpcr, bg, -1)    !cart vers cryst
    iloop0:do i=1,im 
       r2min =100.0
       koo = ielat(i)                          ! Numero de la cellule
       ! pour chaque cel. voisine
       do i1 = 0, ncelvois
          ko1=ncel(koo,i1)
          !              write(6,*)i,idp,koo,i1,ko1,natocr(ko1)
          do i2 = 1, natocr(ko1) !atomes dans la cel dans la conf. init.
             j = lastcr(i2,ko1)

             c1 = xp(1,i)-xpcr(1,j)
             c2 = xp(2,i)-xpcr(2,j)
             c3 = xp(3,i)-xpcr(3,j)
             if (c1>0.5) c1 = c1-1.
             if (c1<(-0.5)) c1 = c1+1.
             if (c2>0.5) c2 = c2-1.
             if (c2<(-0.5)) c2 = c2+1.
             if (c3>0.5) c3 = c3-1.
             if (c3<(-0.5)) c3 = c3+1.
             cv(1,1) = c1
             cv(1,2) = c2
             cv(1,3) = c3
             call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
             r = sqrt(cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3))
             if(r.lt.r2min) then ! j est pour l'instant le site le plus proche de i
!                write(6,*)r,i,j
                r2min=r
                indws(i)=j
             end if
          end do
       end do

!       if (i.ne.indws(i))write(6,*)'WS',i,indws(i),r2min
       j=indws(i)
       natsit(j)=natsit(j)+1
       indatsit(j,natsit(j))=i


    end do iloop0

    iloop20:  do j=1,imcr
!       write(6,*)j,natsit(j)
       if (natsit(j).eq.0)then
          nvac=nvac+1
          indvac(nvac)=j
       elseif(natsit(j).gt.1) then
          do i=1,natsit(j)
             nint=nint+1
             indint(nint)=indatsit(j,i)
          end do
       else
          if (ityp(j).ne.ityp(indatsit(j,1))) then
             nas=nas+1
             indas(nas)=indatsit(j,1)
          end if
       end if
    end do iloop20

    call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
    call cryst_to_cart (imm, xpcr, at, 1)     !cryst vers cart


    write(6,*)'IT = ',it,' nombres de lacunes ',nvac
    write(6,*)'IT = ',it,'nombres d_interstitiels ',nint
    write(6,*)'IT = ',it,'nombres d_antisites ',nanti



    !----------------------------------------------
    if (lpdef) then
       !     if (nvac.ne.0) then
       write(71,*)nvac+2,'IT = ',it,' lacunes'
       write(71,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3),1d8*at(2,3),&
            1d8*at(3,3)     

       write(71,113)'H', 0. ,0. ,0. 
       write(71,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)
       do ivac=1,nvac
          write(71, 113) ty(ityp(indvac(ivac))), xpcr(1,indvac(ivac))*1D+8, xpcr(2,&
               indvac(ivac))*1D+8, xpcr(3,indvac(ivac))*1D+8, indvac(ivac)
       end do
       write(73,*)nvac+2,'IT = ',it,' lacunes'
       write(73,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3),1d8*at(2,3),&
            1d8*at(3,3)     

       write(73,113)'H', 0. ,0. ,0. 
       write(73,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)
       do ivac=1,nvac
          write(73, 113) ty(ityp(indvac(ivac))), xpcr(1,indvac(ivac))*1D+8, xpcr(2,&
               indvac(ivac))*1D+8, xpcr(3,indvac(ivac))*1D+8, indvac(ivac)
       end do
       !     end if

       !     if (nint.ne.0) then
       write(71,*)nint+2,'IT = ',it,' interstitiels'
       write(71,'(9F11.5)'),1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3),1d8*at(2,3),&
            1d8*at(3,3)     

       write(71,113)'H', 0. ,0. ,0. 
       write(71,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)
       do iint=1,nint
          write(71, 113) ty(ityp(indint(iint))), xp(1,indint(iint))*1D+8, xp(2,&
               indint(iint))*1D+8, xp(3,indint(iint))*1D+8, indint(iint)
       end do
       write(74,*)nint+2,'IT = ',it,' interstitiels'
       write(74,'(9F11.5)'),1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3),1d8*at(2,3),&
            1d8*at(3,3)     

       write(74,113)'H', 0. ,0. ,0. 
       write(74,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)
       do iint=1,nint
          write(74, 113) ty(ityp(indint(iint))), xp(1,indint(iint))*1D+8, xp(2,&
               indint(iint))*1D+8, xp(3,indint(iint))*1D+8, indint(iint)
       end do

       if (lc15) then

          write(172,*)nint+nvac+2,'IT = ',it,' lacunes et int'
          write(172,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
               1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3),1d8*at(2,3),&
               1d8*at(3,3)     

          write(172,113)'H', 0. ,0. ,0. 
          write(172,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)
          do ivac=1,nvac
             write(172, 113) 'V ', xpcr(1,indvac(ivac))*1D+8, xpcr(2,&
                  indvac(ivac))*1D+8, xpcr(3,indvac(ivac))*1D+8, indvac(ivac)
          end do
          do iint=1,nint
             write(172, 113) 'I ', xp(1,indint(iint))*1D+8, xp(2,&
                  indint(iint))*1D+8, xp(3,indint(iint))*1D+8, indint(iint)
          end do
       end if



       !     end if
       !     if (nanti.ne.0) then
       write(71,*)nanti+2,'IT = ',it,' antisites'
       write(71,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3), & 
            & 1d8*at(2,3),1d8*at(3,3)             

       write(71,113)'H', 0. ,0. ,0. 
       write(71,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)

       do ias=1,nanti
          write(71, 113) ty(ityp(indas(ias))), xp(1,indas(ias))*1D+8, xp(2,&
               indas(ias))*1D+8, xp(3,indas(ias))*1D+8, indas(ias)
       end do
       write(75,*)nanti+2,'IT = ',it,' antisites'
       write(75,'(9F11.5)')1d8*at(1,1),1d8*at(2,1),1d8*at(3,1),&
            1d8*at(1,2),1d8*at(2,2),1d8*at(3,2),1d8*at(1,3), & 
            & 1d8*at(2,3),1d8*at(3,3)             

       write(75,113)'H', 0. ,0. ,0. 
       write(75,113)'H', 1d8*zl(1),1d8*zl(2),1d8*zl(3)

       do ias=1,nanti
          write(75, 113) ty(ityp(indas(ias))), xp(1,indas(ias))*1D+8, xp(2,&
               indas(ias))*1D+8, xp(3,indas(ias))*1D+8, indas(ias)
       end do
       !     end if



    end if
    if(lvac) then
       deallocate (indvac) ; deallocate (indint) ; deallocate (indas) 
    end if

    !  if (nvac==0) then
    !     write(6,*)'plus de lacunes : STOP'
    !     write(6,*)'it = ', it, ' timel= ', timel
    !     stop
    !  end if

    close(2000); close(2001)

113 format(a2,1x,3(f10.4,1x),1x,1x,i6)

    deallocate (lastcr)
    deallocate(ielatcr)
    deallocate(natocr)



  end subroutine ws







end module posana
