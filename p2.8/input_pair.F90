
! **********************************************************************
subroutine input_pair
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  !   Version du 3dec. 2001
  ! **********************************************************************
  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, j, k, ic, l,  nprns
  integer::lupotin=95
  real(double) :: ror, pmr, dipr, rof1m, &
       rof2m, rom, dipm, pmm, a_factorm,r8m
  character ::  fnampotin*80
  real(double)::rue
  integer :: ntypr ! nb detype de ce potentiel
  ! lecture des paires
  integer::iti,itp,itt,lw,lw2
  real(double)::qr,cmr,catomr,ecrue,x,rmd,rm2d,xd,fcd
  character :: tyr*3
  integer  :: nb_paire_a_lire, lect_paire,tt1,tt2,igr

  ! lecture des termes a trois corps
  !  logical, dimension (:),pointer :: lue_trip
  integer  :: n3c,npg,npd
  real(double) :: lambr,gamgr,gamdr,agcr,adcr,cangler, dmr,amr,rmr,rumaxa

  ! watanabe
  real(double)::Awatr,Bwatr,pwatr,qwatr,rawatr ! variable de lecture pour pot. watanabe
  integer :: num_3c


  !Stillinger Weber Vashista
  real(double)::capHijlu,capDijlu,capWijlu,c3cr
  integer:: ietaijlu
  !-----------------------------------------------
  !   E x t e r n a l   F u n c t i o n s
  !-----------------------------------------------
  !      real(double) , external :: distmin, calcvol
  !-----------------------------------------------


  namelist /ewald/ rue, alpha, precis, ncouc3, ncoucx, ncoucy, ncoucz,&
       kpmex, kpmey, kpmez, lopt,ecrue

  !EWALD
  !  rumax=0.0
  ecrue=0.0
  r3cm=0.0
  rue = 0.0
  alpha = 0.0
  precis = 0.0
  ncouc3 = 0
  ncoucx = 0
  ncoucy = 0
  ncoucz = 0
  kpmex = 0
  kpmey = 0
  kpmez = 0
  lopt=.FALSE.

  if (rang==0) write (6, *)
  select case (ipotentiel)
  case(0)
     if (rang==0) write (6, *) ' -------- POTENTIEL BMH -----------------'
  case(1)
     if (rang==0) write (6, *) '----------- POTENTIEL BUCKINGHAM --------------'
  case(2)
     if (rang==0) write (6, *) '----------- POTENTIEL WATANABE --------------'
  case(3)
     if (rang==0) write (6, *) '----------- POTENTIEL BUCKINGHAM +R8 --------------'
  case(4)
     if (rang==0) write (6, *) '----------- POTENTIEL UO2 -----------'
  case(6)
     if (rang==0) write (6, *) '----------- POTENTIEL Stillinger Weber a la Vashista -----------'
  case(7)
     if (rang==0) write (6, *) '----------- POTENTIEL PAIRE TABULE -----------'
  case(5)
     if (rang==0) write (6, *) '----------- POTENTIEL BMH+Morse -----------'
  case default
     write (6, *) rang,'Bienvenue dans le cote obscur de la force : pas de potentiel ?'
     call arret_ndm
  end select

  ! -------------------------------------------------------------------
  !    ouverture du fichier potentiel.potin


  select case (ipotentiel)
  case(0)
     fnampotin = 'bmh.potin'
  case(1)
     fnampotin = 'buckingham.potin'
  case(2)
     fnampotin = 'watanabe.potin'
  case(3)
     fnampotin = 'buck8.potin'
  case(4)
     fnampotin = 'uo2.potin'
  case(5)
     fnampotin = 'buckmorse.potin'
  case(6)
     fnampotin = 'SWV.potin'
  case(7)
     fnampotin = 'pair_tab.potin'
  case default
     write (6, *) rang, 'Bienvenue dans le cote obscur de la force :pas de potentiel ?'
     call arret_ndm
  end select


  lupotin = 95
  open(unit=lupotin, file=fnampotin, status='old')

  ! *********** Lecture des donnees de potentiel.potin *************

  !   Lectures communes a tous buckingham et bmh .potin
  select case (ipotentiel)
  case(0,1,3,4,5,7)

     read (lupotin, *) iewald, l3c
#if(PARA)
     if (iewald==2) then
	iewald=1
        write(6,*)'IEWALD MIS A 1'
     endif
#endif

     if (iewald==0) then
        if (rang==0) write (6, *) '-*-*-*-* PAS DE SOMMATION D-EWALD *-*-*-*-'
     elseif (iewald==1) then
        if (rang==0) then
           write (6, *) '-*-*-*-*-* SOMMATION D-EWALD CLASSIQUE *-*-*-*-*-'
           if (npotentiel.gt.1) write(6,*)'FONCTIONNEMENT NON GARANTI!!!'
        end if
     elseif (iewald==2) then
        if (rang==0) write (6, *) '-*-*-*-*-* SOMMATION D-EWALD METHODE PME *-*-*-*-*-'
        if (npotentiel.gt.1) write(6,*)'FONCTIONNEMENT NON GARANTI!!!'
     else
        write (6, *) rang, 'Valeur de iewald erronee : iewald=',iewald
        call arret_ndm
     endif
     if (l3c) then
        if (rang==0) write (6, *) '-*-*-*-*-* TERMES  A 3 CORPS *-*-*-*-*-'
     else
        if (rang==0) write (6, *) '-*-*-*-* PAS DE TERMES  A 3 CORPS *-*-*-*-'
     endif

     ! initialisations de ipo3c

     read (lupotin, nml=ewald)            ! lecture de la namelist ewald

     ! MPI




     ! Fin de la lecture commune aux fichiers .potin


     select case (ipotentiel)
     case(7)
        !ntyp
        if (npotentiel .gt.1)then
           read(lupotin,*) ntypr
           write(6,*)'ntypr pour ce pot',ntypr
        else
           read(lupotin,*) ntyp
           npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
           call  alloc_typ
        end if

        !types
        if (npotentiel .gt.1)then
           if (iewald==0) then
              if (rang==0) write (6, *) 'CM, masse,type, NUMERO DU TYPE D ATOME'
           else
              if (rang==0) write (6, *) 'CHARGE,CM, masse,type, NUMERO DU TYPE D ATOME'
           end if
           do i = 1, ntypr
              if (iewald==0)then
                 read (lupotin,  *) cmr,catomr,tyr,iti
              else
                 read (lupotin,  *) qr,cmr,catomr,tyr,iti
              end if
              if(lue_typ(iti).EQV..true.)then
                 if (rang==0)write(6,*) 'type',iti,'deja lu ; verification de la cohérence'
                 if (cmr*umass.ne.cm(iti))then
                    if (rang==0)write(6,*) 'pb avec cm'
                    call endrun
                 end if
                 if (tyr.ne.ty(iti))then
                    if (rang==0)write(6,*) 'pb avec ty'
                    call endrun
                 end if
                 if(iewald.ne.0)then
                    if (tyr.ne.ty(iti))then
                       if (rang==0)write(6,*) 'pb avec ty'
                       call endrun
                    end if
                 end if
              endif
              if (associated(typ_and_pot))typ_and_pot(iti,ipotentiel)=.true.
              write(6,*)cmr,umass
              cm(iti)=cmr*umass;catom(iti)=catomr;ty(iti)=tyr
              if(iewald.ne.0) q(iti)=qr
              lue_typ(iti)=.true.
              if (rang/=0) cycle
              if (iewald==0)then
                 write (6, '(E12.3,F9.3,A5,I4)') cm(iti),catom(iti),ty(iti),iti
              else
                 write (6, '(2E12.3,F9.3,A5,I4)') q(iti),cm(iti),catom(iti),ty(iti),iti
              end if
           end do
        else
           if (iewald==0) then
              if (rang==0) write (6, *) 'CM, masse,type'
           else
              if (rang==0) write (6, *) 'CHARGE,CM, masse,type'
           end if
           do i = 1, ntyp
              if (iewald==0)then
                 read (lupotin,  *) cm(i),catom(i),ty(i)
              else
                 read (lupotin,  *) q(i),cm(i),catom(i),ty(i)
              end if
              cm(i)=cm(i)*umass
              if (rang/=0) cycle
              write (6, '(I4,2F9.3,A5)') i,cm(i),catom(i),ty(i)

           end do
           rue_pair(:)=rue*1d-8
        end if
        
        !paires

        read(lupotin,*)nb_paire_a_lire,ngr
        allocate (pot_pair_tab(0:ngr,0:4,nb_paire_a_lire))
        allocate (ipo_2_pair_tab(npair))
        if (rang==0) write(6,*)'nb de paires grille',  nb_paire_a_lire, ngr
        do lect_paire=1,nb_paire_a_lire
           read(lupotin,*) tt1,tt2
           l=ipo(tt1,tt2)
           ipo_2_pair_tab(l)=lect_paire
           if(lue_paire(l)) then
              write(6,*) rang,'paire l lue deux fois ', l,tt1,tt2
              call arret_ndm
           endif
           lue_paire(l)=.TRUE. ; typ_pot_pair(l)=ipotentiel       
           rue_pair(l)=rue*1d-8
           if (rang==0) write(6,*)'paire l active  ipotentiel: ',l, ipotentiel
           do igr=1,ngr
              read(lupotin,*)pot_pair_tab(igr,0,lect_paire),pot_pair_tab(igr,1,lect_paire)
           end do

           if (ecrue.ne.0)then
                 rm2d=rue-2*ecrue
                 rmd=rue-ecrue
 !                write(6,*)'rmd',rue,ecrue,rmd
              do igr=1,ngr
                 if (pot_pair_tab(igr,0,lect_paire).ge.rm2d) then
                    xd=10.*(pot_pair_tab(igr,0,lect_paire)-rmd)/ecrue
                    fcd=1./(1.+exp(xd))
                    pot_pair_tab(igr,1,lect_paire)=pot_pair_tab(igr,1,lect_paire)*fcd
!                    write(6,*)'dp',pot_pair_tab(igr,0,lect_paire),xd,fcd
                 end if
              end do
           end if
           pot_pair_tab(0,0,lect_paire)=0.
           pot_pair_tab(0,1,lect_paire)=pot_pair_tab(1,1,lect_paire)
           pot_pair_tab(:,0,lect_paire)=pot_pair_tab(:,0,lect_paire)*1d-8
           pot_pair_tab(:,1,lect_paire)=pot_pair_tab(:,1,lect_paire)*ev2erg

!           do igr=1,ngr
!              lw=340+l
!              lw2=360+l
!              write(lw,'(2D15.6)') 1d8*pot_pair_tab(igr,0,lect_paire),pot_pair_tab(igr,1,lect_paire)
!           end do
           call cspline (ngr,pot_pair_tab(:,0,lect_paire) ,pot_pair_tab(:,1,lect_paire)&
                & ,pot_pair_tab(:,2,lect_paire),pot_pair_tab(:,3,lect_paire),pot_pair_tab(:,4,lect_paire))
!           do igr=1,ngr
!              lw=340+l
!              write(lw2,'(5D15.6)') 1d8*pot_pair_tab(igr,0,lect_paire),pot_pair_tab(igr,1,lect_paire)&
!                   & ,pot_pair_tab(igr,2,lect_paire),pot_pair_tab(igr,3,lect_paire),pot_pair_tab(igr,4,lect_paire)
!           end do
        
   
        end do
        

     case(0)

        ! MPI
        ! nombres de types implicite:
        !        read(ntyp)
        ntyp=10 ; npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2

        call  alloc_typ
        rue_pair(:)=rue*1d-8
        read (lupotin, *) (cm(i),i=1,ntyp)         ! masses
        read (lupotin, *) (catom(i),i=1,ntyp)      ! numeros atomiques
        cm(:ntyp) = cm(:ntyp)*umass

        ! initialisations
        h2sm(:ntyp) = tstep**2/cm(:ntyp)/two
        usdh = 1/(two*tstep)

        if (rang==0) write (6, *) 'type ;charge ; rayon ; bm ; shell ; type'
        do i = 1, ntyp
           read (lupotin, *) q(i), ray(i), bm(i), shel(i),ty(i)
           if (rang/=0) cycle
           write (6, '(I2,4F8.4,a4)') i, q(i), ray(i), bm(i), shel(i),ty(i)
        end do


        ! lecture des caracteristiques des paires
        ! 1. paires standards
        read (lupotin, *) rom, dipm, pmm, rof1m, rof2m
        rom = rom*1.0D-8                        ! conversion A --> cm
        dipm = dipm*1.6021892D-60               ! conversion eV.A^6 --> erg.cm^6
        rof1m = rof1m*1.0D-8                    ! conversion A --> cm
        rof2m = rof2m*1.0D-8                    ! conversion A --> cm
        ro(:npair) = rom
        dip(:npair) = dipm
        pm(:npair) = pmm
        roff1(:npair) = rof1m
        roff2(:npair) = rof2m
        lu_roff_pair(:)=.true.
        ! 2. paires non standards
        read (lupotin, *) nprns
        if (rang==0) write(6,*)'nprns ',nprns
        do i = 1, nprns
           !            if (rang==0) write(6,*) i
           read (lupotin, *) l, ror, dipr, pmr, rof1m, rof2m
           ror = ror*1.0D-8                     ! conversion A --> cm
           dipr = dipr*1.6021892D-60            ! conversion eV.A^6 --> erg.cm^6
           rof1m = rof1m*1.0D-8                 ! conversion A --> cm
           rof2m = rof2m*1.0D-8                 ! conversion A --> cm
           ro(l) = ror
           dip(l) = dipr
           pm(l) = pmr
           roff1(l) = rof1m
           roff2(l) = rof2m

        end do


        ! ++++++++ Fin de la lecture specifique du fichier bmh.potin ++++++++++

        ! ++++++++ Lecture specifique du fichier buckingham.potin ++++++++++++
     case(1,3,5)
        ! MPI
        if (npotentiel .gt.1)then
           read(lupotin,*) ntypr

        else
           read(lupotin,*) ntyp
           npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
           call  alloc_typ
        end if
        if (npotentiel .gt.1)then
           if (rang==0) write (6, *) 'numero, charge, CM, masse,type, NUMERO DU TYPE D ATOME'
           do i = 1, ntypr
              read (lupotin,  *) qr,cmr,catomr,tyr,iti
              if(lue_typ(iti).EQV..true.)then
                 if (rang==0)write(6,*) 'type',iti,'deja lu ; verification de la cohérence'
                 if (cmr*umass.ne.cm(iti))then
                    if (rang==0)write(6,*) 'pb avec cm'
                    call endrun
                 end if
                 if (tyr.ne.ty(iti))then
                    if (rang==0)write(6,*) 'pb avec ty'
                    call endrun
                 end if
              endif
              q(iti)=qr;cm(iti)=cmr*umass;catom(iti)=catomr;ty(iti)=tyr
              lue_typ(iti)=.true.
              if (rang/=0) cycle
              write (6, '(I4,3F9.3,A5)') iti, q(iti),cm(iti),catom(iti),ty(iti)
           end do


        else
           if (rang==0) write (6, *) 'numero, charge, CM, masse,type'
           do i = 1, ntyp
              read (lupotin,  *) q(i),cm(i),catom(i),ty(i)
              if (rang/=0) cycle
              write (6, '(I4,3F9.3,A5)') i, q(i),cm(i),catom(i),ty(i)
              cm(i)=cm(i)*umass
           end do
           rue_pair(:)=rue*1d-8
        end if

        ! initialisations
        if(ipotentiel.ne.(5)) then
           read(lupotin,*)nb_paire_a_lire
           if (rang==0) write(6,*)'nb de paires ',  nb_paire_a_lire
           do lect_paire=1,nb_paire_a_lire
              if (ipotentiel==3) then
                 read(lupotin,*) tt1,tt2, a_factorm, rom, dipm, rof1m, rof2m, r8m
              else
                 read(lupotin,*) tt1,tt2, a_factorm, rom, dipm, rof1m, rof2m
              endif
              l=ipo(tt1,tt2)
              if(lue_paire(l)) then
                 write(6,*) rang,'paire l lue deux fois ', l,tt1,tt2
                 call arret_ndm
              endif
              lue_paire(l)=.TRUE. ; typ_pot_pair(l)=ipotentiel       
              rue_pair(l)=rue*1d-8
              if (rang==0) write(6,*)'paire l active  ipotentiel: ',l, ipotentiel

              !        conversions d'unites
              a_factorm = a_factorm*ecgs           ! conversion eV --> erg
              rom = rom*1.0D-8                     ! conversion A --> cm
              dipm = dipm*1.6021892D-60            ! conversion eV.A^6 --> erg.cm^6
              rof1m = rof1m*1.0D-8                 ! conversion A --> cm
              rof2m = rof2m*1.0D-8                 ! conversion A --> cm

              a_factor(l) = a_factorm
              ro(l) = rom
              dip(l) = dipm
              roff1(l) = rof1m
              roff2(l) = rof2m
              lu_roff_pair(l)=.true.
              if (ipotentiel==3) then
                 r8m=r8m*1.6021892D-60*1.0D-16
                 r8p(l)=r8m
              endif
           end do
        else
           !        if (ipotentiel==5) then ! terme Morse
           Dmorse(:)=0. ; amorse(:)=0. ; remorse(:)=2.0d-8
           read(lupotin,*)nb_paire_a_lire
           if (rang==0) write(6,*)'nb de paires MORSE',  nb_paire_a_lire
           do lect_paire=1,nb_paire_a_lire
              read(lupotin,*) tt1,tt2, dmr,amr,rmr
              l=ipo(tt1,tt2)
              if(lue_paire(l)) then
                 if (rang==0)write(6,*)'paire l lue deux fois ', l,tt1,tt2
                 call arret_ndm
              endif
              lue_paire(l)=.TRUE. ; typ_pot_pair(l)=ipotentiel           
              rue_pair(l)=rue*1d-8
              dmorse(l)=dmr*ecgs
              amorse(l)=amr*1.0d8
              remorse(l)=rmr*1.0d-8
              if (rang==0) write(6,*)'paire l active  ipotentiel: ',l, ipotentiel
           end do
        end if


        ! ++++++++ Fin de la lecture specifique du fichier buckingham.potin ++++++++++

        ! ++++++++ Lecture specifique du fichier uo2.potin ++++++++++++
     case(4) 

        if (npotentiel .gt.1)then
           read(lupotin,*) ntypr

        else
           read(lupotin,*) ntyp
           npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
           call  alloc_typ
        end if
        if (npotentiel .gt.1)then
           if (rang==0) write (6, *) 'numero, charge, CM, masse,type, NUMERO DU TYPE D ATOME'
           do i = 1, ntypr
              read (lupotin,  *) qr,cmr,catomr,tyr,iti
              if(lue_typ(iti).EQV..true.)then
                 if (rang==0)write(6,*) 'type',iti,'deja lu ; verification de la cohérence'
                 call arret_ndm
              end if
              if (cmr.ne.cm(iti))then
                 if (rang==0)write(6,*) 'pb avec cm'
                 stop
              end if
              if (tyr.ne.ty(iti))then
                 if (rang==0)write(6,*) 'pb avec ty'
                 stop
              end if
              q(iti)=qr;cm(iti)=cmr*umass;catom(iti)=catomr;ty(iti)=tyr
              lue_typ(iti)=.true.
              if (rang/=0) cycle
              write (6, '(I4,3F9.3,A5)') iti, q(iti),cm(iti),catom(iti),ty(iti)
           end do

        else
           if (rang==0) write (6, *) 'numero, charge, CM, masse,type'
           do i = 1, ntyp
              read (lupotin,  *) q(i),cm(i),catom(i),ty(i)
              cm(i)=cm(i)*umass
              if (rang/=0) cycle
              write (6, '(I4,3F9.3,A5)') i, q(i),cm(i),catom(i),ty(i)
           end do

        end if

        ro(:)=0.29d-8 ; dip(:)=0.d0 ; a_factor(:)=0.d0
        roff1(:)=0.9d-8 ; roff2(:)=1.0d-8; r8p=0.d0 ;  lu_roff_pair(:)=.true.
        lue_paire(:)=.false.

        read(lupotin,*)nb_paire_a_lire
        if (rang==0) write(6,*)'nb de paires ',  nb_paire_a_lire
        do lect_paire=1,nb_paire_a_lire
           read(lupotin,*) tt1,tt2, a_factorm, rom, dipm, rof1m, rof2m
           l=ipo(tt1,tt2)
           if(lue_paire(l)) then
              write(6,*) rang,'paire l lue deux fois ', l,tt1,tt2
              call arret_ndm
           endif
           lue_paire(l)=.TRUE.
           rue_pair(l)=rue *1d-8   
           !            if (rang==0) write(6,*)'paire l active  : ',l

           !        conversions d'unites
           a_factorm = a_factorm*ecgs    ! conversion eV --> erg
           rom = rom*1.0D-8              ! conversion A --> cm
           dipm = dipm*1.6021892D-60     ! conversion eV.A^6 --> erg.cm^6
           rof1m = rof1m*1.0D-8          ! conversion A --> cm
           rof2m = rof2m*1.0D-8          ! conversion A --> cm

           a_factor(l) = a_factorm
           ro(l) = rom
           dip(l) = dipm
           roff1(l) = rof1m
           roff2(l) = rof2m
           lu_roff_pair(l)=.true.
        end do

        read(lupotin,*) poly5(:)
        read(lupotin,*) poly3(:)

        if (rang==0) write(6,*)'Specificite de l interaction O-O '
        do i = 1, 6
           poly5(i) = poly5(i)*ecgs*(1.d8)**(i-1)
           if (rang==0) write(6,*)'Polynome de degre 5 a',i,' = ', poly5(i)
        enddo
        do i = 1, 4
           poly3(i) = poly3(i)*ecgs*(1.d8)**(i-1)
           if (rang==0) write(6,*)'Polynome de degre 3 b',i,' = ', poly3(i)
        enddo

        read(lupotin,*) rbp5, rp5p3, rp3c
        rbp5   = rbp5*1.0D-8     ! conversion A --> cm
        rp5p3 = rp5p3*1.0D-8     ! conversion A --> cm
        rp3c   = rp3c*1.0D-8     ! conversion A --> cm



        ! ++++++++ Fin de lecture du fichier uo2.potin ++++++++++++++

     end select

     rumax=max(rumax,maxval(rue_pair(:)))
     !  Lecture des parametres des forces a 3 corps

     r3cm=0.
     if (l3c) then
        !        allocate(lue_trip(ntrip))
        !        lue_trip(:)=.false.
        l3ctyp(:)=.false.
        l3cpair(:)=.false.
        read(lupotin,*)n3c     !nb de triplets actifs
        if (rang==0) write(6,*)'n3c ',n3c
        do l=1,n3c
           read(lupotin,*)ic,i,j, lambr,gamgr,gamdr,agcr,adcr,cangler
           if (rang==0) write(6,*)l
           gamgr=gamgr*1.0d-8 ; gamdr=gamdr*1.0d-8 ! A -> cm
           adcr=adcr*1.0d-8 ; agcr=agcr*1.0d-8 ! A -> cm
           lambr=lambr*ecgs                    ! eV -> erg
           k=ipo3c(ic,i,j)
           if (rang==0) write(6,*)'k ntrp ', k,ntrip
           if (rang==0) write(6,*) 'triplet lu :',k,ic,i,j
           if (lue_trip(k)) then
              write(6,*) rang,'triplet deja lu ', ic,i,j,k
              call arret_ndm
           endif
           lue_trip(k)=.TRUE.
!           if (rang==0) write(6,*)l
           l3ctyp(ic)=.TRUE. ; l3ctyp(i)=.TRUE. ; l3ctyp(j)=.TRUE.

           npg=ipo(ic,min(i,j)) ; npd=ipo(ic,max(i,j))
           l3cpair(npg)=.TRUE. ; l3cpair(npd)=.TRUE.
           lamb(k)=lambr
           gam(k,npg)=gamgr
           gam(k,npd)=gamdr
           coup3c(k,npg)=agcr
           coup3c(k,npd)=adcr
           cangle(k)=cangler
           !            if (rang==0) write(6,*) 'r3cm=',r3cm,' agcr=',agcr,' adcr=',adcr
           if(agcr.gt.r3cm) r3cm=agcr
           if(adcr.gt.r3cm) r3cm=adcr
           if (rang==0) write(6,*) 'r3cm=',r3cm
        enddo
        c3c(:)=0.
     endif !l3c


  case(2)
     ! ici lecture de watanabe .potin
     ! il y a un terme a 3 corps :
     iewald= 0
     l3c= .TRUE.
     rumaxa=0

     if (rang==0) write(6,*)'-*-*-*-*-* TERMES  A 3 CORPS *-*-*-*-*-'

     if (npotentiel .gt.1)then
        write(6,*)'ipotentiel==2 et Npotentiel> 1 stop'
        call endrun

     else
        read(lupotin,*) ntyp
        npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
        call  alloc_typ
        allocate(lue_trip(ntrip))
     end if


     !     read(lupotin,*)ntyp
     !     npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
     !     call  alloc_typ
     !     
     do  i=1,ntyp   ! ce qui suit ne sert n'intervient pas dans wat.
        q(i)=0.0    ! mais evite des erreurs d'execution
     enddo
     eta= 1.0
     ty(1)='Si ' ; ty(2)='O  '

     read (lupotin, *) (cm(i),i=1,ntyp)         ! masses
     read (lupotin, *) (catom(i),i=1,ntyp)      ! numeros atomiques
     cm(:ntyp) = cm(:ntyp)*umass

     ! initialisations
     !     where(cm.ne.0.0)  h2sm = tstep**2/cm/two
     !     usdh = 1/(two*tstep)


     read(lupotin,*) epswat   ! lu directement en ergs
     read(lupotin,*) sigmawat ! lu en A
     read(lupotin,*) nprns,n3c    ! nbre de paires/triplets effectifs
     !        conversions d'unites
     sigmawat=sigmawat*1.0d-8  ! conversion A --> cm

     ! initialisation des paires
     do l=1,npair
        Awat(l)= 0.0
        Bwat(l)=0.0
        pwat(l)=0.0
        qwat(l)=0.0
        rawat(l)=0.0
        rawat2(l) = 0.
     enddo

     ! initialisation des triplets
     do ic=1,ntrip
        lamb(ic)= 0.0
        cangle(ic)=0.0
        do l=1,npair
           gam(ic,l)=0.0
           coup3c(ic,l)=0.0
           coup3c2(ic,l)=0.0
        enddo
     enddo


     ! lecture des caracteristiques des paires effectives
     read(lupotin,*)
     do l=1,nprns
        read(lupotin,*)num_paire,Awatr,Bwatr,pwatr,qwatr,rawatr
        Awat(num_paire)=Awatr*epswat ! unites reduites -> ergs
        Bwat(num_paire)=Bwatr
        pwat(num_paire)=pwatr
        qwat(num_paire)=qwatr
        rawat(num_paire)=rawatr*sigmawat ! unites reduites -> cm
        rawat2(num_paire) = rawat(num_paire)**2
        if(rawatr.gt.rumaxa) rumaxa=rawatr
     enddo

     ! Lecture des parametres de g pour partie a 2 corps
     read(lupotin,*)gm1,gm2,gm3,gm4,gm5,gR,gD
     gR=gR*sigmawat       ! unites reduites -> cm
     gD=gD*sigmawat       ! unites reduites -> cm

     ! lecture des caracteristiques des inter. a 3 corps
     read(lupotin,*)
     lue_trip(:)=.false.
     l3ctyp(:)=.false.
     l3cpair(:)=.false.

     do l=1,n3c
        read(lupotin,*)ic,i,j,lambr,gamgr,gamdr,agcr,adcr,cangler
        num_3c=ipo3c(ic,i,j)
        if (lue_trip(num_3c)) then
           if (rang==0) write(6,*)'triplet deja lu ', ic,i,j,num_3c
           call arret_ndm
        endif
        lue_trip(num_3c)=.TRUE.
        l3ctyp(ic)=.TRUE. ; l3ctyp(i)=.TRUE. ; l3ctyp(j)=.TRUE.
        npg=ipo(ic,min(i,j)) ; npd=ipo(ic,max(i,j))
        l3cpair(npg)=.TRUE. ; l3cpair(npd)=.TRUE.

        !            read(lupotin,*)num_3c,npg,npd,lambr,gamgr,gamdr,agcr,adcr,cangler
        lamb(num_3c)=lambr*epswat            ! serie de conversions pour
        gam(num_3c,npg)=gamgr*sigmawat       ! travailler avec r en cm et
        gam(num_3c,npd)=gamdr*sigmawat       ! energie en ergs dans calfo3c
        coup3c(num_3c,npg)=agcr*sigmawat
        coup3c(num_3c,npd)=adcr*sigmawat
        coup3c2(num_3c,npg)=(agcr*sigmawat)**2
        coup3c2(num_3c,npd)=(adcr*sigmawat)**2
        cangle(num_3c)=cangler
        if(agcr.gt.r3cm) r3cm=agcr
        if(adcr.gt.r3cm) r3cm=adcr
     enddo
     ! rayon de coupure
     rumaxa=max(rumax,r3cm)
     rumax=max(rumax,rumaxa*sigmawat)       ! U reduite -> cm
     r3cm=r3cm*sigmawat         ! U reduite -> cm
     r3cm2=r3cm**2

     rue_pair(:) = rumax
     if (rang==0) write(6,*)'rayon de coupure max = ', rumax

  case(6)

     if (npotentiel .gt.1)then
        write(6,*)'ipotentiel==2 et Npotentiel> 1 stop'
        call endrun
     end if
     iewald=0
     l3c=.true.

     read(lupotin,*) ntyp
     npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
     call  alloc_typ

     if (rang==0) write (6, *) 'numero, "charge", CM, masse,type'
     do i = 1, ntyp
        read (lupotin,*)q(i), cm(i),catom(i),ty(i)
        if (rang==0) write (6, '(I4,3F9.3,A5)') i, q(i),cm(i),catom(i),ty(i)
     end do

     ! initialisations
     cm(:ntyp) = cm(:ntyp)*umass
     !     h2sm(:ntyp) = tstep**2/cm(:ntyp)/two
     !     usdh = 1/(two*tstep)

     read(lupotin,*)rue, lambda,xsi
     write(6,'(A,3F12.5)')'rue, lambda,xsi',rue, lambda,xsi
     rue=rue*1d-8; lambda=lambda*1d-8; xsi=xsi*1d-8

     rue_pair(:)=rue
     lue_paire(:npair)=.false.

     read(lupotin,*)nb_paire_a_lire
     if (rang==0) write(6,*)'nb de paires ',  nb_paire_a_lire
     do lect_paire=1,nb_paire_a_lire
        read(lupotin,*) tt1,tt2,ietaijlu,capHijlu,capDijlu,capWijlu
        l=ipo(tt1,tt2)
        if(lue_paire(l)) then
           write(6,*) rang,'paire l lue deux fois ', l,tt1,tt2
           call arret_ndm
        endif
        lue_paire(l)=.TRUE.

        !conversions
        ietaij(l)=ietaijlu
        capHij(l)=capHijlu*ev2erg*ang2cm**ietaijlu
        capDij(l)=capDijlu*e2on4pieps0*ang2cm**3
        capWij(l)=capWijlu*ev2erg*ang2cm**6

     end do

     read(lupotin,*) lambr,gamgr,agcr,cangler,c3cr
     gamgr=gamgr*1.0d-8  ! A -> cm
     agcr=agcr*1.0d-8 ! A -> cm
     lambr=lambr*ecgs                    ! eV -> erg

     lamb(:)=lambr
     gam(:,:)=gamgr
     cangle(:)=cangler
     coup3c(:,:)=agcr
     coup3c2(:,:)=agcr**2
     l3ctyp(:)=.true.
     l3cpair(:)=.true.
     c3c(:)=c3cr
     r3cm=agcr

     rumax=max(rumax,rue)     

  case default
     write (6, *) rang, 'Bienvenue dans le cote obscur de la force :pas de potentiel ?'
     call arret_ndm
  end select
!if(associated (typ_and_pot).eqv..false.), i.e. si npotentiel==1 
  if(associated (typ_and_pot).eqv..false.) then
     allocate (typ_and_pot(ntyp,npotmax))
     typ_and_pot(:,:)=.false.
     typ_and_pot(1:ntyp,ipotentiel)=.true.
  end if
     

  ! ********** Fin de lecture des donnees du fichier potentiel.potin ********
  close(lupotin)
  return
end subroutine input_pair
