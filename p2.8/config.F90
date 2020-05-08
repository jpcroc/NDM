module config_mod
        USE period_mod,only:period
        USE divid_mod, only:divid
#ifdef PARA
        USE coord_to_cell_mod,only:coord_to_cell
        USE decoupage_mod,only: decoupage
#endif
  USE gen_com_m, ONLY:at,bg,zls2,tstep,oldtstep,tmean,timel,nox,noy,noz,im,imm,&
  &it,itmax,ldesinteg,lperiod,pmean,zl,xpspr,nzl,normat,cell_debx,cell_deby,cell_debz,&
  &cell_finx,cell_finy,cell_finz,low_limit
  USE var_pot, ONLY:alpha,na,ntyp,rumax
  USE recips_mod,only: recips
  use cryst_to_cart_mod,only:cryst_to_cart
  USE dynalloccell,only:deallocateall
 USE arret_ndm_mod,only: arret_ndm
!USE caltabi_mod,only: caltabi

  implicit none
  contains
subroutine config
!********************************************************************
!             CONSTRUCTION DE LA BOITE DE SIMULATION
!********************************************************************

  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  USE tab_imm_m
  USE suivinonpbc
#ifdef PARA
  USE mod_para
#endif

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
  integer :: i, j, k, ia, ib, ic, icell, iti, icintype, icintypemod&
       , lucin, lugin, imcell, la, lb, lc, typmax, typmin, npoin, natyp, typ
  integer:: indpoint1, indpointdes
  integer :: passe, nb_passes,ncore,lenfn2
#ifndef PARA
  integer :: nprocs
#endif
  integer ,     dimension(:),   allocatable :: itypc
  real(double), dimension(:,:), allocatable :: xc
  integer,      dimension(:), allocatable   :: num_at_buff
  integer, dimension(:),allocatable     :: ibuffer
  real(double), dimension(:,:),allocatable    :: buffer
  real(double),dimension(:,:),allocatable :: tmpxc
  character :: extension*2

#ifdef PARA
  integer,      dimension(ntyp)         :: na_loc

  integer  :: pointeur_loc, i_loc
  integer  :: numcell, numproc
  integer  :: i_glob
  integer  :: cellx,celly,cellz
#endif
  real(double) :: rue_init
  real(double) :: rumax_init
  real(double) :: alpha_init
  !      integer , dimension(imm,ntyp) :: fv
  !         integer , dimension(6000,10) :: fv    !Truc_bizarre_jmd
  real(double), dimension(3) :: rr
  real(double) :: tirax, tiray, tiraz, x1, x2, x3, a1, a2, a3, c1, c2, c3, r2,xpici,cpp
  real(double) :: rsep2
  real(double) :: avemass !average mass of atoms
  character :: fnamcin*80, fnamgin*80
  !              real(double) drand
  !              external drand


  !-----------------------------------------------
  !
  !
  !

  !-----------------------------------------------------
  ! READING FROM THE CONFIGURATION FILE
  !---------------------------------------------------
  allocate (ibuffer(imm_glob))
  allocate (buffer(3,imm_glob))

  if (rang==0) then
     write(6,*)
     write(6,*)' *-*-*-*-*-*CONSTRUCTION DE LA BOITE*-*-*-*-*-*-'
     write(6,*)
  endif

  if (igen.ge.1) then

#ifdef PARA
     ! En parallele, la lecture du fichier de position se fait en passes
     !  - la premiere pour lire toutes les positions et determiner le
     !    meilleur equilibrage/decoupage
     !  - la deuxieme pour lire uniquement les positions propres au
     !    processeur
     nb_passes = 2
#else
     nb_passes = 1
#endif

     do passe=1,nb_passes

        if (passe==2) close(lucin)

        if(rang==0)               write(6,*)'********** reading configuration from file********'

        ! open du fichier .cin
        lucin = 93
        if (lrestart) then
           fnamcin = fnam(1:lenfnam)//'.cout'
        else
           fnamcin = fnam(1:lenfnam)//'.cin'
        end if
        open(unit=lucin, file=fnamcin, form='unformatted', status='unknown', err=456)

        read (lucin, err=456) icintype

        if (rang==0) write (6, *) 'config type de fichier .cin : ', icintype
        if (icintype>3.or.icintype<0) then
           write (6, *) rang, 'wrong icintype'
           call arret_ndm
        endif

        icintypemod = mod(icintype,2)

        if (lrestart.and.icintypemod==0) then
           write (6, *) rang, 'not possible to restart from this file'
           call arret_ndm
        endif
        !at(vect123,xyz)
!        if (icintype>=2) then
           read (lucin, err=456) at
           if(dilat(1).ne.0.0)then
              do i=1,3
                 at(i,:)=at(i,:)*dilat(i)
              end do
           end if
           call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
           do ic = 1, 3
              normat(ic) = 0
              normat(ic) = normat(ic)+sum(at(:,ic)**2)
              normat(ic) = sqrt(normat(ic))
              zl(ic) = normat(ic)
              normat(ic)=0
              normat(ic)=sqrt(sum(bg(:,ic)**2))
              nzl(ic)=1/zl(ic)


           end do
           !                    write(6,*)at
           zls2 = zl/2.0

!        else
!           read (lucin, err=456) zl                      !size of the box
!           if(dilat(1).ne.0.0)then
!              zl(:)=zl(:)*dilat(:)
!           end if
!           if (rang==0) write (6, *) 'zl ', zl
!           at(1,1)=zl(1)
!           at(2,2)=zl(2)
!           at(3,3)=zl(3)
!           at(1,2)=zero
!           at(1,3)=zero
!           at(2,1)=zero
!           at(2,3)=zero
!           at(3,1)=zero
!           at(3,2)=zero
!           call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
!           nzl(:)=zl(:)
!
!           zls2 = zl/2.0
!
!        endif                                   !icintype=2

	! Il est important de conserver rumax et alpha identique a
 	! chaque appel a la routine divid, on sauvegarde donc la valeur
 	! initiale pour la remettre en sortie
	rumax_init=rumax
        alpha_init=alpha
         if(rang==0)write(6,*)'call divid 0'
        call divid(0)
	rumax=rumax_init
	alpha = alpha_init

#if defined DECOUP || defined PARA
	if (passe==1) then
#ifdef DECOUP 
	   open(123, file='decoup.dat', status='old')
	   read (123, *) nprocs,ncore
	   close(123)
#endif
 	   call  decoupage(nprocs,ncore)
	   allocate(num_at_buff(imm))
        endif
#endif
#ifdef DECOUP
        ! Dans ce cas, pas la peine d'aller plus loin dans l'initialisation
	return	
#endif

        read (lucin, err=456) im_glob                         !number of atoms in the box
        if (im>imm_glob) then
           if(rang==0)                    write (6, *) 'im > imM', im, imm
           call arret_ndm
        endif
#ifdef PARA
	if (passe==1) im = im_glob
#else
 	im = im_glob
#endif

  !                    write(6,*)im
#ifdef PARA
        if (passe == 1) then
           ! lecture muette
           read (lucin, err=456) ibuffer   !ityp
        else
           ! lecture du tableau puis dispatch sur les procs
           read (lucin, err=456) ibuffer   !ityp
	   do i_loc=1,im
	      ityp(i_loc)=ibuffer(num_at_buff(i_loc))
	   enddo
        endif
#else
        !crc 24.11.08        read (lucin, err=456) ityp                       !types
        read (lucin, err=456) ibuffer                       !types
#endif

        if (rang==0) write (6, *) 'types'
        na(:ntyp) = 0


#ifdef PARA
        if (passe.ne.1) then
           do i=1,im
              na(ityp(i))=na(ityp(i))+1
           enddo
           ! On somme les valeurs locales
           na_loc = na
           call MPI_ALLREDUCE(na_loc,na,ntyp,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
	endif
#else
        do i=1,im
           !crc 24.11.08           na(ityp(i))=na(ityp(i))+1
           na(ibuffer(i))=na(ibuffer(i))+1
        enddo
#endif


        !  lecture des positions

#ifdef PARA
        if (passe == 1) then
           ! lecture du tableau et chaque proc conserve ce qui le concerne
           read (lucin, err=456) buffer    ! xp
           if (fmt_cin==1) then
              read (lucin, err=456) ibuffer   ! num_at_glob
           endif

           if((ldesinteg).and.(ides.ne.1))then
              indpoint1=0; indpointdes=0
              do i=1,im_glob
                 if (ibuffer(i)==1) indpoint1=i
                 if (ibuffer(i)==ides) indpointdes=i
              end do
              write(6,*)'point1 pointdes', indpoint1,indpointdes
              ibuffer(indpoint1)=ides
              ibuffer(indpointdes)=1
           end if

           im = 0
           do i=1,im_glob
	      ! Dans les buffer lus, on ne garde que les atomes locaux
              call coord_to_cell(buffer(:,i),numcell)
              numproc=proc_cell(numcell)
              if (numproc == myid) then
                 im = im + 1
                 xp(:,im) = buffer(:,i)
		 if (fmt_cin==0) then
                    num_at_glob(im)=i
		    num_at_buff(im)=i
		 else
		    num_at_glob(im)=ibuffer(i)
      ! num_at_buff permet de stocker l'indice dans le buffer/fichier
      ! du imieme atome local pour repositionner les atomes lors de 
      ! la lecture des autres tableaux
		    num_at_buff(im)=i
		 endif
              endif
           enddo
        else
           ! lecture muette : les coordonnees ont ete lues lors de la premiere passe
           read (lucin, err=456) buffer
           if (fmt_cin==1) then
              read (lucin, err=456) ibuffer   ! num_at_glob
           endif
        endif
#else
        !crc24.11.08        read (lucin, err=456) xp
        read (lucin, err=456) buffer
        if (rang==0)write(6,*)'fmt_cin',fmt_cin
        formcin:select case (fmt_cin)
        case (0) formcin
           do i=1,im
              num_at_glob(i) = i
           enddo
        case(1) formcin
           read (lucin, err=456) num_at_glob
        case default  formcin
           if (rang.eq.0) write(6,*) 'precisez le format fmt_cin'
           call arret_ndm
        end select formcin


#endif

     enddo   ! boucle sur les 2 passes de lecture




     !crc 24.11.08

#ifndef PARA
     if((ldesinteg).and.(ides.ne.1))then
        indpoint1=0; indpointdes=0
        do i=1,im
           if (num_at_glob(i)==1) indpoint1=i
           if (num_at_glob(i)==ides) indpointdes=i
        end do
        write(6,*)'point1 pointdes', indpoint1,indpointdes
        num_at_glob(indpoint1)=ides
        num_at_glob(indpointdes)=1
     end if

     do i=1,im
        !        write(6,*)'i, num_at_glob',i,num_at_glob(i),buffer(1,i)
        ityp(num_at_glob(i))=ibuffer(i)
        xp(:,num_at_glob(i))=buffer(:,i)
     end do


#endif
     !crc 24.11.08




     if (icintypemod==1) then
#ifdef PARA
        read (lucin, err=456) buffer                ! (xpp) former positions
	do i_loc=1,im
	   xpp(:,i_loc)=buffer(:,num_at_buff(i_loc))
	enddo

        read (lucin, err=456) buffer                ! (vp) velocities
        do i_loc=1,im
           vp(:,i_loc)=buffer(:,num_at_buff(i_loc))
        enddo

        read (lucin, err=456) buffer                ! (ax) original positions
        do i_loc=1,im
           ax(:,i_loc)=buffer(:,num_at_buff(i_loc))
        enddo

#else
        !crc 24.11.08
        !        read (lucin, err=456) xpp                     !former positions
        !        read (lucin, err=456) vp                      !velocities
        !        read (lucin, err=456) ax                      !original positions

        read (lucin, err=456) buffer                     !former positions
        do i=1,im
           xpp(:,num_at_glob(i))=buffer(:,i)
        end do

        read (lucin, err=456) buffer                     !former positions
        do i=1,im
           vp(:,num_at_glob(i))=buffer(:,i)
        end do

        read (lucin, err=456) buffer                     !former positions
        do i=1,im
           ax(:,num_at_glob(i))=buffer(:,i)
        end do

        do i=1,im
           num_at_glob(i)=i
        end do
        if((ldesinteg).and.(xpspr(1)==-1000))xpspr(:)=xp(:,1)

#endif
        !            lvpread = .TRUE.
        read (lucin, err=456) oldtstep
        !  Si l'option de redemarrage (lrestart) n'est pas activee
        !  alors les positions d'origine ax deviennent les xp du fichier .cin
        if (.not.lrestart) then
           ax(:,:im) = xp(:,:im)
	   if (lsuivinonpbc) axnonpbc(:,:im)=ax(:,:im)
        endif

     else                                    ! si icintypemod=0
        ax(:,:im) = xp(:,:im)
        if (lsuivinonpbc) axnonpbc(:,:im)=ax(:,:im)
        lvpread=.false.
     endif

     if (lrestart) then
        read (lucin, err=456) tmean, pmean, it, timel
        if ((lrestart).and.(nitmax.ge.0)) itmax=it+nitmax
        tstep = oldtstep

        if (rang==0) then

           write (6, *) 'restart parameters'
           write (6, *) 'it =', it, ' time =', timel
           write (6, *) 'pmean', pmean, ' tmean =', tmean
           write (6, *) 'tstep', tstep
        endif                                ! fin rang=0
        usdh = 1.0/(two*tstep)
     endif

     close (lucin)



#ifdef PARA
     deallocate(num_at_buff)
#endif
     !     do i=1,im
     !        write(6,*)'i, num_at_glob',i,num_at_glob(i),xp(1,i)
     !        ityp(num_at_glob(i))=ibuffer(i)
     !        xp(:,num_at_glob(i))=buffer(:,i)
     !     end do

     !-----------------------------------------------------
     ! BUILDING OF THE CRISTAL FROM .GIN FILE
     !-----------------------------------------------------
  else    ! igen.eq.0


     lvpread=.false.


     ! open fichier .gin
     lugin = 92
     fnamgin = fnam(1:lenfnam)//'.gin'
     open(unit=lugin, file=fnamgin, status='unknown')

     !  si coordonnees reduites
     if ( .not.lalea) then

        if (rang==0) write (6, *) '**********construction du reseau************'

        !                                                !number of cells in 3 directions
        read (lugin, *) lat(1), lat(2), lat(3)
        if (rang==0) write (6, *) 'repetition de mailles', lat


        !     **** coordonnes des vecteurs de maille en A dans une base orthonormee ****
        !                                                !a
        read (lugin, *) at(1,1), at(2,1), at(3,1)
        !b
        read (lugin, *) at(1,2), at(2,2), at(3,2)
        !                                                !c
        read (lugin, *) at(1,3), at(2,3), at(3,3)


        do ic = 1, 3
           normat(ic) = 0
           at(:,ic) = at(:,ic)*1.0D-8*lat(ic)
           normat(ic) = normat(ic)+sum(at(:,ic)**2)
           normat(ic) = sqrt(normat(ic))
           zl(ic) = normat(ic)
           zls2(ic)=zl(ic)/2.
        end do
        la = lat(1)
        lb = lat(2)
        lc = lat(3)
        !        end if

	! Il est important de conserver rue et alpha identique a
 	! chaque appel a la routine divid, on sauvegarde donc la valeur
 	! initiale pour la remettre en sortie
	rumax_init=rumax
	alpha_init=alpha
        call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
	do ic=1,3
           normat(ic)=sqrt(sum(bg(:,ic)**2))
           nzl(ic)=1.0/normat(ic)

           !           if(rang==0)  write(6,)'nzl',nzl(ic)*1d8
        enddo

        call divid(0)
	rumax = rumax_init
	alpha = alpha_init

#ifdef DECOUP
	open(123, file='decoup.dat', status='old')
	   read (123, *) nprocs,ncore
	close(123)
#endif
#if defined DECOUP || defined PARA
 	   call  decoupage(nprocs,ncore)
#endif
#ifdef DECOUP
 ! Dans ce cas, pas la peine d'aller plus loin dans l'initialisation
	return	
#endif
        if (rang==0)	write (6,*) ' Lecture imcell' 
        read (lugin, *) imcell               !number of atoms in UC
        if (rang==0) write (6, *) 'atomes par maille:', imcell
        if (imcell>imm_glob) then
           if(rang==0)               write (6, *) 'trop d_atomes dans la cel. unite'
           call arret_ndm
        endif
	allocate(xc(imcell,3))
	allocate(itypc(imcell))
        if (lsuivinonpbc) ALLOCATE(tmpxc(imcell,3))

        im_glob=la*lb*lc*imcell
        if (im_glob>imm_glob) then
           write (6, *) rang,'imm trop petit'
           call arret_ndm
        endif
        if (rang==0) print *,'Nombre de mailles creees : ',la,lb,lc,la*lb*lc

        do i = 1, imcell
           !            write(6,*)i,imcell
           !                                                !coordonnes reduites des atomes
           read (lugin, *) xc(i,1), xc(i,2), xc(i,3), itypc(i)
           !               write(6,*) xc(i,1), xc(i,2), xc(i,3), itypc(i)

        end do
        !            write(6,*) xc
        !        do i=1,imcell
        !           xc(i,1:3)=xc(i,1:3)+0.000135165
        !        end do
        ! im = la*lb*lc*imcell    ! OB doublon avec 20 lignes plus haut


        if (lperiod.EQV..true.)then
           ! NEVER but NEVER rewrite this sequence. In not true for coordinates |x| > 2 
           !           do i=1,imcell
           !              where(xc(i,:).ge.1.0) 
           !                 xc(i,:)=xc(i,:)-1.0
           !              end where
           !              where(xc(i,:).lt.0.0) 
           !                 xc(i,:)=xc(i,:)+1.0
           !              end where
           !           end do
           ! This sequence is coorect:
           if (lsuivinonpbc) then
              do i=1,imcell
                 tmpsuivi (1:3,i) = xc(i,1:3)
              end do
	   end if
           do i=1,imcell
              WHERE ( (xc(i,:).LT.0.d0).OR.(xc(i,:).GE.1.d0) )
                 xc(i,:)  = xc(i,:)  - Dble(Floor(xc(i,:)))
              END WHERE
           end do

        end if

        i  = 0
        im = 0
#ifdef PARA
        i_glob = 0
#endif
        do ia = 1,la
           do ib = 1,lb
              do ic = 1,lc
                 do icell = 1, imcell
                    i  = i + 1
                    im = im + 1
                    xp(1,i) = (xc(icell,1)+float(ia-1))/float(la)
                    xp(2,i) = (xc(icell,2)+float(ib-1))/float(lb)
                    xp(3,i) = (xc(icell,3)+float(ic-1))/float(lc)
#ifdef PARA
                    do k=1,3
                       xpici=xp(k,i)
                       if ( (xpici < 0.d0 ).OR.( xpici >= 1.d0 ) ) then
                          if ( (xpici > -low_limit).and.(xpici<0.d0) ) then
                             xp(k,i)=zero
                          else
                             cpp  = Dble(Floor(xp(ic,i)))
                             xp (k,i) = xpici     - cpp
                          end if
                       end if
                    end do
#endif
		    if (lsuivinonpbc) then
                       xpnonpbc(1,i) = (tmpsuivi(1,icell)+float(ia-1))/float(la)
                       xpnonpbc(2,i) = (tmpsuivi(2,icell)+float(ib-1))/float(lb)
                       xpnonpbc(3,i) = (tmpsuivi(3,icell)+float(ic-1))/float(lc)
		    end if
                    ityp(i) = itypc(icell)
#ifdef PARA
                    i_glob = i_glob + 1
                    num_at_glob(i)=i_glob
                    ! On teste si c'est un atome local pour le prendre
                    ! en compte ou le retirer
                    call coord_to_cellcoord(xp(1,i),xp(2,i),xp(3,i),cellx,celly,cellz)
                    if (cellx<cell_debx.or.cellx>cell_finx .or. &
                         celly<cell_deby.or.celly>cell_finy .or. &
                         cellz<cell_debz.or.cellz>cell_finz) then
                       ! l'atome n'est pas local, on l'elimine du processeur courant
                       i=i-1
                       im=im-1
                    endif
#else
                    num_at_glob(i)=i
#endif
                 end do
              end do
           end do
        end do

        deallocate(xc)
        deallocate(itypc)

        na=0
        do i=1,im
           na(ityp(i))=na(ityp(i))+1
           ax(:,i)=xp(:,i)
        enddo
#ifdef PARA
        ! On somme les valeurs locales
        na_loc = na
        call MPI_ALLREDUCE(na_loc,na,ntyp,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

        call cryst_to_cart (imm, xp, at, 1)  !cryst vers cart
        ax(:,:im) = xp(:,:im)
	if ((lperiod).and.(lsuivinonpbc)) then
           call cryst_to_cart (imm, xpnonpbc, at, 1)  !cryst vers cart
           axnonpbc(:,:im) = xpnonpbc (:,:im)
	end if

 ! génération de verre
     else if ( lalea) then
        ! Cas ou on tire les positions aleatoires


        read (lugin, *) la, lb, lc           !number of cells in 3 directions
        read (lugin, *) rr(1), rr(2), rr(3)  !size of the unit cell

        rr = rr*1.0D-8
        zl(1) = float(la)*rr(1)
        zl(2) = float(lb)*rr(2)
        zl(3) = float(lc)*rr(3)
        zls2(1) = zl(1)*half
        zls2(2) = zl(2)*half
        zls2(3) = zl(3)*half
        at(1,1)=zl(1)
        at(2,2)=zl(2)
        at(3,3)=zl(3)
        at(1,2)=zero
        at(1,3)=zero
        at(2,1)=zero
        at(2,3)=zero
        at(3,1)=zero
        at(3,2)=zero
        call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))

        ! Il est important de conserver rue et alpha identique a
        ! chaque appel a la routine divid, on sauvegarde donc la valeur
        ! initiale pour la remettre en sortie
        rumax_init=rumax
        alpha_init=alpha
        call divid(0)
        rumax=rumax_init
        alpha = alpha_init

#ifdef DECOUP
        open(123, file='decoup.dat', status='old')
        read (123, *) nprocs,ncore
        close(123)
#endif
#if defined DECOUP || defined PARA
 	   call  decoupage(nprocs,ncore)
#endif
#ifdef DECOUP
        ! Dans ce cas, pas la peine d'aller plus loin dans l'initialisation
        return	
#endif

        read (lugin, *) imcell               !number of atoms in UC
        if (rang==0) write (6, *) 'atomes par maille: ', imcell


        if (imcell>imm) then
           if(rang==0)               write (6, *) 'too many atoms in the unit cell'
           call arret_ndm
        endif
        im = la*lb*lc*imcell
        if (im>imm) then
           if(rang==0)               write (6, *) 'im > imm', im, imm
           call arret_ndm
        endif

        na(:ntyp) = 0
        read (lugin, *) (na(i),i=1,ntyp)

        ! --- Debut du tirage aleatoire des positions initiales ---
        typ=0
        natyp=0
        rsep2=rsep*rsep
        i=0
197     continue
        call random_number(tirax)
        call random_number(tiray)
        call random_number(tiraz)
        !            tirax=drand()
        !            tiray=drand()
        !            tiraz=drand()
        x1=tirax*zl(1)-zls2(1)
        x2=tiray*zl(2)-zls2(2)
        x3=tiraz*zl(3)-zls2(3)
        a1=-dsign(zl(1),x1)
        a2=-dsign(zl(2),x2)
        a3=-dsign(zl(3),x3)
        npoin=1
        do 198 j=1,i-1
           if (npoin.eq.0) goto 198
           c1=x1-xp(1,j)
           c2=x2-xp(2,j)
           c3=x3-xp(3,j)
           if (dabs(c1).gt.zls2(1)) c1=c1+a1
           if (dabs(c2).gt.zls2(2)) c2=c2+a2
           if (dabs(c3).gt.zls2(3)) c3=c3+a3
           r2=c1*c1+c2*c2+c3*c3
           if (r2.lt.rsep2) then
              npoin=0
           endif
198        continue
           if (npoin.eq.1) then
              i=i+1
              xp(1,i)=x1
              xp(2,i)=x2
              xp(3,i)=x3
              do while (i.gt.natyp)
                 typ=typ+1
                 natyp=natyp+na(typ)
              end do
              ityp(i)=typ
           endif
           if (i.ne.im) goto 197

           !            fv(:im,:ntyp) = 0

           !            do j = 1, ntyp
           !               where (ityp(:im)==j) fv(:im,j) = 1
           !            end do

           !            do j = 1, ntyp
           !               na(j) = na(j)+sum(fv(:im,j))
           !            end do

           do i=1,im
              na(ityp(i))=na(ityp(i))+1
           enddo

           ax(:,:im) = xp(:,:im)
           if(lsuivinonpbc) axnonpbc(:,:) = ax(:,:)

        endif      !Fin du if general pour  lalea

        ! fin de la construction du cristal

        close(lugin)

     endif

     ! ----------------------------------------------------------
     !  CONDITIONS PERIODIQUES : REMETTRE LES ATOMES DANS BOITE
     ! ----------------------------------------------------------


     if (lperiod.EQV..true.) call period (imm,xp,xpp,ax)

     ! SUMMARY
#ifdef LAMMPS_VERSION

     if((ipotentiel==-10).or.(ipotentiel==-11)) then
         if(rang==0)write(6,*)'write configuration to conf.lmp'
        call config2data (imm,im,xp,ityp,at,ntyp)
     end if
#endif     

     if (rang==0) then

        write(6,*)
        write (6, *) '-------- boite de simulation ------'
        write (6, *) 'nombre d atomes =', im_glob
        !      write(6,*)'taille de la boite ZL ', zl(1),zl(2),zl(3)
        write (6, '(A,3F11.4)') 'taille de la boite ZL ', 1D+08*zl(1), 1D+08*&
             zl(2), 1D+08*zl(3)
        do i=1,3
           write(6,'(A,I2,3F15.6)')'vecteur ',i, (at(ic,i)*1.0d8,ic=1,3)
        end do
        do iti = 1, ntyp
           if (na(iti)==0) cycle
           write (6, *) na(iti), ' atomes de type', iti
        end do
        !           write (6, *) '----------------------------------'

        !      write(6,*)'sortie de config.f'

     endif                                  ! fin rang=0



     if (llangevin.eqv..true.) then
        allocate(Gl(3,imm))
     end if
     deallocate (ibuffer)
     deallocate (buffer)
     write(6,*)



     return

456  print *,'Erreur dans la lecture du fichier .cin, verifier son format&
          & et fmt_cin ATTENTION A BIG_ENDIAN !! SI COMMPILE BIG_ENDIAN NE LIT PLUS QUE CA'


   end subroutine config

   subroutine  coord_to_cellcoord(coordx,coordy,coordz,cellx,celly,cellz)
     !-----------------------------------------------
     !   M o d u l e s
     !-----------------------------------------------
     USE T_kind_param_m, ONLY:  double
     USE gen_com_m, ONLY:

     implicit none

     real(double) :: coordx,coordy,coordz
     integer      :: cellx,celly,cellz

     real(double), dimension(3,1) :: coord_tab
     real(double) :: aux, auy, auz



     coord_tab(1,1)=coordx
     coord_tab(2,1)=coordy
     coord_tab(3,1)=coordz

     aux = coord_tab(1,1)*nox
     auy = coord_tab(2,1)*noy
     auz = coord_tab(3,1)*noz
     cellx = int(aux)+1
     celly = int(auy)+1
     cellz = int(auz)+1


  return
end subroutine coord_to_cellcoord


#ifdef LAMMPS_VERSION

subroutine config2data (imm,im,xp,ityp,at,ntyp)
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:,ONLY : position_conversion_lammps
  USE var_pot, ONLY:q,ipotentiel
  USE Mat_utils_mod,only: Mat_utils
  implicit none
  integer,intent(in)::imm,im,ntyp
  real(double),intent(in)::xp(3,imm),at(3,3)
  integer,intent(in)::ityp(imm)

  logical lrotated,upper
  real(double)::xhi,yhi,zhi,xy,xz,yz,xlo,ylo,zlo
  real(double), dimension(3,3)::at_lammps,passage,passage_inv
  integer::ic,i
  real(double), dimension(3) :: tmp_coord_i,new_tmp_coord_i

!  if ((at(2,1).ne.0).or.(at(3,1).ne.0).or.(at(3,2).ne.0))then
     write(6,*)
     write(6,*)'LAMMPS BOX BUILT FROM NDM'
     write(6,*)'transforms NDM box to lamps shaped box if needed'
     write(6,*)
!     stop
 ! endif

  open(63,file='conf.lmp',status='unknown')
  write(63,*)
  
  call is_upper_triangular(at,upper)
  write(6,*)'upper ? ',upper
  if (.not.upper) then
     call convert_cell (at,at_lammps,passage)
     call matinv(passage, passage_inv)
     do ic=1,3
        write(6,*)passage(:,ic)
     end do
     write(6,*)
     do ic=1,3
        write(6,*)passage_inv(:,ic)
     end do
  else
     at_lammps=at
     passage(:,:)=0
     do ic=1,3
        passage(ic,ic)=1
     end do
     passage_inv(:,:)=passage(:,:)
  end if


  !building  lammps header
  xlo = 0.d0
  ylo = 0.d0
  zlo = 0.d0
  xhi = at_lammps(1,1)/position_conversion_lammps
  yhi = at_lammps(2,2)/position_conversion_lammps
  zhi = at_lammps(3,3)/position_conversion_lammps
  xy = at_lammps(1,2)/position_conversion_lammps
  xz = at_lammps(1,3)/position_conversion_lammps
  yz = at_lammps(2,3)/position_conversion_lammps

!  xhi=at(1,1)/position_conversion_lammps
!  yhi=dsqrt(at(1,2)**2+at(2,2)**2)/position_conversion_lammps
!  zhi=dsqrt(at(1,3)**2+at(2,3)**2+at(3,3)**2)/position_conversion_lammps
!  xz=at(1,1)*at(1,3)/(xhi*zhi*position_conversion_lammps*position_conversion_lammps)
!  xy=at(1,1)*at(1,2)/(xhi*yhi*position_conversion_lammps*position_conversion_lammps)
!  yz=(at(1,2)*at(1,3)+at(2,2)*at(2,3))/(yhi*zhi*position_conversion_lammps*position_conversion_lammps)

  if(IM.lt.10)then
     write(63,"(I1,A)") IM,' atoms'
  elseif((IM.lt.100).and.(IM.gt.10))then
     write(63,"(I2,A)") IM,' atoms'
  elseif((IM.lt.1000).and.(IM.gt.100))then
     write(63,"(I3,A)") IM,' atoms' 
  elseif((IM.lt.10000).and.(IM.gt.1000))then
     write(63,"(I4,A)") IM,' atoms'       
  elseif((IM.lt.100000).and.(IM.gt.10000))then
     write(63,"(I5,A)") IM,' atoms'         
  elseif((IM.lt.1000000).and.(IM.gt.100000))then
     write(63,"(I6,A)") IM,' atoms'         
  elseif((IM.lt.10000000).and.(IM.gt.1000000))then
     write(63,"(I7,A)") IM,' atoms'        
  elseif((IM.lt.100000000).and.(IM.gt.10000000))then
     write(63,"(I8,A)") IM,' atoms'         
  else 
     stop 'ADD FORMAT'
  endif
  write(6,*)
  write(6,*) " a = (xhi-xlo,0,0); b = (xy,yhi-ylo,0); c = (xz,yz,zhi-zlo). "
  write(6,'(f22.16,a,f22.16,a)')xlo,' ',xhi,' xlo xhi'
  write(6,'(f22.16,a,f22.16,a)')ylo,' ',yhi,' ylo yhi'
  write(6,'(f22.16,a,f22.16,a)')zlo,' ',zhi,' zlo zhi'
  write(6,'(f22.16,a,f22.16,a,f22.16,a)')xy,' ',xz,' ',yz,' xy xz yz'
  write(6,*)


  write(63,"(I1,A)") ntyp, ' atom types' 
  write(63,*) 
  write(63,'(f22.16,a,f22.16,a)')xlo,' ',xhi,' xlo xhi'
  write(63,'(f22.16,a,f22.16,a)')ylo,' ',yhi,' ylo yhi'
  write(63,'(f22.16,a,f22.16,a)')zlo,' ',zhi,' zlo zhi'
  write(63,'(f22.16,a,f22.16,a,f22.16,a)')xy,' ',xz,' ',yz,' xy xz yz'
  write(63,*)
  write(63,"(A)")'Atoms'
  write(63,*)
  select case (ipotentiel)
  case(-10)
     do i=1,im
        tmp_coord_i = xp(:,i)
        new_tmp_coord_i = matmul(passage,tmp_coord_i)/position_conversion_lammps
        write(63,'(I8,a,I2,a,f22.15,a,f22.15,a,f22.15)')i,' ',ityp(i),' ',new_tmp_coord_i(1),' '&
             &,new_tmp_coord_i(2),' ',new_tmp_coord_i(3)

!        write(63,"(I8,I6,F21.12,F20.12,F20.12)") i,ityp(i),&
!             & xp(1,i)/position_conversion_lammps,xp(2,i)/position_conversion_lammps,xp(3,i)/position_conversion_lammps
     end do
  case(-11)
     do i=1,im
        tmp_coord_i = xp(:,i)
        new_tmp_coord_i = matmul(passage,tmp_coord_i)/position_conversion_lammps
        write(63,'(I8,a,I2,a,F20.12,a,f22.15,a,f22.15,a,f22.15)')i,' ',ityp(i),' ',q(ityp(i)),&
             &' ',new_tmp_coord_i(1),' ',new_tmp_coord_i(2),' ',new_tmp_coord_i(3)
!        write(63,"(I8,I6,F20.12,F20.12,F20.12,F20.12)") i,ityp(i),&
!             & q(ityp(i)), xp(1,i)/position_conversion_lammps,xp(2,i)/position_conversion_lammps,xp(3,i)/position_conversion_lammps
     end do
  end select
     close (63)
end subroutine config2data

!---------------------------------------------------
subroutine convert_cell(mat_ini,new_mat,transform)
  USE Mat_utils_mod,only: Mat_utils


  implicit none
  real(kind(0.d0)), dimension(3,3), intent(in) :: mat_ini
  real(kind(0.d0)), dimension(3,3), intent(out) :: new_mat, transform

  !internal
  real(kind(0.d0)), dimension(3,3) :: transit_cell,inv_mat_ini
  real(kind(0.d0)), dimension(3) :: A,B,C, Ahat,AxBhat
  real(kind(0.d0)) :: volume
  logical :: upper, right
  integer :: i

  !matrix is already transpose
  transit_cell = mat_ini
  !write(*,*) 'transpose matrix is :',transit_cell

  call is_upper_triangular(transit_cell,upper)
  if (.not.upper) then
     ! rotate bases into triangular matrix
     new_mat(:,:) = 0.d0
     A = transit_cell(:,1)
     B = transit_cell(:,2)
     C = transit_cell(:,3)
     call right_hand_basis(A,B,C,right)

     if (.not.right) then
         write(*,*)"WARNING: your reper is not right handed."
        write(*,*)"WARNING: This is a critical issue. The LAMMPS results are wrong !!!!!"
        stop
     end if

     new_mat(1,1) = norme(A)
     Ahat = A / norme(A)
     AxBhat = cross_product(A, B) / norme(cross_product(A, B))
     new_mat(1,2) = dot_product(B, Ahat)
     new_mat(2,2) = norme(cross_product(Ahat, B))
     new_mat(1,3) = dot_product(C,Ahat)
     new_mat(2,3) = dot_product(C,cross_product(AxBhat, Ahat))
     new_mat(3,3) = abs(dot_product(C, AxBhat))
     !create and save the transformation for coordinates
     !volume = matdet(mat_ini)
     !trans = np.array([np.cross(B, C), np.cross(C, A), np.cross(A, B)])
     !trans = trans / volume
     !coord_transform = np.dot(tri_mat , trans)
     call matinv(mat_ini,inv_mat_ini)
     transform = matmul(new_mat,inv_mat_ini)

  else
     new_mat = mat_ini
     transform(:,:) = 0.d0
     do i=1,3
        transform(i,i) = 1.d0
     enddo

  endif
  return
end subroutine convert_cell

#endif
end module
