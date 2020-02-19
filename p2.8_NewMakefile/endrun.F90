module endrun_mod
        use analyse_mod
        use adf_mod
        use spebc_fin_mod
        use desinteg_insert_mod
        use arret_ndm_mod
        use sauvegarde_mod
        use calfo_mod  
        use rdf_mod
        use rasmol_mod
        implicit none
        contains
! ****************************************************************
subroutine endrun
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
#ifdef PARA
  use mod_para
#endif
#if defined ML && defined PARAML
 use time_measure
#endif
  use posana
  USE cfg_module
  use elec_cell, only:  sauveelec
  !       version MPI du 07 f if (associated(eatom)) eatom(:)=0

  ! ****************************************************************

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

  integer :: i,j, n, nAux_real
  CHARACTER(len=100) :: out_file
  REAL(kind(0.d0)), dimension(:,:), allocatable :: aux_real
  CHARACTER(len=20), dimension(:), allocatable :: aux_title
#ifdef PARA
  integer :: iproc
  real(double), allocatable :: xp_loc(:,:),eatom_loc(:)
  integer, allocatable      :: ityp_loc(:)
  integer, allocatable      :: num_at_glob_loc(:)
  integer :: im_loc
  integer :: proc_source
#endif
  !-----------------------------------------------
  !
  !
  !
  
  if (lPkbar) then
     unitP=1.0d-9
     cunitP='kbar'
  else
     unitP=1.0
     cunitP='d/cm2'
  endif

  ! Un dernier calcul des forces pour la route
  IF (iteTemp.GE.0) iteTemp=1
  IF (iteSigma.GE.0) iteSigma=1

  !flag_fin = .true. !*!
  if(ibound.ne.0) Call spebc_fin (.true.) !*!

  CALL calfo


  ! MPI
  if (rang==0) then
     if (lprtfat)then
        open(unit=10, file='xifi.dat', status='unknown')
        write(10,'(i0)') im
        write (10, '(6g20.8)') (xp(1:3,i)*angst,fp(1:3,i)*erg2eV/angst,i=1,im)
        close(10)
     end if
  end if
  if (lprteat)then
     if (rang==0)open(unit=10, file='xiei.dat', status='unknown')
#ifdef PARA
     ! Le processeur maitre recoit les information des autres processeurs pour les ecrire sur fichier
     if (myid==0) then
        ! Copie des tableaux xp,num_at_glob et ityp locaux 
        allocate(xp_loc(3,imm))
        allocate(ityp_loc(imm))
        allocate(num_at_glob_loc(imm))
        allocate(eatom_loc(imm))
        xp_loc = xp
        ityp_loc = ityp
        num_at_glob_loc = num_at_glob
        im_loc = im
        eatom_loc=eatom
        ! Boucle sur les processeurs
        do iproc=0,nprocs-1
           ! Pour le processeur maitre il n'y a rien a faire
           ! reception des donnees des autres processeurs
           if (iproc.ne.0) then
              call MPI_RECV(im,               1,    MPI_INTEGER,      MPI_ANY_SOURCE, 10001, MPI_COMM_WORLD, status, ierr)
              proc_source = status(MPI_SOURCE)
              call MPI_RECV(xp(1:3,1:im),     3*im, NDM_MPI_REAL_DOUBLE, proc_source, 10002, MPI_COMM_WORLD, status, ierr)
              call MPI_RECV(ityp(1:im),       im,   MPI_INTEGER,         proc_source, 10003, MPI_COMM_WORLD, status, ierr)
              call MPI_RECV(num_at_glob(1:im),im,   MPI_INTEGER,         proc_source, 10004, MPI_COMM_WORLD, status, ierr)
              call MPI_RECV(eatom(1:im),       im,  NDM_MPI_REAL_DOUBLE ,proc_source, 10005, MPI_COMM_WORLD, status, ierr)
           endif
           write (10, '(i6,i3,4g20.8)') (num_at_glob(i),ityp(i),(xp(j,i)*angst,j=1,3),eatom(i)*erg2eV,i=1,im)
           write(6,*)'ZERO',iproc
        enddo  ! fin de boucle sur les processeurs
        ! Le processeur maitre recupere ses donnees locales
        xp = xp_loc
        ityp = ityp_loc
        num_at_glob = num_at_glob_loc
        im = im_loc
        eatom=eatom_loc
        deallocate(xp_loc)
        deallocate(ityp_loc)
        deallocate(eatom_loc)

     else ! Les autres processeurs envoient leurs donnees locales
        call MPI_SEND(im,               1,   MPI_INTEGER,        0,10001,MPI_COMM_WORLD,ierr)
        call MPI_SEND(xp(1:3,1:im),     3*im,NDM_MPI_REAL_DOUBLE,0,10002,MPI_COMM_WORLD,ierr)
        call MPI_SEND(ityp(1:im),       im,  MPI_INTEGER,        0,10003,MPI_COMM_WORLD,ierr)
        call MPI_SEND(num_at_glob(1:im),im,  MPI_INTEGER,        0,10004,MPI_COMM_WORLD,ierr)
        call MPI_SEND(eatom(1:im),im,  NDM_MPI_REAL_DOUBLE ,     0,10005,MPI_COMM_WORLD,ierr)
        write(6,*)'NZ',myid
     endif

#else
     if (lprteattotm.EQV..true.) then

        if (associated(free))then
           do i=1,im
              if (free(i).EQV..true.) write (10, '(i6,i3,4g20.8)') i,ityp(i),(xp(j,i)*angst,j=1,3), & 
              eatomtotm(i)*erg2eV-eatref(ityp(i))
           end do
        else
           do i=1,im
              write (10,'(i6,i3,4g20.8)') i,ityp(i),(xp(j,i)*angst,j=1,3),eatomtotm(i)*erg2eV-eatref(ityp(i))
           end do
        end if

     else
        if (associated(free))then
           do i=1,im
              if (free(i).EQV..true.) write (10, '(i6,i3,4g20.8)') i,ityp(i),(xp(j,i)*angst,j=1,3),eatom(i)*erg2eV
           end do
        else
           do i=1,im
              write (10, '(i6,i3,4g20.8)') i,ityp(i),(xp(j,i)*angst,j=1,3),eatom(i)*erg2eV
           end do
        end if
     end if
#endif


     close(10)

     if (lposmoy.eqv..true.) then
       open(919, file=fnam(1:lenfnam)//'.MOY.mol', form='formatted', &
            status='unknown')
       write (919, '(I9,A)') im_glob, ' POSMOY et ENERGIES '
       at=at*1.d8
       write (919,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
       at=at/1.d8

       do i=1,im
          write (919, 136) ty(ityp(i)),posmoyx(1,i)*1D+08,posmoyx(2,i)*1D+08,posmoyx(3,i)*1D+08,& 
          eatomtotm(i)*erg2eV-eatref(ityp(i)) ,num_at_glob(i)
       end do
    end if

136 format(A,3f10.4,D14.5,I9)

  end if



  if (rang==0) then
#if defined ML && defined PARAML
     write (6, *) 'ML: neighbours  time',  temps_neigh
     write (6, *) 'ML: energy      time',  temps_energy
     write (6, *) 'ML: force       time',  temps_force
     write (6, *) 'ML: stress      time',  temps_stress
     write (6, *) 'ML: descriptors time',  temps_descripteurs
#endif

     write (6, *)
     write (6, *)

     write (6, *) '####### END OF RUN  ######## = ', it, '  time = ', timel
  endif
#ifdef PARA
  temps_dmloop=MPI_Wtime() - temps_dmloop_deb
#endif
  IF (iteSauv.GE.0) then
     call sauvegarde     ! Modif E. Clouet: sauvegarde seulement si voulu
     if (l2T.and.rang==0) call sauveelec
  end IF

  if (.not.linstantrdf) then
     if (iterdf>=0) call rdf
  endif
  if (.not.linstantfda) then
     if (iteangle>=0) call adf
  endif
  if ((dmtype==2).or.(dmtype==3)) then
       it=0
  end if
  call analyse
  if ((ldesinteg.EQV..true.).and.(itdes==nstepdes))call desinteg_insert
  if (iterasmol.GE.0) call rasmol (it)
  if (.not.parallele.and.iteanapos>=0) call anapos (it)

  ! Ecriture d'un fichier atomeye
  if (itecfg.GE.0) then

        WRITE(out_file,'(2a,i0,a)') fnam(1:lenfnam),'.', it, '.cfg'
        OPEN(file=out_file, unit=60, action='write')

        if (dmtype==17)  CALL redefine_ty()

        IF (lPrtEat.OR.lPrtSigat) THEN       ! Energy and/or stress per atom
                nAux_real=0
                IF (lPrtEat)   nAux_real = nAux_real + 1
                IF (lPrtSigat) nAux_real = nAux_real + 6
                IF (Allocated(aux_real)) DeAllocate(aux_real)
                Allocate(aux_real(nAux_real,1:im))
                IF (Allocated(aux_title)) DeAllocate(aux_title)
                Allocate(aux_title(nAux_real))
                n=0
                IF (lPrtEat) THEN
                        aux_title(n+1)="Energy per atom (eV)"
                        aux_real(n+1,1:im)=Eatom(1:im)*erg2eV
                        n = n+1
                END IF
                IF (lPrtSigat) THEN
                        aux_title(n+1) = 'Stress Sxx'
                        aux_title(n+2) = '       Syy'
                        aux_title(n+3) = '       Szz'
                        aux_title(n+4) = '       Syz'
                        aux_title(n+5) = '       Sxz'
                        aux_title(n+6) = '       Sxy (' // cunitP // ')'
                        aux_real(n+1,1:im) = sigat(1,1,1:im)*unitP
                        aux_real(n+2,1:im) = sigat(2,2,1:im)*unitP
                        aux_real(n+3,1:im) = sigat(3,3,1:im)*unitP
                        aux_real(n+4,1:im) = 0.5d0*( sigat(2,3,1:im) + sigat(3,2,1:im) )*unitP
                        aux_real(n+5,1:im) = 0.5d0*( sigat(1,3,1:im) + sigat(3,1,1:im) )*unitP
                        aux_real(n+6,1:im) = 0.5d0*( sigat(1,2,1:im) + sigat(2,1,1:im) )*unitP
                END IF
                CALL WriteCfg(xp, ityp, im, at, 60, nAux_real=nAux_real, aux_real=aux_real, aux_title=aux_title)
                DEALLOCATE(aux_real, aux_title)
        ELSE
                CALL WriteCfg(xp, ityp, im, at, 60)
        END IF
        CLOSE(60)
  endif
  call arret_ndm


  stop
  return
end subroutine endrun
end module
