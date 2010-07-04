! *******************************************************************
subroutine calcdepla2
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
#if(PARA)
  use mod_mpi
#endif
  !
  !
  !       version du 09 decembre 2003
  ! *******************************************************************

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
  integer :: ndeplatot
  integer , dimension(ntyp) :: ndepla
  integer :: i, iti
  integer , dimension(imm) :: indic
  integer :: lufilm2it, lutampon
  real(double), dimension(ntyp) :: dr2
  real(double) :: dri2, a1, a2, a3, c1, c2, c3, racdri2
  real(double), dimension(1,3) :: cv

  character :: fnamtampon*10, fnamfilm2it*20, extension*10
#if(PARA)
  integer,      allocatable :: ityp_depla(:)
  integer,      allocatable :: indic_depla(:)
  real(double), allocatable :: xp_depla(:,:)
  real(double), dimension(ntyp) :: dr2_glob
  integer , dimension(ntyp) :: ndepla_glob
  integer :: ndeplatot_glob
  integer :: ndeplatot_tmp
  integer :: proc_source
#endif

  !-----------------------------------------------
  !
  ! local variables
  !
  !
  !

  lufilm2it = 169                            ! index fichiers positions pour une iteration


  !      write(6,*)'entree dans calcdepla2'
  if(rang==0)      write (6, *)
  if(rang==0)      write (6, *) '----------- Deplacements (2) ----------------'
  !       write(6,*)'tdepla2',tdepla2
  ndeplatot = 0
  dr2(:ntyp) = 0.0
  ndepla(:ntyp) = 0
     call cryst_to_cart (imm, xp, bg, -1)    !cart vers cryst
     call cryst_to_cart (imm, ax, bg, -1)    !cart vers cryst

     do i = 1, im
        c1 = ax(1,i)-xp(1,i)
        c2 = ax(2,i)-xp(2,i)
        c3 = ax(3,i)-xp(3,i)
        !            if (c1>0.5) c1 = c1-1.               ! conditions periodiques
        !            if (c1<(-0.5)) c1 = c1+1.
        !            if (c2>0.5) c2 = c2-1.
        !            if (c2<(-0.5)) c2 = c2+1.
        !            if (c3>0.5) c3 = c3-1.
        !            if (c3<(-0.5)) c3 = c3+1.
        cv(1,1) = c1
        cv(1,2) = c2
        cv(1,3) = c3
        call cryst_to_cart (1, cv, at, 1)    !cryst vers cart sur cv
        dri2 = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)
        dr2(ityp(i)) = dr2(ityp(i))+dri2/nad(ityp(i))
        racdri2 = sqrt(dri2)
        if (racdri2<tdepla2) cycle
        ndepla(ityp(i)) = ndepla(ityp(i))+1
        ndeplatot = ndeplatot+1
        indic(ndeplatot) = i
     end do
     call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
     call cryst_to_cart (imm, ax, at, 1)     !cryst vers cart



#if(PARA)
  call MPI_ALLREDUCE(dr2,dr2_glob,ntyp,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  dr2 = dr2_glob
  call MPI_ALLREDUCE(ndepla,ndepla_glob,ntyp,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  ndepla = ndepla_glob
  call MPI_ALLREDUCE(ndeplatot,ndeplatot_glob,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

  ! Allocation des tableaux d'emission/reception
  if (myid==0) then
     allocate(ityp_depla(ndeplatot_glob))
     allocate(xp_depla(1:3,ndeplatot_glob))
     allocate(indic_depla(ndeplatot_glob))
  else
     allocate(ityp_depla(ndeplatot))
     allocate(xp_depla(1:3,ndeplatot))
     allocate(indic_depla(ndeplatot))
  endif
  ! Preparation des tableaux avec les donnees locales
  do i=1,ndeplatot
     ityp_depla(i)=ityp(indic(i))
     xp_depla(:,i)=xp(:,indic(i))
     indic_depla(i)=num_at_glob(indic(i))
  enddo

  ! Recuperation par l'ensemble des procs des diffÃ©rents deplacements

  if (myid==0) then
     do i=1,nprocs-1
        call MPI_RECV(ndeplatot_tmp,1,MPI_INTEGER,MPI_ANY_SOURCE,14001,MPI_COMM_WORLD,status,ierr)
        if (ndeplatot_tmp.ne.0) then
           proc_source = status(MPI_SOURCE)
           call MPI_RECV(ityp_depla(ndeplatot+1),ndeplatot_tmp,MPI_INTEGER,proc_source,14002,MPI_COMM_WORLD,status,ierr)
           call MPI_RECV(xp_depla(:,ndeplatot+1),3*ndeplatot_tmp,NDM_MPI_REAL_DOUBLE,proc_source,14003,MPI_COMM_WORLD,status,ierr)
           call MPI_RECV(indic_depla(ndeplatot+1),ndeplatot_tmp,MPI_INTEGER,proc_source,14004,MPI_COMM_WORLD,status,ierr)
           ndeplatot = ndeplatot + ndeplatot_tmp
        endif
     enddo
  else
     call MPI_SEND(ndeplatot,1,MPI_INTEGER,0,13001,MPI_COMM_WORLD,ierr)
     if (ndeplatot.ne.0) then
        call MPI_SEND(ityp_depla,ndeplatot,MPI_INTEGER,0,14002,MPI_COMM_WORLD,ierr)
        call MPI_SEND(xp_depla(:,1:ndeplatot),3*ndeplatot,NDM_MPI_REAL_DOUBLE,0,14003,MPI_COMM_WORLD,ierr)
        call MPI_SEND(indic_depla,ndeplatot,MPI_INTEGER,0,14004,MPI_COMM_WORLD,ierr)
     endif
  endif
#endif

  if(rang==0)then
     write (6, '(A,I5,A,D10.3)') '*  ITERATION  = ', it, '  time = ', timel

     write (6, *)
     do iti = 1, ntyp
        if (nad(iti)==0) cycle
        write (6, *) 'deplacement moyen des atomes de type ', iti, ' = ', &
             sqrt(dr2(iti)*1D+16)
        write (6, *) 'nombre d-atomes de type ', iti, ' deplaces = ', ndepla(&
             iti)
     end do
  end if

  ! ***** Ecriture positions atomiques dans fichiers differents *****
  ! ***** 1 iteration -> 1 fichier pour traitement images animees gif *****

  ! conversion entier-->alphanumerique par transfert du nombre
  ! de l'iteration vers fichier tampon relu sous format caractere.

  lutampon = 17
  fnamtampon = 'tampon'
  if (rang==0)then
     open(unit=lutampon, file=fnamtampon, form='formatted', status='unknown')
     if (it<=9) write (17, '(I1)') it
     if (it<=99.and.it>9) write (17, '(I2)') it
     if (it<=999.and.it>99) write (17, '(I3)') it
     if (it<=9999.and.it>999) write (17, '(I4)') it
     if (it<=99999.and.it>9999) write (17, '(I5)') it
     if (it<=999999.and.it>99999) write (17, '(I6)') it
     if (it<=9999999.and.it>99999) write (17, '(I7)') it
     if (it<=99999999.and.it>999999) write (17, '(I8)') it
     if (it<=999999999.and.it>9999999) write (17, '(I9)') it
!    if (it<=9999999999.and.it>99999999) write (17, '(I10)') it
!    if (it>=9999999999) then
     if (it>=999999999) then
        write (6, *) 'probleme de format dans calcdepla2.f'
        stop
     endif
     rewind 17
     read (17, '(A10)') extension


     !    ouverture d'un fichier film2it.(iteration) pour sauvegarde
     !    des positions toutes les itedepla iterations.

     fnamfilm2it = 'film2it.'//extension
     open(unit=lufilm2it, file=fnamfilm2it, form='formatted', status='unknown'&
          )
     ! Ecriture des types et coordonnees des atomes deplaces de plus de
     ! tdepla2 angstroems dans le fichier film2it.extension
     write (lufilm2it, '(I7,A,I7,A,D10.3)')  ndeplatot+2, ' IT =', it, ' Time = ', timel
     at=at*1.d8
     write (lufilm2it,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
     at=at/1.d8
 
     write (lufilm2it, 114)  zls2(1)*1d8,zls2(2)*1d8,zls2(3)*1d8
     write (lufilm2it, 114) -zls2(1)*1d8,-zls2(2)*1d8,-zls2(3)*1d8
114  format('H ',1x,3(f10.4,1x))

     !          write (lufilm2it, *) ndeplatot
     !         write (lufilm2it, *) ' IT', it, ' time ', timel


#if(PARA)
     do i = 1, ndeplatot
        !       write(6,*)i,indic(i)
        write (lufilm2it, 113) ty(ityp_depla(i)), xp_depla(1,i)*1D+8, xp_depla(2,&
             i)*1D+8, xp_depla(3,i)*1D+8, indic_depla(i)
     end do
#else
     do i = 1, ndeplatot
        !       write(6,*)i,indic(i)
        write (lufilm2it, 113) ty(ityp(indic(i))), xp(1,indic(i))*1D+8, xp(2,&
             indic(i))*1D+8, xp(3,indic(i))*1D+8, indic(i)
     end do
#endif
     close(lufilm2it)



     close(17)
  end if
  ! ***** Fin ecriture positions dans plusieurs fichiers *****


113 format(a2,1x,3(f10.4,1x),1x,'!',1x,i6)

  return
end subroutine calcdepla2
