! *********************************************************************
subroutine deftimestep
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
#if(PARA)
  use mod_mpi
#endif
  !         version paraseq du 21 fevrier 2001
  ! *********************************************************************

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
  integer :: i, iti, ic, expos, imax
  real(double) :: tifac1, tifac2, lts, tseuil, vmax2
  real(double), dimension(imm) :: vpmod2
  real(double) :: tmaxv, tmod, vpmod
  real(double) :: tv1
#if(PARA)
  real(double), dimension(3) :: max_loc,max_glob,max_typ
  integer :: ityp_max
#endif


  !-----------------------------------------------

  ! changement de pas en temps.
  ! le pas en temps optimal est le plus grand tel que
  ! le deplacement maximal entre deux iteration
  ! soit inferieur a 0.5d-10 cm =0.005 A.
  ! ce pas vaut tseuil=1.0d-10/1,0*vmax
  ! pour ne pas tout melanger on ne prend que des pas en temps
  ! egaux a 2.0 ou 5.0 ou 10 * 10 **-qqch


  if (it.le.2) return
  vmax2 = 0.0
  imax = 0
  vpmod2(:im) = vp(1,:im)**2+vp(2,:im)**2+vp(3,:im)**2

  do i = 1, im
     if (vpmod2(i)<=vmax2) cycle
     vmax2 = vpmod2(i)
     imax = i
  end do

#if(PARA)
  max_loc(1)=vmax2
  max_loc(2)=myid
  max_loc(3)=ityp(imax)
  call MPI_ALLREDUCE(max_loc,max_glob,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_WORLD,ierr)
  vmax2 = max_glob(1)
  ityp_max=int(max_loc(3))

!  if (int(max_glob(2)).eq.myid) then
!     call MPI_SEND(ityp(imax),1,MPI_INTEGER,0,1000,MPI_COMM_WORLD,ierr)
!     call MPI_SEND(num_at_glob(imax),1,MPI_INTEGER,0,1001,MPI_COMM_WORLD,ierr)
!  elseif (myid.eq.0) then
!     call MPI_RECV(ityp_max,1,MPI_INTEGER,int(max_glob(2)),1000,MPI_COMM_WORLD,status,ierr)
!     call MPI_RECV(imax,1,MPI_INTEGER,int(max_glob(2)),1001,MPI_COMM_WORLD,status,ierr)
     tmaxv = 1./3./bk*cm(ityp_max)*vmax2
!  endif

!  call MPI_BCAST(tmaxv,1,NDM_MPI_REAL_DOUBLE,0,MPI_COMM_WORLD,ierr)

#else
  tmaxv = 1./3./bk*cm(ityp(imax))*vmax2
#endif
 
  vmax = sqrt(vmax2)
  if (vmax==0)return

  if (rang==0) then
     write (6, '(A,I5,A,D14.5)') '*****  ITERATION  = ', it, '  time = ', &
          timel
     write (6, *) 'Vitesse maximale sur I=', it, imax, vmax, tmaxv
  endif                                      ! fin rang=0

  if (lcasca) then
     !          VPMOD2=VP(1,IKO)**2+VP(2,IKO)**2+VP(3,IKO)**2
     tmod = 1./3./bk*cm(ityp(iko))*vpmod2(iko)
     vpmod = sqrt(vpmod2(iko))
     if (rang==0) write (6, *) 'Vitesse du projectile =', it, iko, vpmod, &
          tmod
  endif

  ! *** Technique anti-bug ! ***
  ! -> tseuil a diminuer pour eviter les derives en energies et temperature
  tseuil = 2.0D-10/(1.0D0*vmax)
  tv1=tstep*vmax !variable servant pour imposer une hysteresis
  !      write(6,*)'tseuil ',tseuil
  lts = log10(tseuil)
  !      write(6,*)'lts ',lts
  expos = 1-int(lts)
  !      write(6,*)'expos ',expos

#ifdef NEC
  ! NEC
  tifac1=exp10(lts+expos)
#else
  ! HP, DEC
  tifac1=10**(lts+expos)
#endif

  if (tifac1<1.0) then
     write (6, *) rang,'sthing wrong deftimestep 1.0'
     call arret_ndm
  else if (tifac1<2.0) then
     tifac2 = float(1)
  else if (tifac1<5.0) then
     tifac2 = float(2)
  else if (tifac1<=10.0) then
     tifac2 = float(5)
  else
     write (6, *) rang,'sthing wrong deftimestep 1.0'
     call arret_ndm
  endif
  !      write(6,*)'tifac2 ',tifac2
  oldtstep = tstep
  tstep = tifac2

  !
  do i = 1, expos
     tstep = tstep/float(10)
  end do


   if (dmtype==4) then

  ! Limite a ne pas depasser pour tstep
     if (tstep<=2D-15.and.(tv1.lt.2D-10.OR.&
          tv1.gt.2.5D-10)) then
        if (tstep/=oldtstep) then
           if (rang==0) then
              write (6, *) '*_*_*_*_ changement de pas en temps *_*_*_'
              write (6, *) ' iteration ', it, 'ancien pas en temps', oldtstep
              write (6, *) 'nouveau tstep ', tstep
           endif 
	end if                                  ! rang=0

     else                                       ! cad si tstep >= 2.10-15s
        tstep = oldtstep
        if (rang==0) write (6, *) 'tstep maintenu'
     endif

  end if



  if (dmtype==1) then

  ! Limite a ne pas depasser pour tstep
     if (tstep<=2D-15.and.(tv1.lt.2D-10.OR.&
          tv1.gt.2.5D-10)) then
        if (tstep/=oldtstep) then
           if (rang==0) then
              write (6, *) '*_*_*_*_ changement de pas en temps *_*_*_'
              write (6, *) ' iteration ', it, 'ancien pas en temps', oldtstep
              write (6, *) 'nouveau tstep ', tstep
           endif                                ! rang=0
           usdh = 1/(two*tstep)
           if (it==0) then
              fp(:,:im) = 0.D0
           endif
           do i = 1, im
              xp(:,i) = xpp(:,i)+tstep*vp(:,i)+tstep**2/cm(ityp(i))/two*fp(:,i)
           end do
           if (lperiod) call period
        endif
     else                                       ! cad si tstep >= 2.10-15s
        tstep = oldtstep
        if (rang==0) write (6, *) 'tstep maintenu',tstep
     endif
  end if

  if (dmtype==2) then
     if (tstep.ne.oldtstep) then
        if (rang==0) then
           write (6, *) '*_*_*_*_ changement de pas en temps *_*_*_'
           write (6, *) ' iteration ', it, 'ancien pas en temps', oldtstep
           write (6, *) 'nouveau tstep ', tstep
        endif                                ! rang=0
        usdh = 1/(two*tstep)
        if (it==0) then
           fp(:,:im) = 0.D0
        endif
        do i = 1, im
           xp(:,i) = xpp(:,i)+tstep*vp(:,i)+tstep**2/cm(ityp(i))/two*fp(:,i)
        end do
         if (lperiod) call period

     else                                       ! cad si tstep >= 2.10-15s
        tstep = oldtstep
        if (rang==0) write (6, *) 'tstep maintenu'
     endif
  end if

  !     write(6,*)'sortie deftimestep'
  return
end subroutine deftimestep
