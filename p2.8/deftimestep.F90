module deftimestep_mod
  USE arret_ndm_mod,only:arret_ndm
  USE gen_com_m, ONLY:bk,depmaxts,dmtype,iko,it,itetimestep,lcasca,lperiod,oldtstep,&
       &rang,timel,tsmin,tstep,two,usdh,vmax,l2T,lspaceNDM
  use atomconfig, only : atom_config_d
  USE boxconfig,only:box_config,periodbox
        implicit none
        contains
! *********************************************************************
subroutine deftimestep(atcf,box)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  USE var_pot, ONLY:cm
  USE elec_cell, ONLY:  etstep, necycle, necyclemin
  USE arret_ndm_mod,only: arret_ndm

#ifdef PARA
  USE Tpara,only:myidsp,nprocspace,comm_space
#else
  USE Tpara,only:nprocspace
#endif

  implicit none
  class (atom_config_d)::atcf
  type(box_config)::box


    integer :: i, iti, ic, expos, imax,ikoloc
  real(double) :: tifac1, tifac2, lts, tseuil, vmax2,depmaxts2
  real(double), dimension(:),allocatable :: vpmod2
  real(double) :: tmaxv, tmod, vpmod
  real(double) :: tv1
#ifdef PARA
  real(double), dimension(3) :: max_loc
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
  allocate (vpmod2(atcf%imm))
      depmaxts2=depmaxts*1.125
!  if (it.le.2) return
  vmax2 = 0.0
  imax = 0
  vpmod2(:atcf%im) = atcf%vp(1,:atcf%im)**2+atcf%vp(2,:atcf%im)**2+atcf%vp(3,:atcf%im)**2

  do i = 1, atcf%im
     if (vpmod2(i)<=vmax2) cycle
     vmax2 = vpmod2(i)
     imax = i
  end do

#ifdef PARA
if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
     max_loc(1)=vmax2
     max_loc(2)=imax+0.5
     !  max_loc(3)=0.5+ityp(imax)
     call comm_space%maxloc(max_loc)
     vmax2 = max_loc(1)
     ityp_max=atcf%ityp(int(max_loc(2)))
     !  ityp_max=int(max_glob(3))
     
     
     
     tmaxv = 1./3./bk*cm(ityp_max)*vmax2
     !tmaxv=0
  else
     tmaxv = 1./3./bk*cm(atcf%ityp(imax))*vmax2
  end if
#else
  tmaxv = 1./3./bk*cm(atcf%ityp(imax))*vmax2
#endif
 
  vmax = sqrt(vmax2)

  if (itetimestep.ne.1) then
     if (rang==0) then
        write (6, '(A,I5,A,D14.5)') '*****  ITERATION  = ', it, '  time = ', &
             timel
        write (6, *) 'Vitesse maximale sur I=', it, imax, vmax, tmaxv
     endif                                      ! fin rang=0
  end if
  if (vmax==0)return
  if (lcasca) then
#if PARA
  ikoloc=0
  do i=1,atcf%im
     if (atcf%num_at_glob(i)==iko) ikoloc=i
  end do
  
#else
ikoloc=iko
#endif

  if (ikoloc.gt.0) then
     tmod = 1./3./bk*cm(atcf%ityp(ikoloc))*vpmod2(ikoloc)
     vpmod = sqrt(vpmod2(ikoloc))
     if (rang==0) write (6, *) 'Vitesse du projectile =', it, iko, vpmod, &
          tmod
  endif
endif	

  ! *** Technique anti-bug ! ***
  ! -> tseuil a diminuer pour eviter les derives en energies et temperature
  tseuil = depmaxts/(1.0D0*vmax)
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
     if ((tstep<=tsmin).and.((tv1.lt.depmaxts).OR.(&
          tv1.gt.depmaxts2))) then
        if (tstep/=oldtstep) then
           if (rang==0) then
              write (6, *) '*_*_*_*_ changement de pas en temps *_*_*_'
              write (6, *) ' iteration ', it, 'ancien pas en temps', oldtstep
              write (6, *) 'nouveau tstep ', tstep
           endif 
!           if(llangevin.eqv..true.) then 
!              gamlg=gamlang/tstep
!           end if
	end if                                  ! rang=0

     else                                       ! cad si tstep >= 2.10-15s
        tstep = oldtstep
        if (rang==0) write (6, *) 'tstep maintenu',tstep
     endif

  end if



  if (dmtype==1) then

  ! Limite a ne pas depasser pour tstep
     if (tstep<=tsmin.and.(tv1.lt.depmaxts.OR.&
          tv1.gt.depmaxts2)) then
        if (tstep/=oldtstep) then
           if (rang==0) then
              write (6, *) '*_*_*_*_ changement de pas en temps *_*_*_'
              write (6, *) ' iteration ', it, 'ancien pas en temps', oldtstep
              write (6, *) 'nouveau tstep ', tstep
           endif                                ! rang=0
           usdh = 1/(two*tstep)
           if (it==0) then
              atcf%fp(:,:atcf%im) = 0.D0
           endif
           do i = 1, atcf%im
              atcf%xp(:,i) = atcf%xpp(:,i)+tstep*atcf%vp(:,i)+tstep**2/cm(atcf%ityp(i))/two*atcf%fp(:,i)
           end do
           call periodbox (box,atcf)
        endif
     else                                       ! cad si tstep >= 2.10-15s
        tstep = oldtstep
        if (rang==0) write (6, *) 'tstep maintenu',tstep
     endif
  end if

  if (dmtype==2) then
     if (tstep<=tsmin.and.(tv1.lt.depmaxts.OR.&
          tv1.gt.depmaxts2)) then
        if (tstep.ne.oldtstep) then
           if (rang==0) then
              write (6, *) '*_*_*_*_ changement de pas en temps *_*_*_'
              write (6, *) ' iteration ', it, 'ancien pas en temps', oldtstep
              write (6, *) 'nouveau tstep ', tstep
           endif                                ! rang=0
           usdh = 1/(two*tstep)
           if (it==0) then
              atcf%fp(:,:atcf%im) = 0.D0
           endif
           do i = 1,atcf%im
              atcf%xp(:,i) = atcf%xpp(:,i)+tstep*atcf%vp(:,i)+tstep**2/cm(atcf%ityp(i))/two*atcf%fp(:,i)
           end do
           call periodbox (box,atcf)
           
        else                                       ! cad si tstep >= 2.10-15s
           tstep = oldtstep
           if (rang==0) write (6, *) 'tstep maintenu',tstep
        endif
     end if
  end if
    ! electronic timestep
  if (l2T)then
     etstep=tstep/necyclemin
     if (etstep.gt.6d-17)then
        etstep=2d-16
        necycle=int(tstep/etstep)
        write(6,*)'chgt etstep',etstep,necycle
     end if
  end if
  !     write(6,*)'sortie deftimestep'
  return
end subroutine deftimestep
end module deftimestep_mod
