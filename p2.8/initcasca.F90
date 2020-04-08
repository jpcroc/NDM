module initcasca_mod
  USE cryst_to_cart_mod
  USE period_mod
  USE gen_com_m, ONLY:at,bg,depmaxts,dmtype,ecgs,eko,iko,im,im_glob,imm,lderive,lperiod,ltranche,&
       &oldtstep,parallele,rang,tsmin,tstep,two,usdh,vmax,xko,xx0,yko,yy0,zko,zz0
  implicit none
contains
  ! *************** initialisation de la cascade **********************
  subroutine initcasca
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE var_pot, ONLY:cm
    USE tab_imm_m
    USE elec_cell, ONLY : necycle,etstep,necyclemin
    ! *******************************************************************
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
    integer :: i, ic, i1,ikoloc
    real(double) :: z1, z2, z3, t1, t2, t3, znorm, aux1

    integer :: iti, expos, imax
    real(double) :: tifac1, tifac2, lts, tseuil, vmax2
    real(double), dimension(imm) :: vpmod2
    real(double) :: masstot, vpi(3)

    integer :: seed_size,isl
    integer, dimension(:), allocatable :: iseedt

#ifdef PARA
    real(double), dimension(3) :: max_loc,max_glob,max_typ
    real(double), dimension(1) :: max_typl,max_typG
    integer :: ityp_max
#endif

#if PARA
    ikoloc=0
    do i=1,im
       if (num_at_glob(i)==iko) ikoloc=i
    end do

#else
    ikoloc=iko
#endif

    !-----------------------------------------------
    ! --- Translation de l'atome IKO au centre de la boite de simulation ---
    if (ikoloc.gt.0) then
       write (6, *) 'initialisation de la cascade'
       write(6,*)'ATOME ',iko, '  TYPE ',ityp(ikoloc)
       if (iko>im_glob) then
          write (6, *) 'wrong input cascade iko eko ', iko, eko
          stop
       endif
    endif                                      ! rang=0

    if (xx0.ge.0) then
#if PARA
#else

       znorm=sqrt(xx0**2+yy0**2+zz0**2)
       if (znorm==0) then
          xx0=0.25 ; yy0=0.25 ; zz0=0.25
       end if

       call cryst_to_cart (imm, xp, bg, -1)    !cart vers cryst
       call cryst_to_cart (imm, ax, bg, -1)    !cart vers cryst
       !debug write(6,*)xp(1,ikoloc),xx0
       if (ltranche) then
          t1 = 0.
          t2 = 0.
          t3 = 0.
       else
          t1 = xp(1,ikoloc)-xx0
          t2 = xp(2,ikoloc)-yy0
          t3 = xp(3,ikoloc)-zz0
       endif

       !      write(6,*)'t1 t2 t3 zl', t1,t2,t3,zl(1),zl(2),zl(3)
       xpp(1,:im) = xpp(1,:im)-t1
       xpp(2,:im) = xpp(2,:im)-t2
       xpp(3,:im) = xpp(3,:im)-t3
       xp(1,:im) = xp(1,:im)-t1
       xp(2,:im) = xp(2,:im)-t2
       xp(3,:im) = xp(3,:im)-t3
       ax(1,:im) = ax(1,:im)-t1
       ax(2,:im) = ax(2,:im)-t2
       ax(3,:im) = ax(3,:im)-t3

       call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
       call cryst_to_cart (imm, ax, at, 1)     !cryst vers cart
       if (lperiod)       call period 
#endif

       !  write(6,*)xp(1,iko)
       !                                                !Conditions periodiques

       ! --- Fin de la translation ---
    end if 	! xx0>0
    znorm = sqrt(xko**2+yko**2+zko**2)

    if (znorm==0)then
       if (parallele) then
          write(6,*) 'tirage al�atoire projectile pas programm�'
          stop
       endif
       write(6,*) 'tirage al�atoire xko '
       call random_seed(size=seed_size)
       allocate(iseedt(seed_size))
       iseedt = 0
       !        if (iseed==0) then
       call system_clock (count=isl)
       write(6,*)'ISLxko',isl
       iseedt(1)=isl
       !        else
       !           iseedt(1)=iseed
       !        end if

       call    random_seed (put=iseedt)
       deallocate(iseedt)
       !        do i=1,100
       call random_number(xko)
       call random_number(yko)
       call random_number(zko)
       znorm = sqrt(xko**2+yko**2+zko**2)
    end if

    z1 = xko/znorm
    z2 = yko/znorm
    z3 = zko/znorm

    ! MPI
    if (rang==0) then
       write (6, 576) z1, z2, z3
576    format('Direction du projectile ',3(f8.4,1x))
    endif

577 format('Positions initiales du projectile (A) ',3(f8.4,1x))
#if PARA
    ikoloc=0
    do i=1,im
       if (num_at_glob(i)==iko) ikoloc=i
    end do

    if (ikoloc.gt.0) then
       write(6,*)'IKO dans processeur', rang,ikoloc
       write (6, 577) xp(1:3,ikoloc)*1.0d8     

       aux1 = sqrt(eko*ecgs*2./cm(ityp(ikoloc)))    !Vitesse en cgs
       vp(1,ikoloc) = vp(1,ikoloc)+z1*aux1
       vp(2,ikoloc) = vp(2,ikoloc)+z2*aux1
       vp(3,ikoloc) = vp(3,ikoloc)+z3*aux1
575    continue
       xpp(1,ikoloc) = xp(1,ikoloc)-vp(1,ikoloc)*tstep
       xpp(2,ikoloc) = xp(2,ikoloc)-vp(2,ikoloc)*tstep
       xpp(3,ikoloc) = xp(3,ikoloc)-vp(3,ikoloc)*tstep
       !                                                !Conditions periodiques

    end if

#else
    write (6, 577) xp(1:3,iko)*1.0d8     
    write(6,'("Positions initiales du projectile (CRYST)",3(f8.4,1x))') xx0,yy0,zz0
    aux1 = sqrt(eko*ecgs*2./cm(ityp(iko)))    !Vitesse en cgs
    vp(1,iko) = vp(1,iko)+z1*aux1
    vp(2,iko) = vp(2,iko)+z2*aux1
    vp(3,iko) = vp(3,iko)+z3*aux1
575 continue
    xpp(1,iko) = xp(1,iko)-vp(1,iko)*tstep
    xpp(2,iko) = xp(2,iko)-vp(2,iko)*tstep
    xpp(3,iko) = xp(3,iko)-vp(3,iko)*tstep

    !                                                !Conditions periodiques


    ! correction de la derive par ajout d'une impulsion inverse sur les autres atomes
    if (lderive) then
       masstot=0
       do i=1,im
          if(i==iko)cycle
          masstot=masstot+cm(ityp(i))
       end do
       vpi(:)=-(vp(:,iko)*cm(ityp(iko))/masstot)
       do i=1,im
          if(i==iko)cycle
          vp(:,i)=vp(:,i)+vpi(:)
          xpp(:,i)=xpp(:,i)-vpi(:)*tstep
       end do
    end if


    if (lperiod)       call period 

#endif

    ! --- Modification de la vitesse de l'atome accelere ---



    !if (parallele)  return






    ! Choix du pas en temps initial selon le vmax
    !goto 121
    vmax2 = 0.0
    imax = 0
    vpmod2(:im) = vp(1,:im)**2+vp(2,:im)**2+vp(3,:im)**2

    do i = 1, im
       if (vpmod2(i)<=vmax2) cycle
       vmax2 = vpmod2(i)
       imax = i
    end do

    vmax = sqrt(vmax2)

#ifdef PARA
    max_loc(1)=vmax
    max_loc(2)=rang
    max_loc(3)=0.5+ityp(imax)
    call MPI_ALLREDUCE(max_loc,max_glob,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_WORLD,ierr)
    vmax = max_glob(1)
    ityp_max=int(max_glob(3))

#else
    if (rang==0) then
       write (6, *) 'Vitesse maximale sur I=', imax, vmax
    endif

#endif




    ! -> tseuil a diminuer pour eviter les derives en energies et temperature
    tseuil = depmaxts/(1.0D0*vmax)
    lts = log10(tseuil)
    expos = 1-int(lts)

#ifdef NEC
    ! NEC
    tifac1=exp10(lts+expos)
#else
    ! HP, DEC
    tifac1=10**(lts+expos)
#endif 
    ! HP, DEC

    if (tifac1<1.0) then
       write (6, *) 'sthing wrong deftimestep 1.0'
       stop
    else if (tifac1<2.0) then
       tifac2 = float(1)
    else if (tifac1<5.0) then
       tifac2 = float(2)
    else if (tifac1<=10.0) then
       tifac2 = float(5)
    else
       write (6, *) 'sthing wrong deftimestep 1.0'
       stop
    endif
    oldtstep = tstep
    tstep = tifac2

    do i = 1, expos
       tstep = tstep/float(10)
    end do

    if (tstep>tsmin) tstep=tsmin


    if (rang==0) then
       write (6, *) 'Nouveau tstep : ', tstep, '   Ancien tstep :',oldtstep
    endif                                ! rang=0


    usdh = 1/(two*tstep)

    if (l2T)then
       etstep=tstep/necyclemin
       if (etstep.gt.6d-17)then
          etstep=2d-16
          necycle=int(tstep/etstep)
          write(6,*)'chgt etstep',etstep,necycle
       end if
    end if

    ! redefinition des positions atomiques suite au changement de pas de temps
    if (dmtype==1) then
       if (tstep/=oldtstep) then
          do i = 1, im
             xpp(1,i) = xp(1,i)-vp(1,i)*tstep
             xpp(2,i) = xp(2,i)-vp(2,i)*tstep
             xpp(3,i) = xp(3,i)-vp(3,i)*tstep
          end do
       endif
    endif
121 continue
    return
  end subroutine initcasca
end module initcasca_mod
