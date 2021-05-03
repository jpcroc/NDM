module initcasca_mod
  USE cryst_to_cart_mod,only: cryst_to_cart

  USE gen_com_m, ONLY:depmaxts,dmtype,ecgs,eko,iko,lderive,lperiod,&
       &oldtstep,parallele,rang,tsmin,tstep,two,usdh,vmax,xko,xx0,yko,yy0,zko,zz0,l2T,&
       lspacendm,im_glob,imm_glob
!  use temp_com,only:at,bg,im,imm
  use constrconf_mod,only:repartition

    USE T_kind_param_m, ONLY:  double
    use atomconfig,only:atom_config_e,atom_config_d
    USE cellconfig, only:cell_config,caltabtC
    use boxconfig,only:box_config,periodbox
    USE var_pot, ONLY:cm
    USE elec_cell, ONLY : necycle,etstep,necyclemin
    ! *******************************************************************
#ifdef PARA
    USE Tpara,only:nprocspace,comm_space
    use paraconfig,only:para_config,initparapuresp
    USE parautils,only:initcomp
#else
    use Tpara,only:nprocspace


#endif


  implicit none
contains
  ! *************** initialisation de la cascade **********************
  subroutine initcasca (atcf,celndm,boxndm)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------


    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    class(atom_config_d),target::atcf
    type(box_config):: boxndm
    type(cell_config),target::celndm
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: i, ic, i1,ikoloc
    real(double) :: z1, z2, z3, t1, t2, t3, znorm, aux1

    integer :: iti, expos, imax
    real(double) :: tifac1, tifac2, lts, tseuil, vmax2
    real(double), dimension(imm_glob) :: vpmod2
    real(double) :: masstot, vpi(3)

    integer :: seed_size,isl
    integer, dimension(:), allocatable :: iseedt


#ifdef PARA
    real(double), dimension(2) :: max_loc
    integer :: ityp_max
    type(atom_config_e)::atcfcasc
    type(para_config),target::Cpara
    type(cell_config),target::celcasc

#else
    type(atom_config_e),pointer::atcfcasc
    type(cell_config),pointer::celcasc
#endif    


    select type (atcf)
       type is (atom_config_e) 

    

#ifdef PARA
    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
       call atcfcasc%init(im_glob,imm_glob)
       call initparapuresp(Cpara,rang,comm_space,nprocspace)
       call initcomp(atcfcasc,celcasc,atcf,celndm,boxndm,Cpara,lperiod)
    else
       call initparapuresp(cpara,rang,comm_space,nprocspace)
       atcfcasc=atcf
       celcasc=celndm
    end if

#else
    atcfcasc=>atcf
    celcasc=>celndm
#endif    

!!$    
!!$#if PARA
!!$    ikoloc=0
!!$    do i=1,im
!!$       if (num_at_glob(i)==iko) ikoloc=i
!!$    end do
!!$
!!$#else
!!$    ikoloc=iko
!!$#endif

    !-----------------------------------------------
    ! --- Translation de l'atome IKO au centre de la boite de simulation ---
    !    if (ikoloc.gt.0) then
    if (rang==0) then
       
       write (6, *) 'initialisation de la cascade'
       write(6,*)'ATOM ',iko, '  TYPE ',atcfcasc%ityp(iko), ' energy=',eko
       if (iko>im_glob) then
          write (6, *) 'wrong input cascade iko eko ', iko, eko
          stop
       endif
!    endif                                      ! rang=0

    if (xx0.ge.0) then

       znorm=sqrt(xx0**2+yy0**2+zz0**2)
       if (znorm==0) then
          xx0=0.25 ; yy0=0.25 ; zz0=0.25
       end if

       call cryst_to_cart (atcfcasc%imm, atcfcasc%xp, boxndm%bg, -1)    !cart vers cryst
       if (.not.atcf%lax) then
          write(6,*) 'no ax and casca stop'
          stop
       end if
       call cryst_to_cart (atcfcasc%imm, atcfcasc%ax, boxndm%bg, -1)    !cart vers cryst
       !debug write(6,*)xp(1,ikoloc),xx0
          t1 = atcfcasc%xp(1,iko)-xx0
          t2 = atcfcasc%xp(2,iko)-yy0
          t3 = atcfcasc%xp(3,iko)-zz0

       !      write(6,*)'t1 t2 t3 zl', t1,t2,t3,zl(1),zl(2),zl(3)
       atcfcasc%xpp(1,:atcfcasc%im) = atcfcasc%xpp(1,:atcfcasc%im)-t1
       atcfcasc%xpp(2,:atcfcasc%im) = atcfcasc%xpp(2,:atcfcasc%im)-t2
       atcfcasc%xpp(3,:atcfcasc%im) = atcfcasc%xpp(3,:atcfcasc%im)-t3
       atcfcasc%xp(1,:atcfcasc%im) = atcfcasc%xp(1,:atcfcasc%im)-t1
       atcfcasc%xp(2,:atcfcasc%im) = atcfcasc%xp(2,:atcfcasc%im)-t2
       atcfcasc%xp(3,:atcfcasc%im) = atcfcasc%xp(3,:atcfcasc%im)-t3
       atcfcasc%ax(1,:atcfcasc%im) = atcfcasc%ax(1,:atcfcasc%im)-t1
       atcfcasc%ax(2,:atcfcasc%im) = atcfcasc%ax(2,:atcfcasc%im)-t2
       atcfcasc%ax(3,:atcfcasc%im) = atcfcasc%ax(3,:atcfcasc%im)-t3

       call cryst_to_cart (atcfcasc%imm, atcfcasc%xp, boxndm%at, 1)     !cryst vers cart
       call cryst_to_cart (atcfcasc%imm, atcfcasc%ax, boxndm%at, 1)     !cryst vers cart
       if (lperiod)       call periodbox  (boxndm,atcfcasc)

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

    if (rang==0) then
       write (6, 576) z1, z2, z3
576    format('Direction du projectile ',3(f8.4,1x))
    endif

577 format('Positions initiales du projectile (A) ',3(f8.4,1x))
    write (6, 577) atcfcasc%xp(1:3,iko)*1.0d8     
    write(6,'("Positions initiales du projectile (CRYST)",3(f8.4,1x))') xx0,yy0,zz0
    aux1 = sqrt(eko*ecgs*2./cm(atcfcasc%ityp(iko)))    !Vitesse en cgs
    atcfcasc%vp(1,iko) = atcfcasc%vp(1,iko)+z1*aux1
    atcfcasc%vp(2,iko) = atcfcasc%vp(2,iko)+z2*aux1
    atcfcasc%vp(3,iko) = atcfcasc%vp(3,iko)+z3*aux1
575 continue
    atcfcasc%xpp(1,iko) = atcfcasc%xp(1,iko)-atcfcasc%vp(1,iko)*tstep
    atcfcasc%xpp(2,iko) = atcfcasc%xp(2,iko)-atcfcasc%vp(2,iko)*tstep
    atcfcasc%xpp(3,iko) = atcfcasc%xp(3,iko)-atcfcasc%vp(3,iko)*tstep

    !                                                !Conditions periodiques


    ! correction de la derive par ajout d'une impulsion inverse sur les autres atomes
    if (lderive) then
       masstot=0
       do i=1,atcfcasc%im
          if(i==iko)cycle
          masstot=masstot+cm(atcfcasc%ityp(i))
       end do
       vpi(:)=-(atcfcasc%vp(:,iko)*cm(atcfcasc%ityp(iko))/masstot)
       do i=1,atcfcasc%im
          if(i==iko)cycle
          atcfcasc%vp(:,i)=atcfcasc%vp(:,i)+vpi(:)
          atcfcasc%xpp(:,i)=atcfcasc%xpp(:,i)-vpi(:)*tstep
       end do
    end if


    if (lperiod)       call periodbox  (boxndm,atcfcasc)



    
    ! --- Modification de la vitesse de l'atome accelere ---



    !if (parallele)  return






    ! Choix du pas en temps initial selon le vmax
    !goto 121
    vmax2 = 0.0
    imax = 0
    vpmod2(:atcfcasc%im) = atcfcasc%vp(1,:atcfcasc%im)**2+atcfcasc%vp(2,:atcfcasc%im)**2+atcfcasc%vp(3,:atcfcasc%im)**2

    do i = 1, atcfcasc%im
       if (vpmod2(i)<=vmax2) cycle
       vmax2 = vpmod2(i)
       imax = i
    end do

    vmax = sqrt(vmax2)

           if (rang==0) then
              write (6, *) 'Vitesse maximale sur I=', imax, vmax
           endif
           




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
          do i = 1, atcfcasc%im
             atcfcasc%xpp(1,i) = atcfcasc%xp(1,i)-atcfcasc%vp(1,i)*tstep
             atcfcasc%xpp(2,i) = atcfcasc%xp(2,i)-atcfcasc%vp(2,i)*tstep
             atcfcasc%xpp(3,i) = atcfcasc%xp(3,i)-atcfcasc%vp(3,i)*tstep
          end do
       endif
    endif
121 continue

 end if
#ifdef PARA
 call atcfcasc%send2all(0,Cpara%mpi_image)
 call repartition(atcfcasc,atcf,boxndm,celndm)
#endif
 call caltabtC(celndm,atcf,lperiod,boxndm)
 
! renvoi vers les autres procs
end select

 return
  end subroutine initcasca
end module initcasca_mod
