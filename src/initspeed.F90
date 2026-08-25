module initspeed_mod

   USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m, ONLY:  double
  USE Mat_utils_mod,only: MatInv
  USE tempinstT_mod,only: tempinstT
  USE arret_ndm_mod,only: arret_ndm
  USE gen_com_m, only:uwrt,lwrt,pi,debyetemp,dmtype,hbar,iseed,lcalcjq,lperiod,ltpcel,&
       &lvpread,oldtstep,one,tempdeplainit,tinit,tstep,iseed,&
       bk,lspacendm,mdcg_noise
  use neb_module, only : neb_noise_scale,mdcg_noise_scale
  USE var_pot, ONLY:ntyp,cm
#ifdef PARA
  USE Tpara,only:COMM_space,nprocs,nprocspace,myidsp
#else
    USE Tpara,only:nprocspace,myidsp
#endif

  USE atomconfig,only:atom_config,atom_config_d,atom_config_e
  use boxconfig,only:box_config,periodbox
  implicit none
  logical:: lwrts
contains
  ! *********************************************************************
  subroutine bruit_xp (bruitmd,im)
    USE T_kind_param_m, ONLY:  double

    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
!    class(atom_config)::atcf
    real(double), allocatable::bruitmd(:,:)
    integer::im
    
    integer    :: ia,seed_size
    integer, dimension(:),allocatable :: iseedt
    real(double)  :: zr1,zr2,zr3,zr4,totalbruit

    if (.not.allocated(bruitmd))allocate(bruitmd(3,im))
    call random_seed(size=seed_size)
    allocate(iseedt(seed_size))
    call system_clock (iseed)

    iseedt(:)=iseed

    call random_seed(put=iseedt)
    totalbruit=0.d0
    bruitmd(1:3,1:im)=0.d0
    do ia=1,im
       call random_number(zr1)
       call random_number(zr2)
       call random_number(zr3)
       call random_number(zr4)
       if(zr1.eq.0.d0) zr1=0.000000001d0
       if(zr2.eq.0.d0) zr2=0.000000001d0
       if(zr3.eq.0.d0) zr3=0.000000001d0
       if(zr4.eq.0.d0) zr4=0.000000001d0

       bruitmd(1,ia)=sqrt((-log(zr1)))*cos(2.0*pi*zr3)
       bruitmd(2,ia)=sqrt((-log(zr1)))*sin(2.0*pi*zr3)
       bruitmd(3,ia)=sqrt((-log(zr2)))*cos(2.0*pi*zr4)
       totalbruit=totalbruit + bruitmd(1,ia)**2 + bruitmd(2,ia)**2 + bruitmd(3,ia)**2
    end do

    if (lwrt)  write(uwrt,*) 'ISEED for MD, NORM of the noise ',iseed, neb_noise_scale, totalbruit
!    bruitmd(1:3,1:im) = bruitmd(1:3,1:im) * mdcg_noise_scale * xp(1:3,1:im) / (sqrt(totalbruit))
    bruitmd(1:3,1:im) = bruitmd(1:3,1:im) * mdcg_noise_scale  / (sqrt(totalbruit))
    bruitmd=bruitmd*1d-8
  end subroutine bruit_xp




  ! *********************************************************************
  subroutine initspeed(atcf,boxndm,latcomp,lprt,tinitr,lwrtsR)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE var_pot, ONLY:
    implicit none
    logical,optional::lwrtsR
    !-----------------------------------------------
    class(atom_config_d)::atcf
    !    type(cell_config):: celndm
    class(box_config),intent(in)::boxndm
    logical,optional,intent(in):: latcomp
    logical,optional::lprt
    logical::latc=.false.,lprint=.true.
    integer :: i, ic, ia
    integer, dimension(:), allocatable :: iseedt
    ! ym      real(double), dimension(nce) :: tempc
    !  real(double), dimension(noxyz) :: tempc
    real(double) :: vv, v0, v1, z1, z2, z3, z4
    real(double) :: tempsauv ,totmass
    real(double) :: rx, ry, rz, r2x, r2y, r2z, r2
    real(double) :: prx, pry, prz, px, py, pz, vrx, vry, vrz,ka
    real(double) :: omegax, omegay, omegaz
    real(double), dimension(3) :: vt1,  scom, pav,kinx
    real(double), dimension(3,3) :: ainer, aineri
    real(double), dimension(3,ntyp) :: vav
    real(double),allocatable::bruitmd(:,:)
    integer :: seed_size
    integer::iti,imtot
    real(double)::sd,grnd,theta,fhi ! ,decx(2)

    real(double),optional::tinitr
    real(double)::tinit0
    if(present(lwrtSR)) then
       lwrts=lwrtsR
    else
       lwrts=lwrt
    end if

    if(present(tinitr)) then
       tinit0=tinitr
    else
       tinit0=tinit
    end if

    if (present(lprt))lprint=lprt
!    if (lwrts) write(uwrt,*) 'PARA-T entree initspeed',iseed,lvpread
    if(present(latcomp))latc=latcomp
!    if (lwrts) write(uwrt,*)
!    write(uwrt,*)'ISEED initspeed',iseed
    select case (dmtype)
    case(3,30,5,11,31,32,33,21,22,2)
       atcf%vp = 0.0
       if (mdcg_noise==0) then 
          atcf%vp=0.0;
          select type (atcf)
          class is (atom_config_e)
             if (atcf%lxpp)atcf%xpp=atcf%xp
          end select
       else
          atcf%vp=0.0 
          call bruit_xp(bruitmd,atcf%im)
          atcf%xp(1:3,1:atcf%im) = atcf%xp(1:3,1:atcf%im) + bruitmd(1:3,1:atcf%im)
       end if
       goto 66
    end select
1   continue
    if (lvpread) then
       tempsauv=tempinstT(atcf)
      if (myidsp==0) write(uwrt,*)'tempsauv ',tempsauv
       select type (atcf)
       class is (atom_config_e)
                   if (atcf%lxpp)  atcf%xpp(1:3,1:atcf%im) = atcf%xp(1:3,1:atcf%im)-atcf%vp(1:3,1:atcf%im)*tstep
       end select
       !     vp(:,:im)=vp(:,:im)*tstep/oldtstep

       if (tinit0<0) then
          ! velocities are read from file and not modified
       if ((lwrts).and.(lprint)) write (uwrt,*) 'pas de chgt des vitesses= '
          !        return
       else
          ! velocities are read from file and rescaled
                 if ((lwrts).and.(lprint)) write (uwrt,*) 'scaling read velocities at TINIT = ', &
               tinit0, 'K'
          !   if (lwrts) write(uwrt,*)'tempsauv ',tempsauv
          if (tempsauv.le.1.) then
             lvpread=.false. ; goto 1
          end if
          vv = sqrt(tinit0/tempsauv)
          atcf%vp(:,:atcf%im) = atcf%vp(:,:atcf%im)*vv
          
          select type (atcf)
          class is (atom_config_e)
                   if (atcf%lxpp)              atcf%xpp(:,:atcf%im) = atcf%xp(:,:atcf%im)-atcf%vp(:,:atcf%im)*tstep
          end select
       endif

    else
!       write(6,*)'TINIT0',tinit0

       if (tinit0<=0) then
          ! velocities are not read and no starting temperature is given
             atcf%vp(1,:atcf%im) = 0.0
             atcf%vp(2,:atcf%im) = 0.0
             atcf%vp(3,:atcf%im) = 0.0
       select type (atcf)
       class is (atom_config_e)
          if (atcf%lxpp) then
             atcf%xpp(1,:atcf%im) = atcf%xp(1,:atcf%im)
             atcf%xpp(2,:atcf%im) = atcf%xp(2,:atcf%im)
             atcf%xpp(3,:atcf%im) = atcf%xp(3,:atcf%im)
          end if
       end select
       if ((lwrts).and.(lprint))  write (uwrt,*) 'ZERO VELOCITY '
       else
          !  a starting temperature is given
                 if ((lwrts).and.(lprint))  write (uwrt,*) 'random velocities at TINIT = ', tinit0, &
               'K'
                 call random_seed(size=seed_size)
!          if (lwrts)write(uwrt,*)'seed_size',seed_size
          allocate(iseedt(seed_size))

          iseedt(:)=iseed
          call    random_seed (put=iseedt)
          deallocate(iseedt)

          v0 = sqrt(2.D0*bk*tinit0)
          vt1(:)=0.0
          do i = 1, atcf%im
             !******************************************
             call random_number(z1)
             call random_number(z2)
             call random_number(z3)
             call random_number(z4)
             if(z1.eq.0.d0) z1=0.000000001d0
             if(z2.eq.0.d0) z2=0.000000001d0
             if(z3.eq.0.d0) z3=0.000000001d0
             if(z4.eq.0.d0) z4=0.000000001d0
             v1 = one/sqrt(cm(atcf%ityp(i)))
             atcf%vp(1,i) = v1*v0*sqrt((-log(z1)))*cos(2.0*pi*z2)
             atcf%vp(2,i) = v1*v0*sqrt((-log(z1)))*sin(2.0*pi*z3)
             atcf%vp(3,i) = v1*v0*sqrt((-log(z2)))*cos(2.0*pi*z4)
             theta=acos(1-2*z3)
             fhi=2*pi*z4
             ! vp(1,i) = v1*v0*sqrt((-log(z1)))*sin(theta)*cos(fhi)
             ! vp(2,i) = v1*v0*sqrt((-log(z1)))*sin(theta)*sin(fhi)
             ! vp(3,i) = v1*v0*sqrt((-log(z2)))*cos(theta)


             !             endif
          end do
          kinx(:)=0.d0
          do ic=1,3
             do i=1,atcf%im
                kinx(ic)=kinx(ic)+0.5*atcf%vp(ic,i)*atcf%vp(ic,i)*cm(atcf%ityp(i))
                !write(*,*) ic, i,kinx(ic),  cm(ityp(i)), vp(ic,i)
             end do
             ka=0.5*bk*tinit0 
!!$             imtot=atcf%im
!!$#ifdef PARA
!!$             if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
!!$                call comm_space%sum(imtot)
!!$                call comm_space%sum(kinx)
!!$             end if
!!$#endif
!!$             kinx(ic)=kinx(ic)+0.5*atcf%vp(ic,i)*atcf%vp(ic,i)*cm(atcf%ityp(i))/imtot
             !if (lwrts) write(uwrt,*)'dir ',ic,' ka Ktinit ', kinx(ic),ka
          end do
          !***************************************************************
          !       Make total momentum zero
          !         Calculate the global average momemtum & center-of-mass


          totmass = 0.d0
          scom = 0.d0
          pav  = 0.d0

          do i = 1, atcf%im
             ic = atcf%ityp(i)
             totmass = totmass+cm(ic)
             do ia = 1, 3
                scom(ia) = scom(ia)+atcf%xp(ia,i)*cm(ic)
                pav (ia) = pav(ia) +atcf%vp(ia,i)*cm(ic)
             enddo
          enddo
          imtot=atcf%im
          
#ifdef PARA
          if (.not.latc) then 
             if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
                call comm_space%sum(totmass)
                call comm_space%sum(scom)
                call comm_space%sum(pav)
                call comm_space%sum(imtot)
             end if
          end if
#endif          

     
          do ia = 1, 3
             scom(ia) = scom(ia)/totmass
             pav (ia) = pav (ia)/float(imtot)
          enddo

          if (lwrts) then
!             write (uwrt,*) 'center of mass = ',scom(1)*1d8,&
!                  scom(2)*1d8,scom(3)*1d8
!             write (uwrt,*) 'Momentum/atom  = ',pav(1),pav(2),pav(3)
          endif

          !         Shift velocities to make the total momemtum zero
          do ic = 1, ntyp
             do ia = 1, 3
                vav(ia,ic) = pav(ia)/cm(ic)
             enddo
          enddo

          do i = 1, atcf%im
             do ia = 1, 3
                atcf%vp(ia,i) = atcf%vp(ia,i)-vav(ia,atcf%ityp(i))
             enddo
          enddo

          if (lcalcjq) then
             ka=0.5*bk*tinit0
             kinx=0.
             do ic=1,3
                do i=1,atcf%im
                   kinx(ic)=kinx(ic)+0.5*atcf%vp(ic,i)*atcf%vp(ic,i)*cm(atcf%ityp(i))
                end do
             imtot=atcf%im
#ifdef PARA
             if (.not.latc) then 
                if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
                   call comm_space%sum(imtot)
                   call comm_space%sum(kinx)
                end if
             end if
#endif
             kinx(ic)=kinx(ic)+0.5*atcf%vp(ic,i)*atcf%vp(ic,i)*cm(atcf%ityp(i))/imtot
                if (lwrts)write(uwrt,*)'dir ',ic,' ka Ktinit ', kinx(ic),ka
                do i=1,atcf%im
                   atcf%vp(ic,i)= atcf%vp(ic,i)*dsqrt(ka/kinx(ic))
                end do
             end do
          else
             ainer = 0.d0
             prx = 0.d0
             pry = 0.d0
             prz = 0.d0

             do i = 1, atcf%im
                ic = atcf%ityp(i)
                rx = atcf%xp(1,i)-scom(1)
                ry = atcf%xp(2,i)-scom(2)
                rz = atcf%xp(3,i)-scom(3)
                r2x = rx*rx
                r2y = ry*ry
                r2z = rz*rz
                r2  = r2x+r2y+r2z
                ainer(1,1) = ainer(1,1)+cm(ic)*(r2-r2x)
                ainer(2,2) = ainer(2,2)+cm(ic)*(r2-r2y)
                ainer(3,3) = ainer(3,3)+cm(ic)*(r2-r2z)
                ainer(2,3) = ainer(2,3)-cm(ic)*ry*rz
                ainer(3,1) = ainer(3,1)-cm(ic)*rz*rx
                ainer(1,2) = ainer(1,2)-cm(ic)*rx*ry
                px  = cm(ic)*atcf%vp(1,i)
                py  = cm(ic)*atcf%vp(2,i)
                pz  = cm(ic)*atcf%vp(3,i)
                prx = prx+ry*pz-rz*py
                pry = pry+rz*px-rx*pz
                prz = prz+rx*py-ry*px
             enddo

             ainer(3,2) = ainer(2,3)
             ainer(1,3) = ainer(3,1)
             ainer(2,1) = ainer(1,2)

#ifdef PARA
             if (.not.latc) then 
                if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
                   call comm_space%sum(ainer)
                   call comm_space%sum(prx)
                   call comm_space%sum(pry)
                   call comm_space%sum(prz)  
                end if
             end if
#endif          

!             if (lwrts) then
!                write(uwrt,997) (ainer(1,ib),ib=1,3),prx
!                write(uwrt,997) (ainer(2,ib),ib=1,3),pry
!                write(uwrt,997) (ainer(3,ib),ib=1,3),prz
!             endif
997          format('Inertia/anglm = ',3e12.4,5x,e12.4)

             !         calculate  angular velocity
             call matinv(ainer,aineri)
             omegax = aineri(1,1)*prx+aineri(1,2)*pry+aineri(1,3)*prz
             omegay = aineri(2,1)*prx+aineri(2,2)*pry+aineri(2,3)*prz
             omegaz = aineri(3,1)*prx+aineri(3,2)*pry+aineri(3,3)*prz

             !         shift velocities to make the angular momentum zero
             do i = 1, atcf%im
                rx = atcf%xp(1,i)-scom(1)
                ry = atcf%xp(2,i)-scom(2)
                rz = atcf%xp(3,i)-scom(3)
                vrx = omegay*rz-omegaz*ry
                vry = omegaz*rx-omegax*rz
                vrz = omegax*ry-omegay*rx
                atcf%vp(1,i) = atcf%vp(1,i)-vrx
                atcf%vp(2,i) = atcf%vp(2,i)-vry
                atcf%vp(3,i) = atcf%vp(3,i)-vrz
             enddo

             !      Old velocities translation
             !            atcf%vp(1,:im) = vp(1,:im)-vt1(1)/float(im)
             !            vp(2,:im) = vp(2,:im)-vt1(2)/float(im)
             !            vp(3,:im) = vp(3,:im)-vt1(3)/float(im)

             !      Initialisation of the previous position for Verlet

          end if
       end if
       tempsauv=tempinstT(atcf)
       if (tempsauv.ne.0)       then
          vv = sqrt(tinit0/tempsauv)
          atcf%vp(:,:atcf%im) = atcf%vp(:,:atcf%im)*vv
          select type (atcf)
          class is (atom_config_e)
                   if (atcf%lxpp)           atcf%xpp(:,:atcf%im) = atcf%xp(:,:atcf%im)-atcf%vp(:,:atcf%im)*tstep
          end select

       end if
    endif
    !     write(uwrt,*)'sortie initspeed'
    select type (atcf)
    class is (atom_config_e)
       if (atcf%lxpp) then 
          atcf%xpp(:,:atcf%im) = atcf%xp(:,:atcf%im)-atcf%vp(:,:atcf%im)*tstep
       end if
    end select
    
    tempsauv=tempinstT(atcf)
    if (lwrts) write(uwrt,*)'temperature fin initspeed ',tempsauv

    if (mdcg_noise==1) then
       call bruit_xp (bruitmd,atcf%im)
       do i=1,atcf%im
          ! Variables = cartesian coordinates (in cm)
          atcf%xp(1:3,i)= atcf%xp(1:3,i)+bruitmd(1:3,i)
       end do
    end if
    if (tempdeplainit.gt.0)then

       if ((lwrts).and.(lprint)) then
          write(uwrt,*)
          write(uwrt,*)'-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
          write(uwrt,*)'depla init Tempdeplainit',tempdeplainit,'debyetemp= ',debyetemp

          do iti=1,ntyp
             sd=sqrt((3*tempdeplainit*hbar**2)/(bk*cm(iti)*debyetemp**2))
             write(uwrt,*)'sd2 de iti',sd*sd,iti
          end do
       end if
       !       decx=0
       do i=1,atcf%im
          sd= sqrt((3*tempdeplainit*hbar**2)/(bk*cm(atcf%ityp(i))*debyetemp**2))

          do ic=1,3
             call gaussianrand(grnd)
             !                        write(uwrt,*)grnd
             atcf%xp(ic,i)=atcf%xp(ic,i)+sd*grnd
             select type (atcf)
             class is (atom_config_e)
                if (atcf%lxpp) then
                   atcf%xpp(ic,i)=atcf%xpp(ic,i)+sd*grnd
                end if
             end select!             decx(ityp(i))=decx(ityp(i))+(sd*grnd)**2
          end do
       end do
       call periodbox (boxndm,atcf)
    end if
66  continue


!! Les envois sont à faire dans l'éventuelle routine appellante
!!$#ifdef PARA
!!$       if ((latc).and.(nprocspace.gt.1)) then ! les procs masters myidsp=0 ont toutes les positions., Il faut passer aux autres procs les nouvelles atcf
!!$          call atcf%send2all(0,comm_space)
!!$       end if
!!$#endif

    return




  end subroutine initspeed

  subroutine gaussianrand(gr)
    USE T_kind_param_m
    USE gen_com_m, only:uwrt,lwrt,pi
    implicit none
    real(double),intent(out)::gr

    real(double):: z1,z2

1   continue
    call random_number(z1)  
    call random_number(z2)  
    gr=sqrt(-2*log(z1))*cos(2*pi*z2)


  end subroutine  gaussianrand


  subroutine init_speed_1at(vp, temp,xp,iti,xpp)
    real(double),intent(in)::temp,xp(3)
    integer,intent(in)::iti
    real(double),intent(out)::vp(3)
    real(double),intent(out), optional ::xpp(3)
    
    real(double) :: v0, v1, z1, z2, z3, z4

    v0 = sqrt(2.D0*bk*temp)
    call random_number(z1)
    call random_number(z2)
    call random_number(z3)
    call random_number(z4)
    if(z1.eq.0.d0) z1=0.000000001d0
    if(z2.eq.0.d0) z2=0.000000001d0
    if(z3.eq.0.d0) z3=0.000000001d0
    if(z4.eq.0.d0) z4=0.000000001d0

    v1 = one/sqrt(cm(iti))
    vp(1) = v1*v0*sqrt((-log(z1)))*cos(2.0*pi*z3)
    vp(2) = v1*v0*sqrt((-log(z1)))*sin(2.0*pi*z3)
    vp(3) = v1*v0*sqrt((-log(z2)))*cos(2.0*pi*z4)
    if (present(xpp) )xpp(:) =xp(:)-vp(:)*tstep
  end subroutine init_speed_1at

    
end module initspeed_mod
