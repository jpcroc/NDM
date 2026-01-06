module setcell
  USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m, ONLY:  double
  !  USE read_val,only:nox,noy,noz
  USE arret_ndm_mod,only: arret_ndm
  USE gen_com_m, ONLY:ldemitab,nvat,pi,rang,lrctest,ltpcel,lspacendm,lperiod
  USE var_pot, ONLY:lpotentiel,rue_pot,ipotentiel !ngrid,r3cm,r3cm2,rumax,q,na,rue_pot,lpotentiel,rue_pair,ntyp,csive
  USE recips_mod,only:recips,calcvol,distmin
  USE atomconfig,only: atom_config
  USE boxconfig,only:box_config
  USE cellconfig,only:cell_config
  use cryst_to_cart_mod,only:cryst_to_cart
#ifdef PARA
  use Tpara,only:nprocspace
#endif
  implicit none
contains

  subroutine setnox(boxsn,celsn,rum,lverbose,noxr,noyr,nozr)

    class(box_config),intent(in)::boxsn
    type(cell_config)::celsn
    real(double),intent(in)::rum
    integer,optional::noxr,noyr,nozr
    integer::nox(3)
    integer::izonr,ic
    logical,intent(in),optional::lverbose
    logical::lverb=.true.
    real(double)::zlmin,zlm2,ronl(3)
    if (present(lverbose)) lverb=lverbose
    nox(:)=0
    if (present(noxr))nox(1)=noxr
    if (present(noyr))nox(2)=noyr
    if (present(nozr))nox(3)=nozr

    zlmin = distmin(boxsn%at(:,1),boxsn%at(:,2))
    zlm2 = distmin(boxsn%at(:,1),boxsn%at(:,3))
    zlmin = min(zlmin,zlm2)
    zlm2 = distmin(boxsn%at(:,2),boxsn%at(:,3))
    zlmin = min(zlmin,zlm2)
    zlmin=zlmin*2

    !    if (lpotentiel(10).eqv..true.)      rut=max(rut,2*rue_pot(10))
    !    if (lpotentiel(20).eqv..true.)      rut=max(rut,2*rue_pot(20))
    !     write(6,*)'BIP',rumax,rut,rue_pot(10)
    !  end if
    !    if (lpotentiel(11).eqv..true.) rut=max(rut,2*rue_pot(11))
    !    if (lpotentiel(12).eqv..true.) rut=max(rut,2*rue_pot(12))
    RonL(:)=rum/boxsn%nzl(:)
    !    write(6,*)'IZONR',izonr
    ! MPI
    if ((rang==0).and.(lverb)) write (6, *) 'nox,noy,noz dans .din =', nox(1),nox(2),nox(3)
    do ic=1,3
       celsn%ismall(ic)=.true.
       if (Ronl(ic).Gt.0.5) then  !small direction
          nox(ic)=1+2*int(2*rum/boxsn%nzl(ic))
#ifdef PARA
          if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
             write(6,*)'trop petite boite pour para'
             call arret_ndm
          end if
#endif

       else
          celsn%ismall(ic)=.false.
          if (nox(ic).le.0) then
             
             nox(ic) = int(boxsn%nzl(ic)/rum)
!             if ((rang==0).and.(lverb)) THEN
!                write (6,*) 'nox (',ic,") calcules a partir de ru= Int( ", boxsn%nzl(1),'/',rum,') = ', nox
                !          WRITE(6,'(2(a,g12.4),a,i0)') '  nox = Int( ', boxsn%nzl(1),'/',rum,') = ', nox
!             end if
          else
             if (rang==0)write (6,'(A,I3,A,I6)') 'nox (',ic,") dans din =", nox(ic)
          end if

!!$       IF (nox.LT.3) nox=1

          celsn%celsize(ic) = boxsn%zl(ic)/float(nox(ic))

          celsn%celsize(ic) = boxsn%zl((ic))/float(nox(ic))


       end if
    end do
    if((rang==0).and.(lverb)) THEN
       write (6,'(a)') 'nox noy noz and ghost cells from ru'
       do ic=1,3
          if( celsn%ismall(ic)) then
             WRITE(6,'(a,i3,a,i5,a,g12.4)') ' GHOST DIRECTION',ic,' nox = ', nox(ic), ', =1+2*int(2*rum/boxsn%nzl(:))'&
                  &,rum/boxsn%nzl(ic)

          else
             WRITE(6,'(a,i2,a,i3,a,g12.4,a,g12.4,a)') '  nox in direction ',ic,'=',nox(ic),' if not specified =Int( '&
                  &, boxsn%nzl(1),'/',rum,') '
             !          WRITE(6,'(2(a,g12.4),a,i0)') '  noy = Int( ', boxsn%nzl(2),'/',rum,') = ', noy
             !          WRITE(6,'(2(a,g12.4),a,i0)') '  noz = Int( ', boxsn%nzl(3),'/',rum,') = ', noz
          END IF
       end do
    end if

    call celsn%init(boxsn,nox(1),nox(2),nox(3),ltpc=ltpcel)
    if ((rang==0).and.(lverb)) write(6,'(A,3G15.7)') 'celsizes ',celsn%celsize(:)

  end subroutine setnox


  subroutine setcellconf(celscf,atcf,boxcf,rumax,lverbose)
    type(cell_config)::celscf
    class(atom_config)::atcf
    class(box_config),intent(in)::boxcf
    real(double)::rumax

    integer::natperc,izonr2,nvois,nvperat
    real(double)::rm2,zlm2,zlmin,voluperat,rvois
    logical,intent(in),optional::lverbose
    logical::lverb=.true.
    if (present(lverbose)) lverb=lverbose
    zlmin = distmin(boxcf%at(:,1),boxcf%at(:,2))
    zlm2 = distmin(boxcf%at(:,1),boxcf%at(:,3))
    zlmin = min(zlmin,zlm2)
    zlm2 = distmin(boxcf%at(:,2),boxcf%at(:,3))
    zlmin = min(zlmin,zlm2)
    zlmin=zlmin*2


    !    IF (natperc.LE.0) THEN        ! MODIF Clouet
    !    write(6,*)'TTTTTTTTTTTTTTTTTTTUUUUUUUUUUUUUUUUUUUUUUUUUUUTTTTTTTTTTTTTTT'
    !    write(6,*)celscf%noxyz
    call setnatperc(celscf,atcf,boxcf,natperc)
    !    write(6,*)'natpercN',natperc,celscf%noxyz
    !    natperc= INT(atcf%im_glob/celscf%noxyz)
    !    write(6,*)'natperc0',natperc
    nvat=3*natperc
    natperc=max(int(2*natperc),20)     ! MODIF Clouet
    !    ELSE                          ! MODIF Clouet
    !       nvat=10*natperc       ! MODIF Clouet
    !    END IF                        ! MODIF Clouet


    if ((rang==0).and.(lverb)) &
         write(6,*) 'natperc im/noxyz', natperc, atcf%im_glob/celscf%noxyz
    celscf%natperc=natperc
    if (allocated(celscf%atincel))deallocate(celscf%atincel)
    allocate(celscf%atincel(celscf%natperc,celscf%noxyz))
    celscf%atincel=0
    if ((rang==0).and.(lverb))  write(6,*)'ltabvois',atcf%ltabvois

    if (atcf%ltabvois) then
       rvois=atcf%rvois
       if (rumax>rvois) then
          write (6, *) rang,' rvois trop petit rvois rumax ', rvois, rumax
          call arret_ndm
       endif
       !crc        rm2=max(rumax,2*rvois)
       rm2=max(rumax,rvois)
       izonr2 = int(zlmin/rm2)
       if (izonr2<1) then
          write (6, *) rang,'trop petite boite pour rvois !!!'
          !cosboite   call arret_ndm

          call arret_ndm
       endif
       !write(*,*) 'DEBUG IN DIVID volu, im', volu, im
       voluperat=boxcf%volu/atcf%im_glob
       nvperat=4*Pi*(rvois+1.0d-8)**3/(3*voluperat)
       if (ldemitab) then
          nvois=max(Int(0.8*nvperat*atcf%im_glob),100)
          nvat=max(Int(nvperat*1.3),10)
       else
          nvois=max(Int(1.5*nvperat*atcf%im_glob),100)
          nvat=max(Int(nvperat*1.3),10)
       end if

       if(rang==0)         write (6, *) 'Nvois= ', nvois,atcf%im_glob,nvperat,rvois,boxcf%volu,voluperat
       atcf%nvois=nvois
       if(allocated(atcf%indi))deallocate(atcf%indi)
       allocate(atcf%indi(nvois))
       !       allocate(indi2(nvois))
       if (.not.allocated(atcf%iwmax))allocate(atcf%iwmax(atcf%imm))
    end if
  end subroutine setcellconf

  subroutine setnatperc(celcf,atcf,boxcf,natperc)
    use Tpara,only:mpi_communicator,comm_space
    USE notperiod_mod,only: notperiod
    type(cell_config)::celcf
    class(atom_config)::atcf
    class(box_config),intent(in)::boxcf
    integer,intent(out)::natperc
    integer, allocatable :: natdscel(:)
    real(double), dimension(:,:), allocatable :: xpnp !
    integer::iml,i,kxyz(3),koo,midnox(3),ic
    real(double)::auxyz(3)
    midnox(:)=(celcf%nox(:)+1)/2

    iml=atcf%im

    if (celcf%noxyz==1) then
       natperc=atcf%imm
    else
       ALLOCATE(xpnp(3,iml))
       allocate (natdscel(celcf%noxyz))
       natdscel(:)=0
       call notperiod(iml,atcf%xp,xpnp,boxcf%at,boxcf%bg,lperiod)       
       call cryst_to_cart (iml, xpnp, boxcf%bg, -1) ! cart vers cryst
       do i = 1, iml
          !     if  ((it.ge.1000).and.(i.lt.20)) write(6,'(I5,3G15.7)')i, xpnp(1,i),xpnp(2,i),xpnp(3,i)
          do ic=1,3
             if  (celcf%ismall(ic)) then
                kxyz(ic)=midnox(ic)-1
             else
                auxyz(ic) = xpnp(ic,i)*celcf%nox(ic)
                kxyz(ic) = int(auxyz(ic))
                kxyz(ic) = Modulo(kxyz(ic),celcf%nox(ic))
             end if
          end do
          koo = 1+kxyz(1)+celcf%nox(1)*(kxyz(2)+celcf%nox(2)*kxyz(3))
          
!!$          aux = xpnp(1,i)*celcf%nox(1)
!!$          auy = xpnp(2,i)*celcf%nox(2)
!!$          auz = xpnp(3,i)*celcf%nox(3)
!!$          kx = int(aux)
!!$          ky = int(auy)
!!$          kz = int(auz)
!!$          kx = Modulo(kx,celcf%nox(1))
!!$          ky = Modulo(ky,celcf%nox(2))
!!$          kz = Modulo(kz,celcf%nox(3))


          IF ( (koo.GT.celcf%noxyz).OR.(koo.LT.0) ) THEN
             WRITE(0,'(a,i0,a,3g20.12)') &
                  'Problem with atom ', i, ', x,y,z = ', atcf%xp(1:3,i)
             WRITE(0,'(2(a,i0))') ' koo = ', koo, ' - noxyz = ', celcf%noxyz
             STOP
          END IF
          natdscel(koo)=natdscel(koo)+1


       end do
!       write(6,*)'NATPERCA',natdscel
#ifdef PARA
       if (lspacendm) then
          call comm_space%sum(natdscel)
       end if

#endif
       natperc=maxval(natdscel)
!       write(6,*)'NATPERCB',natperc
    end if
  end subroutine setnatperc

end module setcell
