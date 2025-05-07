module setcell
  USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m, ONLY:  double
!  USE read_val,only:nox,noy,noz
  USE arret_ndm_mod,only: arret_ndm
  USE gen_com_m, ONLY:ldemitab,nvat,pi,rang,lrctest,ltpcel,lspacendm
  USE var_pot, ONLY:lpotentiel,rue_pot,ipotentiel !ngrid,r3cm,r3cm2,rumax,q,na,rue_pot,lpotentiel,rue_pair,ntyp,csive
  USE recips_mod,only:recips,calcvol,distmin
  USE atomconfig,only: atom_config
  USE boxconfig,only:box_config
  USE cellconfig,only:cell_config
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
    integer::nox(3),ic !noy,noz,ic
    integer::izonr
    logical,intent(in),optional::lverbose
    logical::lverb=.true.
    real(double)::zlmin,zlm2,ronz(3)
    if (present(lverbose)) lverb=lverbose
    nox(3)=0!;noy=0;noz=0
    if (present(noxr))nox(1)=noxr
    if (present(noyr))nox(2)=noyr
    if (present(nozr))nox(3)=nozr
    
    zlmin = distmin(boxsn%at(:,1),boxsn%at(:,2))
    zlm2 = distmin(boxsn%at(:,1),boxsn%at(:,3))
    zlmin = min(zlmin,zlm2)
    zlm2 = distmin(boxsn%at(:,2),boxsn%at(:,3))
    zlmin = min(zlmin,zlm2)
    zlmin=zlmin*2

    izonr = int(zlmin/rum)
    ronz(:)=rum/boxsn%nzl(:)

    if ((rang==0).and.(lverb)) write (6, *) 'nox,noy,noz dans .din =', nox(:)
    do ic=1,3
       if (ronz(ic).gt.0.5)then !small box along this direction
          celsn%ismall(ic)=.true.
          celsn%celsize(ic) = boxsn%zl(1)
!          boxsn%ismall=.true.
!          celsn%ngx(ic)=
          nox(ic) =1+2*int(ronz(ic))+1
          celsn%celsize(ic) = boxsn%zl(ic)
       else
          if (nox(ic).le.0) then 
             nox(ic) = max(1,int(boxsn%nzl(ic)/rum))
             celsn%celsize(ic) = boxsn%zl(ic)/float(nox(1))
          end if
       end if
    end do
    
    if((rang==0).and.(lverb)) THEN
       write (6,'(a)') 'nox noy noz and ghost cells from ru'
       do ic=1,3
          if( celsn%ismall(ic)) then
             WRITE(6,'(a,i3,a,i5,a,g12.4)') ' GHOST DIRECTION',ic,' nox = ', nox(ic), ', =1+2*rum/boxsn%nzl(:))',rum/boxsn%nzl(:)
          else
             WRITE(6,'(a,i3,a,g12.4,a,g12.4,a)') '  nox = ',nox(ic),' if not specified =Int( ', boxsn%nzl(1),'/',rum,') '
             !          WRITE(6,'(2(a,g12.4),a,i0)') '  noy = Int( ', boxsn%nzl(2),'/',rum,') = ', noy
             !          WRITE(6,'(2(a,g12.4),a,i0)') '  noz = Int( ', boxsn%nzl(3),'/',rum,') = ', noz
          END IF
       end do
    end if


    call celsn%init(boxsn,nox(1),nox(2),nox(3),ltpc=ltpcel)
    if ((rang==0).and.(lverb)) write(6,'(A,3G15.7)') 'celsizes ',celsn%celsize(:)
    ! nox noy et noz sont determines

    !    celsn%noxyz = nox*noy*noz

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
    natperc= INT(atcf%im_glob/celscf%noxyz)
    nvat=3*natperc
    natperc=max(int(3*natperc),20)     ! MODIF Clouet
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
end module setcell
