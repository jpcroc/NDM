module setcell
#ifndef ML
#endif
  USE T_kind_param_m, ONLY:  double
  USE read_val,only:nox,noy,noz,rvois
  USE arret_ndm_mod,only: arret_ndm
  USE gen_com_m, ONLY:lconstrtot,ldemitab,lconstrtot,nvat,pi,rang,lrctest,ltpcel,lspacendm
    USE var_pot, ONLY:lpotentiel,rue_pot !ngrid,r3cm,r3cm2,rumax,q,na,rue_pot,lpotentiel,rue_pair,ntyp,csive
  USE recips_mod,only:recips,calcvol,distmin
  USE atomconfig,only: atom_config
  USE boxconfig,only:box_config
  USE cellconfig,only:cell_config
#ifdef PARA
  use Tpara,only:nprocspace
#endif
  implicit none
contains

  subroutine setnox(boxsn,celsn,rum,lverbose)

    type(box_config),intent(in)::boxsn
    type(cell_config)::celsn
    real(double),intent(in)::rum
    integer::izonr,ic
    logical,intent(in),optional::lverbose
    logical::lverb=.true.
    real(double)::rut,zlmin,zlm2
    if (present(lverbose)) lverb=lverbose
    
    zlmin = distmin(boxsn%at(1,1),boxsn%at(1,2))
    zlm2 = distmin(boxsn%at(1,1),boxsn%at(1,3))
    zlmin = min(zlmin,zlm2)
    zlm2 = distmin(boxsn%at(1,2),boxsn%at(1,3))
    zlmin = min(zlmin,zlm2)
    zlmin=zlmin*2

    rut=rum
    if (lpotentiel(10).eqv..true.)      rut=max(rut,2*rue_pot(10))
    if (lpotentiel(20).eqv..true.)      rut=max(rut,2*rue_pot(20))
    !     write(6,*)'BIP',rumax,rut,rue_pot(10)
    !  end if
    if (lpotentiel(11).eqv..true.) rut=max(rut,2*rue_pot(11))
    if (lpotentiel(12).eqv..true.) rut=max(rut,2*rue_pot(12))
     izonr = int(zlmin/rut)
    ! MPI
    if ((rang==0).and.(lverb)) write (6, *) 'izonr,zlmin,rut', izonr, zlmin*1d8, rut*1d8
    if (izonr<2) then
       write (6, *) 'trop petite boite !!!'
       !cosboite  stop
       if (lrctest) then
          write (6, *) 'STOP ; supprimer avec lrctest=.false. dans din'
          stop
       endif
    end if
     if ((rang==0).and.(lverb)) write (6, *) 'nox,noy,noz dans .din =', nox, noy, noz
    
    if (nox<=0.or.noy<=0.or.noz<=0) then
       ! détermination de nox noy noz qui ne sont pas donnes dans .din
       !
        if ((rang==0).and.(lverb))write (6, *) 'calcul de nox noy noz !!!'
       ! ==== MODIF CLOUET 2 ====================
       if (izonr<3) then
#ifdef PARA
if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
             write(6,*)'trop petite boite pour para'
             call arret_ndm
          end if
#endif
             
           if ((rang==0).and.(lverb)) then
             WRITE(6,'(a)') "Boite trop petite: le nombre de cellules est fixe a son minimum"
          endif
       endif
       nox = int(boxsn%nzl(1)/rum)
       noy = int(boxsn%nzl(2)/rum)
       noz = int(boxsn%nzl(3)/rum)
         if ((rang==0).and.(lverb)) THEN
          write (6,'(a)') 'nox noy noz calcules a partir de ru'
          WRITE(6,'(2(a,g12.4),a,i0)') '  nox = Int( ', boxsn%nzl(1),'/',rum,') = ', nox
          WRITE(6,'(2(a,g12.4),a,i0)') '  noy = Int( ', boxsn%nzl(2),'/',rum,') = ', noy
          WRITE(6,'(2(a,g12.4),a,i0)') '  noz = Int( ', boxsn%nzl(3),'/',rum,') = ', noz
       END IF
       IF (nox.LT.3) nox=1
       IF (noy.LT.3) noy=1
       IF (noz.LT.3) noz=1
       if ((rang==0).and.(lverb)) write (6,'(a,3(i0,1x))') 'nox noy noz apres correction = '&
            , nox, noy, noz

!       IF ( (nox.LE.3).AND.(noy.LE.3).AND.(noz.LE.3) ) THEN
!          nox=1 ; noy=1 ; noz=1
!          !             ltabvois=.TRUE.
          !             lconstrtot=.TRUE.
!          if (rang==0) write(6,*)'!!!!!!!!!!Envisager ltabvois = true !!!!!!!!!!!!!!'
!       END IF
       ! ==== FIN MODIF CLOUET 2 ================
!       celsn%nox=nox
!       celsn%noy=noy
!       celsn%noz=noz
       celsn%celsize(1) = boxsn%zl(1)/float(nox)
       celsn%celsize(2) = boxsn%zl(2)/float(noy)
       celsn%celsize(3) = boxsn%zl(3)/float(noz)

    else

       ! *** nox noy noz sont donnes dans.din ***

       if (nox==2.or.noy==2.or.noz==2) then
          write (6, *) rang,'wrong noxyz stop'
          call arret_ndm
       endif
!       celsn%nox=nox;celsn%noy=noy;celsn%noz=noz

       celsn%celsize(1) = boxsn%zl(1)/float(nox)
       celsn%celsize(2) = boxsn%zl(2)/float(noy)
       celsn%celsize(3) = boxsn%zl(3)/float(noz)


    endif
    call celsn%init(boxsn,nox,noy,noz,ltpc=ltpcel)
     if ((rang==0).and.(lverb)) write(6,'(A,3G15.7)') 'celsizes ',celsn%celsize(:)
    ! nox noy et noz sont determines

!    celsn%noxyz = nox*noy*noz
    
  end subroutine setnox


  subroutine setcellconf(celscf,atcf,boxcf,im_glob,rumax,lverbose)
    type(cell_config)::celscf
    class(atom_config)::atcf
    type(box_config),intent(in)::boxcf
    integer,intent(in)::im_glob
    real(double)::rumax

    integer::natperc,izonr2,nvois,nvperat
    real(double)::rm2,zlm2,zlmin,voluperat
    logical,intent(in),optional::lverbose
    logical::lverb=.true.
    if (present(lverbose)) lverb=lverbose
    zlmin = distmin(boxcf%at(1,1),boxcf%at(1,2))
    zlm2 = distmin(boxcf%at(1,1),boxcf%at(1,3))
    zlmin = min(zlmin,zlm2)
    zlm2 = distmin(boxcf%at(1,2),boxcf%at(1,3))
    zlmin = min(zlmin,zlm2)
    zlmin=zlmin*2

    
!    IF (natperc.LE.0) THEN        ! MODIF Clouet
!    write(6,*)'TTTTTTTTTTTTTTTTTTTUUUUUUUUUUUUUUUUUUUUUUUUUUUTTTTTTTTTTTTTTT'
!    write(6,*)celscf%noxyz
!    write(6,*)im_glob
    natperc= INT(im_glob/celscf%noxyz)
       nvat=3*natperc
       natperc=max(int(2*natperc),10)     ! MODIF Clouet
!    ELSE                          ! MODIF Clouet
!       nvat=10*natperc       ! MODIF Clouet
!    END IF                        ! MODIF Clouet


    if ((rang==0).and.(lverb)) &
         write(6,*) 'natperc im/noxyz', natperc, im_glob/celscf%noxyz
    celscf%natperc=natperc
    if (allocated(celscf%atincel))deallocate(celscf%atincel)
    allocate(celscf%atincel(celscf%natperc,celscf%noxyz))
    celscf%atincel=0
  if ((rang==0).and.(lverb))  write(6,*)'ltabvois,lconstrtot',atcf%ltabvois,lconstrtot

    if (atcf%ltabvois) then
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
!       if(.not.lconstrtot)rumax=rvois
       !write(*,*) 'DEBUG IN DIVID volu, im', volu, im
       voluperat=boxcf%volu/im_glob
       nvperat=4*Pi*(rvois+1.0d-8)**3/(3*voluperat)
       if (ldemitab) then
          nvois=max(Int(0.8*nvperat*im_glob),100)
          nvat=max(Int(nvperat*1.3),10)
       else
          nvois=max(Int(1.5*nvperat*im_glob),100)
          nvat=max(Int(nvperat*1.3),10)
       end if

       if(rang==0)         write (6, *) 'Nvois= ', nvois,im_glob,nvperat,rvois,boxcf%volu,voluperat
       atcf%nvois=nvois
       if(allocated(atcf%indi))deallocate(atcf%indi)
       allocate(atcf%indi(nvois))
!       allocate(indi2(nvois))
       if (.not.allocated(atcf%iwmax))allocate(atcf%iwmax(atcf%imm))
    end if
  end subroutine setcellconf
end module setcell
