module setnoxsimple_mod
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:ldemitab,pi,rang
  USE var_pot, ONLY:lpotentiel,rue_pot !ngrid,r3cm,r3cm2,rumax,q,na,rue_pot,lpotentiel,rue_pair,ntyp,csive
  USE recips_mod,only:recips,calcvol,distmin
  USE atomconfig,only: atom_config
  USE boxconfig,only:box_config
  USE cellconfig,only:cell_config
  USE arret_ndm_mod,only:arret_ndm

#ifdef PARA
  use Tpara,only:nprocspace
#endif
  implicit none
contains

  subroutine setnoxsimple(atsn,boxsn,celsn,rum,noxr,noyr,nozr)

    type(box_config),intent(in)::boxsn
    type(cell_config)::celsn
    class(atom_config)::atsn
    real(double),intent(in)::rum
    integer,optional::noxr,noyr,nozr
    integer::izonr2,natperc,nox,noy,noz,nvois,nvperat
    real(double)::zlmin,zlm2,voluperat,rvois
    zlmin = distmin(boxsn%at(:,1),boxsn%at(:,2))
    zlm2 = distmin(boxsn%at(:,1),boxsn%at(:,3))
    zlmin = min(zlmin,zlm2)
    zlm2 = distmin(boxsn%at(:,2),boxsn%at(:,3))
    zlmin = min(zlmin,zlm2)
    zlmin=zlmin*2
    nox=0;noy=0;noz=0
    if (present(noxr))nox=noxr
    if (present(noyr))noy=noyr
    if (present(nozr))noz=nozr

       ! détermination de nox noy noz qui ne sont pas donnes dans .din
       !
    if (nox<=0.or.noy<=0.or.noz<=0) then
       nox = int(boxsn%nzl(1)/rum)
       noy = int(boxsn%nzl(2)/rum)
       noz = int(boxsn%nzl(3)/rum)
    end if
!       IF (nox.LT.3) nox=1
!       IF (noy.LT.3) noy=1
!       IF (noz.LT.3) noz=1
       celsn%celsize(1) = boxsn%zl(1)/float(nox)
       celsn%celsize(2) = boxsn%zl(2)/float(noy)
       celsn%celsize(3) = boxsn%zl(3)/float(noz)
!       write(6,*)'setnoxsimple',rum, boxsn%zl(1),nox,noy,noz
    call celsn%init(boxsn,nox,noy,noz)
    
    natperc= INT(atsn%im/celsn%noxyz)
!    write(6,*)'setnoxsimple',nox,noy,noz,natperc
    natperc=max(int(2*natperc),10)     ! MODIF Clouet
    celsn%natperc=natperc
    if (allocated(celsn%atincel))deallocate(celsn%atincel)
    allocate(celsn%atincel(celsn%natperc,celsn%noxyz))
    celsn%atincel=0
    if (atsn%ltabvois) then
       rvois=atsn%rvois
       izonr2 = int(zlmin/rvois)
       !write(*,*) 'DEBUG IN DIVID volu, im', volu, im
       voluperat=boxsn%volu/atsn%im_glob
       nvperat=4*Pi*(rvois+1.0d-8)**3/(3*voluperat)
       if (ldemitab) then
          nvois=max(Int(0.8*nvperat*atsn%im_glob),100)
       else
          nvois=max(Int(1.5*nvperat*atsn%im_glob),100)
       end if

       atsn%nvois=nvois
       if(allocated(atsn%indi))deallocate(atsn%indi)
       allocate(atsn%indi(nvois))
       if (.not.allocated(atsn%iwmax))allocate(atsn%iwmax(atsn%imm))
    end if
  end subroutine setnoxsimple

  subroutine setcellsimple(atcf,box,celsp,nox,noy,noz)
    class (atom_config)::atcf
    class (box_config)::box
    class (cell_config):: celsp
    integer,intent(in)::nox,noy,noz
    integer::natpc
    if (nox<=0.or.noy<=0.or.noz<=0) then
       write(6,*)'nox noy noz must be specified in the call of setcellsimple'
       call arret_ndm
    end if
    celsp%celsize(1) = box%zl(1)/float(nox)
    celsp%celsize(2) = box%zl(2)/float(noy)
    celsp%celsize(3) = box%zl(3)/float(noz)
    natpc=int(float(atcf%im_glob)/(nox*noy*noz))
    natpc=max(int(2*natpc),10)
    call celsp%init(box,nox,noy,noz,natpc)
  end subroutine setcellsimple

    
end module setnoxsimple_mod
