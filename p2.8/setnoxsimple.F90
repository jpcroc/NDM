module setnoxsimple_mod
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:ldemitab,pi,rang
  USE var_pot, ONLY:lpotentiel,rue_pot !ngrid,r3cm,r3cm2,rumax,q,na,rue_pot,lpotentiel,rue_pair,ntyp,csive
  USE recips_mod,only:recips,calcvol,distmin
  USE atomconfig,only: atom_config
  USE boxconfig,only:box_config
  USE cellconfig,only:cell_config
#ifdef PARA
  use Tpara,only:nprocspace
#endif
 USE read_val,only:rvois
  implicit none
contains

  subroutine setnoxsimple(atsn,boxsn,celsn,rum)

    type(box_config),intent(in)::boxsn
    type(cell_config)::celsn
    class(atom_config)::atsn
    real(double),intent(in)::rum
    integer::izonr,ic,izonr2,natperc,nox,noy,noz,nvois,nvperat
    real(double)::zlmin,zlm2,voluperat
    
    zlmin = distmin(boxsn%at(1,1),boxsn%at(1,2))
    zlm2 = distmin(boxsn%at(1,1),boxsn%at(1,3))
    zlmin = min(zlmin,zlm2)
    zlm2 = distmin(boxsn%at(1,2),boxsn%at(1,3))
    zlmin = min(zlmin,zlm2)
    zlmin=zlmin*2

       ! détermination de nox noy noz qui ne sont pas donnes dans .din
       !
       nox = int(boxsn%nzl(1)/rum)
       noy = int(boxsn%nzl(2)/rum)
       noz = int(boxsn%nzl(3)/rum)
       IF (nox.LT.3) nox=1
       IF (noy.LT.3) noy=1
       IF (noz.LT.3) noz=1
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
!       allocate(indi2(nvois))
       if (.not.allocated(atsn%iwmax))allocate(atsn%iwmax(atsn%imm))
    end if
  end subroutine setnoxsimple
end module setnoxsimple_mod
