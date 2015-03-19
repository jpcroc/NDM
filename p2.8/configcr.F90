!********************************************************************
!             CONSTRUCTION DE LA BOITE DE SIMULATION
!********************************************************************

subroutine configcr(xpcr,ityp,lrescale,itypcr)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use posana, ONLY : imcr
  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(double)  :: xpcr(3,imm)
  integer, intent(in)  :: ityp(imm)
  logical:: lrescale    ! rescale des positions de départ sur la forme de la boite d'arrivée
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, j, k, ia, ib, ic, icell, iti, icintype, icintypemod&
       , lucin, lugin, imcell, la, lb, lc, typmax, typmin, npoin, natyp, typ
       

  integer :: itypcr(imm)

  real(double), dimension(3) :: rr
  real(double), dimension(imm,3) :: xc
  real(double) :: tirax, tiray, tiraz, x1, x2, x3, a1, a2, a3, c1, c2, c3, r2
  real(double) :: rsep2,atcr(3,3),zlcr(3),bgcr(3,3)
  character :: fnamcin*80, fnamgin*80
  !              real(double) drand


  !  if(rang==0)               write(6,*)'********** reading configuration from file********'

  ! open du fichier .cin

  !APARA

  lucin = 93
  close (lucin)
  fnamcin = fnam(1:lenfnam)//'.crcin'
  open(unit=lucin, file=fnamcin, form='unformatted', status='unknown')

  read (lucin) icintype

  if (rang==0) write (6, *) 'type de fichier .cin : ', icintype
  if (icintype>3.or.icintype<0) then
     if(rang==0)                    write (6, *) 'wrong icintype'
     stop
  endif

  icintypemod = mod(icintype,2)

  if (icintype>=2) then
     read (lucin) atcr

     if (.not.lrescale) then
        if(ANY(at.ne.atcr)) then
           write(6,*) 'at <> atcr stop '
           write(6,*)at
           write(6,*)
           write(6,*)atcr

           stop
        end if
     end if
  else
     read (lucin) zlcr                      !size of the box
     if (.not.lrescale) then
        if (any(zlcr.ne.zl)) then               
           write(6,*) 'zl <> zlcr stop '
           stop
        end if
     end if
  endif                                   !icintype=2

  read (lucin) imcr                         !number of atoms in the box
  if (im.ne.imcr) then
     if(rang==0)                    write (6, *) 'ATTENTION im <> imcr', im, imcr
    ! stop
  endif

  read (lucin) itypcr                       !types
  do i=1,min(im,imcr)
     if (ityp(i).ne.itypcr(i)) write(6,*)i,ityp(i),itypcr(i)
  end do
  if (any(ityp.ne.itypcr)) then
     if(rang==0)                    write (6, *) 'ityp <> itypcr'
!     stop
  endif

  read (lucin) xpcr
  write(6,*)'lu'
  if (icintype>=2) then
     if (lrescale) then
        call recips (atcr(1,1), atcr(1,2), atcr(1,3), bgcr(1,1), bgcr(1,2), bgcr(1,3))           
        call cryst_to_cart (imm, xpcr, bgcr, -1)    !cart vers cryst           
        where (xpcr(:,:imcr)<0.0)
           xpcr(:,:imcr) = xpcr(:,:imcr)+1.
        end where
        where (xpcr(:,:imcr)>=1.0)
           xpcr(:,:imcr) = xpcr(:,:imcr)-1.0
        end where
        call cryst_to_cart (imm, xpcr, at, 1)    !cryst vers cart

     else
        call cryst_to_cart (imm, xpcr, bg, -1)    !cart vers cryst           
!        where (xpcr(:,:imcr)<0.0)
!           xpcr(:,:imcr) = xpcr(:,:imcr)+1.
!        end where
!        where (xpcr(:,:imcr)>=1.0)
!           xpcr(:,:imcr) = xpcr(:,:imcr)-1.0
!        end where
!        call cryst_to_cart (imm, xpcr, at, 1)    !cryst vers cart

     end if
     !        write(6,*)at
     !        write(6,*)atcr
     !        write(6,*)
     !        write(6,*)bg
     !        write(6,*)bgcr

  else
     if (lrescale) then
        do i=1,imcr
           do ic=1,3
              xpcr(ic,i)=xpcr(ic,i)*zl(ic)/zlcr(ic)
           end do
        end do
     end if
  end if
  return
end subroutine configcr

