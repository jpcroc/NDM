! ***********************************************************
!    Sous-programme initialisation de la FFT
!                Version du 10/12/2001
! ***********************************************************
subroutine fftsetup (n1,n2,n3,ntable,table)
  implicit none
  integer n1,n2,n3,ntable
  real*8 table(ntable,3)
  call cffti (n1,table(1,1))
  call cffti (n2,table(1,2))
  call cffti (n3,table(1,3))
  return
end subroutine fftsetup


subroutine cffti (n,wsave)
  implicit none
  integer n,iw1,iw2
  real*8 wsave(*)

  if (n .eq. 1) return
  iw1 = n+n+1
  iw2 = iw1+n+n
  call cffti1 (n,wsave(iw1),wsave(iw2))
  return
end subroutine cffti


subroutine cffti1 (n,wa,ifac)
  implicit none
  integer i,j,ii,n,ifac(*)
  integer ib,ido,idot,ip,ipm
  integer i1,k1,l1,l2,ld
  integer nl,nf,nq,nr
  integer ntry,ntryh(4)
  real*8 arg,argh,argld
  real*8 fi,wa(*)
  real*8 pi
  data ntryh / 3, 4, 2, 5 /

  pi=3.141592654D0
  !
  nl = n
  nf = 0
  j = 0
10 continue
  j = j + 1
  if (j .le. 4) then
     ntry = ntryh(j)
  else
     ntry = ntry + 2
  end if
20 continue
  nq = nl / ntry
  nr = nl - ntry*nq
  if (nr .ne. 0)  goto 10
  nf = nf + 1
  ifac(nf+2) = ntry
  nl = nq
  if (ntry .eq. 2) then
     if (nf .ne. 1) then
        do i = 2, nf
           ib = nf - i + 2
           ifac(ib+2) = ifac(ib+1)
        end do
        ifac(3) = 2
     end if
  end if
  if (nl .ne. 1)  goto 20
  ifac(1) = n
  ifac(2) = nf
  argh = 2.0d0 * pi / dble(n)
  i = 2
  l1 = 1
  do k1 = 1, nf
     ip = ifac(k1+2)
     ld = 0
     l2 = l1 * ip
     ido = n / l2
     idot = ido + ido + 2
     ipm = ip - 1
     do j = 1, ipm
        i1 = i
        wa(i-1) = 1.0d0
        wa(i) = 0.0d0
        ld = ld + l1
        fi = 0.0d0
        argld = dble(ld) * argh
        do ii = 4, idot, 2
           i = i + 2
           fi = fi + 1.0d0
           arg = fi * argld
           wa(i-1) = cos(arg)
           wa(i) = sin(arg)
        end do
        if (ip .gt. 5) then
           wa(i1-1) = wa(i-1)
           wa(i1) = wa(i)
        end if
     end do
     l1 = l2
  end do
  return
end subroutine cffti1

!     ##########################################################
!     ##                                                      ##
!     ##  subroutine fftfront  --  3-D FFT forward transform  ##
!     ##                                                      ##
!     ##########################################################
subroutine fftfront (n1,n2,n3,kpmex,kpmey,kpmez,kpme,ntable,table,w)
  implicit none
  integer i,j,k
  integer n1,n2,n3
  integer ntable
  integer kpme,kpmex,kpmey,kpmez
  real*8 table(ntable,3)
  double complex work(kpme)
  double complex w(kpmex,kpmey,kpmez)
  !
  !     forward transform along X, Y and then Z
  !
  do k = 1, n3
     do j = 1, n2
        do i = 1, n1
           work(i) = w(i,j,k)
        end do
        call cfftf (n1,work,table(1,1))
        do i = 1, n1
           w(i,j,k) = work(i)
        end do
     end do
  end do
  do k = 1, n3
     do i = 1, n1
        do j = 1, n2
           work(j) = w(i,j,k)
        end do
        call cfftf (n2,work,table(1,2))
        do j = 1, n2
           w(i,j,k) = work(j)
        end do
     end do
  end do
  do i = 1, n1
     do j = 1, n2
        do k = 1, n3
           work(k) = w(i,j,k)
        end do
        call cfftf (n3,work,table(1,3))
        do k = 1, n3
           w(i,j,k) = work(k)
        end do
     end do
  end do
  return
end subroutine fftfront

! ***********************************************************
!    Sous-programme de calcul de la FFT inverse
!                Version du 10/12/2001
! ***********************************************************
subroutine fftback (n1,n2,n3, &
     kpmex,kpmey,kpmez,kpme,ntable,table,w)
  implicit none
  integer i,j,k
  integer n1,n2,n3
  integer ntable
  integer kpme,kpmex,kpmey,kpmez
  real*8 table(ntable,3)
  double complex work(kpme)
  double complex w(kpmex,kpmey,kpmez)
  !
  !     backward transform along X, Y and then Z
  !
  do k = 1, n3
     do j = 1, n2
        do i = 1, n1
           work(i) = w(i,j,k)
        end do
        call cfftb (n1,work,table(1,1))
        do i = 1, n1
           w(i,j,k) = work(i)
        end do
     end do
  end do
  do k = 1, n3
     do i = 1, n1
        do j = 1, n2
           work(j) = w(i,j,k)
        end do
        call cfftb (n2,work,table(1,2))
        do j = 1, n2
           w(i,j,k) = work(j)
        end do
     end do
  end do
  do i = 1, n1
     do j = 1, n2
        do k = 1, n3
           work(k) = w(i,j,k)
        end do
        call cfftb (n3,work,table(1,3))
        do k = 1, n3
           w(i,j,k) = work(k)
        end do
     end do
  end do
  return
end subroutine fftback

! ***********************************************************
!             Sous-programme utilis par fftback
! ***********************************************************
subroutine cfftb (n,c,wsave)
  implicit none
  integer n,iw1,iw2
  real*8 c(*),wsave(*)
  !
  !
  if (n .eq. 1)  return
  iw1 = n + n + 1
  iw2 = iw1 + n + n
  call cfftb1 (n,c,wsave,wsave(iw1),wsave(iw2))
  return
end subroutine cfftb

! ***********************************************************
!             Sous-programme utilis par cfftb
! ***********************************************************
subroutine cfftb1 (n,c,ch,wa,ifac)
  implicit none
  integer i,k1,l1,l2
  integer na,nac,nf,n2
  integer ido,idot,idl1,ip
  integer iw,ix2,ix3,ix4
  integer n,ifac(*)
  real*8 ch(*),c(*),wa(*)
  !
  !
  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1
  do k1 = 1, nf
     ip = ifac(k1+2)
     l2 = ip * l1
     ido = n / l2
     idot = ido + ido
     idl1 = idot * l1
     if (ip .eq. 5) then
        ix2 = iw + idot
        ix3 = ix2 + idot
        ix4 = ix3 + idot
        if (na .eq. 0) then
           call passb5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
        else
           call passb5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
        end if
        na = 1 - na
     else if (ip .eq. 4) then
        ix2 = iw + idot
        ix3 = ix2 + idot
        if (na .eq. 0) then
           call passb4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
        else
           call passb4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
        end if
        na = 1 - na
     else if (ip .eq. 3) then
        ix2 = iw + idot
        if (na .eq. 0) then
           call passb3 (idot,l1,c,ch,wa(iw),wa(ix2))
        else
           call passb3 (idot,l1,ch,c,wa(iw),wa(ix2))
        end if
        na = 1 - na
     else if (ip .eq. 2) then
        if (na .eq. 0) then
           call passb2 (idot,l1,c,ch,wa(iw))
        else
           call passb2 (idot,l1,ch,c,wa(iw))
        end if
        na = 1 - na
     else
        if (na .eq. 0) then
           call passb (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
        else
           call passb (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
        end if
        if (nac .ne. 0)  na = 1 - na
     end if
     l1 = l2
     iw = iw + (ip-1)*idot
  end do
  if (na .ne. 0) then
     n2 = n + n
     do i = 1, n2
        c(i) = ch(i)
     end do
  end if
  return
end subroutine cfftb1

! ***********************************************************
!             Sous-programme utilis par cfftb1
! ***********************************************************
subroutine passb (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
  implicit none
  integer nac,ido,ip,l1,idl1
  integer i,j,k,l,ik,jc,lc
  integer idj,idl,idp,idij,idlj
  integer inc,idot,nt,ipp2,ipph
  real*8 wai,war,wa(*)
  real*8 cc(ido,ip,l1),c1(ido,l1,ip)
  real*8 ch(ido,l1,ip),ch2(idl1,ip)
  real*8 c2(idl1,ip)
  !
  !
  idot = ido / 2
  nt = ip * idl1
  ipp2 = ip + 2
  ipph = (ip+1) / 2
  idp = ip * ido
  if (ido .ge. l1) then
     do j = 2, ipph
        jc = ipp2 - j
        do k = 1, l1
           do i = 1, ido
              ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
              ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
           end do
        end do
     end do
     do k = 1, l1
        do i = 1, ido
           ch(i,k,1) = cc(i,1,k)
        end do
     end do
  else
     do j = 2, ipph
        jc = ipp2 - j
        do i = 1, ido
           do k = 1, l1
              ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
              ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
           end do
        end do
     end do
     do i = 1, ido
        do k = 1, l1
           ch(i,k,1) = cc(i,1,k)
        end do
     end do
  end if
  idl = 2 - ido
  inc = 0
  do l = 2, ipph
     lc = ipp2 - l
     idl = idl + ido
     do ik = 1, idl1
        c2(ik,l) = ch2(ik,1) + wa(idl-1)*ch2(ik,2)
        c2(ik,lc) = wa(idl) * ch2(ik,ip)
     end do
     idlj = idl
     inc = inc + ido
     do j = 3, ipph
        jc = ipp2 - j
        idlj = idlj + inc
        if (idlj .gt. idp)  idlj = idlj - idp
        war = wa(idlj-1)
        wai = wa(idlj)
        do ik = 1, idl1
           c2(ik,l) = c2(ik,l) + war*ch2(ik,j)
           c2(ik,lc) = c2(ik,lc) + wai*ch2(ik,jc)
        end do
     end do
  end do
  do j = 2, ipph
     do ik = 1, idl1
        ch2(ik,1) = ch2(ik,1) + ch2(ik,j)
     end do
  end do
  do j = 2, ipph
     jc = ipp2 - j
     do ik = 2, idl1, 2
        ch2(ik-1,j) = c2(ik-1,j) - c2(ik,jc)
        ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
        ch2(ik,j) = c2(ik,j) + c2(ik-1,jc)
        ch2(ik,jc) = c2(ik,j) - c2(ik-1,jc)
     end do
  end do
  nac = 1
  if (ido .ne. 2) then
     nac = 0
     do ik = 1, idl1
        c2(ik,1) = ch2(ik,1)
     end do
     do j = 2, ip
        do k = 1, l1
           c1(1,k,j) = ch(1,k,j)
           c1(2,k,j) = ch(2,k,j)
        end do
     end do
     if (idot .le. l1) then
        idij = 0
        do j = 2, ip
           idij = idij + 2
           do i = 4, ido, 2
              idij = idij + 2
              do k = 1, l1
                 c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j) &
                      - wa(idij)*ch(i,k,j)
                 c1(i,k,j) = wa(idij-1)*ch(i,k,j) &
                      + wa(idij)*ch(i-1,k,j)
              end do
           end do
        end do
     else
        idj = 2 - ido
        do j = 2, ip
           idj = idj + ido
           do k = 1, l1
              idij = idj
              do i = 4, ido, 2
                 idij = idij + 2
                 c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j) &
                      - wa(idij)*ch(i,k,j)
                 c1(i,k,j) = wa(idij-1)*ch(i,k,j) &
                      + wa(idij)*ch(i-1,k,j)
              end do
           end do
        end do
     end if
  end if
  return
end subroutine passb

! ***********************************************************
!             Sous-programme utilis par cfftb1
! ***********************************************************
subroutine passb2 (ido,l1,cc,ch,wa1)
  implicit none
  integer i,k,ido,l1
  real*8 ti2,tr2
  real*8 cc(ido,2,l1),ch(ido,l1,2),wa1(*)
  !
  !
  if (ido .le. 2) then
     do k = 1, l1
        ch(1,k,1) = cc(1,1,k) + cc(1,2,k)
        ch(1,k,2) = cc(1,1,k) - cc(1,2,k)
        ch(2,k,1) = cc(2,1,k) + cc(2,2,k)
        ch(2,k,2) = cc(2,1,k) - cc(2,2,k)
     end do
  else
     do k = 1, l1
        do i = 2, ido, 2
           ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
           tr2 = cc(i-1,1,k) - cc(i-1,2,k)
           ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
           ti2 = cc(i,1,k) - cc(i,2,k)
           ch(i,k,2) = wa1(i-1)*ti2 + wa1(i)*tr2
           ch(i-1,k,2) = wa1(i-1)*tr2 - wa1(i)*ti2
        end do
     end do
  end if
  return
end subroutine passb2

! ***********************************************************
!             Sous-programme utilis par cfftb1
! ***********************************************************
subroutine passb3 (ido,l1,cc,ch,wa1,wa2)
  implicit none
  integer i,k,ido,l1
  real*8 ti2,tr2,taur,taui
  real*8 ci2,ci3,cr2,cr3
  real*8 di2,di3,dr2,dr3
  real*8 cc(ido,3,l1),ch(ido,l1,3)
  real*8 wa1(*),wa2(*)
  data taur  / -0.5d0 /
  data taui  / 0.866025403784439d0 /
  !
  !
  if (ido .eq. 2) then
     do k = 1, l1
        tr2 = cc(1,2,k) + cc(1,3,k)
        cr2 = cc(1,1,k) + taur*tr2
        ch(1,k,1) = cc(1,1,k) + tr2
        ti2 = cc(2,2,k) + cc(2,3,k)
        ci2 = cc(2,1,k) + taur*ti2
        ch(2,k,1) = cc(2,1,k) + ti2
        cr3 = taui * (cc(1,2,k)-cc(1,3,k))
        ci3 = taui * (cc(2,2,k)-cc(2,3,k))
        ch(1,k,2) = cr2 - ci3
        ch(1,k,3) = cr2 + ci3
        ch(2,k,2) = ci2 + cr3
        ch(2,k,3) = ci2 - cr3
     end do
  else
     do k = 1, l1
        do i = 2, ido, 2
           tr2 = cc(i-1,2,k) + cc(i-1,3,k)
           cr2 = cc(i-1,1,k) + taur*tr2
           ch(i-1,k,1) = cc(i-1,1,k) + tr2
           ti2 = cc(i,2,k) + cc(i,3,k)
           ci2 = cc(i,1,k) + taur*ti2
           ch(i,k,1) = cc(i,1,k) + ti2
           cr3 = taui * (cc(i-1,2,k)-cc(i-1,3,k))
           ci3 = taui * (cc(i,2,k)-cc(i,3,k))
           dr2 = cr2 - ci3
           dr3 = cr2 + ci3
           di2 = ci2 + cr3
           di3 = ci2 - cr3
           ch(i,k,2) = wa1(i-1)*di2 + wa1(i)*dr2
           ch(i-1,k,2) = wa1(i-1)*dr2 - wa1(i)*di2
           ch(i,k,3) = wa2(i-1)*di3 + wa2(i)*dr3
           ch(i-1,k,3) = wa2(i-1)*dr3 - wa2(i)*di3
        end do
     end do
  end if
  return
end subroutine passb3

! ***********************************************************
!             Sous-programme utilis par cfftb1
! ***********************************************************
subroutine passb4 (ido,l1,cc,ch,wa1,wa2,wa3)
  implicit none
  integer i,k,ido,l1
  real*8 ci2,ci3,ci4
  real*8 cr2,cr3,cr4
  real*8 ti1,ti2,ti3,ti4
  real*8 tr1,tr2,tr3,tr4
  real*8 cc(ido,4,l1),ch(ido,l1,4)
  real*8 wa1(*),wa2(*),wa3(*)
  !
  !
  if (ido .eq. 2) then
     do k = 1, l1
        ti1 = cc(2,1,k) - cc(2,3,k)
        ti2 = cc(2,1,k) + cc(2,3,k)
        tr4 = cc(2,4,k) - cc(2,2,k)
        ti3 = cc(2,2,k) + cc(2,4,k)
        tr1 = cc(1,1,k) - cc(1,3,k)
        tr2 = cc(1,1,k) + cc(1,3,k)
        ti4 = cc(1,2,k) - cc(1,4,k)
        tr3 = cc(1,2,k) + cc(1,4,k)
        ch(1,k,1) = tr2 + tr3
        ch(1,k,3) = tr2 - tr3
        ch(2,k,1) = ti2 + ti3
        ch(2,k,3) = ti2 - ti3
        ch(1,k,2) = tr1 + tr4
        ch(1,k,4) = tr1 - tr4
        ch(2,k,2) = ti1 + ti4
        ch(2,k,4) = ti1 - ti4
     end do
  else
     do k = 1, l1
        do i = 2, ido, 2
           ti1 = cc(i,1,k) - cc(i,3,k)
           ti2 = cc(i,1,k) + cc(i,3,k)
           ti3 = cc(i,2,k) + cc(i,4,k)
           tr4 = cc(i,4,k) - cc(i,2,k)
           tr1 = cc(i-1,1,k) - cc(i-1,3,k)
           tr2 = cc(i-1,1,k) + cc(i-1,3,k)
           ti4 = cc(i-1,2,k) - cc(i-1,4,k)
           tr3 = cc(i-1,2,k) + cc(i-1,4,k)
           ch(i-1,k,1) = tr2 + tr3
           cr3 = tr2 - tr3
           ch(i,k,1) = ti2 + ti3
           ci3 = ti2 - ti3
           cr2 = tr1 + tr4
           cr4 = tr1 - tr4
           ci2 = ti1 + ti4
           ci4 = ti1 - ti4
           ch(i-1,k,2) = wa1(i-1)*cr2 - wa1(i)*ci2
           ch(i,k,2) = wa1(i-1)*ci2 + wa1(i)*cr2
           ch(i-1,k,3) = wa2(i-1)*cr3 - wa2(i)*ci3
           ch(i,k,3) = wa2(i-1)*ci3 + wa2(i)*cr3
           ch(i-1,k,4) = wa3(i-1)*cr4 - wa3(i)*ci4
           ch(i,k,4) = wa3(i-1)*ci4 + wa3(i)*cr4
        end do
     end do
  end if
  return
end subroutine passb4

! ***********************************************************
!             Sous-programme utilis par cfftb1
! ***********************************************************
subroutine passb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
  implicit none
  integer i,k,ido,l1
  real*8 ci2,ci3,ci4,ci5
  real*8 cr2,cr3,cr4,cr5
  real*8 di2,di3,di4,di5
  real*8 dr2,dr3,dr4,dr5
  real*8 ti2,ti3,ti4,ti5
  real*8 tr2,tr3,tr4,tr5
  real*8 tr11,ti11,tr12,ti12
  real*8 cc(ido,5,l1),ch(ido,l1,5)
  real*8 wa1(*),wa2(*),wa3(*),wa4(*)
  data tr11  /  0.309016994374947d0 /
  data ti11  /  0.951056516295154d0 /
  data tr12  / -0.809016994374947d0 /
  data ti12  /  0.587785252292473d0 /
  !
  !
  if (ido .eq. 2) then
     do k = 1, l1
        ti5 = cc(2,2,k) - cc(2,5,k)
        ti2 = cc(2,2,k) + cc(2,5,k)
        ti4 = cc(2,3,k) - cc(2,4,k)
        ti3 = cc(2,3,k) + cc(2,4,k)
        tr5 = cc(1,2,k) - cc(1,5,k)
        tr2 = cc(1,2,k) + cc(1,5,k)
        tr4 = cc(1,3,k) - cc(1,4,k)
        tr3 = cc(1,3,k) + cc(1,4,k)
        ch(1,k,1) = cc(1,1,k) + tr2 + tr3
        ch(2,k,1) = cc(2,1,k) + ti2 + ti3
        cr2 = cc(1,1,k) + tr11*tr2 + tr12*tr3
        ci2 = cc(2,1,k) + tr11*ti2 + tr12*ti3
        cr3 = cc(1,1,k) + tr12*tr2 + tr11*tr3
        ci3 = cc(2,1,k) + tr12*ti2 + tr11*ti3
        cr5 = ti11*tr5 + ti12*tr4
        ci5 = ti11*ti5 + ti12*ti4
        cr4 = ti12*tr5 - ti11*tr4
        ci4 = ti12*ti5 - ti11*ti4
        ch(1,k,2) = cr2 - ci5
        ch(1,k,5) = cr2 + ci5
        ch(2,k,2) = ci2 + cr5
        ch(2,k,3) = ci3 + cr4
        ch(1,k,3) = cr3 - ci4
        ch(1,k,4) = cr3 + ci4
        ch(2,k,4) = ci3 - cr4
        ch(2,k,5) = ci2 - cr5
     end do
  else
     do k = 1, l1
        do i = 2, ido, 2
           ti5 = cc(i,2,k) - cc(i,5,k)
           ti2 = cc(i,2,k) + cc(i,5,k)
           ti4 = cc(i,3,k) - cc(i,4,k)
           ti3 = cc(i,3,k) + cc(i,4,k)
           tr5 = cc(i-1,2,k) - cc(i-1,5,k)
           tr2 = cc(i-1,2,k) + cc(i-1,5,k)
           tr4 = cc(i-1,3,k) - cc(i-1,4,k)
           tr3 = cc(i-1,3,k) + cc(i-1,4,k)
           ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
           ch(i,k,1) = cc(i,1,k) + ti2 + ti3
           cr2 = cc(i-1,1,k) + tr11*tr2 + tr12*tr3
           ci2 = cc(i,1,k) + tr11*ti2 + tr12*ti3
           cr3 = cc(i-1,1,k) + tr12*tr2 + tr11*tr3
           ci3 = cc(i,1,k) + tr12*ti2 + tr11*ti3
           cr5 = ti11*tr5 + ti12*tr4
           ci5 = ti11*ti5 + ti12*ti4
           cr4 = ti12*tr5 - ti11*tr4
           ci4 = ti12*ti5 - ti11*ti4
           dr3 = cr3 - ci4
           dr4 = cr3 + ci4
           di3 = ci3 + cr4
           di4 = ci3 - cr4
           dr5 = cr2 + ci5
           dr2 = cr2 - ci5
           di5 = ci2 - cr5
           di2 = ci2 + cr5
           ch(i-1,k,2) = wa1(i-1)*dr2 - wa1(i)*di2
           ch(i,k,2) = wa1(i-1)*di2 + wa1(i)*dr2
           ch(i-1,k,3) = wa2(i-1)*dr3 - wa2(i)*di3
           ch(i,k,3) = wa2(i-1)*di3 + wa2(i)*dr3
           ch(i-1,k,4) = wa3(i-1)*dr4 - wa3(i)*di4
           ch(i,k,4) = wa3(i-1)*di4 + wa3(i)*dr4
           ch(i-1,k,5) = wa4(i-1)*dr5 - wa4(i)*di5
           ch(i,k,5) = wa4(i-1)*di5 + wa4(i)*dr5
        end do
     end do
  end if
  return
end subroutine passb5

! ***********************************************************
!             Sous-programme utilis pour la FFT directe
! ***********************************************************
subroutine cfftf (n,c,wsave)
  implicit none
  integer i,n,iw1,iw2
  real*8 wsave(*)
  double complex c(*)
  IF (N .EQ. 1) RETURN
  IW1 = N+N+1
  IW2 = IW1+N+N
  CALL CFFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
  return
end subroutine cfftf

! ***********************************************************
!       Sous-programme utilis pour la subroutine cfftf
! ***********************************************************
subroutine cfftf1 (N,C,CH,WA,IFAC)
  implicit none
  integer i,k1,l1,l2
  integer na,nac,nf,n2
  integer ido,idot,idl1,ip
  integer iw,ix2,ix3,ix4
  integer n,ifac(*)
  real*8 ch(*),c(*),wa(*)
  NF = IFAC(2)
  NA = 0
  L1 = 1
  IW = 1

  DO 116 K1=1,NF
     IP = IFAC(K1+2)
     L2 = IP*L1
     IDO = N/L2
     IDOT = IDO+IDO
     IDL1 = IDOT*L1
     IF (IP .NE. 4) GO TO 103
     IX2 = IW+IDOT
     IX3 = IX2+IDOT
     IF (NA .NE. 0) GO TO 101
     call passf4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
     GO TO 102
101  CONTINUE
     call passf4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
102  NA = 1-NA
     GO TO 115
103  IF (IP .NE. 2) GO TO 106
     IF (NA .NE. 0) GO TO 104
     call passf2 (IDOT,L1,C,CH,WA(IW))
     GO TO 105
104  CONTINUE
     call passf2 (IDOT,L1,CH,C,WA(IW))
105  NA = 1-NA
     GO TO 115
106  IF (IP .NE. 3) GO TO 109
     IX2 = IW+IDOT
     IF (NA .NE. 0) GO TO 107
     call passf3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
     GO TO 108
107  CONTINUE
     call passf3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
108  NA = 1-NA
     GO TO 115
109  IF (IP .NE. 5) GO TO 112
     IX2 = IW+IDOT
     IX3 = IX2+IDOT
     IX4 = IX3+IDOT
     IF (NA .NE. 0) GO TO 110
     call passf5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
     GO TO 111
110  CONTINUE
     call passf5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
111  NA = 1-NA
     GO TO 115
112  IF (NA .NE. 0) GO TO 113
     call passf (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
     GO TO 114
113  CONTINUE
     call passf (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
114  IF (NAC .NE. 0) NA = 1-NA
115  L1 = L2
     IW = IW+(IP-1)*IDOT
116  CONTINUE
     IF (NA .EQ. 0) RETURN
     N2 = N+N
     DO 117 I=1,N2
        C(I) = CH(I)
117     CONTINUE
        return
      end subroutine cfftf1

      ! ***********************************************************
      !       Sous-programme utilis pour la subroutine cfftf1
      ! ***********************************************************
      subroutine passf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
        implicit none
        integer i,k,ido,l1
        real*8 ci2,ci3,ci4,ci5
        real*8 cr2,cr3,cr4,cr5
        real*8 di2,di3,di4,di5
        real*8 dr2,dr3,dr4,dr5
        real*8 ti2,ti3,ti4,ti5
        real*8 tr2,tr3,tr4,tr5
        real*8 tr11,ti11,tr12,ti12
        real*8 cc(ido,5,l1),ch(ido,l1,5)
        real*8 wa1(*),wa2(*),wa3(*),wa4(*)
        data tr11  /  0.309016994374947d0 /
        data ti11  / -0.951056516295154d0 /
        data tr12  / -0.809016994374947d0 /
        data ti12  / -0.587785252292473d0 /


        if (ido .eq. 2) then
           do k = 1, l1
              ti5 = cc(2,2,k) - cc(2,5,k)
              ti2 = cc(2,2,k) + cc(2,5,k)
              ti4 = cc(2,3,k) - cc(2,4,k)
              ti3 = cc(2,3,k) + cc(2,4,k)
              tr5 = cc(1,2,k) - cc(1,5,k)
              tr2 = cc(1,2,k) + cc(1,5,k)
              tr4 = cc(1,3,k) - cc(1,4,k)
              tr3 = cc(1,3,k) + cc(1,4,k)
              ch(1,k,1) = cc(1,1,k) + tr2 + tr3
              ch(2,k,1) = cc(2,1,k) + ti2 + ti3
              cr2 = cc(1,1,k) + tr11*tr2 + tr12*tr3
              ci2 = cc(2,1,k) + tr11*ti2 + tr12*ti3
              cr3 = cc(1,1,k) + tr12*tr2 + tr11*tr3
              ci3 = cc(2,1,k) + tr12*ti2 + tr11*ti3
              cr5 = ti11*tr5 + ti12*tr4
              ci5 = ti11*ti5 + ti12*ti4
              cr4 = ti12*tr5 - ti11*tr4
              ci4 = ti12*ti5 - ti11*ti4
              ch(1,k,2) = cr2 - ci5
              ch(1,k,5) = cr2 + ci5
              ch(2,k,2) = ci2 + cr5
              ch(2,k,3) = ci3 + cr4
              ch(1,k,3) = cr3 - ci4
              ch(1,k,4) = cr3 + ci4
              ch(2,k,4) = ci3 - cr4
              ch(2,k,5) = ci2 - cr5
           enddo
        else
           do k = 1, l1
              do i = 2, ido, 2
                 ti5 = cc(i,2,k) - cc(i,5,k)
                 ti2 = cc(i,2,k) + cc(i,5,k)
                 ti4 = cc(i,3,k) - cc(i,4,k)
                 ti3 = cc(i,3,k) + cc(i,4,k)
                 tr5 = cc(i-1,2,k) - cc(i-1,5,k)
                 tr2 = cc(i-1,2,k) + cc(i-1,5,k)
                 tr4 = cc(i-1,3,k) - cc(i-1,4,k)
                 tr3 = cc(i-1,3,k) + cc(i-1,4,k)
                 ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
                 ch(i,k,1) = cc(i,1,k) + ti2 + ti3
                 cr2 = cc(i-1,1,k) + tr11*tr2 + tr12*tr3
                 ci2 = cc(i,1,k) + tr11*ti2 + tr12*ti3
                 cr3 = cc(i-1,1,k) + tr12*tr2 + tr11*tr3
                 ci3 = cc(i,1,k) + tr12*ti2 + tr11*ti3
                 cr5 = ti11*tr5 + ti12*tr4
                 ci5 = ti11*ti5 + ti12*ti4
                 cr4 = ti12*tr5 - ti11*tr4
                 ci4 = ti12*ti5 - ti11*ti4
                 dr3 = cr3 - ci4
                 dr4 = cr3 + ci4
                 di3 = ci3 + cr4
                 di4 = ci3 - cr4
                 dr5 = cr2 + ci5
                 dr2 = cr2 - ci5
                 di5 = ci2 - cr5
                 di2 = ci2 + cr5
                 ch(i-1,k,2) = wa1(i-1)*dr2 + wa1(i)*di2
                 ch(i,k,2) = wa1(i-1)*di2 - wa1(i)*dr2
                 ch(i-1,k,3) = wa2(i-1)*dr3 + wa2(i)*di3
                 ch(i,k,3) = wa2(i-1)*di3 - wa2(i)*dr3
                 ch(i-1,k,4) = wa3(i-1)*dr4 + wa3(i)*di4
                 ch(i,k,4) = wa3(i-1)*di4 - wa3(i)*dr4
                 ch(i-1,k,5) = wa4(i-1)*dr5 + wa4(i)*di5
                 ch(i,k,5) = wa4(i-1)*di5 - wa4(i)*dr5
              enddo
           enddo
        endif
        return
      end subroutine passf5

      ! ***********************************************************
      !       Sous-programme utilis pour la subroutine cfftf1
      ! ***********************************************************
      subroutine passf4 (ido,l1,cc,ch,wa1,wa2,wa3)
        implicit none
        integer i,k,ido,l1
        real*8 ci2,ci3,ci4
        real*8 cr2,cr3,cr4
        real*8 ti1,ti2,ti3,ti4
        real*8 tr1,tr2,tr3,tr4
        real*8 cc(ido,4,l1),ch(ido,l1,4)
        real*8 wa1(*),wa2(*),wa3(*)

        if (ido .eq. 2) then
           do k = 1, l1
              ti1 = cc(2,1,k) - cc(2,3,k)
              ti2 = cc(2,1,k) + cc(2,3,k)
              tr4 = cc(2,2,k) - cc(2,4,k)
              ti3 = cc(2,2,k) + cc(2,4,k)
              tr1 = cc(1,1,k) - cc(1,3,k)
              tr2 = cc(1,1,k) + cc(1,3,k)
              ti4 = cc(1,4,k) - cc(1,2,k)
              tr3 = cc(1,2,k) + cc(1,4,k)
              ch(1,k,1) = tr2 + tr3
              ch(1,k,3) = tr2 - tr3
              ch(2,k,1) = ti2 + ti3
              ch(2,k,3) = ti2 - ti3
              ch(1,k,2) = tr1 + tr4
              ch(1,k,4) = tr1 - tr4
              ch(2,k,2) = ti1 + ti4
              ch(2,k,4) = ti1 - ti4
           end do
        else
           do k = 1, l1
              do i = 2, ido, 2
                 ti1 = cc(i,1,k) - cc(i,3,k)
                 ti2 = cc(i,1,k) + cc(i,3,k)
                 ti3 = cc(i,2,k) + cc(i,4,k)
                 tr4 = cc(i,2,k) - cc(i,4,k)
                 tr1 = cc(i-1,1,k) - cc(i-1,3,k)
                 tr2 = cc(i-1,1,k) + cc(i-1,3,k)
                 ti4 = cc(i-1,4,k) - cc(i-1,2,k)
                 tr3 = cc(i-1,2,k) + cc(i-1,4,k)
                 ch(i-1,k,1) = tr2 + tr3
                 cr3 = tr2 - tr3
                 ch(i,k,1) = ti2 + ti3
                 ci3 = ti2 - ti3
                 cr2 = tr1 + tr4
                 cr4 = tr1 - tr4
                 ci2 = ti1 + ti4
                 ci4 = ti1 - ti4
                 ch(i-1,k,2) = wa1(i-1)*cr2 + wa1(i)*ci2
                 ch(i,k,2) = wa1(i-1)*ci2 - wa1(i)*cr2
                 ch(i-1,k,3) = wa2(i-1)*cr3 + wa2(i)*ci3
                 ch(i,k,3) = wa2(i-1)*ci3 - wa2(i)*cr3
                 ch(i-1,k,4) = wa3(i-1)*cr4 + wa3(i)*ci4
                 ch(i,k,4) = wa3(i-1)*ci4 - wa3(i)*cr4
              end do
           end do
        end if
        return
      end subroutine passf4

      ! ***********************************************************
      !       Sous-programme utilis pour la subroutine cfftf1
      ! ***********************************************************
      subroutine passf3 (ido,l1,cc,ch,wa1,wa2)
        implicit none
        integer i,k,ido,l1
        real*8 ti2,tr2,taur,taui
        real*8 ci2,ci3,cr2,cr3
        real*8 di2,di3,dr2,dr3
        real*8 cc(ido,3,l1),ch(ido,l1,3)
        real*8 wa1(*),wa2(*)
        data taur  / -0.5d0 /
        data taui  / -0.866025403784439d0 /

        if (ido .eq. 2) then
           do k = 1, l1
              tr2 = cc(1,2,k) + cc(1,3,k)
              cr2 = cc(1,1,k) + taur*tr2
              ch(1,k,1) = cc(1,1,k) + tr2
              ti2 = cc(2,2,k) + cc(2,3,k)
              ci2 = cc(2,1,k) + taur*ti2
              ch(2,k,1) = cc(2,1,k) + ti2
              cr3 = taui * (cc(1,2,k)-cc(1,3,k))
              ci3 = taui * (cc(2,2,k)-cc(2,3,k))
              ch(1,k,2) = cr2 - ci3
              ch(1,k,3) = cr2 + ci3
              ch(2,k,2) = ci2 + cr3
              ch(2,k,3) = ci2 - cr3
           enddo
        else
           do k = 1, l1
              do i = 2, ido, 2
                 tr2 = cc(i-1,2,k) + cc(i-1,3,k)
                 cr2 = cc(i-1,1,k) + taur*tr2
                 ch(i-1,k,1) = cc(i-1,1,k) + tr2
                 ti2 = cc(i,2,k) + cc(i,3,k)
                 ci2 = cc(i,1,k) + taur*ti2
                 ch(i,k,1) = cc(i,1,k) + ti2
                 cr3 = taui * (cc(i-1,2,k)-cc(i-1,3,k))
                 ci3 = taui * (cc(i,2,k)-cc(i,3,k))
                 dr2 = cr2 - ci3
                 dr3 = cr2 + ci3
                 di2 = ci2 + cr3
                 di3 = ci2 - cr3
                 ch(i,k,2) = wa1(i-1)*di2 - wa1(i)*dr2
                 ch(i-1,k,2) = wa1(i-1)*dr2 + wa1(i)*di2
                 ch(i,k,3) = wa2(i-1)*di3 - wa2(i)*dr3
                 ch(i-1,k,3) = wa2(i-1)*dr3 + wa2(i)*di3
              enddo
           enddo
        endif
        return
      end subroutine passf3

      ! ***********************************************************
      !       Sous-programme utilis pour la subroutine cfftf1
      ! ***********************************************************
      subroutine passf2 (ido,l1,cc,ch,wa1)
        implicit none
        integer i,k,ido,l1
        real*8 ti2,tr2
        real*8 cc(ido,2,l1),ch(ido,l1,2),wa1(*)

        if (ido .le. 2) then
           do k = 1, l1
              ch(1,k,1) = cc(1,1,k) + cc(1,2,k)
              ch(1,k,2) = cc(1,1,k) - cc(1,2,k)
              ch(2,k,1) = cc(2,1,k) + cc(2,2,k)
              ch(2,k,2) = cc(2,1,k) - cc(2,2,k)
           enddo
        else
           do k = 1, l1
              do i = 2, ido, 2
                 ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
                 tr2 = cc(i-1,1,k) - cc(i-1,2,k)
                 ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
                 ti2 = cc(i,1,k) - cc(i,2,k)
                 ch(i,k,2) = wa1(i-1)*ti2 - wa1(i)*tr2
                 ch(i-1,k,2) = wa1(i-1)*tr2 + wa1(i)*ti2
              enddo
           enddo
        endif
        return
      end subroutine passf2

      ! ***********************************************************
      !       Sous-programme utilis pour la subroutine cfftf1
      ! ***********************************************************
      subroutine passf (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
        implicit none
        integer nac,ido,ip,l1,idl1
        integer i,j,k,l,ik,jc,lc
        integer idj,idl,idp,idij,idlj
        integer inc,idot,nt,ipp2,ipph
        real*8 wai,war,wa(*)
        real*8 cc(ido,ip,l1),c1(ido,l1,ip)
        real*8 ch(ido,l1,ip),ch2(idl1,ip)
        real*8 c2(idl1,ip)

        idot = ido / 2
        nt = ip * idl1
        ipp2 = ip + 2
        ipph = (ip+1) / 2
        idp = ip * ido
        if (ido .ge. l1) then
           do j = 2, ipph
              jc = ipp2 - j
              do k = 1, l1
                 do i = 1, ido
                    ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
                    ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
                 enddo
              enddo
           enddo
           do k = 1, l1
              do i = 1, ido
                 ch(i,k,1) = cc(i,1,k)
              enddo
           enddo
        else
           do j = 2, ipph
              jc = ipp2 - j
              do i = 1, ido
                 do k = 1, l1
                    ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
                    ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
                 enddo
              enddo
           enddo
           do i = 1, ido
              do k = 1, l1
                 ch(i,k,1) = cc(i,1,k)
              enddo
           enddo
        endif
        idl = 2 - ido
        inc = 0
        do l = 2, ipph
           lc = ipp2 - l
           idl = idl + ido
           do ik = 1, idl1
              c2(ik,l) = ch2(ik,1) + wa(idl-1)*ch2(ik,2)
              c2(ik,lc) = -wa(idl) * ch2(ik,ip)
           enddo
           idlj = idl
           inc = inc + ido
           do j = 3, ipph
              jc = ipp2 - j
              idlj = idlj + inc
              if (idlj .gt. idp)  idlj = idlj - idp
              war = wa(idlj-1)
              wai = wa(idlj)
              do ik = 1, idl1
                 c2(ik,l) = c2(ik,l) + war*ch2(ik,j)
                 c2(ik,lc) = c2(ik,lc) - wai*ch2(ik,jc)
              enddo
           enddo
        enddo
        do j = 2, ipph
           do ik = 1, idl1
              ch2(ik,1) = ch2(ik,1) + ch2(ik,j)
           enddo
        enddo
        do j = 2, ipph
           jc = ipp2 - j
           do ik = 2, idl1, 2
              ch2(ik-1,j) = c2(ik-1,j) - c2(ik,jc)
              ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
              ch2(ik,j) = c2(ik,j) + c2(ik-1,jc)
              ch2(ik,jc) = c2(ik,j) - c2(ik-1,jc)
           enddo
        enddo
        nac = 1
        if (ido .ne. 2) then
           nac = 0
           do ik = 1, idl1
              c2(ik,1) = ch2(ik,1)
           enddo
           do j = 2, ip
              do k = 1, l1
                 c1(1,k,j) = ch(1,k,j)
                 c1(2,k,j) = ch(2,k,j)
              enddo
           enddo
           if (idot .le. l1) then
              idij = 0
              do j = 2, ip
                 idij = idij + 2
                 do i = 4, ido, 2
                    idij = idij + 2
                    do k = 1, l1
                       c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+ wa(idij)*ch(i,k,j)
                       c1(i,k,j) = wa(idij-1)*ch(i,k,j) - wa(idij)*ch(i-1,k,j)
                    enddo
                 enddo
              enddo
           else
              idj = 2 - ido
              do j = 2, ip
                 idj = idj + ido
                 do k = 1, l1
                    idij = idj
                    do i = 4, ido, 2
                       idij = idij + 2
                       c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j) + wa(idij)*ch(i,k,j)
                       c1(i,k,j) = wa(idij-1)*ch(i,k,j) - wa(idij)*ch(i-1,k,j)
                    enddo
                 enddo
              enddo
           endif
        endif
        return
      end subroutine passf
