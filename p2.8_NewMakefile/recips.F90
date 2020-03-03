module recips_mod
        implicit none
        contains

!
!---------------------------------------------------------------------
subroutine recips(a1, a2, a3, b1, b2, b3)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  !---------------------------------------------------------------------
  !
  !   This routine generates the reciprocal lattice vectors b1,b2,b3
  !   given the real space vectors a1,a2,a3. The b's are units of 2 pi/a.
  !
  !
  !     first the input variables
  !
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(double) , intent(in) :: a1(3)
  real(double) , intent(in) :: a2(3)
  real(double) , intent(in) :: a3(3)
  real(double) , intent(out) :: b1(3)
  real(double) , intent(out) :: b2(3)
  real(double) , intent(out) :: b3(3)
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: iperm, i, j, k, l, ipol
  real(double) :: den, s
  !-----------------------------------------------
  !
  !   then the local variables
  !
  !
  !    first we compute the denominator
  !
  den = 0
  i = 1
  j = 2
  k = 3
  s = 1.D0
  do iperm = 1, 3
     den = den+s*a1(i)*a2(j)*a3(k)
     l = i
     i = j
     j = k
     k = l
  end do
  i = 2
  j = 1
  k = 3
  s = -s
  do while(s<0.D0)
     do iperm = 1, 3
        den = den+s*a1(i)*a2(j)*a3(k)
        l = i
        i = j
        j = k
        k = l
     end do
     i = 2
     j = 1
     k = 3
     s = -s
  end do
  !
  !    here we compute the reciprocal vectors
  !
  i = 1
  j = 2
  k = 3
  do ipol = 1, 3
     b1(ipol) = (a2(j)*a3(k)-a2(k)*a3(j))/den
     b2(ipol) = (a3(j)*a1(k)-a3(k)*a1(j))/den
     b3(ipol) = (a1(j)*a2(k)-a1(k)*a2(j))/den
     l = i
     i = j
     j = k
     k = l
  end do
  return
end subroutine recips



!*******************************************************************
real(kind(0.0d0)) function distmin (a, b)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(double) , intent(in) :: a(3)
  real(double) , intent(in) :: b(3)
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  real(double), dimension(3) :: vp
  real(double) :: nvp, na, nb
  !-----------------------------------------------
  na = sqrt(a(1)**2+a(2)**2+a(3)**2)
  nb = sqrt(b(1)**2+b(2)**2+b(3)**2)
  vp(1) = a(2)*b(3)-a(3)*b(2)
  vp(2) = a(3)*b(1)-a(1)*b(3)
  vp(3) = a(1)*b(2)-a(2)*b(1)
  nvp = sqrt(vp(1)**2+vp(2)**2+vp(3)**2)
  distmin = nvp/max(na,nb)
  return
end function distmin



!*********************************************************************
real(kind(0.0d0)) function calcvol (a1, a2, a3)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(double) , intent(in) :: a1(3)
  real(double) , intent(in) :: a2(3)
  real(double) , intent(in) :: a3(3)
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, j, k, l, s, iperm
  !-----------------------------------------------
  calcvol = 0.0
  i = 1
  j = 2
  k = 3
  s = 1.D0
  do iperm = 1, 3
     calcvol = calcvol+s*a1(i)*a2(j)*a3(k)
     l = i
     i = j
     j = k
     k = l
  end do
  i = 2
  j = 1
  k = 3
  s = -s
  do while(s<0.D0)
     do iperm = 1, 3
        calcvol = calcvol+s*a1(i)*a2(j)*a3(k)
        l = i
        i = j
        j = k
        k = l
     end do
     i = 2
     j = 1
     k = 3
     s = -s
  end do
  calcvol=dabs(calcvol)
  return
end function calcvol
end module
