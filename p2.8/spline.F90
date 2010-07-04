

subroutine cspline(n, x, y, b, c, d)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer , intent(in) :: n
  real(double) , intent(in) :: x(n)
  real(double) , intent(in) :: y(n)
  real(double) , intent(inout) :: b(n)
  real(double) , intent(inout) :: c(n)
  real(double) , intent(inout) :: d(n)
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: nm1, ib, i
  real(double) :: t
  !-----------------------------------------------
  !
  !  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
  !  for a cubic interpolating spline
  !
  !    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
  !
  !    for  x(i) .le. x .le. x(i+1)
  !
  !  input..
  !
  !    n = the number of data points or knots (n.ge.2)
  !    x = the abscissas of the knots in strictly increasing order
  !    y = the ordinates of the knots
  !
  !  output..
  !
  !    b, c, d  = arrays of spline coefficients as defined above.
  !
  !  using  p  to denote differentiation,
  !
  !    y(i) = s(x(i))
  !    b(i) = sp(x(i))
  !    c(i) = spp(x(i))/2
  !    d(i) = sppp(x(i))/6  (derivative from the right)
  !
  !  the accompanying function subprogram  seval  can be used
  !  to evaluate the spline.
  !
  !
  !
! write(6,*)'spline'
  nm1 = n-1
  if (n<2) return
  if (n>=3) then
     !
     !  set up tridiagonal system
     !
     !  b = diagonal, d = offdiagonal, c = right hand side.
     !
     d(1) = x(2)-x(1)
     c(2) = (y(2)-y(1))/d(1)
     d(2:nm1) = x(3:nm1+1)-x(2:nm1)
     b(2:nm1) = 2.*(d(:nm1-1)+d(2:nm1))
     c(3:nm1+1) = (y(3:nm1+1)-y(2:nm1))/d(2:nm1)
     c(2:nm1) = c(3:nm1+1)-c(2:nm1)
     !
     !  end conditions.  third derivatives at  x(1)  and  x(n)
     !  obtained from divided differences
     !
     b(1) = -d(1)
     b(n) = -d(n-1)
     c(1) = 0.
     c(n) = 0.
     if (n/=3) then
        c(1) = c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
        c(n) = c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
        c(1) = c(1)*d(1)**2/(x(4)-x(1))
        c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
     endif
     !
     !  forward elimination
     !
     do i = 2, n
        t = d(i-1)/b(i-1)
        b(i) = b(i)-t*d(i-1)
        c(i) = c(i)-t*c(i-1)
     end do
     !
     !  back substitution
     !
     c(n) = c(n)/b(n)
     do ib = 1, nm1
        i = n-ib
        c(i) = (c(i)-d(i)*c(i+1))/b(i)
     end do
     !
     !  c(i) is now the sigma(i) of the text
     !
     !  compute polynomial coefficients
     !
     b(n) = (y(n)-y(nm1))/d(nm1)+d(nm1)*(c(nm1)+2.*c(n))
     b(:nm1) = (y(2:nm1+1)-y(:nm1))/d(:nm1)-d(:nm1)*(c(2:nm1+1)+2.*c(:nm1))
     d(:nm1) = (c(2:nm1+1)-c(:nm1))/d(:nm1)
     c(:nm1) = 3.*c(:nm1)
     c(n) = 3.*c(n)
     d(n) = d(n-1)
     return
  endif
  !
  b(1) = (y(2)-y(1))/(x(2)-x(1))
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end subroutine cspline


