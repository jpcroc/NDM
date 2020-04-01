! ***********************************************************
!    Sous-programme preparation des tableaux bsarray () et
!                   array() pour la PME
!                Version du 10/12/2001
! ***********************************************************
subroutine moduli
  !----------------------------------------------
  !   M o d u l e s
  !----------------------------------------------
  USE T_kind_param_m
  use gen_com_m
  use var_pot
  implicit none
  integer i
  real(double) w,bsarray(kpme)
  real(double) array(maxorder)
  real(double) darray(maxorder)
  !
  !     compute and load the moduli values
  !
  w = 0.0d0
  call bspline1 (w,maxorder,array,darray)
  do i = 1, kpme
     bsarray(i) = 0.0d0
  end do
  do i = 2, maxorder+1
     bsarray(i) = array(i-1)
  end do
  call dftmod (bsmod1,bsarray,kpmex)
  call dftmod (bsmod2,bsarray,kpmey)
  call dftmod (bsmod3,bsarray,kpmez)

  return
end subroutine moduli

! ***********************************************************
!    Calcul des coefficients pour l'approximation Cardinal
!                    B_spline
!                Version du 10/12/2001
! ***********************************************************
subroutine bspline1 (x,n,c,d)
  USE T_kind_param_m
  implicit none
  integer i,k,n
  real(double) :: x,denom
  real(double) :: c(n),d(n)
  !
  !     initialize the B-spline as the linear case
  !
  c(n) = 0.0d0
  c(2) = x
  c(1) = 1.0d0 - x
  !
  !     compute standard B-spline recursion to order "n-1"
  !
  do k = 3, n-1
     denom = 1.0d0 / dble(k-1)
     c(k) = x * c(k-1) * denom
     do i = 1, k-2
        c(k-i) = ((x+dble(i))*c(k-i-1) &
             + (dble(k-i)-x)*c(k-i)) * denom
     end do
     c(1) = (1.0d0-x) * c(1) * denom
  end do
  !
  !     get the derivative from "n-1"-th order coefficients
  !
  d(1) = -c(1)
  do i = 2, n
     d(i) = c(i-1) - c(i)
  end do
  !
  !     use one final recursion to get to "n"-th order
  !
  denom = 1.0d0 / dble(n-1)
  c(n) = x * c(n-1) * denom
  do i = 1, n-2
     c(n-i) = ((x+dble(i))*c(n-i-1) &
          + (dble(n-i)-x)*c(n-i)) * denom
  end do
  c(1) = (1.0d0-x) * c(1) * denom
  return
end subroutine bspline1


! ***********************************************************
!    Calcul des coefficients pour l'approximation Cardinal
!                         B_spline
!                    Version decembre 2001
! ***********************************************************
subroutine bspline (nbatom,x,n,c,d)
  USE T_kind_param_m
  implicit none
  integer i,k,n,i1,nbatom
  real(double) :: x(nbatom),denom
  real(double) :: c(n,nbatom),d(n,nbatom)
  !
  !     initialize the B-spline as the linear case
  !

  do i=1,nbatom
     c(n,i) = 0.0d0
     c(2,i) = x(i)
     c(1,i) = 1.0d0 - x(i)
  enddo
  !
  !     compute standard B-spline recursion to order "n-1"
  !
  do k = 3, n-1
     denom = 1.0d0 / dble(k-1)
     do i1=1,nbatom
        c(k,i1) = x(i1) * c(k-1,i1) * denom
        do i = 1, k-2
           c(k-i,i1) = ((x(i1)+dble(i))*c(k-i-1,i1) &
                + (dble(k-i)-x(i1))*c(k-i,i1)) * denom
        end do
        c(1,i1) = (1.0d0-x(i1)) * c(1,i1) * denom
     end do
  end do
  !
  !     get the derivative from "n-1"-th order coefficients
  !
  do i1=1,nbatom
     d(1,i1) = -c(1,i1)
     do i = 2, n
        d(i,i1) = c(i-1,i1) - c(i,i1)
     end do
  enddo
  !
  !     use one final recursion to get to "n"-th order
  !
  denom = 1.0d0 / dble(n-1)
  do i1=1,nbatom
     c(n,i1) = x(i1) * c(n-1,i1) * denom
     do i = 1, n-2
        c(n-i,i1) = ((x(i1)+dble(i))*c(n-i-1,i1) &
             + (dble(n-i)-x(i1))*c(n-i,i1)) * denom
     end do
     c(1,i1) = (1.0d0-x(i1)) * c(1,i1) * denom
  enddo
  return
end subroutine bspline


! ***************************************************
!           Initialisation des tableaux bsmod()
!                Version de decembre 2001
! ***************************************************
subroutine dftmod (bsmod,bsarray,nfft)
  USE T_kind_param_m
  implicit none
  integer i,j,nfft
  real(double) :: eps,factor
  real(double) :: arg,sum1,sum2
  real(double) :: bsmod(*),bsarray(*)
  real(double) :: pi
  !
  !     get the modulus of the discrete Fourier transform
  !
  eps = 1.0d-7
  pi=3.141592654d0

  factor = 2.0d0 * pi / dble(nfft)
  do i = 1, nfft
     sum1 = 0.0d0
     sum2 = 0.0d0
     do j = 1, nfft
        arg = factor * dble((i-1)*(j-1))
        sum1 = sum1 + bsarray(j)*cos(arg)
        sum2 = sum2 + bsarray(j)*sin(arg)
     end do
     bsmod(i) = sum1**2 + sum2**2
  end do
  do i = 1, nfft
     if (bsmod(i) .lt. eps) &
          bsmod(i) = 0.5d0 * (bsmod(i-1)+bsmod(i+1))
  end do
  return
end subroutine dftmod
