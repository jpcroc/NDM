module angular_functions

use ml_in_ndm_module, only : j_max,jj_max, rangml, r_factorial
implicit none

public :: spherical_harm,grad_spherical_harm
!public :: cg
!public :: factorial
contains

   subroutine spherical_4d(x, rr, Umm, dUmm)
   ! Input:
   !    x(3) -> x, y, z ;
   !    rr=sqrt(x**2 + y**2 + z**2)
   ! Output:
   !    Spherical 4D and derivatives
   !    Umm, dUmm
      use ml_in_ndm_module, only : r0,jj_max, one_pi,r_cut, inv_r0
      implicit none
      real(kind=kind(1.d0)),dimension(3),intent(in) :: x
      real(kind=kind(1.d0)),intent(in) :: rr
      double complex,dimension(-jj_max:jj_max,-jj_max:jj_max,0:jj_max), intent(out) :: Umm
      double complex,dimension(-jj_max:jj_max,-jj_max:jj_max,0:jj_max,3), intent(out):: dUmm

      integer :: j,m1,m2
      real(kind=kind(1.d0)) :: z0,  th0, l0, l0i, dz, dil0
      real(kind=kind(1.d0)),   dimension(3) :: drr
      double complex :: cplx_l0i, cplx_i, z_plus,z_minus,x_plus,x_minus, tjm1, tjm2
      double complex,dimension(3) :: dz_plus,dz_minus,dx_plus,dx_minus

      cplx_i=cmplx(0.d0,1.d0, kind=kind(1.d0))
      r0 = r_cut/inv_r0
      th0=rr/r0
      l0= rr/dsin(th0)
      l0i= dsin(th0)/rr
      z0=rr/dtan(th0)
      !debug if ((dabs((rr  - 2.4865502))< 1.d-5) .and.&
      !debug     (dabs((x(1)- 1.5997230))< 1.d-5) .and.&
      !debug     (dabs((x(2)- 1.3288085))< 1.d-5) .and.&
      !debug     (dabs((x(3)- 1.3631164))< 1.d-5) ) then
      !debug  write(*,*) 'aaa', r_cut, inv_r0, r0, th0, l0
      !debug end if

      ! z_+/- is cos(theta0) +/- i sin(theta0) * cos(theta)
      z_plus= cmplx(z0, x(3), kind=kind(1.d0))/l0
      z_minus=cmplx(z0,-x(3), kind=kind(1.d0))/l0
      ! x +/- is sin(theta) * sin(theta0) * exp(+/- i \phi)
      x_plus= cmplx(x(1), x(2), kind=kind(1.d0))/l0
      x_minus=cmplx(x(1),-x(2), kind=kind(1.d0))/l0
      !debug write(*,'("x_plus, x_minus, z_plus, z_minux, rr", 5G15.6)') x_plus, x_minus, z_plus, z_minus, rr


      !Here are the derivatives ...
      drr(1:3)=x(1:3)/rr

      !error dz  = (1.d0/tan(th0) - th0/tan(th0)**2)
      dz  = (1.d0/dtan(th0) - th0/dsin(th0)**2)
      dil0= (dcos(th0)/r0 - l0i)/rr
      cplx_l0i=cmplx(0.d0,l0i, kind=kind(1.d0))

      dz_plus(1:3)= ( cmplx(z0,x(3), kind=kind(1.d0))*dil0 + dz*l0i )*drr(1:3)
      dz_plus(3) = dz_plus(3) + cplx_l0i

      dz_minus(1:3)= ( cmplx(z0,-x(3), kind=kind(1.d0))*dil0 + dz*l0i )*drr(1:3)
      dz_minus(3) = dz_minus(3) - cplx_l0i


      dx_plus(1:3) = cmplx(x(1),x(2), kind=kind(1.d0))*dil0*drr(1:3)
      dx_plus(1) = dx_plus(1) + cmplx(l0i,0.d0, kind=kind(1.d0))
      dx_plus(2) = dx_plus(2) + cplx_l0i

      dx_minus(1:3) = cmplx(x(1),-x(2), kind=kind(1.d0))*dil0*drr(1:3)
      dx_minus(1) = dx_minus(1) + cmplx(l0i,0.d0, kind=kind(1.d0))
      dx_minus(2) = dx_minus(2) - cplx_l0i

      Umm(:,:,:)=(0.d0,0.d0)
      Umm(0,0,0)=(1.d0,0.d0)
      do j=1,jj_max
        do m1=-j,j,2
          ! what is that baby ?
          Umm(m1,m1,j)=(1.d0,0.d0)
        end do
      end do

      dUmm(:,:,:,:)=(0.d0,0.d0)

      do j=1,jj_max
         do m1=-j,j,2
            do m2=-j,j,2
               if (m1.eq.j) then
                 if (m2.eq.-j) then
                   tjm2=cplx_i*dsqrt(dble(j-m2)/dble(j+m1))
                   Umm(m2,m1,j) = - tjm2*x_plus*Umm(m2+1,m1-1,j-1)
                   dUmm(m2,m1,j,:) = - tjm2*( dx_plus(:)*Umm(m2+1,m1-1,j-1)+x_plus*dUmm(m2+1,m1-1,j-1,:) )
                 else if (m2.eq.j) then
                   tjm1=dsqrt(dble(j+m2)/dble(j+m1))
                   Umm(m2,m1,j) = tjm1*z_minus*Umm(m2-1,m1-1,j-1)
                   dUmm(m2,m1,j,:) = tjm1*(dz_minus(:)*Umm(m2-1,m1-1,j-1)+z_minus*dUmm(m2-1,m1-1,j-1,:))
                 else
                   tjm1=dsqrt(dble(j+m2)/dble(j+m1))
                   tjm2=cplx_i*dsqrt(dble(j-m2)/dble(j+m1))
                   Umm(m2,m1,j) = tjm1*z_minus*Umm(m2-1,m1-1,j-1) - tjm2*x_plus*Umm(m2+1,m1-1,j-1)
                   dUmm(m2,m1,j,:) = tjm1*(dz_minus(:)*Umm(m2-1,m1-1,j-1) + z_minus*dUmm(m2-1,m1-1,j-1,:))   &
                                    -tjm2*(dx_plus(:)*Umm(m2+1,m1-1,j-1)+x_plus*dUmm(m2+1,m1-1,j-1,:))
                 end if
               else if (m1.eq.-j) then
                 if (m2.eq.j) then
                   tjm1=cplx_i*dsqrt(dble(j+m2)/dble(j-m1))
                   Umm(m2,m1,j) =  - tjm1*x_minus*Umm(m2-1,m1+1,j-1)
                   !write(*,*) 'ddd', tjm1, x_minus, Umm(m2-1,m1+1,j-1)
                   dUmm(m2,m1,j,:) =  - tjm1*(dx_minus(:)*Umm(m2-1,m1+1,j-1)+x_minus*dUmm(m2-1,m1+1,j-1,:))
                 else if (m2.eq.-j) then
                   tjm2= dsqrt(dble(j-m2)/dble(j-m1))
                   Umm(m2,m1,j) = tjm2*z_plus*Umm(m2+1,m1+1,j-1)
                   !if (j==1) write(*,*)  'ddd', Umm(m2,m1,j)
                   dUmm(m2,m1,j,:) = tjm2*(dz_plus(:)*Umm(m2+1,m1+1,j-1)+z_plus*dUmm(m2+1,m1+1,j-1,:))
                 else

                   tjm1=cplx_i*dsqrt(dble(j+m2)/dble(j-m1))
                   tjm2= dsqrt(dble(j-m2)/dble(j-m1))
                   Umm(m2,m1,j) = tjm2*z_plus*Umm(m2+1,m1+1,j-1) - tjm1*x_minus*Umm(m2-1,m1+1,j-1)

                   dUmm(m2,m1,j,:) = tjm2*(dz_plus(:)*Umm(m2+1,m1+1,j-1)+z_plus*dUmm(m2+1,m1+1,j-1,:)) &
                                 - tjm1*(dx_minus(:)*Umm(m2-1,m1+1,j-1)+x_minus*dUmm(m2-1,m1+1,j-1,:))

                 end if
               else
                 if (m2.eq.-j) then
                   tjm2=cplx_i*dsqrt(dble(j-m2)/dble(j+m1))
                   Umm(m2,m1,j) = - tjm2*x_plus*Umm(m2+1,m1-1,j-1)
                   dUmm(m2,m1,j,:) = - tjm2*( dx_plus(:)*Umm(m2+1,m1-1,j-1)+x_plus*dUmm(m2+1,m1-1,j-1,:) )
                 else if (m2.eq.j) then
                   tjm1=dsqrt(dble(j+m2)/dble(j+m1))
                   Umm(m2,m1,j) = tjm1*z_minus*Umm(m2-1,m1-1,j-1)
                   dUmm(m2,m1,j,:) = tjm1*(dz_minus(:)*Umm(m2-1,m1-1,j-1)+z_minus*dUmm(m2-1,m1-1,j-1,:))
                 else
                   tjm1=dsqrt(dble(j+m2)/dble(j+m1))
                   tjm2=cplx_i*dsqrt(dble(j-m2)/dble(j+m1))
                   Umm(m2,m1,j) = tjm1*z_minus*Umm(m2-1,m1-1,j-1) - tjm2*x_plus*Umm(m2+1,m1-1,j-1)
                   dUmm(m2,m1,j,:) = tjm1*(dz_minus(:)*Umm(m2-1,m1-1,j-1) + z_minus*dUmm(m2-1,m1-1,j-1,:))   &
                                    -tjm2*(dx_plus(:)*Umm(m2+1,m1-1,j-1)+x_plus*dUmm(m2+1,m1-1,j-1,:))
                 end if
               end if
               !debug
               !write(*,'(3i4,2G15.6)') j, m2,m1, Umm(m2,m1,j)
            enddo
         enddo
      enddo
   return

   end subroutine spherical_4d



   function solid_harm(l,m,x)
   ! From Math Wolfram World
   ! http://functions.wolfram.com/PDF/SphericalHarmonicY.pdf

      implicit none

      integer,intent(in) :: l,m
      double precision,dimension(3),intent(in) :: x
      integer :: p,q,s
      double complex :: solid_harm

      solid_harm = (0.d0,0.d0)

      do p=0,l
         q=p-m
         s=l-p-q
         !
         if ((q.ge.0).and.(s.ge.0)) then
            solid_harm = solid_harm + (dcmplx(-0.5d0*x(1),-0.5d0*x(2))**p) * (dcmplx(0.5d0*x(1),-0.5d0*x(2))**q) * (x(3)**s) /&
                                      (r_factorial(p)*r_factorial(q)*r_factorial(s))
         endif
         !
      enddo

      solid_harm = solid_harm * dsqrt(r_factorial(l+m)*r_factorial(l-m))

   end function solid_harm


   function spherical_harm(l,m,x)
   ! From Math Wolfram World
   ! http://functions.wolfram.com/PDF/SphericalHarmonicY.pdf
      use ml_in_ndm_module, only : one_pi
      implicit none
      integer,intent(in) :: l,m
      double precision,dimension(3),intent(in) :: x
      double complex :: spherical_harm

      spherical_harm = solid_harm(l,m,x) * dsqrt((2.d0*dble(l)+1.d0)/(4.d0*one_pi)) * (SUM(x(1:3)**2))**(-0.5d0*l)

   end function spherical_harm


   function grad_spherical_harm(l,m,x)
   ! From Math Wolfram World
   ! http://functions.wolfram.com/PDF/SphericalHarmonicY.pdf
   use ml_in_ndm_module, only : one_pi
   implicit none

   integer,intent(in) :: l,m
   double precision,dimension(3),intent(in) :: x
   double complex,dimension(3) :: grad_spherical_harm
   integer :: p,q,s

   grad_spherical_harm(:)=(0.d0,0.d0)
   do p=0,l
      q=p-m
      s=l-p-q

      if ((p.ge.1).and.(q.ge.0).and.(s.ge.0)) then

         grad_spherical_harm(1) = grad_spherical_harm(1) - 0.5d0 * (dcmplx(-0.5d0*x(1),-0.5d0*x(2))**(p-1)) * (dcmplx(0.5d0*x(1),-0.5d0*x(2))**q) * (x(3)**s) /&
                                                           (r_factorial(p-1)*r_factorial(q)*r_factorial(s))

         grad_spherical_harm(2) = grad_spherical_harm(2) - dcmplx(0.d0,0.5d0) * (dcmplx(-0.5d0*x(1),-0.5d0*x(2))**(p-1)) * (dcmplx(0.5d0*x(1),-0.5d0*x(2))**q) * (x(3)**s) /&
                                                           (r_factorial(p-1)*r_factorial(q)*r_factorial(s))
      endif

      if ((p.ge.0).and.(q.ge.1).and.(s.ge.0)) then

         grad_spherical_harm(1) = grad_spherical_harm(1) + 0.5d0 * (dcmplx(-0.5d0*x(1),-0.5d0*x(2))**p) * (dcmplx(0.5d0*x(1),-0.5d0*x(2))**(q-1)) * (x(3)**s) /&
                                                           (r_factorial(p)*r_factorial(q-1)*r_factorial(s))

         grad_spherical_harm(2) = grad_spherical_harm(2) - dcmplx(0.d0,0.5d0) * (dcmplx(-0.5d0*x(1),-0.5d0*x(2))**p) * (dcmplx(0.5d0*x(1),-0.5d0*x(2))**(q-1)) * (x(3)**s) /&
                                                           (r_factorial(p)*r_factorial(q-1)*r_factorial(s))

      endif

      if ((p.ge.0).and.(q.ge.0).and.(s.ge.1)) then

         grad_spherical_harm(3) = grad_spherical_harm(3) + (dcmplx(-0.5d0*x(1),-0.5d0*x(2))**p) * (dcmplx(0.5d0*x(1),-0.5d0*x(2))**q) * (x(3)**(s-1)) /&
                                                           (r_factorial(p)*r_factorial(q)*r_factorial(s-1))

      endif

   enddo
   grad_spherical_harm(:) = grad_spherical_harm(:) * dsqrt((2.d0*l+1.d0)*r_factorial(l+m)*r_factorial(l-m)/(4.d0*one_pi)) * (sum(x(1:3)**2))**(-0.5d0*l)
   grad_spherical_harm(1:3) = grad_spherical_harm(1:3) - dble(l) * spherical_harm(l,m,x) * x(1:3)/sum(x(1:3)**2)
   !debug write(*,*) grad_spherical_harm(1:3)

   end function grad_spherical_harm


end module angular_functions



function cg(j1,m1,j2,m2,j,m,ntype)
   implicit none

   integer,intent(in) :: j1,j2,j,m1,m2,m,ntype
   double precision :: cg
   double precision ::  wigner3j

   cg = (-1.d0)**((j1-j2+m)/ntype) * dsqrt(2.d0*dble(j)/dble(ntype)+1.d0) * wigner3j(j1,m1,j2,m2,j,-m,ntype)

end function cg


function wigner3j(j1,m1,j2,m2,j,m,ntype)
   implicit none

   integer,intent(in) :: j1,j2,j,m1,m2,m, ntype
   double precision ::  wigner3j
   integer :: t,tmin,tmax
   double precision :: factorial, triangle,coef,sum_coef

   coef = (-1.d0)**((j1-j2-m)/ntype) * dsqrt(factorial((j1+m1)/ntype)*factorial((j1-m1)/ntype)* &
                                                     factorial((j2+m2)/ntype)*factorial((j2-m2)/ntype)* &
                                                     factorial((j+m)/ntype)*factorial((j-m)/ntype))

   triangle = dsqrt(factorial((j1+j2-j)/ntype)*factorial((j1-j2+j)/ntype)* &
               factorial((-j1+j2+j)/ntype)/factorial((j1+j2+j)/ntype+1))

   sum_coef=0.d0
   tmin = max((j2-j-m1)/ntype,(j1+m2-j)/ntype,0)
   tmax = min((j1+j2-j)/ntype,(j1-m1)/ntype,(j2+m2)/ntype)
   do t=tmin,tmax
      sum_coef = sum_coef + (-1.d0)**t / (                   &
      factorial(t)*factorial((j1+j2-j)/ntype-t)*factorial((j1-m1)/ntype-t)* &
                   factorial((j2+m2)/ntype-t)*factorial((j-j2+m1)/ntype+t)*factorial((j-j1-m2)/ntype+t) )
   enddo

   wigner3j = coef * triangle * sum_coef

end function wigner3j


function factorial (n) result(res)

   use ml_in_ndm_module, only : rangml
   implicit none

   integer,intent(in) :: n
   double precision :: res
   integer :: i

   if ((n<0).and.(rangml==0)) then
      stop "Argument of factorial function should be positive"
   elseif (n==0) then
      res=1
   endif

   res=1.d0
   do i=2,n
      res = res*i
   enddo

end function factorial




subroutine compute_cg_vector(j_max,ntype)
use ml_in_ndm_module, only : rangml, cg_vector

implicit none
double precision, intent(in) :: j_max
integer, intent(in) ::  ntype
double precision :: cg

integer :: j, j1, j2, m, m1,m2, j0_max, j1_max, j2_max

if ((ntype==1) .and. dabs(int(j_max)-j_max).gt.0.1d0) then
      if (rangml==0) then
        write(*,*) 'ML: Warning: probably you intend to use some SO3 or bi-SO3 -  related descriptor'
        write(*,*) 'ML: Warning: in that case l_max can be only integer, now l_max', j_max
      end if
      stop
end if

j0_max=int(j_max*ntype)
j1_max=int(j_max*ntype)
j2_max=int(j_max*ntype)

if (allocated(cg_vector)) deallocate(cg_vector) ; allocate(cg_vector( 0:j1_max,-j1_max:j1_max, &
                                                                         0:j2_max,-j2_max:j2_max, &
                                                                         0:j0_max,-j0_max:j0_max ))
 cg_vector(:,:,:,:,:,:)=0.d0
do j1=0,j1_max
   do j2=0,j2_max
      do j =0,j0_max
         do m1=-j1, j1,ntype
            do m2=-j2, j2,ntype
               do m=-j, j,ntype
                   if (m-m1-m2/=0) cycle
                   if (j1+j2-j <0) cycle
                   if (j-j1+j2 <0) cycle
                   if (j-j2+j1 <0) cycle
                   if (ntype==2) then
                         if (mod(j+j1+j2,ntype)==1) cycle
                   end if
                   cg_vector(j1,m1,j2,m2,j,m) =  cg(j1,m1,j2,m2,j,m, ntype)
                   !debug write(12,'(6i5, e20.10)') j1,j2,j,  m1,m2, m, cg(j1,m1,j2,m2,j,m, ntype)
               end do
            end do
         end do
      end do
   end do
end do
!stop 'test cg in compute cg vector'
return
end subroutine  compute_cg_vector
