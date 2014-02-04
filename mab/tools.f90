! This file contains the following functions or subroutines:
!  genere_bruit2, genrand, FACT and timestamp      
MODULE random_mab

  ! Random number generator (from "Numerical Recipes").
  ! Returns a uniform random deviate between 0.0 and 1.0.
  ! Set idum to any negative value to initialize or  
  ! reinitialize the sequence.                      

  ! Shared variables
  save
  integer :: idum, inext, inextp
  integer :: iff = 0
  integer, dimension(55) :: ma
end module random_mab

subroutine genere_bruit2 (sig,gau)
   use T_kind_param_m, ONLY : double
   use gen_com_m,      ONLY : pi,im 
   implicit none

   real(double), dimension(6,im+1)::gau
   real(double), dimension(3,im)::sig
   real(double) :: u1,u2,b1
   real(double):: deriv(6)
   integer   :: ic,i
 
   gau(:,:)=0   

!write(*,*) 'sig' ,sig(:,1)

   do i=1,im+1
      do ic=1,6
         call random_number(u1)
         call random_number(u2)
         b1=sqrt(-2.*log(u1))*cos(2.*pi*u2)
         gau(ic,i) = b1
      enddo
   enddo

 
!write(*,*) 'gau', gau(:,:)


   gau(1:3,1:im) = gau(1:3,1:im)*sig(1:3,1:im)
   gau(4:6,1:im) = gau(4:6,1:im)*sig(1:3,1:im)
!write(*,*) 'gau2', gau(1:3,:)

   do ic=1,6
      deriv(ic)    = sum(gau(ic,1:im))/dble(im)
      !write(*,*) deriv(ic)
      gau(ic,1:im) = gau(ic,1:im)-deriv(ic)
   enddo

  end subroutine genere_bruit2




subroutine genere_bruit(gau)! bruit gaussien

   use T_kind_param_m, ONLY : double
   use gen_com_m,      ONLY : pi,im 
   implicit none

   real(double), dimension(6,im+1), intent(out)::gau
   real(double) :: u1,u2,b1

   integer   :: ic,i!,iatom,i
   
   gau(:,:)=0 

     do i=1,im+1
      do ic=1,6
         call random_number(u1)
         call random_number(u2)
         b1=sqrt(-2.*log(u1))*cos(2.*pi*u2)
         gau(ic,i) = b1
      enddo
    enddo
  return      
 
end subroutine genere_bruit



subroutine genere_bruit_one_value(value)! bruit gaussien

   use T_kind_param_m, ONLY : double
   use gen_com_m,      ONLY : pi,im 
   implicit none

   real(double),intent(out) :: value
   real(double) :: u1,u2,b1
   
   value=0 

         call random_number(u1)
         call random_number(u2)
         b1=sqrt(-2.*log(u1))*cos(2.*pi*u2)
         value = b1
   return
 
end subroutine genere_bruit_one_value



real(8) function ran3()
     use T_kind_param_m, ONLY : double 
     use random_mab
      implicit none

      integer, parameter :: mbig  = 1000000000
      integer, parameter :: mseed = 161803398
      integer, parameter :: mz=0
      real(double), parameter :: fac=1.d0/mbig

      integer :: i,mj, mk, ii, k

      ! Any large mbig, and any smaller (but still large) mseed can be
      !  substituted for the above values.

      if(idum.lt.0.or.iff.eq.0)then
           iff=1
           mj=mseed-iabs(idum)
           mj=mod(mj,mbig)
           ma(55)=mj
           mk=1
           do i=1,54
             ii=mod(21*i,55)
             ma(ii)=mk
             mk=mj-mk
             if(mk.lt.mz)mk=mk+mbig
             mj=ma(ii)
           enddo
           do k=1,4
             do i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.mz)ma(i)=ma(i)+mbig
             enddo
           enddo
           inext=0
           inextp=31
           idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
end function ran3




function genrand()
  use T_kind_param_m, ONLY : double 
 real (double) :: genrand 
  real (double) :: x
  call random_number(x)
  if (x.eq.1.0) x=0.99999999999
  genrand=x

end function genrand





!double precision 
        FUNCTION FACT(N)

	INTEGER I, N
        real*8 FACT
	FACT = 1.d0
	IF(N.LT.2)RETURN
	DO 100 I = 2, N
	   FACT = FACT*I
 100	CONTINUE
	RETURN
	END


! Fermi-Dirac function
  function FerDir (x,R,delta)
  use T_kind_param_m, ONLY : double
  real(double) :: FerDir  
  real(double) :: R, delta
  real(double) :: x,y,x1
  
  x1=(x-R)/delta
  y=1.d0/(1.d0+dexp(x1))

  FerDir=y
return
end

! derivative of the function 1-FD(x) which is:  d(1-FD(x))/dx=-d(FD(x))/dx, where FD is the Fermi-Dirac function
  function dFerDir(x,R,delta)
  use T_kind_param_m, ONLY : double
  real(double) :: dFerDir 
  real(double) :: R, delta
  real(double) :: x,y,x1
  
  x1=(x-R)/delta
  if (dabs(x)>=R) then
    y=0.25d0
   else 
    y=dexp(x1)/(1.d0+dexp(x1))**2
  end if 

  dFerDir=-y/delta
return
end

function Fermi_manuel(x,a,b_min,b_max)
 use T_kind_param_m, ONLY : double
real(double):: Fermi_manuel
real(double):: x,b_min,b_max,a
Fermi_manuel=1.d0/(1.d0+dexp(-a*(x-b_min)))+1.d0/(1.d0+dexp(-a*(b_max-x)))-1.d0

end function

Function gaussien_pdf(x,mu,sigma_gaussien)
use T_kind_param_m, ONLY : double
use gen_com_m, ONLY: pi
real(double):: gaussien_pdf,mu,sigma_gaussien,x,var
var=sigma_gaussien**2
gaussien_pdf=1.d0/dsqrt(2.d0*pi*var)*dexp(-(x-mu)**2/(2.d0*var))
end function



subroutine test_vacancy_position

 USE T_kind_param_m, ONLY:  double
 USE gen_com_m, ONLY: im,imm
 USE tab_imm_m
 USE mab_in_ndm_module, ONLY: normxlac, xlacf,xbarini,xbar,test_end,rtestlac
 implicit none
 real(double) :: rtemp

  test_end=.false.
  rtemp=dsqrt(DOT_PRODUCT(xp(:,7)-xlacf(:)-xbar(:)+xbarini(:),xp(:,7)-xlacf(:)-xbar(:)+xbarini(:)))

  !write(6,*) xp(:,7)
  !write(6,*)  rtemp, rtestlac

  if (rtemp <= rtestlac) test_end=.true.




end subroutine test_vacancy_position

 
 subroutine timestamp ( )

!***************************************************************************    **80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!



  implicit none
  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestamp 


Subroutine control_angular_momenta(p,q)

USE T_kind_param_m, ONLY:  double
USE mab_in_ndm_module, only: m_i  
use gen_com_m, only:im 


implicit none
real(double) :: rx, ry, rz, r2x, r2y, r2z, r2
real(double) :: prx, pry, prz, px, py, pz, vrx, vry, vrz
real(double) :: omegax, omegay, omegaz
real(double), dimension(3) ::   scom
real(double), dimension(3,3) :: ainer, aineri
real(double), dimension(3,im) :: p,q
integer ::         i,ic,ib,rang

ainer = 0.d0
prx = 0.d0
pry = 0.d0
prz = 0.d0


! provisoire: doit tenir compte des masses

do ic=1,3
 scom(ic) = sum(q(ic,1:im))/dble(im)
enddo

!write(*,*) 'p : ', p(1:3,1:im) 

do i = 1, im
                           
 rx = q(1,i)-scom(1)
 ry = q(2,i)-scom(2)
 rz = q(3,i)-scom(3)
             
 r2x = rx*rx
 r2y = ry*ry
 r2z = rz*rz
 r2  = r2x+r2y+r2z
 
 ainer(1,1) = ainer(1,1)+m_i(1,i)*(r2-r2x)
 ainer(2,2) = ainer(2,2)+m_i(1,i)*(r2-r2y)
 ainer(3,3) = ainer(3,3)+m_i(1,i)*(r2-r2z)
 ainer(2,3) = ainer(2,3)-m_i(1,i)*ry*rz
 ainer(3,1) = ainer(3,1)-m_i(1,i)*rz*rx
 ainer(1,2) = ainer(1,2)-m_i(1,i)*rx*ry
 px  = p(1,i)
 py  = p(2,i)
 pz  = p(3,i)
 prx = prx+ry*pz-rz*py
 pry = pry+rz*px-rx*pz
 prz = prz+rx*py-ry*px
enddo

ainer(3,2) = ainer(2,3)
ainer(1,3) = ainer(3,1)
ainer(2,1) = ainer(1,2)

rang=1
if (rang==0) then
       write(6,*) 
       write(6,997) (ainer(1,ib),ib=1,3),prx
       write(6,997) (ainer(2,ib),ib=1,3),pry
       write(6,997) (ainer(3,ib),ib=1,3),prz
endif
997    format('Inertia/anglm = ',3e15.6,4x,e14.6e3)

     !     calculate  angular velocity

call matinv(ainer,aineri)
omegax = aineri(1,1)*prx+aineri(1,2)*pry+aineri(1,3)*prz
omegay = aineri(2,1)*prx+aineri(2,2)*pry+aineri(2,3)*prz
omegaz = aineri(3,1)*prx+aineri(3,2)*pry+aineri(3,3)*prz

           !         shift velocities to make the angular momentum zero
do i = 1, im
 rx = q(1,i)-scom(1)
 ry = q(2,i)-scom(2)
 rz = q(3,i)-scom(3)
 vrx = omegay*rz-omegaz*ry
 vry = omegaz*rx-omegax*rz
 vrz = omegax*ry-omegay*rx

 p(1,i) = p(1,i)-vrx*m_i(1,i)
 p(2,i) = p(2,i)-vry*m_i(1,i)
 p(3,i) = p(3,i)-vrz*m_i(1,i)
enddo

end subroutine control_angular_momenta


subroutine init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
end subroutine



