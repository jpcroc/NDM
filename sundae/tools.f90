! This file contains the following functions or subroutines:
!  OU_controle, genere_bruit2,control_angular_momenta, genrand, FACT and timestamp, init_random_seed
     
SUBROUTINE init_random_seed(k)
 INTEGER :: i, k, n, clock
 INTEGER, DIMENSION(:), ALLOCATABLE :: seed

 CALL RANDOM_SEED(size = n)
 ALLOCATE(seed(n))
        
 CALL SYSTEM_CLOCK(COUNT=clock)
         
 seed = k + 37 * (/ (i - 1, i = 1, n) /) !+ clock
 CALL RANDOM_SEED(PUT = seed)
        
 DEALLOCATE(seed)
END SUBROUTINE


Subroutine OU_control(p,q,a_sto,v0)
 USE T_kind_param_m, ONLY:  double 
 use gen_com_m, ONLY : im
real(double), dimension(3,im) :: p,q,a
real(double):: z1,z2,z3,z4,v0,a_sto,r
integer:: i,ic

do i=1,im 
 do ic=1,3
 r=2.d0
 do while (r.ge.1.d0)
   call random_number(z1)
   call random_number(z2)
   z3=2.0*z1-1.0
   z4=2.0*z2-1.0
   r=z3**2+z4**2
 enddo
 a(ic,i)=z3*sqrt(-2.0*log(r)/(r))
 enddo
enddo


a(1,1:im) = a(1,1:im) - sum(a(1,1:im))/dble(im)
a(2,1:im) = a(2,1:im) - sum(a(2,1:im))/dble(im)
a(3,1:im) = a(3,1:im) - sum(a(3,1:im))/dble(im)

call control_angular_momenta(a,q)


p(1,1:im) = p(1,1:im) - sum(p(1,1:im))/dble(im)
p(2,1:im) = p(2,1:im) - sum(p(2,1:im))/dble(im)
p(3,1:im) = p(3,1:im) - sum(p(3,1:im))/dble(im)
call control_angular_momenta(p,q)

!p = p*a_sto + a*v0*sqrt(2.d0-2.d0*a_sto**2)  pour van brutzel 
p = p*a_sto + a*v0*sqrt(1.d0-1.d0*a_sto**2)


end subroutine OU_control





subroutine genere_bruit2 (sig,gau)
   use T_kind_param_m, ONLY : double
   use gen_com_m,      ONLY : pi,im 
   implicit none

   real(double), dimension(6,im+1),intent(out)::gau
   real(double), dimension(3,im)  ,intent(in) ::sig
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



Subroutine control_angular_momenta(p,q)
use T_kind_param_m, ONLY:  double
use sundae_module, ONLY: m_i,rang
use gen_com_m,      ONLY : im 
implicit none
real(double), dimension(3,im) :: p,q
real(double) :: rx, ry, rz, r2x, r2y, r2z, r2
real(double) :: prx, pry, prz, px, py, pz, vrx, vry, vrz
real(double) :: omegax, omegay, omegaz
real(double), dimension(3) ::   scom
real(double), dimension(3,3) :: ainer, aineri
integer i,ic,ib

ainer = 0.d0
prx = 0.d0
pry = 0.d0
prz = 0.d0

!rang =0

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


real(8) Function ran3()
     use T_kind_param_m, ONLY : double 
     use random_art
      implicit none

      integer, parameter :: mbig  = 1000000000
      integer, parameter :: mseed = 161803398
      integer, parameter :: mz=0
      real(double), parameter :: fac=1./mbig

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
end 
