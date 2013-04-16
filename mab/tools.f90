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
