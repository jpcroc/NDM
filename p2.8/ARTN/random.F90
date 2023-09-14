!> @file
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> ART Module random
!!    Random number generator (from "Numerical Recipes").
!!    Returns a uniform random deviate between 0.0 and 1.0.
!!    Set idum to any negative value to initialize or
!!    reinitialize the sequence.

! $Revision$
! $Date: 2015-12-14 11:01:46 -0500 (Lun, 14 déc 2015) $
! $Id: random.f90 1534 2015-12-14 15:59:38Z mickael $

module random
  use gen_com_m,only:rang
  implicit none

  ! Shared variables
  integer :: idum, inext, inextp
  integer :: iff = 0
  integer, dimension(55) :: ma
contains



!> ART Function ran3
real(kind=8) function ran3(idec)

  implicit none
  integer,optional::idec
  integer::idc
  !Local variables
  integer :: mbig=1000000000
  integer :: mseed=161803398
  integer :: mz=0
  real(kind=8) :: fac

  integer :: i, mj, mk, ii, k
  !_______________________


  ! Any large mbig, and any smaller (but still large) mseed can be
  !  substituted for the above values.
  idc=0
  if (present(idec))then
     idc=idec
  else
     idc=rang*13
  end if
  mbig=mbig+idc
  mseed=mseed+idc
  fac=1./mbig
  if ( idum < 0 .or. iff == 0 ) then
     iff = 1
     mj = mseed - iabs(idum)
     mj = mod( mj, mbig )
     ma(55) = mj
     mk = 1

     do i = 1, 54
        ii = mod( 21*i, 55 )
        ma(ii) = mk
        mk = mj - mk
        if ( mk < mz ) mk = mk + mbig
        mj = ma(ii)
     end do

     do k = 1, 4
        do i=1,55
           ma(i) = ma(i) - ma(1 + mod( i+30, 55 ))
           if ( ma(i) < mz ) ma(i) = ma(i) + mbig
        end do
     end do

     inext = 0
     inextp = 31
     idum = 1

  end if

  inext = inext + 1
  if ( inext == 56 ) inext = 1
  inextp = inextp + 1
  if ( inextp == 56 ) inextp = 1
  mj = ma(inext) - ma(inextp)
  if ( mj < mz ) mj = mj + mbig
  ma(inext) = mj

  ran3 = mj * fac

END FUNCTION ran3
END MODULE random
