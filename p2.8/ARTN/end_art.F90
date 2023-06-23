!> @file
!! @author
!!    Copyright (C) Normand Mousseau, June 2001
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> ART end_art
!!    Procedure called to properly end an art simulation

! $Revision$
! $Date: 2015-12-14 11:01:46 -0500 (Lun, 14 déc 2015) $
! $Id: art.f90 1534 2015-12-14 15:59:38Z mickael $


!> ART end_art
module end_art_mod
  use read_parameters_mod,only:timestamp
  USE arret_ndm_mod,only:arret_ndm

contains
subroutine end_art( )

  use defs
  implicit none

  integer      :: ierror, ierr
  real(kind=8) :: t2  ! cputime


  if ( iproc == 0 ) then
     open( unit = FLOG, file = LOGFILE, status = 'unknown',&
          & action = 'write', position = 'append', iostat = ierror )
     write(FLOG,*) '********************** '
     write(FLOG,*)    '  A bientot !'
     write(FLOG,*) '********************** '
     call CPU_TIME( t2)
     write(FLOG,"(' CPU_TIME: ', f12.4, ' seg')") t2-t1
     call timestamp('End')
     close(FLOG)
  end if

  write(*,*) "Master exiting from end_art"
  call arret_ndm

END SUBROUTINE end_art
end module end_art_mod
