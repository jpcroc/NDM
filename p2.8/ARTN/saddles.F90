!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> ART Module saddles

! $Revision$
! $Date: 2015-12-14 11:01:46 -0500 (Lun, 14 déc 2015) $
! $Id: find_saddle.f90 1534 2015-12-14 15:59:38Z mickael $


module saddles

  implicit none

  integer :: natom_displaced    ! # of local atoms displaced
  integer :: KTER_MIN
  integer :: MAXKTER
  integer :: MAXIPERP, MAXKPERP,SMOOTH_DIR_CHANGE
  integer, dimension(:), allocatable  :: atom_displaced  ! Id of local atoms displaced
  real(kind=8), dimension(:), allocatable, target :: dr

  real(kind=8) :: INITSTEPSIZE
  real(kind=8) :: LOCAL_CUTOFF
  real(kind=8) :: INCREMENT,BASIN_FACTOR
  real(kind=8) :: FTHRESHOLD
  real(kind=8) :: EIGEN_THRESH
  real(kind=8) :: EXITTHRESH
  
  character(len=20) :: TYPE_EVENTS

  real(kind=8) :: delta_e  ! current_energy - ref_energy
  integer :: m_perp   ! Accepted perpendicular iterations.

  integer :: try      ! Total # of iterationes in a given hyperplane
  real(kind=8) :: ftot     ! Norm of the total force...
  real(kind=8) :: fpar     ! Parallel projection of force
  real(kind=8) :: fperp    ! Norm of the perpendicular force
  real(kind=8) :: delr     ! Distance between two configurations.
  integer :: npart    ! number of particles having moved by more
                      ! than a THRESHOLD
  logical :: switchDIIS
  logical :: end_activation
  integer :: nsteps_after_eigen_min
  real(kind=8) :: eigen_min

  !____DEV
  real(kind=8) :: coord_length
  integer      :: coord_number
  logical      :: cw_try_again ! for clean_wf

  ! Variables for GUESS_DIRECTION
  character(len=20)                       :: GUESSFILE
  real(kind=8), dimension(:), allocatable :: g_pos    ! Working positions of the atoms at the presumed direction
  real(kind=8)                            :: guess_noise

END MODULE saddles
