!> @file
!! @author
!!    Copyright (C) Normand Mousseau, June 2001
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> @brief ART module saddle_push
!! @detail
!!    Push the system over the saddle point


module saddle_push
  use calcforce_mod,only:calcforce
  use defs
  use boundary_cond_mod,only:boundary_cond
  implicit none

  integer, parameter, private :: push_messages = 4
  integer           , private :: push_log
  integer                     :: push_error
  real(dp)          , private :: push
  real(dp)          , private :: a1, b1

contains

  !> @brief
  !!    The system will be pushed towards a new basin. The displacement is
  !!    made along the direction of negative curvature, and by default away
  !!    from the initial minimum (invert = 1).
  !! @detail
  !!    Push the system over the saddle point:
  !!       pos = pos_saddle + push * negative_eigenvec
  !!    At second order, the change of energy dE_2 is given by
  !!       dE_2 := -push*force(sad) \dot negative_eigenvec(sad) +
  !!                1/2*push^2*negative_eigenval(sad)
  !!    push is decreased until following conditions are met simultaneously:
  !!       a- error := abs[ (E(pos)-E(saddle)-dE_2) / dE_2 ] < tolerance
  !!       b- E(pos)-E(saddle) < 0
  !!    push_log = 1 : Energy could not be decreased. Warn user and come back
  !!                   to the initial minimum.
  !!    push_log = 2 : Tolerance cannot be reached. Warn user but push using
  !!                   the push factor that minimized the error
  !!                   while decreasing energy.
  !!    push_log = 3 : Both "a" and "b" conditions were fulfilled.
  !!    push_log > 4 : 4 is added to push_log if negative_eigenvec is almost
  !!                   perpendicular to pos_sad_-pos_min_. Warn user.
  subroutine push_at_saddle(pos_, pos_sad_, pos_min_, force_, &
       projection_, eigenvalue_, E_sad, invert)


    real(dp), dimension(:), intent(inout) :: pos_
    real(dp), dimension(:), intent(in)    :: pos_sad_
    real(dp), dimension(:), intent(in)    :: pos_min_
    real(dp), dimension(:), intent(in)    :: force_
    real(dp), dimension(:), intent(in)    :: projection_
    real(dp),               intent(in)    :: eigenvalue_
    real(dp),               intent(in)    :: E_sad
    integer, optional,      intent(in)    :: invert

    real(dp), dimension(:), allocatable :: dpos, tmp_pos, tmp_force
    real(dp)                            :: dpos_norm
    real(dp)                            :: push_min
    real(dp)                            :: err, err_min
    real(dp)                            :: E_push
    real(dp)                            :: dE_push, dE_2
    integer                             :: i


    ! dpos      = (  pos_sad_ - pos_min_  ) / dpos_norm
    ! dpos_norm = || pos_sad_ - pos_min_ ||

    allocate(dpos(VECSIZE))
!    call boundary_cond(dpos, pos_sad_, pos_min_)
    dpos_norm = sqrt(dot_product(dpos, dpos))
    dpos = dpos/dpos_norm

    a1 = dot_product(dpos,   projection_)
    if (present(invert)) a1 = a1*invert
    b1 = dot_product(force_, projection_)
    deallocate(dpos)

    allocate(tmp_pos  (VECSIZE))

    push_log = 1
    push     = PUSH_OVER
    push_min = push
    err_min  = 1d+10

    allocate(tmp_force(VECSIZE))
    do i = 1, 15
       tmp_pos = pos_sad_ + sign(push,a1) * projection_
       call calcforce(NATOMS, tmp_pos, box, tmp_force, E_push, evalf_number)
       dE_2     = push*(-b1 + 0.5*push*eigenvalue_)
       dE_push  = E_push - E_sad
       err      = abs(dE_push/dE_2 - 1)
       if (err < err_min .and. dE_push < 0) then
          err_min  = err
          push_min = push
          push_log = 2
       end if
       if (err < PUSH_OVER_RELATIVE_ERROR .and. dE_push < 0) then
          push_log = 3
          exit
       else
          push = push * 0.5
       end if
    end do

    select case(push_log)
    case(1)
       push = PUSH_OVER
       pos_ = pos_sad_ + sign(PUSH_OVER,a1) * dpos_norm * projection_
    case(2)
       push = push_min
       pos_ = pos_sad_ + sign(push_min,a1)   * projection_       
!       pos_ = pos_sad_ + sign(push_min,a1) * dpos_norm   * projection_
    case(3)
       pos_ = tmp_pos
!       pos_ = pos_sad_ + sign(PUSH_OVER,a1) * dpos_norm * projection_
    end select
    
    deallocate(tmp_pos)

    if ( abs(a1) < 0.1d0 ) push_log = push_log + push_messages


  end subroutine push_at_saddle


  !> @brief Notify user about the outcome of the push
  !! @detail
  !!    If unitid is present, assume the associated file is already opened
  function saddle_push_log(unitid) result(ierror)

    integer, optional, intent(in) :: unitid

    integer :: ierror
    integer :: logunitid

    if (present(unitid)) then
       logunitid = unitid
    else
       logunitid = FLOG
    end if

    if (logunitid == FLOG) &
         open(unit = FLOG, file = LOGFILE, status = 'unknown', &
         action = 'write', position = 'append', iostat = ierror)

    if (push_log >= push_messages) then
       write(unit6P,*) 'BART : WARNING -- negative_eigenvec is almost &
            perpendicular to pos_sad-pos_min, their dot product equals ', a1
       write(logunitid,*)     '# WARNING -- negative_eigenvec is almost &
            perpendicular to pos_sad-pos_min, their dot product equals ', a1
    end if

    push_log = mod(push_log, push_messages)
    select case(push_log)
    case (1)
       write(unit6P,*)    'BART: PUSH -- failure, push increases energy, &
            with push = ', push
       write(logunitid,*) '# PUSH -- failure, push increases energy, &
            with push = ', push
    case(2)
       write(unit6P,*)    'BART: PUSH -- half success, energy decreases without &
            reaching second order approximation, with push = ', push
       write(logunitid,*) '# PUSH -- half success, energy decreases without &
            reaching second order approximation, with push = ', push
    case(3)
       write(unit6P,*)    'BART: PUSH -- success, with push = ', push
       write(logunitid,*) '# PUSH -- success, with push = ', push
    end select

    if (logunitid == FLOG) close(FLOG)

    if (push_log == 1) then
       push_error = 1
    else
       push_error = 0
    end if

  end function saddle_push_log


end module saddle_push
