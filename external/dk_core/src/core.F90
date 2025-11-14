!******************************************************************************!
!                                   core
!______________________________________________________________________________!
!> Basic parameters for the libcore library.
!
!  Part of the Core library version 2.0.
!  Written by Paul Fossati, <paul.fossati@gmail.com>
!  Copyright (c) 2009-2018 Paul Fossati
!
!------------------------------------------------------------------------------!
! Redistribution and use in source and binary forms, with or without           !
! modification, are permitted provided that the following conditions are met:  !
!                                                                              !
!     * Redistributions of source code must retain the above copyright notice, !
!       this list of conditions and the following disclaimer.                  !
!                                                                              !
!     * Redistributions in binary form must reproduce the above copyright      !
!       notice, this list of conditions and the following disclaimer in the    !
!       documentation and/or other materials provided with the distribution.   !
!                                                                              !
!     * The name of the author may not be used to endorse or promote products  !
!      derived from this software without specific prior written permission    !
!      from the author.                                                        !
!                                                                              !
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"  !
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    !
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE   !
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE     !
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR          !
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF         !
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS     !
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      !
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)      !
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE   !
! POSSIBILITY OF SUCH DAMAGE.                                                  !
!******************************************************************************!
module core
    use dk_parameters

    implicit none (external, type)
    public
    save

    character(*), parameter:: CORE_NAME = PACKAGE_NAME
    character(*), parameter:: CORE_STRING = PACKAGE_STRING
    character(*), parameter:: CORE_VERSION = PACKAGE_VERSION

!******************************************************************************!
!  Misc.
!******************************************************************************!
    logical :: testing = .false. !< If set to .true., failed assertions do not stop the program

!******************************************************************************!
!  Generic interface to the assertion subroutines
!******************************************************************************!
    interface assert
        module procedure assertExpression
    end interface
    public :: assert

!******************************************************************************!
!  Generic interface to the comparison functions
!******************************************************************************!
    interface operator(.approx.)
#ifdef HAS_REAL128
        module procedure &
            almostEqual_sp_sp, almostEqual_sp_dp, almostEqual_sp_qp, &
            almostEqual_dp_sp, almostEqual_dp_dp, almostEqual_dp_qp, &
            almostEqual_qp_sp, almostEqual_qp_dp, almostEqual_qp_qp
#else
        module procedure &
            almostEqual_sp_sp, almostEqual_sp_dp, &
            almostEqual_dp_sp, almostEqual_dp_dp
#endif
    end interface
    public :: operator(.approx.)

!******************************************************************************!
!  Generic interface to the comparison functions using custom tolerances.
!******************************************************************************!
    interface approx
        module procedure approx_sp_sp, approx_dp_dp
    end interface
    public :: approx

contains
!******************************************************************************!
!> Compare two single-precision real numbers using custom tolerances.
!******************************************************************************!
    elemental function approx_sp_sp(x, y, maxdiff, maxreldiff) result(approx)
        real(sp), intent(in) :: x
        real(sp), intent(in) :: y
        real(sp), intent(in) :: maxdiff
        real(sp), intent(in) :: maxreldiff
        logical :: approx
!------
        include 'almostEqual.inc'
    end function
!******************************************************************************!
!> Compare two double-precision real numbers using custom tolerances.
!******************************************************************************!
    elemental function approx_dp_dp(x, y, maxdiff, maxreldiff) result(approx)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y
        real(dp), intent(in) :: maxdiff
        real(dp), intent(in) :: maxreldiff
        logical :: approx
!------
        include 'almostEqual.inc'
    end function
!******************************************************************************!
!> Compare two single-precision real numbers
!******************************************************************************!
    elemental function almostEqual_sp_sp(x, y) result(approx)
        real(sp), intent(in) :: x
        real(sp), intent(in) :: y
        logical :: approx
!------
        include 'almostEqual_wrapper.inc'
    end function
!******************************************************************************!
!> Compare single- and double-precision real numbers
!******************************************************************************!
    elemental function almostEqual_sp_dp(x, y) result(approx)
        real(sp), intent(in) :: x
        real(dp), intent(in) :: y
        logical :: approx
!------
        include 'almostEqual_wrapper.inc'
    end function
!******************************************************************************!
!> Compare two double-precision real numbers
!******************************************************************************!
    elemental function almostEqual_dp_dp(x, y) result(approx)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y
        logical :: approx
!------
        include 'almostEqual_wrapper.inc'
    end function
!******************************************************************************!
!> Compare double- and single-precision real numbers
!******************************************************************************!
    elemental function almostEqual_dp_sp(x, y) result(approx)
        real(dp), intent(in) :: x
        real(sp), intent(in) :: y
        logical :: approx
!------
        include 'almostEqual_wrapper.inc'
    end function
#ifdef HAS_REAL128
!******************************************************************************!
!> Compare single- and quadruple-precision real numbers
!******************************************************************************!
    elemental function almostEqual_sp_qp(x, y) result(approx)
        real(sp), intent(in) :: x
        real(qp), intent(in) :: y
        logical :: approx
!------
        include 'almostEqual_wrapper.inc'
    end function
!******************************************************************************!
!> Compare double- and quadruple-precision real numbers
!******************************************************************************!
    elemental function almostEqual_dp_qp(x, y) result(approx)
        real(dp), intent(in) :: x
        real(qp), intent(in) :: y
        logical :: approx
!------
        include 'almostEqual_wrapper.inc'
    end function
!******************************************************************************!
!> Compare two quadruple-precision real numbers
!******************************************************************************!
    elemental function almostEqual_qp_qp(x, y) result(approx)
        real(qp), intent(in) :: x
        real(qp), intent(in) :: y
        logical :: approx
!------
        include 'almostEqual_wrapper.inc'
    end function
!******************************************************************************!
!> Compare quadruple- and single-precision real numbers
!******************************************************************************!
    elemental function almostEqual_qp_sp(x, y) result(approx)
        real(qp), intent(in) :: x
        real(sp), intent(in) :: y
        logical :: approx
!------
        include 'almostEqual_wrapper.inc'
    end function
!******************************************************************************!
!> Compare quadruple- and double-precision real numbers
!******************************************************************************!
    elemental function almostEqual_qp_dp(x, y) result(approx)
        real(qp), intent(in) :: x
        real(dp), intent(in) :: y
        logical :: approx
!------
        include 'almostEqual_wrapper.inc'
    end function
#endif
!******************************************************************************!
!> Stop the program if a condition is not true
!******************************************************************************!
    subroutine assertExpression(expression, message)
        logical, intent(in) :: expression !< Expression to test
        character(*), intent(in), optional :: message !< message to print if the expression is false
!------

        ! Do nothing if the expression is true
        if (expression) return

        ! Print a custom message if provided
        if (present(message)) then
            write(*,*) "Assertion failed: " // message
        else
            write(*,*) "Assertion failed"
        end if

        ! Do not stop during test runs
        if (.not.testing) then
            error stop
        end if
    end subroutine
end module
