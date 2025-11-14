!******************************************************************************!
!                              com_parameters
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
 module dk_parameters
    use iso_fortran_env

    implicit none (external, type)
    public

!******************************************************************************!
!  Real numbers kinds
!******************************************************************************!
    integer, parameter :: sp = kind(0.0e0)                             !< Kind of single precision real variables
    integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp)) !< Kind of double precision real variables
    integer, parameter :: qp = max(selected_real_kind(2*precision(1.0_dp)),dp) !< Kind of quadruple precision real variables
    integer, parameter :: wp = dp                            !< Default kind for real variables
    integer, parameter :: e  = dp                            !< Alias of dp
    integer, parameter :: s  = sp                            !< Alias of sp
    integer, parameter :: q  = qp                            !< Alias of qp

!******************************************************************************!
!  Integer numbers kinds
!******************************************************************************!
    integer, parameter :: long = int64 !< Kind of long integers

!******************************************************************************!
!  String representation of boolean constants
!******************************************************************************!
    character(*), parameter :: YES   = "yes"                 !< True value for character strings
    character(*), parameter :: TRUE  = YES                   !< Alternative name for `YES`
    character(*), parameter :: NO    = "no"                  !< False value for character strings
    character(*), parameter :: FALSE = NO                    !< Alternative name for `NO`
    character(*), parameter :: UNDEF = "undef"               !< Undefined character string

!******************************************************************************!
!  Mathematical constants
!******************************************************************************!
    real(e), parameter :: ee           = 2.7182818284590452353602874713526625_e !< e
    real(e), parameter :: log2ee       = 1.4426950408889634073599246810018921_e !< log_2(e)
    real(e), parameter :: log10ee      = 0.4342944819032518276511289189166051_e !< log_10(e)
    real(e), parameter :: ln2          = 0.6931471805599453094172321214581766_e !< ln(2)
    real(e), parameter :: ln10         = 2.3025850929940456840179914546843642_e !< ln(10)
    real(e), parameter :: pi           = 3.1415926535897932384626433832795029_e !< Pi
    real(e), parameter :: twopi        = 2 * pi !< Two Pi
    real(e), parameter :: sqrtPi       = 1.7724538509055160272981674833411451_e !< Square root of Pi
    real(e), parameter :: piSq         = 9.8696044010893586188344909998761508_e !< Pi squared
    real(e), parameter :: halfPi       = 1.5707963267948966192313216916397514_e !< Pi / 2
    real(e), parameter :: quartPi      = 0.7853981633974483096156608458198757_e !< Pi / 4
    real(e), parameter :: invPi        = 0.3183098861837906715377675267450287_e !< One over Pi
    real(e), parameter :: twoInvPi     = 0.6366197723675813430755350534900574_e !< Two over Pi
    real(e), parameter :: twoInvSqrtPi = 1.1283791670955125738961589031215452_e !< Two over the square root of Pi
    real(e), parameter :: sqrt2        = 1.4142135623730950488016887242096981_e !< Square root of 2
    real(e), parameter :: invSqrt2     = 0.7071067811865475244008443621048490_e !< One over the square root of 2
    real(e), parameter :: sqrt3        = 1.7320508075688772935274463415058723_e !< Square root of 3
    real(e), parameter :: invSqrt3     = 0.5773502691896257645091487805019574_e !< One over the square root of 3

contains
end module
