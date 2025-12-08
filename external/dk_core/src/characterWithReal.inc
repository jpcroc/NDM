!******************************************************************************!
!                             characterWithReal
!______________________________________________________________________________!
!> Content of a subroutine to create character variables representing a
!> real value. This file is included in the ext_character module.
!
!  Written by Paul Fossati, <paul.fossati@gmail.com>
!  Copyright (c) 2009-2017 Paul Fossati
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
        integer :: stat, n, m, total, i
        character(:), allocatable :: format, tmp
        real(kind(number)) :: infinity
        integer :: decimalPlaces

        infinity = huge(number)

        ! Special cases
        if (number == 0) then
            this = "0"
            return
        else if (number >= infinity) then
            this = "Infinity"
            return
        else if (number <= -infinity) then
            this = "-Infinity"
            return
        else if (number /= number) then
            this = "NaN"
            return
        end if

        total = 0

        ! Get the number of positions on the right side of the decimal point
        decimalPlaces = precision(number)

        ! Calculate the number of positions on the left side of the decimal point
        n = floor(log10(abs(number)))
        if (n < 0) then
            n = n - 1
        else
            n = n + 1
        end if

        ! Simplify the number if it is an integer
        if (number < huge(1) .and. number > -huge(1)) then ! Avoid overflow during conversion
            if (nint(number) == number) then
                ! the variable is an integer
                this = nint(number)
                return
            end if
        end if

        if ((n < -maxexponent(number)) .or. (n > maxexponent(number))) then
            ! If the exponent is larger than the maximum exponent for this number kind, then there is a problem somewhere (most probably an optimising compiler)
            this = "NaN"

        else if (n < -2 .or. n > 2 ) then
            ! Scientific notation for small and large numbers
            m = ceiling(log10(abs(n*1.0_e)))
            total = 4 + decimalPlaces + m ! first digit + dot + 6 digits + E + sign + exponent
            if (number < 0) total = total + 1 ! Add one position for the - sign
            allocate(character(13) :: format)
            write(format,'(a,i0,a,i0,a,i0,a)') "(ES", total, ".", decimalPlaces,"E", m, ")"

            ! Print the number
            allocate(character(total) :: this)
            write(this,format,iostat=stat) number
            if (stat /= 0) this = "internal i/o error"

        else
            ! Standard case
            total = abs(n) + 1 + decimalPlaces
            if (number < 0) total = total + 1 ! Add one position for the - sign
            allocate(character(10) :: format)
            write(format,'(a,i0,a,i0,a)') "(F", total, ".", decimalPlaces, ")"

            ! Print the number
            allocate(character(total) :: tmp)
            write(tmp,format,iostat=stat) number
            if (stat /= 0) tmp = "internal i/o error"

            ! Try to remove trailing 0s
            ! Also remove leading spaces which can happen when NaN/=NaN is .false. because of optimisation flags
            i = index(tmp, "0", .true.)
            do while (i == len_trim(tmp))
                tmp(i:i) = ""
                i = index(tmp, "0", .true.)
            end do
            this = trim(adjustl(tmp))
        end if
