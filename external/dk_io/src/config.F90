!******************************************************************************!
!                               config
!______________________________________________________________________________!
!> Compile-time constants for the dk_io library. This module is automatically
!> generated, do not modify it.
!
!  Written by Paul Fossati, <paul.fossati@cea.fr>
!  Copyright (c) 2023-2025 CEA
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
module dk_io_config
    implicit none (external, type)
    public
    save

    character(*), parameter:: DK_IO_STRING = PACKAGE_STRING
    character(*), parameter:: DK_IO_VERSION = PACKAGE_VERSION
    integer, parameter:: DK_IO_VERSION_MAJOR = 0
    integer, parameter:: DK_IO_VERSION_MINOR = 1
    character(*), parameter:: DK_IO_SOURCE_PATH = "/Users/paul/Codes/coredynamics/dk_core"
    character(*), parameter:: DK_IO_TEST_PATH = "/Users/paul/Codes/coredynamics/dk_core/tests"
contains
    subroutine dk_io() bind(C, name="dk_io")
    end subroutine
end module
