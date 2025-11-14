!******************************************************************************!
!                              bzlib.f90
!------------------------------------------------------------------------------!
!> Fortran binding for functions in the `bzlib.h` header from the bzip2 file
!> compression library.
!
!  Written by Paul Fossati, <paul.fossati@cea.fr>
!  Copyright (c) 2024 CEA
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
module bzlib
    use iso_c_binding
    
    implicit none(external, type)
    public

    integer(c_int), parameter :: BZ_RUN              = 0
    integer(c_int), parameter :: BZ_FLUSH            = 1
    integer(c_int), parameter :: BZ_FINISH           = 2

    integer(c_int), parameter :: BZ_OK               = 0
    integer(c_int), parameter :: BZ_RUN_OK           = 1
    integer(c_int), parameter :: BZ_FLUSH_OK         = 2
    integer(c_int), parameter :: BZ_FINISH_OK        = 3
    integer(c_int), parameter :: BZ_STREAM_END       = 4
    integer(c_int), parameter :: BZ_SEQUENCE_ERROR   = -1
    integer(c_int), parameter :: BZ_PARAM_ERROR      = -2
    integer(c_int), parameter :: BZ_MEM_ERROR        = -3
    integer(c_int), parameter :: BZ_DATA_ERROR       = -4
    integer(c_int), parameter :: BZ_DATA_ERROR_MAGIC = -5
    integer(c_int), parameter :: BZ_IO_ERROR         = -6
    integer(c_int), parameter :: BZ_UNEXPECTED_EOF   = -7
    integer(c_int), parameter :: BZ_OUTBUFF_FULL     = -8
    integer(c_int), parameter :: BZ_CONFIG_ERROR     = -9

    integer(c_int), parameter :: BZ_MAX_UNUSED = 5000

    ! High-level procedures
    interface

        function BZ2_bzReadOpen(bzerror, f, verbosity, small, unused, nUnused) bind(C, name="BZ2_bzReadOpen")
            use iso_c_binding
            type(c_ptr) :: BZ2_bzReadOpen
            integer(c_int), intent(out) :: bzerror
            type (c_ptr), value :: f
            integer(c_int), value :: verbosity
            integer(c_int), value :: small
            type(c_ptr), value :: unused
            integer(c_int), value :: nUnused
        end function

        function BZ2_bzRead(bzerror, b, buf, len) bind(C, name="BZ2_bzRead")
            use iso_c_binding
            integer(c_int) :: BZ2_bzRead
            integer(c_int) :: bzerror
            type(c_ptr), value :: b
            character(kind=c_char), dimension(*) :: buf
            integer(c_int), value :: len
        end function

        subroutine BZ2_bzReadGetUnused(bzerror, b, unused, nUnused) bind(C, name="BZ2_bzReadGetUnused")
            use iso_c_binding
            integer(c_int), intent(out) :: bzerror
            type(c_ptr), value :: b
            type(c_ptr), value :: unused
            integer(c_int), intent(out) :: nUnused
        end subroutine

        subroutine BZ2_bzReadClose(bzerror, b) bind(c, name="BZ2_bzReadClose")
            use iso_c_binding
            integer(c_int), intent(inout) :: bzerror
            type(c_ptr), value :: b
        end subroutine

        function BZ2_bzWriteOpen(bzerror, f, blockSize100k, verbosity, workFactor) bind(C, name="BZ2_bzWriteOpen")
            use iso_c_binding
            type(c_ptr) :: BZ2_bzWriteOpen
            integer(c_int), intent(out) :: bzerror
            type (c_ptr), value :: f
            integer(c_int), value :: blockSize100k
            integer(c_int), value :: verbosity
            integer(c_int), value :: workFactor
        end function

        subroutine BZ2_bzWrite(bzerror, b, buf, len) bind(C, name="BZ2_bzWrite")
            use iso_c_binding
            integer(c_int) :: bzerror
            type(c_ptr), value :: b
            character(kind=c_char), dimension(*) :: buf
            integer(c_int), value :: len
        end subroutine

        subroutine BZ2_bzWriteClose(bzerror, b, abandon, nbytes_in, nbytes_out) bind(c, name="BZ2_bzWriteClose")
            use iso_c_binding
            integer(c_int), intent(out) :: bzerror
            type(c_ptr), value :: b
            integer(c_int), value :: abandon
            integer(c_int), intent(out) :: nbytes_in
            integer(c_int), intent(out) :: nbytes_out
        end subroutine
    end interface

    ! Utility functions
    interface
        function BZ2_bzBuffToBuffCompress(dest, destLen, source, sourceLen, blockSize100k, verbosity, workFactor) bind(C, name="BZ2_bzBuffToBuffCompress")
            use iso_c_binding
            integer(c_int) :: BZ2_bzBuffToBuffCompress
            character(1,c_char), dimension(*), intent(out) :: dest
            integer(c_int), intent(out) :: destLen !< This is an unsigned int in the C library; it could be negative
            character(1,c_char), dimension(*), intent(in) :: source
            integer(c_int), value :: sourceLen !< This is an unsigned int in the C library
            integer(c_int), value :: blockSize100k
            integer(c_int), value :: verbosity
            integer(c_int), value :: workFactor
        end function

        function BZ2_bzBuffToBuffDecompress(dest, destLen, source, sourceLen, small, verbosity) bind(C, name="BZ2_bzBuffToBuffDecompress")
            use iso_c_binding
            integer(c_int) :: BZ2_bzBuffToBuffDecompress
            character(1,c_char), dimension(*), intent(out) :: dest
            integer(c_int), intent(out) :: destLen !< This is an unsigned int in the C library; it could be negative
            character(1,c_char), dimension(*), intent(in) :: source
            integer(c_int), value :: sourceLen !< This is an unsigned int in the C library
            integer(c_int), value :: small
            integer(c_int), value :: verbosity
        end function
    end interface
    
    ! zlib compatibility functions
    interface
        function bzlibVersion() bind(c, name="BZ2_bzlibVersion")
            use iso_c_binding
            type(c_ptr) :: bzlibVersion
        end function
        
        function bzopen(path, mode) bind(C, name="BZ2_bzopen")
            use iso_c_binding
            type (c_ptr) :: bzopen
            character(kind=c_char), dimension(*) :: path, mode
        end function

        function bzdopen(fd, mode) bind(C, name="BZ2_bzdopen")
            use iso_c_binding
            type (c_ptr) :: bzdopen
            integer(c_int), value :: fd
            character(kind=c_char), dimension(*) :: mode
        end function

        function bzread(file, buf, len) bind(C, name="BZ2_bzread")
            use iso_c_binding
            integer(c_int) :: bzread
            type(c_ptr), value :: file
            character(kind=c_char), dimension(*) :: buf
            integer(c_int), value :: len
        end function

        function bzwrite(file, buf, len) bind(C, name="BZ2_bzwrite")
            use iso_c_binding
            integer(c_int) ::  bzwrite
            type(c_ptr), value :: file
            character(kind=c_char), dimension(*) :: buf
            integer(c_int), value :: len
        end function

        function bzflush(file) bind(C, name="BZ2_bzflush")
            use iso_c_binding
            integer(c_int) :: bzflush
            type (c_ptr), value :: file
        end function

        function bzclose(file) bind(C, name="BZ2_bzclose")
            use iso_c_binding
            integer(c_int) :: bzclose
            type (c_ptr), value :: file
        end function

        function BZ2_bzerror(b, errnum) bind(C, name="BZ2_bzerror")
            use iso_c_binding
            type(c_ptr) :: BZ2_bzerror
            type (c_ptr), value :: b
            integer(c_int) :: errnum
        end function
    end interface
        
end module
