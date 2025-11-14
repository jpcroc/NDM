!******************************************************************************!
!                         dk_posix_io module
!------------------------------------------------------------------------------!
!> Bindings to C functions from the `stdio.h` header to make the included
!> functions callable from Fortran programs.
!
!  Part of the dk_io library version 0.1.
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
module dk_posix_io
    use iso_c_binding

    use ext_character, only: c_f_string
    
    implicit none(external, type)
    public

! Interfaces to functions defined in stdio.h
    interface
        function feof(stream) bind(C, name="feof")
            use iso_c_binding
            integer(c_int) :: feof
            type(c_ptr), value :: stream
        end function

        function ferror(stream) bind(C, name="ferror")
            use iso_c_binding
            integer(c_int) :: ferror
            type(c_ptr), value :: stream
        end function

        function fflush(stream) bind(C, name="fflush")
            use iso_c_binding
            integer(c_int) :: fflush
            type(c_ptr), value :: stream
        end function
        
        function fopen(name, mode) bind(C, name="fopen")
            use iso_c_binding
            type(c_ptr) :: fopen
            character(1,c_char), dimension(*) :: name
            character(1,c_char), dimension(*) :: mode
        end function

        function fgets(str, n, stream) bind(C, name="fgets")
            use iso_c_binding
            type(c_ptr) :: fgets
            character(c_char), dimension(*), intent(out) :: str
            integer(c_int), intent(in), value :: n
            type(c_ptr), value :: stream
        end function

        function fputs(str, stream) bind(C, name="fputs")
            use iso_c_binding
            integer(c_int) :: fputs
            character(1,c_char), dimension(*), intent(in) :: str
            type(c_ptr), value :: stream
        end function

        function fread(buffer, blocSize, blocCount, stream) bind(C, name="fread")
            use iso_c_binding
            integer(c_size_t) :: fread
            character(kind=c_char), dimension(*) :: buffer
            integer(c_size_t), value :: blocSize
            integer(c_size_t), value :: blocCount
            type(c_ptr), value :: stream
        end function

        function fwrite(buffer, blocSize, blocCount, stream) bind(C, name="fwrite")
            use iso_c_binding
            integer(c_size_t) :: fwrite
            character(kind=c_char), dimension(*) :: buffer
            integer(c_size_t), value :: blocSize
            integer(c_size_t), value :: blocCount
            type(c_ptr), value :: stream
        end function

        function fclose(stream) bind(C, name="fclose")
            use iso_c_binding
            integer(c_int) :: fclose
            type(c_ptr), value :: stream
        end function

        function fseek(stream, offset, whence) bind(C, name="fseek")
            use iso_c_binding
            integer(c_int) :: fseek
            type(c_ptr), value :: stream
            integer(c_long), value :: offset
            integer(c_int), value :: whence
        end function
        
        function ftell(stream) bind(C, name="ftell")
            use iso_c_binding
            integer(c_long) :: ftell
            type(c_ptr), value :: stream
        end function
        
        function strerror(code) bind(C, name="strerror")
            use iso_c_binding
            type(c_ptr) :: strerror
            integer(c_int), value :: code
        end function

        subroutine rewind(stream) bind(C, name="rewind")
            use iso_c_binding
            type(c_ptr), value :: stream
        end subroutine

        function fileno(stream) bind(C, name="fileno")
            use iso_c_binding
            integer(c_int) :: fileno
            type(c_ptr), value :: stream
        end function
    end interface

! Interfaces to the wrappers in io_wrappers.c
    interface
        function errno() bind(C, name="errno_wrapper")
            use iso_c_binding
            integer(c_int) :: errno
        end function

        function seek_cur() bind(C, name="seek_cur")
            use iso_c_binding
            integer(c_int) :: seek_cur
        end function

        function seek_end() bind(C, name="seek_end")
            use iso_c_binding
            integer(c_int) :: seek_end
        end function

        function seek_set() bind(C, name="seek_set")
            use iso_c_binding
            integer(c_int) :: seek_set
        end function
    end interface
contains
!******************************************************************************!
!> Get the message corresponding to the latest error condition from the C
!> library.
!******************************************************************************!
    function error_message() result(message)
        character(:), allocatable :: message
!------
        type(c_ptr) :: str
        integer(c_int) :: errnum
!------
        errnum = errno()
        if (errnum == 0) then
            message = ""
        else
            str = strerror(errnum)
            call c_f_string(str, message)
        end if
    end function
end module
