!******************************************************************************!
!                         dk_filereader module
!------------------------------------------------------------------------------!
!> `[[FileReader(type)]]` derived type and type-bound procedures. This is an a
!> bstract type serving as a template for objects implementing operations on
!> different file types.
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
module dk_filereader
    use iso_c_binding, only: c_char

    use dk_exception, only: FException

    implicit none(external, type)

    ! integer, parameter, public :: BUFFER_LENGTH = 16000
    ! integer, parameter, public :: BUFFER_LENGTH = 262154 ! 256kB buffer
    integer, parameter, public :: BUFFER_LENGTH = 2097152 ! 2MB buffer
    
    type, abstract, public :: FileReader
        private
    contains
        procedure(reader_open_file), deferred, pass(this), public :: open_file
        procedure(reader_close_file), deferred, pass(this), public :: close_file
        procedure(reader_write), deferred, pass(this), public :: write
        procedure(reader_read), deferred, pass(this), public :: read
        procedure(reader_rewind), deferred, pass(this), public :: rewind
    end type

    abstract interface
        subroutine reader_open_file(this, filename, mode, stat)
            import FileReader, FException
            class(FileReader), intent(inout), target :: this
            character(*), intent(in) :: filename
            character(*), intent(in) :: mode
            type(FException), intent(out) :: stat
        end subroutine

        subroutine reader_close_file(this)
            import FileReader, FException
            class(FileReader), intent(inout) :: this
        end subroutine

        subroutine reader_read(this, buffer, n, stat)
            import FileReader, FException, c_char
            class(FileReader), intent(inout) :: this
            character(:,c_char), allocatable, intent(inout) :: buffer
            integer, intent(out) :: n
            type(FException), intent(inout) :: stat
        end subroutine

        subroutine reader_write(this, buffer, n, stat)
            import FileReader, FException, c_char
            class(FileReader), intent(inout) :: this
            character(*,c_char), intent(in) :: buffer
            integer, intent(in) :: n
            type(FException), intent(out) :: stat
        end subroutine

        subroutine reader_rewind(this)
            import FileReader
            class(FileReader), intent(inout) :: this
        end subroutine
    end interface
end module
