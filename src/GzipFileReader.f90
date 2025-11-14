!******************************************************************************!
!                        dk_gzipfilereader module
!------------------------------------------------------------------------------!
!> `[[GzipFileReader(type)]]` derived type and type-bound procedures, to add
!> support for gzip-compressed files to `[[FileObject(type)]]`.
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
module dk_gzipfilereader
    use iso_c_binding, only: c_ptr, c_char, c_int, c_associated, C_NULL_PTR, C_NULL_CHAR

    use zlib
    
    use ext_character, only: c_f_string, bold
    use dk_exception, only: FException, FExceptionDescription, EVENT_LEVEL_LOG
    use dk_filereader, only: FileReader, BUFFER_LENGTH
    use dk_posix_io, only: feof, ferror, errno, strerror, error_message

    implicit none(external, type)
    private

    type, extends(FileReader), public :: GzipFileReader
        private
        character(:), allocatable :: filename
        type(c_ptr) :: stream = C_NULL_PTR !< Pointer to the C stream
        logical :: read_mode = .false.  !< Indicates whether the file has been opened for reading
        logical :: write_mode = .false. !< Indicates whether the file has been opened for writing
    contains
        procedure, public :: open_file
        procedure, public :: close_file
        procedure, public :: read
        procedure, public :: write
        procedure, public :: rewind

        final :: finalise
    end type

contains
!******************************************************************************!
!> Open an gzip-compressed file
!******************************************************************************!
    subroutine open_file(this, filename, mode, stat)
        class(GzipFileReader), intent(inout), target :: this
        character(*), intent(in) :: filename
        character(*), intent(in) :: mode
        type(FException), intent(out) :: stat
!------
        character(:), allocatable :: message
        type(c_ptr) :: c_message
        integer :: i, err
!------

        ! Check whether the file is already open
        if (c_associated(this%stream)) then

            ! If yes, raise an error
            call stat%raise(FExceptionDescription("File reader already in use", "when trying to open file " // bold(filename) // ":"))
            return
        end if

        select case(mode)

        ! Open the file for reading
        case ("r", "read")

            ! Open the file
            this%stream = gzopen(filename//C_NULL_CHAR, &
                "rb"//C_NULL_CHAR)

            ! Set the stream's buffer size to match that of the file reader
            err = gzbuffer(this%stream, BUFFER_LENGTH)

            !TODO: check errors

            ! Error check
            if (.not.c_associated(this%stream)) then
                c_message = strerror(errno())
                call c_f_string(c_message, message)
                call stat%raise(FExceptionDescription(message, "when trying to open file " // bold(filename) // " for reading:"))
                return
            end if

            ! Set the mode
            this%read_mode = .true.
            this%write_mode = .false.

            this%filename = filename

        ! Open the file for writing
        case ("w", "write")

            ! Open the file
            this%stream = gzopen(filename//C_NULL_CHAR, &
                "wb"//C_NULL_CHAR)

            ! Set the stream's buffer size to match that of the file reader
            err = gzbuffer(this%stream, BUFFER_LENGTH)

            !TODO: check errors

            ! Error check
            if (.not.c_associated(this%stream)) then
                c_message = strerror(errno())
                call c_f_string(c_message, message)
                call stat%raise(FExceptionDescription(message, "when trying to open file " // bold(filename) // " for writing:"))
                return
            end if

            ! Set the mode
            this%read_mode = .false.
            this%write_mode = .true.

            this%filename = filename
            
        ! Open the file for writing
        case ("a", "append")

            ! Open the file
            this%stream = gzopen(filename//C_NULL_CHAR, &
                "ab"//C_NULL_CHAR)

            ! Set the stream's buffer size to match that of the file reader
            err = gzbuffer(this%stream, BUFFER_LENGTH)

            !TODO: check errors

            ! Error check
            if (.not.c_associated(this%stream)) then
                c_message = strerror(errno())
                call c_f_string(c_message, message)
                call stat%raise(FExceptionDescription(message, "when trying to open existing file " // bold(filename) // " for writing:"))
                return
            end if

            ! Set the mode
            this%read_mode = .false.
            this%write_mode = .true.

            this%filename = filename
            
        ! Handle unknown modes
        case default
            call stat%raise(FExceptionDescription("Invalid mode", "when trying to open file " // bold(filename) // ":"))
            return
        end select

        ! Set the buffer
        i = gzbuffer(this%stream, 4194304_c_int)
    end subroutine
!******************************************************************************!
!> Close a gzip-compressed file. 
!******************************************************************************!
    subroutine close_file(this)
        class(GzipFileReader), intent(inout) :: this
!------
        integer :: err
!------
        if (c_associated(this%stream)) then
            err = gzclose(this%stream)
            this%stream = C_NULL_PTR
            !TODO: check errors
        end if
     end subroutine
!******************************************************************************!
!> Read from a gzip-compressed file.
!******************************************************************************!
    subroutine read(this, buffer, n, stat)
        class(GzipFileReader), intent(inout) :: this !< File reader object
        character(:,c_char), allocatable, intent(inout) :: buffer
        integer, intent(out) :: n !< Number of characters read
        type(FException), intent(inout) :: stat !< Error status
!------
        integer(c_int) :: errnum
        type(c_ptr) :: c_message
        character(:), allocatable :: message
!------

        ! Allocate the buffer if needed
        if (.not.allocated(buffer)) allocate(character(BUFFER_LENGTH,c_char) :: buffer)

        ! Stop here if the file is not open
        if (.not.c_associated(this%stream)) then
            call stat%raise("File not open.")
            return
        end if

        ! Stop here if the file was not opened for reading
        if (.not. this%read_mode) then
            call stat%raise(FExceptionDescription("File not open for reading.", "when trying to read from file " // bold(this%filename) // ":"))
            return
        end if

        ! Stop here if we are already at the end of the file before reading
        if (gzeof(this%stream) /= 0) then
            call stat%raise(FExceptionDescription("Trying to read past end of file.", "when trying to read from file " // bold(this%filename) // ":"))
            return
        end if

        ! Fill the buffer from the file
        n = gzread(this%stream, buffer, len(buffer, c_int))
        
        ! Error handling
        c_message = gzerror(this%stream, errnum)
        call c_f_string(c_message, message)
        select case(errnum)
        case (Z_ERRNO)
            call stat%raise("Z_ERRNO")
        case (Z_DATA_ERROR)
            call stat%raise("Data integrity error")
        case (Z_STREAM_ERROR)
            call stat%raise("Z_STREAM_ERROR")
        case (Z_NEED_DICT)
            call stat%raise("Z_NEED_DICT")
        case (Z_MEM_ERROR)
            call stat%raise("Z_MEM_ERROR")
        case (Z_BUF_ERROR)
            call stat%raise("Unexpected end of file.")
        case (Z_OK)
            if (gzeof(this%stream) /= 0) call stat%raise("End of file.", EVENT_LEVEL_LOG)
        case default
            write(*,*) ":: ", message, errnum
        end select

     end subroutine
!******************************************************************************!
!> Write to a gzip-compressed file.
!******************************************************************************!
    subroutine write(this, buffer, n, stat)
        class(GzipFileReader), intent(inout) :: this
        character(*,c_char), intent(in) :: buffer
        integer, intent(in) :: n
        type(FException), intent(out) :: stat
!------
        integer :: size_written
        integer(c_int) :: errnum
        type(c_ptr) :: c_message
        character(:), allocatable :: message
!------

        ! Stop here if the file is not open
        if (.not.c_associated(this%stream)) then
            call stat%raise("File not open.")
            return
        end if

        ! Stop here if the file was not opened for writing
        if (.not. this%write_mode) then
            call stat%raise(FExceptionDescription("File not open for writing.", "when trying to write to file " // bold(this%filename) // ":"))
            return
        end if

        ! Write the buffer to the file
        size_written = gzwrite(this%stream, buffer, int(len(buffer), c_int))

        c_message = gzerror(this%stream, errnum)
        call c_f_string(c_message, message)
        select case(errnum)
        case (Z_ERRNO)
            call stat%raise("Z_ERRNO")
        case (Z_DATA_ERROR)
            call stat%raise("Data integrity error")
        case (Z_STREAM_ERROR)
            call stat%raise("Z_STREAM_ERROR")
        case (Z_NEED_DICT)
            call stat%raise("Z_NEED_DICT")
        case (Z_MEM_ERROR)
            call stat%raise("Z_MEM_ERROR")
        case (Z_BUF_ERROR)
            call stat%raise("Unexpected end of file.")
        case (Z_OK)
        case default
            write(*,*) ":: ", message, errnum
        end select


        ! Check what happened if the whole buffer could not be written
        if (size_written /= len(buffer)) then
            !if (gzerror(this%stream) /= 0) then
                ! Possible values: Z_ERRNO (check errno), Z_STREAM_ERROR, Z_BUF_ERROR, Z_MEM_ERROR
            !end if
            call stat%raise(error_message())
        end if
    end subroutine
!******************************************************************************!
!> Move the file position indicator to the beginning of the file.
!******************************************************************************!
    subroutine rewind(this)
        class(GzipFileReader), intent(inout) :: this
!------
        integer :: err
!------
        if (c_associated(this%stream)) then
            err = gzrewind(this%stream)
            !TODO: error handling
        end if
    end subroutine
!******************************************************************************!
!> Final subroutine, to make sure that the file is closed when the reader
!> disappears.
!******************************************************************************!
    elemental impure subroutine finalise(this)
        type(GzipFileReader), intent(inout) :: this
!------
        call close_file(this)
    end subroutine
end module
