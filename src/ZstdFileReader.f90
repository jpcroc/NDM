!******************************************************************************!
!                        dk_zstdfilereader module
!------------------------------------------------------------------------------!
!> `[[ZstdFileReader(type)]]` derived type and type-bound procedures, to add
!> support for zstd-compressed files to `[[FileObject(type)]]`.
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
module dk_zstdfilereader
    use iso_c_binding, only: c_ptr, c_char, c_int, c_long, c_size_t, c_associated, &
        C_NULL_PTR, C_NULL_CHAR, c_loc

    use zstd

    use ext_character, only: c_f_string, bold, hex_string
    use dk_exception, only: FException, FExceptionDescription, EVENT_LEVEL_ERROR
    use dk_filereader, only: FileReader, BUFFER_LENGTH
    use dk_posix_io, only: feof, ferror, fopen, fclose, fwrite, errno, &
        strerror, fread, error_message, rewind_c => rewind

    use event_levels

    implicit none(external, type)
    private

    integer, parameter :: MODE_UNDEF = 0
    integer, parameter :: MODE_READ = 1
    integer, parameter :: MODE_WRITE = 2
    integer, parameter :: MODE_APPEND = 3

    type, extends(FileReader), public :: ZstdFileReader
        private
        character(:,c_char), pointer :: in_buffer_ => null()
        character(:,c_char), pointer :: out_buffer_ => null()
        type(zstd_inbuffer) :: in_buffer
        type(zstd_outbuffer) :: out_buffer
        type(c_ptr) :: dctx = C_NULL_PTR
        integer :: buffer_bottom = 1
        integer :: buffer_top = 0

        type(c_ptr) :: stream = C_NULL_PTR
        character(:), allocatable :: filename
        integer :: mode = MODE_UNDEF
    contains
        procedure, public :: open_file
        procedure, public :: close_file
        procedure, public :: read
        procedure, public :: write
        procedure, public :: rewind => rewind_

        final :: finalise
    end type
contains
!******************************************************************************!
!> Open a compressed zstandard file for reading
!******************************************************************************!
    subroutine open_file(this, filename, mode, stat)
        class(ZstdFileReader), intent(inout), target :: this
        character(*), intent(in) :: filename
        character(*), intent(in) :: mode
        type(FException), intent(out) :: stat
!------
        character(:), allocatable :: message
        type(c_ptr) :: c_message
        integer(c_long) :: ret
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
            this%stream = fopen(filename//C_NULL_CHAR, &
                "rb"//C_NULL_CHAR)

            ! Allocate the internal buffers
            if (associated(this%in_buffer_)) deallocate(this%in_buffer_)
            allocate(character(ZSTD_DStreamInSize(),c_char) :: this%in_buffer_)
            this%in_buffer = zstd_inbuffer(src=c_loc(this%in_buffer_), pos=len(this%in_buffer_), size=len(this%in_buffer_))
            if (associated(this%out_buffer_)) deallocate(this%out_buffer_)
            allocate(character(ZSTD_DStreamOutSize(),c_char) :: this%out_buffer_)
            this%out_buffer = zstd_outbuffer(dst=c_loc(this%out_buffer_), pos=0, size=len(this%out_buffer_))
            if (c_associated(this%dctx)) ret = ZSTD_freeDCtx(this%dctx)
            this%dctx = ZSTD_createDCtx()
            this%buffer_bottom = 1

            ! Error check
            if (.not.c_associated(this%stream)) then
                c_message = strerror(errno())
                call c_f_string(c_message, message)
                call stat%raise(FExceptionDescription(message, "when trying to open file " // bold(filename) // " for reading:"))
                return
            end if

            ! Set the mode
            this%mode = MODE_READ

            this%filename = filename

        ! Open the file for writing
        case ("w", "write")

            ! Open the file
            this%stream = fopen(filename//C_NULL_CHAR, &
                "wb"//C_NULL_CHAR)

            !TODO: check errors

            ! Error check
            if (.not.c_associated(this%stream)) then
                c_message = strerror(errno())
                call c_f_string(c_message, message)
                call stat%raise(FExceptionDescription(message, "when trying to open file " // bold(filename) // " for writing:"))
                return
            end if

            ! Set the mode
            this%mode = MODE_WRITE

            this%filename = filename

        ! Open the file for writing
        case ("a", "append")

            ! Open the file
            this%stream = fopen(filename//C_NULL_CHAR, &
                "ab"//C_NULL_CHAR)

            !TODO: check errors

            ! Error check
            if (.not.c_associated(this%stream)) then
                c_message = strerror(errno())
                call c_f_string(c_message, message)
                call stat%raise(FExceptionDescription(message, "when trying to open existing file " // bold(filename) // " for writing:"))
                return
            end if

            ! Set the mode
            this%mode = MODE_WRITE

            this%filename = filename

        ! Handle unknown modes
        case default
            call stat%raise(FExceptionDescription("Invalid mode", "when trying to open file " // bold(filename) // ":"))
            return
        end select
    end subroutine
!******************************************************************************!
!> Close a compressed zstandard file. 
!******************************************************************************!
    subroutine close_file(this)
        class(ZstdFileReader), intent(inout) :: this
!------
        integer :: err
!------
        if (c_associated(this%stream)) then
            err = fclose(this%stream)
            this%stream = C_NULL_PTR
            !TODO: check errors
        end if
     end subroutine
!******************************************************************************!
!> Read from a compressed zstandard file.
!>
!> The output buffer has a set size.
!******************************************************************************!
    subroutine read(this, buffer, n, stat)
        class(ZstdFileReader), intent(inout) :: this !< File reader object
        character(:,c_char), allocatable, intent(inout) :: buffer
        integer, intent(out) :: n !< Number of characters read
        type(FException), intent(inout) :: stat !< Error status
!------
        integer(c_size_t) :: size_read, ret
        type(c_ptr) :: c_message
        character(:,c_char), allocatable :: message
!------
        
        ! Allocate the buffer if needed
        if (.not.allocated(buffer)) allocate(character(len(this%in_buffer_),c_char) :: buffer)

        ! Stop here if the file is not open
        if (.not.c_associated(this%stream)) then
            call stat%raise("File not open.")
            return
        end if

        ! Do not bother if we are already at the end of file
        if (feof(this%stream) /= 0 .and. this%in_buffer%pos == 0) then
            call stat%raise("Trying to read past end of file.")
            return
        end if

        ! Stop here if the file was not opened for reading
        if (this%mode /= MODE_READ) then
            call stat%raise(FExceptionDescription("File not open for reading", "when trying to read from file " // bold(this%filename) // ":"))
            return
        end if

        this%out_buffer%pos = 0
        do while (this%out_buffer%pos < this%out_buffer%size)

            ! Check whether we need to refill the input buffer
            if (this%in_buffer%pos == this%in_buffer%size) then

                size_read = fread(this%in_buffer_, 1_c_size_t, len(this%in_buffer_, kind=c_size_t), this%stream)

                this%in_buffer%src = c_loc(this%in_buffer_)
                this%in_buffer%pos = 0
                this%in_buffer%size = size_read
            end if

            ret = zstd_decompressstream(this%dctx, this%out_buffer, this%in_buffer)
            if (ret /= 0) then
                if (zstd_iserror(ret)) then
                    c_message = ZSTD_getErrorName(ret)
                    call c_f_string(c_message, message)
                    call stat%raise(trim(message), EVENT_LEVEL_ERROR)
                end if
                exit
            end if
        end do

        buffer(:) = this%out_buffer_(:this%out_buffer%pos)
        n = int(this%out_buffer%pos, kind=kind(n))

        if (feof(this%stream) == 1 .and. stat==0) then
            call stat%raise("End of file.", EVENT_LEVEL_LOG)
        end if
    end subroutine
!******************************************************************************!
!> Write to a compressed zstandard file.
!******************************************************************************!
    subroutine write(this, buffer, n, stat)
        class(ZstdFileReader), intent(inout) :: this
        character(*,c_char), intent(in) :: buffer
        integer, intent(in) :: n
        type(FException), intent(out) :: stat
!------
        integer(c_size_t) :: size_written, compressed_size, actual_size
        character(:,c_char), allocatable, target :: compressed_buffer
!------

        ! Stop here if the file is not open
        if (.not.c_associated(this%stream)) then
            call stat%raise("File not open.")
            return
        end if

        if (this%mode /= MODE_WRITE) then
            call stat%raise(FExceptionDescription("File not open for writing", "when trying to write to file " // bold(this%filename) // ":"))
            return
        end if

        ! Get the size of the compressed buffer
        compressed_size = ZSTD_compressBound(len(buffer, kind=c_size_t))
        allocate(character(compressed_size) :: compressed_buffer)
        compressed_buffer(:) = ""

        ! Compress the buffer
        actual_size = ZSTD_compress(compressed_buffer, len(compressed_buffer, kind=c_size_t), buffer, len(buffer, kind=c_size_t), 2)

        ! Write the compressed buffer to the file
        size_written = fwrite(compressed_buffer, 1_c_size_t, actual_size, this%stream)

        ! Check what happened if the whole buffer could not be written
        if (size_written /= actual_size) then
            ! This line is ignore when testing, because it is difficult to make fwrite error reliably
            call stat%raise(error_message()) ! LCOV_EXCL_LINE
        end if
    end subroutine
!******************************************************************************!
!> Move the file position indicator to the beginning of the file.
!******************************************************************************!
    subroutine rewind_(this)
        class(ZstdFileReader), intent(inout) :: this
!------
        type(FException) :: stat
!------
        call close_file(this)
        if (this%mode == MODE_READ) then
            call open_file(this, this%filename, "r", stat)
        end if
        if (this%mode == MODE_WRITE) then
            call open_file(this, this%filename, "w", stat)
        end if
     end subroutine
!******************************************************************************!
!> Final subroutine, to make sure that the file is closed when the reader
!> disappears.
!******************************************************************************!
    elemental impure subroutine finalise(this)
        type(ZstdFileReader), intent(inout) :: this
!------
        integer(c_size_t) :: ret
!------
        call close_file(this)
        if (associated(this%in_buffer_)) deallocate(this%in_buffer_)
        if (associated(this%out_buffer_)) deallocate(this%out_buffer_)
        if (c_associated(this%dctx)) ret = ZSTD_freeDCtx(this%dctx)
    end subroutine
end module
