!******************************************************************************!
!                         dk_bzip2filereader module
!------------------------------------------------------------------------------!
!> `[[Bzip2FileReader(type)]]` derived type and type-bound procedures, to add
!> support for bzip2-compressed files to `[[FileObject(type)]]`.
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
module dk_bzip2filereader
    use iso_c_binding, only: c_ptr, c_char, c_int, c_associated, c_loc, C_NULL_PTR, C_NULL_CHAR

    use bzlib, only: BZ2_bzReadGetUnused, BZ2_bzReadOpen, BZ2_bzReadClose, bz2_bzwrite, &
        BZ2_bzread, bz2_bzwriteopen, BZ2_bzWriteClose, BZ_CONFIG_ERROR, BZ_PARAM_ERROR, BZ_IO_ERROR, &
        BZ_MEM_ERROR, BZ_CONFIG_ERROR, BZ_SEQUENCE_ERROR, BZ_DATA_ERROR, BZ_DATA_ERROR_MAGIC, BZ_OK, &
        BZ_STREAM_END, BZ_SEQUENCE_ERROR, BZ_UNEXPECTED_EOF, BZ_MAX_UNUSED
    use ext_character, only: c_f_string, bold, operator(//)
    use dk_exception, only: FException, FExceptionDescription
    use dk_filereader, only: FileReader, BUFFER_LENGTH
    use dk_posix_io, only: feof, ferror, errno, strerror, error_message, rewind_c => rewind, fopen, fclose

    implicit none(external, type)
    private

    type, extends(FileReader), public :: Bzip2FileReader
        private
        type(c_ptr) :: file = C_NULL_PTR !< Pointer to the C file handle
        type(c_ptr) :: stream = C_NULL_PTR !< Pointer to the BZFILE data stream
        logical :: read_mode = .false.  !< Indicates whether the file has been opened for reading
        logical :: write_mode = .false. !< Indicates whether the file has been opened for writing
        character(:), allocatable :: filename
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
!> Open a bzip2-compressed file.
!******************************************************************************!
    subroutine open_file(this, filename, mode, stat)
        class(Bzip2FileReader), intent(inout), target :: this
        character(*), intent(in) :: filename
        character(*), intent(in) :: mode
        type(FException), intent(out) :: stat
!------
        character(:), allocatable :: message
        type(c_ptr) :: c_message
        integer(c_int) :: bzerror
!------
        this%filename = filename

        ! Make sure the file is not already in use
        if (c_associated(this%stream)) then
            call stat%raise(FExceptionDescription("File reader already in use", "when trying to open file " // bold(filename) // ":"))
            return
        end if

        select case(mode)

        ! Open the file for reading
        case ("r", "read")

            ! Call the C library function
            this%file = fopen(filename//C_NULL_CHAR, &
                "rb"//C_NULL_CHAR)

            ! Error handling
            if (.not.c_associated(this%file)) then
                c_message = strerror(errno())
                call c_f_string(c_message, message)
                call stat%raise(FExceptionDescription(message, "when trying to open file " // bold(filename) // " for reading:"))
                return
            end if

            ! Initialise the data structures for decompression
            this%stream =  BZ2_bzReadOpen(bzerror, this%file, &
                verbosity=0_c_int, &
                small=0_c_int, &
                unused=c_null_ptr, &
                nUnused=0_c_int)

            ! Error handling
            select case(bzerror)
            case (BZ_CONFIG_ERROR)
                call stat%raise("The libbz2 library was mis-compiled and is unusable")
                return
            case (BZ_PARAM_ERROR)
                call stat%raise("Parameter error")
                return
            case (BZ_IO_ERROR)
                call stat%raise("IO error")
                return
            case (BZ_MEM_ERROR)
                call stat%raise("Not enough memory")
                return
            end select

            ! Set the mode
            this%read_mode = .true.
            this%write_mode = .false.

        ! Open the file for writing
        case ("w", "write")

            ! Call the C library function
            this%file = fopen(filename//C_NULL_CHAR, &
                "wb"//C_NULL_CHAR)

            ! Error handling
            if (.not.c_associated(this%file)) then
                c_message = strerror(errno())
                call c_f_string(c_message, message)
                call stat%raise(FExceptionDescription(message, "when trying to open file " // bold(filename) // " for writing:"))
                return
            end if

            ! Initialise the data structures for compression
            this%stream = bz2_bzwriteopen(bzerror, this%file, 5_c_int, &
                verbosity=0_c_int, &
                workfactor=0_c_int)

            ! Error handling
            select case(bzerror)
            case (BZ_CONFIG_ERROR)
                call stat%raise("The libbz2 library was mis-compiled and is unusable")
                return
            case (BZ_PARAM_ERROR)
                call stat%raise("Parameter error")
                return
            case (BZ_IO_ERROR)
                call stat%raise("IO error")
                return
            case (BZ_MEM_ERROR)
                call stat%raise("Not enough memory")
                return
            end select

            ! Set the mode
            this%read_mode = .false.
            this%write_mode = .true.

        ! Open the file for appending
        case ("a", "append")

            ! Call the C library function
            this%file = fopen(filename//C_NULL_CHAR, &
                "ab"//C_NULL_CHAR)

            ! Error handling
            if (.not.c_associated(this%file)) then
                c_message = strerror(errno())
                call c_f_string(c_message, message)
                call stat%raise(FExceptionDescription(message, "when trying to open file " // bold(filename) // " for writing:"))
                return
            end if

            ! Initialise the data structures for compression
            this%stream = bz2_bzwriteopen(bzerror, this%file, 5_c_int, &
                verbosity=0_c_int, &
                workfactor=0_c_int)

            ! Error handling
            select case(bzerror)
            case (BZ_CONFIG_ERROR)
                call stat%raise("The libbz2 library was mis-compiled and is unusable")
                return
            case (BZ_PARAM_ERROR)
                call stat%raise("Parameter error")
                return
            case (BZ_IO_ERROR)
                call stat%raise("IO error")
                return
            case (BZ_MEM_ERROR)
                call stat%raise("Not enough memory")
                return
            end select

            ! Set the mode
            this%read_mode = .false.
            this%write_mode = .true.

        ! Handle unknown modes
        case default
            call stat%raise(FExceptionDescription("Invalid mode", "when trying to open file " // bold(filename) // ":"))
            return
        end select
    end subroutine
!******************************************************************************!
!> Close a bzip2-compressed file.
!******************************************************************************!
    subroutine close_file(this)
        class(Bzip2FileReader), intent(inout) :: this
!------
        integer(c_int) :: err, in, out
!------
        ! Close the bzlib stream
        if (c_associated(this%stream)) then
            if (this%read_mode) then
                 call BZ2_bzReadClose(err, this%stream)
            end if
            if (this%write_mode) then
                call BZ2_bzWriteClose(err, this%stream, 0_c_int, in, out)
            end if
            this%stream = C_NULL_PTR
        end if

        ! Close the file descriptor
        if (c_associated(this%file)) then
            err = fclose(this%file)
            this%file = C_NULL_PTR
        end if
     end subroutine
!******************************************************************************!
!> Read from a bzip2-compressed file.
!******************************************************************************!
    subroutine read(this, buffer, n, stat)
        class(Bzip2FileReader), intent(inout) :: this !< File reader object
        character(:,c_char), allocatable, intent(inout) :: buffer
        integer, intent(out) :: n !< Number of characters read
        type(FException), intent(inout) :: stat !< Error status
!------
        integer(c_int) :: bzerror, err, nUnused
        character(:,c_char), allocatable, target :: unused
        type(c_ptr), target :: unused_
!------
        ! Do not try to read an un-opened file
        if (.not. c_associated(this%stream)) then
            call stat%raise(FExceptionDescription("File not open.", &
                "when trying to read from file " // bold(this%filename) // ":"))
            return
        end if

        ! Allocate the buffer if needed
        if (.not.allocated(buffer)) allocate(character(BUFFER_LENGTH,c_char) :: buffer)

        ! Fill the buffer from the file
        n = BZ2_bzread(bzerror, this%stream, buffer, len(buffer, kind=c_int))

        ! Avoid reading past the end of file
        ! Are there cases where this%file could not be associated at this point?
        if (feof(this%file) /= 0 .and. bzerror == BZ_UNEXPECTED_EOF) then
            call stat%raise(FExceptionDescription("Trying to read past end of file.", &
                "when trying to read from file " // bold(this%filename) // ":"))
            return
        end if

        ! Make sure we can keep reading if there is another stream after the one that just ended
        if (bzerror == BZ_STREAM_END) then

            ! Get any compressed memory read after the end of stream
            allocate(character(BZ_MAX_UNUSED,c_char) :: unused)
            unused_ = c_loc(unused)
            call BZ2_bzReadGetUnused(err, this%stream, c_loc(unused_), nUnused)

            ! Close the stream
            call BZ2_bzReadClose(err, this%stream)

            ! Re-open the stream to pick up where the previous one left off
            this%stream =  BZ2_bzReadOpen(err, this%file, &
                verbosity=0_c_int, &
                small=0_c_int, &
                unused=unused_, &
                nUnused=nUnused)
        end if

        ! Handle errors
        select case(bzerror)
        case (BZ_PARAM_ERROR)
            ! This should never happen: buffer cannot be unallocated, so its length cannot be negative
            call stat%raise("Internal error in bzip2filereader%read: BZ_PARAM_ERROR") ! LCOV_EXCL_LINE
        case (BZ_SEQUENCE_ERROR)
            call stat%raise(FExceptionDescription("The stream was opened for compression", "when trying to read from file " // bold(this%filename) // ":"))
        case (BZ_IO_ERROR)
            call stat%raise(FExceptionDescription("Could not read from file", "when trying to read from file " // bold(this%filename) // ":"))
        case (BZ_UNEXPECTED_EOF)
            call stat%raise(FExceptionDescription("Unexpected end of file.", "when trying to read from file " // bold(this%filename) // ":"))
        case (BZ_DATA_ERROR)
            call stat%raise(FExceptionDescription("Data integrity error", "when trying to read from file " // bold(this%filename) // ":"))
        case (BZ_DATA_ERROR_MAGIC)
            call stat%raise(FExceptionDescription("The file does not contain valid bzip2-compressed data", "when trying to read from file " // bold(this%filename) // ":"))
        case (BZ_MEM_ERROR)
            call stat%raise(FExceptionDescription("Not enough memory available", "when trying to read from file " // bold(this%filename) // ":"))
        case (BZ_OK)
        case (BZ_STREAM_END)
        case default
            ! This should never happen: all situations should be handled by the other cases
            call stat%raise("Internal error in bzip2filereader%read: bzerror="//bzerror) ! LCOV_EXCL_LINE
        end select
     end subroutine
!******************************************************************************!
!> Write to a bzip2-compressed file.
!******************************************************************************!
    subroutine write(this, buffer, n, stat)
        class(Bzip2FileReader), intent(inout) :: this
        character(*,c_char), intent(in) :: buffer
        integer, intent(in) :: n
        type(FException), intent(out) :: stat
!------
        integer(c_int) :: bzerror
!------
        ! Write the buffer to the file
        call bz2_bzwrite(bzerror, this%stream, buffer, int(len(buffer), c_int))

        ! Handle errors
        select case(bzerror)
        case (BZ_PARAM_ERROR)
            if (.not.c_associated(this%stream)) then
                call stat%raise(FExceptionDescription("File not open.", "when trying to write to file " // bold(this%filename) // ":"))
            else
                ! This should never happen: buffer cannot be unallocated, so its length cannot be negative
                call stat%raise("Internal error in bzip2filereader%write: BZ_PARAM_ERROR") ! LCOV_EXCL_LINE
            end if
        case (BZ_SEQUENCE_ERROR)
            call stat%raise(FExceptionDescription("The stream was opened for decompression", "when trying to write to file " // bold(this%filename) // ":"))
        case (BZ_IO_ERROR)
            call stat%raise(FExceptionDescription("Could not write to file", "when trying to write to file " // bold(this%filename) // ":"))
        end select
    end subroutine
!******************************************************************************!
!> Move the file position indicator to the beginning of the file.
!******************************************************************************!
    subroutine rewind_(this)
        class(Bzip2FileReader), intent(inout) :: this
!------
        type(FException) :: stat_
!------
        ! There is no handy rewind() function in libbz2; we have to emulate it by closing the file and reopening it
        call close_file(this)
        if (this%read_mode) then
            call open_file(this, this%filename, "r", stat_)
        else if (this%write_mode) then
            call open_file(this, this%filename, "w", stat_)
        end if
    end subroutine
!******************************************************************************!
!> Final subroutine, to make sure that the file is closed when the reader
!> disappears.
!******************************************************************************!
    elemental impure subroutine finalise(this)
        type(Bzip2FileReader), intent(inout) :: this
!------
        call close_file(this)
    end subroutine
end module
