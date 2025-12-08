!******************************************************************************!
!                         dk_textfilereader module
!------------------------------------------------------------------------------!
!> `[[TextFileReader(type)]]` derived type and type-bound procedures, to add
!> support for uncompressed files to `[[FileObject(type)]]`.
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
module dk_textfilereader
    use iso_c_binding, only: c_ptr, c_char, c_size_t, c_associated, C_NULL_PTR, C_NULL_CHAR, c_int, c_long

    use ext_character, only: c_f_string, bold
    use dk_exception, only: FException, FExceptionDescription
    use dk_filereader, only: FileReader, BUFFER_LENGTH
    use dk_posix_io, only: feof, ferror, fflush, fopen, fclose, fwrite, errno, &
        strerror, fgets, fread, error_message, rewind_c => rewind, fileno, fseek, ftell

    implicit none(external, type)
    private

    type, extends(FileReader), public :: TextFileReader
        private
        type(c_ptr) :: stream = C_NULL_PTR !< Pointer to the C stream
        logical :: read_mode = .false.  !< Indicates whether the file has been opened for reading
        logical :: write_mode = .false. !< Indicates whether the file has been opened for writing
    contains
        procedure, public :: open_file
        procedure, public :: close_file
        procedure, public :: read
        procedure, public :: write
        procedure, public :: rewind => rewind_
        
        final :: finalise
    end type

    interface
        function posix_fadvise(fd, offset, len, advice) bind(C, name="posix_fadvise")
            use iso_c_binding
            integer(c_int) :: posix_fadvise
            integer(c_int), value :: fd
            integer(c_size_t), value :: offset
            integer(c_size_t), value :: len
            integer(c_int), value :: advice
        end function
    end interface

    integer(c_int), parameter :: POSIX_FADV_NORMAL      = 0
    integer(c_int), parameter :: POSIX_FADV_SEQUENTIAL  = 1
    integer(c_int), parameter :: POSIX_FADV_RANDOM      = 2
    integer(c_int), parameter :: POSIX_FADV_WILLNEED    = 3
    integer(c_int), parameter :: POSIX_FADV_DONTNEED    = 4
    integer(c_int), parameter :: POSIX_FADV_NOREUSE     = 5

contains
!******************************************************************************!
!> Open an uncompressed text file for reading
!******************************************************************************!
    subroutine open_file(this, filename, mode, stat)
        class(TextFileReader), intent(inout), target :: this
        character(*), intent(in) :: filename
        character(*), intent(in) :: mode
        type(FException), intent(out) :: stat
!------
        character(:), allocatable :: message
        type(c_ptr) :: c_message
        integer(c_int) :: ret
        integer(c_long) :: length
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

            ! Call the C library function
            this%stream = fopen(filename//C_NULL_CHAR, &
                "r"//C_NULL_CHAR)
            !TODO: make sure the error check is sufficient
                
            ! Error check
            if (.not.c_associated(this%stream)) then
                c_message = strerror(errno())
                call c_f_string(c_message, message)
                call stat%raise(FExceptionDescription(message, "when trying to open file " // bold(filename) // " for reading:"))
                return
            end if

            ret = fseek(this%stream, 0_c_long, 2_c_int)
            length = ftell(this%stream)
            call rewind_c(this%stream)
            ret = posix_fadvise(fileno(this%stream), 0_c_size_t, length, ior(POSIX_FADV_RANDOM, POSIX_FADV_WILLNEED))

            ! Set the mode
            this%read_mode = .true.
            this%write_mode = .false.

        ! Open the file for writing
        case ("w", "write")

            ! Call the C library function
            this%stream = fopen(filename//C_NULL_CHAR, &
                "w"//C_NULL_CHAR)
            !TODO: make sure the error check is sufficient
    
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

        ! Open an existing file for writing
        case ("a", "append")

            ! Call the C library function
            this%stream = fopen(filename//C_NULL_CHAR, &
                "a"//C_NULL_CHAR)
            !TODO: make sure the error check is sufficient

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

        ! Handle unknown modes
        case default
            call stat%raise(FExceptionDescription("Invalid mode", "when trying to open file " // bold(filename) // ":"))
            return
        end select
    end subroutine
!******************************************************************************!
!> Close an uncompressed text file. 
!******************************************************************************!
    subroutine close_file(this)
        class(TextFileReader), intent(inout) :: this
!------
        integer :: err
!------
        if (c_associated(this%stream)) then
            err = fclose(this%stream)
            this%stream = C_NULL_PTR
            !TODO: error check
        end if
     end subroutine
!******************************************************************************!
!> Read from an uncompressed text file.
!******************************************************************************!
    subroutine read(this, buffer, n, stat)
        class(TextFileReader), intent(inout) :: this !< File reader object
        character(:,c_char), allocatable, intent(inout) :: buffer
        integer, intent(out) :: n !< Number of characters read
        type(FException), intent(inout) :: stat !< Error status
!------
        integer(c_size_t) :: size_read, bloc_size, bloc_count
!------

        ! Allocate the buffer if needed
        if (allocated(buffer)) then
            if (len(buffer) /= BUFFER_LENGTH) deallocate(buffer)
        end if
        if (.not.allocated(buffer)) allocate(character(BUFFER_LENGTH,c_char) :: buffer)

        ! Stop here if the file is not open
        if (.not.c_associated(this%stream)) then
            n = 0
            call stat%raise("File not open.")
            return
        end if

        ! Stop here if the end of file was reached already
        if (feof(this%stream) /= 0) then
            n = 0
            call stat%raise("Trying to read past end of file.")
            return
        end if

        bloc_size = 1
        bloc_count = BUFFER_LENGTH

        ! Fill the buffer from the file
        size_read = fread(buffer, bloc_size, bloc_count, this%stream)

        ! Set the size to zero if an error occured
        if (ferror(this%stream)/=0 .and. feof(this%stream)==0) then
            size_read = -1
        end if
        
        ! Check what happened if the buffer was not filled
        if (size_read < bloc_count) then
        
            ! There was an error
            if (ferror(this%stream) /= 0) then
                call stat%raise(error_message())
            end if

            ! We reached the end of file
            if (feof(this%stream) /= 0) then
                call stat%raise("End of file.")
            end if
        end if

        n = int(size_read, kind(n))
        if (n /= size_read) then
            !TODO: handle long lines
        end if
     end subroutine
!******************************************************************************!
!> Write to an uncompressed text file.
!******************************************************************************!
    subroutine write(this, buffer, n, stat)
        class(TextFileReader), intent(inout) :: this
        character(*,c_char), intent(in) :: buffer
        integer, intent(in) :: n
        type(FException), intent(out) :: stat
!------
        integer(c_size_t) :: size_written, bloc_size, bloc_count
!------
        ! Stop here if the file is not open
        if (.not.c_associated(this%stream)) then
            call stat%raise("File not open.")
            return
        end if

        bloc_size = 1
        bloc_count = n

        ! Write the buffer to the file
        size_written = fwrite(buffer//new_line(buffer), bloc_size, bloc_count, this%stream)
        !TODO: make sure the error check is sufficient
        
        ! Check what happened if the whole buffer could not be written
        if (size_written /= len(buffer)) then
            call stat%raise(error_message())
        end if
    end subroutine
!******************************************************************************!
!> Set the file position indicator to the beginning of the file.
!******************************************************************************!
    subroutine rewind_(this)
        class(TextFileReader), intent(inout) :: this

        if (c_associated(this%stream)) call rewind_c(this%stream)
    end subroutine
!******************************************************************************!
!> Final subroutine, to make sure that the file is closed when the reader
!> disappears.
!******************************************************************************!
    elemental impure subroutine finalise(this)
        type(TextFileReader), intent(inout) :: this
!------
        if (c_associated(this%stream)) call close_file(this)
    end subroutine
end module
