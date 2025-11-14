!******************************************************************************!
!                         dk_mmapfilereader module
!------------------------------------------------------------------------------!
!> `[[MmapFileReader(type)]]` derived type and type-bound procedures, to add
!> support for uncompressed files to `[[FileObject(type)]]` using mmap.
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
module dk_mmapfilereader
    use iso_c_binding, only: c_ptr, c_char, c_size_t, c_associated, C_NULL_PTR, C_NULL_CHAR, c_int, c_long

    use ext_character, only: c_f_string, bold, c_f_buffer
    use dk_exception, only: FException, FExceptionDescription
    use dk_filereader, only: FileReader, BUFFER_LENGTH
    use dk_posix_io, only: feof, ferror, fflush, fopen, fclose, fwrite, errno, strerror, fgets, fread, error_message, rewind_c => rewind, fclose, fileno, fseek, ftell

    use mman

    implicit none(external, type)
    private

    type, extends(FileReader), public :: MmapFileReader
        private
        type(c_ptr) :: stream = C_NULL_PTR !< Pointer to the C stream
        type(c_ptr) :: map_start = C_NULL_PTR
        integer(c_size_t) :: map_size = 0
        integer(c_long) :: file_pos = 0
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
    function mmap(addr, length, prot, flags, fd, offset) bind(C, name="mmap") result(res)
        use iso_c_binding
        type(c_ptr), intent(in), value :: addr
        integer(c_size_t), intent(in), value :: length
        integer(c_int), intent(in), value :: prot
        integer(c_int), intent(in), value :: flags
        integer(c_int), intent(in), value :: fd
        integer(c_long), intent(in), value :: offset
        type(c_ptr) :: res
    end function mmap
end interface
contains
!******************************************************************************!
!> Open an uncompressed mmap-ed text file for reading
!******************************************************************************!
    subroutine open_file(this, filename, mode, stat)
        class(MmapFileReader), intent(inout), target :: this
        character(*), intent(in) :: filename
        character(*), intent(in) :: mode
        type(FException), intent(out) :: stat
!------
        character(:), allocatable :: message
        type(c_ptr) :: c_message, file
        integer(c_int) :: fd, ret
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
            file = fopen(filename//C_NULL_CHAR, &
                "r"//C_NULL_CHAR)

            if (.not.c_associated(file)) then
                c_message = strerror(errno())
                call c_f_string(c_message, message)
                call stat%raise(FExceptionDescription(message, "when trying to open file " // bold(filename) // " for reading:"))
                return
            end if

            fd = fileno(file)
            ret = fseek(file, 0_c_long, 2_c_int)
            this%map_size = ftell(file)
            call rewind_c(file)
            this%map_start = mmap(C_NULL_PTR, this%map_size, PROT_READ, MAP_PRIVATE, fd, 0_c_size_t)
            this%file_pos = 0
            ret = fclose(file)

            !TODO: make sure the error check is sufficient
            !TODO: use proper error checks for mmap

            ! Error check
            if (.not.c_associated(this%map_start)) then
                c_message = strerror(errno())
                call c_f_string(c_message, message)
                call stat%raise(FExceptionDescription(message, "when trying to open file " // bold(filename) // " for reading:"))
                return
            end if

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
!> Close an uncompressed mmap-ed text file.
!******************************************************************************!
    subroutine close_file(this)
        class(MmapFileReader), intent(inout) :: this
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
!> Read from an uncompressed mmap-ed text file.
!******************************************************************************!
    subroutine read(this, buffer, n, stat)
        use pointer, only: c_ptr_add
        use ext_character, only: c_f_string
        class(MmapFileReader), intent(inout) :: this !< File reader object
        character(:,c_char), allocatable, intent(inout) :: buffer
        integer, intent(out) :: n !< Number of characters read
        type(FException), intent(inout) :: stat !< Error status
!------
        type(c_ptr) :: buffptr
!------

        ! Allocate the buffer if needed
        if (.not.allocated(buffer)) allocate(character(BUFFER_LENGTH,c_char) :: buffer)

        ! Stop here if the mapping was not set properly
        if (.not.c_associated(this%map_start)) then
            call stat%raise("File not open.")
            return
        end if

        ! Stop here if all the mapping was read already
        if (this%file_pos >= this%map_size) then
            call stat%raise("Trying to read past end of file.")
            return
        end if

        ! Return the end of the file if there is not enough remaining data to fill the buffer
        if (this%file_pos+len(buffer) > this%map_size) then
            buffptr = c_ptr_add(this%map_start, this%file_pos)
            call c_f_buffer(buffptr, buffer, this%map_size-this%file_pos)
            n = int(this%map_size-this%file_pos)
            this%file_pos = this%map_size
            return
        end if

        buffptr = c_ptr_add(this%map_start, this%file_pos)
        call c_f_buffer(buffptr, buffer, len(buffer, kind=c_size_t))
        n = len(buffer)
        this%file_pos = this%file_pos + len(buffer)
     end subroutine
!******************************************************************************!
!> Write to an uncompressed mmap-ed text file.
!******************************************************************************!
    subroutine write(this, buffer, n, stat)
        class(MmapFileReader), intent(inout) :: this
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
        class(MmapFileReader), intent(inout) :: this

        if (this%read_mode) then
            this%file_pos = 0
        else
            if (c_associated(this%stream)) call rewind_c(this%stream)
        end if
    end subroutine
!******************************************************************************!
!> Final subroutine, to make sure that the file is closed when the reader
!> disappears.
!******************************************************************************!
    elemental impure subroutine finalise(this)
        type(MmapFileReader), intent(inout) :: this
!------
        if (c_associated(this%stream)) call close_file(this)
    end subroutine
end module
