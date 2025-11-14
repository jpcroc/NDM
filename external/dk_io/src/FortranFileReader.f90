!******************************************************************************!
!                         dk_fortranfilereader module
!------------------------------------------------------------------------------!
!> `[[FortranFileReader(type)]]` derived type and type-bound procedures, to add
!> support for uncompressed files to `[[FileObject(type)]]` using native
!> Fortran i/o facilities.
!
!  Part of the dk_io library version 0.1.
!  Written by Paul Fossati, <paul.fossati@cea.fr>
!  Copyright (c) 2024-2025 CEA
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
module dk_fortranfilereader
    use iso_c_binding, only: c_ptr, c_char, c_size_t, c_associated, C_NULL_PTR, C_NULL_CHAR

    use ext_character, only: c_f_string, bold
    use dk_exception, only: FException, FExceptionDescription
    use dk_filereader, only: FileReader, BUFFER_LENGTH

    implicit none(external, type)
    private

    type, extends(FileReader), public :: FortranFileReader
        private
        integer :: unit = 0
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
contains
!******************************************************************************!
!> Open an uncompressed text file for reading
!******************************************************************************!
    subroutine open_file(this, filename, mode, stat)
        class(FortranFileReader), intent(inout), target :: this
        character(*), intent(in) :: filename
        character(*), intent(in) :: mode
        type(FException), intent(out) :: stat
!------
        character(1000) :: iomsg
        integer :: iostat, i
!------

        ! Check whether the file is already open
        if (this%unit /= 0) then

            ! If yes, raise an error
            call stat%raise(FExceptionDescription("File reader already in use", "when trying to open file " // bold(filename) // ":"))
            return
        end if

        select case(mode)

        ! Open the file for reading
        case ("r", "read")

            ! Open the file
            open(file=filename, newunit=this%unit, action="read", status="old", access="stream", iostat=iostat, iomsg=iomsg)

            ! Error check
            if (iostat /= 0) then
                i = index(iomsg, ":")
                if (i > 0) iomsg = iomsg(i+2:)
                call stat%raise(FExceptionDescription(trim(iomsg), "when trying to open file " // bold(filename) // " for reading:"))
                return
            end if

            ! Set the mode
            this%read_mode = .true.
            this%write_mode = .false.

        ! Open the file for writing
        case ("w", "write")

            ! Open the file
            open(file=filename, newunit=this%unit, action="write", status="replace", access="stream", iostat=iostat, iomsg=iomsg)

            ! Error check
            if (iostat /= 0) then
                i = index(iomsg, ":")
                if (i > 0) iomsg = iomsg(i+2:)
                call stat%raise(FExceptionDescription(trim(iomsg), "when trying to open file " // bold(filename) // " for writing:"))
                return
            end if

            ! Set the mode
            this%read_mode = .false.
            this%write_mode = .true.

        ! Open the file for writing
        case ("a", "append")

            ! Open the file
            open(file=filename, newunit=this%unit, action="write", status="old", position="append", access="stream", iostat=iostat, iomsg=iomsg)

            ! Error check
            if (iostat /= 0) then
                i = index(iomsg, ":")
                if (i > 0) iomsg = iomsg(i+2:)
                call stat%raise(FExceptionDescription(trim(iomsg), "when trying to open existing file " // bold(filename) // " for writing:"))
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
        class(FortranFileReader), intent(inout) :: this
!------
        integer :: iostat
!------
        if (this%unit /= 0) then
            close(this%unit, iostat=iostat)
            this%unit = 0
            
            if (iostat /= 0) then
                !TODO: error check
            end if
        end if
     end subroutine
!******************************************************************************!
!> Read from an uncompressed text file.
!******************************************************************************!
    subroutine read(this, buffer, n, stat)
        use iso_fortran_env, only: int64

        class(FortranFileReader), intent(inout) :: this !< File reader object
        character(:,c_char), allocatable, intent(inout) :: buffer
        integer, intent(out) :: n !< Number of characters read
        type(FException), intent(inout) :: stat !< Error status
!------
        integer :: iostat
        integer(int64) :: pos
!------

        ! Allocate the buffer if needed
        if (.not.allocated(buffer)) allocate(character(BUFFER_LENGTH,c_char) :: buffer)

        ! Stop here if the file is not open
        if (this%unit == 0) then
            call stat%raise("File not open.")
            return
        end if

        inquire(this%unit, pos=pos)
        buffer(:) = "  "
        read(this%unit, iostat=iostat) buffer

        if (iostat < 0) then
            ! We reached the end of file and need to re-read carefully

            ! Reposition the file pointer
            iostat = 0
            read(this%unit, pos=pos, iostat=iostat) buffer(1:1)

            if (iostat /= 0) then
                call stat%raise("Trying to read past end of file.")
                n = 0
                return
            end if

            ! Read to the end of the file character by character
            do n=2, len(buffer)
                read(this%unit, iostat=iostat) buffer(n:n)
                if (iostat /= 0) then
                    exit
                end if
            end do
        else if (iostat > 0) then
            call stat%raise("Error when reading file")
            n = 0
        else
             n = len(buffer)
        end if
     end subroutine
!******************************************************************************!
!> Write to an uncompressed text file.
!******************************************************************************!
    subroutine write(this, buffer, n, stat)
        class(FortranFileReader), intent(inout) :: this
        character(*,c_char), intent(in) :: buffer
        integer, intent(in) :: n
        type(FException), intent(out) :: stat
!------
!------
        ! Stop here if the file is not open
        if (this%unit == 0) then
            call stat%raise("File not open.")
            return
        end if

        write(this%unit) buffer(:n)
        !TODO: error check
        
    end subroutine
!******************************************************************************!
!> Set the file position indicator to the beginning of the file.
!******************************************************************************!
    subroutine rewind_(this)
        class(FortranFileReader), intent(inout) :: this

        if (this%unit /= 0) then
            rewind(this%unit)
        end if
    end subroutine
!******************************************************************************!
!> Final subroutine, to make sure that the file is closed when the reader
!> disappears.
!******************************************************************************!
    elemental impure subroutine finalise(this)
        type(FortranFileReader), intent(inout) :: this
!------
        if (this%unit /= 0) then
            close(this%unit)
        end if
    end subroutine
end module
