!******************************************************************************!
!                         dk_fileobject module
!------------------------------------------------------------------------------!
!> Contains the derived type [[FileObject(type)]] and its type-bound procedures,
!> to read and  write files in various compression formats. Also contains the
!> [[FileError(type)]] derived type for error reporting.
!
!  Part of the dk_io library version 0.1.
!  Written by Paul Fossati, <paul.fossati@cea.fr>
!  Copyright (c) 2024-2005 CEA
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
module dk_fileobject
    use iso_fortran_env, only: int8, IOSTAT_END

    use ext_character, only: bold
    use dk_exception, only: FExceptionDescription, FException
    use event_levels
    use ext_character
    use dk_filereader, only: FileReader, BUFFER_LENGTH
    use dk_textfilereader,  only: TextFileReader
    use dk_fortranfilereader, only: FortranFileReader
    use dk_mmapfilereader, only: MmapFileReader
#ifdef HAVE_LIBZ
    use dk_gzipfilereader,  only: GzipFileReader
#endif
#ifdef HAVE_LIBBZ2
    use dk_bzip2filereader, only: Bzip2FileReader
#endif
#ifdef HAVE_LIBZSTD
    use dk_zstdfilereader, only: ZstdFileReader
#endif

    implicit none(external, type)
    private

    ! Modes for open files
    character(*), parameter :: MODE_UNSET  = "UNDEF"
    character(*), parameter :: MODE_READ   = "read"
    character(*), parameter :: MODE_WRITE  = "write"
    character(*), parameter :: MODE_APPEND = "append"

    ! Magic numbers for various compression formats
    character(*), parameter, private :: MAGIC_BZIP2 = "42 5A 68"
    character(*), parameter, private :: MAGIC_GZIP  = "1F 8B 08"
    character(*), parameter, private :: MAGIC_ZSTD  = "28 B5 2F FD"

!******************************************************************************!
!> Derived type holding the description of a file object.
!******************************************************************************!
    type, public :: FileObject
        private

        character(:), allocatable, public :: name !< Full file name, including extensions
        integer, public :: current_line = 0 !< Counter indicating the number of the current line in the file

        class(FileReader), allocatable :: reader

        character(:), allocatable :: path !< Path of the directory containing the file
        character(:), allocatable :: ext  !< File extension (including the compression extension if any)

        character(:), allocatable :: buffer !< Internal buffer
        integer :: buffer_top = 1 !< Cursor to the available space in the buffer
        integer :: buffer_bottom = 1

        character(8) :: mode = MODE_UNSET !< File mode, indicating whether the file was opened for reading, writing, or appending

        logical :: exists = .false. !< Whether the file exists
        logical :: opened = .false. !< Whether the file was opened
    contains
        procedure, public :: init_fortran
        procedure, public :: init_mmap
        procedure, public :: init => init_file
        procedure, public :: close
        procedure, public :: full_name
        procedure, public :: line_number
        procedure, public :: read_line
        procedure, public :: read_nonempty_line
        procedure, public :: write_line
        procedure, public :: rewind !< Move the file cursor back to the beginning of the file
        procedure, public :: compression_format !< Get the compression format of the file

        final :: finalise !< Clean up when the file object is destroyed
    end type

!******************************************************************************!
!> I/O error description.
!******************************************************************************!
    type, public, extends(FExceptionDescription) :: FileError
        private
        character(:), allocatable :: filename !< Name of the file related to the error condition
        integer :: line = 0 !< Number of the line where the error happened
        character(:), allocatable :: content !< Context
    contains
        procedure, public :: string => fileError_string
    end type

!******************************************************************************!
!> I/O error constructor
!******************************************************************************!
    interface FileError
        module procedure :: FileError_init
    end interface
    public :: fileError_init

contains
    pure function compression_format(this)
        class(FIleObject), intent(in) :: this
        character(:), allocatable :: compression_format

        if (allocated(this%reader)) then
            select type(r => this%reader)
            class is (FortranFileReader)
                compression_format = "ASCII, Fortran reader"
            class is (MMapFileReader)
                compression_format = "ASCII, mmap reader"
            class is (TextFileReader)
                compression_format = "ASCII"
#ifdef HAVE_LIBZSTD
            class is (Bzip2FileReader)
                compression_format = "bzip2"
#endif
#ifdef HAVE_LIBZSTD
            class is (GzipFileReader)
                compression_format = "gzip"
#endif
#ifdef HAVE_LIBZSTD
            class is (ZstdFileReader)
                compression_format = "zstandard"
#endif
            end select
        end if
    end function
!******************************************************************************!
!> Make sure the file is closed properly before a file object is destroyed.
!******************************************************************************!
    elemental impure subroutine finalise(this)
        type(FileObject), intent(inout) :: this
!------
        type(FException) :: stat_
!------
        ! Close the file properly before the object is destroyed
        call this%close(stat_)
        call stat_%discard
    end subroutine

    subroutine close(this, stat)
        class(FileObject), intent(inout) :: this
        type(FException), intent(out), optional :: stat
!------
        type(FException) :: stat_
!------

        if (.not.allocated(this%reader)) return

        ! Write the rest of the buffer if needed before the file object is destroyed
        if (this%mode == MODE_WRITE .or. this%mode == MODE_APPEND) then
            call this%reader%write(this%buffer(this%buffer_bottom:this%buffer_top-1), (this%buffer_top-this%buffer_bottom), stat_)
            if (stat_ /= 0 .and. stat_ /= "File not open.") then
                write(*,*) "Error in final WRITE"
                call stat_%report
            end if
            this%buffer_bottom = 1
            this%buffer_top = 1
            if (stat_ == "File not open.") call stat_%discard
        end if

        call this%reader%close_file
        this%current_line = 0
        if (present(stat)) call stat%transfer(stat_)
    end subroutine
!******************************************************************************!
!> Initialise a file object to read or write an uncompressed text file using
!> the native Fortran facilities.
!******************************************************************************!
    subroutine init_fortran(this, path, mode, stat)
        class(FileObject), intent(out) :: this
        character(*), intent(in) :: path
        character(*), intent(in) :: mode
        type(FException), intent(out), optional :: stat
!------
        integer :: i
        type(FException) :: istat
!------

        ! Get the file extension
        i = index(path, ".", back=.true.)
        if (i > 0) then
            this%ext = path(i+1:)
            if (this%ext=="gz" .or. this%ext=="bz2" .or. this%ext=="zstd") then
                i = index(path(:i-1), ".", back=.true.)
                if (i > 0) this%ext = path(i+1:)
            end if
        end if

        ! Get the file name and path
        i = index(path, "/", back=.true.)
        if (i > 0) then
            this%path = path(:i-1)
            this%name = path(i+1:)
        else
            this%name = path
        end if

        ! Initialise the internal file reader object
        allocate(FortranFileReader :: this%reader)

        ! Set the mode
        select case(mode)
        case ("r", "R", "read", "READ")
            this%mode = MODE_READ
            call this%reader%open_file(this%full_name(), "r", istat)
        case ("w", "W", "write", "WRITE")
            this%mode = MODE_WRITE
            allocate(character(BUFFER_LENGTH) :: this%buffer)
            call this%reader%open_file(this%full_name(), "w", istat)
        case ("a", "A", "append", "APPEND")
            this%mode = MODE_APPEND
            allocate(character(BUFFER_LENGTH) :: this%buffer)
            call this%reader%open_file(this%full_name(), "a", istat)
        case default
            call istat%raise("wrong file mode: " // this%mode // ".")
        end select
    end subroutine
!******************************************************************************!
!> Initialise a file object to read or write an uncompressed text file using
!> mmap().
!******************************************************************************!
    subroutine init_mmap(this, path, mode, stat)
        class(FileObject), intent(out) :: this
        character(*), intent(in) :: path
        character(*), intent(in) :: mode
        type(FException), intent(out), optional :: stat
!------
        integer :: i
        type(FException) :: istat
!------

        ! Get the file extension
        i = index(path, ".", back=.true.)
        if (i > 0) then
            this%ext = path(i+1:)
            if (this%ext=="gz" .or. this%ext=="bz2" .or. this%ext=="zstd") then
                i = index(path(:i-1), ".", back=.true.)
                if (i > 0) this%ext = path(i+1:)
            end if
        end if

        ! Get the file name and path
        i = index(path, "/", back=.true.)
        if (i > 0) then
            this%path = path(:i-1)
            this%name = path(i+1:)
        else
            this%name = path
        end if

        ! Initialise the internal file reader object
        allocate(MmapFileReader :: this%reader)

        ! Set the mode
        select case(mode)
        case ("r", "R", "read", "READ")
            this%mode = MODE_READ
            call this%reader%open_file(this%full_name(), "r", istat)
        case ("w", "W", "write", "WRITE")
            this%mode = MODE_WRITE
            allocate(character(BUFFER_LENGTH) :: this%buffer)
            call this%reader%open_file(this%full_name(), "w", istat)
        case ("a", "A", "append", "APPEND")
            this%mode = MODE_APPEND
            allocate(character(BUFFER_LENGTH) :: this%buffer)
            call this%reader%open_file(this%full_name(), "a", istat)
        case default
            call istat%raise("wrong file mode: " // this%mode // ".")
        end select
    end subroutine
!******************************************************************************!
!> Initialise a file object.
!******************************************************************************!
    subroutine init_file(this, path, mode, stat)
        class(FileObject), intent(out) :: this
        character(*), intent(in) :: path
        character(*), intent(in) :: mode
        type(FException), intent(out), optional :: stat
!------
        integer :: i
        type(FException) :: istat
!------

        ! Get the file extension
        i = index(path, ".", back=.true.)
        if (i > 0) then
            this%ext = path(i+1:)
            if (this%ext=="gz" .or. this%ext=="bz2" .or. this%ext=="zst") then
                i = index(path(:i-1), ".", back=.true.)
                if (i > 0) this%ext = path(i+1:)
            end if
        end if

        ! Get the file name and path
        i = index(path, "/", back=.true.)
        if (i > 0) then
            this%path = path(:i-1)
            this%name = path(i+1:)
        else
            this%name = path
        end if

        ! Initialise the current line counter
        this%current_line = 0

        ! Set the mode and initialise the internal file reader object
        select case(mode)
        case ("r", "R", "read", "READ")
            this%mode = MODE_READ
            call open_read(this, istat)
        case ("w", "W", "write", "WRITE")
            this%mode = MODE_WRITE
            allocate(character(BUFFER_LENGTH) :: this%buffer)
            call open_write(this, istat)
        case ("a", "A", "append", "APPEND")
            this%mode = MODE_APPEND
            allocate(character(BUFFER_LENGTH) :: this%buffer)
            call open_append(this, istat)
        case default
            call istat%raise("wrong file mode: " // this%mode // ".")
        end select

        ! Initialise the current line counter
        this%current_line = 0

        ! Report any issue
        if (present(stat)) call stat%transfer(istat)
    end subroutine
!******************************************************************************!
!> Open a file for reading.
!******************************************************************************!
    subroutine open_read(this, stat)
        class(FileObject), intent(inout) :: this
        type(FException), intent(out) :: stat
!------
        integer :: unit, iostat
        character(1000) :: iomsg
        integer(int8), dimension(6) :: magicBytes
        character(17) :: magicString
        logical :: exists
!------
        ! Check whether the file exists
        inquire(file=this%full_name(), exist=exists)

        ! Stop here if the file does not exist
        if (.not.exists) then
            call stat%raise(FExceptionDescription("No such file or directory.", &
                "could not open file " // bold(this%full_name()) // " for reading:"))
            return
        end if

        ! The file exists, try to open it
        open(newunit=unit, &
                file=this%full_name(), &
              action="read", &
              status="old", &
              access="stream", &
                form="unformatted", &
              iostat=iostat, &
               iomsg=iomsg)

        ! Check whether the file could be opened
        if (iostat /= 0) then
            call stat%raise(FExceptionDescription( &
                trim(iomsg), &
                "could not open file " // this%full_name() // ":"))
            close(unit, iostat=iostat)
            return
        end if

        ! Try to read the magic bytes from the file
        read(unit, iostat=iostat) magicBytes
        close(unit)

        if (iostat == 0) then
            ! Compare the magic string to some known compression formats
            write(magicString,'(z2.2,5(1x,z2.2))') magicBytes
            if (magicString(:len(MAGIC_BZIP2)) == MAGIC_BZIP2) then
                ! The file looks like a bzip2 archive
#ifdef HAVE_LIBBZ2
                allocate(Bzip2FileReader :: this%reader)
#else
                call stat%raise("The dk_io library was not compiled with support for bzip2 compressed files.")
#endif
            else if (magicString(:len(MAGIC_GZIP)) == MAGIC_GZIP) then
                ! The file looks like a gzip archive
#ifdef HAVE_LIBZ
                allocate(GzipFileReader :: this%reader)
#else
                call stat%raise("The dk_io library was not compiled with support for gzip compressed files.")
#endif
            else if (magicString(:len(MAGIC_ZSTD)) == MAGIC_ZSTD) then
                ! The file looks like a zstd archive
#ifdef HAVE_LIBZSTD
                allocate(ZstdFileReader :: this%reader)
#else
                call stat%raise("The dk_io library was not compiled with support for zstd compressed files.")
#endif
            else
                ! Assume that the file is not compressed
                allocate(TextFileReader :: this%reader)
            end if
        else if (iostat == IOSTAT_END) then
            ! Assume that the file is not compressed
            allocate(TextFileReader :: this%reader)
        else
            ! LCOV_EXCL_START
            call stat%raise(FExceptionDescription( &
                "Unknown error, please report the circumstances to fix this", &
                "could not read from file " // this%full_name() // ":"))
            close(unit, iostat=iostat)
            return
            ! LCOV_EXCL_STOP
        end if

        ! Initialise the file reader
        if (allocated(this%reader)) call this%reader%open_file(this%full_name(), "r", stat)
    end subroutine
!******************************************************************************!
!> Open a file for writing.
!******************************************************************************!
    subroutine open_write(this, stat)
        class(FileObject), intent(inout) :: this
        type(FException), intent(out) :: stat
!------
        integer :: i
        logical :: exists
!------
        if (allocated(this%path)) then
#ifdef INQUIRE_DIR
            inquire(INQUIRE_DIR=this%path, exist=exists)
#else
            error stop "unsupported compiler"
#endif

            ! Stop here if the directory in which the file is to be written does not exist
            if (.not.exists) then
                call stat%raise(FExceptionDescription("No such file or directory.", &
                    "could not open file " // bold(this%full_name()) // " for writing:"))
                return
            end if
        end if

        ! Guess the compression format from the file extension and set the file reader accordingly
        if (allocated(this%ext)) then
            i = index(this%ext, ".", back=.true.)
            select case(lower_case(this%ext(i+1:)))
            case ("bz2")
#ifdef HAVE_LIBBZ2
                allocate(Bzip2FileReader :: this%reader)
#else
                call stat%raise("The dk_io library was not compiled with support for bzip2 compressed files.")
#endif
            case ("gz")
#ifdef HAVE_LIBZ
                allocate(GzipFileReader :: this%reader)
#else
                call stat%raise("The dk_io library was not compiled with support for gzip compressed files.")
#endif
            case ("zst", "zstd")
#ifdef HAVE_LIBZSTD
                allocate(ZstdFileReader :: this%reader)
#else
                call stat%raise("The dk_io library was not compiled with support for zstd compressed files.")
#endif
            case default
                allocate(TextFileReader :: this%reader)
            end select
        else
            allocate(TextFileReader :: this%reader)
        end if

        ! Initialise the file reader
        if (allocated(this%reader)) then
            call this%reader%open_file(this%full_name(), "w", stat)
            call stat%report
        end if
    end subroutine
!******************************************************************************!
!> Open an existing file for appending.
!******************************************************************************!
    subroutine open_append(this, stat)
        class(FileObject), intent(inout) :: this
        type(FException), intent(out) :: stat
!------
        integer :: unit, iostat, i
        character(1000) :: iomsg
        integer(int8), dimension(6) :: magicBytes
        character(17) :: magicString
        logical :: exists
!------
        if (allocated(this%path)) then

#ifdef INQUIRE_DIR
            inquire(INQUIRE_DIR=this%path, exist=exists)
#else
            error stop "unsupported compiler"
#endif

            ! Stop here if the directory in which the file is to be written does not exist
            if (.not.exists) then
                call stat%raise(FExceptionDescription("No such file or directory.", &
                    "could not open file " // bold(this%full_name()) // " for writing:"))
                return
            end if
        end if

        inquire(file=this%full_name(), exist=exists)

        ! Check whether the file could be opened
        if (exists) then

            ! If the file exists, we guess the compression format from its magic bytes

            open(newunit=unit, &
                    file=this%full_name(), &
                action="read", &
                status="old", &
                access="stream", &
                    form="unformatted", &
                iostat=iostat, &
                iomsg=iomsg)

            ! Check whether the file could be opened
            if (iostat /= 0) then
                call stat%raise(FExceptionDescription( &
                    trim(iomsg), &
                    "could not open file " // this%full_name() // ":"))
                close(unit, iostat=iostat)
                return
            end if

            ! Try to read the magic bytes from the file
            read(unit, iostat=iostat) magicBytes
            close(unit)

            ! Compare the magic string to some known compression formats
            write(magicString,'(z2.2,5(1x,z2.2))') magicBytes
            if (magicString(:len(MAGIC_BZIP2)) == MAGIC_BZIP2) then
                ! The file looks like a bzip2 archive
#ifdef HAVE_LIBBZ2
                allocate(Bzip2FileReader :: this%reader)
#else
                call stat%raise("The dk_io library was not compiled with support for bzip2 compressed files.")
#endif
            else if (magicString(:len(MAGIC_GZIP)) == MAGIC_GZIP) then
                ! The file looks like a gzip archive
#ifdef HAVE_LIBZ
                allocate(GzipFileReader :: this%reader)
#else
                call stat%raise("The dk_io library was not compiled with support for gzip compressed files.")
#endif
            else if (magicString(:len(MAGIC_ZSTD)) == MAGIC_ZSTD) then
                ! The file looks like a zstd archive
#ifdef HAVE_LIBZSTD
                allocate(ZstdFileReader :: this%reader)
#else
                call stat%raise("The dk_io library was not compiled with support for zstd compressed files.")
#endif
            else
                ! Assume that the file is not compressed
                allocate(TextFileReader :: this%reader)
            end if

        else
            ! If the file does not exist, we guess the compression format from the file extension

            i = index(this%ext, ".", back=.true.)
            select case(lower_case(this%ext(i+1:)))
            case ("bz2")
#ifdef HAVE_LIBBZ2
                allocate(Bzip2FileReader :: this%reader)
#else
                call stat%raise("The dk_io library was not compiled with support for bzip2 compressed files.")
#endif
            case ("gz")
#ifdef HAVE_LIBZ
                allocate(GzipFileReader :: this%reader)
#else
                call stat%raise("The dk_io library was not compiled with support for gzip compressed files.")
#endif
            case ("zip")
                call stat%raise("The zip compression format is not supported yet")
            case ("7z")
                call stat%raise("The 7-zip compression format is not supported yet")
            case ("zstd", "zst")
#ifdef HAVE_LIBZSTD
                allocate(ZstdFileReader :: this%reader)
#else
                call stat%raise("The dk_io library was not compiled with support for zstd compressed files.")
#endif
            case default
                allocate(TextFileReader :: this%reader)
            end select
        end if

        ! Initialise the file reader
        if (allocated(this%reader)) then
            call this%reader%open_file(this%full_name(), "a", stat)
            call stat%report
        end if
    end subroutine
!******************************************************************************!
!> Get the full name of a file
!******************************************************************************!
    pure function full_name(this)
        class(FileObject), intent(in) :: this
        character(:), allocatable :: full_name
!------
        if (allocated(this%path) .and. allocated(this%name)) then
            full_name = this%path // "/" // this%name
        else if (.not.allocated(this%path) .and. allocated(this%name)) then
            full_name = this%name
        else if (.not.allocated(this%path) .and. .not.allocated(this%name)) then
            full_name = "<unnamed file>"
        end if
    end function
!******************************************************************************!
!> Get the current line number
!******************************************************************************!
    pure function line_number(this)
        class(FileObject), intent(in) :: this
        integer :: line_number
!------
        line_number = this%current_line
    end function
!******************************************************************************!
!> Read a line from a file into a character variable. This subroutine returns
!> the first non-empty line that is found, if any. The arguments have the same
!> semantics as in [[read_line(subroutine)]].
!>
!> @warning
!> No attempt is made to check whether the comment markers are in multi-line
!> character strings. This subroutine should not be used if the file being
!> read can have those.
!> @endwarning
!******************************************************************************!
    subroutine read_nonempty_line(this, line, comment_chars, trim, stat)
        class(FileObject), intent(inout) :: this !< File to read
        character(:), allocatable, intent(inout) :: line !< Variable containing the next non-empty and uncommented line in the file
        character(*), intent(in), optional :: comment_chars(:) !< Array of characters marking comments in the file
        logical, intent(in), optional :: trim !< Flag indicating whether trailing whitespace is acceptable or not
        type(FException), intent(out), optional :: stat
!------
        logical :: is_trimmed
        type(FException) :: istat
!------
        if (present(trim)) then
            is_trimmed = trim
        else
            is_trimmed = .false. ! Default behavior is not to trim
        end if

        call read_line(this, line, stat=istat)
        do while (is_comment_or_empty(line, comment_chars) .and. istat == 0)
            call read_line(this, line, stat=istat)
        end do

        if (present(stat)) call stat%transfer(istat)
    end subroutine
    function is_comment_or_empty(str, comment_chars) result(flag)
        character(*), intent(in) :: str
        character, optional, intent(in) :: comment_chars(:)
        logical :: flag
!------
        integer :: i
!------
        flag = .false.

        if (len_trim(str) == 0) then
            flag = .true.
        else if (present(comment_chars)) then
            do i = 1, size(comment_chars)
                if (index(str, comment_chars(i)) > 0) then
                    if (str(:index(str, comment_chars(i))-1) == " ") then
                        flag = .true.
                        return
                    end if
                end if
            end do
        end if
    end function
!******************************************************************************!
!> Read a line from a file into a character variable. The line variable is only
!> reallocated if needed to improve performance. If the variable is longer than
!> the actual line read from file, it is padded with spaces. A `trim` optional
!> argument can be passed to force the line to be reallocated and ensure that
!> it does not have any padding, at a performance cost.
!******************************************************************************!
    subroutine read_line(this, line, trim, stat)
        class(FileObject), intent(inout) :: this !< File to read
        character(:), allocatable, intent(inout) :: line !< Variable containing the next line in the file
        logical, intent(in), optional :: trim !< Flag indicating whether trailing whitespace is acceptable or not
        type(FException), intent(out), optional :: stat
!------
        integer :: i, n, m, length
        type(FException) :: stat_
!------
        body: block

            length = 0

            ! Check whether the file was open
            if (.not.allocated(this%reader)) then
                call stat_%raise(FExceptionDescription("File not open.", &
                    "could not read from file " // bold(this%full_name()) // ":"))
                exit body
            end if

            ! Check whether the file was open for reading
            if (.not.this%mode == MODE_READ) then
                call stat_%raise(FExceptionDescription("File not open for reading.", &
                    "could not read from file " // bold(this%full_name()) // ":"))
                exit body
            end if

            ! Deallocate the line variable if it should be trimmed for the output
            if (present(trim)) then
                if (trim .and. allocated(line)) deallocate(line)
            end if

            ! Check whether a full line is in the buffer
            i = 0
            do m=this%buffer_bottom, this%buffer_top-1
                if (this%buffer(m:m) == new_line(this%buffer)) then
                    i = m
                    exit
                end if
            end do

            ! Read data into the buffer if it does not contain a full line
            do while (i == 0 .and. stat_ == 0)

                ! Move the remaining data in the internal buffer to the output string
                if (this%buffer_top > this%buffer_bottom) then
                    call concatenate_strings(line, length, this%buffer(this%buffer_bottom:this%buffer_top-1))
                end if
                this%buffer_bottom = 1
                this%buffer_top = 1

                ! Re-fill the buffer
                call this%reader%read(this%buffer, n, stat_)
                this%buffer_top = this%buffer_bottom + n

                ! Check whether a full line is in the buffer
                i = 0
                do m=this%buffer_bottom, this%buffer_top-1
                    if (this%buffer(m:m) == new_line(this%buffer)) then
                        i = m
                        exit
                    end if
                end do

                ! Handle errors and end of file
            end do

            if (i > 0) then
                call concatenate_strings(line, length, this%buffer(this%buffer_bottom:i-1))
                this%buffer_bottom = i+1

                ! Increment the line counter
                this%current_line = this%current_line + 1
            else
                ! We went to the end of the file without seeing a new line
                call concatenate_strings(line, length, this%buffer(this%buffer_bottom:this%buffer_top-1))
                this%buffer_bottom = this%buffer_top+1

                ! Increment the line counter
                this%current_line = this%current_line + 1
            end if

            if (stat_ == "End of file.") then
                call stat_%discard
            end if
        end block body
        if (allocated(line)) then
            if (length < len(line)) line(length+1:) = ""
        end if

        if (present(stat)) call stat%transfer(stat_)
    contains
        subroutine concatenate_strings(string1, length, string2)
            character(:), allocatable, intent(inout) :: string1
            integer, intent(inout) :: length
            character(*), intent(in) :: string2
!------
            if (.not.allocated(string1)) then
                string1 = string2
                length = len(string2)
                return
            end if
            if (length + len(string2) > len(string1)) then
                call realloc_string(string1, length + len(string2))
            end if
            string1(length+1:length+len(string2)) = string2
            length = length+len(string2)
        end subroutine

        subroutine realloc_string(string, length)
            character(:), allocatable, intent(inout) :: string
            integer, intent(in) :: length
!------
            character(:), allocatable :: tmp
!------
            if (allocated(string)) then
                call move_alloc(from=string, to=tmp)
                allocate(character(length) :: string)
                string(:len(tmp)) = tmp
            else
                allocate(character(length) :: string)
                return
            end if
        end subroutine
    end subroutine
!******************************************************************************!
!> Write a character variable as a line in a file.
!******************************************************************************!
    subroutine write_line(this, line, stat)
        class(FileObject), intent(inout) :: this
        character(*), intent(in) :: line
        type(FException), intent(out), optional :: stat
!------
        type(FException) :: stat_
        integer :: available, remaining, i
!------
        body: block

            ! Check whether the file was open
            if (.not.allocated(this%reader)) then
                call stat%raise(FExceptionDescription("File not open.", &
                    "could not write to file " // bold(this%full_name()) // ":"))
                exit body
            end if

            ! Check whether the file was open for reading
            if (this%mode/=MODE_WRITE .and. this%mode/=MODE_APPEND) then
                call stat%raise(FExceptionDescription("File not open for writing.", &
                    "could not write to file " // bold(this%full_name()) // ":"))
                exit body
            end if

            this%current_line = this%current_line + 1

            ! call write_to_buffer(this, line, stat_)
            i = 0
            remaining = len(line)+1-i ! Remaining length of the line
            available = len(this%buffer)-this%buffer_top+1 ! Space remaining in the buffer
            do while(remaining > available)
                this%buffer(this%buffer_top:) = line(i+1:i+available)
                call this%reader%write(this%buffer, len(this%buffer), stat_)
                if (stat_ /= 0) exit body
                i = i+available
                this%buffer_top = 1
                remaining = len(line)+1-i
                available = len(this%buffer)-this%buffer_top+1
            end do

            ! if (stat_ == 0) then
                this%buffer(this%buffer_top:) = line(i+1:)
                this%buffer_top = this%buffer_top + len(line(i+1:))
                this%buffer(this%buffer_top:this%buffer_top) = new_line(this%buffer)
                this%buffer_top = this%buffer_top + 1
            ! end if
        end block body

        ! Report any issue
        if (present(stat)) call stat%transfer(stat_)
    end subroutine
!******************************************************************************!
!> Move the file pointer to the begining of the file.
!******************************************************************************!
    subroutine rewind(this, stat)
        class(FileObject), intent(inout) :: this
        type(FException), intent(out), optional :: stat
!------
        type(FException) :: stat_
!------
        if (allocated(this%reader)) then

            ! Rewind the underlying file pointer
            call this%reader%rewind!(stat_)

            ! Reset the buffer and the line counter
            this%buffer_bottom = 1
            this%buffer_top = 1
            this%current_line = 0

            if (stat_ /= 0) then
                !TODO: error handling
            end if
        end if

        if (present(stat)) call stat%transfer(stat_)
    end subroutine
!******************************************************************************!
!> Set up a file error object based on a file and a user-provided message.
!******************************************************************************!
    function fileError_init(file, message, context, level) result(this)
        class(FileObject), intent(in) :: file
        character(*), intent(in) :: message
        integer, intent(in), optional :: level
        character(*), intent(in), optional :: context
        type(FileError) :: this
!------
        character(:), allocatable :: icontext
        integer :: ilevel
!------
        if (allocated(file%name) .and. allocated(file%path)) then
            this%filename = file%path // "/" // file%name
        else if (allocated(file%name)) then
            this%filename = file%name
        else
            this%filename = "<unnamed file>"
        end if
        this%line = file%line_number()

        ! Get the context (file name and possibly line number)
        if (this%line == 0) then
            icontext = bold(this%filename // ":")
        else
            icontext = bold(this%filename // ":" // this%line // ":")
        end if

        if (present(context)) icontext = icontext // NEW_LINE("") // context

        ! Set the message level
        if (present(level)) then
            ilevel = level
        else
            ilevel = EVENT_LEVEL_ERROR
        end if

        ! Call the super-class' constructor
        this%FExceptionDescription = &
            FExceptionDescription(message=message, &
                                  context=icontext, &
                                    level=ilevel)
    end function
!******************************************************************************!
!> Get a short description of a file error object.
!******************************************************************************!
    pure function fileError_string(this) result(string)
        class(FileError), intent(in) :: this
        character(:), allocatable :: string
!------
        if (this%line == 0) then
            string = this%filename // ": " // this%FExceptionDescription%string()
        else
            string = this%filename // ":" // this%line // ": " // this%FExceptionDescription%string()
        end if
    end function
end module
