!******************************************************************************!
!                             cfgconvert.f90
!------------------------------------------------------------------------------!
!> Convert structure files between the formats supported by the dk_io library.
!
!  Written by Paul Fossati, <paul.fossati@cea.fr>
!  Copyright (c) 2025 CEA
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
program cfgconvert
    use iso_fortran_env, only: ERROR_UNIT

    use dk_structure_io, only: read_structure, write_structure, TAG_LENGTH
    use dk_exception,    only: FException
    use dk_datatable,    only: DataTable
    use dk_datacolumn,   only: DataColumn

    implicit none(external, type)

    integer, parameter :: e = kind(1.0d0)
    character(:), allocatable :: input_file, output_file
    integer :: i
    character(1000) :: buffer
    real(e), dimension(3,3) :: box
    real(e), dimension(:,:), allocatable :: positions
    character(TAG_LENGTH), dimension(:), allocatable :: tags
    type(FException) :: stat
    logical :: show_auxiliaries
    type(DataTable) :: auxiliaries
    character(:), allocatable :: input_format, output_format
    integer :: verbosity = 0

    i = 0
    show_auxiliaries = .false.
    do while (i < command_argument_count())
        i = i + 1
        call get_command_argument(i, buffer)

        select case(buffer)
        case ("-h", "--help")
            call print_help
            stop
        case ("-v", "--verbose")
            verbosity = verbosity + 1
        case ("-q", "--quiet")
            verbosity = verbosity - 1
        case ("-i", "--info")
            call print_info
            stop
        case ("-f", "--input-format")
            if (i == command_argument_count()) then
                call print_help
                error stop
            end if
            i = i + 1
            call get_command_argument(i, buffer)
            input_format = trim(buffer)
        case ("-t", "--output-format")
            if (i == command_argument_count()) then
                call print_help
                error stop
            end if
            i = i + 1
            call get_command_argument(i, buffer)
            output_format = trim(buffer)
        case default
            if (allocated(input_file)) then
                output_file = trim(buffer)
            else
                input_file = trim(buffer)
            end if
        end select
    end do

    if (.not.allocated(input_file)) then
        call print_help
        write(*,*)
        write(ERROR_UNIT, '(a)') "Please provide an input file"
        stop 1, quiet=.true.
    end if

    if (.not.allocated(output_file)) then
        call print_help
        write(*,*)
        write(ERROR_UNIT, '(a)') "Please provide an output file"
        stop 1, quiet=.true.
    end if

!------ Read the input file
    if (verbosity>0) write(*,*) "Input file: "//input_file
    if (allocated(input_format)) then
        if (verbosity>0) write(*,*) "Input file format: "//input_format
        call read_structure(input_file, box, positions, tags, aux=auxiliaries, format=input_format)
    else
        call read_structure(input_file, box, positions, tags, aux=auxiliaries)
    end if

!------ Write the output file
    if (verbosity>0) write(*,*) "Output file: "//output_file
    if (allocated(output_format)) then
        if (verbosity>0) write(*,*) "Input file format: "//output_format
        call write_structure(output_file, box, positions, tags, aux=auxiliaries, format=output_format, stat=stat)
    else
        call write_structure(output_file, box, positions, tags, aux=auxiliaries, stat=stat)
    end if

    if (stat /= 0) then
        call stat%report
    end if

!------ Cleanup
    if (allocated(tags)) deallocate(tags)
    if (allocated(positions)) deallocate(positions)
    if (allocated(input_file)) deallocate(input_file)
    if (allocated(output_file)) deallocate(output_file)
    if (allocated(input_format)) deallocate(input_format)
    if (allocated(output_format)) deallocate(output_format)
contains
    subroutine print_help
        write(*,'(a)') "Usage: cfgconvert [options] input output"
        write(*,'(a)') ""
        write(*,'(a)') "Available options:"
        write(*,'(a)') "-f FORMAT, --input-format FORMAT   specify the format of the input file"
        write(*,'(a)') "-h, --help       print this help message and exit"
        write(*,'(a)') "-i, --info       print detailed information about the program and the libraries it uses"
        write(*,'(a)') "-q, --quiet      suppress unnecessary output during processing"
        write(*,'(a)') "-t FORMAT, --output-format FORMAT  specify the format of the output file"
        write(*,'(a)') "-v, --verbose    print additional output during processing"
    end subroutine

    subroutine print_info
#ifdef HAVE_LIBZ
        use zlib, only: zlibVersion
#endif
#ifdef HAVE_LIBBZ2
        use bzlib, only: bzlibVersion
#endif
#ifdef HAVE_SPGLIB
        use spglib_f08, only: spg_get_version
#endif
        use iso_c_binding, only: c_ptr
        use iso_fortran_env, only: compiler_version, compiler_options
        use ext_character

        interface
            function ZSTD_versionString() bind(C, name="ZSTD_versionString")
                use iso_c_binding, only: c_ptr
                type(c_ptr) :: ZSTD_versionString
            end function
        end interface

        type(c_ptr) :: c_string
        character(:), allocatable :: string

        write(*,'(a)') "DynamicsKit:cfgconvert v." // PACKAGE_VERSION
        write(*,'(a)') "Copyright (C) 2025 CEA"

        write(*,*)
        write(*,'(a)') "Enabled libraries:"
#ifdef HAVE_LIBZ
        c_string = zlibVersion()
        call c_f_string(c_string, string)
        write(*,'(a)') " zlib v." // string
#else
        write(*,'(a)') " zlib support disabled"
#endif
#ifdef HAVE_LIBBZ2
        c_string = bzlibVersion()
        call c_f_string(c_string, string)
        write(*,'(a)') " libbz2 v." // string
#else
        write(*,'(a)') " libbz2 support disabled"
#endif
#ifdef HAVE_LIBZSTD
        c_string = ZSTD_versionString()
        call c_f_string(c_string, string)
        write(*,'(a)') " libbzstd v." // string
#else
        write(*,'(a)') " libbzstd support disabled"
#endif
#ifdef HAVE_SPGLIB
        write(*,'(a)') " spglib v." // trim(spg_get_version())
#else
        write(*,'(a)') " spglib support disabled"
#endif
        if (verbosity > 0) then
            write(*,*)
            write(*,'(a)') "Compiled with "//trim(compiler_version())
            write(*,'(a)') "Build options: "//trim(compiler_options())
        end if
    end subroutine
end program
