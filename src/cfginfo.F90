!******************************************************************************!
!                             cfginfo.f90
!------------------------------------------------------------------------------!
!> cfginfo program to print high-level informations about the atomistic
!> structures in a given file.
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
program cfginfo
    use dk_structure_io, only: read_structure, TAG_LENGTH, parameters_with_box
    use dk_exception,    only: FException
    use dk_datatable,    only: DataTable
    use dk_datacolumn,   only: DataColumn
    use ext_character,   only: operator(//)

    implicit none(external, type)

    integer, parameter :: e = kind(1.0d0)
    character(:), allocatable :: filename
    integer :: i
    character(1000) :: buffer
    real(e), dimension(3,3) :: box
    real(e), dimension(:,:), allocatable :: positions
    character(TAG_LENGTH), dimension(:), allocatable :: tags
    type(FException) :: stat
    logical :: show_auxiliaries, show_symmetry
    type(DataTable) :: auxiliaries
    type(DataColumn), pointer :: col
    character(:), allocatable :: format
    integer :: verbosity = 0
    real(e), dimension(6) :: parameters

    i = 0
    show_auxiliaries = .false.
    show_symmetry = .false.
    do while (i< command_argument_count())
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
        case ("-f", "--format")
            i = i + 1
            call get_command_argument(i, buffer)
            format = trim(buffer)
        case ("-a", "--show-auxiliaries")
            show_auxiliaries = .true.
        case ("-s", "--show-symmetry")
            show_symmetry = .true.
        case default
            filename = trim(buffer)
        end select
    end do

    if (.not.allocated(filename)) then
        call print_help
        write(*,*)
        error stop "Please provide an input file"
    end if

    write(*,*) "File name: ", filename

    block
        if (allocated(format)) then
            call read_structure(filename, box, positions, tags, aux=auxiliaries, format=format)
        else
            call read_structure(filename, box, positions, tags, aux=auxiliaries)
        end if
    end block

    if (stat /= 0) then
        call stat%report
    end if

    write(*,*) "Cell shape:"
    write(*,*) box(:,1)
    write(*,*) box(:,2)
    write(*,*) box(:,3)

    parameters = parameters_with_box(box)
    write(*,*)
    write(*,*) "Cell parameters:"
    write(*,*) "a", parameters(1)
    write(*,*) "b", parameters(2)
    write(*,*) "c", parameters(3)
    write(*,*) "α", parameters(4)
    write(*,*) "β", parameters(5)
    write(*,*) "γ", parameters(6)

    write(*,*)
    write(*,*) "Particles:", size(tags)

    if (show_auxiliaries) then
        write(*,*)
        write(*,*) "Auxiliary fields:"
        do i=1, auxiliaries%num_columns()
            col => auxiliaries%column(i)
            if (col%name() /= "") call col%print_summary
        end do
    end if

    if (show_symmetry) then
        write(*,*)
#ifdef HAVE_SPGLIB
        block
            use spglib_f08

            type(SpglibDataset) :: dset
            integer, dimension(:), allocatable :: types
            integer :: i, j, n
            character(:), allocatable :: string

            allocate(types(size(tags)))
            types = 0
            n = 0
            do i=1, size(tags)
                if (types(i) > 0) cycle
                n = n + 1
                types(i) = n
                do j=i+1, size(tags)
                    if (tags(i)==tags(j)) types(j) = n
                end do
            end do

            write(*,*) "Crystal structure:"
            dset = SpglibDataset(box, positions, types, size(positions, 2), 1.0e-2_e)
            write(*,*) " Space group: ", dset%international_symbol // " (" // dset%spacegroup_number // ")"
            write(*,*) " Wyckoff sites: "
            do i=minval(dset%wyckoffs), maxval(dset%wyckoffs)
                if (count(dset%wyckoffs==i) == 0) cycle
                write(*,*) count(dset%wyckoffs==i) // achar(97+i)
                string = ""
                n = 0
                do j=1, size(dset%wyckoffs)
                    if (dset%wyckoffs(j)/=i) cycle
                    string = string // "  (" // &
                            coordinate_string(dset%std_positions(1,j)) // " , " // &
                            coordinate_string(dset%std_positions(2,j))  // " , " // &
                            coordinate_string(dset%std_positions(3,j)) // ")"
                    n = n + 1
                    if (n == 3) then
                        write(*,*) string
                        string = ""
                        n = 0
                    end if
!                     if (dset%wyckoffs(j)==i) then
!                         write(*,*) "  (" // &
!                             coordinate_string(dset%std_positions(1,j)) // " , " // &
!                             coordinate_string(dset%std_positions(2,j))  // " , " // &
!                             coordinate_string(dset%std_positions(3,j)) // ")"
!                     end if
                end do

                if (string /= "") write(*,*) string
            end do
        end block
#else
        write(*,*) "The symmetry cannot be shown because support for the spglib library was not enabled."
#endif
    end if
contains
    function coordinate_string(x) result(string)
        real(e), intent(in) :: x
        character(:), allocatable :: string

        if (x == 0.25) then
            string = "1/4"
        else if (x == 0.5) then
            string = "1/2"
        else if (x == 0.75) then
            string = "3/4"
        else
            string = "" // x
        end if
    end function

    subroutine print_help
        write(*,*) "Usage: cfginfo [options] file"
        write(*,'(a)') ""
        write(*,'(a)') "Available options:"
        write(*,'(a)') "-a, --show-auxiliaries"
        write(*,'(a)') "-f, --format"
        write(*,'(a)') "-h, --help       print this help message and exit"
        write(*,'(a)') "-i, --info       print detailed information about the program and the libraries it uses"
        write(*,'(a)') "-q, --quiet      suppress unnecessary output during processing"
        write(*,'(a)') "-s, --show-symmetry"
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

        write(*,'(a)') "DynamicsKit:cfginfo v." // PACKAGE_VERSION
        write(*,'(a)') "Copyright (C) 2025 CEA"

        write(*,*)
        write(*,'(a)') "Supported libraries:"
#ifdef HAVE_LIBZ
        c_string = zlibVersion()
        call c_f_string(c_string, string)
        write(*,'(a)') " zlib v. " // string
#else
        write(*,'(a)') " zlib support disabled"
#endif
#ifdef HAVE_LIBBZ2
        c_string = bzlibVersion()
        call c_f_string(c_string, string)
        write(*,'(a)') " libbz2 v. " // string
#else
        write(*,'(a)') " libbz2 support disabled"
#endif
#ifdef HAVE_LIBZSTD
        c_string = ZSTD_versionString()
        call c_f_string(c_string, string)
        write(*,'(a)') " libbzstd v. " // string
#else
        write(*,'(a)') " libbzstd support disabled"
#endif
#ifdef HAVE_SPGLIB
        write(*,'(a)') " spglib v. " // trim(spg_get_version())
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
