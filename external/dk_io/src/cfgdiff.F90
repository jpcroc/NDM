!******************************************************************************!
!                             cfgdiff.f90
!------------------------------------------------------------------------------!
!> Compare the content of structure files supported by the dk_io library.
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
program cfgdiff
    use iso_fortran_env, only: ERROR_UNIT

    use dk_structure_io, only: read_structure, write_structure, TAG_LENGTH
    use dk_exception,    only: FException
    use dk_datatable,    only: DataTable
    use dk_datacolumn,   only: DataColumn
    use ext_character,   only: operator(//)

    implicit none(external, type)

    integer, parameter :: e = kind(1.0d0)
    character(:), allocatable :: input_file, output_file
    integer :: i
    character(1000) :: buffer
    real(e), dimension(3,3) :: box, box2
    real(e), dimension(:,:), allocatable :: positions, positions2
    real(e), dimension(:,:), allocatable :: velocities, velocities2
    real(e), dimension(:,:), allocatable :: forces, forces2
    real(e), dimension(:,:), allocatable :: data, data2
    character(TAG_LENGTH), dimension(:), allocatable :: tags, tags2
    logical :: auxiliary_mode, velocity_mode, force_mode
    type(DataTable) :: auxiliaries, auxiliaries2
    integer :: verbosity = 0

    i = 0
    auxiliary_mode = .false.
    velocity_mode = .false.
    force_mode = .false.
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
        case ("-a", "--auxiliaries")
            auxiliary_mode = .true.
        case ("-l", "--velocities")
            velocity_mode = .true.
        case ("-r", "--forces")
            force_mode = .true.
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
        write(ERROR_UNIT, '(a)') "Please provide two input files"
        stop 1, quiet=.true.
    end if

    if (.not.allocated(output_file)) then
        call print_help
        write(*,*)
        write(ERROR_UNIT, '(a)') "Please provide two input files"
        stop 1, quiet=.true.
    end if

!------ Read the input file
    if (verbosity>0) write(*,*) "Files to compare: "//input_file//", "//output_file
    call read_structure(input_file, box, positions, tags, forces=forces, velocities=velocities, aux=auxiliaries)
    call read_structure(output_file, box2, positions2, tags2, forces=forces2, velocities=velocities2, aux=auxiliaries2)

!------ Set the first data array
    if (velocity_mode .and. force_mode) then
        allocate(data(9,size(positions, 2)))
        data = 0
        data(1:3,:) = positions(:,:)
        if (allocated(velocities)) data(4:6,:) = velocities(:,:)
        if (allocated(forces)) data(7:9,:) = forces(:,:)

    else if (velocity_mode .and. .not.force_mode) then
        allocate(data(6,size(positions, 2)))
        data = 0
        data(1:3,:) = positions(:,:)
        if (allocated(velocities)) data(4:6,:) = velocities(:,:)

    else if (.not.velocity_mode .and. force_mode) then
        allocate(data(6,size(positions, 2)))
        data = 0
        data(1:3,:) = positions(:,:)
        if (allocated(forces)) data(4:6,:) = forces(:,:)

    else if (.not.velocity_mode .and. .not.force_mode) then
        allocate(data(3,size(positions, 2)))
        data = 0
        data(1:3,:) = positions(:,:)
    end if

!------ Set the second data array
    if (velocity_mode .and. force_mode) then
        allocate(data2(9,size(positions2, 2)))
        data2(1:3,:) = positions2(:,:)
        if (allocated(velocities2)) data2(4:6,:) = velocities2(:,:)
        if (allocated(forces2)) data2(7:9,:) = forces2(:,:)

    else if (velocity_mode .and. .not.force_mode) then
        allocate(data2(6,size(positions2, 2)))
        data2(1:3,:) = positions2(:,:)
        if (allocated(velocities2)) data2(4:6,:) = velocities2(:,:)

    else if (.not.velocity_mode .and. force_mode) then
        allocate(data2(6,size(positions2, 2)))
        data2(1:3,:) = positions2(:,:)
        if (allocated(forces2)) data2(4:6,:) = forces2(:,:)

    else if (.not.velocity_mode .and. .not.force_mode) then
        allocate(data2(3,size(positions2, 2)))
        data2(1:3,:) = positions2(:,:)
    end if

!------ Compare the data arrays
    call compare_structures

contains
    subroutine compare_structures

        integer :: i, j, a, n, differences

        i=1
        j=1
        differences = 0
        do while(i<=size(tags) .and. j<=size(tags2))

            ! Process the next pair of atoms
            if (all(data(:,i)==data2(:,j))) then
                i = i + 1
                j = j + 1
                cycle
            end if

            ! Look for an atom equivalent to i in the second structure
            a = j
            do while (any(data(:,i) /= data2(:,a)))
                a = a + 1
                if (a>size(data2, 2)) exit
            end do
            ! Here, either data(:,i) == data2(:,j), or we ran out of atoms in structure 2

            if (all(data(:,i) == data2(:,a))) then

                ! If we found one, then all atoms between j and a were added in the second structure
                do n=j, a-1
                    write(*,*) "+ ", n, tags2(n), data2(:,n)
                    differences = differences + 1
                end do
                j = a
                cycle
            else
                ! If we did not find one, then the atom i was deleted
                write(*,*) "- ", i, tags(i), data(:,i)
                differences = differences + 1
                i = i + 1
            end if
        end do

        if (i<=size(tags)) then
            do while (i<=size(tags))
                write(*,*) "- ", i, tags(i), data(:,i)
                differences = differences + 1
                i = i + 1
            end do
        end if

        if (i>size(tags)) then
            do n=j, size(tags2)
                write(*,*) "+ ", n, tags2(n), data2(:,n)
                differences = differences + 1
            end do
        end if

        if (differences > 0) stop differences, quiet = .true.
    end subroutine

    subroutine print_help
        write(*,'(a)') "Usage: cfgdiff [options] file1 file2"
        write(*,'(a)') ""
        write(*,'(a)') "Available options:"
        write(*,'(a)') "-a, --auxiliaries  consider auxiliary fields when comparing structures"
        write(*,'(a)') "-f, --forces     consider atomic forces when comparing structures"
        write(*,'(a)') "-h, --help       print this help message and exit"
        write(*,'(a)') "-i, --info       print detailed information about the program and the libraries it uses"
        write(*,'(a)') "-l, --velocities consider atomic velocities when comparing structures"
        write(*,'(a)') "-q, --quiet      suppress unnecessary output during processing"
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

        write(*,'(a)') "DynamicsKit:cfgdiff v." // PACKAGE_VERSION
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
