!******************************************************************************!
!                        dk_structure_io module
!------------------------------------------------------------------------------!
!> Procedures to read and write structure files. The subroutines in this module
!> are high-level and call specific subroutines depending on the type of
!> structure file.
!>
!> ## Supported file formats
!>
!> This table summarises the different features supported when reading the supported file formats.
!>
!>  Format  | Positions | Velocities | Forces | Spin    | Masses | Charges | Auxiliaries | Multiple frames | Symmetry
!> ---------|-----------|------------|--------|---------|--------|---------|-------------|-----------------|----------
!>  Abinit  | 🕒        | ❌         | ❌     | 🕒     | ❌    | ❌      | ❌         |❌              | ❌
!>  Atomeye | ✔        | ✔          | ❌ [0] | ❌ [0] | 🕒    | ❌      | ✔          |❌              | ❌
!>  CASTEP  | ✔        | ✔          | ❌     | 🕒     | ❌    | ❌      | ❌          |❌              | 🕒 [1]
!>  CIF     | ✔        | ❌         | ❌     | ❌     | ❌    | ❌      | ❌          |🕒              | ✔
!>  DL_POLY | ✔        | ✔          | ✔     | ❌     | ❌    | ❌      | ❌          |🕒              | ❌
!>  GULP    | ✔        | ❌         | ❌     | ❌     | ❌    | ❌      | ❌          |🕒              | ✔
!>  LAMMPS  | ✔        | 🕒         | ❌     | ❌     | ❌    | ❌      | ❌          |❌              | ❌
!>  VASP    | ✔        | 🕒         | ❌     | 🕒     | ❌    | ❌      | ❌          |❌              | ❌
!>  XYZ     | ✔        | 🕒         | 🕒 [0] | 🕒 [0] | 🕒    | 🕒      | ✔          |🕒              | ❌
!>
!> [0]: Forces and magnetic moments can be read from the auxiliary fields.
!>
!> [1]: The full structure can be build from the irreductible crystal sites if the SYMMETRY_OPS block is present.
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
module dk_structure_io
    use iso_c_binding, only: long=>c_long

    use iso_c_binding, only: c_double, C_NULL_CHAR, c_size_t, c_char

    use ext_character, only: lower_case

    use dk_datatable,  only: DataTable
    use dk_exception,  only: FException
    use dk_structure,  only: Structure, TAG_LENGTH
    use dk_parameters, only: e, pi, qp
    use dk_fileobject, only: FileObject, FileError

    implicit none(external, type)
    private

    type, public :: Element
        character(:), allocatable :: tag
        real(e) :: mass = 0.0_e
        real(e) :: charge = 0.0_e
    end type

    character(*), parameter, public :: UNDEF = "UNDEF"

    character(5), parameter :: COMPRESSION_NONE  = "ascii"
    character(5), parameter :: COMPRESSION_GZIP  = "gzip"
    character(5), parameter :: COMPRESSION_BZIP2 = "bzip2"
    character(5), parameter :: COMPRESSION_ZIP   = "zip"
    character(5), parameter :: COMPRESSION_ZSTD  = "zstd"

    character(5), dimension(*), parameter, public :: COMPRESSION_FORMATS = [ &
        COMPRESSION_NONE, &
        COMPRESSION_GZIP, &
        COMPRESSION_BZIP2, &
        COMPRESSION_ZIP, &
        COMPRESSION_ZSTD &
    ]

    character(10), parameter :: STRUCTURE_UNKNOWN    = "unknown"
    character(10), parameter :: STRUCTURE_ABINIT     = "abinit"
    character(10), parameter :: STRUCTURE_CASTEP     = "castep"
    character(10), parameter :: STRUCTURE_CASTEP_OUT = "castep-out"
    character(10), parameter :: STRUCTURE_CIF        = "cif"
    character(10), parameter :: STRUCTURE_CFG        = "cfg"
    character(10), parameter :: STRUCTURE_DLPOLY     = "dlpoly"
    character(10), parameter :: STRUCTURE_GULP       = "gulp"
    character(10), parameter :: STRUCTURE_GULP_OUT   = "gulp-out"
    character(10), parameter :: STRUCTURE_LAMMPS     = "lammps"
    character(10), parameter :: STRUCTURE_VASP       = "vasp"
    character(10), parameter :: STRUCTURE_XFG        = "xfg"
    character(10), parameter :: STRUCTURE_XSF        = "xsf"
    character(10), parameter :: STRUCTURE_XYZ        = "xyz"

    public :: TAG_LENGTH

    character(10), dimension(*), parameter, public :: STRUCTURE_FORMATS = [ &
        STRUCTURE_ABINIT, &
        STRUCTURE_CASTEP, &
        STRUCTURE_CASTEP_OUT, &
        STRUCTURE_CIF, &
        STRUCTURE_CFG, &
        STRUCTURE_DLPOLY, &
        STRUCTURE_GULP, &
        STRUCTURE_GULP_OUT, &
        STRUCTURE_LAMMPS, &
        STRUCTURE_VASP, &
        STRUCTURE_XFG, &
        STRUCTURE_XYZ &
    ]

    type, public :: FileInfo
        character(:), allocatable :: data_format
        character(:), allocatable :: compression_format
        integer :: num_frames = 0
    end type

    public :: read_structure
    public :: write_structure
    public :: count_frames_in_file
    public :: get_file_format
    public :: get_test_structure
    public :: get_triclinic_test_structure

    public :: parameters_with_box
    public :: box_with_parameters
    public :: get_mark

    interface
        module subroutine read_abinit(input_file, box, positions, tags, info, stat)
            type(FileObject), intent(inout) :: input_file
            real(e), dimension(3,3), intent(out) :: box
            real(e), dimension(:,:), allocatable, intent(out) :: positions
            character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
            type(FileInfo), intent(out), optional :: info
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine write_abinit(file, box, positions, tags, masses, include, elements, stat)
            character(*), intent(in) :: file
            real(e), dimension(3,3), intent(in) :: box
            real(e), dimension(:,:), intent(in) :: positions
            character(TAG_LENGTH), dimension(:), intent(in) :: tags
            real(e), dimension(:), intent(in), optional :: masses
            logical, dimension(:), intent(in), optional :: include
            type(Element), dimension(:), intent(in), optional :: elements
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine read_castep(input_file, box, positions, tags, velocities, info, stat)
            type(FileObject), intent(inout) :: input_file
            real(e), dimension(3,3), intent(out) :: box
            real(e), dimension(:,:), allocatable, intent(out) :: positions
            character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
            real(e), dimension(:,:), allocatable, intent(out), optional :: velocities
            type(FileInfo), intent(out), optional :: info
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine write_castep(filename, box, positions, tags, masses, velocities, include, elements, stat)
            character(*), intent(in) :: filename
            real(e), dimension(3,3), intent(in) :: box
            real(e), dimension(:,:), intent(in) :: positions
            character(TAG_LENGTH), dimension(:), intent(in) :: tags
            real(e), dimension(:), intent(in), optional :: masses
            real(e), dimension(:,:), intent(in), optional :: velocities
            logical, dimension(:), intent(in), optional :: include
            type(Element), dimension(:), intent(in), optional :: elements
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine read_cif(input_file, box, positions, tags, info, stat)
            type(FileObject), intent(inout) :: input_file
            real(e), dimension(3,3), intent(out) :: box
            real(e), dimension(:,:), allocatable, intent(out) :: positions
            character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
            type(FileInfo), intent(out), optional :: info
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine write_cif(filename, box, positions, tags, include, scale, stat)
            character(*), intent(in) :: filename
            real(e), dimension(3,3), intent(in) :: box
            real(e), dimension(:,:), intent(in) :: positions
            character(TAG_LENGTH), dimension(:), intent(in) :: tags
            logical, dimension(:), intent(in), optional :: include
            real(e), intent(in), optional :: scale
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine read_dlpoly(input_file, box, positions, tags, velocities, forces, info, stat)
            type(FileObject), intent(inout) :: input_file
            real(e), dimension(3,3), intent(out) :: box
            real(e), dimension(:,:), allocatable, intent(out) :: positions
            character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
            real(e), dimension(:,:), allocatable, intent(out), optional :: velocities
            real(e), dimension(:,:), allocatable, intent(out), optional :: forces
            type(FileInfo), intent(out), optional :: info
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine write_dlpoly(filename, box, positions, tags, include, elements, velocities, forces, stat)
            character(*), intent(in) :: filename
            real(e), dimension(3,3), intent(in) :: box
            real(e), dimension(:,:), intent(in) :: positions
            character(TAG_LENGTH), dimension(:), intent(in) :: tags
            logical, dimension(:), intent(in), optional :: include
            type(Element), dimension(:), intent(in), optional :: elements
            real(e), dimension(:,:), intent(in), optional :: velocities
            real(e), dimension(:,:), intent(in), optional :: forces
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine read_gulp(input_file, box, positions, tags, info, stat)
            type(FileObject), intent(inout) :: input_file
            real(e), dimension(3,3), intent(out) :: box
            real(e), dimension(:,:), allocatable, intent(out) :: positions
            character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
            type(FileInfo), intent(out), optional :: info
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine write_gulp(filename, box, positions, tags, masses, include, elements, stat)
            character(*), intent(in) :: filename
            real(e), dimension(3,3), intent(in) :: box
            real(e), dimension(:,:), intent(in) :: positions
            character(TAG_LENGTH), dimension(:), intent(in) :: tags
            real(e), dimension(:), intent(in), optional :: masses
            logical, dimension(:), intent(in), optional :: include
            type(Element), dimension(:), intent(in), optional :: elements
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine read_xfg(input_file, box, positions, tags, masses, velocities, forces, aux, info, stat)
            type(FileObject), intent(inout) :: input_file
            real(e), dimension(3,3), intent(out) :: box
            real(e), dimension(:,:), allocatable, intent(out) :: positions
            character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
            real(e), dimension(:), allocatable, intent(out), optional :: masses
            real(e), dimension(:,:), allocatable, intent(out), optional :: velocities
            real(e), dimension(:,:), allocatable, intent(out), optional :: forces
            type(DataTable), intent(out), optional :: aux
            type(FileInfo), intent(out), optional :: info
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine write_xfg(filename, box, positions, tags, masses, velocities, forces, aux, include, scale, stat)
            character(*), intent(in) :: filename
            real(e), dimension(3,3), intent(in) :: box
            real(e), dimension(:,:), intent(in) :: positions
            character(TAG_LENGTH), dimension(:), intent(in) :: tags
            real(e), dimension(:), intent(in), optional :: masses
            real(e), dimension(:,:), intent(in), optional :: velocities
            real(e), dimension(:,:), intent(in), optional :: forces
            type(DataTable), intent(in), optional :: aux
            logical, dimension(:), intent(in), optional :: include
            real(e), intent(in), optional :: scale
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine read_xyz(input_file, box, positions, tags, masses, velocities, forces, aux, info, stat)
            type(FileObject), intent(inout) :: input_file
            real(e), dimension(3,3), intent(out) :: box
            real(e), dimension(:,:), allocatable, intent(out) :: positions
            character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
            real(e), dimension(:), allocatable, intent(out), optional :: masses
            real(e), dimension(:,:), allocatable, intent(out), optional :: velocities
            real(e), dimension(:,:), allocatable, intent(out), optional :: forces
            type(DataTable), intent(out), optional :: aux
            type(FileInfo), intent(out), optional :: info
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine write_xyz(filename, box, positions, tags, masses, velocities, forces, aux, include, scale, stat)
            character(*), intent(in) :: filename
            real(e), dimension(3,3), intent(in) :: box
            real(e), dimension(:,:), intent(in) :: positions
            character(TAG_LENGTH), dimension(:), intent(in) :: tags
            real(e), dimension(:), intent(in), optional :: masses
            real(e), dimension(:,:), intent(in), optional :: velocities
            real(e), dimension(:,:), intent(in), optional :: forces
            type(DataTable), intent(in), optional :: aux
            logical, dimension(:), intent(in), optional :: include
            real(e), intent(in), optional :: scale
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine read_lammps(input_file, box, positions, tags, velocities, info, stat)
            type(FileObject), intent(inout) :: input_file
            real(e), dimension(3,3), intent(out) :: box
            real(e), dimension(:,:), allocatable, intent(out) :: positions
            character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
            real(e), dimension(:,:), allocatable, intent(out), optional :: velocities
            type(FileInfo), intent(out), optional :: info
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine write_lammps(file, box, positions, tags, masses, velocities, forces, charges, aux, include, elements, stat)
            character(*), intent(in) :: file
            real(e), dimension(3,3), intent(in) :: box
            real(e), dimension(:,:), intent(in) :: positions
            character(TAG_LENGTH), dimension(:), intent(in) :: tags
            real(e), dimension(:), intent(in), optional :: masses
            real(e), dimension(:,:), intent(in), optional :: velocities
            real(e), dimension(:,:), intent(in), optional :: forces
            real(e), dimension(:), intent(in), optional :: charges
            type(DataTable), intent(in), optional :: aux
            logical, dimension(:), intent(in), optional :: include
            type(Element), dimension(:), intent(in), optional :: elements
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine read_vasp(input_file, box, positions, tags, velocities, info, stat)
            type(FileObject), intent(inout) :: input_file
            real(e), dimension(3,3), intent(out) :: box
            real(e), dimension(:,:), allocatable, intent(out) :: positions
            character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
            real(e), dimension(:,:), allocatable, intent(out), optional :: velocities
            type(FileInfo), intent(out), optional :: info
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine write_vasp(filename, box, positions, tags, masses, include, elements, stat)
            character(*), intent(in) :: filename
            real(e), dimension(3,3), intent(in) :: box
            real(e), dimension(:,:), intent(in) :: positions
            character(TAG_LENGTH), dimension(:), intent(in) :: tags
            real(e), dimension(:), intent(in), optional :: masses
            logical, dimension(:), intent(in), optional :: include
            type(Element), dimension(:), intent(in), optional :: elements
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine read_xsf(input_file, box, positions, tags, velocities, info, stat)
            type(FileObject), intent(inout) :: input_file
            real(e), dimension(3,3), intent(out) :: box
            real(e), dimension(:,:), allocatable, intent(out) :: positions
            character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
            real(e), dimension(:,:), allocatable, intent(out), optional :: velocities
            type(FileInfo), intent(out), optional :: info
            type(FException), intent(out), optional :: stat
        end subroutine

        module subroutine write_xsf(filename, box, positions, tags, masses, include, elements, stat)
            character(*), intent(in) :: filename
            real(e), dimension(3,3), intent(in) :: box
            real(e), dimension(:,:), intent(in) :: positions
            character(TAG_LENGTH), dimension(:), intent(in) :: tags
            real(e), dimension(:), intent(in), optional :: masses
            logical, dimension(:), intent(in), optional :: include
            type(Element), dimension(:), intent(in), optional :: elements
            type(FException), intent(out), optional :: stat
        end subroutine
    end interface

    interface
        function dk_read_array(str, x, n, pos) bind(C, name="dk_read_array")
            use iso_c_binding, only: c_int, c_double, c_size_t, c_char
            import
            integer(c_int) :: dk_read_array
            character(c_char), dimension(*), intent(inout) :: str
            real(c_double), dimension(*), intent(out) :: x
            integer(c_size_t), value :: n
            integer(c_size_t), intent(out) :: pos
        end function
        function writeArray(str, x, n, buffsize) bind(C, name="dk_write_array")
            use iso_c_binding, only: c_int, c_double, c_size_t, c_char
            import
            integer(c_int) :: writeArray
            character(len=1,kind=c_char), dimension(*), intent(inout) :: str
            real(c_double), dimension(*), intent(in) :: x
            integer(c_size_t), value :: n
            integer(c_int), value :: buffsize
        end function
    end interface

    public :: read_array
contains
!******************************************************************************!
!> Read numerical values from a character variable.
!******************************************************************************!
    subroutine read_array(line, array, n, pos)
        character(*), intent(inout) :: line !< Character variable to parse
        real(e), dimension(:), intent(inout) :: array !< Data read from the string
        integer, intent(out) :: n !< Number of elements that could be read
        integer, intent(out), optional :: pos
!------
        integer :: l
        character(:,c_char), allocatable :: tmp
        integer(c_size_t) :: pos_
!------
        ! Make sure the string is null-terminated
        l = len_trim(line)

        ! Return if the line is empty
        if (l==0) return

        ! Make sure the line is NULL-terminated
        if (line(l:l) == C_NULL_CHAR) then
            n = dk_read_array(line, array, size(array, kind=c_double), pos_)
        else if (len(line) > len_trim(line)) then
            line(l+1:l+1) = C_NULL_CHAR
            n = dk_read_array(line, array, size(array, kind=c_double), pos_)
        else
            allocate(character(len(line)+1,c_char) :: tmp)
            tmp(:len(line)) = line
            tmp(len(tmp):) = C_NULL_CHAR
            n = dk_read_array(tmp, array, size(array, kind=c_double), pos_)
        end if

        if (present(pos)) then
            pos = pos_
        end if
    end subroutine
!******************************************************************************!
!> Get lattice parameters from a box matrix.
!******************************************************************************!
    pure function parameters_with_box(box) result(p)
        real(e), dimension(3,3), intent(in) :: box !< Box matrix
        real(e), dimension(6) :: p !< Parameters calculated from the box
!------
        p(1) = sqrt(dot_product(box(1,:), box(1,:)))
        p(2) = sqrt(dot_product(box(2,:), box(2,:)))
        p(3) = sqrt(dot_product(box(3,:), box(3,:)))
        p(4) = acos(dot_product(box(2,:), box(3,:)) / (p(2)*p(3))) * 180 / PI
        p(5) = acos(dot_product(box(1,:), box(3,:)) / (p(1)*p(3))) * 180 / PI
        p(6) = acos(dot_product(box(1,:), box(2,:)) / (p(1)*p(2))) * 180 / PI
    end function
!******************************************************************************!
!> Get a box matrix from lattice parameters.
!******************************************************************************!
    pure function box_with_parameters(parameters) result(box)
        real(e), dimension(6), intent(in) :: parameters !< Parameters
        real(e), dimension(3,3) :: box !< Box calculated from the parameters
!------
        real(e) :: a, b, c, alpha, beta, gamma
!------
        a = parameters(1)
        b = parameters(2)
        c = parameters(3)
        alpha = parameters(4)! * PI / 180
        beta = parameters(5)! * PI / 180
        gamma = parameters(6)! * PI / 180

        box = 0

        ! First vector
        box(1,1) = a

        ! Second vector
        box(2,1) = b*cosd(gamma)
        box(2,2) = sqrt(b**2-box(2,1)**2)

        ! Third vector
        box(3,1) = c*cosd(beta)
        box(3,2) = (b*c*cosd(alpha) - box(2,1)*box(3,1)) / box(2,2)
        box(3,3) = sqrt(c**2-box(3,1)**2-box(3,2)**2)
    end function
!******************************************************************************!
!> Get a comment line to print in output files.
!******************************************************************************!
    function get_mark()
        use ext_character
        use dk_io_config, only: DK_IO_VERSION
        character(:), allocatable :: get_mark
        integer, dimension(4) :: time
        integer, dimension(3) :: date
        integer, dimension(8) :: values

        call date_and_time(values=values)
        date(1) = values(3)
        date(2) = values(2)
        date(3) = values(1)
        time(1) = values(5)
        time(2) = values(6)
        get_mark = "#Generated by DynamicsKit:dk_io v. " // trim(DK_IO_VERSION) // " - " // &
            date(1) // "/" // date(2) // "/" // date(3) // " at " // time(1) // ":" // time(2)
    end function
!******************************************************************************!
!> Get a simple structure for testing purposes.
!******************************************************************************!
    subroutine get_test_structure(box, positions, tags, velocities, forces, masses, aux)
        real(e), dimension(3,3), intent(out) :: box !< Simulation box matrix
        real(e), dimension(:,:), allocatable, intent(out) :: positions !< List of atomic positions
        character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags !< List of atomic tags
        real(e), dimension(:,:), allocatable, intent(out), optional :: velocities !< List of atomic velocities
        real(e), dimension(:,:), allocatable, intent(out), optional :: forces !< List of atomic forces
        real(e), dimension(:),   allocatable, intent(out), optional :: masses !< List of atomic masses
        type(DataTable), intent(out), optional :: aux !< Table with auxiliary values
!------
        integer :: i
!------
        ! Set the structure manually
        box = 0
        box(1,1) = 5
        box(2,2) = 5
        box(3,3) = 5

        ! Set the positions
        allocate(positions(3,12))
        positions(:,1)  = [0.0_e, 0.0_e, 0.0_e]
        positions(:,2)  = [0.0_e, 0.5_e, 0.5_e]
        positions(:,3)  = [0.5_e, 0.0_e, 0.5_e]
        positions(:,4)  = [0.5_e, 0.5_e, 0.0_e]
        positions(:,5)  = [0.25_e, 0.25_e, 0.25_e]
        positions(:,6)  = [0.75_e, 0.25_e, 0.25_e]
        positions(:,7)  = [0.25_e, 0.75_e, 0.25_e]
        positions(:,8)  = [0.25_e, 0.25_e, 0.75_e]
        positions(:,9)  = [0.25_e, 0.75_e, 0.75_e]
        positions(:,10) = [0.75_e, 0.25_e, 0.75_e]
        positions(:,11) = [0.75_e, 0.75_e, 0.25_e]
        positions(:,12) = [0.75_e, 0.75_e, 0.75_e]

        ! Set the tags
        allocate(tags(12))
        tags(1:4) = "Ca"
        tags(5:12) = "F"

        ! Set the velocities
        if (present(velocities)) then
            allocate(velocities(3,12))
            do i=1, 12
                velocities(1,i) = sin(i*0.2d0)
                velocities(2,i) = cos(i*0.2d0)
                velocities(3,i) = tan(i*0.2d0)
            end do
        end if

        ! Set the forces
        if (present(forces)) then
            allocate(forces(3,12))
            do i=1, 12
                forces(1,i) = sinh(i*0.2d0)
                forces(2,i) = cosh(i*0.2d0)
                forces(3,i) = tanh(i*0.2d0)
            end do
        end if

        ! Set the masses
        if (present(masses)) then
            allocate(masses(12))
            masses(1:4) = 1.0d0
            masses(5:12) = 2.0d0
        end if

        ! Set the auxiliary fields
        if (present(aux)) then
            call aux%add_column("Auxiliary")
            do i=1, 12
                call aux%append("Auxiliary", sin(i*0.2d0))
            end do
        end if
    end subroutine
!******************************************************************************!
!> Get a triclinic structure for testing purposes.
!******************************************************************************!
    subroutine get_triclinic_test_structure(box, positions, tags, velocities, forces, masses, aux)
        real(e), dimension(3,3), intent(out) :: box !< Simulation box matrix
        real(e), dimension(:,:), allocatable, intent(out) :: positions !< List of atomic positions
        character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags !< List of atomic tags
        real(e), dimension(:,:), allocatable, intent(out), optional :: velocities !< List of atomic velocities
        real(e), dimension(:,:), allocatable, intent(out), optional :: forces !< List of atomic forces
        real(e), dimension(:),   allocatable, intent(out), optional :: masses !< List of atomic masses
        type(DataTable), intent(out), optional :: aux !< Table with auxiliary values
!------
        integer :: i
!------

        ! Set the cell shape
        box = box_with_parameters([5.860_e, 5.104_e, 5.836_e, 93.82_e, 107.31_e, 115.81_e])

        ! Set the positions
        allocate(positions(3,10))
        positions(:,1)  = [2.172e-1_e, 2.158e-1_e, 7.256e-1_e]
        positions(:,2)  = [7.828e-1_e, 7.842e-1_e, 2.744e-1_e]
        positions(:,3)  = [4.873e-1_e, 9.385e-1_e, 3.037e-1_e]
        positions(:,4)  = [5.127e-1_e, 6.150e-2_e, 6.963e-1_e]
        positions(:,5)  = [4.530e-2_e, 6.001e-1_e, 1.957e-1_e]
        positions(:,6)  = [9.547e-1_e, 3.100e-1_e, 8.043e-1_e]
        positions(:,7)  = [4.198e-1_e, 2.990e-1_e, 1.522e-1_e]
        positions(:,8)  = [5.802e-1_e, 7.010e-1_e, 8.478e-1_e]
        positions(:,9)  = [1.193e-1_e, 2.469e-1_e, 3.303e-1_e]
        positions(:,10) = [8.807e-1_e, 7.531e-1_e, 6.697e-1_e]

        ! Set the tags
        allocate(tags(10))
        tags(1:2) = "Mn"
        tags(3:10) = "P"

        ! Set the velocities
        if (present(velocities)) then
            allocate(velocities, mold=positions)
            do i=1, size(velocities, 2)
                velocities(1,i) = sin(i*0.2d0)
                velocities(2,i) = cos(i*0.2d0)
                velocities(3,i) = tan(i*0.2d0)
            end do
        end if

        ! Set the forces
        if (present(forces)) then
            allocate(forces, mold=positions)
            do i=1, size(forces, 2)
                forces(1,i) = sinh(i*0.2d0)
                forces(2,i) = cosh(i*0.2d0)
                forces(3,i) = tanh(i*0.2d0)
            end do
        end if

        ! Set the masses
        if (present(masses)) then
            allocate(masses(size(positions, 2)))
            masses(1:4) = 1.0d0
            masses(5:12) = 2.0d0
        end if

        ! Set the auxiliary fields
        if (present(aux)) then
            call aux%add_column("Auxiliary")
            do i=1, size(positions, 2)
                call aux%append("Auxiliary", sin(i*0.2d0))
            end do
        end if
    end subroutine
!******************************************************************************!
!> Read a structure from a file.
!******************************************************************************!
    subroutine read_structure(filename, box, positions, tags, format, compression, options, &
        masses, velocities, forces, aux, file_info, stat, msg, except)

        character(*), intent(in) :: filename
        real(e), dimension(3,3), intent(out) :: box
        real(e), dimension(:,:), allocatable, intent(out) :: positions
        character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
        character(*), intent(in), optional :: format
        character(*), intent(in), optional :: compression
        real(e), dimension(:), allocatable, intent(out), optional :: masses
        real(e), dimension(:,:), allocatable, intent(out), optional :: velocities
        real(e), dimension(:,:), allocatable, intent(out), optional :: forces
        type(DataTable), intent(out), optional :: aux
        integer, intent(in), optional :: options
        type(FileInfo), intent(out), optional :: file_info
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: msg
        type(FException), intent(out), optional :: except
!------
        type(FException) :: istat
        character(:), allocatable :: format_
        type(FileObject) :: input_file
        type(FileInfo) :: info
!------
        body: block
            if (present(format)) then
                format_ = format
            else
                if (filename == "CONFIG" .or. filename == "REVCON" .or. filename == "HISTORY") then
                    format_ = STRUCTURE_DLPOLY
                else if (filename == "POSCAR" .or. filename == "CONTCAR") then
                    format_ = STRUCTURE_VASP
                end if
            end if

            ! Try to open the file
            call input_file%init(filename, mode="read", stat=istat)
            if (istat /= 0) then
                exit body
            end if

            if (allocated(format_)) then
                select case(lower_case(format_))
                case ("abinit")
                    call read_abinit(input_file, box, positions, tags, info=info, stat=istat)
                case ("castep-output")
                case ("cell", "castep")
                    call read_castep(input_file, box, positions, tags, velocities, info=info, stat=istat)
                case ("cif")
                    call read_cif(input_file, box, positions, tags, info=info, stat=istat)
                case ("history", "config", "revcon", "dlpoly", "dl_poly")
                    call read_dlpoly(input_file, box, positions, tags, velocities, forces, info=info, stat=istat)
                case ("xfg", "cfg")
                    call read_xfg(input_file, box, positions, tags, masses=masses, &
                        velocities=velocities, forces=forces, aux=aux, info=info, stat=istat)
                case ("gin", "gulp")
                    call read_gulp(input_file, box, positions, tags, info=info, stat=istat)
                case ("gout", "gulp-output")
                case ("lmp", "lconfig", "lammps")
                    call read_lammps(input_file, box, positions, tags, info=info, stat=istat)
                case ("xyz")
                    call read_xyz(input_file, box, positions, tags, masses=masses, &
                        velocities=velocities, forces=forces, aux=aux, info=info, stat=istat)
                case ("xsf")
                    call read_xsf(input_file, box, positions, tags, info=info, stat=istat)
                case ("poscar", "contcar", "vasp")
                    call read_vasp(input_file, box, positions, tags, velocities, info=info, stat=istat)
                case default

                end select
            else

                ! First attempt: extended CFG
                call read_xfg(input_file, box, positions, tags, masses=masses, &
                    velocities=velocities, forces=forces, aux=aux, info=info, stat=istat)

                ! Second attempt: XYZ
                if (istat /= 0) then
                    call istat%discard
                    call input_file%rewind
                    call read_xyz(input_file, box, positions, tags, masses=masses, &
                        velocities=velocities, forces=forces, aux=aux, info=info, stat=istat)
                end if

                ! Third attempt: LAMMPS
                if (istat /= 0) then
                    call istat%discard
                    call input_file%rewind
                    call read_lammps(input_file, box, positions, tags, info=info, stat=istat)
                end if

                ! Fourth attempt: VASP
                if (istat /= 0) then
                    call istat%discard
                    call input_file%rewind
                    call read_vasp(input_file, box, positions, tags, velocities, info=info, stat=istat)
                end if

                ! Fifth attempt: CIF
                if (istat /= 0) then
                    call istat%discard
                    call input_file%rewind
                    call read_cif(input_file, box, positions, tags, info=info, stat=istat)
                end if

                ! Sixth attempt: DL_POLY
                if (istat /= 0) then
                    call istat%discard
                    call input_file%rewind
                    call read_dlpoly(input_file, box, positions, tags, &
                        velocities=velocities, forces=forces, info=info, stat=istat)
                end if

                ! Seventh attempt: GULP
                if (istat /= 0) then
                    call istat%discard
                    call input_file%rewind
                    call read_gulp(input_file, box, positions, tags, info=info, stat=istat)
                end if

                ! Eighth attempt: CASTEP
                if (istat /= 0) then
                    call istat%discard
                    call input_file%rewind
                    call read_castep(input_file, box, positions, tags, &
                    velocities=velocities, info=info, stat=istat)
                end if

                ! Ninth attempt: XSF
                if (istat /= 0) then
                    call istat%discard
                    call input_file%rewind
                    call read_xsf(input_file, box, positions, tags, info=info, stat=istat)
                end if

                ! Tenth attempt: Abinit
                if (istat /= 0) then
                    call istat%discard
                    call input_file%rewind
                    call read_abinit(input_file, box, positions, tags, info=info, stat=istat)
                end if

                ! No format worked
                if (istat /= 0) then
                    call istat%discard
                    call istat%raise("Could not determine file format.")
                end if
            end if

!------ Make sure the velocities array is allocated if it is present
            if (present(velocities)) then
                if (.not.allocated(velocities)) then
                    allocate(velocities, mold=positions)
                    velocities = 0
                end if
            end if

!------ Make sure the forces array is allocated if it is present
            if (present(forces)) then
                if (.not.allocated(forces)) then
                    allocate(forces, mold=positions)
                    forces = 0
                end if
            end if
        end block body

        if (present(file_info)) file_info = info

        ! Report any error to the caller if provided with one of the optional arguments
        if (present(msg)) then
            if (istat==0) then
                msg = ""
            else
                msg = istat%string()
            end if
        end if
        if (present(stat)) then
            stat = 0
            if (istat /= 0) stat = -1
        end if
        if (present(except)) call except%transfer(istat)
        if (present(stat) .or. present(msg) .or. present(except)) call istat%discard
    end subroutine
!******************************************************************************!
!> Write a structure to a file.
!******************************************************************************!
    subroutine write_structure(filename, box, positions, tags, format, options, compression, velocities, forces, masses, aux, stat)
        character(*), intent(in) :: filename
        real(e), dimension(3,3), intent(in) :: box
        real(e), dimension(:,:), intent(in) :: positions
        character(TAG_LENGTH), dimension(:), intent(in) :: tags
        character(*), intent(in), optional :: format
        character(*), intent(in), optional :: compression
        real(e), dimension(:,:), intent(in), optional :: velocities
        real(e), dimension(:,:), intent(in), optional :: forces
        real(e), dimension(:), intent(in), optional :: masses
        type(DataTable), intent(in), optional :: aux
        integer, intent(in), optional :: options
        type(FException), intent(out), optional :: stat
!------
        character(:), allocatable :: format_, compression_, extension, basename, fext
        integer :: i
        type(FException) :: istat
!------
        body: block

!------ Set the compression format
            if (present(compression)) then
                ! Use the format provided by the caller if there is one
                compression_ = compression
            else
                ! Otherwise, try to guess from the file extension
                compression_ = COMPRESSION_NONE
                i = index(filename, ".", back=.true.)
                if (i > 0) then
                    select case(lower_case(filename(i+1:)))
                    case ("gz", "gzip")
                        compression_ = COMPRESSION_GZIP
                        basename = filename(:i-1)
                    case ("bz2", "bzip2")
                        compression_ = COMPRESSION_BZIP2
                        basename = filename(:i-1)
                    case ("zstd", "zst")
                        compression_ = COMPRESSION_ZSTD
                        basename = filename(:i-1)
                    end select
                end if
            end if
            if (.not.allocated(basename)) basename = filename

            ! Sanity check: is the compression format valid?
            if (compression_/=COMPRESSION_NONE .and. &
                compression_/=COMPRESSION_GZIP .and. &
                compression_/=COMPRESSION_BZIP2 .and. &
                compression_/=COMPRESSION_ZSTD) then

                call istat%raise("Unknown compression format.")
                exit body
            end if

!------ Set the file format
            ! Use the format provided by the caller if there is one
            if (present(format)) then
                if (format /= UNDEF .and. format /= "" .and. format /= "*") then
                    format_ = format
                end if
            end if

            ! Otherwise, try to guess from the file name
            if (.not.allocated(format_)) then
                if (filename == "CONFIG" .or. filename == "REVCON" .or. filename == "HISTORY") then
                    format_ = STRUCTURE_DLPOLY
                else if (filename == "POSCAR" .or. filename == "CONTCAR") then
                    format_ = STRUCTURE_VASP
                end if
            end if

            ! Otherwise, try to guess from the file extension
            if (.not.allocated(format_)) then
                i = index(basename, ".", back=.true.)
                if (i > 0) then
                    extension = basename(i+1:)
                end if

                if (allocated(extension)) then
                    i = index(extension, ".")
                    if (i>0) then
                        fext = lower_case(extension(:i-1))
                    else
                        fext = lower_case(extension)
                    end if

                    select case(lower_case(fext))
                    case ("abinit")
                        format_ = STRUCTURE_ABINIT
                    case ("cell", "castep")
                        format_ = STRUCTURE_CASTEP
                    case ("castep-out")
                        format_ = STRUCTURE_CASTEP_OUT
                    case ("cif")
                        format_ = STRUCTURE_CIF
                    case ("config", "revcon", "dlpoly", "dl_poly", "history")
                        format_ = STRUCTURE_DLPOLY
                    case ("xfg", "cfg")
                        format_ = STRUCTURE_XFG
                    case ("gin", "gulp", "res")
                        format_ = STRUCTURE_GULP
                    case ("gout", "gulp-output")
                        format_ = STRUCTURE_GULP_OUT
                    case ("lmp", "lconfig", "lammps")
                        format_ = STRUCTURE_LAMMPS
                    case ("xsf")
                        format_ = STRUCTURE_XSF
                    case ("xyz")
                        format_ = STRUCTURE_XYZ
                    case ("poscar", "contcar", "vasp")
                        format_ = STRUCTURE_VASP
                    end select
                else
                    select case (lower_case(filename))
                        case ("config", "revcon", "HISTORY")
                            format_ = STRUCTURE_DLPOLY
                        case ("poscar", "contcar")
                            format_ = STRUCTURE_VASP
                    end select
                end if
            end if

            ! Use extended cfg as the default file format if all else fails
            if (.not.allocated(format_)) format_ = STRUCTURE_XFG

            ! Call the right writer
            select case(format_)
            case (STRUCTURE_ABINIT)
            case (STRUCTURE_CASTEP)
                call write_castep(filename, box, positions, tags, velocities=velocities, &
                    stat=istat)
            case (STRUCTURE_CIF)
                call write_cif(filename, box, positions, tags, stat=istat)
            case (STRUCTURE_CFG)
            case (STRUCTURE_DLPOLY)
                call write_dlpoly(filename, box, positions, tags, velocities=velocities, &
                    forces=forces, stat=istat)
            case (STRUCTURE_GULP)
                call write_gulp(filename, box, positions, tags, stat=istat)
            case (STRUCTURE_LAMMPS)
                call write_lammps(filename, box, positions, tags, stat=istat)
            case (STRUCTURE_VASP, "poscar")
                call write_vasp(filename, box, positions, tags, stat=istat)
            case (STRUCTURE_XFG)
                call write_xfg(filename, box, positions, tags, velocities=velocities, &
                    masses=masses, forces=forces, aux=aux, stat=istat)
            case (STRUCTURE_XSF)
                call write_xsf(filename, box, positions, tags, &
                    masses=masses, stat=istat)
            case (STRUCTURE_XYZ)
                call write_xyz(filename, box, positions, tags, velocities=velocities, &
                    masses=masses, forces=forces, aux=aux, stat=istat)
            case default
                call istat%raise("Unknown file format: " // format_)
            end select
        end block body

        if (present(stat)) call stat%transfer(istat)
    end subroutine
!******************************************************************************!
!> Count the number of frames defined in a given file.
!>
!> @todo
!> Not implemented yet.
!> @endtodo
!******************************************************************************!
    subroutine count_frames_in_file(filename, frames_count, format, compression, stat, msg)
        character(*), intent(in) :: filename
        integer, intent(out) :: frames_count
        character(*), intent(in), optional :: format
        character(*), intent(in), optional :: compression
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: msg
!------

    end subroutine
!******************************************************************************!
!> Try to determine the format of a given file. This is not very efficient as
!> it tries to read the file using the different formats successively.
!******************************************************************************!
    subroutine get_file_format(filename, format, stat)
        character(*), intent(in) :: filename !< Name of the file to test
        character(:), allocatable, intent(out) :: format !< Detected file format
        type(FException), intent(out), optional :: stat !< Status
!------
        type(FException) :: istat
        real(e), dimension(3,3) :: box
        real(e), dimension(:,:), allocatable :: positions
        character(TAG_LENGTH), dimension(:), allocatable :: tags
        integer :: err
        type(FileInfo) :: info
!------

        body: block
            call read_structure(filename, box, positions, tags, file_info=info, stat=err)
            if (err /= 0) call istat%raise("Unknown file format.")
            format = info%data_format
        end block body

        if (present(stat)) call stat%transfer(istat)
    end subroutine
end module
