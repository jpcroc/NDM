!******************************************************************************!
!                    dk_structure_io_lammps submodule
!------------------------------------------------------------------------------!
!> Specific procedures to read and write LAMMPS structure files.
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
submodule (dk_structure_io) dk_structure_io_lammps
    use dk_math, only: cross_product, det, matrix_inverse

    use dk_datatable, only: DataTable
    use dk_string, only: FString

    use dk_exception, only: FException, FExceptionDescription, EVENT_LEVEL_LOG
    use ext_character

    implicit none(external, type)

contains
!******************************************************************************!
!> Read a structure from a LAMMPS file.
!******************************************************************************!
    module subroutine read_lammps(input_file, box, positions, tags, velocities, info, stat)
        type(FileObject), intent(inout) :: input_file
        real(e), dimension(3,3), intent(out) :: box
        real(e), dimension(:,:), allocatable, intent(out) :: positions
        character(TAG_LENGTH), dimension(:), allocatable, intent(out) :: tags
        real(e), dimension(:,:), allocatable, intent(out), optional :: velocities
        type(FileInfo), intent(out), optional :: info
        type(FException), intent(out), optional :: stat
!------
        integer :: i, n, iostat, j, atoms_count
        logical :: atoms
        character(:), allocatable :: format
        character(:), allocatable :: buffer, keyword, line
        type(FException) :: stat_
        real(e), dimension(3) :: range
        integer, dimension(3) :: intDum
        real(e), dimension(7) :: realDum
        type(FException) :: istat
        real(e), dimension(3,3) :: ibox
!------

        body: block

            n = 0
            atoms = .false.

!------ Read the header

            ! Ignore the first line
            call input_file%read_line(line, stat=istat)
            if (istat /= 0) exit body

            box = 0

            do while(stat_ == 0)
                call input_file%read_line(line, stat=stat_)
                if (adjustl(line) == "Atoms" .or. adjustl(line) == "Velocities" .or. adjustl(line) == "Masses") exit

                ! Work on a copy of the line to keep an unmodified version for error messages
                buffer = adjustl(line)

                ! Remove comments
                if (index(buffer, "#") > 0) buffer(index(buffer,"#"):) = ""

                ! Ignore empty lines
                if (buffer == "") cycle

                ! Get the keyword and values
                if (allocated(keyword)) deallocate(keyword)
                allocate(character(100) :: keyword)
                read(buffer,*, iostat=iostat) range(1), range(2), range(3), keyword
                if (iostat /= 0) read(buffer,*, iostat=iostat) range(1), range(2), keyword
                if (iostat /= 0) read(buffer,*, iostat=iostat) range(1), keyword
                if (iostat /= 0) then
                    call istat%raise("Error in LAMMPS file: " // input_file%name // &
                        ": syntax error in file header: " // trim(buffer))
                    exit body
                end if

                keyword = trim(buffer(index(buffer,trim(keyword)):))

                select case(keyword)
                case("xlo xhi")
                    read(buffer(:index(buffer, "xlo xhi")-1),*,iostat=iostat) range(1), range(2)
                    if (iostat == 0) then
                        box(1,1) = range(2) - range(1)
                    else
                        call istat%raise("Error in LAMMPS file: " // input_file%name // &
                            ": invalid box definition along the x direction: " // trim(buffer))
                        exit body
                    end if
                case("ylo yhi")
                    read(buffer(:index(buffer, "ylo yhi")-1),*,iostat=iostat) range(1), range(2)
                    if (iostat == 0) then
                        box(2,2) = range(2) - range(1)
                    else
                        call istat%raise("Error in LAMMPS file: " // input_file%name // &
                            ": invalid box definition along the y direction: " // trim(buffer))
                        exit body
                    end if
                case("zlo zhi")
                    read(buffer(:index(buffer, "zlo zhi")-1),*,iostat=iostat) range(1), range(2)
                    if (iostat == 0) then
                        box(3,3) = range(2) - range(1)
                    else
                        call istat%raise("Error in LAMMPS file: " // input_file%name // &
                            ": invalid box definition along the z direction: " // trim(buffer))
                        exit body
                    end if
                case("xy xz yz")
                    read(buffer(:index(buffer, "xy xz yz")),*,iostat=iostat) range(1), range(2), range(3)
                    box(1,2) = range(1)
                    box(1,3) = range(2)
                    box(2,3) = range(3)
                case("atoms")
                    read(buffer(:index(buffer, "atoms")-1),*,iostat=iostat) atoms_count
                    if (iostat /= 0) then
                        call istat%raise("the data do not follow LAMMPS format: " // &
                            "wrong number of atoms: " // trim(buffer) // ".")
                        exit body
                    end if
                    atoms = .true.
                case("bonds")
                case("angles")
                case("dihedrals")
                case("impropers")
                case("atom types")
                case("bond types")
                case("angle types")
                case("dihedral types")
                case("improper types")
                case("extra bond per atom")
                case("extra angle per atom")
                case("extra dihedral per atom")
                case("extra improper per atom")
                case("extra special per atom")
                case("ellipsoids")
                case("lines")
                case("triangles")
                case("bodies")
                case("")
                case default
                    call istat%raise("the data do not follow LAMMPS format: " // &
                        "unknown keyword """ // trim(keyword) // """ in header.")
                    exit body
                end select
            end do

            ! Check that the number of atoms is defined
            if (.not.atoms) then
                call istat%raise("the data do not follow LAMMPS format: " // &
                    "missing atoms keyword.")
                exit body
            end if
            allocate(tags(atoms_count))
            allocate(positions(3,atoms_count))
            if (present(velocities)) allocate(velocities(3,atoms_count))

            ! Check that the box matrix is not singular
            if (det(box) < 1.0e-5_e) then
                call istat%raise("Error in LAMMPS file: " // input_file%name // &
                    ": box matrix is singular.")
            end if
            ibox = matrix_inverse(box)

!------ Read the atomic positions
            ! Look for the atoms keyword
            do while(line /= "Atoms")
                call input_file%read_line(line, stat=stat_)
                if (stat_ /= 0) then
                    call istat%raise("Error in LAMMPS file: " // input_file%name // &
                    ": unexpected end of file before atoms section.")
                    exit body
                end if
            end do

            call input_file%read_line(line, stat=stat_)
            if (line /= "") then
                call istat%raise("Syntax error in LAMMPS file: " // input_file%name // &
                    ": expected an empty line after the Atoms keyword.")
                exit body
            end if

            ! Guess data format from first atom
            format = UNDEF
            detectFormat: block
                call input_file%read_line(buffer, stat=stat_)
                do while(buffer == "" .and. stat_ == 0)
                    call input_file%read_line(buffer, stat=stat_)
                end do
                read(buffer,*,iostat=iostat) intDum(1:3), realDum(1:4)
                if (iostat == 0) then
                    format = "full"
                    exit detectFormat
                end if
                read(buffer,*,iostat=iostat) intDum(1:3), realDum(1:3)
                if (iostat == 0) then
                    format = "angle"
                    exit detectFormat
                end if
                read(buffer,*,iostat=iostat) intDum(1:2), realDum(1:7)
                if (iostat == 0) then
                    format = "dipole"
                    exit detectFormat
                end if
                read(buffer,*,iostat=iostat) intDum(1:2), realDum(1:6)
                if (iostat == 0) then
                    format = "electron"
                    exit detectFormat
                end if
                read(buffer,*,iostat=iostat) intDum(1:2), realDum(1:4)
                if (iostat == 0) then
                    format = "charge"
                    exit detectFormat
                end if
                iostat = 0
                read(buffer,*,iostat=iostat) intDum(1:2), realDum(1:3)
                if (iostat == 0) then
                    format = "atomic"
                    exit detectFormat
                end if

                call istat%raise("Error in LAMMPS file: " // input_file%name // &
                ": could not detect atoms data format.")
                exit body
            end block detectFormat
            write(*,*) "LAMMPS file format: " // format

            ! Read atoms data
            do i=1, atoms_count
                select case(format)
                    case("angle")
                        read(buffer,*,iostat=iostat) intDum(1:3), realDum(1:3)
                        j = intDum(1)
                        positions(:,j) = matmul(ibox, realDum(1:3))
                        tags(j) = "" // intDum(3)
                    case("atomic")
                        read(buffer,*,iostat=iostat) intDum(1:2), realDum(1:3)
                        j = intDum(1)
                        positions(:,j) = matmul(ibox, realDum(1:3))
                        tags(j) = "" // intDum(2)
                    case("charge")
                        read(buffer,*,iostat=iostat) intDum(1:2), realDum(1:4)
                        j = intDum(1)
                        positions(:,j) = matmul(ibox, realDum(2:4))
                        tags(j) = "" // intDum(2)
                    case("dipole")
                        read(buffer,*,iostat=iostat) intDum(1:2), realDum(1:7)
                        j = intDum(1)
                        positions(:,j) = matmul(ibox, realDum(2:4))
                        tags(j) = "" // intDum(2)
                    case("electron")
                        read(buffer,*,iostat=iostat) intDum(1:2), realDum(1:6)
                        j = intDum(1)
                        positions(:,j) = matmul(ibox, realDum(4:6))
                        tags(j) = "" // intDum(2)
                    case("full")
                        read(buffer,*,iostat=iostat) intDum(1:3), realDum(1:4)
                        j = intDum(1)
                        positions(:,j) = matmul(ibox, realDum(2:4))
                        tags(j) = "" // intDum(3)
                    case default

                end select
                if (i < atoms_count) call input_file%read_line(buffer, stat=stat_)
            end do

!------ Read the atomic velocities
            if (present(velocities)) then

                ! Look for velocities keyword
                do while(lower_case(buffer) /= "velocities")
                    call input_file%read_line(buffer, stat=istat)
                    if (istat /= 0) then
                        exit body
                    end if
                end do
            end if

            if (present(info)) then
                info%data_format = trim(STRUCTURE_LAMMPS)
                info%num_frames = 1
                info%compression_format = input_file%compression_format()
            end if
        end block body

        ! Cleanup
        call stat_%discard
        if (present(stat)) call stat%transfer(istat)
    contains
        pure function is_alphabetic(char)
            character(1), intent(in) :: char
            logical :: is_alphabetic
!------
            is_alphabetic = (iachar(char) >= iachar('a') .and. iachar(char) <= iachar('z')) .or. &
                (iachar(char) >= iachar('A') .and. iachar(char) <= iachar('Z'))
        end function
    end subroutine
!******************************************************************************!
!> Write a structure to a LAMMPS data file.
!******************************************************************************!
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
!------
        integer :: a, i, atoms_count, nelements
        logical, dimension(:), allocatable :: mask
        character(1000) :: buffer
        logical :: has_velocities
        real(e) :: ax, bx, by, cx, cy, cz
        real(e) :: lx, ly, lz, xy, xz, yz
        real(e), dimension(6) :: parameters
        real(e), dimension(3) :: x0, x
        real(e), dimension(3,3) :: T, S
        integer, dimension(:), allocatable :: moleculeId, element_id
        type(Element), dimension(:), allocatable :: elements_
        type(FileObject) :: output_file
        type(FException) :: istat
!------

        body: block

!------ Prepare the file
            call output_file%init(file, "w", istat)

            ! Return if the file could not be opened
            if (istat /= 0) exit body

            ! Set whether we need to save atomic velocities
            has_velocities = present(velocities)

!------ Build the inclusion mask
            if (present(include)) then
                if (size(include) /= size(positions, 2)) then
                    call istat%raise("Wrong mask size")
                    return
                end if
                allocate(mask, source=include)
            else
                allocate(mask(size(positions, 2)))
                mask = .true.
            end if
            atoms_count = count(mask)

            ! Set the element id of each atom
            allocate(element_id(size(tags)))
            element_id = -42
            if (present(elements)) then
                nelements = size(elements)
                do i=1, size(tags)
                    do a=1, nelements
                        if (tags(i) == elements(a)%tag) then
                            element_id(i) = a
                            exit
                        end if
                    end do
                end do
                allocate(elements_, source=elements)
            else
                allocate(elements_(atoms_count))
                nelements = 0
                atom_loop: do i=1, size(tags)
                    do a=1, nelements
                        if (tags(i) == elements_(a)%tag) then
                            element_id(i) = a
                            cycle atom_loop
                        end if
                    end do
                    nelements = nelements + 1
                    elements_(nelements)%tag = tags(i)
                    if (present(masses)) elements_(nelements)%mass = masses(i)
                    if (present(charges)) elements_(nelements)%charge = charges(i)
                    element_id(i) = nelements
                end do atom_loop
            end if

!------ Write the file header
            call output_file%write_line(get_mark())
            call output_file%write_line("")
            call output_file%write_line(atoms_count // " atoms")
            call output_file%write_line(nelements // " atom types")
            call output_file%write_line("0 bonds")
            call output_file%write_line("0 bond types")
            call output_file%write_line("")

            parameters = parameters_with_box(box)
            lx = parameters(1)
            xy = parameters(2)*cos(parameters(6)*PI/180)
            xz = parameters(3)*cos(parameters(5)*PI/180)
            ly = sqrt(parameters(2)**2 - xy**2)
            yz = (parameters(2)*parameters(3)*cos(parameters(4)*PI/180) - xy * xz) / ly
            lz = sqrt(parameters(3)**2 - xz**2 - yz**2)
            write(buffer,'(f0.6,1x,f0.6,1x,a)') 0.0_e, lx, "xlo xhi"
            call output_file%write_line(trim(buffer))
            write(buffer,'(f0.6,1x,f0.6,1x,a)') 0.0_e, ly, "ylo yhi"
            call output_file%write_line(trim(buffer))
            write(buffer,'(f0.6,1x,f0.6,1x,a)') 0.0_e, lz, "zlo zhi"
            call output_file%write_line(trim(buffer))
            write(buffer,'(f0.6,1x,f0.6,1x,f0.6,1x,a)') xy, xz, yz, "xy xz yz"
            call output_file%write_line(trim(buffer))
            call output_file%write_line("")

            S = 0
            ax = parameters(1)
            bx = parameters(2) * cos(parameters(6)*PI/180)
            by = parameters(2) * sin(parameters(6)*PI/180)
            cx = parameters(3) * cos(parameters(5)*PI/180)
            cy = (dot_product(box(:,2), box(:,3)) - bx*cx) / by
            cz = sqrt(parameters(3)**2 - cx**2 - cy**2)

!------ Write the atoms section
            call output_file%write_line("")
            call output_file%write_line("Atoms")
            call output_file%write_line("")
            T(1,:) = cross_product(box(:,2), box(:,3))
            T(2,:) = cross_product(box(:,3), box(:,1))
            T(3,:) = cross_product(box(:,1), box(:,2))
            T  = T / det(box)
            S = 0
            S(1,1) = ax
            S(1,2) = bx
            S(2,2) = by
            S(1,3) = cx
            S(2,3) = cy
            S(3,3) = cz
            T = matmul(S, T)
            do i=1, atoms_count
                if (.not.mask(i)) cycle

                x0 = matmul(box, positions(:,i))
                x = matmul(T, x0)

                ! Using the atom_style full:
                ! atom-ID molecule-ID atom-type q x y z
                if (allocated(moleculeId)) then
                    write(buffer,*) i, moleculeId(i), element_id(i), elements_(element_id(i))%charge, x
                    call output_file%write_line(trim(buffer))
                else
                    write(buffer,*) i, i, element_id(i), elements_(element_id(i))%charge, x
                    call output_file%write_line(trim(buffer))
                end if
            end do

!------ Write the velocities section
            if (has_velocities) then
                call output_file%write_line("")
                call output_file%write_line("Velocities")
                call output_file%write_line("")
                do i=1, atoms_count
                    write(buffer,*) i, velocities(:,i)
                    call output_file%write_line(trim(buffer))
                end do
            end if

!------ Write the bonds section
!         if (associated(shell)) then
!             write(unit,*)
!             write(unit,*) "Bonds"
!             write(unit,*)
!             n = 0
!             do i=1, natoms
!                 if (shell(i) == 0) cycle
!
!                 n = n + 1
!                 key = model%coreKey(config%tag(i))
!                 write(unit,*) n, key, i, nint(shell(i))
!             end do
!         end if
        end block body

        call output_file%close

        if (present(stat)) call stat%transfer(istat)
    end subroutine

end submodule
