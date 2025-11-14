

module dk_structurereader
    use ext_character, only: split

    use dk_structure, only: Structure
    use dk_exception, only: FException
    use dk_structure_io, only: read_structure, FileInfo

    implicit none(external, type)

    type, public :: StructureReader
        integer :: files_count = 0
        integer :: num_frames = 0
        character(:), allocatable :: source
        character(:), dimension(:), allocatable :: files
        integer, dimension(:), allocatable :: frames_count
        integer, dimension(:), allocatable :: read_frames
        character(12), dimension(:), allocatable :: format
    contains
        procedure, public :: init
        procedure, public :: get_next_frame
    end type

contains
    subroutine init(this, source)
        class(StructureReader), intent(out) :: this
        character(*), intent(in) :: source
!------
        integer :: f
!------
        call split(source, " ", this%files, .true.)

        ! Check whether we have a glob
        if (size(this%files) == 1) then
            ! Do nothing for now and assume the source is a single file
        end if

        ! Set up the reader object
        this%files_count = size(this%files)
        allocate(this%frames_count(this%files_count))
        this%frames_count = 0
        allocate(this%format(this%files_count))
        allocate(this%read_frames(this%files_count))
        this%read_frames = 0

        ! Process the list of files
        do f=1, this%files_count

        end do
    end subroutine
!******************************************************************************!
!> Get the next frame in a sequence. The frame will be the next one in the
!> current file, or the first one from the next file in the sequence.
!******************************************************************************!
    subroutine get_next_frame(this, struct, stat)
        class(StructureReader), intent(inout) :: this
        type(Structure), intent(out) :: struct
        type(FException), intent(out), optional :: stat
!------
        integer :: i
        type(FException) :: istat
        type(FileInfo) :: info
!------
        body: block
            do i=1, this%files_count
                if (this%read_frames(i) < this%frames_count(i)) then
                    call read_structure(this%files(i), struct%box, struct%r, struct%tag, &
                        velocities=struct%v, forces=struct%f, aux=struct%aux, file_info=info, except=istat)
                    this%read_frames(i) = this%read_frames(i)+ 1
                    this%format(i) = info%data_format
                    exit body
                end if
            end do

            ! If we reach this point, it means that we read all the frames in the sequence
            call istat%raise("Reached end of sequence")
        end block body

        if (present(stat)) call stat%transfer(istat)
    end subroutine
end module
