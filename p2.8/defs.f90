module defs
! This module defines all variables used accross the program ART01
!
!  Copyright N. Mousseau, May 2001

  implicit none
  save

  real(8)::ENERGYDIFFERENCE

  real(8), parameter :: version = 1.003   ! Numero de version
  real(8) :: TEMPERATURE        ! Temperature in eV
  integer :: NATOMS             ! Number of atoms in the system
!  integer :: MAXNEI             ! Maximum number of nearest neighbours
  integer :: VECSIZE     ! Length of the force and position vectors

  integer :: NUMBER_EVENTS      ! Total number of events in this run
  logical :: NEW_EVENT          ! Total number of events in this run

  integer, parameter :: FCONF        = 1         ! Units for printing/reading
  integer, parameter :: FCOUNTER     = 2        
  integer, parameter :: FLIST        = 3         
  integer, parameter :: FLOG         = 4         
  integer, parameter :: FREFCONFIG   = 11
  integer, parameter :: FSTARTCONF   = 12
  integer, parameter :: FRESTART     = 9        
  integer, parameter :: XYZ          = 13  
  ! Name of the file storing the current configurations
  character(len=15) :: conf_initial, conf_saddle, conf_final
  
  integer, dimension(:), allocatable         :: type                ! Atomic type
  real(8), dimension(:), allocatable, target  :: force      ! Working forces on the atoms
  real(8), dimension(:), allocatable, target  :: pos        ! Working positions of the atoms
  real(8), dimension(:), allocatable, target  :: posref     ! Reference position
  real(8), dimension(:), allocatable, target  :: direction_restart  
  real(8), dimension(:), allocatable :: initial_direction  ! Initial move for leaving harmonic  well
  character(len=5), dimension(:), allocatable :: Atom

  real(8), dimension(:), pointer :: x, y, z   ! Pointers for working position
  real(8), dimension(:), pointer :: xref, yref, zref   ! Pointers for reference position
  real(8), dimension(:), pointer :: fx, fy, fz   ! Pointers for working force


  real(8) :: boxref(3)                         ! Reference boxsize
  real(8) :: box(3)                            ! Working boxsize

  real(8) :: scalaref                       ! Reference volume scaling
  real(8) :: scala                          ! Working volume scaling
  real(8) :: fscala                         ! Working forces on volume

  logical :: restart                        ! State of restart (true or false)
  integer :: state_restart                  ! start of restart (1 - harmonic well, 2 - 
                                            ! activation, 3 - relaxation)

  integer :: preferred_atom                 ! Atom at the center of the event
  real(8) :: radius_initial_deformation     ! Radius of the initial deformation

  integer :: ievent                         ! actual number of events
  integer :: iter_restart                   ! iteraction number of restart
  integer :: ievent_restart                 ! Event number at restart
  integer :: mincounter                     ! Counter for output files
  integer :: refcounter                     ! Id of reference file

  logical :: print_details                  ! Print or not the details of activation
  integer :: kprint, iprint, mprint         ! Printing interval for details of activation

  integer :: evalf_number                   ! Number of force evalutions

  real(8) :: total_energy, ref_energy       ! Energies

  character(len=20) :: LOGFILE
  character(len=20) :: EVENTSLIST
  character(len=20) :: REFCONFIG
  character(len=3)  :: FINAL
  character(len=3)  :: SADDLE 
  character(len=11) :: COUNTER
  character(len=11) :: RESTARTFILE
  character(len=3)  :: bulk_structure
  real (8) :: art_ann,art_Rinst
  integer  :: art_ninst
end module defs
