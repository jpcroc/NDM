!> @file
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> ART Module defs
!! This module defines almost all variables used accross the program ART


! $Revision$
! $Date: 2015-12-14 11:01:46 -0500 (Lun, 14 déc 2015) $
! $Id: defs.f90 1534 2015-12-14 15:59:38Z mickael $


module defs

#ifdef PARA
  use mpi
#endif

  implicit none

  integer, parameter :: dp = kind(0.d0) ! double precision
  real(kind=8), parameter :: VERSION_NUMBER  = 3.0    ! Version of the code
  character(len=20), parameter :: BIGREVNO ="1.6-dev.12-538" ! Bigdft version

  integer::ivisuart
!  logical lwgin rather from gen_com_m
  integer::unit6P

  
  real(kind=8) :: t1                 ! Initial Date (cputime)

  real(kind=8) :: my_gnrm            ! We use a higher convergence criterion for the
                                     ! wavefunction optimization in the Lanczos procedure

  integer      :: iproc, nproc       ! MPI proc identificators
  integer      :: INFLECTION
  logical      :: setup_initial      ! The tests done for knowing the best set of parameters
                                     ! for doing lanczos. Four times Lanczos for the minimum, and
                                     ! four for the inflection point.

  real(kind=8) :: TEMPERATURE        ! Temperature in eV
  integer      :: NATOMS             ! Number of atoms in the system
  integer      :: MAXNEI             ! Maximum number of nearest neighbours
  integer      :: VECSIZE            ! Length of the force and position vectors
  integer      :: nbr_quantum        !Number of quantum atoms for BSW force_calc. These are the first atoms in the input file
  integer      :: nbr_quantum_trash  !Number of quantum atoms where the forces will be wrong (to close to the H bonds). The are the very first atoms in the file
                                     !Thus nbr_quantum = nbr_quantum_trash + nbr_quantum_good
  integer      :: nbr_to_fit         !number of atoms where we fit SW to the hybrid computation
  logical      :: passivate          !Do we passivate the cluster with H ?
  logical, dimension(:), allocatable :: should_fit  !should we fit this atom or not ?

!  character(len=3)    :: energy_type  !To choose the energy calc routine ! replaced by NDM choices
  integer             :: DIM          !Spatial dimension

  integer      :: NUMBER_EVENTS      ! Total number of events in this run
  logical      :: NEW_EVENT          ! Total number of events in this run


  ! LOCAL FORCE
  logical :: LOCAL_FORCE             ! If true, uses local forces to compute generic events
  logical :: LOCAL_LANCZOS           ! If true, makes only a local deformation for the initial random vector in Lanczos
  logical :: GLOBAL_CONVERGENCE      ! If true, use global forces for final convergence to saddle of specific events
  logical :: use_local_forces        ! If true, use local forces in the force calculation
  real(kind=8)  :: INNER_REGION      ! maximum radius for inner region
  real(kind=8)  :: OUTER_REGION      ! maximum radius of outer region
  real(kind=8)  :: local_ref_energy  ! reference local energy to ensure total energy is computed properly
  real(kind=8)  :: global_ref_energy ! rererence global energy for the same goal

  ! Units for printing/reading

  integer, parameter :: FCONF       = 1
  integer, parameter :: FCOUNTER    = 2
  integer, parameter :: FLIST       = 3
  integer, parameter :: FLOG        = 4
  integer, parameter :: FREFCONFIG  = 11
  integer, parameter :: FSTARTCONF  = 12
  integer, parameter :: FRESTART    = 9
  integer, parameter :: XYZ         = 13
  integer, parameter :: CRESTART    = 14
  integer, parameter :: ASCII       = 15

  ! Name of the file storing the current configurations
  character(len=20) :: conf_initial, conf_saddle, conf_final

  ! Parameter used to count nb of atoms that moved
  real(dp) :: NPART_DR_THRESHOLD ! In atomic distance units
                                 ! for instance Angstroms

  ! Name of the files storing the path of steepest descent
  character(len=30) :: conf_sdp

  integer,      dimension(:), allocatable          :: typat    ! Atomic type
  integer,      dimension(:), allocatable          :: constr   ! Constraint over atoms
  real(kind=8), dimension(:), allocatable, target  :: force    ! Working forces on the atoms
  real(kind=8), dimension(:), allocatable, target  :: pos      ! Working positions of the atoms
  real(kind=8), dimension(:), allocatable, target  :: posref   ! Reference position

  logical      :: write_xyz
  ! restart
  logical      :: write_restart_file
  logical      :: restart                   ! State of restart (true or false)
  real(kind=8), dimension(:),  allocatable, target :: direction_restart
  real(kind=8), dimension(:,:),allocatable, target :: diis_forces_restart
  real(kind=8), dimension(:,:),allocatable, target :: diis_pos_restart
  real(kind=8), dimension(:),  allocatable         :: diis_norm_restart

  integer      :: state_restart             ! start of restart (1 - harmonic well, 2 -
                                            ! lanczos, 3 - relaxation, 4-DIIS)
  integer      :: iter_restart              ! Lanczos iteraction number of restart
  integer      :: ievent_restart            ! Event number at restart
  integer      :: nsteps_after_eigen_min_r
  integer      :: maxter_r
  integer      :: atp                       ! No. of attempts for a given event
  real(kind=8) :: eigen_min_r
  real(kind=8) :: eigenvalue_r
  !__________________

  real(kind=8),     dimension(:), allocatable :: initial_direction  ! Initial move for leaving harmonic  well
  character(len=5), dimension(:), allocatable :: Atom
  character(len=5), dimension(5)              :: type_name  ! replaced by NDM
!CROCTODO replace type,ame by typat
  real(kind=8), dimension(:), pointer :: x, y, z           ! Pointers for working position
  real(kind=8), dimension(:), pointer :: xref, yref, zref  ! Pointers for reference position
  real(kind=8), dimension(:), pointer :: fx, fy, fz        ! Pointers for working force

  real(kind=8) :: PUSH_OVER                                ! Fraction of displacement for pushing of saddle point.
  real(kind=8) :: PUSH_OVER_RELATIVE_ERROR ! Adjust push over saddle, until
                                           ! relative error between energy
                                           ! variation and its 2nd order
                                           ! approximation becomes lower
                                           ! than this threshold
                                           ! - Optional

  character(len=1)          :: boundary
  real(kind=8),dimension(3) :: boxref         ! Reference boxsize
  real(kind=8),dimension(3) :: box,small_box    ! Working boxsize

  real(kind=8), dimension(3,3) :: cell			! Triclinic cell parameters defined following LAMMPS
  real(kind=8), dimension(3,3) :: invcell		! Matrix going from triclinic to orthorombic coordinates

  real(kind=8) :: scalaref                    ! Reference volume scaling
  real(kind=8) :: scala                       ! Working volume scaling
  real(kind=8) :: fscala                      ! Working forces on volume

  real(kind=8) :: sym_break_dist              ! Distance atoms are pushed to
                                              ! break symmetry

  integer      :: preferred_atom              ! Atom at the center of the event
  real(kind=8) :: radius_initial_deformation  ! Radius of the initial deformation

  integer      :: pas                         ! Counter: number of steps in the
                                              ! event

  integer      :: ievent                      ! actual number of events
  integer      :: mincounter                  ! Counter for output files
  integer      :: refcounter                  ! Id of reference file
  integer      :: evalf_number                ! Number of force evalutions
  logical      :: WRITE_REJECTED_EVENT        ! If true, also write failed events

  real(kind=8) :: total_energy, ref_energy  ! Energies

  real(kind=8) :: NORM_CRITERIUM_FIRE       ! Norm 2 criterium for stopping minimization
  real(kind=8) :: FMAX_CRITERIUM_FIRE       ! Norm inf criterium for stopping minimization
  real(kind=8) :: DT_MAX_FIRE               ! Maximum step for FIRE
  integer      :: MAX_PERP_ITER_FIRE        ! Maximum number of relaxation steps in the perp direction during activation
  integer      :: MAX_ITER_FIRE             ! Maximum number of relaxation step during global minimization
  integer      :: MAX_PERP_INCR_FIRE        ! Increment in the number of relaxation step during activation
  integer      :: MAX_LANCZOS_STEPS         ! Maximum number of Lanczos steps during activation
  logical      :: ITERATIVE                 ! Iterative use of Lanczos & DIIS
  logical      :: USE_DIIS                  ! Use DIIS for final convergence to saddle
  logical      :: DIIS_CHECK_EIGENVEC       ! Check whether the metastable point is a saddle
  integer      :: DIIS_MEMORY               ! Number of steps kept in memory
  integer      :: MAXPAS
  real(kind=8) :: DIIS_FORCE_THRESHOLD      ! Force Threshold to call the algorithm
  real(kind=8) :: DIIS_STEP                 ! Step used to update positions in DIIS

  real(kind=8) :: factor_diis               ! new position pos(i) is accepted if the norm
                                            ! of ( pos - previous_pos(i-1)) is
                                            ! lower than factor_diis*INCREMENT
  integer      :: maxdiis                   ! max allowed diis steps per call
  logical      :: SAVE_CONF_INT             ! Save the configuration at every step?
  !__________________

  ! Development :
  ! if delta_e < delta_thr .and. delr < delr_thr => end_activation = .true.

  logical     :: CHECK_CONNECTIVITY         ! Checks if reversing displacement from saddle leads to initial minimum

  logical      :: calc_proj                 ! if it is true; it calculates the projection only
                                            ! at every two steps but after 4 steps above of an
                                            ! inflection in the eigenvalue, and if the last
                                            ! a1  >0.9d0

  real(kind=8) :: delta_thr                 ! default = 0.0
  real(kind=8) :: delr_thr                  ! default = 0.0
  logical      :: clean_wf                  ! we are not going to use previously calculated wave
                                            ! function at some key points. This is for charged systems.

  ! for dual_search:
  ! for new projection vector. It will be start with a random displacement localized
  ! only around a given atom.
  ! Forces in the other molecule will be zero
  ! This is only for Type_of_events == local or list_local
  logical      :: dual_search
  integer      :: central_atom
  real(kind=8) :: size_system
  integer      :: type_sel
  !__________________

  character(len=20) :: LOGFILE
  character(len=20) :: EVENTSLIST
  character(len=64) :: REFCONFIG
  character(len=20) :: INPUT_REFINE_SADDLE
  character(len=11) :: RESTARTFILE
  character(len=3)  :: FINAL                ! is the prefix for minima files
  character(len=3)  :: SADDLE
  character(len=3)  :: CHECK
  character(len=11) :: COUNTER



  character(len=20) :: eventtype
  character(len=11),parameter :: pos_units = "angstroemd0"

  !-------------------------------------
 character(len=40) :: INPUT_LAMMPS_FILE !input file need by lammps
 logical :: firsttime_lammps  !trigger used to initialize lammps
 real(8) :: energy_conversion, position_conversion ! Factor for converting energy and position between force input and that used in the code
! character(len = 7):: UNITS_CONVERSION ! Flag for energy conversion !replaced by NDM
 integer :: orig_group,orig_comm ! MPI param
 integer,dimension(:),allocatable :: all_rank
 integer :: lammps_ranksize
 integer :: ranks(1)
 integer,parameter :: STOP_WORKER_TAG = 2
 integer,parameter :: LAMMPS_WORKER_TAG = 10
 integer,parameter :: MASTER_TAG = 11

 real(kind=8),allocatable:: atdisp(:)
 logical,allocatable ::ldisp(:)
  !-------------------------------------

END MODULE defs



