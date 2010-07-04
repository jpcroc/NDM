
subroutine read_parameters()

  use defs
  use lanczos_defs
  use saddles
  use random_art

  implicit none  
  integer :: ierror
  real(8) :: ran3

  integer, dimension(8) :: value
  character(8) :: date
  character(10) :: time
  character(5) :: zone
  character(len=40) :: temporary 
  character(len=10) :: eventtype

  ! We first read the parameters defining the run


  call getenv('EVENT_TYPE', temporary)
  if (temporary .eq. '') then
     new_event = .true.
     eventtype = 'NEW'
  else if (temporary .eq. 'NEW') then
     new_event = .true.
     eventtype = 'NEW'
  else if (temporary .eq. 'REFINE_SADDLE') then
     new_event = .false.
     eventtype = 'REFINE_SADDLE'
  else 
     write(*,*) 'Error: event_types permitted are NEW and REFINE_SADDLE'
     stop
  endif

  call getenv('Temperature', temporary)
  if (temporary .eq. '') then
     write(*,*) 'Error: Metropolis temperature is not defined'
     stop
  else
     read(temporary,*) temperature
  endif

  call getenv('NATOMS', temporary)
  if (temporary .eq. '') then
     write(*,*) 'Error: NATOMS is not defined'
     stop
  else
     read(temporary,*) natoms
  endif

  ! Initalise the random number generator
  call getenv('Random_seed', temporary)
  if (temporary .eq. '') then
    call date_and_time(date,time,zone,value)
    idum = -1 * mod( (1000 * value(7) + value(8)), 1024)
  else
     read(temporary,*) idum
  endif

!  call getenv('MAXNEI', temporary)
!  if (temporary .eq. '') then
!     maxnei = natoms
!  else
!     read(temporary,*) maxnei
!  endif

  call getenv('MAX_NUMBER_EVENTS', temporary)
  if (temporary .eq. '') then
     number_events = 100
  else
     read(temporary,*) number_events
  endif



  ! File names
  call getenv('LOGFILE', temporary)
  if (temporary .eq. '') then
     LOGFILE   = 'log.file'
  else
     read(temporary,*) LOGFILE
  endif

  call getenv('EVENTSLIST', temporary)
  if (temporary .eq. '') then
     EVENTSLIST   = 'events.list'
  else
     read(temporary,*) EVENTSLIST
  endif

  call getenv('REFCONFIG', temporary)
  if (temporary .eq. '') then
     REFCONFIG   = 'refconfig.dat'
  else
     read(temporary,*) REFCONFIG
  endif

  call getenv('FINAL', temporary)
  if (temporary .eq. '') then
     FINAL  = 'min'
  else
     read(temporary,*) FINAL
  endif

  call getenv('SADDLE', temporary)
  if (temporary .eq. '') then
     SADDLE  = 'sad'
  else
     read(temporary,*) SADDLE
  endif

  call getenv('FILECOUNTER', temporary)
  if (temporary .eq. '') then
     COUNTER  = 'filecounter'
  else
     read(temporary,*) COUNTER
  endif

  call getenv('RESTART_FILE', temporary)
  if (temporary .eq. '') then
     RESTARTFILE = 'restart.dat'
  else
     read(temporary,*) RESTARTFILE
  endif

  call getenv('Print_details', temporary)
  if (temporary .eq. '') then
     print_details = .true.
  else
     read(temporary,*) print_details
  endif

  call getenv('kprint', temporary)
  if (temporary .eq. '') then
     kprint = 5
  else
     read(temporary,*) kprint
  endif

  call getenv('mprint', temporary)
  if (temporary .eq. '') then
     mprint = 100
  else
     read(temporary,*) mprint
  endif

  call getenv('iprint', temporary)
  if (temporary .eq. '') then
     iprint = 10
  else
     read(temporary,*) iprint
  endif


  call getenv('bulk_structure', temporary)
  if (temporary .eq. '') then
     bulk_structure='bcc'
  else
     read(temporary,*) bulk_structure
  endif

  call getenv('number_defects', temporary)
  if (temporary .eq. '') then
     art_ninst= 1
  else
     read(temporary,*) art_ninst
  endif


  ! Reading details for activation
  ! Read type of events - local or global
  call getenv('Type_of_Events', TYPE_EVENTS)
  if ( (TYPE_EVENTS .ne. 'global') .and.&
       (TYPE_EVENTS .ne. 'local' ) .and.&
       (TYPE_EVENTS .ne. 'energy' ) .and.&
       (TYPE_EVENTS .ne. 'defect') ) then
     write(*,*) 'Error : only global or local or defect type of events are accepted - provided: ', TYPE_EVENTS
     stop
  endif

  if ( (TYPE_EVENTS .eq. 'local')  .or.  &
       (TYPE_EVENTS .eq. 'energy' ) .or.&
       (TYPE_EVENTS .eq. 'defect') ) then
     call getenv('Radius_Initial_Deformation', temporary)
     if (temporary .eq. '') then
        write(*,*) 'Error: Radius_Initial_Deformation must be defined when TYPE_EVENTS is local'
        stop
     else
        read(temporary,*) LOCAL_CUTOFF
     endif

     call getenv('Central_Atom',temporary)
     if (temporary .eq. '') then
       preferred_atom = -1
     else
       read(temporary,*) preferred_atom
     endif
  endif

  ! Info regarding initial displacement
  call getenv('Initial_Step_Size', temporary)
  if (temporary .eq. '') then
     INITSTEPSIZE = 0.001             ! Size of initial displacement in Ang.
  else
     read(temporary,*) INITSTEPSIZE
  endif

  ! Info regarding energy difference to classify as possible initial move center
  call getenv('energy_difference', temporary)
  if (temporary .eq. '') then
     ENERGYDIFFERENCE = 0.1             ! energy difference in eV
  else
     read(temporary,*) ENERGYDIFFERENCE
  endif

  ! Minimum number of steps in kter-loop before we call lanczos
  call getenv('Min_Number_KSteps',temporary)
  if (temporary .eq. '') then
     KTER_MIN = 1
  else
     read(temporary,*) KTER_MIN
  endif

  ! Maximum number of iteration in the activation
  call getenv('Max_Iter_Activation',temporary)
  if (temporary .eq. '') then
     MAXITER = 1000
  else
     read(temporary,*) MAXITER
  endif

  ! Maximum number of iteration in the activation
  call getenv('Max_Iter_Basin',temporary)
  if (temporary .eq. '') then
     MAXKTER = 100
  else
     read(temporary,*) MAXKTER
  endif


  ! Number of relaxation perpendicular moves - basin and activation
  call getenv('Max_Perp_Moves_Basin',temporary)
  if (temporary .eq. '') then
     MAXKPERP = 2
  else
     read(temporary,*) MAXKPERP
  endif

  call getenv('Max_Perp_Moves_Activ',temporary)
  if (temporary .eq. '') then
     MAXIPERP = 12
  else
     read(temporary,*) MAXIPERP
  endif

  ! Increment size - overall scaling (in angstroems)
  call getenv('Increment_Size',temporary)
  if (temporary .eq. '') then
     INCREMENT = 0.01
  else
     read(temporary,*) INCREMENT
  endif

  ! Eigenvalue threshold
  call getenv('Eigenvalue_Threshold',temporary)
  if (temporary .eq. '') then
     write(*,*) 'Error : No eigenvalue threshold provided  (Eigenvalue_Threshold)'
     stop
  else
     read(temporary,*) EIGEN_THRESH
  endif

  ! Force threshold for the perpendicular relaxation
  call getenv('Force_Threshold_Perp_Rel',temporary)
  if (temporary .eq. '') then
     FTHRESHOLD = 1.0
  else
     read(temporary,*) FTHRESHOLD
  endif
  FTHRESH2 = FTHRESHOLD * FTHRESHOLD

  ! Force threshold for the convergence at the saddle point
  call getenv('Exit_Force_Threshold',temporary)
  if (temporary .eq. '') then
     EXITTHRESH = 1.0
  else
     read(temporary,*) EXITTHRESH
  endif

  ! Force threshold for the convergence at the saddle point
  call getenv('Number_Lanczos_Vectors',temporary)
  if (temporary .eq. '') then
     NVECTOR_LANCZOS = 16
  else
     read(temporary,*) NVECTOR_LANCZOS
  endif

  call getenv('Lanczos_step',temporary)
  if (temporary .eq. '') then
     lanczos_step = 0.001
  else
     read(temporary,*) lanczos_step
  endif

  ! We set up other related parameters
  VECSIZE = 3*natoms

  ! And we allocate the vectors
  allocate(type(natoms))
  allocate(force(VECSIZE))
  allocate(pos(VECSIZE))
  allocate(posref(VECSIZE))
  allocate(direction_restart(VECSIZE))
  allocate(initial_direction(VECSIZE))
  allocate(Atom(VECSIZE))

  ! Vectors for lanczos
  allocate(old_projection(VECSIZE))
  allocate(projection(VECSIZE))
  allocate(first_projection(VECSIZE))

  x => pos(1:NATOMS)
  y => pos(NATOMS+1:2*NATOMS)
  z => pos(2*NATOMS+1:3*NATOMS)

  xref => posref(1:NATOMS)
  yref => posref(NATOMS+1:2*NATOMS)
  zref => posref(2*NATOMS+1:3*NATOMS)

  fx => force(1:NATOMS)
  fy => force(NATOMS+1:2*NATOMS)
  fz => force(2*NATOMS+1:3*NATOMS)


  ! We write down the various parameters for the simulation
  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='rewind',iostat=ierror)  
  write(FLOG,*)            '*************************'
  write(FLOG,'(A15,F8.3)') ' Version     : ', version
  write(FLOG,*)            '*************************'
  write(flog,'(A39,A16  )')  ' Event type                           : ', eventtype
  write(flog,'(A39,f16.4)')  ' Temperature                          : ', temperature
  write(flog,'(A39,I16  )')  ' Number of atoms                      : ', natoms     
  write(flog,'(A39,I16  )')  ' Number of events                     : ', number_events
!  write(flog,'(A39,I16  )')  ' Maximum number of neighbours         : ', maxnei    
  write(flog,'(A39,I16  )')  ' Initial value for the seed           : ',idum

  write(flog,*) ' '
  write(flog,*) 'Selection of the event '
  write(flog,*) '********************** '
  write(flog,'(1X,A39,A12)')   ' - Type of events                    : ', TYPE_EVENTS
  if (TYPE_EVENTS .eq. 'local') then 
    write(flog,'(1X,A39,F12.4)')   ' - Radius of deformation (local ev.) : ', LOCAL_CUTOFF
    if (preferred_atom .gt. 0) then 
      write(flog,'(1X,A39,I12)')   ' - Central atom for events           : ', preferred_atom
    else
      write(flog,'(1X,A39,A12)')   ' - Central atom for events           :   none (all equal) '
    endif
  endif
  write(flog,*) ' '
  write(flog,*) 'Activation parameters '
  write(flog,*) '********************* '
  write(flog,'(1X,A39,F12.4)') ' - Eigenvalue threshold              : ', EIGEN_THRESH
  write(flog,'(1X,A39,F12.4)') ' - Total force threshold (saddle)    : ', EXITTHRESH

  write(flog,'(1X,A39,F12.4)') ' - Initial step size                 : ', INITSTEPSIZE
  write(flog,'(1X,A39,F12.4)') ' - Increment size                    : ', INCREMENT
  write(flog,'(1X,A51,I8)')    ' - Min. number of ksteps before calling lanczos  : ', KTER_MIN
  write(flog,'(1X,A51,I8)')    ' - Maximum number of iteractions (basin -kter)   : ', MAXKTER
  write(flog,'(1X,A51,I8)')    ' - Maximum number of iteractions (activation)    : ', MAXITER
  write(flog,'(1X,A51,I8)')    ' - Maximum number of perpendicular moves (basin) : ', MAXKPERP
  write(flog,'(1X,A51,I8)')    ' - Maximum number of perpendicular moves (activ) : ', MAXIPERP
  write(flog,'(1X,A51,F8.4)')  ' - Force threshold for perpendicular relaxation  : ', FTHRESHOLD

  write(flog,*) ' '
  write(flog,'(1X,A51,I8   )') ' - Number of vectors computed by Lanzcos         : ', NVECTOR_LANCZOS
  write(flog,'(1X,A51,F8.4)')  ' - Lanczos step (for numerical derivative)       : ', lanczos_step

  if (print_details) then 
     write(flog,*) ' '
     write(flog,*) 'Printing of details '
     write(flog,*) '********************* '
     write(flog,'(A39,A16  )')  ' Printing details                     : ',print_details
     write(flog,'(A39,I16  )')  ' Printing inverval - leaving the basin: ',kprint
     write(flog,'(A39,I16  )')  ' Printing inverval - activating       : ',iprint
     write(flog,'(A39,I16  )')  ' Printing inverval - minimization     : ',mprint
  endif


  write(flog,*) ' '
  write(flog,*) 'Input / Output '
  write(flog,*) '********************* '
  write(flog,'(A39,A16  )')  ' Name of log file                     : ', logfile
  write(flog,'(A39,A16  )')  ' Liste of events                      : ', eventslist 
  write(flog,'(A39,A16  )')  ' Reference configuration              : ', refconfig   
  write(flog,'(A39,A16  )')  ' Restart file                         : ', restartfile
  write(flog,'(A39,A16  )')  ' Prefix for minima (file)             : ', FINAL      
  write(flog,'(A39,A16  )')  ' Prefix for saddle points             : ', SADDLE
  write(flog,'(A39,A16  )')  ' File with filecounter                : ', counter

  close(flog)
end subroutine
