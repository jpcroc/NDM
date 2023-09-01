!> @file
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> ART read_parameters
!! Read the parameters defining the simulation

! $Revision$
! $Date: 2015-12-14 11:01:46 -0500 (Lun, 14 déc 2015) $
! $Id: read_parameters.f90 1534 2015-12-14 15:59:38Z mickael $
module read_parameters_mod
use ndm2art2ndm,only:parapath
use arret_ndm_mod,only:arret_ndm
use newunit_mod,only:newunit
use gen_com_m,only:lprteat

contains


  subroutine read_parameters( )

    use defs
    use lanczos_defs
    use saddles
    use generate_local_region_mod,only: 
    !nat_inner, nat_outer, nat_local
    !integer, dimension(:), allocatable :: inner_list, outer_list
    !real(8)                            :: r2_inner, r2_outer
    !real(8), dimension(:), allocatable :: mask_inner_region
    implicit none
    character(len=20) :: TYPE_of_EVENTS
    !Local variables
    integer::unitart,ierr
    Character(len=40)  :: temporary
    integer :: Number_Lanczos_Vectors
    integer :: activation_maxiter,Max_Perp_Moves_Basin,Min_Number_KSteps,Max_Iter_Basin,Lanczos_SCLoop,&
         &Number_Lanczos_Vectors_h,Number_Lanczos_Vectors_c,Max_Perp_Moves_Activ,MAX_DIIS,&
         &Type_selected,max_number_events
    logical :: exists_already,calc_of_projection,Lanczos_of_minimum,DIIS_Check_Eigenvector,Clean_wavefunct,Dual_system
    character(len=20):: event_type
    real*8::initial_step_size,increment_size,Exit_Force_Threshold,Force_Threshold_Perp_Rel,Eigenvalue_Threshold,Lanczos_collinear,&
         &delta_disp_Lanczos,Prefactor_Push_Over_Saddle,Relative_To_Second_Order_Error,&
         &delta_threshold,delr_threshold,DIIS_Step_size
    integer::MAX_REL_STEPS_PERP_FIRE,MAX_REL_STEPS_INCREMENT_FIRE,MAX_RELAXATION_STEP_FIRE

  character(len=11) :: FILECOUNTER,RESTART_FILE

    namelist /inputart/setup_initial,guess_noise,GUESSFILE,dim,check_connectivity,local_lanczos,INNER_REGION,&
         & OUTER_REGION ,NORM_CRITERIUM_FIRE,FMAX_CRITERIUM_FIRE,dT_MAX_FIRE, MAX_RELAXATION_STEP_FIRE,MAX_REL_STEPS_PERP_FIRE,&
         &MAX_REL_STEPS_INCREMENT_FIRE,&
         & MAX_LANCZOS_STEPS,temperature,max_Number_Events, TYPE_of_EVENTS,Radius_Initial_Deformation,Central_Atom,sym_break_dist,&
         & initial_step_size,INCREMENT_size,Exit_Force_Threshold,Force_Threshold_Perp_Rel,&
         &BASIN_FACTOR,Max_Perp_Moves_Basin, SMOOTH_DIR_CHANGE&
         &, Min_Number_KSteps,Eigenvalue_Threshold,Max_Iter_Basin ,calc_of_projection,Lanczos_of_minimum,Lanczos_SCLoop,&
         &Lanczos_collinear,Number_Lanczos_Vectors,Number_Lanczos_Vectors_H,Number_Lanczos_Vectors_C,&
         &delta_disp_Lanczos,Max_Perp_Moves_Activ, Prefactor_Push_Over_Saddle, Relative_To_Second_Order_Error,delta_threshold,&
         &delr_threshold,USE_DIIS,ITERATIVE ,ivisuart,&
         &INFLECTION,DIIS_FORCE_THRESHOLD,DIIS_MEMORY, DIIS_Step_size,FACTOR_DIIS,max_diis,DIIS_Check_Eigenvector,&
         & fileCOUNTER, NPART_DR_THRESHOLD,LOGFILE,EVENTSLIST,SAVE_CONF_INT, write_restart_file,WRITE_REJECTED_EVENT,&
         &RESTART_FILE , write_xyz ,REFCONFIG,FINAL,SADDLE,CHECK,Coord_length,coord_number,&
         &Type_selected,Dual_system,size_system,event_type,activation_maxiter
    size_system=-1.
    Dual_system=.false.
    Type_selected=0
    coord_number=-1
    coord_length=-1.
    ivisuart=40
    maxnei = 1
    clean_wf = .false.
    CHECK  = 'con'
    SADDLE  = 'sad'
    FINAL  = 'min'
    REFCONFIG   = 'refconfig.dat'
    write_xyz = .false.
    RESTART_FILE = 'restart.dat'
    WRITE_REJECTED_EVENT = .true.
    write_restart_file = .false.
    SAVE_CONF_INT = .false.
    EVENTSLIST   = 'events.list'
    LOGFILE   = 'log.file'
    fileCOUNTER  = 'filecounter'
    NPART_DR_THRESHOLD=0.1
    DIIS_Check_Eigenvector = .true.
    FMAX_CRITERIUM_FIRE = 0.020d0
    max_diis = 100
    FACTOR_DIIS=5.0
    DIIS_Step_size = 0.01d0
    DIIS_MEMORY=12
    DIIS_FORCE_THRESHOLD = 1.0d0
    INFLECTION=-1
    ITERATIVE = .false.
    USE_DIIS = .false.
    delr_threshold = 0.0d0
    delta_threshold = 0.0d0
    Relative_To_Second_Order_Error = 1d+10
    Prefactor_Push_Over_Saddle = 0.15
    Max_Perp_Moves_Activ=12
delta_disp_Lanczos=0.01
    Number_Lanczos_Vectors = 16
    Number_Lanczos_Vectors_H = -1
    Number_Lanczos_Vectors_C = -1

    Lanczos_collinear= 0.7d0
    Lanczos_SCLoop = 1
    Lanczos_of_minimum = .False.
    calc_of_projection = .false.
    Max_Iter_Basin = 100
    Eigenvalue_Threshold=-1000.0
    Min_Number_KSteps=1
    SMOOTH_DIR_CHANGE = 3
    Max_Perp_Moves_Basin=2
    BASIN_FACTOR=1.
    Force_Threshold_Perp_Rel=-1
    Exit_Force_Threshold=0.1 
    INCREMENT_size=0.01
    initial_step_size=0.001
    activation_maxiter=100
    sym_break_dist=0
    Central_Atom=-1
    TYPE_of_EVENTS='undefined'
    Radius_Initial_Deformation=-1
    max_Number_Events=100
    temperature=-1.
    MAX_LANCZOS_STEPS = 280
    MAX_REL_STEPS_PERP_FIRE = 25
    MAX_REL_STEPS_INCREMENT_FIRE = 3
    MAX_RELAXATION_STEP_FIRE = 500
    DT_MAX_FIRE = 0.15d0
    FMAX_CRITERIUM_FIRE = 0.020d0
    NORM_CRITERIUM_FIRE = 0.008d0
    setup_initial=.false.
    eventtype='NEW'
    event_type='NEW'
    guess_noise=0
    GUESSFILE = 'initdir'
    dim=3
    CHECK_CONNECTIVITY = .false.
    LOCAL_LANCZOS = .false.
    local_force=.false. ! IS NOT in the namelist
    global_convergence=.true. ! not in the namelist
    inner_region=6.0
    OUTER_REGION = 6.0

    
    call newunit(unitart)
    open(unit=unitart, file='art.in', status='unknown')
    read (unitart, nml=inputart)

    !SECTION____________________________ ATOMS  ( Obsolete )

!!$  call getenv('NATOMS', temporary)
!!$  if (temporary .eq. '') then
!!$     write(unit6P,*) 'Error: NATOMS is not defined'
!!$     stop
!!$  else
!!$     read(temporary,*) natoms
!!$  end if
    !CRC natoms set to atart%im in ndm2art
    !!__________________
    ! We now get the types - define up to 5
!!$  call getenv('type1',temporary)
!!$  if (temporary .eq. '') then
!!$     write(unit6P,*) 'Error, must at least define 1 type of atoms -  use "setenv type1 Si", for example'
!!$     stop
!!$  else
!!$     read(temporary,*) type_name(1)
!!$  end if
!!$
!!$  call getenv('type2',temporary)
!!$  if (temporary .eq. '') then
!!$     type_name(2) = ''
!!$  else
!!$     read(temporary,*) type_name(2)
!!$  end if
!!$
!!$  call getenv('type3',temporary)
!!$  if (temporary .eq. '') then
!!$     type_name(3) = ''
!!$  else
!!$     read(temporary,*) type_name(3)
!!$  end if
!!$
!!$  call getenv('type4',temporary)
!!$  if (temporary .eq. '') then
!!$     type_name(4) = ''
!!$  else
!!$     read(temporary,*) type_name(4)
!!$  end if
!!$
!!$  call getenv('type5',temporary)
!!$  if (temporary .eq. '') then
!!$     type_name(5) = ''
!!$  else
!!$     read(temporary,*) type_name(5)
!!$  end if
!!$
!!$  call getenv('type6',temporary)
!!$  if (temporary .ne. '') then
!!$     write(unit6P,*) 'Error: The code can only handle 5 atomic types, change read_parameters.f90, to allow for more.'
!!$     stop
!!$  end if
    !crc type_name is used only to define Atom(i) (in retart !) which is the character type of the atom, useless in NDM

    !SECTION_____________________________ ART

    ! Hide Option, by default .false.
    ! Only Lanczos analysis for the minimum and the inflection point.



!!$  call getenv('Setup_Initial',temporary)
!!$  if (temporary == '') then
!!$     setup_initial = .false.
!!$  else
!!$     read(temporary,*)  setup_initial
!!$  end if

    if (event_type.ne.'NEW') then
       write(unit6P,*)'event_type<>new not coded'
       stop
    end if

!!$  if (eventtype=='undefined') then ! EVENTYPE ALWAYS NEW for Emp Pot
!!$     
!!$  !!__________________
!!$!  call getenv('EVENT_TYPE', temporary)
!!$  selectcase( eventtype )
!!$     case ( 'NEW' )
          new_event = .true.
          eventtype = 'NEW'
!!$     case ( 'REFINE_SADDLE' )
!!$          new_event = .false.
!!$          eventtype = 'REFINE_SADDLE'
!!$     case ( 'REFINE_AND_RELAX' )
!!$          new_event = .false.
!!$          eventtype = 'REFINE_AND_RELAX'
!!$     case default
!!$          write(unit6P,*) 'Error: eventtypes permitted are:'
!!$          write(unit6P,*) 'NEW, REFINE_SADDLE, REFINE_AND_RELAX'
!!$          stop
!!$  end select

    !!__________________
    !  if ( eventtype == 'GUESS_DIRECTION' ) then

    !!__________________HIDE
    !     call getenv('Noise',temporary)
!!$     if (temporary == '') then
!!$        guess_noise = 0.0d0
!!$     else
!!$        read(temporary,*) guess_noise
!!$     end if
    !!__________________NAME OF FILE
    !     call getenv('INITDIR', temporary)
    !     if (temporary .eq. '') then
    !        GUESSFILE = 'initdir'
    !     else
    !        read(temporary,*) GUESSFILE
    !     end if
    GUESSFILE = adjustl(GUESSFILE)

    !  end if

    !!__________________
    ! type of energy and force calculation
!!$  call getenv('ENERGY_CALC', temporary) !energy_type à supprimer
!!$  if (temporary .eq. '') then
!!$     write(unit6P,*) "Error: energy calculation type is not defined: ENERGY_CALC "
!!$     write(unit6P,*) " choose: DFT, SWP (Stillinger-Weber Si 'dia'), SWA (Stillinger-Weber amorphous Si) or LAM (LAMMPS) "
!!$     stop
!!$  else
!!$     read(temporary,*) energy_type
!!$  endif

    !!__________________
    ! spatial dimension
!!$  call getenv('DIM', temporary)
!!$  if (temporary .eq. '') then
!!$     DIM = 3
!!$  else
!!$     read(temporary,*) DIM
!!$  endif

    !!_________________
    ! Conversion of energy units
!!$  call getenv('UNITS_CONVERSION',temporary)
!!$  if (temporary .eq. '' .or. temporary .eq. 'metal') then
!!$     UNITS_CONVERSION = temporary
!!$     energy_conversion = 1.0d0
!!$     position_conversion = 1.0d0
!!$  else if (temporary == 'real') then
!!$     UNITS_CONVERSION = temporary
!!$     energy_conversion =   4.3363d-2    !=1.0d0/23.0609d0  ! Unit eV/(kcal/mol)
!!$     position_conversion = 1.0d0
!!$  else if (temporary == 'electron') then
!!$     UNITS_CONVERSION = temporary
!!$     energy_conversion = 27.2107d0      ! Unit eV/(Ha)
!!$     position_conversion = 1.889726d0   ! Unit Bohr/Angstrom
!!$  else if (temporary == 'si') then
!!$     UNITS_CONVERSION = temporary
!!$     energy_conversion   = 6.24181d18   ! Unit eV/J
!!$     position_conversion = 1.0d-10      ! Unit m/Angstrom
!!$  else if(temporary == 'kcal2ev') then
!!$     UNITS_CONVERSION = temporary
!!$     energy_conversion = 1.0d0/23.0609d0
!!$     position_conversion = 1.0d0
!!$  else
!!$     write(unit6P,*) "Error: UNITS_CONVERSION is only defined for none or metal or real or electron or si"
!!$     write(unit6P,*) "You must modify read_param.f90 for more choices"
!!$     write(unit6P,*) "or use none as argument"
!!$     stop
!!$  endif

!!$  ! Use a local force calculation for a first evaluation of barriers
!!$  call getenv ('LOCAL_FORCE',temporary)
!!$  if (temporary == '') then
!!$     LOCAL_FORCE = .false.
!!$  else
!!$     read(temporary,*) LOCAL_FORCE
!!$  endif
!!$
!!$  ! Use global convergence after local convergence at saddle points 
!!$  call getenv ('GLOBAL_CONVERGENCE',temporary)
!!$  if (temporary == '') then
!!$     GLOBAL_CONVERGENCE = .true.
!!$  else
!!$     if (LOCAL_FORCE) then
!!$        read(temporary,*) GLOBAL_CONVERGENCE
!!$     else
!!$        GLOBAL_CONVERGENCE = .true.
!!$     endif
!!$  endif
!!$  

    ! Use check connectivity to validate the reversibility of events
!!$  call getenv ('CHECK_CONNECTIVITY',temporary)
!!$  if (temporary == '') then
!!$     CHECK_CONNECTIVITY = .false.
!!$  else
!!$     read(temporary,*) CHECK_CONNECTIVITY
!!$  endif

    ! Use a local calculation for LANCZOS, useful for very large systems
    if (LOCAL_FORCE .and. .not. LOCAL_LANCZOS) then
       write(unit6P,*) 'Error: when LOCAL_FORCES is true, LOCAL_LANCZOS should be true also.'
       write(unit6P,*) 'Error: The code will automatically set LOCAL_LANCZOS to true.'
       LOCAL_LANCZOS = .true.
    endif


!!$  if (LOCAL_FORCE .or. LOCAL_LANCZOS) then
!!$     call getenv('INNER_REGION_RADIUS',temporary)
!!$     if (temporary == '') then
!!$        INNER_REGION = 6.0
!!$     else
!!$        read(temporary,*) INNER_REGION
!!$     endif
!!$  endif
!!$
!!$  if (LOCAL_FORCE .or. LOCAL_LANCZOS) then
!!$     call getenv('OUTER_REGION_WIDTH',temporary)
!!$     if (temporary == '') then
!!$        OUTER_REGION = 6.0
!!$     else
!!$        read(temporary,*) OUTER_REGION
!!$     endif
!!$  endif



    !!_________________ FIRE minimisation

    ! Norm 2 criterium for stopping minimization
!!$    call getenv('NORM_CRITERIUM_FIRE', temporary)
!!$    if (temporary .eq. '') then
!!$       NORM_CRITERIUM_FIRE = 0.008d0
!!$    else
!!$       read(temporary,*) NORM_CRITERIUM_FIRE
!!$    endif
!!$
!!$    ! Norm inf criterium for stopping minimization
!!$    call getenv('FMAX_CRITERIUM_FIRE', temporary)
!!$    if (temporary .eq. '') then
!!$       FMAX_CRITERIUM_FIRE = 0.020d0
!!$    else
!!$       read(temporary,*) FMAX_CRITERIUM_FIRE
!!$    endif
!!$
!!$    ! Maximum time step
!!$    call getenv('DT_MAX_FIRE', temporary)
!!$    if (temporary .eq. '') then
!!$       DT_MAX_FIRE = 0.15d0
!!$    else
!!$       read(temporary,*) DT_MAX_FIRE
!!$    endif
!!$
!!$    ! Maximum number of relaxation step during global minimization
!!$    call getenv('MAX_RELAXATION_STEP_FIRE', temporary)
!!$    if (temporary .eq. '') then
       MAX_ITER_FIRE = MAX_RELAXATION_STEP_FIRE
!!$    else
!!$       read(temporary,*) MAX_ITER_FIRE
!!$    endif
!!$
!!$    ! Maximum number of relaxation steps in the perp direction during activation
!!$    call getenv('MAX_REL_STEPS_PERP_FIRE', temporary)
!!$    if (temporary .eq. '') then
!!$       MAX_PERP_ITER_FIRE = 25
!!$    else
!!$       read(temporary,*) MAX_PERP_ITER_FIRE
!!$    endif
!!$
     MAX_PERP_ITER_FIRE =MAX_REL_STEPS_PERP_FIRE
!!$    ! Increment in the number of relaxation step during activation
!!$    call getenv('MAX_REL_STEPS_INCREMENT_FIRE', temporary)
!!$    if (temporary .eq. '') then
       MAX_PERP_INCR_FIRE = MAX_REL_STEPS_INCREMENT_FIRE
!!$    else
!!$       read(temporary,*) MAX_PERP_INCR_FIRE
!!$    endif
!!$
!!$    ! Maximum number of Lanczos steps during activation
!!$    call getenv('MAX_LANCZOS_STEPS', temporary)
!!$    if (temporary .eq. '') then
!!$       MAX_LANCZOS_STEPS = 280
!!$    else
!!$       read(temporary,*) MAX_LANCZOS_STEPS
!!$    endif
!!$
!!$

    !!-----------------
    ! name of the input file from lammps
!!$  call getenv('INPUT_LAMMPS_FILE', temporary)
!!$  if (temporary .eq. '' .and. energy_type == 'LAM') then
!!$     inquire(file = 'in.lammps', exist = exists_already)
!!$     if (.not. exists_already) then
!!$        write(unit6P,*) "ERROR, INPUT_LAMMPS_FILE missing and default 'in.lammps' not present."
!!$        write(unit6P,*) "The program will stop"
!!$        stop
!!$     end if
!!$     INPUT_LAMMPS_FILE  = 'in.lammps'
!!$  else if (energy_type == 'LAM') then
!!$     read(temporary,*) INPUT_LAMMPS_FILE
!!$     inquire(file = INPUT_LAMMPS_FILE, exist = exists_already)
!!$     if (.not. exists_already) then
!!$        write(unit6P,*) "ERROR, INPUT_LAMMPS_FILE with name ", temporary," is missing."
!!$        write(unit6P,*) "The program will stop"
!!$        stop
!!$     end if
!!$  endif

    !sts/__________________
    ! Fictive temperature, if negative always reject the event
    !  call getenv('Temperature', temporary)
!!$    if (temperature.eq.-1.) then
!!$       if (rang==0) write(unit6P,*) 'Error: Metropolis temperature is not defined'
!!$       stop
!!$  else
!!$     read(temporary,*) temperature
!!$    end if

    !!__________________
    ! Maximum number of events
!!$  call getenv('Max_Number_Events', temporary)
!!$  if (temporary .eq. '') then
!!$     number_events = 100
!!$  else
!!$     read(temporary,*) number_events
!!$  end if
    !CROCTODO  if (number_events.lt.nb de taches para) the stop
    number_events=max_number_events


    
    !!__________________
    ! Read type of events
    ! Activation: global, local, list_local, list, local_coord
    !  call getenv('Type_of_Events', TYPE_EVENTS)
    type_events=type_of_events
    if ( (TYPE_EVENTS .ne. 'energy') .and. (TYPE_EVENTS .ne. 'global') &
         &.and. (TYPE_EVENTS .ne. 'local') .and. &
         & (TYPE_EVENTS .ne. 'list')   .and. (TYPE_EVENTS .ne. 'list_local') .and. &
         & (TYPE_EVENTS .ne. 'local_coord') ) then
       write(unit6P,*) 'Error : wrong type of  events : ',&
            & TYPE_EVENTS
       stop
    end if

    if (TYPE_EVENTS .eq. 'global' .or. TYPE_EVENTS .eq. 'list' ) then
!!$     if (LOCAL_FORCE) then
!!$        write(unit6P,*) 'Error: EVENTS MUST BE LOCAL or LIST_LOCAL WHEN USING LOCAL_FORCE'
!!$        stop 
!!$     endif
    endif
    !!________________
    local_cutoff=Radius_Initial_Deformation
    ! Cutoff for local and list_local (in angstroems)
    if (TYPE_EVENTS .eq. 'local' .or. TYPE_EVENTS .eq. 'list_local' ) then
       !     call getenv('Radius_Initial_Deformation', temporary)
       if (local_cutoff==-1) then
          write(unit6P,*) 'Error: local_cutoffmust be defined when TYPE_EVENTS is local'
          stop
!!$     else
!!$        read(temporary,*) LOCAL_CUTOFF
       end if
    end if
    if (TYPE_EVENTS .eq. 'energy' ) then
       if (.not.lprteat) then
          write(unit6P,*)'event type energy: set lprteat=.true. in *.din'
          stop
       end if
    end if

    !!__________________
    ! Number of the atom around which the initial move takes place
  if (TYPE_EVENTS .eq. 'local' ) then
!!$     call getenv('Central_Atom',temporary)
!!$     if (temporary .eq. '') then
        preferred_atom = Central_Atom
!!$     else
!!$        read(temporary,*) preferred_atom
!!$     end if
  end if

    !!__________________
    ! Breaks the symmetry of the crystal by randomly displacing
    ! all atoms by this distance
!!$  call getenv('sym_break_dist',temporary)
!!$  if (temporary .eq. '') then
!!$     sym_break_dist = 0.0
!!$  else
!!$     read(temporary,*) sym_break_dist
!!$  end if

    !!__________________
    ! Maximum number of iteraction for reaching the saddle point
!!$  call getenv('Activation_MaxIter',temporary)
!!$  if (temporary .eq. '') then
!!$     MAXPAS = 100
!!$  else
!!$     read(temporary,*) MAXPAS
!!$  end if
    maxpas=activation_maxiter
    !!__________________
    ! Info regarding initial displacement
!!$  call getenv('Initial_Step_Size', temporary)
!!$  if (temporary .eq. '') then
!!$     INITSTEPSIZE = 0.001             ! Size of initial displacement in Ang.
!!$  else
!!$     read(temporary,*) INITSTEPSIZE
!!$  end if
!!$
    INITSTEPSIZE=initial_step_size
!!$  !!__________________
!!$  ! Increment size - overall scaling (in angstroems)
!!$  call getenv('Increment_Size',temporary)
!!$  if (temporary .eq. '') then
!!$     INCREMENT = 0.01
!!$  else
!!$     read(temporary,*) INCREMENT
!!$  end if
    increment=increment_size
    !!__________________
    ! Force threshold for the convergence at the saddle point
!!$  call getenv('Exit_Force_Threshold',temporary)
!!$  if (temporary .eq. '') then
!!$     EXITTHRESH = 0.1
!!$  else
!!$     read(temporary,*) EXITTHRESH
!!$  end if
    EXITTHRESH=Exit_Force_Threshold
    !!__________________
    ! Force threshold for the perpendicular relaxation
    !  call getenv('Force_Threshold_Perp_Rel',temporary)
    FTHRESHOLD=Force_Threshold_Perp_Rel
    if (FTHRESHOLD==-1) then
       FTHRESHOLD = EXITTHRESH !1.0
    else
       !     read(temporary,*) FTHRESHOLD
       if(FTHRESHOLD>EXITTHRESH) then
          FTHRESHOLD = EXITTHRESH
          write(unit6P,*) 'Warning: Force_Threshold_Perp_Rel should be smaller or equal to Exit_Force_Threshold'
          write(unit6P,*) 'The code has automatically set Force_Threshold_Perp_Rel=Exit_Force_Threshold'
       endif
    end if

    !SECTION_____________________________ HARMONIC WELL

    ! Size of the parallel displacement in for leavign the basin with respect to
    ! the overall increment
!!$  call getenv('Basin_Factor',temporary)
!!$  if (temporary .eq. '') then
!!$     BASIN_FACTOR = 1.0d0
!!$  else
!!$     read(temporary,*) BASIN_FACTOR
!!$  end if

    !!__________________
    ! Maximum number of perpendicular steps leaving basin
!!$  call getenv('Max_Perp_Moves_Basin',temporary)
!!$  if (temporary .eq. '') then
!!$     MAXKPERP = 2
!!$  else
!!$     read(temporary,*) MAXKPERP
!!$  end if
    maxkperp=Max_Perp_Moves_Basin
    ! NUmber of steps for smoothing the transition from the harmonic basin to the activation regime
!!$  call getenv('Smooth_Dir_Change',temporary)
!!$  if (temporary == '') then
!!$          SMOOTH_DIR_CHANGE = 3
!!$  else
!!$          read(temporary,*) SMOOTH_DIR_CHANGE
!!$  endif

    !!__________________
    ! Minimum number of steps in kter-loop before we call lanczos
!!$  call getenv('Min_Number_KSteps',temporary)
!!$  if (temporary .eq. '') then
     KTER_MIN = Min_Number_KSteps
!!$  else
!!$     read(temporary,*) KTER_MIN
!!$  end if

    !!__________________
    ! Eigenvalue threshold for leaving basin
    !  call getenv('Eigenvalue_Threshold',temporary)
     eigen_thresh=Eigenvalue_Threshold
     if (EIGEN_THRESH==-1000.) then
       write(unit6P,*) 'Error : No eigenvalue threshold provided  (Eigenvalue_Threshold)'
       stop
    end if
!!$  else
!!$     read(temporary,*) EIGEN_THRESH
!!$    end if

    !!__________________
    ! Maximum number of iteration in the activation
!!$  call getenv('Max_Iter_Basin',temporary)
!!$  if (temporary .eq. '') then
     MAXKTER = Max_Iter_Basin
!!$  else
!!$     read(temporary,*) MAXKTER
!!$  end if

    !SECTION_____________________________ LANCZOS

    !!_________________HIDE
    ! inflection in the eigenvalue
    ! if it is true; it calculates the projection only at every two steps
    ! but after 4 steps above of an inflection in the eigenvalue, and if the
    ! last a1 >0.9d0
!!$  call getenv('calc_of_projection',temporary)
!!$  if (temporary == '') then
     calc_proj = calc_of_projection
!!$  else
!!$     read(temporary,*) calc_proj
!!$  end if

    !!__________________
    ! Calculation of the Hessian for each minimum
!!$  call getenv('Lanczos_of_minimum',temporary)
!!$  if (temporary .eq. '') then
     LANCZOS_MIN = Lanczos_of_minimum
!!$  else
!!$     read(temporary,*)  LANCZOS_MIN
!!$  end if
!!$
!!$  !!__________________
!!$  ! Number of iterations in the lanczos Self consistent loop
!!$  call getenv('Lanczos_SCLoop',temporary)
!!$  if (temporary .eq. '') then
     LANCZOS_SCL = Lanczos_SCLoop
!!$  else
!!$     read(temporary,*) LANCZOS_SCL
!!$  end if

    !!__________________
    ! The convergence criteria in the Lanczos Self consistent loop
!!$  call getenv('Lanczos_collinear',temporary)
!!$  if (temporary .eq. '') then
     collinear_factor= Lanczos_collinear
!!$  else
!!$     read(temporary,*) collinear_factor
!!$  end if
!!$
!!$  !!__________________
!!$  ! Number of vectors included in lanczos procedure in the Harmonic well
!!$  call getenv('Number_Lanczos_Vectors',temporary)
!!$  if (temporary .eq. '') then
     NVECTOR_LANCZOS = Number_Lanczos_Vectors
     NVECTOR_LANCZOS_h = Number_Lanczos_Vectors_h
     NVECTOR_LANCZOS_c = Number_Lanczos_Vectors_c
!!$  else
!!$     read(temporary,*) NVECTOR_LANCZOS
!!$  end if

    if (NVECTOR_LANCZOS_H==-1)NVECTOR_LANCZOS_H = NVECTOR_LANCZOS
    if(NVECTOR_LANCZOS_C ==-1)NVECTOR_LANCZOS_C = NVECTOR_LANCZOS_H

    !!__________________
    ! Number of vectors included in lanczos procedure in the Harmonic well
!!$  call getenv('Number_Lanczos_Vectors_H',temporary)
!!$  if (temporary .eq. '') then
!!$     NVECTOR_LANCZOS_H = NVECTOR_LANCZOS
!!$  else
!!$     read(temporary,*) NVECTOR_LANCZOS_H
!!$  end if
!!$
!!$  !!__________________
!!$  ! Number of vectors included in lanczos procedure in convergence
!!$  call getenv('Number_Lanczos_Vectors_C',temporary)
!!$  if (temporary .eq. '') then
!!$     NVECTOR_LANCZOS_C = NVECTOR_LANCZOS_H
!!$  else
!!$     read(temporary,*) NVECTOR_LANCZOS_C
!!$  end if

    !!__________________
    ! The step of the numerical derivative of forces for the Hessian (in Ang)
!!$  call getenv('delta_disp_Lanczos',temporary)
!!$  if (temporary .eq. '') then
     DEL_LANCZOS =  delta_disp_Lanczos
!!$  else
!!$     read(temporary,*) DEL_LANCZOS
!!$  end if

    !SECTION____________________________ CONVERGENCE

    ! Maximum number of perpendicular steps during activation
!!$  call getenv('Max_Perp_Moves_Activ',temporary)
!!$  if (temporary .eq. '') then
     MAXIPERP = Max_Perp_Moves_Activ
!!$  else
!!$     read(temporary,*) MAXIPERP
!!$  end if

    !!__________________
    ! The prefactor for pushing over the saddle point, fraction of distance from
    ! initial minimum to saddle point
!!$  call getenv('Prefactor_Push_Over_Saddle',temporary)
!!$  if (temporary .eq. '') then
     PUSH_OVER = Prefactor_Push_Over_Saddle
!!$  else
!!$     read(temporary,*) PUSH_OVER
!!$  end if

    !!__________________
    ! Adjust push over saddle, until relative error between energy
    ! variation and its 2nd order approximation becomes lower
    ! than this threshold. When this threshold is set to a high value,
    ! push is performed as soon as it decreases energy. When set to a ultra low
    ! value, the push that minimizes the error while decreasing energy is
    ! selected.
!!$  call getenv('Relative_To_Second_Order_Error',temporary)
!!$  if (temporary .eq. '') then
     PUSH_OVER_RELATIVE_ERROR = Relative_To_Second_Order_Error
!!$  else
!!$     read(temporary,*) PUSH_OVER_RELATIVE_ERROR
!!$  end if

    !!__________________
    ! if delta_e < delta_thr .and. delr < delr_thr => end_activation = .true.
    ! Set up them to zero if you dont want to use these criteria.
!!$  call getenv('delta_threshold',temporary)
!!$  if (temporary == '') then
     delta_thr = delta_threshold
!!$  else
!!$     read(temporary,*) delta_thr
!!$  end if

!!$  call getenv('delr_threshold',temporary)
!!$  if (temporary == '') then
     delr_thr = delr_threshold
!!$  else
!!$     read(temporary,*) delr_thr
!!$  end if

    !SECTION_____________________________ DIIS

    ! Do we use DIIS for refining and converging to saddle
!!$  call getenv('Use_DIIS',temporary)
!!$  if (temporary .eq. '') then
!!$     USE_DIIS = .false.
!!$  else
!!$     read(temporary,*) USE_DIIS
!!$  end if

    !!__________________
    if (USE_DIIS) then
       ! Iterative use of Lanczos & DIIS
!!$     call getenv('Iterative',temporary)
!!$     if (temporary .eq. '') then
!!$        ITERATIVE = .false.
!!$     else
!!$        read(temporary,*) ITERATIVE
!!$     end if
       if (inflection==-1)  INFLECTION = MAXPAS

!!$     !!__________________
!!$     ! Number of Lanczos steps after an inflection in the eigenvalue
!!$     call getenv('Inflection', temporary)
!!$     if (temporary == '') then
!!$        INFLECTION = MAXPAS
!!$     else
!!$        read(temporary,*) INFLECTION
!!$     end if

       !!__________________
       ! Force threshold for call DIIS
!!$     call getenv('DIIS_Force_Threshold',temporary)
!!$     if (temporary .eq. '') then
!!$        DIIS_FORCE_THRESHOLD = 1.0d0
!!$     else
!!$        read(temporary,*) DIIS_FORCE_THRESHOLD
!!$     end if

       !!__________________
       ! Number of vectors kepts in memory for algorithm
!!$     call getenv('DIIS_Memory',temporary)
!!$     if (temporary .eq. '') then
!!$        DIIS_MEMORY = 12
!!$     else
!!$        read(temporary,*) DIIS_MEMORY
!!$     end if

       !!__________________
       ! prefactor multiplying forces
!!$     call getenv('DIIS_Step_size',temporary)
!!$     if (temporary .eq. '') then
        DIIS_STEP = DIIS_Step_size
!!$     else
!!$        read(temporary,*) DIIS_STEP
!!$     end if

       !!__________________
       ! times Increment_Size, max allowed diis step size

!!$     call getenv('FACTOR_DIIS',temporary)
!!$     if (temporary .eq. '') then
!!$        factor_diis = 5.0
!!$     else
!!$        read(temporary,*) factor_diis
!!$     end if
!!$
!!$     !!__________________
!!$     ! max diis iterations per call
!!$     call getenv('MAX_DIIS',temporary)
!!$     if (temporary .eq. '') then
        maxdiis = MAX_DIIS
!!$     else
!!$        read(temporary,*) maxdiis
!!$     end if

       !!__________________
       ! Check that the final state is indeed a saddle
!!$     call getenv('DIIS_Check_Eigenvector',temporary)
!!$     if (temporary .eq. '') then
        DIIS_CHECK_EIGENVEC = DIIS_Check_Eigenvector
!!$     else
!!$        read(temporary,*) DIIS_CHECK_EIGENVEC
!!$     end if

    else  ! Is this necessary ??
       DIIS_FORCE_THRESHOLD = 1.0d0
       DIIS_STEP            = 1.0d0
       DIIS_MEMORY          = 1
       DIIS_CHECK_EIGENVEC  = .false.
    end if


    !SECTION______________________ INPUT/OUTPUT

    !!__________________
    ! Set displacement threshold used to count
    ! how many atoms are involved in the event
!!$  call getenv('NPART_DR_THRESHOLD', temporary)
!!$  if (temporary .eq. '') then
!!$     NPART_DR_THRESHOLD = 0.1_dp
!!$  else
!!$     read(temporary,*) NPART_DR_THRESHOLD
!!$  end if


    ! File tracking  the file (event) number - facultative
!!$  call getenv('FILECOUNTER', temporary)
!!$  if (temporary .eq. '') then
     COUNTER  = FILECOUNTER
!!$  else
!!$     read(temporary,*) COUNTER
!!$  end if

    !!__________________
    ! General output for message
!!$  call getenv('LOGFILE', temporary)
!!$  if (temporary .eq. '') then
!!$     LOGFILE   = 'log.file'
!!$  else
!!$     read(temporary,*) LOGFILE
!!$  end if

    !!__________________
    ! list of events with success or failure
!!$  call getenv('EVENTSLIST', temporary)
!!$  if (temporary .eq. '') then
!!$     EVENTSLIST   = 'events.list'
!!$  else
!!$     read(temporary,*) EVENTSLIST
!!$  end if

    !!__________________
    ! Save the configuration at every step?
!!$  call getenv('Save_Conf_Int', temporary)
!!$  if (temporary .eq. '') then
!!$     SAVE_CONF_INT = .false.
!!$  else
!!$     read(temporary,*) SAVE_CONF_INT
!!$  end if

    !!__________________
    ! RESTART: It is useful only for ab-initio
!!$  call getenv('Write_restart_file', temporary)
!!$  if (temporary .eq. '') then
!!$     write_restart_file = .false.
!!$  else
!!$     read(temporary,*) write_restart_file
!!$  endif



    !!__________________
    ! If WRITE_REJECTED_EVENT is false, then we only write accepted ones
!!$  call getenv('WRITE_REJECTED_EVENT', temporary)
!!$  if (temporary .eq. '') then
!!$    WRITE_REJECTED_EVENT = .true.
!!$  else
!!$    read(temporary,*) WRITE_REJECTED_EVENT
!!$  endif

    !!__________________
    ! current data for restarting event
!!$  call getenv('RESTART_FILE', temporary)
!!$  if (temporary .eq. '') then
     RESTARTFILE = RESTART_FILE
!!$  else
!!$     read(temporary,*) RESTARTFILE
!!$  end if
!!$
!!$  !!__________________
!!$  ! Writes min. and sad. configurations in .xyz format.
!!$  call getenv('Write_xyz', temporary)
!!$  if (temporary .eq. '') then
!!$     write_xyz = .false.
!!$  else
!!$     read(temporary,*) write_xyz
!!$  endif

    !!__________________
    ! Reference configuration for refine saddle. Without ext.
!!$  call getenv('REFCONFIG', temporary)
!!$  if (temporary .eq. '') then
!!$     REFCONFIG   = 'refconfig.dat'
!!$  else
!!$     read(temporary,*) REFCONFIG
!!$  end if

    !!__________________HIDE: name for minima files
!!$  call getenv('FINAL', temporary)
!!$  if (temporary .eq. '') then
!!$     FINAL  = 'min'
!!$  else
!!$     read(temporary,*) FINAL
!!$  end if
    FINAL = adjustl(FINAL)

    !!__________________HIDE: name for saddle file
!!$  call getenv('SADDLE', temporary)
!!$  if (temporary .eq. '') then
!!$     SADDLE  = 'sad'
!!$  else
!!$     read(temporary,*) SADDLE
!!$  end if


    !!__________________HIDE: name for saddle file
!!$  call getenv('CHECK', temporary)
!!$  if (temporary .eq. '') then
!!$     CHECK  = 'con'
!!$  else
!!$     read(temporary,*) SADDLE
!!$  end if


    !_____________________________HIDE: Clean_wavefunct
    ! Hide Option, by default .false.
    ! we do not use the previously calculated wave function at some key
    ! points. This is for charged systems.
!!$  call getenv('Clean_wavefunct',temporary)
!!$  if (temporary == '') then
!!$     clean_wf = Clean_wavefunct
!!$  else
!!$     read(temporary,*) clean_wf
!!$     if (.not.( energy_type=="BSW" .or. energy_type=="OTF" .or. &
!!$                energy_type=="BAY" .or. energy_type=="BIG" )     &
!!$          .and. clean_wf ) then
!!$        write(unit6P,*) "Error : Clean_wavefunct option is only for bigdft"
!!$        stop
!!$     end if
!!$  end if
    !_____________________________HIDE: QM/MM
!!$  call getenv('MAXNEI', temporary)
!!$  if (temporary .eq. '') then
!!$     if (energy_type=='SWP') then
!!$        maxnei = min(natoms,100)
!!$     else
!!$        maxnei = 1
!!$     endif
!!$  else
!!$     read(temporary,*) maxnei
!!$  end if
!!$
!!$   ! number of quantum atoms
!!$  call getenv('NBR_QUNT', temporary)
!!$  if (temporary .eq. '' .and. energy_type .ne. "BSW" &
!!$       &  .and. energy_type .ne. "OTF" .and. energy_type .ne. "BAY") then
!!$     nbr_quantum = natoms
!!$  elseif ((energy_type == "BSW" .or. energy_type == "OTF" .or. energy_type == "BAY"&
!!$           &) .and. temporary .eq. "" ) then
!!$     write(unit6P,*) "Error : you have not given the number of quantum atoms"
!!$     write(unit6P,*) "The program will stop"
!!$     stop
!!$  elseif (temporary .ne. "") then
!!$     read(temporary,*) nbr_quantum
!!$  endif
!!$
!!$   ! number of quantum atoms to trash (buffer zone)
!!$  call getenv('NBR_QUNT_BUF', temporary)
!!$  if (temporary .ne. '' .and. energy_type .ne. "BSW" .and. energy_type .ne. "OTF" &
!!$            & .and. energy_type .ne. "BAY") then
!!$     write(unit6P,*) "Error: number of quantum atoms only usefull for BSW energy_calc "
!!$     write(unit6P,*) " we will stop  "
!!$     stop
!!$  elseif ( (energy_type == "BSW" .or. energy_type == "OTF" .or. energy_type == "BAY" &
!!$              &) .and. temporary .eq. "" ) then
!!$     write(unit6P,*) "Error : you have not given the number of quantum atoms"
!!$     write(unit6P,*) "The program will stop"
!!$     stop
!!$  elseif (temporary .eq. '') then
!!$     nbr_quantum_trash = 0
!!$  else
!!$     read(temporary,*) nbr_quantum_trash
!!$  endif
!!$
!!$  call getenv('PASSIVATE', temporary)
!!$  if ( (temporary .eq. ".true.") .and. energy_type .ne. "BSW" .and. energy_type .ne. "OTF" &
!!$              & .and. energy_type .ne. "BAY") then
!!$     write(unit6P,*) "Error: should only passivate if BSW energy_calc "
!!$     write(unit6P,*) " we will stop  "
!!$     write(unit6P,*) temporary
!!$     stop
!!$  elseif (temporary .eq. "" .or. temporary .ne. ".true." ) then
!!$     passivate = .false.
!!$  else
!!$     read(temporary,*) passivate
!!$  endif

    !_____________________________HIDE: local_coord
    if (TYPE_EVENTS == 'local_coord') then

       !     call getenv('Radius_Initial_Deformation', temporary)
       if (local_cutoff==-1.) then
          write(unit6P,*) 'Error: Radius_Initial_Deformation must be defined when TYPE_EVENTS is local'
          stop
!!$     else
!!$        read(temporary,*) LOCAL_CUTOFF
       end if

       !     call getenv('Coord_radius', temporary)
       !     if (temporary .eq. '') then
       if (coord_length==-1.)then
          write(unit6P,*) 'Error: Coord_length must be defined when TYPE_EVENTS is local_coord'
          stop
!!$     else
!!$        read(temporary,*) coord_length
       end if
       if (coord_number==-1)then
          write(unit6P,*) 'Error: Coord_number must be defined when TYPE_EVENTS is local_coord'
          stop
       end if

!!$     call getenv('Coord_number',temporary)
!!$     if (temporary .eq. '') then
!!$        write(unit6P,*) 'Error: Coord_number must be defined when TYPE_EVENTS is local_coord'
!!$        stop
!!$     else
!!$        read(temporary,*) coord_number
!!$     end if

       ! type of atom
!!$     call getenv('Type_selected',temporary)
!!$     if (temporary .eq. '') then
        type_sel = Type_selected
!!$     else
!!$        read(temporary,*) type_sel
!!$     end if

    end if
    !_____________________________HIDE: Dual_system
    dual_search=Dual_system
!!$  call getenv('Dual_system',temporary)
!!$  if (temporary == '') then
!!$     dual_search = .false.
!!$  else
    !     read(temporary,*) dual_search
    if ( .not. (TYPE_EVENTS=='local' .or. TYPE_EVENTS=='list_local') &
         .and.  dual_search ) then
       write(unit6P,*) 'Error: Dual_system is defined only for local or list_local'
       write(unit6P,*)  dual_search, TYPE_EVENTS
       stop
    end if
!!$  end if

    if ( dual_search ) then
       !     if (size_system==-1.)then
       !     call getenv('Size_system', temporary)
       !    if (temporary .eq. '') then
       !        write(unit6P,*) 'Error: size_system must be defined if dual_search'
       !        stop
       !     else
       !     read(temporary,*) size_system
       if ( size_system <= LOCAL_CUTOFF ) then
          write(unit6P,*) 'Error: Radius_Initial_Deformation > Size_system'
          stop
       end if
       !     end if
    end if
    !_____________________________

    ! We set up other related parameters
    vecsize = 3*natoms

    ! And we allocate the vectors
    !FUCKING ALLOCATION!
    allocate(typat(natoms))
    allocate(constr(natoms))
    allocate(force(vecsize))
    allocate(pos(vecsize))
    allocate(posref(vecsize))
    allocate(initial_direction(vecsize))
    allocate(Atom(vecsize))
    allocate(old_projection(VECSIZE))
    allocate(projection(VECSIZE))

    if (type_events=='energy')then
       allocate(eatom(natoms))
    end if
    
    force(:) = 0.0d0
    pos(:) = 0.0d0
    posref(:) = 0.0d0
    initial_direction(:) = 0.0d0
    old_projection(:) = 0.0d0
    projection(:) = 0.0d0

    x => pos(1:NATOMS)
    y => pos(NATOMS+1:2*NATOMS)
    z => pos(2*NATOMS+1:3*NATOMS)

    xref => posref(1:NATOMS)
    yref => posref(NATOMS+1:2*NATOMS)
    zref => posref(2*NATOMS+1:3*NATOMS)

    fx => force(1:NATOMS)
    fy => force(NATOMS+1:2*NATOMS)
    fz => force(2*NATOMS+1:3*NATOMS)

    ! If LANCZOS_MIN, we check how it changes the energy of the system by applying the projection.
    ! This is done at the end of LANCZOS procedure, but this is done only for the minima.
    ! (see subroutine check_min )
    IN_MINIMUN = .False.

!!$  if (LOCAL_FORCE) then
!!$     allocate(inner_list(natoms))
!!$     allocate(outer_list(natoms))
!!$     allocate(mask_inner_region(vecsize))
!!$     r2_inner = INNER_REGION*INNER_REGION
!!$     r2_outer = (INNER_REGION+OUTER_REGION)*(INNER_REGION+OUTER_REGION)
!!$  endif

  END SUBROUTINE read_parameters


  subroutine write_parameters( )

    use defs
    use lanczos_defs
    use saddles
    use random
    implicit none

    !Local variables
    integer :: ierror, i
    integer, dimension(8) :: values
    character(len=20) :: fname  ! same len than LOGFILE in defs.f90
    character(len=4)  :: ext
    character(len=20)  :: fmt_underheader,fmt_header,fmt_int,fmt_real,fmt_log,fmt_str,fmt_ints

    logical :: exists_already
    !_______________________
    ! Initialization of seed in random.
    ! WARNING: There will be an idum per
    ! iproc. Suggested by Laurent Karim Beland,
    ! UdeM 2011
    call date_and_time(values=values)

    if ( .not. setup_initial ) then
       idum = -1 * mod( (1000 * values(7) + values(8))+iproc, 1024)
    else
       idum = 0
    end if

    if ( iproc == 0 ) then

       ! We check what is the last 'LOGFILE.number' in the directory. LOGFILE will be LOGFILE.number+1

       i = 1
       do

          write(ext,'(i4)')  i
          fname = trim(LOGFILE)// "." //adjustl(ext)
          fname = adjustl(fname)
          inquire( file = fname, exist = exists_already )

          if ( .not. exists_already ) then
             logfile = fname
             exit
          end if
          i = i + 1

       end do

       ! We write down the various parameters for the simulation

       ! Declare specific format for KFLOG writes.
       fmt_header = '(T55,A)'
       fmt_underheader = '(T50,A)'
       fmt_str    = '(A,T75,A30)'
       fmt_int    = '(A,T75,I30)'
       fmt_log    = '(A,T75,L30)'
       fmt_real   = '(A,T75,F30.4)'
       fmt_ints   = '(A,T75,3I10)'

       open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='rewind',iostat=ierror)
       write(flog,*) '****************************** '
       write(flog,*) 'WELCOME TO  ART nouveau (with NDM :-)!) '
       write(flog,*) '****************************** '
       call timestamp ('Start')
       write(flog,'(1X,A39,f12.3)')  ' - Version number of  ART            : ', VERSION_NUMBER
       write(flog,*) ''
       write(flog,'(1X,A39,A16  )')  ' - Event type                        : ', trim(eventtype)
       write(flog,'(1X,A39,f12.4)')  ' - Temperature                       : ', temperature
       write(flog,'(1X,A39,I12  )')  ' - Number of atoms                   : ', natoms
       write(flog,'(1X,A39,I12  )')  ' - Number of events                  : ', number_events
       write(flog,'(1X,A39,I12  )')  ' - Maximum number of neighbours      : ', maxnei

       write(flog,'(1X,A39  )')  ' - Atomic types                      : '
       do i =1, 5
          if (type_name(i) .ne. '') then
             write(flog,'(1X,A31,I3,A5,a12  )')  ' Type ',i,': ', trim(type_name(i))
          end if
       end do

       write(flog,*) ' '
       write(flog,*) 'Selection of the event '
       write(flog,*) '********************** '
       write(flog,'(1X,A39,A15)')   ' - Type of events                    : ', trim(TYPE_EVENTS)
       if (TYPE_EVENTS .eq. 'local') then
          write(flog,'(1X,A39,f12.4)')   ' - Radius of deformation (local ev.) : ', LOCAL_CUTOFF
          if (preferred_atom .gt. 0) then
             write(flog,'(1X,A39,i12)')   ' - Central atom for events           : ', preferred_atom
          else
             write(flog,'(1X,A39,A12)')   ' - Central atom for events           :   none (all equal) '
          end if
       end if
       write(FLOG,fmt_log) ' Check if events are reversible (CHECK_CONNECTIVITY)  : ', CHECK_CONNECTIVITY


       write(flog,*) ' '
       write(flog,*) 'local force parameters '
       write(flog,*) '*********************** '
!!$  write(FLOG,fmt_log) ' Use of local forces (LOCAL_FORCES)  : ', LOCAL_FORCE
       write(FLOG,fmt_log) ' Use of local Lanczos (LOCAL_LANCZOS) : ', LOCAL_LANCZOS
       write(FLOG,fmt_log) ' Use of global convergence (GLOBAL_CONVERGENCE) : ', GLOBAL_CONVERGENCE
       if (LOCAL_LANCZOS) then
          write(FLOG,fmt_real)  ' Radius of the inner region: ', INNER_REGION
          write(FLOG,fmt_real)  ' Width of the outer region: ', OUTER_REGION
          write(FLOG,fmt_log) ' Use of global convergence (GLOBAL_CONVERGENCE) : ', GLOBAL_CONVERGENCE
       endif



       write(flog,*) ' '
       write(flog,*) 'Minimization parameters '
       write(flog,*) '*********************** '
       write(flog,'(1X,A51,F8.4)')  ' - DT_MAX_FIRE                                  : ', DT_MAX_FIRE
       write(flog,'(1X,A51,I8)')    ' - MAX_RELAXATION_STEP_FIRE                     : ', MAX_ITER_FIRE
       write(flog,'(1X,A51,I8)')    ' - MAX_REL_STEPS_PERP_FIRE                      : ', MAX_PERP_ITER_FIRE
       write(flog,'(1X,A51,I8)')    ' - MAX_REL_STEPS_INCREMENT_FIRE                 : ', MAX_PERP_INCR_FIRE

       write(flog,*) ' '
       write(flog,*) 'Activation parameters '
       write(flog,*) '********************* '
       write(flog,'(1X,A51,F8.4)')  ' - Eigenvalue threshold                          : ', EIGEN_THRESH
       write(flog,'(1X,A51,F8.4)')  ' - Total force threshold (saddle)                : ', EXITTHRESH
       write(flog,'(1X,A51,F8.4)')  ' - Initial step size                             : ', INITSTEPSIZE
       write(flog,'(1X,A51,F8.4)')  ' - Increment size                                : ', INCREMENT
       write(flog,'(1X,A51,F8.4)')  ' - Atomic displacement for breaking the symmetry : ', sym_break_dist
       write(flog,'(1X,A51,2I8)')    ' - Number of vectors computed by Lanzcos         : ', NVECTOR_LANCZOS_H,&
            NVECTOR_LANCZOS_C
       write(flog,'(1X,A51,F8.4)')  ' - Delta displacement for derivative in Lanzcos  : ', DEL_LANCZOS
       write(flog,*) ' '
       write(flog,'(1X,A51,L8)')    ' - Calculation of the Hessian for each minimum   : ', LANCZOS_MIN
       write(flog,'(1X,A51,I8)')    ' - Min. number of ksteps before calling lanczos  : ', KTER_MIN
       write(flog,'(1X,A51,F8.4)')  ' - Factor mult. INCREMENT for leaving basin      : ', BASIN_FACTOR
       write(flog,'(1X,A51,I8)')    ' - Maximum number of iteractions (basin -kter)   : ', MAXKTER
       write(flog,'(1X,A51,I8)')    ' - Maximum number of perpendicular moves (basin) : ', MAXKPERP
       write(flog,'(1X,A81,I8)')    ' - Number of steps for smoothing exit from basin  (Smooth_Dir_Change): ', SMOOTH_DIR_CHANGE

       write(flog,*) ' '
       write(flog,'(1X,A51,I8)')    ' - Maximum number of iteractions (activation)    : ', MAXPAS
       write(flog,'(1X,A51,I8)')    ' - Maximum number of perpendicular moves (activ) : ', MAXIPERP
       write(flog,'(1X,A51,F8.4)')  ' - Force threshold for perpendicular relaxation  : ', FTHRESHOLD
       write(flog,'(1X,A51,F8.4)')  ' - Fraction of displacement over the saddle      : ', PUSH_OVER

       write(flog,*) ' '
       write(flog,'(1X,A39,L12)')    ' Use DIIS for convergence to saddle     : ', USE_DIIS
       if (USE_DIIS) then
          write(flog,'(1X,A40,F12.4)') ' - Total force threshold to call DIIS : ', DIIS_FORCE_THRESHOLD
          write(flog,'(1X,A40,F12.4)') ' - Step size to update positions      : ', DIIS_STEP
          write(flog,'(1X,A40,I12)')   ' - Memory (number of steps)           : ', DIIS_MEMORY
          write(flog,'(1X,A40,L12)')   ' - Check eigenvector at saddle        : ', DIIS_CHECK_EIGENVEC
       end if

       write(flog,*) ' '
       write(flog,*) 'Input / Output '
       write(flog,*) '********************* '
       write(flog,'(1X,A50,A15)')  ' - Name of log file                  : ', trim(LOGFILE)
       write(flog,'(1X,A50,A15)')  ' - Liste of events                   : ', trim(eventslist)
       write(flog,'(1X,A50,A15)')  ' - Reference configuration           : ', trim(refconfig)
       write(flog,'(1X,A50,A15)')  ' - Restart file                      : ', trim(restartfile)
       write(flog,'(1X,A50,A15)')  ' - Prefix for minima (file)          : ', trim(FINAL)
       write(flog,'(1X,A50,A15)')  ' - Prefix for saddle points          : ', trim(SADDLE)
       write(flog,'(1X,A50,L15)')  ' - Write rejected step (WRITE_REJECTED_EVENT) :', WRITE_REJECTED_EVENT 
       write(flog,*) '********************* '
       write(flog,*) ' '
       write(flog,'(1X,A34,I17)') ' - The seed is                  : ', idum
       write(flog,'(1X,A34,L17)') ' - Restart                      : ', restart
       close(flog)

    end if

  END SUBROUTINE write_parameters

  !> ART timestamp
  !! it stamps the time at the beggining and at the final of the simulation
  subroutine timestamp (str)

    use defs, only : FLOG
    implicit none
    !Arguments
    character(len=*), intent(in) :: str
    !Local variables
    integer         :: sec, min, hour, day, month, year
    character(len=1), parameter :: dash = "-"
    character(len=1), parameter :: colon = ":"
    character(len=3), parameter :: prefix = ">> "
    character(len=3)  :: month_str(12) =     &
         (/'JAN','FEB','MAR','APR','MAY','JUN', &
         'JUL','AUG','SEP','OCT','NOV','DEC'/)

    integer :: values(8)

    call date_and_time(values=values)
    year  = values(1)
    month = values(2)
    day   = values(3)
    hour  = values(5)
    min   = values(6)
    sec   = values(7)

    write(FLOG,1000)  prefix, trim(str), colon, day, dash, month_str(month),&
         &                 dash, year, hour, colon, min, colon, sec

1000 format(2a,a1,2x,i2,a1,a3,a1,i4,2x,i2,a1,i2.2,a1,i2.2)

  END SUBROUTINE timestamp


end module read_parameters_mod

!>  ART git module
!!  Write Git commit and other compilation info
module git

  implicit none

contains

  subroutine write_git_commit(fileunit)

    use defs
    use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
         stdout=>output_unit, &
         stderr=>error_unit

    integer, intent(in), optional :: fileunit
    integer :: fileunit_


    if (present(fileunit)) then
       fileunit_ = fileunit
    else
       fileunit_ = stdout
    end if

#ifdef GIT_COMMIT
    if (iproc == 0) then
       write(fileunit_,*)
       write(fileunit_,*) 'COMMIT ID: '//GIT_COMMIT
       write(fileunit_,*) 'DATE     : '//GIT_DATE
       if (GIT_CLEAN.eq.'1') then
          write(fileunit_,*) 'Local files match commit id.'
       else
          write(fileunit_,*) 'Local files do not match commit id.'
       endif
       write(fileunit_,*) '****************************** '
    end if
#endif
    if (iproc == 0) then
       write(fileunit_,*) 'VERSION COMPILED ON :          '
       write(fileunit_,*) __DATE__//' '//__TIME__
       write(fileunit_,*) '****************************** '
       write(fileunit_,*)
    end if

  end subroutine write_git_commit

end module git



!>  ART write_parameters
!! * Define and open the log.file
!! * Write in it the parameters defining the simulation
