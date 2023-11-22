!> @file
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!! Modified by:
!! -EM 2010, see ~/AUTHORS
!! -Laurent Karim Beland, UdeM, 2011. For working with QM/MM !!
!! -EM 2011, see ~/AUTHORS
!!
!> ART initialize
!!   Initialization of art method
!!   It relaxes it into a local minimum without volume optimization
!!   The initial configuration is a "reference configuration", it will be the reference
!!   configuration until a new event is accepted.
!!

! $Revision$
! $Date: 2015-12-14 11:01:46 -0500 (Lun, 14 déc 2015) $
! $Id: initialize.f90 1534 2015-12-14 15:59:38Z mickael $


module initialize_mod
  use calcforce_mod,only:calcforce
  use end_art_mod,only:end_art
  use min_converge_mod,only:min_converge
  use convert_to_chain_mod,only:convert_to_chain
  use write_refconfig_mod,only:write_refconfig
  use storage,only:store
  use min_converge_mod,only:check_min
  use generate_local_region_mod, only: initial_local_region
  use montecarlo_mod,only:lparapath,nparapath
  use ndm2art2ndm,only:parapath,atcfart,celart,boxart
!  use posana,only:anaposart
contains
subroutine initialize()

  use defs
  use lanczos_defs,  only : LANCZOS_MIN
  use saddles,       only : g_pos, GUESSFILE
  implicit none

  !Local variables
  integer :: i, ierror
  character(len=20) :: dummy, fname
  character(len=4)  :: scounter
  logical           :: flag, success

  integer           :: nat_test
  integer, allocatable               :: typ_a(:)    ! Atomic type
  integer, allocatable               :: const_a(:)  ! Constraints
  real(kind=8), allocatable, target  :: pos_a(:)    ! Working positions of the atoms
  real(kind=8), dimension(3)          :: boxref_     ! Reference box from posinp file

  character(len=1)                    :: boundary_b
  integer, allocatable                :: typ_b(:)    ! Atomic type
  integer, allocatable                :: const_b(:)  ! Constraints
  real(kind=8), allocatable, target   :: pos_b(:)    ! Working positions of the atoms
  real(kind=8), dimension(:), pointer :: xb, yb, zb
  real(kind=8), dimension(3) :: boxl

  real(kind=8)                        :: determinant
  !_______________________
  ! assigned value of VECSIZE
  VECSIZE = 3*NATOMS

  ! Set local force to false
  use_local_forces = .false.

  ! Read the counter in order to continue the run where it stopped or for refine
  ! Format:
  ! Counter:     1000

!!$  if (.not. restart) then


     inquire( file = COUNTER, exist = flag )
     if ( flag .and. iproc == 0 ) then
        open(unit=FCOUNTER,file=COUNTER,status='old',action='read',iostat=ierror)
        read(FCOUNTER,'(A12,I6)') dummy, mincounter
        close(FCOUNTER)
     else
        mincounter = 1000
     end if
#ifdef PARA
     if (lparapath) then
        mincounter=mincounter+parapath%image*number_events
     end if
#endif
     
!!$  end if
     refcounter = mincounter
  ! we read the initial/reference configuration


!!$  if ( new_event .or. restart ) then   ! new_event=true with eventtype=new
!!$     ! Read initial atomic file
!!$     ! is it neccesary for restart ?
!!$
!!$     nat_test = NATOMS
!!$     allocate(typ_a(NATOMS))
!!$     allocate(pos_a(3*NATOMS))
!!$     allocate(const_a(NATOMS))
!!$     call init_all_atoms( nat_test, typ_a, pos_a, const_a, boxref_, boundary, nproc, iproc, trim(REFCONFIG) )
!!$  else if ( eventtype == 'REFINE_AND_RELAX' .or. eventtype == 'REFINE_SADDLE' ) then
!!$     ! Read reference atomic file
!!$     write(unit6P,*)'notcoded initialize 1  eventtype <> new'
!!$     stop
!!$     nat_test = NATOMS
!!$     allocate(typ_a(NATOMS))
!!$     allocate(pos_a(3*NATOMS))
!!$     allocate(const_a(NATOMS))
!!$     call init_all_atoms( nat_test, typ_a, pos_a, const_a, boxref_, boundary, nproc, iproc, trim(REFCONFIG) )
!!$  end if

  !assign the data from the atomic file
!!$  if ( .not. restart ) then
!!$     typat(:)   = typ_a(:)
!!$     pos(:)     = pos_a(:)
!!$     constr(:)  = const_a(:)
!!$     boxref(:)  = boxref_(:)

!!$     box = boxref
!!$
!!$     deallocate(typ_a)
!!$     deallocate(pos_a)
!!$     deallocate(const_a)
!!$
!!$     ! We read the position for the presumed saddle point.
!!$     if ( eventtype == "GUESS_DIRECTION" ) then
!!$        write(unit6P,*)'notcoded initialize 2  eventtype <> new'
!!$        stop
!!$
!!$        nat_test = NATOMS
!!$        allocate(typ_b(NATOMS))
!!$        allocate(pos_b(3*NATOMS))
!!$        allocate(const_b(NATOMS))
!!$        call init_all_atoms( nat_test, typ_b, pos_b, const_b, boxref_, boundary_b, nproc, iproc, trim(GUESSFILE) )
!!$
!!$        ! Let's check if it is in the same conditions as posinp.
!!$        write(unit6P,*) "eventtype: ", eventtype
!!$        if ( boundary /= boundary_b ) then
!!$           if ( iproc == 0 ) write(unit6P,*) "GUESS: Different type of boundary conditions"
!!$           call end_art()
!!$        end if
!!$        do i = 1, NATOMS
!!$           if ( typ_b(i) /= typat(i) ) then
!!$              if ( iproc == 0 ) write(unit6P,*) "GUESS: Different type of atoms"
!!$              call end_art()
!!$           end if
!!$           if ( const_b(i) /= constr(i) ) then
!!$              if ( iproc == 0 ) write(unit6P,*) "GUESS: Different type of constraints"
!!$              call end_art()
!!$           end if
!!$        end do
!!$        ! g_pos in module saddles
!!$        allocate(g_pos(vecsize))
!!$        g_pos(:) = 0.0d0
!!$        g_pos(:) = pos_b(:)
!!$
!!$        deallocate(typ_b)
!!$        deallocate(pos_b)
!!$        deallocate(const_b)
!!$     end if
!!$  endif

  ! write
  if ( iproc == 0 ) then
     open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
     write(FLOG,'(1X,A34,I17)') ' - Mincounter                   : ', mincounter
     close(FLOG)
  end if

  ! We rescale the coordinates. For what ??
  scalaref = 1.0d0
  scala = scalaref
!  write(unit6P,*) 'before initialize_potential'
  !call initialize_potential()         ! Initialize Potential (CORE)
  call calcforce( NATOMS, pos, boxref, force, total_energy, evalf_number )
!  write(unit6P,*) 'after initialize_potential this is ok for this test'

  ! for output files
  if ( iproc == 0 ) call convert_to_chain( refcounter, 4, scounter )
  fname = FINAL // scounter
  conf_initial = fname
  ! If this is a new event we relax
  If_ne: if ( new_event .and. (.not. restart) ) then  ! cas standard
     posref = pos                     ! New reference configuration.
     call min_converge( success )     ! Converge the configuration to a local minimum
     write(unit6P,*)'initial minimization success',success
!     call anaposart(atcfart,celart,boxart)
     posref = pos                     ! New reference configuration.
     ref_energy = total_energy

     if ( iproc == 0 ) then
        ! THIS iS NOT NECESSARY IN BIGDFT
        call write_refconfig( )       ! Write reference in REFCONFIG.
        call store( fname )           ! Store the configuration into fname. fname=minxxx

        open( unit = FLOG, file = LOGFILE, status = 'unknown',&
             & action = 'write', position = 'append', iostat = ierror )
        write(unit6P,*) 'BART: Configuration stored in file ',fname
        write(FLOG,'(1X,A34,A17)') ' - Configuration stored in file : ', trim(fname)
        if ( .not. success ) then
           write(FLOG,'(1X,A)') "ERROR: Initial configurations is not a minimum"
           !write(unit6P,*) 'Just before end_art() (initialize line 171)'
           call end_art()
        end if
        close(FLOG)
     end if
     mincounter = mincounter + 1
     ! if dual_search we dont do this check at the
     ! beginning. It is not well defined.
     if ( success .and. .not. dual_search ) then
        if ( LANCZOS_MIN .or. setup_initial ) call check_min( 'M' )
     end if

  else if ( (.not. new_event) .and. (.not. restart) ) then
     write(unit6P,*)'pas_OKARTINIT2'
     stop
!!$     ! once we have the total energy we copy as reference values
!!$     posref = pos
!!$     ref_energy = total_energy
!!$
!!$     if ( iproc == 0 ) then
!!$        call store( fname )           ! Store the reference configuration into fname.
!!$
!!$        open( unit = FLOG, file = LOGFILE, status = 'unknown',&
!!$             & action = 'write', position = 'append', iostat = ierror )
!!$        write(unit6P,*) 'BART: Ref. Configuration stored in file ',fname
!!$        write(FLOG,'(1X,A34,A17)') ' - Configuration stored in file : ', trim(fname)
!!$        close(FLOG)
!!$     end if
!!$     mincounter = mincounter + 1
!!$
!!$     ! now we read the positions of the configuration we want to refine it
!!$     nat_test = NATOMS
!!$     allocate(typ_b(NATOMS))
!!$     allocate(pos_b(3*NATOMS))
!!$     allocate(const_b(NATOMS))
!!$     call init_all_atoms( nat_test, typ_b, pos_b, const_b, boxref_, boundary_b, nproc, iproc, 'posinp' )
!!$
!!$     ! Let's check if it is in the same conditions as the reference.
!!$     if ( boundary /= boundary_b ) then
!!$        if ( iproc == 0 ) write(unit6P,*) "posinp: Different boundary condition"
!!$        call end_art()
!!$     end if
!!$
!!$     do i = 1, NATOMS
!!$        if ( typ_b(i) /= typat(i) ) then
!!$           if ( iproc == 0 ) write(unit6P,*) "posinp: Different type of atoms"
!!$           call end_art()
!!$        end if
!!$
!!$        if ( const_b(i) /= constr(i) ) then
!!$           if ( iproc == 0 ) write(unit6P,*) "posinp: Different type of constraints"
!!$           call end_art()
!!$        end if
!!$
!!$     end do
!!$     ! our real starting point
!!$     pos(:)  = pos_b(:)
!!$     box(:)  = boxref_(:)
!!$
!!$     if ( boundary == 'T' ) then
!!$        write(unit6P,'(3f14.6)') cell
!!$     else
!!$        write(unit6P,*) 'box: ', box
!!$     end if
!!$
!!$     write(unit6P,*) 'positions: ', pos(1), pos(1+NATOMS), pos(1+2*NATOMS)
!!$
!!$     deallocate(typ_b)
!!$     deallocate(pos_b)
!!$     ! we need the total energy for this configuration
!!$     boxl = box * scala               ! we compute at constant volume.
!!$
!!$     call calcforce( NATOMS, pos, boxl, force, total_energy, evalf_number )

  end if If_ne

  if (LOCAL_FORCE)  call initial_local_region()
END SUBROUTINE initialize
end module initialize_mod
