!> @file
!! @author
!!    Copyright (C) Normand Mousseau, June 2001
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> ART module run_art
!! 5 routines to:
!! init_mpi:          initialize mpi processes
!! init_conf:         initialize the system
!! art_search:        find a neighboring minimum
!! report_and_check : report info about the search and
!!                    check that final state is a minimum


module run_art
  use Tpara,only:mpi_world
  use defs
  use arret_ndm_mod,only:arret_ndm
  use montecarlo_mod,only:lparapath,nparapath
  use git
  use ndm2art2ndm,only:ndm2art
  use random
  use end_art_mod,only:end_art
  use find_saddle_mod,only:find_saddle
  use initialize_mod,only:initialize
  use min_converge_mod,only:min_converge
  use read_parameters_mod,only:read_parameters,write_parameters
  use save_restart,only:save_state,restart_states
  use print_event_mod
  use disp,only:displacement
  USE atomconfig,only : atom_config
  USE cellconfig, only:cell_config
  USE boxconfig,only:box_config
  use Tpara,only:para_space_config
  use storage,only:store
  use convert_to_chain_mod,only:convert_to_chain
  use  min_converge_mod,only:check_min
  implicit none

  character(4)  :: scounter
  real(kind=8)  :: saddle_energy
  real(kind=8)  :: delta_e
  real(kind=8)  :: a1, b1, c1, prod
  logical       :: success
  character(8)  :: accept


contains




  subroutine init_conf

!    character(len=*), intent(in), optional :: which
    
    integer                                :: ierror, ierr

    ! Read the various parameters and options defining the run
    call read_parameters( )
#ifdef PARA
    if (lparapath)then 
       if (mod(number_events,nparapath).ne.0) then
          write(unit6P,*)'nparapath/number_events <>0 STOP',nparapath,number_events
          call arret_ndm
       else
          number_events=number_events/nparapath
       end if
    end if
#endif
!    write(unit6P,*)'JP out readp'
    call  ndm2art
!    write(unit6P,*)'JP out ndm2art'
    ! If restartfile exists, then we restart from where we left.
!    inquire ( file = restartfile, exist = restart )
!    if ( restart ) &
!         & call restart_states( state_restart, ievent_restart, iter_restart, atp )


    
!!$    
!!$#ifdef PARA
!!$    call mpi_world%barrier
!!$    !call MPI_Barrier( MPI_COMM_WORLD, ierr )
!!$    call write_parameters( )            ! Write options in LOGFILE.
!!$    call mpi_world%barrier
!!$!call MPI_Barrier( MPI_COMM_WORLD, ierr )
!!$#else
    
    call write_parameters( )            ! Write options in LOGFILE.
!    write(unit6P,*)'JP out writeparam'
!!$#endif

    ! Open the log file and get ready for the simulation
    if ( iproc == 0 ) then              ! Report
       open( unit = FLOG, file = LOGFILE, status = 'unknown',&
            & action = 'write', position = 'append', iostat = ierror )
       if ( restart ) then
          write(FLOG,'(1X,A)') ' Call restart'
       else
          write(FLOG,'(1X,A)') ' Start with new event    '
       end if
       close(FLOG)
    end if


    
!   call initialize_potential()
!ifdef LAMMPS_VERSION
!endif

!!$    if (present(which)) then
!!$       if (which == "potential") then
!!$
!!$          scalaref = 1.0d0
!!$          scala    = scalaref
!!$
!!$          if (iproc == 0) write(unit6P,*) "initialized up to potential setup"
!!$          return
!!$       end if
!!$    else

!!$       write(unit6P,*) "(1) I am the proc :", iproc
       
       call initialize( )         ! Initialize positions and potential

!    write(unit6P,*)'JP out initialiaze'
!!$       write(unit6P,*) "(2) I am the proc :", iproc
       ! If the restart file exists, then we make sure that we do not overwrite the files
!!$       if ( restart ) then
!!$          if ( iproc == 0 ) call convert_to_chain( refcounter, 4, scounter )
!!$          ! Information for FLIST
!!$          conf_initial = trim( FINAL // scounter )
!!$       else
          ! Set up of mincounter & reference configuration,
          ! if new_event then relaxes it into a local minimum.
          evalf_number = 0
          ievent_restart = 1
!!$       end if
!!$       write(unit6P,*) "(3) I am the proc :", iproc

!!$    end if

    return

  end subroutine init_conf


  subroutine art_search(fname)

    use lanczos_defs, only: projection, eigenvalue
    use saddle_push

    implicit none

    character(20), intent(out)              :: fname
    integer                                 :: ierror
    integer                                 :: npart ! Number of atoms participating to the event
    real(kind=8)                            :: delr
    real(dp), dimension(:), allocatable     :: tmp_pos
    real(dp), dimension(:), allocatable     :: pos_saddle,force_saddle, projection_saddle,pos_min
    real(dp)                                :: energy_saddle,eigenvalue_saddle,energy_min
  real(kind=8), dimension(:), allocatable :: del_pos
  real(kind=8) :: difpos
  real(kind=8)  :: a1, b1, c1, prod



    ! If it is a restart event for the activation,
    ! phase 1, 2 or 4, or it is not a restart event
    ! then we call find_saddle.
    if ( .not. ( restart .and. ( state_restart == 3 ) ) ) then
       ! atp is the number of attempts for event.
       if (.not. restart) atp = 1
       do
          if ( iproc == 0 ) call print_event( ievent, temperature )
          call find_saddle( success, saddle_energy )
          if ( success ) then
             exit
          else
             atp = atp + 1
          end if
       end do
    end if

    ! If not a new event, then a convergence to
    ! the saddle point, We do not go further.
    if ( .not. NEW_EVENT .and. ( eventtype == 'REFINE_SADDLE' ) ) call end_art ()
    ! If it is a restart event of type 3, we
    ! are starting at the right point.
    if ( restart .and. ( state_restart == 3 ) ) then
       if ( iproc == 0 ) then        ! Report
          call print_event( ievent, temperature )
          open( unit = FLOG, file = LOGFILE, status = 'unknown',&
               & action = 'write', position = 'append', iostat = ierror )
          write(FLOG,'(1X,A)') ' - Restart event: We restart at the saddle point'
          close(FLOG)
          call convert_to_chain( mincounter, 4, scounter )
       end if
       saddle_energy = total_energy
       conf_saddle = trim( SADDLE // scounter )
       restart = .false.
    else if ( write_restart_file ) then
       ! We save the current state for a possible
       ! restart.
       state_restart = 3
       iter_restart  = 0
       if ( iproc == 0 ) call save_state( state_restart, iter_restart, projection )
    end if


    ! Saddle has been reached, try pushingkthe system towards a new basin.
    ! The displacement is made along the direction of negative curvature,
    ! away from the initial minimum.
    ! Notify user about the outcome of the push.
    ! If pushing could not be made without increasing energy, try finding an
    ! other saddle.

    write(unit6P,*)'JP501',success
    if (CHECK_CONNECTIVITY) then
       !CRC
       STOP
!!$       allocate(pos_saddle(VECSIZE))
!!$         allocate(force_saddle(VECSIZE))
!!$         allocate(projection_saddle(VECSIZE))
!!$         pos_saddle = pos
!!$         energy_saddle = saddle_energy
!!$         force_saddle = force
!!$         projection_saddle = projection
!!$         eigenvalue_saddle = eigenvalue
    endif 

    !CRC version Mousseau emplacé par version MySrC
    
!!$    allocate(tmp_pos(VECSIZE))
!!$    call push_at_saddle(tmp_pos, pos, posref, force, &
!!$                        projection, eigenvalue, saddle_energy)
!!$    pos = tmp_pos
!!$    deallocate(tmp_pos)
!!$    ierror = saddle_push_log()
!!$
!!$    ! If the push went well, converge to the new minimum.
!!$    ! Write the configuration in a min.... file
!!$    if (push_error == -1) then
!!$       success = .false.
!!$       delta_e = total_energy - ref_energy
!!$       fname = FINAL // '_push_failed'
!!$       conf_final = fname
!!$    else
!!$       call min_converge( success )     ! And we converge to the new minimum.
!!$       delta_e = total_energy - ref_energy
!!$       if ( iproc == 0 ) then
!!$          ! We write the configuration in a min.... file.
!!$          call convert_to_chain( mincounter, 4, scounter )
!!$          write(unit6P,*) 'BART: Mincounter is ', mincounter,', scounter is ', scounter
!!$          fname = FINAL // scounter
!!$          conf_final = fname
!!$       end if
!!$    end if

   
!    The next lines are for the displacement from the saddle away from the initial minimum 
!    in order to place the configuration in a new basin. The displacement is made along the 
!    direction of negative curvature away from the initial minimum.
     
     allocate(del_pos(VECSIZE))       ! We compute the displacement.
     call boundary_cond( del_pos, pos, posref )   ! Applies the boundary conditions
     difpos  = sqrt( dot_product(del_pos,del_pos) )
     del_pos = del_pos/difpos 


     a1 = dot_product(del_pos,projection) 
     b1 = dot_product(force,projection)
     c1 = dot_product(del_pos,force)
     deallocate(del_pos)

     if ( abs(a1) < 0.1d0 ) then      ! pushing in the direction of projection (assuming sign ok) 
        prod = 1.0d0
        if ( iproc == 0 ) then 
         write(unit6P,*) 'BART :WARNING'
         write(unit6P,*) 'BART :Projection and displacement vectors almost perpendicular'
         write(unit6P,*) 'BART :to each other. Assuming projection points in right direction'
        end if
     else                             ! just keep the sign of the dot product
        if ( a1 > 0.0 ) then
           prod =  1.0d0
        else
           prod = -1.0d0
        end if 
     end if 

     ! We finally push over the saddle point
     pos = pos + prod * PUSH_OVER * difpos * projection
!     write(unit6P,*)'JP502',iproc
     call min_converge( success )     ! And we converge to the new minimum.
     delta_e = total_energy - ref_energy 
!     if ( iproc == 0 ) then
                                      ! We write the configuration in a min.... file.
        call convert_to_chain( mincounter, 4, scounter )
        write(unit6P,*) 'BART: Mincounter is ', mincounter,', scounter is ', scounter
        fname = FINAL // scounter
        conf_final = fname
        call store( fname )  !fname =minxxx
!     end if
                                      ! Magnitude of the displacement (utils.f90).


!    write(unit6P,*)'JP503',fname
!    call store( fname )

    if (CHECK_CONNECTIVITY) then
       !CRC
       STOP
!!$       allocate(pos_min(VECSIZE))
!!$        pos_min = pos
!!$        energy_min = total_energy
!!$
!!$        pos = pos_saddle 
!!$        saddle_energy =  energy_saddle
!!$        force =  force_saddle 
!!$        projection =  projection_saddle
!!$        eigenvalue =  eigenvalue_saddle 
!!$
!!$        allocate(tmp_pos(VECSIZE))
!!$        call push_at_saddle(tmp_pos, pos, posref, force, &
!!$                        projection, eigenvalue, saddle_energy,-1)
!!$        pos = tmp_pos
!!$        deallocate(tmp_pos)
!!$        ierror = saddle_push_log()
!!$
!!$        ! If the push went well, converge to the new minimum.
!!$        if (push_error == -1) then
!!$           success = .false.
!!$           delta_e = total_energy - ref_energy
!!$        else
!!$           call min_converge( success )     ! And we converge to the new minimum.
!!$           delta_e = total_energy - ref_energy
!!$           call displacement( posref, pos, delr, npart)
!!$           if ( iproc == 0 ) then
!!$              write(unit6P,*) 'BART: CHECK CONNECTIVITY :  Delta_E = ',delta_e, '  delr: ', delr
!!$              open( unit = FLOG, file = LOGFILE, status = 'unknown',&
!!$              & action = 'write', position = 'append', iostat = ierror )
!!$              write(FLOG,"(' ','CONNECTIVITY  |E(con-ini)= ', f9.4,&
!!$            & '  |npart= ', i4,' |delr= ', f8.3)")&
!!$            &  delta_e,npart, delr
!!$            close(FLOG)
!!$
!!$           end if
!!$        end if
!!$
!!$         fname = CHECK // scounter
!!$         call store( fname )
!!$
!!$
!!$         pos = pos_min 
!!$         total_energy = energy_min 
!!$         delta_e = total_energy - ref_energy
!!$
!!$         deallocate(pos_min)
!!$         deallocate(pos_saddle)
!!$         deallocate(force_saddle)
!!$         deallocate(projection_saddle)
    endif


  end subroutine art_search


  subroutine report_and_check(fname,accept)

    use lanczos_defs, only: LANCZOS_MIN

    implicit none

    character(20), intent(in) :: fname
    character(8), intent(in)  :: accept


    integer       :: npart ! Number of atoms participating to the event
    real(kind=8)  :: delr
    real(kind=8)  :: difpos
    integer       :: ierror

    ! Magnitude of the displacement (utils.f90).

    call displacement( posref, pos, delr, npart )

    if ( iproc == 0 ) then           ! Report
       close(FLIST)
       open( unit = FLOG, file = LOGFILE, status = 'unknown',&
            & action = 'write', position = 'append', iostat = ierror )
       write(unit6P,*) 'BART: Configuration stored in file ',fname
       write(FLOG,'(1X,A34,A17)') ' - Configuration stored in file : ', trim(fname)
       write(FLOG,'(1X,A34,(1p,e17.10,0p))')&
            &  ' - Total energy Minimum (eV)    : ', total_energy
       write(FLOG,"(' ','MINIMUM',i5, a9,' |E(fin-ini)= ', f9.4,' |E(fin-sad)= ',&
            & f9.4,' |npart= ', i4,' |delr= ', f8.3,' |evalf=', i6,' |')")&
            & mincounter, adjustr(accept), delta_e,                                 &
            & total_energy - saddle_energy, npart, delr, evalf_number
       write(unit6P,"(' ','BART: MINIMUM',i5, a9,' |E(fin-ini)= ', f9.4,' |E(fin-sad)= ',&
            & f9.4,' |npart= ', i4,' |delr= ', f8.3,' |evalf=', i6,' |',f8.3,3f7.2)")       &
            & mincounter, adjustr(accept), delta_e,                              &
            & total_energy - saddle_energy, npart, delr, evalf_number, difpos,   &
            & a1, b1, c1

       close(FLOG)
    end if
    ! Is a real minimum ?
    ! If dual_search, this check is done only on
    ! one system
    if ( LANCZOS_MIN .and. success ) call check_min( 'M' )

    if ( eventtype == "REFINE_AND_RELAX" .or. eventtype == "GUESS_DIRECTION" ) call end_art()

    if (accept=='ACCEPTED' .or. WRITE_REJECTED_EVENT) mincounter = mincounter + 1

    if ( iproc == 0 ) then
       open(unit=FCOUNTER,file=COUNTER,status='unknown',action='write',iostat=ierror)
       write(FCOUNTER,'(A12,I6)') 'Counter:    ', mincounter
       close(FCOUNTER)
    end if

    evalf_number = 0

  end subroutine report_and_check


end module run_art
