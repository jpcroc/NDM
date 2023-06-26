module saddle_converge_mod
!> @file
!! @author
!!    Copyright (C) Normand Mousseau, June 2001
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> ART saddle_converge
!!   This subroutine brings the configuration to a saddle point. It does that
!!   by first pushing the configuration outside of the harmonic well, using
!!   the initial direction selected in find_saddle. Once outside the harmonic
!!   well, as defined by the appearance of a negative eigenvalue (or
!!   reasonnable size) the configuration follows the direction corresponding
!!   to this eigenvalue until the force components parallel and perpendicular
!!   to the eigendirection become close to zero.

! $Revision$
! $Date: 2015-12-14 11:01:46 -0500 (Lun, 14 déc 2015) $
! $Id: saddle_converge.f90 1534 2015-12-14 15:59:38Z mickael $

 use min_converge_mod,only:check_min
 use save_intermediate_mod,only:save_intermediate


  use art_step_mod,only:allocate_activation,deallocate_activation,apply_diis,apply_lanczos
  use calcforce_mod,only:calcforce
  use end_art_mod,only:end_art
  use save_restart,only:save_state
  use disp,only:displacement
  use force_projection_mod,only:force_projection_art
  use write_step_report,only:write_step
  use lanczos_mod,only:lanczos

contains
  subroutine saddle_converge( ret, saddle_energy )

    use defs
    use saddles
    use lanczos_defs
    use diis_defs
   implicit none

   !Arguments
   integer, intent(out) :: ret
   real(kind=8), intent(out) :: saddle_energy  ! Energy at saddle point.

   !Local variables
   logical :: new_projection              ! For lanczos.
   ! Loop indeces :
   integer :: i, kter, kter_init, liter, diter, step_rejected
   integer :: ierror, ierr                ! File and MPI control.

   real(kind=8) :: step                        ! This is the step in the hyperplane.
   real(kind=8),dimension(3) :: boxl
   real(kind=8) :: a1
   real(kind=8) :: current_energy              ! Accepted energy.
   real(kind=8) :: ftot_b                      ! ftot  for evaluation.
   real(kind=8) :: fpar_b                      ! fpar  for evaluation.
   real(kind=8) :: fperp_b                     ! fperp for evaluation.
   real(kind=8), dimension(VECSIZE) :: pos_b        ! Position for evaluation.
   real(kind=8), dimension(VECSIZE) :: force_b      ! Total Force for evaluation.
   real(kind=8), dimension(VECSIZE) :: perp_force   ! Perpendicular force...
   real(kind=8), dimension(VECSIZE) :: perp_force_b ! ...& for evaluation.
   ! __________________
   ! Write header in log.file.
   if ( iproc == 0 ) then
      open( unit = FLOG, file = LOGFILE, status = 'unknown',&
         & action = 'write', position = 'append', iostat = ierror )
      write(FLOG, '(a23,a8,a8,a12,a12,a12,a11,a7,a6,a5)')&
         &   'E-Eref','m_perp','ftot','fpar','fperp','eigen','delr','npart','evalf','a1'
      write(FLOG, '(a23,a31,a27)' ) '( eV )', '( eV/Ang )', '( eV/Ang**2 )'
      close(FLOG)
   end if
   ! initialization
   saddle_energy = 0.0d0
   boxl = box * scala                  ! We compute at constant volume.
   cw_try_again = .True.               ! For clean_wf

   if ( .not. restart ) pas = 0
   ! If restarts in harmonic well
   if ( restart .and. ( state_restart == 1 ) ) then

      kter_init = iter_restart         ! init value of kter loop.
      initial_direction = direction_restart
      call displacement( posref, pos, delr, npart )
      deallocate(direction_restart)

      restart = .false.
      ! Write
      if ( iproc == 0 ) then
         write(unit6P,*) 'BART: Restart  in harmonic well '
         write(unit6P,*) 'BART: kter : ', kter_init
         write(unit6P,*) 'BART: pos: ', pos(1), pos(2), pos(3)
      end if
   else
      kter_init = 0                    ! init value of kter loop.

   end if
   ! _________
   !                HARMONIC WELL

   If_restart: if ( ( .not. restart ) .and. NEW_EVENT ) then

      eigenvalue = 0.0d0               ! Only for the report.
      a1         = 0.0d0               ! Only for the report.
      step = 0.4*INCREMENT             ! The step in the hyperplane.
      new_projection = .true.          ! We do not use previously computed
      ! lowest direction.
      Do_kter: do kter = kter_init, MAXKTER
         ! Reference energy and force.
!         write(unit6P,*)'JP5',kter,kter_init, MAXKTER
         call calcforce( NATOMS, pos, boxl, force, current_energy, evalf_number )
         ! We now project out the direction of the initial displacement from the
         ! minimum from the force vector so that we can minimize the energy in
         ! the direction perpendicular to the direction of escape from harmonic well

         call force_projection_art( fpar, perp_force, fperp, ftot, force, initial_direction )

         try = 0                       ! Total # of iterationes in a given hyperplane.
         m_perp = 0                    ! Perpendicular iterations accepted.s
         step_rejected = 0             ! Perpendicular steps rejected.

         ! We relax perpendicularly using a simple variable-step steepest descent
         While_perpk: do
!         write(unit6P,*)'JP6'
            pos_b = pos + step * perp_force

            call calcforce( NATOMS, pos_b, boxl, force_b, total_energy, evalf_number )
            ! New force's components.  write(unit6P,*) 'Outside of calcforce_local_lammps'

            call force_projection_art( fpar_b, perp_force_b, fperp_b, ftot_b, &
               &   force_b, initial_direction )

            if ( total_energy < current_energy ) then
               pos            = pos_b
               current_energy = total_energy
               fpar           = fpar_b
               perp_force     = perp_force_b
               fperp          = fperp_b
               ftot           = ftot_b
               force          = force_b

               step = 1.2 * step
               m_perp = m_perp + 1
               step_rejected = 0
            else
               step = 0.6 * step
               step_rejected = step_rejected + 1
            end if
            try = try + 1

            if ( fperp < FTHRESHOLD .or. m_perp >= MAXKPERP &
               &   .or. step_rejected > 5 ) exit While_perpk

         end do While_perpk

         delta_e = current_energy - ref_energy
         ! Magnitude of the displacement (utils.f90).
         call displacement( posref, pos, delr, npart )
!         write(unit6P,*)'JP7'
         if ( SAVE_CONF_INT ) call save_intermediate( 'K' )

         ! We start checking of negative eigenvalues only after a few steps.
!         write(unit6P,*)'JP8 ', kter,kter_min,setup_initial 
         if ( kter == KTER_MIN ) then
            ! First time, twice !!
            do i = 1, 1
!               write(unit6P,*)'JP91'
               call lanczos( NVECTOR_LANCZOS_H, new_projection, a1 )

               new_projection = .false.
               if ( iproc == 0 ) write(unit6P,'(a,3I5,f12.6,f7.2)') &
                  &   'BART COLLINEAR:', pas, kter, i, eigenvalue, a1
            end do

         else if ( setup_initial .and. kter > KTER_MIN + 1 ) then
            call check_min( 'I' )
            call write_step ( 'K', kter, a1, current_energy )
#ifdef MPI_VERSION_ART
!            call MPI_Barrier(MPI_COMM_WORLD, ierr )
#endif
            call end_art( )

         else if ( kter > KTER_MIN  ) then
            ! we get eigen direction for the minimum of this hyperplane.

            call lanczos( NVECTOR_LANCZOS_H, new_projection, a1 )
!            write(unit6P,*)'JP92'
            ! Lanczos call, we start from the
            new_projection = .false.      ! previous direction each time.
         end if

         ! Write
         call write_step ( 'K JP10', kter, a1, current_energy )

!         write(unit6P,*)'JP93',eigenvalue, EIGEN_THRESH ,setup_initial
         ! For restart ( restart.f90 )
         if ( write_restart_file ) then
            state_restart = 1
            total_energy = current_energy
            if ( iproc == 0 ) call save_state( state_restart, kter+1, initial_direction )
         end if
         pas = pas + 1

         ! Is the configuration out of the harmonic basin?
         if ( eigenvalue < EIGEN_THRESH .and. (.not. setup_initial) ) exit Do_kter
         ! If not, we move the configuration along
         ! the initial direction.
         pos = pos + BASIN_FACTOR * INCREMENT * initial_direction

      end do Do_kter

      ! If after MAXKTER iterations we have not found
      ! an inflection the event is killed.
      if ( eigenvalue > EIGEN_THRESH ) then
         saddle_energy = current_energy ! For the report.
         ret = 900000 + kter

         return
      end if
!      write(unit6P,*)'JP200'
      ! The configuration is now out of the harmonic well, we can now bring
      ! it to the saddle point.
      ! First, we must now orient the direction of the eigenvector corresponding to the
      ! negative eigendirection (called projection) such that it points away from minimum.

      fpar = dot_product( force, projection )

      if ( fpar > 0.0d0 ) projection = -1.0d0 * projection

      liter = 1                   ! init value of lanczos loop.
      diter = 1                   ! init value of diis loop.

      switchDIIS= .False.
      ! _________
   else if ( restart .and. &
      &   ( state_restart == 2 .or. state_restart == 4 ) ) then ! Else If_restart

      !           RESTART FROM A PREVIOUS ACTIVATION PROCESS

      call allocate_activation ()
      projection      = direction_restart
      previous_forces = diis_forces_restart
      previous_pos    = diis_pos_restart
      previous_norm   = diis_norm_restart
      maxter     = maxter_r
      eigen_min  = eigen_min_r
      eigenvalue = eigenvalue_r
      nsteps_after_eigen_min = nsteps_after_eigen_min_r
      delta_e = total_energy - ref_energy

      deallocate(direction_restart)
      deallocate(diis_forces_restart)
      deallocate(diis_pos_restart)
      deallocate(diis_norm_restart)

      if ( state_restart == 2 ) then
         liter = iter_restart
         diter = 1
         switchDIIS= .False.           ! We go to apply_lanczos.
         if ( iproc == 0 ) write(unit6P,*) "BART: Restart = 2"
      elseif ( state_restart == 4 ) then
         diter = iter_restart
         liter = 1
         switchDIIS= .True.            ! We go to apply_diis.
         if ( iproc == 0 ) write(unit6P,*) "BART: Restart = 4"
      else                             !DEBUG
         if ( iproc == 0 ) write(unit6P,*) "BART: HOUSTON, we've got a problem"
         call end_art ()
      end if

      call displacement( posref, pos, delr, npart )
      if (iproc==0) write(unit6P,*) "BART: delr npart", delr, npart
      call force_projection_art( fpar, perp_force, fperp, ftot, force, projection )
      ! _________
   else if ( .not. new_event ) then    ! Else If_restart

      !                'REFINE' EVENT
      ! we must compute the eigenvector
      ! with precision.
      new_projection = .true.

      do i = 1, 5

         call lanczos( NVECTOR_LANCZOS_H, new_projection, a1 )

         new_projection = .false.
         if ( iproc == 0 ) write(unit6P,'(a,2I5,f12.6,f7.2)') &
              &   'BART COLLINEAR:', pas, i, eigenvalue, a1
         if ( a1 > collinear_factor ) exit
      end do

      delta_e = total_energy - ref_energy
      call displacement( posref, pos, delr, npart )
      fpar = dot_product(force,projection)
      if( fpar > 0.0d0 ) projection = -1.0d0 * projection
      call force_projection_art( fpar, perp_force, fperp, ftot, force, projection )
      ! This will be just for divide the INCREMENT,
      liter = 10                       ! i.e, a small move.
      ! Write
      call write_step ( 'L', liter, a1, total_energy )

      liter = liter + 1
      diter = 1                        ! init value of diis loop.
      pas   = pas + 1
      switchDIIS= .False.
      ! _________
   else                                ! Else If_restart
      write(unit6P,*) 'BART: Problem with restart and state_restart : '
      write(unit6P,*) 'BART: restart = ', restart, ' state_restart = ', state_restart
      stop
   end if If_restart
   ! _________
   !                ACTIVATION PART
   write(unit6P,*)'JP201 pre ACT'
   if ( .not. restart ) call allocate_activation ()

   ret = 0
   end_activation = .false.

   While_activation: do
!      write(unit6P,*)'JP202',switchDIIS
      if ( .not. switchDIIS ) then
         call apply_lanczos( liter, saddle_energy, ret )
      else
         call apply_diis( diter, saddle_energy, ret )
         if ( diter > 1 ) liter = 1
         diter = 1
      end if

      if ( end_activation ) exit
   end do While_activation

   call deallocate_activation ()
   write(unit6P,*)'JP203endact',end_activation
 END SUBROUTINE saddle_converge



!> ART write_step
end module saddle_converge_mod
