module write_step_report
  use force_projection_mod,only:force_projection_art
  use convert_to_chain_mod,only:convert_to_chain
  use storage,only:store
  use disp,only:displacement
  use save_intermediate_mod,only:save_intermediate
  use lanczos_mod,only:lanczos

contains
  subroutine write_step ( stage, it, a1, energy )

   use defs
   use saddles
   use lanczos_defs
   use lanczos_mod,only:lanczos
   implicit None

   !Arguments
   character(len=1), intent(in) :: stage
   integer,          intent(in) :: it
   real(kind=8),     intent(in) :: a1
   real(kind=8),     intent(in) :: energy

   !Local variables
   integer :: ierror
   character(len=2)  :: etape

   etape = stage//"="

   if ( iproc == 0 ) then
      open( unit = FLOG, file = LOGFILE, status = 'unknown',&
         & action = 'write', position = 'append', iostat = ierror )

      write( FLOG,'(i4,2x,a2,i4,1x,F16.6,i3,i2,4f12.4,1x,1f10.3,i4,i6,f7.2)')  &
         &   pas, etape, it, delta_e, m_perp, try, ftot, fpar, fperp,  &
         &  eigenvalue, delr, npart, evalf_number, a1
      close( FLOG )

      write(unit6P,'(a,i4,2x,a2,i4,1x,(1p,e17.10,0p),1x,2i3,4f12.6,1x,1f10.4,i4,i6)') &
         &   " BART:", pas, etape, it, delta_e,  m_perp, try, ftot, fpar, fperp,   &
         &   eigenvalue, delr, npart, evalf_number
   end if

END SUBROUTINE write_step


!> ART end_report
subroutine end_report ( success, ret, saddle_energy )

   use defs
   use saddles
   use lanczos_defs


   implicit None

   !Arguments
   logical, intent(out)   :: success
   integer, intent(inout) :: ret
   real(kind=8), intent(in) :: saddle_energy

   !Local variables
   integer :: i
   integer :: ierror                   ! File control.
   real(kind=8) :: a1
   real(kind=8), dimension(VECSIZE) :: perp_force
   logical :: new_projection           ! For lanczos.
   character(len=10) :: converg        ! For report.
   character(len=4)  :: scounter
   character(len=20) :: fname
   ! __________________

   select case( ret )

   case( 200000 : 299999 )                ! ftot < EXITTHRESH

      If_diis: if ( USE_DIIS ) then
         fpar  = 0.0d0
         fperp = 0.0d0

         If_check: if ( DIIS_CHECK_EIGENVEC ) then
            ! Lanczos several times.
            new_projection = .false.
            Do_lanc: do i = 1, 4
               call lanczos( NVECTOR_LANCZOS_C, new_projection, a1 )
               new_projection = .false.
               ! Exit of loop.
               if ( eigenvalue < 0.0d0 ) exit Do_lanc
            end do Do_lanc

            if ( eigenvalue >= 0.0d0 ) then
               ret = 600000 + pas
            else                            ! Else of eigenvalue.
               ! New fpar and fperp.
               call force_projection_art( fpar, perp_force, fperp,&
                  &   ftot, force, projection )
               ret = 100000 + pas
            end if

         else                               ! Else of If_check
            eigenvalue = 0.0d0
            ret = 100000 + pas
         end if If_check

      else                                  ! Pure Lanczos
         ! Else of If_diis
         if ( eigenvalue > 0.0d0 ) then
            ret = 600000 + pas
         else
            ret = 200000 + pas
         end if
      end if If_diis

   case( 400000 : 499999 )                ! clean_wf
      ! We accept the event is the new total
      ! force changes only in max 0.1 eV
      if ((ftot -0.1d0) < EXITTHRESH) then
         ! we check always the eigenvalue
         new_projection = .false.
         Do_lanc_c: do i = 1, 4
            call lanczos( NVECTOR_LANCZOS_C, new_projection, a1 )
            new_projection = .false.
            ! Exit of loop.
            if ( eigenvalue < 0.0d0 ) exit Do_lanc_c
         end do Do_lanc_c

         if ( eigenvalue >= 0.0d0 ) then
            ret = 600000 + pas
         else                          ! Else of eigenvalue.
            ! New fpar and fperp.
            call force_projection_art( fpar, perp_force, fperp,&
               &   ftot, force, projection )
            ret = 300000 + pas
         end if
      end if

   end select

   select case ( ret )

   case( 100000 : 399999 )

      converg = 'CONVERGED'
      success = .true.
      ! We write the configuration in a sad.... file
      call convert_to_chain( mincounter, 4, scounter )
      fname = SADDLE // scounter
      if (iproc == 0 ) call store( fname )
      conf_saddle = fname

   case default
      converg = 'FAILED'
      success = .false.
   end select

   delta_e = saddle_energy - ref_energy
   call displacement( posref, pos, delr, npart )
   ! Final report of saddle point
   if ( iproc == 0 ) then
      open( unit = FLOG, file = LOGFILE, status = 'unknown',&
         & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,"(/' ','SADDLE',i5, a10,' |ret ',i6,' |delta energy= '," //  &
         &       "f9.4, ' |force_(tot,par,perp)= ', 3f10.4," // &
         &       "' |eigenval=',f9.4,' |npart= ', i4,' |delr= ', f8.3,' |evalf='," // &
         &       "i6, ' |')")                                          &
         & mincounter, adjustr(converg), ret, delta_e, ftot,  &
         & fpar, fperp, eigenvalue, npart, delr, evalf_number
      write(unit6P,"(/' ','BART: SADDLE',i5, a10,' |ret ',i6,' |delta energy= '," //  &
         &    "f9.4, ' |force_(tot,par,perp)= ', 3f10.4," //       &
         &    "' |eigenval=',f9.4,' |npart= ', i4,' |delr= ', f8.3,' |evalf='," // &
         &    "i6, ' |')")                                          &
         & mincounter, adjustr(converg), ret, delta_e, ftot,  &
         & fpar, fperp, eigenvalue, npart, delr, evalf_number

      if ( success ) then
         ! Write
         write(unit6P,*) 'BART: Configuration stored in file ',fname
         write(FLOG,'(1X,A34,A17)') ' - Configuration stored in file : ', trim(fname)
         write(FLOG,'(1X,A34,(1p,e17.10,0p))') &
            &   ' - Total energy Saddle (eV)     : ', saddle_energy
         write(FLOG,*) ' '
      end if
      close(FLOG)
   end if

END SUBROUTINE end_report
end module write_step_report
