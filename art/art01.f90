subroutine art90
!-----------------------------------------------

      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use tab_imm_m
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use DEFS
      use random_art
      use lanczos_defs
      use art_in_ndm_module
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      implicit none
!      integer  :: ielat(imm)
!      integer  :: iwmax(imm)
!      integer  :: ityp(imm)
!      real(double)  :: xp(3,imm)
!      real(double)  :: xpp(3,imm)
!      real(double)  :: vp(3,imm)
!      real(double)  :: ax(3,imm)
!      real(double)  :: fp(3,imm)




! This is the main program for ART nouveau version 2001
! 
! This version is made to work with SIESTA 2001
!
! Copyright Normand Mousseau, July 2001

  integer :: i, ierror
  integer :: npart             ! Number of atoms participating to the eventt
  integer :: seed
  real(8) :: del_r
  real(8) :: ran3

  logical :: success
  character(len=20) :: fname
  character(4) :: scounter
  character(len=150) :: format_output
  real(8), dimension(:), allocatable, target  :: diff


  write(6,*)
  write(6,*)
  write(6,*)'************ DEBUT DE ART ****************'
  write(6,*)
  write(6,*)

  call read_parameters()
  allocate(diff(VECSIZE))

  ! We now decide whether or not we restart if the restart exists, then we restart

  write(*,*) 'Start with new event'
  call allocate_art()
  call init_art_in_ndm_module(xp, xpp, vp, ax, fp, ielat, iwmax, ityp) 
  call initialize()       ! Read in the initial configuration
  ievent_restart = 1

  do ievent= ievent_restart, NUMBER_EVENTS         ! Main loop over the events
    call print_newevent(ievent,temperature)
    ! We look for a local saddle point b

    ! If it is a restart event for the activation, phase 1 or 2, or it is not a restart event
    ! then we call find_saddle
    if ( .not. (restart .and. (state_restart .eq. 3))  ) then
       do 
          call find_saddle( success )
          if ( success ) exit
       end do
    endif

    ! If not a new event, then a convergence to the saddle point, We do not go further
    if( .not.NEW_EVENT) stop

    ! Push the configuration slightly over the saddle point in order to minimise 
    ! the odds that it falls back into its original state
    diff = pos-posref
    del_r = sqrt(dot_product(diff,diff))*sign(1.0,DOT_PRODUCT(diff,projection))
    
    if (dabs(del_r)>= 0.4 ) del_r=0.4*sign(1.0,del_r)
!    write(*,*) 'GO TO THE VALEY', del_r
!ooo    pos = pos + 0.15 * del_r * projection ! A vectorial operation
   pos = pos + del_r * projection ! A vectorial operation

    ! And we converge to the new minimum.
    call min_converge()

    ! We need to write the configuration in a min.... file
    call convert_to_chain(mincounter,scounter)
    write(*,*) ' Mincounter is : ', mincounter, ' and scounter is ', scounter
    fname =   FINAL // scounter


! Compute the displacement and the number of atoms involved in the event
     
     do i = 1,NATOMS
!         Atom(i) = 'Si'
         Atom(i) = ty(ityp(i))
     enddo

    call displacement(pos,posref,del_r,npart)
    
    ! We now store the configuration into fname
    call store(fname)

    open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
    open(unit=FLIST,file=EVENTSLIST,status='unknown',action='write',position='append',iostat=ierror)
    write(*,*) 'Configuration stored in file ',fname
    conf_final = fname

    mincounter = mincounter+1
    write(6,*)'DIFF ENERGIE FINAL-INITIAL ', total_energy - ref_energy
    ! Now, we accept or reject this move based on a Boltzmann weight
    if(  (total_energy - ref_energy) < -temperature * log(ran3()) ) then
      write(*,*) 'New configuration accepted, mincounter was : ', mincounter-1
      write(FLOG,*) 'New configuration accepted, mincounter was : ', mincounter-1

      ! We write out various information to both screen and file
      format_output = "('Total energy: ', f12.4, ' npart: ',i4, ' delr: ', f9.4, ' Nr. Force Eval: ', i10, ' accepted')  "
      write(*,format_output) total_energy, npart, del_r, evalf_number
      write(FLOG,format_output) total_energy, npart, del_r, evalf_number

      write(FLIST,*) conf_initial, conf_saddle, conf_final,'    accepted'

      ! We now redefine the reference configuration
      scalaref = scala
      posref= pos          ! This is a vectorial copy
      conf_initial = conf_final
      
      ref_energy = total_energy

      ! Update the reference configuration file, which serves as the initial
      ! configuration for events.
      call write_refconfig() 
    else
      write(*,*) 'New configuration rejected, mincounter was : ', mincounter-1
      write(FLOG,*) 'New configuration rejected, mincounter was : ', mincounter-1

      ! We write out various information to both screen and file
      format_output = "('Total energy: ', f12.4, ' npart: ',i4, ' delr: ', f9.4, ' Nr. Force Eval: ', i10, ' rejected')  "
      write(*,format_output) total_energy, npart, del_r, evalf_number
      write(FLOG,format_output) total_energy, npart, del_r, evalf_number

      ! The events is not accepted; we start from the previous refconfig
      if((total_energy - ref_energy) > 1.0d-5)  then
        write(FLIST,*) conf_initial, conf_saddle, conf_final,'    rejected'
      else  
        write(FLIST,*) conf_initial, conf_saddle, conf_final,'    exchanged'
      endif
    endif
    close(FLIST)
    close(FLOG)

    open(unit=FCOUNTER,file=COUNTER,status='unknown',action='write',iostat=ierror)
    write(FCOUNTER,'(A12,I6)') 'Counter:    ', mincounter
    close(FCOUNTER)
  end do
end subroutine art90


! This subroutine prints the initial details for a new events
subroutine print_newevent(ievent_current,temperat)
  use defs
  implicit none
  integer, intent(in) :: ievent_current
  real(8), intent(in) :: temperat
  integer :: ierror;

  write(*,*) 'Simulation : ', ievent_current
  write(*,*) 'Starting from minconf : ', mincounter
  write(*,*) 'Temperature : ', temperat

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(FLOG,*) 'Simulation : ', ievent_current
  write(FLOG,*) 'Starting from minconf : ', mincounter
  write(FLOG,*) 'Temperature : ', temperat
  close(FLOG)

  return
end subroutine print_newevent
