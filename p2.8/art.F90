!> @file
!! @author
!!    Copyright (C) Normand Mousseau, June 2001
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> ART Program art90
!! Main program to use BigDFT with art nouveau method
module art_mod
  use ndm2art2ndm,only:set_pointers_art,lchg,pscart,init_mpi_art2,parapath
  use min_converge_mod,only:min_converge
  use random,only:ran3
  use run_art,only:init_conf,art_search,report_and_check,accept,success
  use write_refconfig_mod,only:write_refconfig
  use storage,only:store
  use end_art_mod,only:end_art
  USE atomconfig,only : atom_config
  USE cellconfig, only:cell_config
  USE boxconfig,only:box_config
  use Tpara,only:para_space_config
  use defs, only : conf_final,conf_initial,conf_saddle,conf_final,eventslist,ievent,ievent_restart,iproc,nproc,&
       &flist,local_force,mincounter,natoms,number_events,pos,posref,ref_energy,refcounter,restart,scala,scalaref,&
       temperature,total_energy,use_local_forces,WRITE_REJECTED_EVENT
contains
  subroutine art90(atdml,celndm,boxndm,psc)



    implicit None

    type(atom_config)::atdml
    type(cell_config)::celndm
    type(box_config)::boxndm
    
    type(para_space_config)::psc

    integer       :: ierror
    real(kind=8)  ::  random_number
!    real(kind=8)  :: ran3, random_number
    character(20) :: fname
    logical       :: local_success

    call init_mpi_art2(atdml,celndm,boxndm,psc,parapath) ! most is done in init_mpi_art
    NATOMS= atdml%im
    restart=.false.
    call init_conf(atdml,celndm,boxndm)

    !initialization of local and workers 

    ! _________
    !                MAIN LOOP OVER THE EVENTS.

    Do_ev: do ievent = ievent_restart, NUMBER_EVENTS

       call art_search(fname)

       ! Now, we accept or reject this move based
       ! on a Boltzmann weight.
       if ( iproc == 0 ) then
          random_number = ran3( )
          open( unit = FLIST, file = EVENTSLIST, status = 'unknown', &
               & action = 'write', position = 'append', iostat = ierror )
       end if

       If_bol: if ( ( (total_energy - ref_energy) < -temperature * log( random_number ) )&
            & .and. ( temperature >= 0.0d0 ) .and. success ) then
          accept = "ACCEPTED"
          if (LOCAL_FORCE) then
             use_local_forces = .false.
             call min_converge( local_success ) 
             write(*,*) "Global minimisation -  success: ", local_success, " - total_energy: ", total_energy
             write(*,*) "Rewrite the minimisation file ", conf_final
             call store( conf_final )
          endif
          if ( iproc == 0 )&
               &  write(FLIST,*) conf_initial, conf_saddle, conf_final,'    accepted'
          ! We now redefine the reference configuration
          scalaref     = scala
          posref       = pos
          conf_initial = conf_final
          ref_energy   = total_energy
          refcounter   = mincounter

          if ( iproc == 0 ) call write_refconfig( )
       else
          ! Else If_bol:
          ! If the event is not accepted we start
          ! from the previous refconfig.

          accept = "REJECTED"
          if (.not.WRITE_REJECTED_EVENT) NUMBER_EVENTS = NUMBER_EVENTS + 1
          if ( iproc == 0 ) then        ! Write

             ! the exchange does not have any meaning without a geometric analysis.
             !if (( total_energy - ref_energy ) > 1.0d-5 )  then
             write(FLIST,*) conf_initial, conf_saddle, conf_final,'    rejected'
             !else
             !    write(FLIST,*) conf_initial, conf_saddle, conf_final,'    exchanged'
             !end if
          end if

       end if If_bol

       call report_and_check(fname,accept)

    end do Do_ev

    call end_art

  END subroutine art90

  
end module art_mod
