!> @file
!!   Subroutine which calls the right for type inside art
!! @author
!!   Written by Laurent Karim Beland, UdeM 2011!!
!!   Copyright (C) 2010-2011 BigDFT group, Normand Mousseau
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS


!> Subroutine which calls the right for type inside art


! $Revision$
! $Date: 2015-12-14 11:01:46 -0500 (Lun, 14 déc 2015) $
! $Id: calcforce.f90 1534 2015-12-14 15:59:38Z mickael $

module calcforce_mod
  use ndm2art2ndm,only:calcforce_ndm
contains
subroutine calcforce(nat, posa, boxl, forca, energy, evalf_number)

  use defs, only : typat, use_local_forces, local_ref_energy, global_ref_energy
  implicit none

  !Arguments
  integer,      intent(in)                            :: nat
  real(kind=8), intent(in),  dimension(3*nat)         :: posa
  real(kind=8), dimension(3), intent(inout)           :: boxl
  real(kind=8), intent(out), dimension(3*nat)         :: forca
  real(kind=8), intent(out)                           :: energy
  integer,      intent(inout)                         :: evalf_number

!!$  interface
!!$     subroutine SWcalcforce(nat,typat,posa,boxl,forca, energy)
!!$       integer,      intent(in)                            :: nat
!!$       integer,      intent(in),dimension(nat)             :: typat
!!$       real(kind=8), intent(in),  dimension(3*nat)         :: posa
!!$       real(kind=8), dimension(3), intent(inout)           :: boxl
!!$       real(kind=8), intent(out), dimension(3*nat)         :: forca
!!$       real(kind=8), intent(out)                           :: energy
!!$     end subroutine SWcalcforce
!!$
!!$#ifdef LAMMPS_VERSION
!!$     subroutine calcforce_lammps(nat,posa,boxl,forca, energy)
!!$       integer,      intent(in)                            :: nat
!!$       real(kind=8), intent(in),  dimension(3*nat)         :: posa
!!$       real(kind=8), dimension(3), intent(inout)           :: boxl
!!$       real(kind=8), intent(out), dimension(3*nat)         :: forca
!!$       real(kind=8), intent(out)                           :: energy
!!$     end subroutine calcforce_lammps
!!$#endif
!!$
!!$  end interface
!!$
!!$  if(use_local_forces) then
!!$#ifdef LAMMPS_VERSION
!!$     if(energy_type == "LAM") then
!!$        call calcforce_local_lammps(nat,posa,boxl,forca,energy)
!!$
!!$        ! We now correct the potential energy for a global scale
!!$        energy = energy - local_ref_energy + global_ref_energy
!!$        return
!!$     endif
!!$#endif
!!$
!!$     write(*,*) "Your choice of potential is not set up for local force calculations, please change this choice and restart"
!!$     stop
!!$
!!$  else
!!$
!!$     if (energy_type == "SWP" .or. energy_type == "SWA")  then
!!$        call SWcalcforce(nat,typat,posa,boxl,forca, energy )
!!$        evalf_number = evalf_number +1
!!$     endif
!!$     ! -----------------------------------------------------
!!$
!!$#ifdef NDM_VERSION
!     if (energy_type=="NDM") then
        call calcforce_NDM(nat,posa,forca, energy)
        evalf_number = evalf_number +1
!     end if
!#endif
!!$     
!!$#ifdef LAMMPS_VERSION
!!$     if (energy_type == "LAM")  then
!!$        call calcforce_lammps(nat,posa,boxl,forca,energy)
!!$        evalf_number = evalf_number +1
!!$     endif
!!$#endif
!!$
!!$  endif
  ! -----------------------------------------------------

END SUBROUTINE calcforce
end module calcforce_mod

! -------------------------------------------------------------------------------------------
!  This routine update the atom list for the local force calculation.
!
!  At the moment, this function only works with lammps.
!
! Copyright  Normand Mousseau (December 2015)
!            Universite de Montreal
! 
!!$subroutine update_local_potential(nat,posart,boxl,tmp_force,temp_energy)
!!$  use defs
!!$  implicit none
!!$
!!$  integer,      intent(in)                            :: nat
!!$  real(kind=8), intent(in),  dimension(3*nat)         :: posart
!!$  real(kind=8), dimension(3), intent(inout)           :: boxl
!!$  real(kind=8), intent(out), dimension(3*nat)         :: tmp_force
!!$  real(kind=8), intent(out)                           :: temp_energy
!!$
!!$  interface
!!$
!!$#ifdef LAMMPS_VERSION
!!$     subroutine update_local_lammps(nat,posa,boxl,forca, energy)
!!$       integer,      intent(in)                            :: nat
!!$       real(kind=8), intent(in),  dimension(3*nat)         :: posa
!!$       real(kind=8), dimension(3), intent(inout)           :: boxl
!!$       real(kind=8), intent(out), dimension(3*nat)         :: forca
!!$       real(kind=8), intent(out)                           :: energy
!!$     end subroutine update_local_lammps
!!$#endif
!!$
!!$  end interface
!!$
!!$  
!!$#ifdef LAMMPS_VERSION
!!$  if(energy_type == "LAM")  then
!!$     write(*,*) 'Call update_local_lammps'
!!$     call update_local_LAMMPS(nat,posart,boxl,tmp_force,temp_energy)
!!$     return
!!$  endif
!!$#endif
!!$
!!$  write(*,*) "Your choice of potential is not set up for local force calculations, please change this choice and restart"
!!$  stop
!!$
!!$end subroutine update_local_potential

