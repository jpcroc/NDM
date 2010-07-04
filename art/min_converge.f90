! Subroutine min_converge
!
! Minimizes the energy at constant volume. It uses a steepest descent
! algorithm which is fast and precise enough for our needs.
!
! This minimization is done with only a minimal knowledge of the physics
! of the problem so that it is portable
!
! The minimization uses a simple steepest descent with variable step size.
!
! This file contains 2 routines
!
! June 2001  Normand Mousseau
!

MODULE minimization
!  This module defines a number of parameters used during the minimization
  implicit none
  save

  integer, parameter :: MAX_ITER = 8000
!orig  real(8), parameter :: FTHRESHOLD = 1.0D-1
!  real(8), parameter :: FTHRESHOLD = 5.0D-3  !This is for SiC
  real(8), parameter :: FTHRESHOLD = 1.0D-3  !This is for SiC
  real(8), parameter :: STEPSIZE = 0.0001  ! Size in angstroems

  real(8), parameter :: FTHRESH2 = FTHRESHOLD * FTHRESHOLD
END MODULE minimization


subroutine min_converge
  use defs
  use minimization
  use art_in_ndm_module
  implicit none

  integer :: iter, i, npart,ierror
  real(8) :: current_energy, ftot,ftot2, step, boxl(3), delr
  real(8), dimension(:),allocatable :: posb, forceb
  ! local variables .... by me
  real(8)                       ::   posprov
  real(double), dimension(VECSIZE) :: vit_art
  real(double), dimension(imm)     :: ftot_max
  integer               :: iatom

  vit_art(:)=0.d0
  it_art=0

  allocate(posb(VECSIZE))
  allocate(forceb(VECSIZE))
  posb(:)=pos(:)

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)

  ! We compute at constant volume

  boxl(:) = box(:) *scala
  call calcforce(NATOMS,type,pos,boxl,force,total_energy)
  evalf_number = evalf_number + 1
  current_energy = total_energy
  write(*,*)  'current energy :', total_energy

  step = STEPSIZE   ! which is in fact Initial_Step_Size
  do iter = 1, MAX_ITER

     do i=1,VECSIZE
        ! 
        if (vit_art(i)*force(i)>0) then
           posprov = 2.d0*posb(i) - pos(i) + tmass_art(i)*force(i)
        else 
           posprov = posb(i) + tmass_art(i) * force(i)
        end if
        !
        vit_art(i) = (posprov-pos(i))*usdh_art
        pos(i) = posb(i)
        posb(i) = posprov
        !
     end do

     call calcforce(NATOMS,type,posb,boxl,forceb,total_energy)
     evalf_number = evalf_number + 1

     do i=1, NATOMS
        ftot_max(i) = forceb(i)**2 + forceb(i+NATOMS) **2 +forceb(i+2*NATOMS)**2
     end do
     ftot2 = MAXVAL(ftot_max(1:NATOMS))


     force(:) = forceb(:)
     current_energy = total_energy
     call displacement(posref, posb, delr,npart)

     if (print_details .and. mod(iter,mprint) == 0 ) then
        write(*, "(' ','it: ',i5,' ener: ', f12.4,' ftot: ', f12.6, ' step: ', &
             & f12.6,' evalf: ',i4,' delr: ', f12.6,' npart: ', i4)") iter, &
             & total_energy, sqrt(ftot2),tstep_art, evalf_number, delr, npart 
        write(FLOG, "(' ','it: ',i5,' ener: ', f12.4,' ftot: ', f12.6, ' step: ', &
             & f12.6,' evalf: ',i4,' delr: ', f12.6,' npart: ', i4)") iter, &
             & total_energy, sqrt(ftot2),tstep_art, evalf_number, delr, npart 
     endif

     if(ftot2 < FTHRESH2) then
        pos(:)=posb(:)
        exit
     end if
  end do

  ftot = sqrt(ftot2)  
  if (ftot < FTHRESHOLD ) then
     write(*,*) 'Minimization successful   ftot : ', ftot
     write(*,*) 'total energy',total_energy
  else
     write(*,*) 'Minimization failed   ftot : ', ftot
  endif

  close(flog)
  deallocate(posb)
  deallocate(forceb)
end subroutine min_converge
