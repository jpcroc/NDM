module print_event_mod
contains

  subroutine print_event( ievent_current, temperat )

  use defs
  implicit none

  !Arguments
  integer, intent(in) :: ievent_current
  real(kind=8), intent(in) :: temperat

  !Local variables
  integer :: ierror

  write(*,*) 'BART: Simulation : ', ievent_current
  write(*,*) 'BART: Attempt    : ', atp
  write(*,*) 'BART: Starting from minconf : ', refcounter
  write(*,*) 'BART: Reference Energy (eV) : ', ref_energy
  write(*,*) 'BART: Temperature : ', temperat

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(FLOG,*) ' __________________________________________________'
  write(FLOG,'(1X,A34,I17)') ' - Simulation                   : ', ievent_current
  write(FLOG,'(1X,A34,I17)') ' - Attempt                      : ', atp
  write(FLOG,'(1X,A34,I17)') ' - Starting from minconf        : ', refcounter
  write(FLOG,'(1X,A34,(1p,e17.10,0p))') ' - Reference Energy (eV)        : ', ref_energy
  write(FLOG,'(1X,A34,F17.6)') ' - Temperature                  : ', temperat
  close(FLOG)

END SUBROUTINE print_event
end module print_event_mod
