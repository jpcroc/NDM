module write_refconfig_mod
contains
  subroutine write_refconfig( )

  use defs
  implicit none

  !Local variables
  integer :: i, ierror
  real(kind=8),dimension(3) :: boxl

  boxl = box * scala                  ! Update the box size

  ! switch replace for unknown
  open(unit=FREFCONFIG,file=REFCONFIG,status='unknown',action='write',iostat=ierror)
  write(FREFCONFIG,*) 'run_id: ', refcounter
  write(FREFCONFIG,*) 'total_energy: ', total_energy
  if (.not. boundary == 'T') then
     write(FREFCONFIG,*) boundary, boxl
  else
     write(FREFCONFIG,*) boundary, cell(:,1)
     write(FREFCONFIG,*) ' ', cell(:,2)
     write(FREFCONFIG,*) ' ', cell(:,3)
  end if

  do i = 1, NATOMS
     write(FREFCONFIG,'(1x,i6, 3(2x,F16.8))') typat(i), x(i), y(i), z(i)
  end do
  close(FREFCONFIG)

END SUBROUTINE write_refconfig
end module write_refconfig_mod
