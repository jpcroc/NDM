module print_proj_mod
  use convert_to_chain_mod,only:convert_to_chain

contains
  !> ART print_proj
subroutine print_proj( repetitions, stage, vector, eigenvalue, stepsize )

  use defs

  implicit none
  integer,          intent(in) :: repetitions
  character(len=1), intent(in) :: stage
  real(kind=8),     intent(in) :: eigenvalue
  real(kind=8),     intent(in), dimension(3*natoms) :: vector
  real(kind=8),     intent(in) :: stepsize

  !Local variables

  integer :: ierror
  integer :: i
  real(kind=8), dimension(3) :: boxl
  real(kind=8), allocatable  :: pc(:,:)
  character(len=40) :: fname
  character(len=4)  :: rcounter

  allocate(pc(3, NATOMS))
  do i = 1, NATOMS, 1
     pc(:, i) = (/ vector(i), vector(natoms + i), vector(2 * natoms + i) /)
  end do

  if ( iproc == 0 ) then

     boxl = box * scala  ! Update the box size

     call convert_to_chain( repetitions, 2, rcounter )
     fname = 'proj_'//trim(rcounter)//'_'//stage//".xyz"
     fname = trim(fname)

     open(unit=XYZ,file=fname,status='unknown',action='write',iostat=ierror)

     write(XYZ,*) NATOMS,  'angstroem'

     if (boundary == 'P') then
        write(XYZ,'(a,3(1x,1p,e24.17,0p),2x,a,2x,F12.6)')'periodic', (boxl(i),i=1,3),'#', eigenvalue
     else if (boundary == 'S') then
        write(XYZ,'(a,3(1x,1p,e24.17,0p),2x,a,2x,F12.6)')'surface',  (boxl(i),i=1,3),'#', eigenvalue
     else if (boundary == 'T') then
        write(XYZ,'(a,3(1x,1p,e24.17,0p))')'triclinic', cell(:,1)
        write(XYZ,'(3(e24.17))') cell(:,2)
        write(XYZ,'(3(e24.17))') cell(:,3)
     else
        write(XYZ,*)'free ',' # ', eigenvalue
     end if

     do i= 1, NATOMS
        write(XYZ,'(1x,A2,3(2x,f16.8),2x,A,2x,I3,3(2x,f12.8))') &
             &                               Atom(i), x(i) + stepsize*pc(1,i), &
             &                                        y(i) + stepsize*pc(2,i), &
             &                                        z(i) + stepsize*pc(3,i), &
             &                               '#',i, pc(1,i), pc(2,i),pc(3,i)
     end do

     close(XYZ)
  end if

  deallocate(pc)

END SUBROUTINE print_proj

end module print_proj_mod
