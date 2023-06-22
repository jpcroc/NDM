module save_intermediate_mod
 use convert_to_chain_mod,only:convert_to_chain
contains
  !> ART save_intermediate
!!    It saves the configuration at every step in xyz format.
!!    The name of the file will look like this example:
!!    p_1001_05_030_K.xyz
!!    1001:: is the 'mincounter' of the event,
!!    05  :: the attempt
!!    030 :: the step
!!    K is the argument 'stage' ( K= basin activation, L=lanczos, D= DIIS )
subroutine save_intermediate( stage )

  use defs
  implicit none

  !Arguments
  character(len=1), intent(in) :: stage

  !Local variables
  integer :: i, ierror
  real(kind=8), dimension(3) :: boxl
  character(len=40) :: fname
  character(len=4)  :: scounter, rcounter, pcounter
  character(len=4), dimension(natoms) :: frzchain

  boxl = box * scala                  ! Update the box size

  ! subroutines in utils.f90
  if ( iproc == 0 ) then
     ! set up of xyz file.
     call convert_to_chain( mincounter, 4, scounter )
     call convert_to_chain( atp       , 2, rcounter )
     call convert_to_chain( pas       , 3, pcounter )
     fname = 'p_'//trim(scounter)//'_'//trim(rcounter)//'_'//trim(pcounter)//'_'//stage//".xyz"
     fname = trim(fname)

     ! If there is a constraint over a given atom, is written in the geometry file.
     do i = 1, NATOMS
        if      ( constr(i)== 0) then
           frzchain(i)='    '
        else if (constr(i) == 1) then
           frzchain(i)='   f'
        else if (constr(i) == 2) then
           frzchain(i)='  fy'
        else if (constr(i) == 3) then
           frzchain(i)=' fxz'
        end if
     end do

     open(unit=XYZ,file=fname,status='unknown',action='write',iostat=ierror)

     write(XYZ,*) NATOMS,  'angstroem'

     if (boundary == 'P') then
        write(XYZ,'(a,3(1x,1p,e24.17,0p))') 'periodic', (boxl(i),i=1,3)
     else if (boundary == 'S') then
        write(XYZ,'(a,3(1x,1p,e24.17,0p))') 'surface',  (boxl(i),i=1,3)
     else if (boundary == 'T') then
        write(XYZ,'(a,3(1x,1p,e24.17,0p))')'triclinic', cell(:,1)
        write(XYZ,'(3(e24.17))') cell(:,2)
        write(XYZ,'(3(e24.17))') cell(:,3)
     else
        write(XYZ,*)'free'
     end if

     do i= 1, NATOMS
        write(XYZ,'(1x,A2,3(2x,f16.8),2x,a4)')   Atom(i), x(i), y(i), z(i), frzchain(i)
     end do

     write(XYZ,*) '# simulation ',scounter,", attempt ", atp,", step ", pas,", stage ", stage
     write(XYZ,'(a,(1p,e17.10,0p))') ' # total energy (eV) : ', total_energy

     close(XYZ)

  end if

END SUBROUTINE save_intermediate
end module save_intermediate_mod
