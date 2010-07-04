! This subroutine writes the atomic positions and others to a "refconfig" file
! which will be used a the reference point until a new events gets accepted
!
!  
! Copyright Normand Mousseau May 2001

subroutine write_refconfig()
   use defs
   implicit none
   integer :: i,ierror
   real(8) :: boxl(3)

   boxl(:) = box(:) * scala  ! Update the box size 

   open(unit=FREFCONFIG,file=REFCONFIG,status='replace',action='write',iostat=ierror)
   write(FREFCONFIG,*) 'run_id: ', mincounter
   write(FREFCONFIG,*) 'total energy : '
   write(FREFCONFIG,*) total_energy
   write(FREFCONFIG,*) boxl(1), boxl(2), boxl(3)
   do i=1, NATOMS
     write(FREFCONFIG,'(1x,i6, 3(2x,F16.8))') type(i),x(i),y(i), z(i)
   enddo
   close(FREFCONFIG)

   return
end subroutine
