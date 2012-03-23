! then relaxes it into a  a local minimum with or without volume optimization
! depending on the compilation flags. 
!
! The initial configuration is a "reference configuration", it will be the reference 
! configuration until a new event is accepted. 
!
! The reference configuration should contain the following information:
!
!   first line:         refcounter 1232      
!   second line:        total energy: 
!   third line:          -867.33454
!   fourth line:        8.34543 8.21345 8.67789
!   next NATOMS lines:  0  0.0437844 0.96894 0.847555
!   
!   The first number is the configuration number associated with the reference 
!   configuration. 
!
!   The second line gives the total energy so that it can be compared with a 
!   direct calculation.
!   
!   The third line indicates the full box size along the x, y and z directions.
!
!   The NATOMS lines are the atomic species and coordinates in Angstroems
!

subroutine initialize()
  use defs
  use gen_com_m
  use tab_imm_m 
  use art_in_ndm_module
  
  implicit none  
  integer :: i, ierror
  real(8) :: ran3
  real(8) :: tmp_local

  character(len=20) :: dummy,fname
  character(len=4) ::scounter
  logical :: flag

  ! Read the atomic positions
!old_art  open(unit=FREFCONFIG,file=REFCONFIG,status='old',action='read',iostat=ierror)
!old_art    read(FREFCONFIG, '(A10,i5)' ) dummy,refcounter
!old_art    read(FREFCONFIG, '(A17)') dummy
!old_art    read(FREFCONFIG,*) ref_energy
!old_art    read(FREFCONFIG,*) boxref(1), boxref(2),boxref(3)
!old_art    read(FREFCONFIG,*) (type(i),x(i),y(i),z(i),i=1,NATOMS)
!old_art    close(FREFCONFIG)

  ! Read the counter in order to continue the run where it stopped
  ! Format:
  ! Counter:     1000

  inquire(file=COUNTER, exist=flag)
  if (flag) then 
    open(unit=FCOUNTER,file=COUNTER,status='old',action='read',iostat=ierror)
    read(FCOUNTER,'(A12,I6)') dummy, mincounter 
    close(FCOUNTER)
  else
    mincounter = 1000
  endif

      do i = 1,NATOMS
!         Atom(i) = 'Si'
         Atom(i) = ty(ityp(i))
     enddo


  write(*,*) 'Mincounter :', mincounter

  ! We rescale the coordinates and define the reference coordinates
  scalaref = 1.0d0
  scala = scalaref
  
!old_art  box(:) = boxref(:)
!old_art  posref = pos


! the box, boxref from ndm into art

   do i=1,3
    boxref(i)=at(i,i) * angst
   end do
   box(:)=boxref(:)
   
   if ( NATOMS /= im ) then
    write(*,*)'initialize(): NATOMS and im:', NATOMS, im
    write(*,*)'initialize(): NATOMS must have the same value as im'
    stop   
   end if  
 
   write(*,*) 'initialize():  art boxref:', boxref(1),boxref(2), boxref(3)
 
    
!   tmp_local = dabs( SUM(at(:,:)**2) - SUM(boxref(:)**2)/angst**2 )
!debug_art   write(*,*) 'initialize(): tmp_local:', tmp_local
!   if (tmp_local  > low_limit ) then
!   !
!    write(*,*) 'initialize(): the box is not rectangular'
!     do i=1,3
!      write(*,'(3f18.7)') at(:,i)*angst
!     end do
!    stop
!   !
!   end if  
   

! the xp from the NDM go into pos,posref of the ART
 
   do i=1,NATOMS
   !
    pos(i)         = xp(1,i) * angst
    pos(NATOMS+i)  = xp(2,i) * angst
    pos(2*NATOMS+i)= xp(3,i) * angst



   !
   type(i)=ityp(i)
   end do
   posref(:) = pos(:)
   
   
    


  total_energy = ref_energy

  ! Initialise the potential if needed 
  write(*,*) 'Initial total energy : ', ref_energy

!  call init_potential()
  open(unit=CSAD,file=CSADDLE)
  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(*,*) 'Initial total energy : ', ref_energy
  write(FLOG,*) 'Initial total energy : ', ref_energy
  close(FLOG)  

  ! If this is a new event, we minimize before starting
  if (new_event) then 
     ! Converge the configuration to a local minimum
     call min_converge()
     open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
     write(*,*) 'Relaxed energy : ', total_energy
     write(FLOG,*) 'Relaxed energy : ', total_energy

     ! Assign the relaxed configuration as the reference and write it to ref-file
     posref = pos
     call write_refconfig()
 
     call convert_to_chain(mincounter,scounter)
     fname =   FINAL // scounter
     conf_initial = fname
     ! We now store the configuration into fname
     call store(fname)

     write(*,*) 'Configuration stored in file ',fname
     write(FLOG,*) 'Configuration stored in file ',fname

     mincounter = mincounter + 1
     ref_energy = total_energy
     close(flog)     
  endif


end subroutine
