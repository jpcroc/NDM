!  Subroutine find_saddle
!
!  This subroutine initiates the random displacement at the start
!  of the ART algorithm. 
!
!  After  random escape direction has been selected, the routine call 
!  saddle_converge which will try to bring the configuration to a saddle point.
!
!  If the convergence fails, find_saddle will restart; if it succeeds, it saves 
!  the event and returns to the main loop.

!
!  The displacement can be either LOCAL and involve a single atom
!  and its neareast neighbours or
!
!  NON-LOCAL and involve ALL the atoms.
!
!  For large cells, it is preferable to use a local initial
!  displacement to prevent the creation of many trajectories
!  at the same time in different sections of the cell.
!
!   Normand Mousseau  19 June 2001
!
module saddles
  use defs
  implicit none
  save
  
  character(len=20) :: TYPE_EVENTS
 
  real(8) :: INITSTEPSIZE 
  real(8) :: LOCAL_CUTOFF

  integer :: KTER_MIN, NVECTOR_LANCZOS
  integer :: MAXITER, MAXKTER
  integer :: MAXIPERP, MAXKPERP,MAXIPERP_ORIG
  real(8) :: INCREMENT 
  real(8) :: FTHRESHOLD, FTHRESH2
  real(8) :: EIGEN_THRESH 
  real(8) :: EXITTHRESH
  real(8) :: MAXSTEPUPHILL


  integer, dimension(:), allocatable  :: atom_displaced     ! Id of local atoms displaced
  integer                     :: natom_displaced    ! # of local atoms displaced
end module saddles

subroutine find_saddle(success)
  use random_art
  use defs
  use saddles
  use lanczos_defs
  implicit none

  logical, intent(out) :: success
  integer :: i, j, k, that, ret, npart, ierror
  real(8) :: fperp, fpar, del_r, saddle_energy
  real(8) :: ran3
  character(len=4)  :: scounter
  character(len=20) :: fname

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)

  if ( NEW_EVENT ) then 

     evalf_number = 0
     allocate(atom_displaced(natoms))       

! We first copy the reference position, half box and scaling to the working ones
     scala = scalaref
     box(:)  = boxref(:)

! And copy the reference position into the working vector
     pos = posref  !Vectorial operation

! These two subroutines modify the vector pos and generate a vector of length 1
! indicating the direction of the random displacement
     if (TYPE_EVENTS .eq. 'global') then
        call global_move()
     else if (TYPE_EVENTS .eq.'defect') then
!        write(*,*) 'there are a lot of defects, defects defects, defects'
        call first_neighbours_distance()
        call local_move()
     else if (TYPE_EVENTS .eq.'lenergy') then
!        write(*,*) 'there are a lot of defects, defects defects, defects'
        call locate_energy()
        call local_move()
     else if (TYPE_EVENTS .eq.'energy') then
        write(*,*) 'there are a lot of defects, defects defects, defects'
                
        call locate_energy()
	call selected_move()
     else if (TYPE_EVENTS .eq.'oenergy') then
        write(*,*) 'There are a lot of defects which go all in the same direction'
        write(*,*) 'Remember....: all in the same direction'
                
        call locate_energy()
	call organized_selected_move()


     else
        call local_move()	
     endif
  endif

  
  ! Now, activate per se
  call saddle_converge(ret, saddle_energy, fpar, fperp)
  
  if (saddle_energy .lt. ref_energy ) then
   write(*,*) 'The energy of the SP is lower than the Minium: REJECTED:'
   ret=60077
  end if   
  ! We compute the displacement (del_r) and number of involved particles (npart)

  call displacement(pos, posref, del_r,npart)
  
  ! We write out various information to both screen and file
  write(*,"(' ','Total energy S: ',f16.4,'  npart: ', i4,'  delr: ',&
   & f12.6,'  fpar: ',f12.6,'  fperp: ',f12.6,'  ret: ',i6,' force eval: ',&
   & i6)") saddle_energy, npart,del_r,fpar,fperp,ret,evalf_number


  write(FLOG,"(' ','Total energy S: ',f16.4,'  npart: ', i4,'  delr: ',&
   & f12.6,'  fpar: ',f12.6,'  fperp: ',f12.6,'  ret: ',i6,' force eval: ',&
   & i6)") saddle_energy, npart,del_r,fpar,fperp,ret,evalf_number

  ! If the activation did not converge, for whatever reason, we restart the 
  ! routine and do not accept the new position

  if ( ( ret < 0) .or. (ret > 30000) ) then
    success = .false.
  else
    success = .true.
    ! We need to write the configuration in a sad.... file  
    call convert_to_chain(mincounter,scounter)
    write(*,*) ' Mincounter is : ', mincounter, ' and scounter is ', scounter
    fname =   SADDLE // scounter
    conf_saddle = fname
    write(6,*)'HAUTEUR DE CE COL (COL-INITIAL) ', saddle_energy - ref_energy  
    ! We now store the configuration into fname
    call store(fname)
  
    write(*,*) 'Configuration stored in file ',fname
    write(FLOG,*) 'Configuration stored in file ',fname
  endif
 
  close(FLOG)

  deallocate(atom_displaced)       

end subroutine find_saddle

! Subroutine local_move
!
! The initial random direction is taken from a restricted space based on 
! the local bonding environment. For this, we need to know the list of neighbours
! and the cut-off of the potential. Other approaches could be used also.

subroutine local_move()
  use defs
  use random_art
  use saddles

  integer                          :: i, j, that, k, i_id, j_id
  real(8)                          :: lcutoff2    ! Cut-off for local moves, squared
  real(8), dimension(VECSIZE), target :: dr
  real(8), dimension(:), pointer   :: dx, dy, dz
  real(8)                          :: dr2, dxi, dyi, dzi, dxij, dyij, dzij
  real(8)                          :: xsum, ysum, zsum, xnorm, ynorm, znorm, norm
  real(8)                          :: ran3, inf, sup, boxl(3)

  boxl(:) = box(:) *scala

  ! We assign a few pointers 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  ! Square the cut-off
  lcutoff2 = LOCAL_CUTOFF * LOCAL_CUTOFF   ! Vectorial operation

  ! Select an atom at random in this case it must be an atom having coordination number 3 
  if (preferred_atom .lt. 0 ) then
    that = int ( NATOMS * ran3() + 1)    ! Between 1 and NATOMS
  else 
    that = preferred_atom
  endif
  write(*,*) 'That = ', that
  write(FLOG,*) 'That = ', that

  !otherwise we put to zero all dx,dy and dz
  ! dr(1:3*NATOMS)=0   
  dr(:) = 0.0d0  ! Vectorial operation
  atom_displaced(:) = 0  !Vectorial operation

  do 
    dx(that) = 0.5d0 - ran3()
    dy(that) = 0.5d0 - ran3()
    dz(that) = 0.5d0 - ran3()
   
    dr2 = dx(that)**2 + dy(that)**2 + dz(that)**2
    if (dr2 < 0.25d0 ) exit  ! Ensure that the random displacement is isotropic
  end do 
  natom_displaced = 1 
  atom_displaced(that) = 1

  i_id = type(that)
  ! Now we also displace all atoms within a cut-off distance, LOCAL_CUTOFF
  xi = x(that)
  yi = y(that)
  zi = z(that)
  do j=1, NATOMS
      j_id = type(j)
      xij = x(j) - xi - boxl(1) * nint((x(j)-xi)/boxl(1))
      yij = y(j) - yi - boxl(2) * nint((y(j)-yi)/boxl(2))
      zij = z(j) - zi - boxl(3) * nint((z(j)-zi)/boxl(3))

      xij = xij
      yij = yij
      zij = zij
      dr2 = xij*xij + yij*yij + zij*zij

      if(dr2 < lcutoff2 ) then  ! Close enough, give a random displacement
        do
          dx(j) = 0.5d0 - ran3()
          dy(j) = 0.5d0 - ran3()
          dz(j) = 0.5d0 - ran3()
   
          dr2 = dx(j)**2 + dy(j)**2 + dz(j)**2
          if (dr2 < 0.25d0 ) exit  ! Ensure that the random displacement is isotropic
        end do 
        natom_displaced = natom_displaced + 1
        atom_displaced(j) = 1
      endif
   end do     

   ! We now center only on the atoms that have been randomly selected
   xsum = sum(dx(:)) / natom_displaced
   ysum = sum(dy(:)) / natom_displaced
   zsum = sum(dz(:)) / natom_displaced

   dx = dx - xsum * atom_displaced
   dy = dy - ysum * atom_displaced
   dz = dz - zsum * atom_displaced

   ! And  normalize the total random displacement effected to the value desired
   ! This renormalizes in angstroems to a displacement INITSTEPSIZE
   norm = 1.0d0 / sqrt( dot_product(dr(:), dr(:)) )

   ! Now, we normalize dr to get the initial_direction (note that this had to be
   ! done after the transfer into box units
   initial_direction  = dr * norm

   dr(:) = dr(:) * norm *INITSTEPSIZE  ! Vectorial operation

   ! Update the position using this random displacement
   pos(:) = pos(:) +  dr(:)  ! Vectorial operation

   write(*,*) 'Number of displaced atoms initially: ',natom_displaced
end subroutine local_move


! Subroutine global_move
!
! The initial random direction is taken from the full 3N-dimensional space

subroutine global_move()
  use defs
  use random_art
  use saddles

  integer :: i
  real(8) :: norm, xnorm, ynorm, znorm
  real(8), dimension(VECSIZE), target :: dr
  real(8), dimension(:), pointer    :: dx, dy, dz
  real(8) :: ran3

  ! We assign a few pointers 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  ! All atoms are displaced
  atom_displaced = 1  ! Vectorial operation
  natom_dispaced = NATOMS
 

  ! Generate a random displacement
  do i=1, VECSIZE
    dr(i) = 0.5d0 - ran3()
  end do

  ! Keep the center of mass fixed 
  call center(dr,VECSIZE)
  
  ! And renormalize the total displacement to the value desired
  norm = 0.0d0
  do i=1, VECSIZE
     norm = norm + dr(i) * dr(i)
  end do

  !  This renormalizes in angstroems to a displacement INITSTEPSIZE
  norm = 1.0d0 / sqrt(norm)

  ! Now, we normalize dr to get the initial_direction (note that this had to be
  ! done after the transfer into box units
  initial_direction  = dr * norm

  
  ! The displacement is now in length units
  dr = dr * norm *INITSTEPSIZE  ! Vectorial operation


  ! Update the position using this random displacement
  pos = pos +  dr  ! Vectorial operation

end subroutine global_move


subroutine selected_move()

  use defs
  use random_art
  use saddles

  integer                          :: i, j, that, k, i_id, j_id
  real(8), dimension(VECSIZE), target :: dr
  real(8), dimension(:), pointer   :: dx, dy, dz
  real(8)                          :: dr2
  real(8)                          :: xsum, ysum, zsum, xnorm, ynorm, znorm, norm
  real(8)                          :: ran3, inf, sup, boxl(3)

  ! We assign a few pointers 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)


 
  dr(:) = 0.0d0          ! Vectorial operation
  atom_displaced(:) = 0  !Vectorial operation
  natom_displaced = 0
  
  do i=1,selected_atoms_max
      j=selected_atoms(i)
      atom_displaced(j)=1
      natom_displaced = natom_displaced + 1
        
            do
          dx(j) = 0.5d0 - ran3()
          dy(j) = 0.5d0 - ran3()
          dz(j) = 0.5d0 - ran3()
   
          dr2 = dx(j)**2 + dy(j)**2 + dz(j)**2
          if (dr2 < 0.25d0 ) exit  ! Ensure that the random displacement is isotropic
        end do 

    
   
  end do
  
  if (natom_displaced /= selected_atoms_max) then
   write(*,*) '   -------------SEVERE--------------------  '
   write(*,*) 'THERE ARE SOME PROBLEMS IN THE SELECTED ATOMS'
   write(*,*) 'natom_displaced ...........:', natom_displaced
   write(*,*) 'selected_atoms_max ........:', selected_atoms_max
   write(*,*) 'stop'
   stop
  end if
   ! We now center only on the atoms that have been randomly selected
   xsum = sum(dx(:)) / natom_displaced
   ysum = sum(dy(:)) / natom_displaced
   zsum = sum(dz(:)) / natom_displaced

   dx = dx - xsum * atom_displaced
   dy = dy - ysum * atom_displaced
   dz = dz - zsum * atom_displaced

   ! And  normalize the total random displacement effected to the value desired
   ! This renormalizes in angstroems to a displacement INITSTEPSIZE
   norm= sqrt( dot_product(dr(:), dr(:)) )
   
   if (norm .eq. 0.0) then
    write(*,*) 'Coane esti nefericit!'
    write(*,*) 'Try Again!'
    stop
   end if
   
   norm = 1.0d0 / norm

   ! Now, we normalize dr to get the initial_direction (note that this had to be
   ! done after the transfer into box units
   initial_direction  = dr * norm

   dr(:) = dr(:) * norm *INITSTEPSIZE  ! Vectorial operation

   ! Update the position using this random displacement
   pos(:) = pos(:) +  dr(:)  ! Vectorial operation

   write(*,*) 'Number of displaced atoms initially: ',natom_displaced
   
   
end subroutine selected_move


subroutine organized_selected_move()

  use defs
  use random_art
  use saddles

  integer                          :: i, j, that, k, i_id, j_id
  real(8), dimension(VECSIZE), target :: dr
  real(8), dimension(:), pointer   :: dx, dy, dz
  real(8)                          :: dr2
  real(8)                          :: xsum, ysum, zsum, xnorm, ynorm, znorm, norm
  real(8)                          :: ran3, inf, sup, boxl(3)
  real(8)                          :: rdir1,rdir2,rdir3

  ! We assign a few pointers 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)


 
  dr(:) = 0.0d0          ! Vectorial operation
  atom_displaced(:) = 0  !Vectorial operation
  natom_displaced = 0
  do 
   rdir1 = 0.5d0 - ran3()
   rdir2 = 0.5d0 - ran3()
   rdir3 = 0.5d0 - ran3()
   dr2=rdir1**2+rdir2**2+rdir3**3 !Ensure that the random displacement is isotropic
   if (dr2 < 0.25d0 ) exit  ! Ensure that the random displacement is isotropic
  end do

  do i=1,selected_atoms_max
      j=selected_atoms(i)
      atom_displaced(j)=1
      natom_displaced = natom_displaced + 1
      dx(j) = rdir1
      dy(j) = rdir2
      dz(j) = rdir3
  end do
  
  if (natom_displaced /= selected_atoms_max) then
   write(*,*) '   -------------SEVERE--------------------  '
   write(*,*) 'THERE ARE SOME PROBLEMS IN THE SELECTED ATOMS'
   write(*,*) 'natom_displaced ...........:', natom_displaced
   write(*,*) 'selected_atoms_max ........:', selected_atoms_max
   write(*,*) 'stop'
   stop
  end if
   ! We now center only on the atoms that have been randomly selected
   xsum = sum(dx(:)) / natom_displaced
   ysum = sum(dy(:)) / natom_displaced
   zsum = sum(dz(:)) / natom_displaced

   dx = dx - xsum * atom_displaced
   dy = dy - ysum * atom_displaced
   dz = dz - zsum * atom_displaced

   ! And  normalize the total random displacement effected to the value desired
   ! This renormalizes in angstroems to a displacement INITSTEPSIZE
   norm= sqrt( dot_product(dr(:), dr(:)) )
   
   if (norm .eq. 0.0) then
    write(*,*) 'Coane esti nefericit! (Trust me deserve romanian translation)'
    write(*,*) 'Try Again!'
    stop
   end if
   
   norm = 1.0d0 / norm

   ! Now, we normalize dr to get the initial_direction (note that this had to be
   ! done after the transfer into box units
   initial_direction  = dr * norm

   dr(:) = dr(:) * norm *INITSTEPSIZE  ! Vectorial operation

   ! Update the position using this random displacement
   pos(:) = pos(:) +  dr(:)  ! Vectorial operation

   write(*,*) 'Number of displaced atoms initially: ',natom_displaced
   
   
end subroutine organized_selected_move
