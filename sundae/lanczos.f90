!Subroutines related to lanczos calculations ....
! - lanczos 
! - center 


subroutine lanczos(N,maxvec,q1s2,new_projection,projection) !!!!!!!!!!!!! SERVE POTENZIALE PER LA LACUNA!!!!!!!!!!!!!!!!!!!!
 USE T_kind_param_m, ONLY:  double 
 use gen_com_m, ONLY : erg2ev,ev2erg
 use sundae_module, ONLY: it_trajectory,nl_iter
 use lanczos_defs
  !use random_art
  !use art_in_ndm_module
  implicit none
  integer :: maxvec
  integer :: N
  logical ::  new_projection
  logical :: lanczos_failed=.false.
  real(8), dimension( 2 * maxvec -1 ) :: scratcha
  real(8), dimension(maxvec) :: diag
  real(8), dimension(maxvec-1) :: offdiag
  real(8), dimension(maxvec, maxvec) :: vector
  real(8), dimension(3*N, maxvec), target :: lanc
  real(8):: sum_forcenew, sum_force
  ! Vectors used to build the matrix for Lanzcos algorithm 
  real(8), dimension(:), pointer :: z0, z1, z2
  ! Projection direction based on lanczos computations of lowest eigenvalues

  integer :: i,k, i_err, ivec, nl_failed, it_art,evalf_number
  real(8) :: a1,a0,b2,b1,increment!, eigenvalue
  real(8) :: excited_energy,c1,norm
  real(8) :: xsum, ysum, zsum, sum2, invsum
  real(8), dimension(3*N) :: pos, newpos,newforce,ref_force
  real(8), dimension(3*N) :: newforcep1,newforcem1,projection
  real (double), dimension (1:3,1:N) :: q1s2  ! vector of position
  real(double), parameter:: cmTOang=1.0d8
  real(double)  :: ran3

!  do i =1,N
!  pos(i)=qx(i-1)*cmTOang
!  pos(i+N)=qy(i-1)*cmTOang
!  pos(i+2*N)=qz(i-1)*cmTOang
!  enddo

  pos(1:N)       = q1s2(1,1:N)*cmTOang
  pos(N+1:2*N)   = q1s2(2,1:N)*cmTOang
  pos(2*N+1:3*N) = q1s2(3,1:N)*cmTOang



  !lanczos_step =0.001 ! in Angstroems  mettre en parametre d'entree.
  !boxl(:) = box(:) * scala
  increment = lanczos_step  ! Increment, convert in box units
  overlap=0.0
  evalf_number=0
  nl_iter=0
  nl_failed=0

  if(.not. new_projection ) then
    old_before_sc_projection = projection*cmTOang*cmTOang
    projection=old_before_sc_projection  ! Vectorial operation
    eigenvalue=eigenvalue/(ev2erg*cmTOang*cmTOang)
  end if

  ! We now take the current position as the reference point and will make 
  ! a displacement in a random direction or using the previous direction as
  ! the starting point.
  it_art=it_trajectory
  !debugC write(*,*) 'lanczos',it_art,it_trajectory,mcmoves
  34 continue
   !d1 call calcforce_lanc(N,pos,ref_force,total_energy,it_art)
   !d1 evalf_number = evalf_number + 1
  z0 => lanc(:,1)

  if(((.not. new_projection ).or. self_consistent).and.(.not. lanczos_failed)) then
    z0 = projection             ! Vectorial operation
    old_projection = projection ! Vectorial operation
    old_eigenvalue = eigenvalue 
  else
    do i=1, 3*N
      z0(i) = 0.5d0 - ran3()
    end do

    z1 => lanc(1,:)

    xsum = 0.0d0
    ysum = 0.0d0
    zsum = 0.0d0
    do i=1, N
      xsum = xsum + z0(i)
      ysum = ysum + z0(i+N)
      zsum = zsum + z0(i + 2*N)
    end do
    xsum = xsum / real(N)
    ysum = ysum / real(N)
    zsum = zsum / real(N)
    do i=1, N
      z0(i)  = z0(i)  - xsum
      z0(i+N )   = z0(i+N)   - ysum
      z0(i+2*N)  = z0(i+2*N) - zsum 
    end do

    if(first_time) then
      old_projection = z0   ! Vectorial operation
      first_time = .false.
    else
      old_projection = projection
    endif 
  endif
  ! We normalize the displacement to 1 total
  sum2 = 0.0d0
  do i=1, 3*N
    sum2 = sum2 + z0(i) * z0(i)
  end do
  invsum = 1.0/sqrt(sum2)
  z0 = z0 * invsum
    newpos = pos + z0 * increment
    call calcforce_lanc(N,newpos,newforcep1,excited_energy,it_art)
    evalf_number = evalf_number + 1

    !d4 newpos = pos +2.d0* z0 * increment
    !d4 call calcforce_lanc(N,newpos,newforcep2,excited_energy,it_art)
    !d4 evalf_number = evalf_number + 1

    newpos = pos - z0 * increment
    call calcforce_lanc(N,newpos,newforcem1,excited_energy,it_art)
    evalf_number = evalf_number + 1

    !d4 newpos = pos - 2.d0*z0 * increment
    !d4 call calcforce_lanc(N,newpos,newforcem2,excited_energy,it_art)
    !d4 evalf_number = evalf_number + 1
     
     !d4 newforce(:)=(newforcem2(:)-8.d0*newforcem1(:)+8.d0*newforcep1(:)-newforcep2(:))/12.d0
     newforce(:)=(newforcep1(:) - newforcem1(:))/2.0d0 
     !d2 newforce(:)= newforcep1(:)-ref_force(:)
    ! We extract lanczos(1)

!write(*,*) '44'!,maxval(newforce(:)), maxval(newforce1(:)), minval(newforce2(:))

  ! We get a0
  a0 = 0.0d0
  do i=1, 3*N
    a0 = a0 + z0(i) * newforce(i)
  end do
  diag(1) = a0
  z1 => lanc(:,2)
  z1 = newforce - a0 * z0    ! Vectorial operation
!write (*,*) z1

!stop
  b1=0.0d0
  do i=1, 3*N
    b1 = b1 + z1(i) * z1(i)
  end do
  offdiag(1) = sqrt(b1)

  invsum = 1.0d0 / sqrt ( b1 )
  z1 = z1 * invsum           ! Vectorial operation

  ! We can now repeat this game for the next vectors
  do ivec = 2, maxvec-1
    z1 => lanc(:,ivec)
    newpos = pos + z1 * increment
    call calcforce_lanc(N,newpos,newforcep1,excited_energy,it_art)
    evalf_number = evalf_number + 1

    !d4 newpos = pos + 2.d0*z1 * increment
    !d4 call calcforce_lanc(N,newpos,newforcep2,excited_energy,it_art)
    !d4 evalf_number = evalf_number + 1

     newpos = pos - z1 * increment
     call calcforce_lanc(N,newpos,newforcem1,excited_energy,it_art)
     evalf_number = evalf_number + 1

    !newpos = pos - 2.d0*z1 * increment
    !call calcforce_lanc(N,newpos,newforcem2,excited_energy,it_art)
    !evalf_number = evalf_number + 1


     !d4 newforce(:)=(newforcem2(:)-8.d0*newforcem1(:)+8.d0*newforcep1(:)-newforcep2(:))/12.d0
     newforce(:)=(newforcep1(:) - newforcem1(:))/2.0d0  
     !d1 newforce(:)=newforcep1(:)-ref_force(:)
 
    a1 = 0.0d0
    do i=1,3*N
      a1 = a1 + z1(i) * newforce(i)
    end do
    diag(ivec) = a1

    b1 = offdiag(ivec-1)
    z0 => lanc(:,ivec-1)
    z2 => lanc(:,ivec+1)
    z2 = newforce - a1*z1 -b1*z0

    b2=0.0d0
    do i=1, 3*N
      b2 = b2 + z2(i)*z2(i)
    end do

    offdiag(ivec) = sqrt(b2)
    
    invsum = 1.0/sqrt(b2)
    z2 = z2 * invsum
  end do
!write(*,*) 'cinque', newpos(1)
  ! We now consider the last line of our matrix
  ivec = maxvec
  z1 => lanc(:,maxvec)
    newpos = pos + z1 * increment
    call calcforce_lanc(N,newpos,newforcep1,excited_energy,it_art)
    evalf_number = evalf_number + 1

    !d4 newpos = pos + 2.d0*z1 * increment
    !d4 call calcforce_lanc(N,newpos,newforcep2,excited_energy,it_art)
    !d4 evalf_number = evalf_number + 1

     newpos = pos - z1 * increment
     call calcforce_lanc(N,newpos,newforcem1,excited_energy,it_art)
     evalf_number = evalf_number + 1

    !d4 newpos = pos - 2.d0*z1 * increment
    !d4 call calcforce_lanc(N,newpos,newforcem2,excited_energy,it_art)
    !d4 evalf_number = evalf_number + 1


      !d4 newforce(:)=(newforcem2(:)-8.d0*newforcem1(:)+8.d0*newforcep1(:)-newforcep2(:))/12.d0
      newforce(:)=(newforcep1(:) - newforcem1(:))/2.0d0  
      !d1 newforce(:)=newforcep1(:)-ref_force(:)

    sum_force = 0.0 
    sum_forcenew = 0.0

    do i = 1 , 3*N
       sum_force = sum_force + ref_force(i)
       sum_forcenew =  sum_forcenew + newforce(i) 
    end do
!   write(*,*) 'the sum of the forces before the move', sum_force,sum_forcenew

!   write(*,*) 'the sum of the forces After the move', sum_forcenew
! newforce = newforce - ref_force
  sum_forcenew = 0.0
  a1 = 0.0d0
  do i=1, 3*N
    a1 = a1 + z1(i) * newforce(i)
    sum_forcenew =  sum_forcenew + newforce(i)
  end do
  !write(*,*) 'the difference between the forces ' , sum_forcenew
!  write(*,*)
  diag(maxvec) = a1

  ! We now have everything we need in order to diagonalise and find the
  ! eigenvectors.

  diag = -1.0d0 * diag
  offdiag = -1.0d0 * offdiag

  ! We now need the routines from Lapack. We define a few values
  i_err = 0

  ! We call the routine for diagonalizing a tridiagonal  matrix
  call dstev('V',maxvec,diag,offdiag,vector,maxvec,scratcha,i_err)

  ! We now reconstruct the eigenvectors in the real space
  ! Of course, we need only the first 5*N elements of vec

  projection = 0.0d0    ! Vectorial operation
  do k=1, maxvec
    z1 => lanc(:,k)
    a1 = vector(k,1)
    projection = projection + a1 * z1   ! Vectorial operation
  end do 
  c1=0.0d0
  do i=1, 3*N
     c1 = c1 + projection(i) * projection(i)
  end do

     norm = 1.0/sqrt(c1)
     projection = projection * norm 

  ! The following lines are probably not needed.
  !newpos = pos + projection * increment   ! Vectorial operation
  !call calcforce_lanc(N,type,newpos,boxl,newforce,excited_energy)
  !evalf_number = evalf_number + 1
  !newforce = newforce - ref_force

  eigenvalue=diag(1)/increment
  do i=1, 4
    eigenvals(i) = diag(i) / increment
  end do
 !write(*,*) 'lambda', eigenvalue

  a1=0.0d0
  b1=0.0d0
  do i=1, 3*N
    a1 = a1 + old_projection(i) * projection(i)
    b1 = b1 + projection(i) * projection(i)
    c1 = c1 + old_projection(i) * old_projection(i)

  end do
!ooo  
!ooo
!ooo!   The condition on the scalar product , we reject the point where we loose the eigen value and try to reduce the step size
!ooo   if(abs(a1)<=0.2) reject = .true. 
!ooo
overlap=a1 
!d if(a1<0.0d0) then
!d    projection = -1.0d0 * projection
!d    overlap=-overlap
!d end if    
!ooo
!ooo  
!d  call center(projection,3*N)
  
! write(*,*) 'eigenv', eigenvalue-old_eigenvalue, a1

!write(*,'("lanczos nliter",i5,f5.2,2f12.3)') nl_iter,overlap, eigenvalue,old_eigenvalue-eigenvalue
! 10-3 eV/A² is the same thing as 6.15 erg/cm². For the reasons we  will take 5.  
if ( dabs((old_eigenvalue-eigenvalue)) .gt. 1.d-3) then  ! mettre lanczos_threshold en parametre d'entree
   self_consistent=.true.
   lanczos_failed=.false.
   nl_iter=nl_iter+1
   nl_failed=nl_failed+1
   if (nl_failed.gt.30) then
      nl_failed=0
      !lanczos_failed=.true.
      write(*,*) 'WARNING: LANCZOS FAILED ... convergence not reached'
   go to 35
   end if
   go to 34
 end if
35 continue 
 self_consistent=.false.
 
  a1=0.0d0
  do i=1, 3*N
    a1 = a1 + old_before_sc_projection(i) * projection(i)
  end do

  eigenvalue=eigenvalue*ev2erg*cmTOang*cmTOang

  eigenvals = eigenvals*ev2erg*cmTOang*cmTOang

  projection=projection/(cmTOang*cmTOang)
  overlap=a1 
  lanczos_iter=nl_iter*maxvec

!   The condition on the scalar product , we reject the point where we loose the eigen value and try to reduce the step size
!debug    if(abs(a1)<=0.2) reject = .true. 

 overlap=a1

!debug if(a1<0.0d0) then
!debug    projection = -1.0d0 * projection
    !overlap=-overlap
!debug end if    

 call center(projection,3*N)

!write(*,'("lanczos nliter",i5,f5.2,2f12.3)') nl_iter,overlap, eigenvalue,old_eigenvalue
!write(*,*)'newproj',projection(1) 

end subroutine lanczos



subroutine center(vector,VECSIZE)
  integer, intent(IN) :: VECSIZE
  real(8), dimension(VECSIZE),intent(inout) :: vector

  integer :: i, natoms
  real(8), dimension(VECSIZE/3) :: x, y, z     ! Pointers for coordinates
  real(8) :: xtotal, ytotal, ztotal

  natoms = VECSIZE / 3

  ! We first set-up pointers for the x, y, z components 
  x(1:natoms)= vector(1:natoms)
  y(1:natoms)= vector(natoms+1:2*natoms)
  z(1:natoms)= vector(2*natoms+1:3*natoms)

  xtotal = 0.0d0
  ytotal = 0.0d0
  ztotal = 0.0d0

  do i = 1, natoms
    xtotal = xtotal + x(i)
    ytotal = ytotal + y(i)
    ztotal = ztotal + z(i)
  enddo 

  xtotal   = xtotal / natoms
  ytotal   = ytotal / natoms
  ztotal   = ztotal / natoms

  do i = 1, natoms
    x(i)   = x(i) - xtotal
    y(i)   = y(i) - ytotal
    z(i)   = z(i) - ztotal
  end do
end subroutine

