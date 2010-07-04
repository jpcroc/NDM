module lanczos_defs
  use defs
  implicit none
  save

  logical :: first_time = .true., reject= .false.
  real(8) :: eigenvalue, old_eigenvalue
  real(8) :: lanczos_step
  real(8), dimension(10) :: eigenvals

  ! Projection direction based on lanczos computations of lowest eigenvalues
  real(8), dimension(:), allocatable :: old_projection, projection, first_projection
end module lanczos_defs
  
subroutine lanczos(maxvec,new_projection)
  use defs
  use lanczos_defs
  use random_art
  implicit none

  integer, intent(in) :: maxvec
  logical, intent(in) ::  new_projection
  integer, dimension(maxvec) :: iscratch
  real(8), dimension( 2 * maxvec -1 ) :: scratcha
  real(8), dimension(maxvec) :: diag
  real(8), dimension(maxvec-1) :: offdiag
  real(8), dimension(maxvec, maxvec) :: vector
  real(8), dimension(VECSIZE, maxvec), TARGET :: lanc, proj
  
  
  real(8):: sum_forcenew, sum_force
  ! Vectors used to build the matrix for Lanzcos algorithm 
  real(8), dimension(:), pointer :: z0, z1, z2
 
  integer :: i,j,k, i_err, scratcha_size,ivec

  real(8) :: a1,a0,b2,b1,increment
  real(8) :: boxl(3), excited_energy,c1,norm
  real(8) :: xsum, ysum, zsum, sum2, invsum
  real(8), dimension(VECSIZE) :: newpos,newforce,ref_force,newpos_,newforce_
  real(8) :: ran3

 
  boxl(:) = box(:) * scala
  increment = lanczos_step  ! Increment, convert in box units

  ! We now take the current position as the reference point and will make 
  ! a displacement in a random direction or using the previous direction as
  ! the starting point.
34 continue
  call calcforce(NATOMS,type,pos,boxl,ref_force,total_energy)
  evalf_number = evalf_number + 1
  z0 => lanc(:,1)
  
  if(.not. new_projection ) then
     z0 = projection             ! Vectorial operation
     old_projection = projection ! Vectorial operation
     old_eigenvalue = eigenvalue 
  else
     
     xsum = 0.0d0
     ysum = 0.0d0
     zsum = 0.0d0
     do i=1, VECSIZE
        z0(i) = 0.5d0 - ran3()
        if (i >=1 .and. i <= NATOMS) then
           xsum = xsum + z0(i)
        else if (i > NATOMS .and. i <= 2*NATOMS) then
           ysum = ysum + z0(i)
        else
           zsum = zsum + z0(i)
        end if   
        xsum = xsum / NATOMS
        ysum = ysum / NATOMS
        zsum = zsum / NATOMS
     end do
     
     sum2 = 0.0d0
     do i=1, NATOMS
        z0(i)           = z0(i)          - xsum
        z0(i+NATOMS )   = z0(i+NATOMS)   - ysum
        z0(i+2*NATOMS)  = z0(i+2*NATOMS) - zsum 
        sum2 = sum2 + z0(i)**2 + z0(i+NATOMS)**2 + z0(i+2*NATOMS)**2
     end do
     invsum = 1.0/sqrt(sum2)
     z0 = z0 * invsum
     
     if(first_time) then
        old_projection = z0   ! Vectorial operation
        first_time = .false.
     else
        old_projection = projection
     endif
  endif
  ! We normalize the displacement to 1 total

  newpos  = pos + z0 * increment   ! Vectorial operation
  !newpos_ = pos - z0 * increment   ! Vectorial operation
  call calcforce(NATOMS,type,newpos,boxl,newforce,excited_energy)
  evalf_number = evalf_number + 1
  !call calcforce(NATOMS,type,newpos_,boxl,newforce_,excited_energy)
  !evalf_number = evalf_number + 1
  
  !newforce = (newforce - newforce_)/2 !/(2. * increment)     ! KIMIYA VERSION

  ! We extract lanczos(1)
   newforce = (newforce - ref_force) !/increment     ! COSMIN VERSION

  ! We get a0

  a0 = 0.0d0
  do i=1, VECSIZE
     a0 = a0 + z0(i) * newforce(i)   ! a_1 = q_1' * r_1 = q_1' * (A*q_1 - b_0*q_0(=0)) 
  end do
  diag(1) = a0
  
  z1 => lanc(:,2)
  z1 =  newforce - a0 * z0    ! r_1 = r_1 - a_1 * q_1  in the lecture notes
  
  b1 = DOT_PRODUCT(z1(:),z1(:))
  offdiag(1) = sqrt(b1)
  
  invsum = 1.0d0 / sqrt ( b1 )
  z1 = z1 * invsum           ! Vectorial operation
  
  ! We can now repeat this game for the next vectors
  do ivec = 2, maxvec-1
     z1 => lanc(:,ivec)
     
     newpos = pos + z1 * increment
     !newpos_ = pos - z1 * increment   ! Vectorial operation
     call calcforce(NATOMS,type,newpos,boxl,newforce,excited_energy)
     evalf_number = evalf_number + 1
     !call calcforce(NATOMS,type,newpos_,boxl,newforce_,excited_energy)
     !evalf_number = evalf_number + 1
     newforce = (newforce - ref_force) !/increment
     
     !newforce = (newforce - newforce_)/2 !/(2. * increment)     ! KIMIYA VERSION
     
     a1 = DOT_PRODUCT(z1(:),newforce(:))
     diag(ivec) = a1
     
     b1 = offdiag(ivec-1)
     z0 => lanc(:,ivec-1)
     z2 => lanc(:,ivec+1)
     z2 = newforce - a1*z1 -b1*z0  !! so far correct !!
     
     b2 = DOT_PRODUCT(z2(:),z2(:))
     offdiag(ivec) = sqrt(b2)
     
     invsum = 1.0/sqrt(b2)
     z2 = z2 * invsum
  end do
  
  ! We now consider the last line of our matrix
  ivec = maxvec
  z1 => lanc(:,maxvec)
  newpos = pos + z1 * increment    ! Vectorial operation
  !newpos_ = pos - z1 * increment   ! Vectorial operation
  call calcforce(NATOMS,type,newpos,boxl,newforce,excited_energy)
  evalf_number = evalf_number + 1
  !call calcforce(NATOMS,type,newpos_,boxl,newforce_,excited_energy)
  !evalf_number = evalf_number + 1
  newforce = (newforce - ref_force)!/increment
      
  !newforce = (newforce - newforce_) !/(2. * increment)     ! KIMIYA VERSION

  a1 = DOT_PRODUCT(z1(:),newforce(:))
  !  write(*,*) 'the difference between the forces ' , sum_forcenew
  !  write(*,*)
  diag(maxvec) = a1
  
  ! We now have everything we need in order to diagonalise and find the
  ! eigenvectors.
  
  diag = -1.0d0 * diag
  offdiag = -1.0d0 * offdiag

  ! We now need the routines from Lapack. We define a few values
  i_err = 0

  ! We call the routine for diagonalizing a tridiagonal  matrix
  call dstev('V', maxvec, diag, offdiag, vector, maxvec, scratcha, i_err)

  ! We now reconstruct the eigenvectors in the real space
  ! Of course, we need only the first maxvec elements of vec

  projection = 0.0d0    ! Vectorial operation
  do k=1, maxvec
    z1 => lanc(:,k)
    a1 = vector(k,1)
    projection = projection + a1 * z1   ! Vectorial operation
  end do 

  c1 = DOT_PRODUCT(projection(:),projection(:))
  
  norm = 1.0/sqrt(c1)
  projection = projection*norm 
  
  ! The following lines are probably not needed.
!  newpos = pos + projection * increment   ! Vectorial operation
!  call calcforce(NATOMS,type,newpos,boxl,newforce,excited_energy)
!  evalf_number = evalf_number + 1
!  newforce = newforce - ref_force

  eigenvalue=diag(1) /(increment)
  do i=1, 4
    eigenvals(i) = diag(i) / (increment)
  end do

  a1=0.0d0
  b1=0.0d0
  do i=1, VECSIZE
     a1 = a1 + old_projection(i) * projection(i)
     b1 = b1 + projection(i) * projection(i)
     c1 = c1 + old_projection(i) * old_projection(i)
  end do

!!$  z1 => lanc(:,1)
!!$  do i = 1,maxvec
!!$     z2 => lanc(:,i)
!!$     write(*,*) 'DOT product of vectors = ', DOT_PRODUCT(z1(:),z2(:))
!!$  end do

  
  
  !   The condition on the scalar product ...
  !   we reject the point where we lose the eigenvalue and try to reduce the step size
  if(abs(a1)<=0.2) reject = .true. 
  
!ooo  if(a1<0.0d0) projection = -1.0d0 * projection 
  
  call center(projection,VECSIZE)
end subroutine lanczos
