MODULE random_art

  ! Random number generator (from "Numerical Recipes").
  ! Returns a uniform random deviate between 0.0 and 1.0.
  ! Set idum to any negative value to initialize or  
  ! reinitialize the sequence.                      

  ! Shared variables
  save
  integer :: idum, inext, inextp
  integer :: iff = 0
  integer, dimension(55) :: ma
end module random_art

module lanczos_defs
  !use defs
  implicit none
  save

  logical :: first_time = .true., reject= .false., self_consistent= .false.
  integer :: lanczos_iter
  real(8) :: eigenvalue, old_eigenvalue
  real(8) :: lanczos_step, overlap
  real(8), dimension(10) :: eigenvals

  ! Projection direction based on lanczos computations of lowest eigenvalues
  real(8), dimension(:), allocatable :: old_projection, first_projection,old_before_sc_projection

end module lanczos_defs