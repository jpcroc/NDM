module diis_defs

  implicit none

  real(kind=8), dimension(:,:), allocatable :: previous_forces
  real(kind=8), dimension(:,:), allocatable :: previous_pos
  real(kind=8), dimension(:),   allocatable :: previous_norm
  integer :: maxter

END MODULE diis_defs
