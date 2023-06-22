module generate_local_region_mod
use defs!, only: NATOMS, cell, invcell,boundary
use triclinic_pbc_positions_mod,only:triclinic_distance
 

integer                            :: nat_inner, nat_outer, nat_local
  integer, dimension(:), allocatable :: inner_list, outer_list
  real(8)                            :: r2_inner, r2_outer
  real(8), dimension(:), allocatable :: mask_inner_region

contains
  
  subroutine generate_local_region(centralatom,posART)

  implicit none


  integer, intent(in)                              :: centralatom
  real(kind=8), dimension(3*NATOMS), intent(in)	   :: posART

  integer :: i,j,k, l, this_label
  integer :: ierror
  integer :: x_central_cell, y_central_cell, z_central_cell
  integer :: xmin, xmax, ymin, ymax, zmin, zmax
  integer :: xlabel, ylabel,zlabel


  integer   :: this_atom, that_atom
  real(8) :: xi,yi,zi,xij,yij,zij, dr2
  real(8), dimension(3)       :: dr, invbox

  if (.not. boundary == 'T') invbox = 1.0d0/box

  this_atom   = centralatom
  xi = posART(this_atom)
  yi = posART(this_atom + NATOMS)
  zi = posART(this_atom + 2*NATOMS)

  ! Initialise the variable
  nat_inner = 0
  nat_outer = 0
  inner_list = 0
  outer_list = 0
  mask_inner_region = 0.0d0

  do i = 1, NATOMS
     if ( boundary == 'P' ) then
        xij = posART(i)          - xi - box(1) * nint((posART(i)         -xi)*invbox(1))
        yij = posART(i+natoms)   - yi - box(2) * nint((posART(i+natoms)  -yi)*invbox(2))
        zij = posART(i+2*natoms) - zi - box(3) * nint((posART(i+2*natoms)-zi)*invbox(3))

     else if ( boundary == 'S' ) then
        ! Be carefull with boundaries if surface
        xij = posART(i)          - xi - box(1) * nint((posART(i)         -xi)*invbox(1))
        zij = posART(i+2*natoms) - zi - box(3) * nint((posART(i+2*natoms)-zi)*invbox(3))
        yij = posART(i+natoms) - yi

     else if ( boundary == 'T' ) then
        call triclinic_distance( xi, posART(i), yi, posART(i+natoms), zi, posART(i+2*natoms), xij, yij, zij )
     end if
     ! Compute norm 
     dr2 = xij*xij+yij*yij+zij*zij

     if (dr2 .le. r2_inner ) then
        nat_inner = nat_inner + 1
        inner_list(nat_inner) = i
        mask_inner_region(i)          = 1.0d0 
        mask_inner_region(i+NATOMS)   = 1.0d0 
        mask_inner_region(i+2*NATOMS) = 1.0d0 
     else if (dr2 .le. r2_outer) then
        nat_outer = nat_outer + 1
        outer_list(nat_outer) =i
     endif
  end do

  ! Get the total number of atoms in the region
  nat_local = nat_inner + nat_outer
  write(*,*) 'nat_inner :', nat_inner, ' nat_outer: ',nat_outer, 'nat_local: ', nat_local
  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(FLOG,*) 'Size of local region:  nat_inner :', nat_inner, ' nat_outer: ',nat_outer, 'nat_local: ', nat_local
  close(FLOG)
end subroutine generate_local_region

subroutine initial_local_region()
  use defs, only: NATOMS, cell, invcell
  
  implicit none
  integer :: i,j,k, ierr
  integer :: max_incell, cells_width, ncells_region
  real(8) :: outer_diameter
  real(8), dimension(3)                      :: mask_box    ! Controls the application along 3 directions



  ! Compute the mask for applying local box. If the diameter is larger than the box size in x,y or z, we do not apply the local
  ! conditions along this direction
  outer_diameter = dsqrt(r2_outer)*2
  if (outer_diameter.ge.box(1)) then
     mask_box(1) = 0.0d0
     write(*,*) 'Box along x is smaller than local box, use real box'
  else
     mask_box(1) = 1.0d0
  endif

  if (outer_diameter.ge.box(2)) then
     mask_box(2) = 0.0d0
     write(*,*) 'Box along y is smaller than local box, use real box'
  else
     mask_box(2) = 1.0d0
  endif

  if (outer_diameter.ge.box(3)) then
     mask_box(3) = 0.0d0
     write(*,*) 'Box along z is smaller than local box, use real box'
  else
     mask_box(3) = 1.0d0
  endif

end subroutine initial_local_region

end module generate_local_region_mod
