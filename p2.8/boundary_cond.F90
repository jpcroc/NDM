module boundary_cond_mod
  use triclinic_pbc_positions_mod,only:triclinic_distance
contains

  !> ART boundary_cond
subroutine boundary_cond ( posr, posa, posb )

  use defs, only : natoms, vecsize, box, boundary
  implicit none

  !Arguments
  real(kind=8), dimension(vecsize), intent(out), target :: posr
  real(kind=8), dimension(vecsize), intent(in), target  :: posa
  real(kind=8), dimension(vecsize), intent(in), target  :: posb

  !Local variables
  integer :: i
  real(kind=8),dimension(3) :: invbox
  real(kind=8), dimension(:), pointer :: xr, yr, zr, xa, ya, za, xb, yb, zb

  ! We first set-up pointers for the x, y, z components for posr, posa and posb

  posr  = 0.0d0
  xr => posr(1:NATOMS)
  yr => posr(NATOMS+1:2*NATOMS)
  zr => posr(2*NATOMS+1:3*NATOMS)

  xa => posa(1:NATOMS)
  ya => posa(NATOMS+1:2*NATOMS)
  za => posa(2*NATOMS+1:3*NATOMS)

  xb => posb(1:NATOMS)
  yb => posb(NATOMS+1:2*NATOMS)
  zb => posb(2*NATOMS+1:3*NATOMS)

  if (.not. boundary == 'T') invbox = 1.0d0/box

  do i = 1, natoms
     if ( boundary == 'P' ) then
        xr(i) = ( xa(i) - xb(i) ) - box(1) * nint(( xa(i) - xb(i) )*invbox(1))
        yr(i) = ( ya(i) - yb(i) ) - box(2) * nint(( ya(i) - yb(i) )*invbox(2))
        zr(i) = ( za(i) - zb(i) ) - box(3) * nint(( za(i) - zb(i) )*invbox(3))

     else if ( boundary == 'S' ) then
        !be carefull with boundaries if surface: In bigdftif the surface is
        !always perpendicular to "y"

        xr(i) = ( xa(i) - xb(i) ) - box(1) * nint(( xa(i) - xb(i) )*invbox(1))
        yr(i) = ( ya(i) - yb(i) )
        zr(i) = ( za(i) - zb(i) ) - box(3) * nint(( za(i) - zb(i) )*invbox(3))

     else if ( boundary == 'F' ) then

        xr(i) = ( xa(i) - xb(i) )
        yr(i) = ( ya(i) - yb(i) )
        zr(i) = ( za(i) - zb(i) )

     else if ( boundary == 'T' ) then
        call triclinic_distance( xb(i), xa(i), yb(i), ya(i), zb(i), za(i), xr(i), yr(i), zr(i) )

     end if
  enddo

END SUBROUTINE boundary_cond
end module boundary_cond_mod
