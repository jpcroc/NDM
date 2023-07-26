module disp
  use triclinic_pbc_positions_mod,only:triclinic_distance
contains

  !> ART displacement
  !!    It computes the distance between two configurations and
  !!    the number of particles having moved by more than a THRESHOLD
  subroutine displacement( posa, posb, delr, npart, idmax )

    use defs
    implicit none

    !Arguments
    real(kind=8), dimension(vecsize), intent(in), target :: posa
    real(kind=8), dimension(vecsize), intent(in), target :: posb
    real(kind=8), intent(out)                            :: delr
    integer,      intent(out)                            :: npart
    integer,      intent(out), optional                  :: idmax

    !Local variables
    integer :: i
    real(kind=8),dimension(3) :: invbox
    real(kind=8), parameter :: THRESHOLD = 0.1d0  ! In Angstroems
    real(kind=8), dimension(:), pointer :: xa, ya, za, xb, yb, zb
    real(kind=8) :: delx, dely, delz, dr, dr2, delr2, dr_

    ! We first set-up pointers for the x, y, z components for posa and posb

    if (.not.allocated(atdisp))allocate(atdisp(natoms))
    if (.not.allocated(ldisp))allocate(ldisp(natoms))
    atdisp(:)=0
    ldisp=.false.
    xa => posa(1:NATOMS)
    ya => posa(NATOMS+1:2*NATOMS)
    za => posa(2*NATOMS+1:3*NATOMS)

    xb => posb(1:NATOMS)
    yb => posb(NATOMS+1:2*NATOMS)
    zb => posb(2*NATOMS+1:3*NATOMS)

    delr2 = 0.0d0
    npart = 0
    dr_   = -1.0d0
    if (present(idmax)) idmax = -1d9

    if (.not. boundary == 'T') invbox = 1.0d0/box

    do i = 1, NATOMS
       if ( boundary == 'P' ) then
          delx = ( xa(i) - xb(i) ) - box(1) * nint(( xa(i) - xb(i) )*invbox(1))
          dely = ( ya(i) - yb(i) ) - box(2) * nint(( ya(i) - yb(i) )*invbox(2))
          delz = ( za(i) - zb(i) ) - box(3) * nint(( za(i) - zb(i) )*invbox(3))

       else if ( boundary == 'S' ) then
          !be carefull with boundaries if surface: In bigdft the surface is
          !always perpendicular to "y"

          delx = ( xa(i) - xb(i) ) - box(1) * nint(( xa(i) - xb(i) )*invbox(1))
          dely = ( ya(i) - yb(i) )
          delz = ( za(i) - zb(i) ) - box(3) * nint(( za(i) - zb(i) )*invbox(3))
       else if ( boundary == 'F' ) then

          delx = ( xa(i) - xb(i) )
          dely = ( ya(i) - yb(i) )
          delz = ( za(i) - zb(i) )

       else if ( boundary == 'T' ) then
          call triclinic_distance( xb(i), xa(i), yb(i), ya(i), zb(i), za(i), delx, dely, delz )

       end if

       dr2   = delx*delx + dely*dely + delz*delz
       delr2 = delr2 + dr2
       dr    = sqrt(dr2)
       atdisp(i)=dr
       if (dr > dr_ .and. present(idmax)) then
          idmax = i
          dr_   = dr
       end if

       ! could comment this part if you are not interested in counting the moved atoms

       if ( dr > NPART_DR_THRESHOLD ) then
          npart = npart + 1
          ldisp(i)=.true.
       end if
    end do

    delr = sqrt(delr2)

  END SUBROUTINE displacement
end module disp
