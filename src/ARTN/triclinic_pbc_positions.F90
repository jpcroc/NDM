module triclinic_pbc_positions_mod
  contains
  subroutine triclinic_pbc_positions( pos_ )


  use defs, only: NATOMS, cell, invcell

  implicit none

  real(kind=8), dimension(3*NATOMS), intent(inout)		:: pos_

  real(kind=8), dimension(NATOMS)						:: trix, triy, triz, orthox, orthoy, orthoz
  integer												:: i


  ! We separate the x, y and z position in the triclinic box
  trix = pos_(1 : NATOMS)
  triy = pos_(NATOMS+1 : 2*NATOMS)
  triz = pos_(2*NATOMS+1 : 3*NATOMS)

  do i=1, NATOMS
     ! We transform the triclinic positions to orthorombic positions
     orthox(i) = invcell(1,1)*trix(i) + invcell(2,1)*triy(i) + invcell(3,1)*triz(i)
     orthoy(i) = invcell(1,2)*trix(i) + invcell(2,2)*triy(i) + invcell(3,2)*triz(i)
     orthoz(i) = invcell(1,3)*trix(i) + invcell(2,3)*triy(i) + invcell(3,3)*triz(i)

     ! Then, we apply the pbc
     orthox(i) = orthox(i) - floor(orthox(i))
     orthoy(i) = orthoy(i) - floor(orthoy(i))
     orthoz(i) = orthoz(i) - floor(orthoz(i))

     ! And then, we transform the orthorombic positions back to the triclinic ones
     trix(i) = cell(1,1)*orthox(i) + cell(2,1)*orthoy(i) + cell(3,1)*orthoz(i)
     triy(i) = cell(1,2)*orthox(i) + cell(2,2)*orthoy(i) + cell(3,2)*orthoz(i)
     triz(i) = cell(1,3)*orthox(i) + cell(2,3)*orthoy(i) + cell(3,3)*orthoz(i)
  end do

  ! We put the triclinic coordinates back in the position matrix
  pos_(1 : NATOMS) = trix
  pos_(NATOMS+1 : 2*NATOMS) = triy
  pos_(2*NATOMS+1 : 3*NATOMS) = triz


END SUBROUTINE triclinic_pbc_positions

! Routine to calculate the position of the atoms in a triclinic box
subroutine triclinic_pbc_positions_local(nlocal,local_list, pos_ )


  use defs, only: NATOMS, cell, invcell

  implicit none

  integer, intent(in)                               :: nlocal
  integer, dimension(3*NATOMS), intent(in)          :: local_list
  real(kind=8), dimension(3*NATOMS), intent(inout)  :: pos_

  real(kind=8), dimension(NATOMS)		          :: trix, triy, triz, orthox, orthoy, orthoz
  integer		                                  :: i, it


  ! We separate the x, y and z position in the triclinic box
  trix = pos_(1 : NATOMS)
  triy = pos_(NATOMS+1 : 2*NATOMS)
  triz = pos_(2*NATOMS+1 : 3*NATOMS)

  do it=1, nlocal
     i = local_list(i)
     ! We transform the triclinic positions to orthorombic positions
     orthox(i) = invcell(1,1)*trix(i) + invcell(2,1)*triy(i) + invcell(3,1)*triz(i)
     orthoy(i) = invcell(1,2)*trix(i) + invcell(2,2)*triy(i) + invcell(3,2)*triz(i)
     orthoz(i) = invcell(1,3)*trix(i) + invcell(2,3)*triy(i) + invcell(3,3)*triz(i)

     ! Then, we apply the pbc
     orthox(i) = orthox(i) - floor(orthox(i))
     orthoy(i) = orthoy(i) - floor(orthoy(i))
     orthoz(i) = orthoz(i) - floor(orthoz(i))

     ! And then, we transform the orthorombic positions back to the triclinic ones
     trix(i) = cell(1,1)*orthox(i) + cell(2,1)*orthoy(i) + cell(3,1)*orthoz(i)
     triy(i) = cell(1,2)*orthox(i) + cell(2,2)*orthoy(i) + cell(3,2)*orthoz(i)
     triz(i) = cell(1,3)*orthox(i) + cell(2,3)*orthoy(i) + cell(3,3)*orthoz(i)
  end do

  ! We put the triclinic coordinates back in the position matrix
  pos_(1 : NATOMS) = trix
  pos_(NATOMS+1 : 2*NATOMS) = triy
  pos_(2*NATOMS+1 : 3*NATOMS) = triz

END SUBROUTINE triclinic_pbc_positions_local



! Routine to calculate the distance between two points in a triclinic box
subroutine triclinic_distance( tx1, tx2, ty1, ty2, tz1, tz2, tx12, ty12, tz12 )


  use defs, only: cell, invcell

  implicit none

  real(kind=8), intent(in)		:: tx1, tx2, ty1, ty2, tz1, tz2		! Triclinic coordinates
  real(kind=8), intent(out)	:: tx12, ty12, tz12					! Triclinic distances

  real(kind=8)					:: ox1, ox2, oy1, oy2, oz1, oz2		! Orthorombic coordinates
  real(kind=8)					:: ox12, oy12, oz12                             ! Orthorombic distances


  ! We transform the triclinic box coordinates to the orthorombic box coordinates
  ox1 = invcell(1,1)*tx1 + invcell(2,1)*ty1 + invcell(3,1)*tz1
  oy1 = invcell(1,2)*tx1 + invcell(2,2)*ty1 + invcell(3,2)*tz1
  oz1 = invcell(1,3)*tx1 + invcell(2,3)*ty1 + invcell(3,3)*tz1

  ox2 = invcell(1,1)*tx2 + invcell(2,1)*ty2 + invcell(3,1)*tz2
  oy2 = invcell(1,2)*tx2 + invcell(2,2)*ty2 + invcell(3,2)*tz2
  oz2 = invcell(1,3)*tx2 + invcell(2,3)*ty2 + invcell(3,3)*tz2

  ! Then, we apply the pbc
  ox12 = ox2 - ox1 - nint( ox2 - ox1 )
  oy12 = oy2 - oy1 - nint( oy2 - oy1 )
  oz12 = oz2 - oz1 - nint( oz2 - oz1 )

  ! And then, we transform the orthorombic box positions back to the triclinic ones
  tx12 = cell(1,1)*ox12 + cell(2,1)*oy12 + cell(3,1)*oz12
  ty12 = cell(1,2)*ox12 + cell(2,2)*oy12 + cell(3,2)*oz12
  tz12 = cell(1,3)*ox12 + cell(2,3)*oy12 + cell(3,3)*oz12


END SUBROUTINE triclinic_distance

end module triclinic_pbc_positions_mod
