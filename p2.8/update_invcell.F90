module update_invcell_mod
contains
  subroutine update_invcell( )


  use defs, only: cell, invcell

  implicit none

  real(kind=8) :: determinant


  ! We start by computing the adjoint
  invcell(1,1) = cell(2,2)*cell(3,3)-cell(2,3)*cell(3,2)
  invcell(2,1) = cell(2,3)*cell(3,1)-cell(2,1)*cell(3,3)
  invcell(3,1) = cell(2,1)*cell(3,2)-cell(2,2)*cell(3,1)
  invcell(1,2) = cell(1,3)*cell(3,2)-cell(1,2)*cell(3,3)
  invcell(2,2) = cell(1,1)*cell(3,3)-cell(1,3)*cell(3,1)
  invcell(3,2) = cell(1,2)*cell(3,1)-cell(1,1)*cell(3,2)
  invcell(1,3) = cell(1,2)*cell(2,3)-cell(1,3)*cell(2,2)
  invcell(2,3) = cell(1,3)*cell(2,1)-cell(1,1)*cell(2,3)
  invcell(3,3) = cell(1,1)*cell(2,2)-cell(1,2)*cell(2,1)

  ! Then, we compute the determinant
  determinant = cell(1,1)* invcell(1,1) + cell(1,2) * invcell(2,1) + cell(1,3) * invcell(3,1)

  ! Finally we get the inverse
  invcell = invcell / determinant


END SUBROUTINE update_invcell
end module update_invcell_mod
