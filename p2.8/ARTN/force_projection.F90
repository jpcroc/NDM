module force_projection_mod
  !> ART force_projection
  !!   It calculates:
  !!    the magnitude the force in the direction of some 'direction' (F_par),
  !!    the force vector perpendicular to that 'direction' (F_perp_V) and its
  !!    magnitude(F_perp), and the norm of the total force (norm_F).
contains
  subroutine force_projection_art ( F_par, F_perp_V, F_perp, norm_F, ref_F, direction )

    use defs , only : VECSIZE
    implicit none

    ! Arguments
    real(kind=8),                     intent(out) :: F_par     ! Norm of parallel force.
    real(kind=8), dimension(VECSIZE), intent(out) :: F_perp_V  ! Perpendicular force Vector.
    real(kind=8),                     intent(out) :: F_perp    ! Norm of F_perp_V.
    real(kind=8),                     intent(out) :: norm_F    ! Norm of the input force.
    real(kind=8), dimension(VECSIZE), intent(in)  :: ref_F     ! Input force.
    real(kind=8), dimension(VECSIZE), intent(in)  :: direction ! Eigen direction.

    ! Internal variables
    real(kind=8) :: F_perp2                  ! F_perp*F_perp.

    F_par  = dot_product( ref_F, direction )
    F_perp_V  = ref_F - F_par * direction

    F_perp2 = dot_product( F_perp_V , F_perp_V )
    F_perp  = sqrt( F_perp2 )
    ! This is = sqrt( dot_product( ref_F, ref_F ) )
    ! i.e, the norm of the force, but cheaper.
    norm_F = sqrt( F_par*F_par + F_perp2 )

  END SUBROUTINE force_projection_art
end module force_projection_mod
