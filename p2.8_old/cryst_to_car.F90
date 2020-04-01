!
!-----------------------------------------------------------------------
subroutine cryst_to_cart(nvec, vec, trmat, iflag)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  !             Version du 01 septembre 2000
  !-----------------------------------------------------------------------
  !
  !     This routine transforms the atomic positions or the k-point
  !     components from crystallographic to carthesian coordinates ( iflag=1)
  !     and viceversa ( iflag=-1 ).
  !     Output carth. coordinates are stored in the input ('vec') array.
  !
  !
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer , intent(in) :: nvec
  integer , intent(in) :: iflag
  real(double) , intent(inout) :: vec(3,nvec)
  real(double) , intent(in) :: trmat(3,3)
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: nv, kpol
  real(double), dimension(3) :: vau
  !-----------------------------------------------
  !
  !     first the dummy variables
  !
  !                          ! if iflag=1:
  !                          !    trmat = at ,  basis of the real-space latt.
  !                          !                  for atoms   or
  !                          !          = bg ,  basis of the rec.-space latt.
  !                          !                  for k-points
  !                          ! if iflag=-1: the opposite
  !
  !    here the local variables
  !
  !
  !
  !     Compute the carth. coordinates of each vectors
  !     (atomic positions or k-points components)
  !
  do nv = 1, nvec
     if (iflag==1) then
        vau = trmat(:,1)*vec(1,nv)+trmat(:,2)*vec(2,nv)+trmat(:,3)*vec(3,nv&
             )
     else
        vau = trmat(1,:)*vec(1,nv)+trmat(2,:)*vec(2,nv)+trmat(3,:)*vec(3,nv&
             )
     endif
     vec(:,nv) = vau

     !     Technique anti-bug
     !     Translation des atomes en bord de boite vers xp=0.0
     !        if (iflag.eq.-1) then
     !         do kpol =1,3
     !         if (vec(kpol,nv).ge.1.0) then
     !          vec(kpol,nv)=0.0
     !         endif
     !         enddo
     !        endif

  end do
  !
  return
end subroutine cryst_to_cart
