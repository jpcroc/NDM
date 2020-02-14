module extrapolation

interface krr_extrapolation
module procedure kernel_extrapolation_one
end interface krr_extrapolation


contains

subroutine kernel_extrapolation_one(dim_extra,x_to_extra,y_extra,error_extra)
use ml_in_ndm_module , ONLY: rangml, krr_error, kernel_type
use temporary_data_cov, ONLY: dim_train,matfor,dim_xdesc,kbfunc,xdesc_train
use def_kernels

implicit none
integer, intent(in) :: dim_extra
real(kind=kind(1.d0)),intent(in)   ::  x_to_extra(dim_xdesc,dim_extra)
real(kind=kind(1.d0)),intent(out)  ::  y_extra(dim_extra), error_extra(dim_extra)
! local variables ...
real(kind=kind(1.d0)) :: k_ij, temp_j(dim_train),temp_chol(dim_train)
integer :: i,j

temp_j(:)=0
error_extra(:)=0
do j=1,dim_extra
  do i=1,dim_train

         select case(kernel_type)
          case (kernel_se)
             k_ij=kernel_square_exp(x_to_extra(:,j),xdesc_train(:,i),dim_xdesc)
          case (kernel_ou)
             k_ij=kernel_ornstein_uhlenbeck (x_to_extra(:,j),xdesc_train(:,i),dim_xdesc)
          case (kernel_mc)
             k_ij=kernel_matern_class(x_to_extra(:,j),xdesc_train(:,i),dim_xdesc)
          case (kernel_so)
             k_ij=kernel_soap(x_to_extra(:,j),xdesc_train(:,i),dim_xdesc)
        end select
             temp_j(i)=k_ij
  end do

  y_extra(j)=DOT_PRODUCT(temp_j(:),kbfunc(:))

! Error evaluation in the kernel method. Implemented only in the serial case.
  if (krr_error) then

    temp_chol(:) = temp_j(:)
#ifdef PARAML
    if (rangml==0) then
#endif
    call cholesky_solve_k(matfor, temp_chol,dim_train)
#ifdef PARAML
    end if
#endif
    error_extra(j) =1.d0-DOT_PRODUCT(temp_j(:),temp_chol(:))
   end if


end do
return
end subroutine kernel_extrapolation_one

subroutine cholesky_solve_k(matfor,b,n)
use ml_in_ndm_module, ONLY  : rangml
implicit none

integer, intent(in)   :: n
real(kind=kind(1.d0)),intent(in)    :: matfor(n,n)
real(kind=kind(1.d0)),intent(inout) :: b(n)
integer :: info

call dtrtrs( 'L', 'N', 'N', n, 1, matfor, n, b, n, info )

!call dtrtrs( uplo, trans, diag, n, nrhs, a, lda, b, ldb, info )
!
!uplo
!CHARACTER*1. Must be 'U' or 'L'.
!Indicates whether A is upper or lower triangular:
!If uplo = 'U', then A is upper triangular.
!If uplo = 'L', then A is lower triangular.
!trans
!CHARACTER*1. Must be 'N' or 'T' or 'C'.
!If trans = 'N', then A*X = B is solved for X.
!If trans = 'T', then AT*X = B is solved for X.
!If trans = 'C', then AH*X = B is solved for X.
!diag
!CHARACTER*1. Must be 'N' or 'U'.
!If diag = 'N', then A is not a unit triangular matrix.
!If diag = 'U', then A is unit triangular: diagonal elements of A are assumed to be 1 and not referenced in the array a.
!n
!INTEGER. The order of A; the number of rows in B; n � 0.
!nrhs
!INTEGER. The number of right-hand sides; nrhs � 0.
!a, b
!REAL for strtrs
!DOUBLE PRECISION for dtrtrs
!COMPLEX for ctrtrs
!DOUBLE COMPLEX for ztrtrs.
!Arrays: a(lda,*), b(ldb,*).
!The array a contains the matrix A.
!The array b contains the matrix B whose columns are the right-hand sides for the systems of equations.
!The second dimension of a must be at least max(1,n), the second dimension of b at least max(1,nrhs).
!lda
!INTEGER. The leading dimension of a; lda � max(1, n).
!ldb
!INTEGER. The leading dimension of b; ldb � max(1, n).
!Output Parameters
!b
!Overwritten by the solution matrix X.
!info
!INTEGER. If info=0, the execution is successful.
!If info = -i, the i-th parameter had an illegal value.


if (info<0) then

 if (rangml==0) write(6,*) ' the i-th parameter had an illegal value.', info
 if (rangml==0) write(6,*) ' the program will stop in cholesky_solve_k'
 stop
end if

return
end subroutine cholesky_solve_k

end module extrapolation
