module opt_marginal_likelihood
use ml_in_ndm_module, ONLY : rangml, lambda_krr,isave_ml
use temporary_data_cov
use def_kernels, ONLY : length_kse,sigma_kse
use set_limits
use extrapolation
use math
use cholesky_mod
use build_kernel_matrix_mod
contains



subroutine min_marginal_likelihood
implicit none

real(kind=kind(1.d0)) :: r_coeff

 call fill_xdesc_yfunc_marginal_likelihood (rangml)
 call set_limit_for_cov(rangml, dim_train,i_final_cov,i_start_cov)
 call build_kernel_matrix(xdesc_train,dim_xdesc,dim_train)
#ifdef PARAML
       if (isave_ml==2) then
         call cholesky_threading()
        else
         call cholesky_scalapack()
       end if
#else
       call cholesky_threading()
#endif
    ! else
    ! end if
    if (allocated(y_extra))  deallocate(y_extra); allocate(y_extra(dim_valid))
    call krr_extrapolation (dim_valid,xdesc_valid, y_extra,error_valid)
    call correlation_coef(yfunc_valid, y_extra, dim_valid, r_coeff)
    if(rangml==0) write(6,*)'ML:  Marginal L ....:',  r_coeff
    if (rangml==0) then

!     do i=1,dim_valid
!       if (rangml==0) write(62,*) i, yfunc_valid(i), y_extra(i), error_valid(i)
!     end do
    end if
   !debug  write(*,*) rangml,im, dim_train, i_start_at, i_final_at,i_final_cov,i_start_cov

   if(rangml==0) write(50,'(3(G16.8))')lambda_krr,length_kse,sigma_kse

return
end subroutine min_marginal_likelihood






end module opt_marginal_likelihood




subroutine fill_xdesc_yfunc_marginal_likelihood (rangml)
use ml_in_ndm_module, only: debug
use temporary_data_cov, ONLY: dim_data, dim_train, dim_valid, dim_xdesc, &
                              xdesc,       yfunc ,      &
                              xdesc_train, yfunc_train, &
                              xdesc_valid, yfunc_valid, error_valid
implicit none
integer, intent(in) :: rangml
integer :: dim_ml_valid, dim_ml_train

if (debug) then
 if (rangml==0) write(6,'("ML: Entering in fill_xdesc_yfunc_marginal_likelihood")')
end if

dim_ml_train=int(2.d0*dim_data/3.d0)
dim_ml_valid=dim_data-dim_ml_train


if (allocated(xdesc_train)) deallocate(xdesc_train); allocate(xdesc_train(dim_xdesc, dim_ml_train))
if (allocated(xdesc_valid)) deallocate(xdesc_valid); allocate(xdesc_valid(dim_xdesc, dim_ml_valid))

if (allocated(yfunc_train)) deallocate(yfunc_train); allocate(yfunc_train(dim_ml_train))
if (allocated(yfunc_valid)) deallocate(yfunc_valid); allocate(yfunc_valid(dim_ml_valid))



xdesc_train(1:dim_xdesc, 1:dim_ml_train)=xdesc(1:dim_xdesc, 1:dim_ml_train)
yfunc_train(1:dim_ml_train) = yfunc(1:dim_ml_train)
xdesc_valid(1:dim_xdesc, 1:dim_ml_valid)=xdesc(1:dim_xdesc, dim_ml_train+1:dim_data)
yfunc_valid(1:dim_ml_valid)=yfunc(dim_ml_train+1:dim_data)


dim_train=dim_ml_train
dim_valid=dim_ml_valid

if (allocated(error_valid)) deallocate(error_valid); allocate(error_valid(dim_valid))


return

end subroutine fill_xdesc_yfunc_marginal_likelihood
