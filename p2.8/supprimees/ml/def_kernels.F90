
module def_kernels
USE T_kind_param_m, ONLY:  double
implicit none


!real(double), dimension(:,:), allocatable :: xdesc ! first  index is the dimension of the descriptor

real(double)            :: length_kse
real(double)            :: sigma_kse
real(double), parameter :: l_kou=1.d0


integer, parameter ::  kernel_se=1, &
                       kernel_ou=2, &
                       kernel_mc=3, &
                       kernel_so=4
integer :: kappa_soap

contains



real (double) function kernel_square_exp(a_local, b_local,n_dimension_desc)
implicit none
integer :: n_dimension_desc
real(double), dimension(n_dimension_desc) :: a_local, b_local

  kernel_square_exp=sigma_kse*dexp(-real(SUM((a_local(:)-b_local(:))*(a_local(:)-b_local(:))))/(2.d0*length_kse))

return
end function kernel_square_exp


real (double) function d_length_kernel_square_exp(a_local, b_local,n_dimension_desc)
implicit none
integer :: n_dimension_desc
real(double), dimension(n_dimension_desc) :: a_local, b_local
real(double) :: temp

temp=SUM((a_local(:)-b_local(:))*(a_local(:)-b_local(:)))/(2.d0*length_kse)

  d_length_kernel_square_exp=sigma_kse*dexp(-temp)-temp/length_kse

return
end function d_length_kernel_square_exp



real(double) function kernel_ornstein_uhlenbeck(a_local, b_local,n_dimension_desc)
implicit none
integer :: n_dimension_desc
real(double), dimension(n_dimension_desc) :: a_local, b_local

  kernel_ornstein_uhlenbeck=dexp(-dsqrt((SUM((a_local(:)-b_local(:))*(a_local(:)-b_local(:)))))/dble(l_kou))

return
end function kernel_ornstein_uhlenbeck


real(double) function kernel_matern_class(a_local, b_local,n_dimension_desc)
implicit none
integer :: n_dimension_desc
real(double), dimension(n_dimension_desc) :: a_local, b_local

  kernel_matern_class=dexp(-dsqrt((SUM((a_local(:)-b_local(:))*(a_local(:)-b_local(:)))))/dble(l_kou))

return
end function kernel_matern_class


real(double) function kernel_soap(a_local, b_local,n_dimension_desc)
implicit none
integer :: n_dimension_desc
real(double), dimension(n_dimension_desc) :: a_local, b_local

  kernel_soap=dot_product(a_local(:),b_local(:))**kappa_soap

return
end function kernel_soap

end module def_kernels
