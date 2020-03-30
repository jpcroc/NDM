module svd_mod
        implicit none
        contains
subroutine svd(entropy,n_temp)

USE T_kind_param_m, ONLY:  double
use ml_in_ndm_module
use temporary_data_cov, ONLY: xdesc, dim_xdesc

implicit none

integer, intent(in) :: n_temp
double precision,intent(out) :: entropy

integer :: dim_s,lwork,info
double precision, dimension(:), allocatable :: work,s,s_out,norm_spec
double precision, dimension(:,:), allocatable :: u,vt


lwork=max(3*min(dim_xdesc,n_temp)+max(dim_xdesc,n_temp),5*min(dim_xdesc,n_temp))
if (allocated(u)) deallocate(u)
allocate(u(dim_xdesc,dim_xdesc))
if (allocated(vt)) deallocate(vt)
allocate(vt(n_temp,n_temp))
if (allocated(s)) deallocate(s)
allocate(s(dim_xdesc))
if (allocated(work)) deallocate(work)
allocate(work(lwork))
call dgesvd('N', 'N', dim_xdesc, n_temp, xdesc, dim_xdesc, s, u, dim_xdesc, vt, n_temp, work, lwork, info)
dim_s=count(s>1d-16)
if (allocated(s_out)) deallocate(s_out)
allocate(s_out(dim_s))
s_out=pack(s,s>1d-16)
if (allocated(norm_spec)) deallocate(norm_spec)
allocate(norm_spec(dim_s))
norm_spec(:)=s(:)**2/sum(s(:)**2)
entropy=-sum(norm_spec(:)*log(norm_spec(:)))/size(norm_spec)

if (debug.and.(rangml==0)) then
write(*,*)'singular values'
write(*,*)s(:)
write(*,*)'singular values out'
write(*,*)s_out(:)
write(*,*)'normalised spectrum'
write(*,*)norm_spec(:)
write(*,*)'entropy'
write(*,*)entropy
endif

deallocate(u,vt,s,work,norm_spec)

return
end subroutine svd
end module
