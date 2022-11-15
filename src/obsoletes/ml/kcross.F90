module   temporary_data_cov
implicit none

 integer :: dim_xdesc
 integer :: dim_extra,  dim_train, dim_valid

! dim_data - dimension of the database
! dim_data_train + dim_data_valid= dim_data
 integer :: dim_data, dim_data_train, dim_data_test
 integer :: dim_data_constraints

 real(kind=kind(1.d0)),dimension(:), allocatable    :: yfunc_average


 integer, parameter :: dim_yfunc=1
 real(kind=kind(1.d0)), dimension(:),   allocatable :: yfunc        ! vector of dimension M, size of the database
 real(kind=kind(1.d0)), dimension(:),   allocatable :: yfunc_train, yfunc_valid ! vector of dimension M, size of the database
 real(kind=kind(1.d0)), dimension(:),   allocatable :: y_test
 real(kind=kind(1.d0)), dimension(:),   allocatable :: y_extra  ! vector of dimension M, size of the database

 !if dim_yfunc > 1 probably the following functions should be used:

 real(kind=kind(1.d0)), dimension(:,:),   allocatable :: yfunc_nd        ! vector of dimension YxM, M = size of the database, Y the componenets of the  y function



 real(kind=kind(1.d0)), dimension(:,:), allocatable :: matfor
 real(kind=kind(1.d0)), dimension(:,:), allocatable :: xdesc ! vetor of dimension (dim_desc, M) desc size x size of the database
 real(kind=kind(1.d0)), dimension(:,:), allocatable :: xdesc_train, xdesc_valid, xdesc_test   ! vetor of dimension (dim_desc, P) desc size x size of the database
 real(kind=kind(1.d0)), dimension(:),   allocatable :: xdesc_average   ! center of data xdesc(:,dim_data)
 real(kind=kind(1.d0)), dimension(:),   allocatable :: error_train, error_valid, error_test
 real(kind=kind(1.d0)), dimension(:),   allocatable :: mlocal,kbfunc,eigen_values
 real(kind=kind(1.d0)), dimension(:),   allocatable :: d_mlocal                   ! the values of the derivatives respect with the length of kernel
 integer, dimension(:), allocatable :: u_local, v_local

 integer :: i_start_cov, i_final_cov,i_local_cov
 integer, dimension(:), allocatable :: neigh_cov,list_cov

 ! MD objects ... for snap
 ! energy MD descriptor per  box with the shape (dim_snap = dim_xdesc+1)
 real(kind=kind(1.d0)), dimension(:), allocatable   :: xdesc_emd
 ! force MD descriptor per atom (d,3,Nat) shape
 real(kind=kind(1.d0)), dimension(:,:,:), allocatable :: xdesc_fmd

 ! stress MD descriptor per atom (d,6) shape
 real(kind=kind(1.d0)), dimension(:,:), allocatable :: xdesc_smd

end module temporary_data_cov




module k_cross_validation

logical :: kcross
integer :: n_kcross
integer, allocatable, dimension(:) :: dim_kcell_valid, dim_kcell_train, imin_kcell, imax_kcell, rdm

contains


subroutine set_kcross(nd_local_data,rangml)
! subroutine which span allocate the dimension of each space i from 1 to n_kcross for the Kcross procedure.
!
! in input:  nd_local_data the dimension of the data nd_local_data
!            rangml        the rang of the core
!            n_kcross      the dimension of the kcross process.
!
! in output output:
! - the dimensions of the validation and trainning for each i between i=1,n_kcross
!   dim_kcell_valid(i), dim_kcell_train(i) for i =1, n_kcross
! - the limits of the 1<= i^th <= n_kcross cell given by the vectors
!   imim_kcell(i) and imax_kcell(i)

implicit none
integer, intent(in) :: nd_local_data, rangml
integer :: i

  if (allocated(rdm)) deallocate(rdm); allocate (rdm(nd_local_data))
  call shuffle(nd_local_data,rdm)

   if (allocated(dim_kcell_train))  deallocate(dim_kcell_train); allocate (dim_kcell_train(n_kcross))
   if (allocated(dim_kcell_valid))  deallocate(dim_kcell_valid); allocate (dim_kcell_valid(n_kcross))
   if (allocated(imin_kcell))       deallocate(imin_kcell)     ; allocate (imin_kcell(n_kcross))
   if (allocated(imax_kcell))       deallocate(imax_kcell)     ; allocate (imax_kcell(n_kcross))

 if ((n_kcross<=0) .or. (n_kcross==1)) then
  if (rangml==0) write(6,*) 'ML: This program is not implemented for n_kcross', n_kcross
  if (rangml==0) write(6,*) 'ML: Legal value for n_kcross should be an integer > 1. Usually 7 is a good value'
  stop
 end if

 if (n_kcross>1) then
   do i=1,n_kcross-1
      dim_kcell_train(i)=nd_local_data-int(nd_local_data/n_kcross)
      dim_kcell_valid(i)=int(nd_local_data/n_kcross)
      imin_kcell(i)=(i-1)*int(nd_local_data/n_kcross) + 1
      imax_kcell(i)=i*int(nd_local_data/n_kcross)
   end do
     dim_kcell_valid(n_kcross)=int(nd_local_data/n_kcross)+mod(nd_local_data,n_kcross)
     dim_kcell_train(n_kcross)=nd_local_data-dim_kcell_valid(n_kcross)
     imin_kcell(n_kcross)=(n_kcross-1)*int(nd_local_data/n_kcross)+1
     imax_kcell(n_kcross)=nd_local_data

 else
     dim_kcell_valid(n_kcross)=0
     dim_kcell_train(n_kcross)=nd_local_data
     imin_kcell(n_kcross)=nd_local_data
     imax_kcell(n_kcross)=nd_local_data
 end if

!debug  do i =1,n_kcross
!debug   write(*,*) i, dim_kcell_train(i), dim_kcell_valid(i), imin_kcell(i), imax_kcell(i)
!debug  end do


return

end subroutine set_kcross

subroutine fill_xdesc_yfunc_kcross(ik,rangml)
use temporary_data_cov, ONLY: dim_data_train, dim_xdesc, &
                              xdesc,       yfunc ,      &
                              xdesc_train, yfunc_train, &
                              xdesc_valid, yfunc_valid, error_valid, &
                              dim_train, dim_valid
implicit none
integer, intent(in) :: ik,rangml
integer :: i,icount1, icount2

  if (allocated(xdesc_train)) deallocate(xdesc_train)
  allocate(xdesc_train(dim_xdesc, dim_kcell_train(ik)))
  if (allocated(xdesc_valid)) deallocate(xdesc_valid)
  allocate(xdesc_valid(dim_xdesc, dim_kcell_valid(ik)))

  if (allocated(yfunc_train)) deallocate(yfunc_train)
  allocate(yfunc_train(dim_kcell_train(ik)))
  if (allocated(yfunc_valid)) deallocate(yfunc_valid)
  allocate(yfunc_valid(dim_kcell_valid(ik)))

  icount1=0
  icount2=0
  do i=1,dim_data_train

    if ( (i>=imin_kcell(ik)) .and. (i<=imax_kcell(ik)) ) then
     icount2=icount2+1
     xdesc_valid(1:dim_xdesc, icount2)=xdesc(1:dim_xdesc, rdm(i))
     yfunc_valid(icount2) = yfunc(rdm(i))
!     xdesc_valid(1:dim_xdesc, icount2)=xdesc(1:dim_xdesc, i)
!     yfunc_valid(icount2) = yfunc(i)
    else
     icount1=icount1+1
     xdesc_train(1:dim_xdesc, icount1)=xdesc(1:dim_xdesc, rdm(i))
     yfunc_train(icount1) = yfunc(rdm(i))
!     xdesc_train(1:dim_xdesc, icount1)=xdesc(1:dim_xdesc, i)
!     yfunc_train(icount1) = yfunc(i)
    end if
  end do !

!debug    write(*,*) 'ml train', ik, icount1, dim_kcell_train(ik)
!debug    write(*,*) 'ml valid', ik, icount2, dim_kcell_valid(ik)

  if (icount1 /= dim_kcell_train(ik) ) then
     if (rangml==0) write(6,*) 'error in kcross_validation train. Wrong partition train-valid', ik, icount1, dim_kcell_train(ik)
     stop
  end if

  if (icount2 /= dim_kcell_valid(ik) ) then
    if (rangml==0) write(6,*) 'error in kcross_validation validation.Wrong partition train-valid', ik, icount2, dim_kcell_valid(ik)
    stop
  end if

  dim_train=dim_kcell_train(ik)
  dim_valid=dim_kcell_valid(ik)

  if (allocated(error_valid)) deallocate(error_valid)
  allocate(error_valid(dim_valid))

!debug   write(*,*) 'im',  dim_train, dim_valid, imin_kcell(ik), imax_kcell(ik)

return

end subroutine fill_xdesc_yfunc_kcross



subroutine fill_xdesc_yfunc_oneshot(rangml)
use temporary_data_cov, ONLY: dim_data_train, dim_xdesc, &
                              xdesc,       yfunc ,      &
                              xdesc_train, yfunc_train, &
                              xdesc_valid, yfunc_valid, error_valid, &
                              dim_train, dim_valid
implicit none
integer, intent(in) :: rangml
integer :: i,icount1, icount2


  if (allocated(rdm)) deallocate(rdm); allocate (rdm(dim_data_train))
  call shuffle(dim_data_train,rdm)



  if (allocated(xdesc_train)) deallocate(xdesc_train) ; allocate(xdesc_train(dim_xdesc, dim_train))
  if (allocated(xdesc_valid)) deallocate(xdesc_valid) ; allocate(xdesc_valid(dim_xdesc, dim_valid))

  if (allocated(yfunc_train)) deallocate(yfunc_train) ; allocate(yfunc_train(dim_train))
  if (allocated(yfunc_valid)) deallocate(yfunc_valid) ; allocate(yfunc_valid(dim_valid))

  icount1=0
  icount2=0
  do i=1,dim_data_train

    if (i<=dim_valid) then
     icount2=icount2+1
     xdesc_valid(1:dim_xdesc, icount2)=xdesc(1:dim_xdesc, rdm(i))
     yfunc_valid(icount2) = yfunc(rdm(i))
    else
     icount1=icount1+1
     xdesc_train(1:dim_xdesc, icount1)=xdesc(1:dim_xdesc, rdm(i))
     yfunc_train(icount1) = yfunc(rdm(i))
    end if

  end do !


  if (icount1 /= dim_train) then
     if (rangml==0) write(6,*) 'ML: < fill_xdesc_yfunc_oneshot > error wrong partition train-valid', icount1, dim_train
     stop
  end if

  if (icount2 /= dim_valid) then
    if (rangml==0) write(6,*) 'ML: < fill_xdesc_yfunc_oneshot > error wrong partition train-valid', icount2, dim_valid
    stop
  end if


  if (allocated(error_valid)) deallocate(error_valid)
  allocate(error_valid(dim_valid))


return

end subroutine fill_xdesc_yfunc_oneshot






subroutine allocate_xdesc (n_frac, rangml)
! take the dim_data points which are splitted into dim_data_train and dim_data_test.
use temporary_data_cov, ONLY: dim_data_train, dim_data_test, dim_data, &
                              dim_xdesc, &
                              xdesc,       xdesc_test


implicit none
real(8), intent(in) :: n_frac
integer, intent(in) :: rangml



 if (n_frac.gt.0d0) then
    dim_data_train=ceiling(dim_data*(1.d0-n_frac))
    dim_data_test=dim_data - dim_data_train

    if (allocated(xdesc_test)) deallocate(xdesc_test); allocate(xdesc_test(dim_xdesc,dim_data_test))
    xdesc_test(:,:)=xdesc(1:dim_xdesc,dim_data_train+1:dim_data)

 else
    dim_data_train=dim_data
    dim_data_test=0
 endif

    if (rangml==0) then
           write(*,'("ML: <allocate_xdesc> dim_data       ",i9)') dim_data
           write(*,'("ML: <allocate_xdesc> dim_train      ",i9)') dim_data_train
           write(*,'("ML: <allocate_xdesc> dim_test       ",i9)') dim_data_test
    end if


return

end subroutine allocate_xdesc



subroutine shuffle(n,a)

implicit none

integer,intent(in) :: n
integer,dimension(n),intent(out) :: a
integer :: i, randpos, temp
real :: r

a=(/ (i, i=1,n) /)

do i = n,2,-1
  call random_number(r)
  randpos = int(r * i) + 1
  temp = a(randpos)
  a(randpos) = a(i)
  a(i) = temp
end do

return
end subroutine shuffle

end module k_cross_validation
