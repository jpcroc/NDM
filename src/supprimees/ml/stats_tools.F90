module math

contains

  FUNCTION matdet(A) result(det)

    implicit none

    REAL(kind(0.d0)), dimension(:,:), intent(in) :: A
    REAL(kind(0.d0)) :: det

    IF ( (size(A,1).NE.3).AND.(size(A,2).NE.3) ) &
        STOP '< MatDet >: size of matrix to invert should be equal to 3'

    det =  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) &
         + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1) &
         - a(1,1)*a(2,3)*a(3,2) - a(1,2)*a(2,1)*a(3,3)

  END FUNCTION matdet


  SUBROUTINE matinv_gen(A, B)

    implicit none

    REAL(kind(0.d0)), dimension(:,:), intent(in) :: A
    REAL(kind(0.d0)), dimension(:,:), intent(out) :: B

    !REAL(kind(0.d0)) :: matdet
    REAL(kind(0.d0)) :: invdet

    ! Variables used with Lapack subroutine
    INTEGER :: i, n, INFO
    INTEGER, dimension(:), allocatable :: IPIV
    REAL(kind(0.d0)), dimension(:,:), allocatable :: tempA, invA
    if (size(A,1).NE.size(A,2)) &
         STOP '< MatInv >: matrix to invert is not square'

    SELECT CASE(size(A,1))
    CASE(1)
       B(1,1)=1.d0/A(1,1)

    CASE(2)
       invdet=1.d0/( A(1,1)*A(2,2)-A(1,2)*A(2,1) )
       B(1,1)=A(2,2)*invdet
       B(2,2)=A(1,1)*invdet
       B(1,2)=-A(2,1)*invdet
       B(2,1)=-A(1,2)*invdet

    CASE(3)
       invdet=1.d0/matdet(A)

       b(1,1) = a(2,2)*a(3,3) - a(2,3)*a(3,2)
       b(2,1) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
       b(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)

       b(1,2) = a(3,2)*a(1,3) - a(3,3)*a(1,2)
       b(2,2) = a(3,3)*a(1,1) - a(3,1)*a(1,3)
       b(3,2) = a(3,1)*a(1,2) - a(3,2)*a(1,1)

       b(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
       b(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
       b(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

       b(1:3,1:3)=b(1:3,1:3)*invdet

    CASE DEFAULT
!!$         STOP '< MatInv >: size of matrix to invert should be less than 3 &
!!$            &(lapack not implemented)'
       ! Use Lapack subroutine DGESV
       n=size(A,1)
       allocate(tempA(n,n))
       tempA(1:n,1:n)=A(1:n,1:n) ! Needed in order to not modify A
       allocate(IPIV(1:n))
       allocate(invA(n,n))
       invA(1:n,1:n)=0.d0
       do i=1,n
          invA(i,i)=1.d0
       enddo
       call DGESV(n, n, tempA, n, IPIV, invA, n, INFO)
       if (INFO.NE.0) then
          write(0,*) 'Result of DGESV subroutine: INFO=', INFO
          write(0,*)' < 0: if INFO = -i, the i-th argument had an illegal value'
          write(0,*)' > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization'
          write(0,*)'has been completed, but the factor U is exactly'
          write(0,*)'singular, so the solution could not be computed.'
          STOP '< MatInv >'
       endif
       B(1:n,1:n)=invA(1:n,1:n)
       deallocate(tempA, invA, IPIV)

    END SELECT

  END SUBROUTINE matinv_gen


subroutine correlation_coef(y1,y2,n,r)
implicit none
integer, intent(in)  :: n
real(kind=kind(1.d0)), intent(in) :: y1(n),y2(n)
real(kind=kind(1.d0)), intent(out):: r

!local
real(kind=kind(1.d0)) :: y1_m, y2_m, sigma_y1, sigma_y2




y1_m = 1.d0/dble(n)*SUM(y1(:))
y2_m = 1.d0/dble(n)*SUM(y2(:))

sigma_y1=dsqrt(SUM((y1(:)-y1_m)**2)/dble(n))
sigma_y2=dsqrt(SUM((y2(:)-y2_m)**2)/dble(n))


r = SUM((y1(:)-y1_m)*(y2(:)-y2_m))/dble(n)/(sigma_y1*sigma_y2)

return
end



!---------------------------------------
subroutine determination_coef(y_measure,y_predicted,n,r)
!---------------------------------------
implicit none
integer, intent(in)               :: n
real(kind=kind(1.d0)), intent(in) :: y_measure(n),y_predicted(n) ! y1-> database, y2-> fitted function
real(kind=kind(1.d0)), intent(out):: r
!local variables
real(kind=kind(1.d0)) :: ym_m, yp_m, sigma_ym, sigma_yp



ym_m = 1.d0/dble(n)*SUM(y_measure(:))
yp_m = 1.d0/dble(n)*SUM(y_predicted(:))

sigma_yp=dsqrt(SUM((y_predicted(:)-yp_m)**2)/dble(n))
sigma_ym=dsqrt(SUM((y_measure(:)-ym_m)**2)/dble(n))


r = SUM((y_predicted(:)-ym_m)**2) / SUM((y_measure(:)-ym_m)**2)

return
end subroutine determination_coef


!---------------------------------------
subroutine rmse_mae(y_measure,y_predicted,n,rmse,mae)
!---------------------------------------
implicit none
integer, intent(in)               :: n
real(kind=kind(1.d0)), intent(in) :: y_measure(n),y_predicted(n) ! y1-> database, y2-> fitted function
real(kind=kind(1.d0)), intent(out):: rmse,mae
!local variables



mae = 1.d0/dble(n)*SUM(dabs(y_measure(:)-y_predicted(:)))
rmse = dsqrt(1.d0/dble(n)*SUM((y_measure(:)-y_predicted(:))**2))


return
end subroutine rmse_mae

subroutine recentrate_database(nd_local_data,dim_yfunc,yfunc,yfunc_average)
implicit none
integer, intent(in)   :: nd_local_data,dim_yfunc
real(kind=kind(1.d0)) :: yfunc(dim_yfunc,nd_local_data), yfunc_average(dim_yfunc)

integer :: i


do i=1,dim_yfunc
   yfunc_average(i)=SUM(yfunc(i,:))/dble(nd_local_data)
end do

do i=1,dim_yfunc
  yfunc(i,:)=yfunc(i,:) - yfunc_average(i)
end do


return

end subroutine recentrate_database


subroutine print_message(iproc, text)
use ml_in_ndm_module, only : rangml
implicit none
integer, intent(in) :: iproc
character(len=150), intent(in) :: text

if (rangml==iproc) write(6,*) text

return

end subroutine print_message



integer function nitems2(line)
character,intent(in):: line*(*)
integer i, n, toks

i = 1;
n = len_trim(line)
toks = 0
nitems2 = 0
do while(i <= n)
   do while(line(i:i) == ' ')
     i = i + 1
     if (n < i) return
   enddo
   toks = toks + 1
   nitems2 = toks
   do
     i = i + 1
     if (n < i) return
     if (line(i:i) == ' ') exit
   enddo
enddo
end function nitems2

 ! number of space-separated items in a line
integer function nitems(line)
implicit none
    character line*(*)
    logical back
    integer length
    integer :: k

    back = .true.
    line=' '//line
    length = len_trim(line)
    k = index(line(1:length), ' ', back)
    !write(*,*) 'nnnnnniiit', k, len_trim(line)
    if (k == 0) then
        nitems = 0
        return
    end if

    nitems = 1
    do
        ! starting with the right most blank space,
        ! look for the next non-space character down
        ! indicating there is another item in the line
        do
            if (k <= 0) exit

            if (line(k:k) == ' ') then
                k = k - 1
                cycle
            else
                nitems = nitems + 1
                exit
            end if

        end do

        ! once a non-space character is found,
        ! skip all adjacent non-space character
        do
            if ( k<=0 ) exit

            if (line(k:k) /= ' ') then
                k = k - 1
                cycle
            end if

            exit

        end do

        if (k <= 0) exit

    end do
end function nitems

integer function itest_there_is_a_number (line)
 implicit none
 character :: line*(*)
 character(1) , dimension(10) :: numbers
 character(80) :: trimline
 integer :: i

 numbers(1:10)=(/ '0', '1', '2','3','4','5','6','7','8','9' /)
 trimline=trim(adjustl(line))
 itest_there_is_a_number=0
 do i=1,10
    if (numbers(i)==trimline(1:1)) itest_there_is_a_number=1

 end do
end function itest_there_is_a_number
end module  math
