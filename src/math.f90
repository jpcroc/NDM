module dk_math
    use dk_parameters, only: e

    implicit none(external, type)
    private

    !> Identity matrix
    real(e), dimension(3,3), parameter, public :: id3 = reshape([1,0,0,0,1,0,0,0,1],[3,3])

    public :: cross_product
    public :: matrix_inverse
    public :: det
    public :: approx_equal
contains
!******************************************************************************!
!> Determine whether two numbers are almost equal, i.e. at least one of their
!> absolute or relative differences is smaller than specific thresholds.
!******************************************************************************!
    elemental function approx_equal(x, y)
        real(e), intent(in) :: x
        real(e), intent(in) :: y
        logical :: approx_equal
!------
        integer, parameter :: k = min(kind(x), kind(y))
        ! real(k), parameter :: maxdiff = 10.0_k**(floor(log10(epsilon(0.0_k))))
        real(k), parameter :: maxdiff = 10.0_k**(-precision(0.0_k))
        real(k), parameter :: maxreldiff = 10.0_k ** (-min(precision(x), precision(y)))
        real(k) :: diff
!------
        if (x==y) then

            ! Do not even bother if the numbers are actually equal
            approx_equal = .true.
        else

            ! Otherwise, check whether the numbers are really close
            diff = real(abs(x-y), k)
            if (diff < maxdiff) then
                approx_equal = .true.
            else
                if (diff <= max(real(abs(x),k), real(abs(y),k)) * maxreldiff) then
                    approx_equal = .true.
                else
                    approx_equal = .false.
                end if
            end if
        end if

    end function
!******************************************************************************!
!> Calculate the cross product of 3-dimensional real vectors.
!******************************************************************************!
    pure function cross_product(u, v) result(w)
        real(e), dimension(3), intent(in) :: u
        real(e), dimension(3), intent(in) :: v
        real(e), dimension(3) :: w
!------
        w(1) = u(2)*v(3) - u(3)*v(2)
        w(2) = u(3)*v(1) - u(1)*v(3)
        w(3) = u(1)*v(2) - u(2)*v(1)
    end function
!******************************************************************************!
!> Calculate the determinant of a 3x3 real matrix.
!******************************************************************************!
     pure function det(M)
        real(e), dimension(3,3), intent(in) :: M    !! 3x3 matrix
        real(e) :: det
!------
        det = M(1,1) * ((M(2,2)*M(3,3)) - (M(2,3)*M(3,2)))
        det = det + M(2,1) * (M(3,2)*M(1,3) - M(3,3)*M(1,2))
        det = det + M(3,1) * ((M(1,2)*M(2,3)) - (M(1,3)*M(2,2)))
        return
    end function
!******************************************************************************!
!> Calculate the inverse of a 3x3 matrix.
!******************************************************************************!
    function matrix_inverse(a)
        real(e) , dimension(3,3), intent(in) :: a
        real(e) , dimension(3,3) :: matrix_inverse
!------
        integer :: i, im, ip, j, jm, jp
        real(e) :: deta, detai
!------

!----- Build the transposed cofactor matrix
        do j=1,3
            jm=mod(j+1,3)+1
            jp=mod(j,  3)+1
            do i=1,3
                im=mod(i+1,3)+1
                ip=mod(i,  3)+1

                matrix_inverse(i,j)=a(jp,ip)*a(jm,im)-a(jp,im)*a(jm,ip)
            end do
        end do

!----- Calculate the determinant
        deta=0.0_e
        do i=1,3
            deta=deta+a(i,1)*matrix_inverse(1,i)
        end do
        detai=1.0_e/deta

!----- Adjust the inverse marix
        do j=1,3
            do i=1,3
                matrix_inverse(i,j)=matrix_inverse(i,j)*detai
            end do
        end do
    end function
end module
