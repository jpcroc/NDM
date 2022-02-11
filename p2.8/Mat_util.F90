module Mat_utils_mod
   USE T_kind_param_m, ONLY:  double
   USE arret_ndm_mod,only:arret_ndm
   implicit none

  interface fillbuffer3D
     module procedure fillbuffer3Dreal
  end interface fillbuffer3D
  interface fillbuffer1D
     module procedure fillbuffer1Dinteger,fillbuffer1Dlogical,fillbuffer1Dreal
  end interface fillbuffer1D
  interface fillbuffer9D
     module procedure fillbuffer9Dreal
  end interface fillbuffer9D
contains
  !**********************************************************************
  subroutine matscl(a,s,b)
    !     Multiplies a scalar, s, to a 3-by-3 matrices a.
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double) , intent(in) :: a(3,3), s
    real(double) , intent(out) :: b(3,3)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: j, i
    !-----------------------------------------------
    do j = 1,3
       do i = 1,3
          b(i,j) = a(i,j)*s
       enddo
    enddo
    return
  end subroutine matscl

  !**********************************************************************
  !subroutine matmul(a,b,c)
  !  !     Multiplies 3-by-3 matrices a & b, and stores the result in c.
  !  !-----------------------------------------------
  !  USE T_kind_param_m, ONLY:  double
  !  implicit none
  !  !-----------------------------------------------
  !  !   D u m m y   A r g u m e n t s
  !  !-----------------------------------------------
  !  real(double) , intent(in) :: a(3,3), b(3,3)
  !  real(double) , intent(out) :: c(3,3)
  !  !-----------------------------------------------
  !  !   L o c a l   V a r i a b l e s
  !  !-----------------------------------------------
  !  integer :: j, i, k
  !  !-----------------------------------------------
  !  do j = 1,3
  !     do i = 1,3
  !        c(i,j) = 0d0
  !        do k = 1,3
  !           c(i,j) = c(i,j)+a(i,k)*b(k,j)
  !        enddo
  !     enddo
  !  enddo
  !  return
  !end subroutine matmul

  !**********************************************************************
  subroutine mattrp(a,at)
    !     Transposes a 3-by-3 matrix, a, and stores the result in at.
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double) , intent(in) :: a(3,3)
    real(double) , intent(out) :: at(3,3)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: j, i
    !-----------------------------------------------
    do j = 1,3
       do i = 1,3
          at(i,j) = a(j,i)
       enddo
    enddo
    return
  end subroutine mattrp

  !----------------------------------------------------------------------
subroutine matinv_gen(A, B)
!use tools
implicit none
real(double), dimension(3,3), intent(in) :: A
real(double), dimension(3,3), intent(out) :: B
real(double) :: invdet
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

return
end subroutine matinv_gen
FUNCTION matdet(A) result(det)
implicit none
real(double), dimension(:,:), intent(in) :: A
real(double) :: det
IF ( (size(A,1).NE.3).AND.(size(A,2).NE.3) ) &
    STOP '< MatDet >: size of matrix to invert should be equal to 3'
det =  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) &
       + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1) &
        - a(1,1)*a(2,3)*a(3,2) - a(1,2)*a(2,1)*a(3,3)
END FUNCTION matdet

  subroutine matinv(a, ai)
    !  Invert a 3-by-3 matrix a, and store the result in ai
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double) , intent(in) :: a(3,3)
    real(double) , intent(out) :: ai(3,3)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: j, jm, jp, i, im, ip
    real(double) :: deta, detai
    !-----------------------------------------------

    !      Transpose cofactor matrix
    do j = 1, 3
       jm = mod(j+1,3)+1
       jp = mod(j,  3)+1
       do i = 1, 3
          im = mod(i+1,3)+1
          ip = mod(i,  3)+1
          ai(i,j) = a(jp,ip)*a(jm,im)-a(jp,im)*a(jm,ip) ! transposed cofactor matrix
       enddo
    enddo

    !      Determinant
    deta = 0.d0
    do i = 1, 3
       deta = deta+a(i,1)*ai(1,i)
    enddo

    detai = 1d0/deta

    !      Inverse matrix
    do j = 1, 3
       do i = 1, 3
          ai(i,j) = ai(i,j)*detai
       enddo
    enddo

    return
  end subroutine matinv

  !**********************************************************************
  subroutine matcof(a,c)
    !--------------------------------------------------------------------
    !  Given a 3-by-3 matrix a, calculates its cofactor matrix and stores
    !  the result in c.
    !--------------------------------------------------------------------
    USE T_kind_param_m, ONLY:  double
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double) , intent(in) :: a(3,3)
    real(double) , intent(out) :: c(3,3)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: j, jm, jp, i, im, ip
    !-----------------------------------------------
    do j = 1,3
       jm = mod(j+1,3)+1
       jp = mod(j,  3)+1
       do i = 1,3
          im = mod(i+1,3)+1
          ip = mod(i,  3)+1
          c(i,j) = a(ip,jp)*a(im,jm)-a(im,jp)*a(ip,jm)
       enddo
    enddo
    return
  end subroutine matcof

  !**********************************************************************
  real(kind(0.0d0)) function detmat(a)
    !  Returns the determinant of a 3-by-3 matrices a.
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double) , intent(in) :: a(3,3)
    !-----------------------------------------------
    detmat = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
         -a(2,1)*(a(1,2)*a(3,3)-a(1,3)*a(3,2))   &
         +a(3,1)*(a(1,2)*a(2,3)-a(1,3)*a(2,2))
    return
  end function detmat

  !**********************************************************************
  subroutine boxmat(at,volu)
    !     Sets up matrices allocated with the MD box.
    !----------------------------------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY: ati
    ! *********************************************************************
    implicit none
    real(double), dimension(3,3) :: at
    real(double)::volu
    
    real(double), dimension(3,3) :: atit, sgm, ah

    !-----Areal tensor, SGM
    call mattrp(at,ah)
    call matcof(ah,sgm)
    !-----MD-box volumeb
    volu=ah(1,1)*sgm(1,1)+ah(2,1)*sgm(2,1)+ah(3,1)*sgm(3,1)
    call matscl(sgm,1d0/volu,atit)
    !-----Inverse MD-box tensor, HI
    call mattrp(atit,ati)

    return
  end subroutine boxmat

  !**********************************************************************
  subroutine eigen3(a,d,v)
    !----------------------------------------------------------------------
    !  Diagonalizes a real symmetric 3x3 matrix.
    !     a(3,3): Input matrix
    !     d(3):   Return eigenvalues
    !     v(3,3): Return eigenvectors--v(*,i)  =  i-th eigenvector
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double)  :: a(3,3),astock(3,3)
    real(double)  :: v(3,3), d(3)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: j, i
    !-----------------------------------------------

    !-----Find eigenvalues
    call jacobi(a,3,3,d,v,i)
    !-----Restore the original matrix
    astock(:,:)=a(:,:)
    do i = 1,3
       do j = i+1,3
          a(i,j) = astock(j,i)
       enddo
    enddo
    return
  end subroutine eigen3

  !**********************************************************************
  subroutine jacobi(a,n,np,d,v,nrot)
    !     Diagonalizes a real symmetric matrix by Jacobi transformation.
    !     From "Numerical Recipes".
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer :: n, np, nrot
    real(double)  :: a(np,np)
    real(double)  :: v(np,np),d(np)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: j, i, iq, ip, nmax
    parameter (nmax = 500)
    real(double) :: sm, tresh, g, h, t, theta, c, s, tau
    real(double) :: b(nmax), z(nmax)
    !-----------------------------------------------
    do ip = 1,n
       do iq = 1,n
          v(ip,iq) = 0.
       enddo
       v(ip,ip) = 1.
    enddo
    do ip = 1,n
       b(ip) = a(ip,ip)
       d(ip) = b(ip)
       z(ip) = 0.
    enddo
    nrot = 0
    do i = 1,50
       sm = 0.
       do ip = 1,n-1
          do iq = ip+1,n
             sm = sm+abs(a(ip,iq))
          enddo
       enddo
       if(sm==0.)return
       if(i<4)then
          tresh = 0.2*sm/n**2
       else
          tresh = 0.
       endif
       do ip = 1,n-1
          do iq = ip+1,n
             g = 100.*abs(a(ip,iq))
             if((i>4).and.(abs(d(ip))+   &
                  g==abs(d(ip))).and.(abs(d(iq))+g==abs(d(iq)))) then
                a(ip,iq) = 0.
             else if(abs(a(ip,iq))>tresh)then
                h = d(iq)-d(ip)
                if(abs(h)+g==abs(h))then
                   t = a(ip,iq)/h
                else
                   theta = 0.5*h/a(ip,iq)
                   t = 1./(abs(theta)+sqrt(1.+theta**2))
                   if(theta<0.)t = -t
                endif
                c = 1./sqrt(1+t**2)
                s = t*c
                tau = s/(1.+c)
                h = t*a(ip,iq)
                z(ip) = z(ip)-h
                z(iq) = z(iq)+h
                d(ip) = d(ip)-h
                d(iq) = d(iq)+h
                a(ip,iq) = 0.
                do j = 1,ip-1
                   g = a(j,ip)
                   h = a(j,iq)
                   a(j,ip) = g-s*(h+g*tau)
                   a(j,iq) = h+s*(g-h*tau)
                enddo
                do j = ip+1,iq-1
                   g = a(ip,j)
                   h = a(j,iq)
                   a(ip,j) = g-s*(h+g*tau)
                   a(j,iq) = h+s*(g-h*tau)
                enddo
                do j = iq+1,n
                   g = a(ip,j)
                   h = a(iq,j)
                   a(ip,j) = g-s*(h+g*tau)
                   a(iq,j) = h+s*(g-h*tau)
                enddo
                do j = 1,n
                   g = v(j,ip)
                   h = v(j,iq)
                   v(j,ip) = g-s*(h+g*tau)
                   v(j,iq) = h+s*(g-h*tau)
                enddo
                nrot = nrot+1
             endif
          enddo
       enddo
       do ip = 1,n
          b(ip) = b(ip)+z(ip)
          d(ip) = b(ip)
          z(ip) = 0.
       enddo
    enddo
    write(6,*) 'too many iterations in jacobi'
    return
  end subroutine jacobi

  !**********************************************************************
  subroutine myrnd(rnd,dseed)
    !----------------------------------------------------------------------
    !  Random-number generator.
    !----------------------------------------------------------------------
    USE T_kind_param_m, ONLY:  double
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double) :: dseed
    real(double) ::rnd
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    real(double) :: d2p31m, d2p31
    data d2p31m/2147483647d0/
    data d2p31 /2147483648d0/
    !-----------------------------------------------

    dseed = dmod(16807d0*dseed,d2p31m)
    rnd   = dseed/d2p31
    return
  end subroutine myrnd

  function cross_product(vecta,vectb) result(vectaxb) ! a intégrer
    USE T_kind_param_m, ONLY:  double
    real(double), dimension(3), intent(in) :: vecta, vectb
    real(double), dimension(3) :: vectaxb

    vectaxb(1) = vecta(2)*vectb(3) - vecta(3)*vectb(2)
    vectaxb(2) = vecta(3)*vectb(1) - vecta(1)*vectb(3)
    vectaxb(3) = vecta(1)*vectb(2) - vecta(2)*vectb(1)
  end function cross_product

  !-----------------------------------------------------
  function norme(vect) result(norme_vect) ! a intégrer
    USE T_kind_param_m, ONLY:  double
    real(double), dimension(3), intent(in) :: vect
    real(double) :: norme_vect
    norme_vect = sqrt(dot_product(vect,vect))
  end function norme


  subroutine right_hand_basis(vect1,vect2,vect3,l_right) ! inutile a inliner
    USE T_kind_param_m, ONLY:  double
    implicit none

    real(double), dimension(3), intent(in) :: vect1, vect2, vect3
    logical, intent(out) :: l_right
    real(double), dimension(3)::cp
    l_right = .true.
    cp=cross_product(vect1,vect2)
    if (dot_product(cp,vect3).lt.0) then
       l_right = .false.
    endif

    return
  end subroutine right_hand_basis


  !---------------------------------------------------
  subroutine is_upper_triangular(M,l_triang)
    USE T_kind_param_m, ONLY:  double
    implicit none
    real(double), dimension(3,3), intent(in)  :: M
    logical, intent(out) :: l_triang
    !write(*,*)'2,1',M(2,1) !debug
    !write(*,*)'3,1',M(3,1) !debug
    !write(*,*)'3,2',M(3,2) !debug
    real(double)::v1,v2,v3
    v1=abs(M(2,1))
    v2=abs(M(3,1))
    v3=abs(M(3,2))
    l_triang = .true.
    if ((v1.gt.1.d-16).or.(v2.gt.1.d-16).or.(v3.gt.1.d-16)) then
       l_triang = .false.
    endif

    return
  end subroutine is_upper_triangular

  subroutine fillbuffer3Dreal(buf,vect,masq)
    real(double),allocatable,intent(in)::vect(:,:)
    real(double),allocatable::buf(:,:)
    logical, allocatable,intent(in)::masq(:)
    integer::i,iloc
    iloc=0
    buf=0
!    write(6,*)'masq',size(masq),size(vect),size(buf)
    do i=1,size(masq)
 !      write(6,*)i,masq(i),vect(:,i),iloc
       if (masq(i))then
          iloc=iloc+1
          buf(:,iloc)=vect(:,i)
       end if
    end do

  end subroutine fillbuffer3Dreal
  
  subroutine fillbuffer9Dreal(buf,vect,masq)
    real(double),allocatable,intent(in)::vect(:,:,:)
    real(double),allocatable::buf(:,:,:)
    logical, allocatable,intent(in)::masq(:)
    integer::i,iloc
    iloc=0
    buf=0
!    write(6,*)'masq',size(masq),size(vect),size(buf)
    do i=1,size(masq)
 !      write(6,*)i,masq(i),vect(:,i),iloc
       if (masq(i))then
          iloc=iloc+1
          buf(:,:,iloc)=vect(:,:,i)
       end if
    end do

  end subroutine fillbuffer9Dreal

  subroutine fillbuffer1Dinteger(ibuf,vect,masq)
    integer,allocatable,intent(in)::vect(:)
    integer,allocatable::ibuf(:)
    logical, allocatable,intent(in)::masq(:)
    integer::i,iloc
    iloc=0
    do i=1,size(masq)
       if (masq(i))then
          iloc=iloc+1
          ibuf(iloc)=vect(i)
       end if
    end do
  end subroutine fillbuffer1Dinteger
  subroutine fillbuffer1Dreal(rbuf,vect,masq)
    real(double),allocatable,intent(in)::vect(:)
    real(double),allocatable::rbuf(:)
    logical, allocatable,intent(in)::masq(:)
    integer::i,iloc
    iloc=0
    do i=1,size(masq)
       if (masq(i))then
          iloc=iloc+1
          rbuf(iloc)=vect(i)
       end if
    end do
  end subroutine fillbuffer1Dreal

  subroutine matmult3(A,b,c)
    real(double),intent(in)::a(3,3),b(3)
    real(double),intent(out)::c(3)
    integer::in
    c(1)=a(1,1)*b(1)+a(1,2)*b(2)+a(1,3)*b(3)
    c(2)=a(2,1)*b(1)+a(2,2)*b(2)+a(2,3)*b(3)
    c(3)=a(3,1)*b(1)+a(3,2)*b(2)+a(3,3)*b(3)
  end subroutine matmult3
  
  subroutine fillbuffer1Dlogical(lbuf,vect,masq)
    logical,allocatable,intent(in)::vect(:)
    logical,allocatable::lbuf(:)
    logical, allocatable,intent(in)::masq(:)
    integer::i,iloc
    iloc=0
    do i=1,size(masq)
       if (masq(i))then
          iloc=iloc+1
          lbuf(iloc)=vect(i)
       end if
    end do
  end subroutine fillbuffer1Dlogical
  !---------------------------------------------------
  subroutine convert_cell(mat_ini,new_mat,transform)



    implicit none
    real(kind(0.d0)), dimension(3,3), intent(in) :: mat_ini
    real(kind(0.d0)), dimension(3,3), intent(out) :: new_mat, transform

    !internal
    real(kind(0.d0)), dimension(3,3) :: transit_cell,inv_mat_ini
    real(kind(0.d0)), dimension(3) :: A,B,C, Ahat,AxBhat,Mtmp(3)
    real(kind(0.d0)) :: volume
    logical :: upper, right
    integer :: i

    !matrix is already transpose
    transit_cell = mat_ini
    !write(*,*) 'transpose matrix is :',transit_cell

    call is_upper_triangular(transit_cell,upper)
    if (.not.upper) then
       ! rotate bases into triangular matrix
       new_mat(:,:) = 0.d0
       A = transit_cell(:,1)
       B = transit_cell(:,2)
       C = transit_cell(:,3)
       call right_hand_basis(A,B,C,right)

       if (.not.right) then
          write(*,*)"WARNING: your reper is not right handed."
          write(*,*)"WARNING: This is a critical issue. The LAMMPS results are wrong !!!!!"
          call arret_ndm
       end if

       new_mat(1,1) = norme(A)
       Ahat = A / norme(A)
       Mtmp=cross_product(A, B)
       AxBhat = cross_product(A, B) / norme(Mtmp)
       new_mat(1,2) = dot_product(B, Ahat)
       Mtmp=cross_product(Ahat, B)
       new_mat(2,2) = norme(Mtmp)
       new_mat(1,3) = dot_product(C,Ahat)
       new_mat(2,3) = dot_product(C,cross_product(AxBhat, Ahat))
       new_mat(3,3) = abs(dot_product(C, AxBhat))
       !create and save the transformation for coordinates
       !volume = matdet(mat_ini)
       !trans = np.array([np.cross(B, C), np.cross(C, A), np.cross(A, B)])
       !trans = trans / volume
       !coord_transform = np.dot(tri_mat , trans)
       call matinv_gen(mat_ini,inv_mat_ini)
       transform = matmul(new_mat,inv_mat_ini)

    else
       new_mat = mat_ini
       transform(:,:) = 0.d0
       do i=1,3
          transform(i,i) = 1.d0
       enddo

    endif
    return
  end subroutine convert_cell


  subroutine reslin33(x,b,a)
    real(double),intent(out)::x(3)
    real(double),intent(in)::b(3),a(3,3)
    real(double)::deta,detab(3),ab1(3,3),ab2(3,3),ab3(3,3),am1(3,3)
!!$    deta=matdet(a)
!!$
!!$    ab1=a
!!$    ab1(:,1)=b(:)
!!$    detab(1)=matdet(ab1)
!!$
!!$    ab2=a
!!$    ab2(:,2)=b(:)
!!$    detab(2)=matdet(ab2)
!!$
!!$    ab3=a
!!$    ab3(:,3)=b(:)
!!$    detab(3)=matdet(ab3)
!!$    
!!$    x(:)=detab(:)/deta
!!$    write(6,*)'X1',x
    call matinv(a,am1)
    call matmult3(am1,b,x)
!    write(6,*)'X2',x
  end subroutine reslin33
end module Mat_utils_mod
