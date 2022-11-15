module tools

contains
!---------------------------------------------------
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

!-----------------------------------------------------
function cross_product(vecta,vectb) result(vectaxb)
real(kind(0.d0)), dimension(3), intent(in) :: vecta, vectb
real(kind(0.d0)), dimension(3) :: vectaxb

vectaxb(1) = vecta(2)*vectb(3) - vecta(3)*vectb(2)
vectaxb(2) = vecta(3)*vectb(1) - vecta(1)*vectb(3)
vectaxb(3) = vecta(1)*vectb(2) - vecta(2)*vectb(1)
end function cross_product

!-----------------------------------------------------
function norme(vect) result(norme_vect)
real(kind(0.d0)), dimension(3), intent(in) :: vect
real(kind(0.d0)) :: norme_vect
norme_vect = sqrt(dot_product(vect,vect))
end function norme

end module tools


!-------------------------------------------------
!-------------------------------------------------
subroutine dump_lammps(imm, im, itype, at_ang, xp_cart_ang)
!use derived_types_ph, only: config_real
use var_pot, only: ntyp
use phondy_in_ndm_module, only: debug_ph, rangph, passage, passage_inv

implicit none
integer, intent(in) :: imm, im,  itype(imm)
real(kind(0.d0)), intent(in) :: xp_cart_ang(3,imm), at_ang(3,3)

!internal
real(kind(0.d0)), dimension(3,3) :: ini_cell, new_cell_matrix, passage_matrix
real(kind(0.d0)), dimension(3) :: tmp_coord_i,new_tmp_coord_i
real(kind(0.d0)) :: xlo, ylo, zlo, xhi, yhi, zhi, xy, xz, yz
integer :: Nat, i


if (allocated(passage))     deallocate(passage)     ; allocate(passage(3,3))
if (allocated(passage_inv)) deallocate(passage_inv) ; allocate(passage_inv(3,3))


Nat = im
ini_cell = at_ang
!ini_cell = config_real(iconf)%cell
if (debug_ph) then
    if (rangph==0) then
        write(*,*) 'here is the big cell :', ini_cell
        write(*,*)'first column should be',ini_cell(:,1)
    end if
end if
call convert_cell(ini_cell,new_cell_matrix,passage_matrix)

if (debug_ph) then
    if (rangph==0) then
        write(*,*)'the passage matrix is :', passage_matrix
        write(*,*)'new_cell_matrix is :', new_cell_matrix
    end if
end if
passage = passage_matrix
call matinv_gen(passage, passage_inv)

!building header lammps
xlo = 0.d0
ylo = 0.d0
zlo = 0.d0
xhi = new_cell_matrix(1,1)
yhi = new_cell_matrix(2,2)
zhi = new_cell_matrix(3,3)
xy = new_cell_matrix(1,2)
xz = new_cell_matrix(1,3)
yz = new_cell_matrix(2,3)

if (rangph==0) then
    open(87, file='in.lmp',status='unknown')
    !write(87,*)'# add your comment ... '
    write(87,'(a)')'# add your comment ... '
    write(87,*)''
    write(87,'(I5,a)')Nat,' atoms'
    write(87,'(I1,a)')ntyp,' atom types'
    write(87,*)''
    write(87,'(f22.16,a,f22.16,a)')xlo,' ',xhi,' xlo xhi'
    write(87,'(f22.16,a,f22.16,a)')ylo,' ',yhi,' ylo yhi'
    write(87,'(f22.16,a,f22.16,a)')zlo,' ',zhi,' zlo zhi'
    write(87,'(f22.16,a,f22.16,a,f22.16,a)')xy,' ',xz,' ',yz,' xy xz yz'
    write(87,*)''
    write(87,'(a)')'Atoms'
    write(87,*)''

    !building coordinates
    do i=1,Nat
        tmp_coord_i = xp_cart_ang(:,i)
        new_tmp_coord_i = matmul(passage_matrix,tmp_coord_i)
        write(87,'(I5,a,I1,a,f22.16,a,f22.16,a,f22.16)')i,' ',itype(i),' ',new_tmp_coord_i(1),' ',new_tmp_coord_i(2),' ',new_tmp_coord_i(3)
    enddo
    close(87,status='keep')
end if

return
end subroutine  dump_lammps

!---------------------------------------------------
subroutine convert_cell(mat_ini,new_mat,transform)
use tools
use phondy_in_ndm_module, only: debug_ph

implicit none
real(kind(0.d0)), dimension(3,3), intent(in) :: mat_ini
real(kind(0.d0)), dimension(3,3), intent(out) :: new_mat, transform

!internal
real(kind(0.d0)), dimension(3,3) :: transit_cell,inv_mat_ini
real(kind(0.d0)), dimension(3) :: A,B,C, Ahat,AxBhat
real(kind(0.d0)) :: volume
logical :: upper, right
integer :: i

!matrix is already transpose
transit_cell = mat_ini
!write(*,*) 'transpose matrix is :',transit_cell

call is_upper_triangular(transit_cell,upper)
if (debug_ph) then
    write(*,*)'triangular up?',upper !debug
end if
if (.not.upper) then
  ! rotate bases into triangular matrix
    new_mat(:,:) = 0.d0
    A = transit_cell(:,1)
    B = transit_cell(:,2)
    C = transit_cell(:,3)
    call right_hand_basis(A,B,C,right)

    if (debug_ph) write(*,*)'direct ?',right !debug
    if (.not.right) then
        if (debug_ph) write(*,*)"WARNING: your reper is not right handed."
        if (debug_ph) write(*,*)"WARNING: This is a critical issue. The LAMMPS results are wrong !!!!!"
        stop
    end if

    new_mat(1,1) = norme(A)
    Ahat = A / norme(A)
    AxBhat = cross_product(A, B) / norme(cross_product(A, B))
    new_mat(1,2) = dot_product(B, Ahat)
    new_mat(2,2) = norme(cross_product(Ahat, B))
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



!---------------------------------------------------
subroutine right_hand_basis(vect1,vect2,vect3,l_right)
use tools

implicit none
real(kind(0.d0)), dimension(3), intent(in) :: vect1, vect2, vect3
logical, intent(out) :: l_right

l_right = .true.
if (dot_product(cross_product(vect1,vect2),vect3).lt.0) then
    l_right = .false.
endif

return
end subroutine right_hand_basis



!---------------------------------------------------
subroutine is_upper_triangular(M,l_triang)
implicit none
real(kind(0.d0)), dimension(3,3), intent(in)  :: M
logical, intent(out) :: l_triang
!write(*,*)'2,1',M(2,1) !debug
!write(*,*)'3,1',M(3,1) !debug
!write(*,*)'3,2',M(3,2) !debug

l_triang = .true.
if (abs(M(2,1).gt.1e-6).or.(abs(M(3,1).gt.1e-6)).or.(abs(M(3,2)).gt.1e-6)) then
    l_triang = .false.
endif

return
end subroutine is_upper_triangular


!-----------------------------------------------------

subroutine matinv_gen(A, B)
use tools
implicit none
REAL(kind(0.d0)), dimension(3,3), intent(in) :: A
REAL(kind(0.d0)), dimension(3,3), intent(out) :: B
REAL(kind(0.d0)) :: invdet
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
