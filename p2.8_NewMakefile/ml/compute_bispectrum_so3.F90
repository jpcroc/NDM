subroutine compute_bispectrum_so3(i_start_at,i_final_at,d_n_neigh, d_kind_neigh, local_bispectrum_so3_out,local_bispectrum_so3_deriv_out, iconf)

USE T_kind_param_m, ONLY:  double
use gen_com_m, ONLY: imm,A2cm,lperiod,bg,at,indi2
use tab_imm_m, ONLY : iwmax2,xp
use angular_functions, only : spherical_harm,grad_spherical_harm
use ml_in_ndm_module, ONLY: imm_neigh,  n_rbf,l_max,mconf, w2_rho,weighted, cg_vector, &
                            lbso3_diag, bisso3_dim, desc_forces, W_pow_so3, coeff_rbf, r_cut, factor_weight_mass, linvisible

use derived_types, only : config_real

implicit none
integer, intent (in) :: i_start_at,i_final_at
integer, dimension(imm),intent(out)  :: d_n_neigh
integer, dimension(imm,imm_neigh), intent(out) :: d_kind_neigh
double precision,dimension(bisso3_dim,imm),intent(out) :: local_bispectrum_so3_out
double precision,dimension(bisso3_dim,imm, 0:imm_neigh, 3),intent(out) :: local_bispectrum_so3_deriv_out
integer, optional :: iconf

logical :: small
real(double), dimension(:,:), allocatable :: xpnp
real(double), dimension(3) :: dxp_ji, ds
integer :: ia,ja,iw,iw1,iw2, ia_n, i_count


integer :: p,l,m
integer :: l1,l2,m1
real(double) :: r2_ji,r_ji


real(double), dimension(n_rbf) :: phi_ji, dphi_ji, rbf_ji, drbf_ji
real(double) :: factor_ia, factor_ja


!  double precision,dimension(n_rbf,ivoismax,imm) :: rbf,drbf
!powso3 double complex,dimension(-l_max:l_max,pow_so3_dim) :: c
!powso3 double complex,dimension(-l_max:l_max,pow_so3_dim,imm_neigh,3) :: dc
double complex,dimension(-l_max:l_max,0:l_max,n_rbf) :: c
double complex,dimension(-l_max:l_max,0:l_max,n_rbf,imm_neigh,3) :: dc
double complex, dimension(imm_neigh,3) :: dsum_bi
double complex :: sum_bi
integer :: l1_min, l1_max


if ((i_start_at==0) .and. (i_final_at==0)) then
  d_n_neigh(:)=0
  d_kind_neigh(:,:)=0
  local_bispectrum_so3_out(:,:)=0.d0
  local_bispectrum_so3_deriv_out(:,:,:,:)=0.d0
  return
end if

small=.false.
if (present(iconf)) then
small = config_real(iconf)%small
end if


ALLOCATE(xpnp(3,imm))
if (lperiod) then
 xpnp(:,:)=xp(:,:)
else
 call notperiod(xp,xpnp)
end if
! call cryst_to_cart (imm, xpnp, bg, -1)



d_n_neigh(:)=0
d_kind_neigh(:,:)=0
local_bispectrum_so3_out(:,:)=0.d0
local_bispectrum_so3_deriv_out(:,:,:,:)=0.d0



!  if (allocated(m_typ)) deallocate(m_typ)
!  allocate(m_typ(ntyp,ntyp))
!  call mconf(m_typ(:,:))

if (i_start_at==1)  iw2=0
if (i_start_at > 1) iw2=iwmax2(i_start_at-1)

do ja=i_start_at,i_final_at
     !begin small box or not 1/
     if (small) then
       iw1=1
       iw2=config_real(iconf)%n_neigh(ja)
     else
       iw1=iw2+1
       iw2=iwmax2(ja)
     end if
     !end   small box or not 1/


     if (linvisible.and.weighted) then
        if( config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ja)) ) cycle
     end if

     if (weighted) then
        !factor=w2_rho(massat(ityp(ja)),massat(ityp(ia)))
        !factor_ja=massat(ityp(ja))
        factor_ja=config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))!/factor_weight_mass
     else
        factor_ja=1.d0
     endif
     !tmp_pow_so3_out(:)=0.d0
     !tmp_pow_so3_deriv_out(:,:,:)=0.d0
     c(:,:,:)=(0.d0,0.d0)
     dc(:,:,:,:,:)=(0.d0,0.d0)
     ia_n=0
     do iw=iw1,iw2
        !begin small box or not 2/
        if (small) then
           ia = config_real(iconf)%kind_neigh(ja,iw)
        else
           ia=indi2(iw)
           if (ja==ia) cycle
        end if
        !end   small box or not 2/

        if (linvisible.and.weighted) then
            if( config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ia)) ) cycle
        end if

        if (weighted) then
           !factor_ia=massat(ityp(ia))
           factor_ia=config_real(iconf)%weight_per_type(config_real(iconf)%itype(ia))!/factor_weight_mass
        else
           factor_ia=1.d0
        endif

        if (small) then
           r_ji = config_real(iconf)%r_ij(ja,iw)
           r2_ji=r_ji**2
           dxp_ji(:) = config_real(iconf)%u_ij(ja,iw,:)
        else
           dxp_ji(1:3) = xpnp(1:3,ia) - xpnp(1:3,ja)
           ds=MatMul(dxp_ji,bg)
           WHERE ( (ds.GT.0.5d0).OR.(ds.LT.-0.5d0) )
              ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
           END WHERE
           dxp_ji(:) = MatMul(at(:,:),ds(:))/A2cm
           r2_ji = Sum( dxp_ji(1:3)**2 )
           r_ji = dsqrt(r2_ji)
        end if

        phi_ji(:)=0.d0
        if (r_ji >= r_cut) cycle
        ia_n = ia_n + 1

           do p=1,n_rbf
             phi_ji(p) = coeff_rbf(p)*(r_cut-r_ji)**(p+2.d0)*factor_ia
             if (desc_forces) dphi_ji(p)=-coeff_rbf(p)*(p+2.d0)*(r_cut-r_ji)**(p+1.d0)*factor_ia
           end do

           rbf_ji(:) =matmul(W_pow_so3(:,:) , phi_ji(:))
           if (desc_forces) drbf_ji(:)=matmul(W_pow_so3(:,:) ,dphi_ji(:))
           do p=1,n_rbf
              do l=0,l_max
                 do m=-l,l
                    c(m,l,p) = c(m,l,p) + rbf_ji(p) * spherical_harm(l,m,dxp_ji(:))
                    if (desc_forces) dc(m,l,p,ia_n,1:3) = dc(m,l,p,ia_n,1:3) +  &
                                                              drbf_ji(p)*(dxp_ji(1:3)/r_ji) * spherical_harm(l,m,dxp_ji(:)) + &
                                                                             rbf_ji(p) * grad_spherical_harm(l,m,dxp_ji(:))
                   !debug if ((ja==1).and.(ia==1)) then
                   !debug write(60, '(5i4, 4e20.10, 7e20.10)') ja, ia, p,l,m, dxp_ji(1:3), r_ji, rbf_ji(p),  grad_spherical_harm(l,m,dxp_ji(1:3))
                   !debug end if
                 enddo !p
              enddo    !l
           enddo       !m

     d_kind_neigh(ja,ia_n)=ia
     end do   !ia end of neighbours iterations ...
     d_n_neigh(ja)=ia_n

        i_count=0
        do p=1,n_rbf
          do l1=0,l_max
            if (lbso3_diag) then
              !GaborLike l2=l1
              l1_min=l1
              l1_max=l1
            else
              !T Aidan like do l2=0,l1
              l1_min=0
              l1_max=l1
            end if

            do l2=l1_min, l1_max
              do l=abs(l1-l2),min(l_max,l1+l2)
                if( mod(l1,2)==1 .and. mod(l2,2)==1 .and. mod(l,2)==1 ) cycle
                  i_count=i_count+1
                  !if (mod(l1+l2+l,2)==1) cycle ! assure invariance par réflexion
                  do m=-l,l
                    sum_bi=(0d0,0d0)
                    dsum_bi(:,:)=(0d0,0d0)
                    do m1=max(-l1,m-l2),min(l1,m+l2)
                      sum_bi = sum_bi + cg_vector(l1,m1,l2,m-m1,l,m)*c(m1,l1,p)*c(m-m1,l2,p)
                      if (desc_forces) then
                        do ia=1,ia_n
                           dsum_bi(ia,1:3) = dsum_bi(ia,1:3) + cg_vector(l1,m1,l2,m-m1,l,m)*(dc(m1,l1,p,ia,1:3)*c(m-m1,l2,p)+c(m1,l1,p)*dc(m-m1,l2,p,ia,1:3))
                        end do
                      end if
                    enddo

                    local_bispectrum_so3_out(i_count,ja) = local_bispectrum_so3_out(i_count,ja) + real(dconjg(c(m,l,p)) * sum_bi)*factor_ja
                    !debug write(60,'(2i5,3e15.7)') ja, i_count, c(m,l,p) * sum_bi, real(c(m,l,p) * sum_bi)

                    if (desc_forces) then
                       do ia=1,ia_n
                         local_bispectrum_so3_deriv_out(i_count,ja,ia,1:3) = local_bispectrum_so3_deriv_out(i_count,ja,ia,1:3) + &
                                                                                  real( dconjg(dc(m,l,p,ia,1:3)) * sum_bi + dconjg(c(m,l,p)) * dsum_bi(ia,1:3) )*factor_ja
                        end do
                    end if
                  enddo !m
                enddo !l
              enddo !l2
          enddo !l1
       enddo !p

       do ia=1,ia_n
         local_bispectrum_so3_deriv_out(:,ja,0,:) = local_bispectrum_so3_deriv_out(:,ja,0,:) - local_bispectrum_so3_deriv_out(:,ja,ia,:)
       end do

  end do !ja


 deallocate(xpnp)

return
end subroutine compute_bispectrum_so3



subroutine gen_dimension_for_bispectrum_so3_all
use ml_in_ndm_module, ONLY:  n_rbf, l_max, bisso3_l, bisso3_l1, bisso3_l2, bisso3_dim
implicit none
integer :: l1, l2, l, p
integer :: i_bi

i_bi=0
do p=1,n_rbf
do l1=0,l_max
  !l2=l1 ! following Gabor  only the diagonal elements are important
  do l2=0, l1  ! in the end we want only the componenets with l1 <= l2 <= l
     do l=abs(l1-l2),min(l_max,l1+l2)
       !if (mod(l1+l2+l,2)==1) cycle ! assure invariance par réflexion
       !if (l < l1) cycle  ! this comes from Thompson
       i_bi = i_bi+1
     enddo
  end do
enddo
enddo

bisso3_dim=i_bi
if (allocated(bisso3_l  )) deallocate(bisso3_l  ) ; allocate (bisso3_l(bisso3_dim))
if (allocated(bisso3_l1 )) deallocate(bisso3_l1 ) ; allocate (bisso3_l1(bisso3_dim))
if (allocated(bisso3_l2 )) deallocate(bisso3_l2 ) ; allocate (bisso3_l2(bisso3_dim))

i_bi=0
do p=1,n_rbf
do l1=0,l_max
     !l2=l1 ! following Gabor only diagonal elements are importants
     do l2=0, l1
        do l=abs(l1-l2),min(l_max,l1+l2)
          if( mod(l1,2)==1 .and. mod(l2,2)==1 .and. mod(l,2)==1 ) cycle
          !if (mod(l1+l2+l,2)==1) cycle ! assure invariance par réflexion
          !if (l < l1) cycle  ! this comes from Thompson
          !if (l < l1) cycle           ! in the end we want only the componenets with l1 <= l2 <= l
          i_bi = i_bi+1
          bisso3_l1 (i_bi) = l1
          bisso3_l2 (i_bi) = l2
          bisso3_l  (i_bi) = l
        end do
     end do
enddo
enddo
return
end subroutine gen_dimension_for_bispectrum_so3_all


subroutine  gen_dimension_for_bispectrum_so3_diagonal
use ml_in_ndm_module, ONLY: n_rbf, l_max, bisso3_l, bisso3_l1, bisso3_l2, bisso3_dim
implicit none
integer :: l1, l2, l, p
integer :: i_bi

i_bi=0
do p=1,n_rbf
do l1=0,l_max
     l2=l1 ! following Gabor  only the diagonal elements are important
     do l=abs(l1-l2),min(l_max,l1+l2)
       if( mod(l1,2)==1 .and. mod(l2,2)==1 .and. mod(l,2)==1 ) cycle
       !if (mod(l1+l2+l,2)==1) cycle ! assure invariance par réflexion
       !if (l < l1) cycle  ! this comes from Thompson
       i_bi = i_bi+1
     enddo
enddo
enddo

bisso3_dim=i_bi
if (allocated(bisso3_l  )) deallocate(bisso3_l  ) ; allocate (bisso3_l(bisso3_dim))
if (allocated(bisso3_l1 )) deallocate(bisso3_l1 ) ; allocate (bisso3_l1(bisso3_dim))
if (allocated(bisso3_l2 )) deallocate(bisso3_l2 ) ; allocate (bisso3_l2(bisso3_dim))

i_bi=0
do l1=0,l_max
     l2=l1 ! following Gabor only diagonal elements are importants
        do l=abs(l1-l2),min(l_max,l1+l2)
          if( mod(l1,2)==1 .and. mod(l2,2)==1 .and. mod(l,2)==1 ) cycle
          i_bi = i_bi+1
          bisso3_l1 (i_bi) = l1
          bisso3_l2 (i_bi) = l2
          bisso3_l  (i_bi) = l
        end do
enddo

return
end subroutine gen_dimension_for_bispectrum_so3_diagonal
