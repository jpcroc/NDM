module compute_pow_so3_mod
        use notperiod_mod
        implicit none
        contains
subroutine compute_pow_so3(i_start_at,i_final_at,d_n_neigh, d_kind_neigh, local_pow_so3_out,local_pow_so3_deriv_out, iconf)

 USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: imm,A2cm,lperiod,bg,at,indi2
  use tab_imm_m, ONLY : iwmax2,xp
  use angular_functions, only : spherical_harm,grad_spherical_harm
  use derived_types, only: config_real
  use ml_in_ndm_module, ONLY: imm_neigh, one_pi,  r_cut,n_rbf,l_max,  w2_rho,weighted, &
                              pow_so3_dim, desc_forces, W_pow_so3, coeff_rbf, factor_weight_mass, linvisible, weighted

 implicit none
  integer, intent (in) :: i_start_at,i_final_at
  !double precision,dimension(0:l_max,n_rbf,n2_typ,imm),intent(out) :: local_pow_so3_out
  !double precision,dimension(0:l_max,n_rbf,3,n2_typ,imm),intent(out) :: local_pow_so3f_out
  integer, dimension(imm),intent(out)  :: d_n_neigh
  integer, dimension(imm,imm_neigh), intent(out) :: d_kind_neigh
  double precision, dimension(pow_so3_dim, imm),intent(out) :: local_pow_so3_out
  double precision, dimension(pow_so3_dim, imm ,0:imm_neigh,3),intent(out) :: local_pow_so3_deriv_out
  integer, optional :: iconf

  logical :: small
  real(double), dimension(:,:), allocatable :: xpnp
  real(double), dimension(3) :: dxp_ji, ds
  integer :: ia,ja,iw,iw1,iw2, ia_n, i_count

  integer :: p,l,m
  real(double) :: r2_ji,r_ji
  complex(kind=8), dimension(-l_max:l_max,pow_so3_dim) :: c
  complex(kind=8), dimension(-l_max:l_max,pow_so3_dim,imm_neigh,3) :: dc

  real(double), dimension(pow_so3_dim) :: tmp_pow_so3_out
  real(double), dimension(pow_so3_dim, 0:imm_neigh,3):: tmp_pow_so3_deriv_out

  real(double), dimension(n_rbf) :: phi_ji, dphi_ji, rbf_ji, drbf_ji
  real(double) :: factor_ia, factor_ja
if ((i_start_at==0) .and. (i_final_at==0)) then
  d_n_neigh(:)=0
  d_kind_neigh(:,:)=0
  local_pow_so3_out(:,:)=0.d0
  local_pow_so3_deriv_out(:,:,:,:)=0.d0
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
local_pow_so3_out(:,:)=0.d0
local_pow_so3_deriv_out(:,:,:,:)=0.d0



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
     tmp_pow_so3_out(:)=0.d0
     tmp_pow_so3_deriv_out(:,:,:)=0.d0
     c(:,:)=(0.d0,0.d0)
     dc(:,:,:,:)=(0.d0,0.d0)
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
           i_count=0
           do p=1,n_rbf
              do l=0,l_max
                 i_count = i_count + 1
                 do m=-l,l
                    c(m,i_count) = c(m,i_count) + rbf_ji(p) * spherical_harm(l,m,dxp_ji(:))
                    if (desc_forces) dc(m,i_count,ia_n,1:3) = dc(m,i_count,ia_n,1:3) +  &
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
        do l=0,l_max
           i_count = i_count + 1

            do m=-l,l
               tmp_pow_so3_out(i_count) =  tmp_pow_so3_out(i_count) + real(dconjg(c(m,i_count)) * c(m,i_count), kind(0.d0))
               do ia=1,d_n_neigh(ja)
                 if (desc_forces) tmp_pow_so3_deriv_out(i_count,ia, 1:3) =  tmp_pow_so3_deriv_out(i_count,ia, 1:3) + 2.d0*real(dconjg(c(m,i_count))*dc(m,i_count,ia,1:3), kind(0.d0) )
              end do
            enddo !m

        enddo    !p
     enddo       !l

     if (desc_forces) then
      do ia=1,d_n_neigh(ja)
        local_pow_so3_deriv_out(1:pow_so3_dim,ja,  ia,1:3) =  tmp_pow_so3_deriv_out(1:pow_so3_dim,ia,1:3)*factor_ja
        local_pow_so3_deriv_out(1:pow_so3_dim,ja,   0,1:3)= local_pow_so3_deriv_out(1:pow_so3_dim,ja,0,1:3) - tmp_pow_so3_deriv_out(1:pow_so3_dim,ia,1:3)
      end do
     end if
     local_pow_so3_out(:,ja) = tmp_pow_so3_out(:)*factor_ja
  end do      !ja


 deallocate(xpnp)

return
end subroutine compute_pow_so3



subroutine init_pow_so3_rbf()
use ml_in_ndm_module, only: n_rbf, W_pow_so3, coeff_rbf, r_cut
use compute_afs_mod
implicit none
real(kind(0.d0)),dimension(n_rbf,n_rbf) :: S,V
real(kind(0.d0)),dimension(n_rbf) :: L
!local
integer :: p, q, nb, ilaenv,lwork

if (allocated(W_pow_so3)) deallocate(W_pow_so3) ; allocate(W_pow_so3(n_rbf,n_rbf))
if (allocated(coeff_rbf)) deallocate(coeff_rbf) ; allocate(coeff_rbf(n_rbf))


do q=1,n_rbf
   do p=1,n_rbf
      S(p,q)=dsqrt((2.d0*dble(p)+5.d0)*(2.d0*dble(q)+5.d0))/(dble(p+q)+5.d0)
   enddo
enddo

 nb = ilaenv( 1, 'DSYTRD', 'L', n_rbf, -1, -1, -1 )
 lwork=(nb+2)*n_rbf
 call diagsym(S,n_rbf,lwork,L)
 V(:,:)=0d0
 do p=1,n_rbf

    if (L(p)==0.d0) then
       V(p,p) = 0.d0
    else
       V(p,p) = L(p)/dabs(L(p)) * dabs(L(p))**(-0.5)
    end if

    coeff_rbf(p) = dsqrt((2.d0*p+5.d0)*r_cut**(-2.d0*p-5.d0))
 enddo

W_pow_so3(:,:)=matmul(S(:,:),matmul(V(:,:),transpose(S(:,:))))
!W_pow_so3(:,:)=matmul( matmul(S(:,:), V(:,:)),TRANSPOSE(S(:,:)))

return
end subroutine init_pow_so3_rbf
end module
