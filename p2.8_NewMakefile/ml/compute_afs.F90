!subroutine compute_afs(i_start_at,i_final_at, d_n_neigh, d_kind_neigh,  local_afs_out,config_desc(iconf)%force, iconf)
subroutine compute_afs(i_start_at,i_final_at, d_n_neigh, d_kind_neigh, iconf)

USE T_kind_param_m, ONLY:  double
use gen_com_m, ONLY: imm,A2cm,lperiod,bg,at,indi2
use tab_imm_m, ONLY : iwmax2,xp
use ml_in_ndm_module, ONLY: n_rbf,n_cheb,afs_dim, tconf,write_desc,w3_rho,weighted,massat, &
                            imm_neigh, r_cut, n_rbf, W_afs, desc_forces, factor_weight_mass, coeff_rbf_afs, linvisible
use derived_types, only : config_real, config_desc

 implicit none
  integer, intent (in) :: i_start_at,i_final_at
integer, dimension(imm),intent(out)  :: d_n_neigh
integer, dimension(imm,imm_neigh), intent(out) :: d_kind_neigh
!double precision,dimension(afs_dim, imm),intent(out)   :: config_desc(iconf)%energy
!double precision,dimension(afs_dim, imm,0:imm_neigh, 3),intent(out) :: config_desc(iconf)%force
integer, optional :: iconf

logical :: small
real(double), dimension(:,:), allocatable :: xpnp
real(double), dimension(3) :: dxp_ji,dxp_jk, dxp_ik, ds
real(double), dimension(3) :: cdxp_ji,cdxp_jk, cdxp_ik
real(double), dimension(n_rbf) :: phi_ji, dphi_ji, phi_jk, dphi_jk, rbf_ji, rbf_jk, drbf_ji, drbf_jk
integer :: ia,ja,ka, ia_n, ka_n, iw,iw1,iw2, i_desc

integer :: iz,p1,p2
double precision :: r2_ji,r2_jk,r_ji,r_jk, r2_ik, r_ik, cos_jik
double precision,dimension(0:n_cheb) :: T_cos_jik,U_cos_jik
double precision,dimension(3) :: dcos_jik_ji,dcos_jik_jk,dT_cos_jik_ji,dT_cos_jik_jk
double precision :: factor_ia, factor_ja, factor_ka, term_local, term_local_3d(3)

namelist /input_ml/ write_desc,weighted,massat

if ((i_start_at==0) .and. (i_final_at==0)) then
  d_n_neigh(:)=0
  d_kind_neigh(:,:)=0
  config_desc(iconf)%energy(:,:)=0.d0
!  config_desc(iconf)%force(:,:,:,:)=0.d0
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
!   call cryst_to_cart (imm, xpnp, bg, -1)


d_n_neigh(:)=0
d_kind_neigh(:,:)=0
config_desc(iconf)%energy(:,:)=0.d0
if (desc_forces) config_desc(iconf)%force(:,:,:,:)=0.d0

!call compute_rbf(i_start_at,i_final_at,xpnp,rbf,drbf)

!  if (allocated(t_typ)) deallocate(t_typ)
!  allocate(t_typ(ntyp,ntyp,ntyp))
!  call tconf(t_typ(:,:,:))

if (i_start_at==1)  iw2=0
if (i_start_at > 1) iw2=iwmax2(i_start_at-1)

!debug write(*,*) i_start_at, i_final_at, rangml
do ja=i_start_at,i_final_at
    if (linvisible.and.weighted) then
        if( config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ja)) ) cycle
    end if
    !begin small box or not 1/
    if (small) then
      iw1=1
      iw2=config_real(iconf)%n_neigh(ja)
    else
      iw1=iw2+1
      iw2=iwmax2(ja)
    end if
    !end   small box or not 1/

    if (weighted) then
        !factor=w2_rho(massat(ityp(ja)),massat(ityp(ia)))
        !factor_ja=massat(ityp(ja))
        factor_ja=config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))!/factor_weight_mass
    else
        factor_ja=1.d0
    endif
    ia_n = 0
     do iw=iw1,iw2
        !begin small box or not 2/
        if (small) then
           ia = config_real(iconf)%kind_neigh(ja,iw)
        else
           ia=indi2(iw)
           if (ja==ia) cycle
        end if


        if (linvisible.and.weighted) then
            if( config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ia)) ) cycle
        end if

        !end   small box or not 2/

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
           dxp_ji = MatMul(at,ds)/A2cm
           r2_ji = Sum( dxp_ji(1:3)**2 )
           r_ji = dsqrt(r2_ji)
        end if
        phi_ji(:)=0.d0
        if (r_ji >= r_cut) cycle
        ia_n = ia_n + 1

        ka_n=0
        do iz=iw1,iw2

           !begin small box or not 2/
           if (small) then
              ka = config_real(iconf)%kind_neigh(ja,iz)
           else
              ka=indi2(iz)
              if ( (ja==ka).or.(ia==ka)) cycle
           end if
           !end   small box or not 2/



           if (linvisible.and.weighted) then
               if( config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ka)) ) cycle
           end if

           if (weighted) then
              !factor_ka=massat(ityp(ka))
              factor_ka=config_real(iconf)%weight_per_type(config_real(iconf)%itype(ka))!/factor_weight_mass
           else
              factor_ka=1.d0
           endif

           if (small) then
              r_jk = config_real(iconf)%r_ij(ja,iz)
              r2_jk=r_jk**2
              dxp_jk(:) = config_real(iconf)%u_ij(ja,iz,:)
              dxp_ik(:) = dxp_jk(:) - dxp_ji(:)
              r2_ik = Sum(dxp_ik(:)**2)
              r_ik = dsqrt(r2_ik)
           else
              dxp_jk(1:3) = xpnp(1:3,ka) - xpnp(1:3,ja)
              ds=MatMul(dxp_jk,bg)
              WHERE ( (ds.GT.0.5d0).OR.(ds.LT.-0.5d0) )
                 ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
              END WHERE
              dxp_jk = MatMul(at,ds)/A2cm
              r2_jk = Sum( dxp_jk(1:3)**2 )
              r_jk = dsqrt(r2_jk)

              dxp_ik(1:3) = xpnp(1:3,ka) - xpnp(1:3,ia)
              ds=MatMul(dxp_ik,bg)
              WHERE ( (ds.GT.0.5d0).OR.(ds.LT.-0.5d0) )
                 ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
              END WHERE
              dxp_ik = MatMul(at,ds)/A2cm
              r2_ik = Sum( dxp_ik(1:3)**2 )
              r_ik = dsqrt(r2_ik)
           end if

           if (r_jk >= r_cut) cycle
           ka_n = ka_n + 1
           !why_not if (r_ik >= r_cut) cycle
           if (r_ik <= 1d-20) cycle

           do p1=1,n_rbf
            !phi_ji(p1) = dsqrt((2.d0*p1+5.d0)*r_cut**(-2.d0*p1-5.d0))*(r_cut-r_ji)**(p1+2.d0)
            phi_ji(p1) = coeff_rbf_afs(p1)*(r_cut-r_ji)**(p1+2.d0)
            !if (desc_forces) dphi_ji(p1)=-dsqrt((2.d0*p1+5.d0)*r_cut**(-2.d0*p1-5.d0))*(p1+2.d0)*(r_cut-r_ji)**(p1+1.d0)
            if (desc_forces) dphi_ji(p1)=-coeff_rbf_afs(p1)*(p1+2.d0)*(r_cut-r_ji)**(p1+1.d0)

            !phi_jk(p1) = dsqrt((2.d0*p1+5.d0)*r_cut**(-2.d0*p1-5.d0))*(r_cut-r_jk)**(p1+2.d0)
            phi_jk(p1) = coeff_rbf_afs(p1)*(r_cut-r_jk)**(p1+2.d0)
            !if (desc_forces) dphi_jk(p1)=-dsqrt((2.d0*p1+5.d0)*r_cut**(-2.d0*p1-5.d0))*(p1+2.d0)*(r_cut-r_jk)**(p1+1.d0)
            if (desc_forces) dphi_jk(p1)=-coeff_rbf_afs(p1)*(p1+2.d0)*(r_cut-r_jk)**(p1+1.d0)
          enddo


           rbf_ji(:) =matmul(W_afs(:,:) , phi_ji(:))
           if (desc_forces) drbf_ji(:)=matmul(W_afs(:,:) ,dphi_ji(:))

           rbf_jk(:) =matmul(W_afs(:,:) , phi_jk(:))
           if (desc_forces) drbf_jk(:)=matmul(W_afs(:,:) ,dphi_jk(:))

           cdxp_ji(:)=dxp_ji(:)/r_ji
           cdxp_jk(:)=dxp_jk(:)/r_jk
           cdxp_ik(:)=dxp_ik(:)/r_ik

           cos_jik=dot_product(cdxp_ji(1:3),cdxp_jk(1:3))

           if (desc_forces) dcos_jik_ji(1:3)=(cdxp_jk(1:3) - cos_jik*cdxp_ji(1:3))/r_ji
           if (desc_forces) dcos_jik_jk(1:3)=(cdxp_ji(1:3) - cos_jik*cdxp_jk(1:3))/r_jk

           T_cos_jik(0)=1.d0
           T_cos_jik(1)=cos_jik
           U_cos_jik(0)=1.d0
           U_cos_jik(1)=2.d0*cos_jik
           if (n_cheb.ge.2) then
              do p2=2,n_cheb
                 T_cos_jik(p2)=2d0*cos_jik*T_cos_jik(p2-1)-T_cos_jik(p2-2)
                 U_cos_jik(p2)=2d0*cos_jik*U_cos_jik(p2-1)-U_cos_jik(p2-2)
              enddo
           endif

           i_desc=0
           do p1=1,n_rbf
           do p2=0,n_cheb

              if (p2==0) then
                 dT_cos_jik_ji(1:3)=0d0
                 dT_cos_jik_jk(1:3)=0d0
              else
                 dT_cos_jik_ji(1:3)=p2*U_cos_jik(p2-1)*dcos_jik_ji(1:3)
                 dT_cos_jik_jk(1:3)=p2*U_cos_jik(p2-1)*dcos_jik_jk(1:3)
              endif

                 i_desc = i_desc + 1
                 !oldWesley config_desc(iconf)%energy(i_desc,ja) = config_desc(iconf)%energy(i_desc,ja) + factor_ia*factor_ka* rbf_ji(p1,ia,ja) * rbf(p1,ka,ja) * T_cos_jik(p2)
                 term_local=rbf_ji(p1) * rbf_jk(p1) * T_cos_jik(p2)/2.d0
                 config_desc(iconf)%energy(i_desc,ja) = config_desc(iconf)%energy(i_desc,ja) + factor_ia*factor_ka* term_local
                 if (weighted) config_desc(iconf)%energy(i_desc+afs_dim,ja) = config_desc(iconf)%energy(i_desc+afs_dim,ja) +  term_local
                 !oldWesley config_desc(iconf)%force(i_desc,ja,ia_n,1:3) = config_desc(iconf)%force(i_desc, ja, ia_n, 1:3) - factor_ia* factor_ka* ( rbf(p1,ia,ja)*rbf(p1,ka,ja)*(dT_cos_jik_ji(1:3)+dT_cos_jik_jk(1:3)) + &
                 !oldWesley                                        (drbf(p1,ia,ja)*cdxp_ji(1:3)*rbf(p1,ka,ja)+rbf(p1,ia,ja)*drbf(p1,ka,ja)*cdxp_jk(1:3)) * T_cos_jik(p2) )

                 !config_desc(iconf)%force(i_desc,ja,ia_n,1:3) = config_desc(iconf)%force(i_desc, ja, ia_n, 1:3) - factor_ia* factor_ka* ( rbf_ji(p1)*rbf_jk(p1)*(dT_cos_jik_ji(1:3)+dT_cos_jik_jk(1:3)) + &
                 !                                          (drbf_ji(p1)*cdxp_ji(1:3)*rbf_jk(p1)+rbf_ji(p1)*drbf_jk(p1)*cdxp_jk(1:3)) * T_cos_jik(p2) )
                 if (desc_forces) then
                        term_local_3d(1:3)= rbf_ji(p1)*rbf_jk(p1)*dT_cos_jik_ji(1:3) + drbf_ji(p1)*cdxp_ji(1:3)*rbf_jk(p1)*T_cos_jik(p2)
                        config_desc(iconf)%force(i_desc,ja,ia_n,1:3) = config_desc(iconf)%force(i_desc, ja, ia_n, 1:3) + factor_ia* factor_ka* term_local_3d(1:3)
                        if (weighted) then
                           config_desc(iconf)%force(i_desc+afs_dim,ja,ia_n,1:3) = config_desc(iconf)%force(i_desc+afs_dim, ja, ia_n, 1:3) + term_local_3d(1:3)
                           !3d config_desc(iconf)%force(i_desc+2*afs_dim,ja,ia_n,1:3) = config_desc(iconf)%force(i_desc+2*afs_dim, ja, ia_n, 1:3) + sqrt(factor_ia)*sqrt(factor_ka)*term_local_3d(1:3)
                        end if
                 end if

                !config_desc(iconf)%force(i_desc,ja,0,1:3) = config_desc(iconf)%force(i_desc, ja, ia_n, 1:3) - factor_ia* factor_ka* 0.5d0* ( &
                !                                           rbf_ji(p1)*rbf_jk(p1)*(dT_cos_jik_ji(1:3) + dT_cos_jik_jk(1:3)) + &
                !                                          drbf_ji(p1)*cdxp_ji(1:3)*rbf_jk(p1)*T_cos_jik(p2) + rbf_ji(p1)*drbf_jk(p1)*cdxp_jk(1:3)* T_cos_jik(p2) )
              enddo
           enddo
        end do ! ka of the  neigh of ja
        d_kind_neigh(ja,ia_n)=ia
        if (desc_forces) config_desc(iconf)%force(:,ja, 0,:) =  config_desc(iconf)%force(:,ja, 0,:) -  config_desc(iconf)%force(:,ja, ia_n,:)
     end do  ! ia_n loop over neigh of ja
  if (weighted) then
    config_desc(iconf)%energy(1:afs_dim,ja) = config_desc(iconf)%energy(1:afs_dim,ja)*factor_ja
    !3d config_desc(iconf)%energy(2*afs_dim+1:3*afs_dim,ja) = config_desc(iconf)%energy(2*afs_dim+1:3*afs_dim,ja)*sqrt(factor_ja)
  end if
  if (desc_forces)  then
    if (weighted) then
      config_desc(iconf)%force(1:afs_dim,ja, :,:) = config_desc(iconf)%force(1:afs_dim,ja, :,:)*factor_ja
      !3d config_desc(iconf)%force(2*afs_dim+1:3*afs_dim,ja, :,:) = config_desc(iconf)%force(2*afs_dim+1:3*afs_dim,ja, :,:)*sqrt(factor_ja)
    end if
  end if
  d_n_neigh(ja)=ia_n
  end do   !ja main loop


 deallocate(xpnp)

return
end subroutine compute_afs

subroutine pre_compute_neighbours_descriptors(i_start_at,i_final_at, d_n_neigh, d_kind_neigh,  iconf)

USE T_kind_param_m, ONLY:  double
use gen_com_m, ONLY: imm,A2cm,lperiod,bg,at,indi2
use tab_imm_m, ONLY : iwmax2,xp
use ml_in_ndm_module, ONLY: n_rbf,n_cheb,afs_dim, tconf,write_desc,w3_rho,weighted,massat, &
                            imm_neigh, r_cut, n_rbf, W_afs, desc_forces, factor_weight_mass, coeff_rbf_afs
use derived_types, only : config_real

implicit none
integer, intent (in) :: i_start_at,i_final_at
integer, dimension(imm),intent(out)  :: d_n_neigh
integer, dimension(imm,imm_neigh), intent(out) :: d_kind_neigh
integer  :: iconf

logical :: small
real(double), dimension(:,:), allocatable :: xpnp
real(double), dimension(3) :: dxp_ji,dxp_jk, dxp_ik, ds
real(double), dimension(3) :: cdxp_ji,cdxp_jk, cdxp_ik
real(double), dimension(n_rbf) :: phi_ji, dphi_ji, phi_jk, dphi_jk, rbf_ji, rbf_jk, drbf_ji, drbf_jk
integer :: ia,ja,ka, ia_n, ka_n, iw,iw1,iw2, i_desc

integer :: iz,p1,p2
double precision :: r2_ji,r2_jk,r_ji,r_jk, r2_ik, r_ik, cos_jik
double precision,dimension(0:n_cheb) :: T_cos_jik,U_cos_jik
double precision,dimension(3) :: dcos_jik_ji,dcos_jik_jk,dT_cos_jik_ji,dT_cos_jik_jk
double precision :: factor_ia, factor_ja, factor_ka

namelist /input_ml/ write_desc,weighted,massat

if ((i_start_at==0) .and. (i_final_at==0)) then
  d_n_neigh(:)=0
  d_kind_neigh(:,:)=0
  return
end if

small = config_real(iconf)%small

ALLOCATE(xpnp(3,imm))
if (lperiod) then
  xpnp(:,:)=xp(:,:)
else
   call notperiod(xp,xpnp)
end if
!   call cryst_to_cart (imm, xpnp, bg, -1)


d_n_neigh(:)=0
d_kind_neigh(:,:)=0


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

    ia_n = 0
     do iw=iw1,iw2
        !begin small box or not 2/
        if (small) then
           ia = config_real(iconf)%kind_neigh(ja,iw)
        else
           ia=indi2(iw)
           if (ja==ia) cycle
        end if
        !end   small box or not 2/

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
           dxp_ji = MatMul(at,ds)/A2cm
           r2_ji = Sum( dxp_ji(1:3)**2 )
           r_ji = dsqrt(r2_ji)
        end if
        phi_ji(:)=0.d0
        if (r_ji >= r_cut) cycle
        ia_n = ia_n + 1
        d_kind_neigh(ja,ia_n)=ia
     end do  ! ia_n loop over neigh of ja

  d_n_neigh(ja)=ia_n
  end do   !ja main loop

 deallocate(xpnp)

return
end subroutine pre_compute_neighbours_descriptors


subroutine init_afs_rbf()
use ml_in_ndm_module, only: afs_dim, n_rbf, n_cheb, W_afs, coeff_rbf_afs, r_cut
implicit none
real(kind(0.d0)),dimension(n_rbf,n_rbf) :: S,V
real(kind(0.d0)),dimension(n_rbf) :: L
!local
integer :: p, q, nb, ilaenv,lwork

if (allocated(W_afs)) deallocate(W_afs) ; allocate(W_afs(n_rbf,n_rbf))
if (allocated(coeff_rbf_afs)) deallocate(coeff_rbf_afs) ; allocate(coeff_rbf_afs(n_rbf))

afs_dim = n_rbf*(n_cheb+1)

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

    coeff_rbf_afs(p) = dsqrt((2.d0*p+5.d0)*r_cut**(-2.d0*p-5.d0))
 enddo

W_afs(:,:)=matmul(S(:,:),matmul(V(:,:),transpose(S(:,:))))
!W_pow_so3(:,:)=matmul( matmul(S(:,:), V(:,:)),TRANSPOSE(S(:,:)))

return
end subroutine init_afs_rbf


subroutine diagsym(A,n_rbf,lwork,L)
! input - the symmetric matrix A
! output - the ortho-normalized vectors that diaonalize A and the eigenvalues L
use ml_in_ndm_module, only: rangml
implicit none
integer,intent(in) :: n_rbf,lwork
double precision,dimension(n_rbf,n_rbf),intent(inout) :: A
double precision,dimension(n_rbf),intent(out) :: L
double precision,dimension(lwork) :: work
integer :: info

call dsyev('V','L',n_rbf,A,n_rbf,L,work,lwork,info)

if (rangml==0) then
  if (.not.(info == 0)) then
   write(6,*) 'ML: WARNING the nomalization of the overlap matrix is WRONG in diagsym compute_afs.F90'
  end if
end if

return
end subroutine diagsym



subroutine reverse_index_of_ghost_atoms(i_start_at, i_final_at, iconf)
use derived_types, only: config_desc, config_real
implicit none
integer, intent(in) :: i_start_at, i_final_at, iconf
integer :: ik,igh, rang_ia, ia
if ( (i_start_at==0).and.(i_final_at==0) ) return

  do ik = i_start_at, i_final_at
    do igh=1,config_desc(iconf)%n_neigh_ghost(ik)
      ia=config_desc(iconf)%kind_neigh_ghost(ik,igh)
      rang_ia = config_desc(iconf)%kind_neigh_proc(ik,igh)
      !call MPI_SEND( data, count, MPI_DOUBLE_PRECISION, rang_ia, tag, MPI_COMM_WORLD, ierr )
    end do
  end do


return
end subroutine reverse_index_of_ghost_atoms



subroutine distribute_ghost_descritors(i_start_at, i_final_at, iconf)
use derived_types, only: config_desc, config_real
implicit none
integer, intent(in) :: i_start_at, i_final_at, iconf
integer :: ik,igh, rang_ia, ia
if ( (i_start_at==0).and.(i_final_at==0) ) return

  do ik = i_start_at, i_final_at
    do igh=1,config_desc(iconf)%n_neigh_ghost(ik)
      ia=config_desc(iconf)%kind_neigh_ghost(ik,igh)
      rang_ia = config_desc(iconf)%kind_neigh_proc(ik,igh)
      !call MPI_SEND( data, count, MPI_DOUBLE_PRECISION, rang_ia, tag, MPI_COMM_WORLD, ierr )
    end do
  end do


return
end subroutine distribute_ghost_descritors
