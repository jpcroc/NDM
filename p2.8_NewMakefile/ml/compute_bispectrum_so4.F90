module compute_bispectrum_so4_mod
        implicit none
        contains
!subroutine compute_bispectrum_so4(i_start_at,i_final_at,d_n_neigh, d_kind_neigh, local_bispectrum_so4_out,local_bispectrum_so4_deriv_out, iconf)
subroutine compute_bispectrum_so4(i_start_at,i_final_at,d_n_neigh, d_kind_neigh, iconf)

use T_kind_param_m, ONLY:  double
use gen_com_m, ONLY: imm, A2cm,lperiod,bg,at,indi2
use tab_imm_m, ONLY : iwmax2,xp
use angular_functions
use ml_in_ndm_module, ONLY: rangml,pi,r_cut,j_max,jj_max, imm_neigh, mconf,   &
                            debug, w2_rho,weighted, cg_vector, bisso4_dim, lbso4_diag, desc_forces, factor_weight_mass,bisso4_dim, linvisible
use derived_types, only : config_real, config_desc
use notperiod_mod
implicit none
integer, intent (in) :: i_start_at,i_final_at
integer, dimension(imm),intent(out)  :: d_n_neigh
integer, dimension(imm,imm_neigh), intent(out) :: d_kind_neigh
!double precision,dimension(bisso4_dim,imm),intent(out) :: config_desc(iconf)%energy
!double precision,dimension(bisso4_dim,imm, 0:imm_neigh, 3),intent(out) :: config_desc(iconf)%force
integer, optional :: iconf

logical :: small
double complex,dimension(-jj_max:jj_max,-jj_max:jj_max,0:jj_max) :: Umm
double complex,dimension(-jj_max:jj_max,-jj_max:jj_max,0:jj_max,3):: dUmm
! cmm derivatives for the componenets for 4D spherical functions.
! cmm (m1,m2,j), with m1,2=-j,j, j=0,2*j_max
! please note that for derivatives wer have the components cmm(m1,m2,j,ia_n,3) with ia_n  running from 0 to max of neighbours of a central atom.
! ia_n =  0 is for central atom
! ia_n != 0 is for any other atom not central
double complex, allocatable, dimension(:,:,:) :: cmm, cmm2
double complex, allocatable, dimension(:,:,:,:,:) :: dcmm, dcmm2




real(double), dimension(:,:), allocatable :: xpnp
real(double), dimension(3) :: dxp_ji, ds
integer :: j,iw,iw1,iw2

integer :: l,l1,l2,m1,m2,m1_1,m1_2
integer :: ia, ja, ix, i_bi, ia_n
double precision :: r2_ji,r_ji
double complex :: sum_bi, sum_bi_w, czero
double complex,dimension(3) :: dsum_bi, dsum_bi_w
double precision :: factor_ia, factor_ja, fcut, dfcut, fcut_w, dfcut_w, cg_local
integer :: l1_min, l1_max


if (allocated(cmm)) deallocate(cmm);    allocate(cmm(-jj_max:jj_max, &
                                                     -jj_max:jj_max, &
                                                      0:jj_max) )

if (weighted) then
    if (allocated(cmm2)) deallocate(cmm2);    allocate(cmm2(-jj_max:jj_max, &
                                                     -jj_max:jj_max, &
                                                      0:jj_max) )
end if


if (desc_forces) then
  if (allocated(dcmm)) deallocate(dcmm) ; allocate(dcmm(-jj_max:jj_max, &
                                                      -jj_max:jj_max, &
                                                       0:jj_max, &
                                                       0:imm_neigh,  &
                                                       1:3)  )
  if (weighted) then
    if (allocated(dcmm2)) deallocate(dcmm2) ; allocate(dcmm2(-jj_max:jj_max, &
                                                      -jj_max:jj_max, &
                                                       0:jj_max, &
                                                       0:imm_neigh,  &
                                                       1:3)  )
  end if

end if

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
czero=cmplx(0.d0,0.d0, kind=kind(1.d0))

ALLOCATE(xpnp(3,imm))
if (lperiod) then
   xpnp(:,:)=xp(:,:)
else
   call notperiod(xp,xpnp)
end if
!call cryst_to_cart (imm, xpnp, bg, -1)

d_n_neigh(:)=0
d_kind_neigh(:,:)=0
config_desc(iconf)%energy(:,:)=0.d0
if (desc_forces) config_desc(iconf)%force(:,:,:,:)=0.d0

!if (allocated(m_typ)) deallocate(m_typ); allocate(m_typ(ntyp,ntyp))
!call mconf(m_typ(:,:))

if (i_start_at==1)  iw2=0
if (i_start_at > 1) iw2=iwmax2(i_start_at-1)

do ja=i_start_at,i_final_at

    !debug if (rangml==1)  write(6,*) config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ja)), config_real(iconf)%Z_per_type(config_real(iconf)%itype(ja))
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

    cmm(:,:,:) = 0.d0
    if (weighted) cmm2(:,:,:)=0.d0
    do j = 0,jj_max
      do m1=-j,j,2
        cmm(m1,m1,j) = 1.d0
        if (weighted) cmm2(m1,m1,j) = 1.d0
      end do
    end do
    if (desc_forces) then
      dcmm(:,:,:,:,:) = 0.d0
      if (weighted)  dcmm2(:,:,:,:,:) = 0.d0
    end if
    ia_n=0

    if (debug) then
      if (mod(ja-1,10)==0) then
        if(rangml==0) write(6,'("in bso4 i_start_at i_final_at  ja:  ",3i9)') i_start_at, i_final_at , ja
      end if
    end if


    if (weighted) then
      factor_ja=config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))!/factor_weight_mass
      !write(44,*) 'factor_ja', factor_ja
    else
     factor_ja=1.d0
    endif




    do iw=iw1,iw2
        !begin small box or not 2/
        if (small) then
           ia = config_real(iconf)%kind_neigh(ja,iw)
           !write(*,*) 'debug kind_neigh', ia
        else
           ia=indi2(iw)
           if (ja==ia) cycle
        end if
        !end   small box or not 2/

        if (weighted) then

            if (linvisible) then
              if ( config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ia)) ) cycle
            end if

            factor_ia=config_real(iconf)%weight_per_type(config_real(iconf)%itype(ia))!/factor_weight_mass
            !write(44,*) 'factor_ia', factor_ia
        else
           factor_ia=1.d0
        endif

        if (small) then
           r_ji = config_real(iconf)%r_ij(ja,iw)
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

        if (r_ji .gt. r_cut) cycle
        ia_n = ia_n + 1
        fcut=0.5d0*(cos(pi*r_ji/r_cut)+1d0)
        dfcut= -0.5d0*pi/r_cut*sin(pi*r_ji/r_cut)
        fcut_w=fcut*factor_ia
        dfcut_w=dfcut_w*factor_ia
        call spherical_4d(dxp_ji, r_ji, Umm,dUmm)
        !debug write(*,'("ddd  ", 4f15.7)') dxp_ji(:),r_ji
        !debug write(*,'("ppp  ", 2f15.7)') Umm(-1,-1,1)
        cmm(:,:,:)= cmm(:,:,:)+Umm(:,:,:)*fcut
        if (weighted) cmm2(:,:,:)= cmm2(:,:,:)+Umm(:,:,:)*fcut_w
        if (desc_forces) then
          do ix=1,3
              dcmm(:,:,:,ia_n,ix) =  (dfcut*dxp_ji(ix)*Umm(:,:,:)/r_ji + fcut*dUmm(:,:,:,ix))
            if (weighted) dcmm2(:,:,:,ia_n,ix) =  (dfcut_w*dxp_ji(ix)*Umm(:,:,:)/r_ji + fcut_w*dUmm(:,:,:,ix))
          end do
        end if
        d_kind_neigh(ja,ia_n)=ia
    end do  !iw
    d_n_neigh(ja)=ia_n
    if (desc_forces) then
      dcmm(:,:,:,0,:) = - sum(dcmm(:,:,:,:,:), dim=4)
      if (weighted) dcmm2(:,:,:,0,:) = - sum(dcmm2(:,:,:,:,:), dim=4)
    end if

    do ia=1,ia_n
      i_bi=0
      do l1=0,jj_max
              if (lbso4_diag) then
                   !GaborLike l2=l1
                   l1_min=l1
                   l1_max=l1
                else
                   !Tlike do l2=0,l1
                   l1_min=0
                   l1_max=l1
              end if
              do l2=l1_min, l1_max
              do l=abs(l1-l2),min(jj_max,l1+l2)
                 if (mod(l1+l2+l,2)==1) cycle
                 if (.not.(lbso4_diag)) then
                      if (l < l1) cycle  ! this comes from Aidan Thompson and SNAP
                 end if
                 i_bi=i_bi+1

                 do m1=-l,l,2
                    do m2=-l,l,2

                        sum_bi=czero
                        sum_bi_w=czero
                        dsum_bi(1:3)=czero
                        dsum_bi_w(1:3)=czero

                        do m1_1=max(-l1,m1-l2),min(l1,m1+l2),2
                          do m1_2=max(-l1,m2-l2),min(l1,m2+l2),2
                            cg_local=cg_vector(l1,m1_1,l2,m1-m1_1,l,m1)*cg_vector(l1,m1_2,l2,m2-m1_2,l,m2)

                            sum_bi = sum_bi +  cg_local*cmm(m1_1,m1_2,l1)*cmm(m1-m1_1, m2-m1_2, l2)
                            if (weighted) then
                              sum_bi_w = sum_bi_w + cg_local*cmm2(m1_1,m1_2,l1)*cmm2(m1-m1_1, m2-m1_2, l2)
                            end if

                            if (desc_forces) then
                              dsum_bi(1:3) = dsum_bi(1:3) + cg_local*&
                                                          (dcmm(m1_1,m1_2,l1,ia, 1:3)*cmm(m1-m1_1,m2-m1_2,l2)+cmm(m1_1,m1_2,l1)*dcmm(m1-m1_1, m2-m1_2,l2,ia, 1:3))

                              if (weighted) then
                                dsum_bi_w(1:3) = dsum_bi_w(1:3) + cg_local*&
                                                          (dcmm2(m1_1,m1_2,l1,ia, 1:3)*cmm2(m1-m1_1,m2-m1_2,l2)+cmm2(m1_1,m1_2,l1)*dcmm2(m1-m1_1, m2-m1_2,l2,ia, 1:3))
                              end if
                            end if


                          enddo
                        enddo

                                     config_desc(iconf)%energy(i_bi,ja) = config_desc(iconf)%energy(i_bi,ja) + real ( conjg(cmm(m1,m2,l)) * sum_bi, kind=kind(0.d0))
                       if (weighted) config_desc(iconf)%energy(i_bi+bisso4_dim,ja) = config_desc(iconf)%energy(i_bi+bisso4_dim,ja) + real ( conjg(cmm2(m1,m2,l)) * sum_bi_w, kind=kind(0.d0))
                       !if (weighted) config_desc(iconf)%energy(i_bi+bisso4_dim,ja) = config_desc(iconf)%energy(i_bi,ja)

                       if (desc_forces) then
                        config_desc(iconf)%force(i_bi,ja, ia,1:3) = config_desc(iconf)%force(i_bi,ja, ia,1:3) - &
                                        real ( ( conjg(dcmm(m1,m2,l,ia, 1:3)) * sum_bi + conjg(cmm(m1,m2,l)) * dsum_bi(1:3) ), kind=kind(0.d0))
                        if (weighted) then
                          config_desc(iconf)%force(i_bi+bisso4_dim,ja, ia,1:3) = config_desc(iconf)%force(i_bi+bisso4_dim,ja, ia,1:3) - &
                                        real ( ( conjg(dcmm2(m1,m2,l,ia, 1:3)) * sum_bi_w + conjg(cmm2(m1,m2,l)) * dsum_bi_w(1:3) ), kind=kind(0.d0))

                        end if
                       end if
                    enddo  !m2
                 enddo     !m1
              end do  !l
              end do
            !Tlike enddo
           enddo !l1

       if (desc_forces) config_desc(iconf)%force(:,ja, 0,:) =  config_desc(iconf)%force(:,ja, 0,:) -  config_desc(iconf)%force(:,ja, ia,:)
    enddo  !ia from ia_m
    !NON NORMALIZED VERSION
    config_desc(iconf)%energy(:,ja) = config_desc(iconf)%energy(:,ja)/dble(ia_n)
    if (weighted) config_desc(iconf)%energy(bisso4_dim+1:2*bisso4_dim,ja) = config_desc(iconf)%energy(bisso4_dim+1:2*bisso4_dim,ja) *factor_ja
    if (desc_forces) then
      if (weighted) config_desc(iconf)%force(bisso4_dim+1:2*bisso4_dim,ja, :,:) =  config_desc(iconf)%force(bisso4_dim+1:2*bisso4_dim,ja, :,:)*factor_ja!/dble(ia_n)
    end if
  enddo !ja
 deallocate(xpnp)

return
end subroutine compute_bispectrum_so4

subroutine gen_dimension_for_bispectrum_so4_all
use ml_in_ndm_module, ONLY:  jj_max, bisso4_l, bisso4_l1, bisso4_l2, bisso4_dim
implicit none
integer :: l1, l2, l
integer :: i_bi

i_bi=0
do l1=0,jj_max
  !l2=l1 ! following Gabor  only the diagonal elements are important
  do l2=0, l1  ! in the end we want only the componenets with l1 <= l2 <= l
     do l=abs(l1-l2),min(jj_max,l1+l2)
       if (mod(l1+l2+l,2)==1) cycle ! assure invariance par réflexion
       if (l < l1) cycle  ! this comes from Thompson
       i_bi = i_bi+1
     enddo
  end do
enddo

bisso4_dim=i_bi
if (allocated(bisso4_l  )) deallocate(bisso4_l  ) ; allocate (bisso4_l(bisso4_dim))
if (allocated(bisso4_l1 )) deallocate(bisso4_l1 ) ; allocate (bisso4_l1(bisso4_dim))
if (allocated(bisso4_l2 )) deallocate(bisso4_l2 ) ; allocate (bisso4_l2(bisso4_dim))

i_bi=0
do l1=0,jj_max
     !l2=l1 ! following Gabor only diagonal elements are importants
     do l2=0, l1
        do l=abs(l1-l2),min(jj_max,l1+l2)
          if (mod(l1+l2+l,2)==1) cycle ! assure invariance par réflexion
          if (l < l1) cycle           ! in the end we want only the componenets with l1 <= l2 <= l
          i_bi = i_bi+1
          bisso4_l1 (i_bi) = l1
          bisso4_l2 (i_bi) = l2
          bisso4_l  (i_bi) = l
        end do
     end do
enddo

return
end subroutine gen_dimension_for_bispectrum_so4_all


subroutine  gen_dimension_for_bispectrum_so4_diagonal
use ml_in_ndm_module, ONLY: jj_max, bisso4_l, bisso4_l1, bisso4_l2, bisso4_dim
implicit none
integer :: l1, l2, l
integer :: i_bi

i_bi=0
do l1=0,jj_max
     l2=l1 ! following Gabor  only the diagonal elements are important
     do l=abs(l1-l2),min(jj_max,l1+l2)
       if (mod(l1+l2+l,2)==1) cycle ! assure invariance par réflexion
       i_bi = i_bi+1
     enddo
enddo

bisso4_dim=i_bi
if (allocated(bisso4_l  )) deallocate(bisso4_l  ) ; allocate (bisso4_l(bisso4_dim))
if (allocated(bisso4_l1 )) deallocate(bisso4_l1 ) ; allocate (bisso4_l1(bisso4_dim))
if (allocated(bisso4_l2 )) deallocate(bisso4_l2 ) ; allocate (bisso4_l2(bisso4_dim))

i_bi=0
do l1=0,jj_max
     l2=l1 ! following Gabor only diagonal elements are importants
        do l=abs(l1-l2),min(jj_max,l1+l2)
          if (mod(l1+l2+l,2)==1) cycle ! assure invariance par réflexion
          i_bi = i_bi+1
          bisso4_l1 (i_bi) = l1
          bisso4_l2 (i_bi) = l2
          bisso4_l  (i_bi) = l
        end do
enddo

return
end subroutine gen_dimension_for_bispectrum_so4_diagonal
end module
