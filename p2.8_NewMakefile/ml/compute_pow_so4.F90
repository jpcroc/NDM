subroutine compute_pow_so4(i_start_at,i_final_at,d_n_neigh, d_kind_neigh, &
                           local_pow_so4_out, local_pow_so4_deriv_out, iconf)

use T_kind_param_m, ONLY:  double
use gen_com_m, ONLY: imm, A2cm,lperiod,bg,at,indi2
use tab_imm_m, ONLY : iwmax2,xp
!use var_pot, only : ntyp
use angular_functions!, only : Umm
use derived_types, only : config_real
use ml_in_ndm_module, ONLY: imm_neigh, one_pi, rangml, r_cut,j_max,jj_max, weighted, factor_weight_mass, linvisible

implicit none
integer, intent (in) :: i_start_at,i_final_at
integer, dimension(imm),intent(out)  :: d_n_neigh
integer, dimension(imm,imm_neigh), intent(out) :: d_kind_neigh
double precision,dimension(0:jj_max,imm),intent(out) :: local_pow_so4_out
double precision,dimension(0:jj_max,imm,0:imm_neigh,3), intent(out) :: local_pow_so4_deriv_out
integer, optional :: iconf

logical :: small
double complex,dimension(-jj_max:jj_max,-jj_max:jj_max,0:jj_max) :: Umm
double complex,dimension(-jj_max:jj_max,-jj_max:jj_max,0:jj_max,3):: dUmm
! cmm derivatives for the componenets for 4D spherical functions.
! cmm (m1,m2,j), with m1,2=-j,j, j=0,2*j_max
! please note that for derivatives wer have the components cmm(m1,m2,j,ia_n,3) with ia_n  running from 0 to max of neighbours of a central atom.
! ia_n =  0 is for central atom
! ia_n != 0 is for any other atom not central
double complex, allocatable, dimension(:,:,:) :: cmm
double complex, allocatable, dimension(:,:,:,:,:) :: dcmm


real(double), dimension(:,:), allocatable :: xpnp
real(double), dimension(3) :: dxp_ji, ds
integer :: ia,ja, ia_n, ix, j,iw,iw1,iw2

integer :: l,m1,m2
double precision :: r2_ji,r_ji
double precision :: factor_ia,factor_ja, fcut, dfcut

if (allocated(cmm)) deallocate(cmm);    allocate(cmm(-jj_max:jj_max, &
                                                     -jj_max:jj_max, &
                                                      0:jj_max) )
if (allocated(dcmm)) deallocate(dcmm) ; allocate(dcmm(-jj_max:jj_max, &
                                                      -jj_max:jj_max, &
                                                       0:jj_max, &
                                                       0:imm_neigh,  &
                                                       1:3)  )


if ((i_start_at==0) .and. (i_final_at==0)) then
  d_n_neigh(:)=0
  d_kind_neigh(:,:)=0
  local_pow_so4_out(:,:)=0.d0
  local_pow_so4_deriv_out(:,:,:,:)=0.d0
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
!call cryst_to_cart (imm, xpnp, bg, -1)

d_n_neigh(:)=0
d_kind_neigh(:,:)=0
local_pow_so4_out(:,:)=0.d0
local_pow_so4_deriv_out(:,:,:,:)=0.d0

!if (allocated(m_typ)) deallocate(m_typ); allocate(m_typ(ntyp,ntyp))
!call mconf(m_typ(:,:))

if (i_start_at==1)  iw2=0
if (i_start_at > 1) iw2=iwmax2(i_start_at-1)

do ja=i_start_at,i_final_at
    !begin small box or not 1/

     if (small) then
      iw1=1
      iw2=config_real(iconf)%n_neigh(ja)
      !debug write(*,'("rangml, debug neigh",4i5)') rangml, ja, iw1,  iw2
      !debug write(*,*) rangml, config_real(iconf)%filename
    else
      iw1=iw2+1
      iw2=iwmax2(ja)
      !write(*,*) 'debug neigh', iw1, iw2
    end if
    !end   small box or not 1/
    !set-up initialization of some variables for descriptors ...


    if (linvisible.and.weighted) then
        if( config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ja)) ) cycle
    end if

    if (weighted) then
        factor_ja=config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))!/factor_weight_mass
    else
        factor_ja=1.d0
    endif


    cmm(:,:,:) = 0.d0
    do j = 0,jj_max
      do m1=-j,j,2
        cmm(m1,m1,j) = 1.d0
      end do
    end do
    dcmm(:,:,:,:,:) = 0.d0
    ia_n=0
    !end set-up of some variables for initialization of descriptors variables
    !write(*,*) 'debuj ja',ja
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

      if (linvisible.and.weighted) then
        if( config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ia)) ) cycle
      end if
      if (weighted) then
        factor_ia=config_real(iconf)%weight_per_type(config_real(iconf)%itype(ia))!/factor_weight_mass
      else
        factor_ia=1.d0
      endif

      if (small) then
        r_ji = config_real(iconf)%r_ij(ja,iw)
        dxp_ji(:) = config_real(iconf)%u_ij(ja,iw,:)
       else
        dxp_ji(1:3) = xpnp(1:3,ia) - xpnp(1:3,ja)
        ds(:) = MatMul( dxp_ji(:), bg(:,:) )
        WHERE ( (ds.GT.0.5d0).OR.(ds.LT.-0.5d0) )
           ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
        END WHERE
        dxp_ji(:) = MatMul(at(:,:),ds(:))
        dxp_ji(:) = dxp_ji(:)/A2cm
        r2_ji = Sum( dxp_ji(1:3)**2 )
        r_ji = dsqrt(r2_ji)
      end if

      if (r_ji >= r_cut) cycle
      ia_n = ia_n + 1
      fcut=0.5d0*(cos(one_pi*r_ji/r_cut)+1d0)*factor_ia
      dfcut= -0.5d0*one_pi/r_cut*sin(one_pi*r_ji/r_cut)*factor_ia
      call spherical_4d(dxp_ji, r_ji, Umm, dUmm)
      cmm(:,:,:)= cmm(:,:,:)+Umm(:,:,:)*fcut
      do ix=1,3
        dcmm(:,:,:,ia_n,ix) =  (dfcut*dxp_ji(ix)*Umm(:,:,:)/r_ji + fcut*dUmm(:,:,:,ix))
      end do
      d_kind_neigh(ja,ia_n)=ia
    end do
    d_n_neigh(ja)=ia_n
    dcmm(:,:,:,0,:) = - sum(dcmm(:,:,:,:,:), dim=4)

    do ia=0,ia_n
      do l=0,jj_max
        do m1=-l,l,2
          do m2=-l,l,2
            local_pow_so4_out(l,ja) = local_pow_so4_out(l,ja) + real(conjg(cmm(m2,m1,l)) * cmm(m2,m1,l))

            local_pow_so4_deriv_out(l,ja,ia,1:3) = local_pow_so4_deriv_out(l,ja,ia, 1:3) + &
                                              2.d0*conjg(dcmm(m2,m1,l,ia,1:3))*cmm(m2,m1,l)*factor_ja
          enddo
        enddo
      enddo
    end do !ia
    !write(*,'("dbg pow", 2i5, 8e20.10)') 1, 1, local_deriv(0:7,1,1,1)
    ! there is a normalisation with the number of neighbours, proposed by Gabor G/ Bartok 2013.
    local_pow_so4_out(0:jj_max,ja) = factor_ja*local_pow_so4_out(0:jj_max,ja)/dble(ia_n+1)
enddo !ja

 deallocate(xpnp)

return
end subroutine compute_pow_so4
