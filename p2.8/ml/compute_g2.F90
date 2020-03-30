module compute_g2_mod
        use notperiod_mod
        implicit none
        contains
subroutine compute_g2(i_start_at,i_final_at, local_d_n_neigh, local_d_kind_neigh,  local_g2_out,local_g2_deriv_out, iconf)

USE T_kind_param_m, ONLY:  double
use gen_com_m, ONLY: imm,A2cm,lperiod,bg,at,indi2
use tab_imm_m, ONLY : iwmax2,xp
use ml_in_ndm_module, ONLY: mconf,w2_rho,weighted, &
                            one_pi, imm_neigh, &
                            g2_eta, g2_rs, &
                            r_cut,  &
                            g2_dim, factor_weight_mass, desc_forces, weighted, linvisible
use derived_types, only: config_real
use notperiod_mod
implicit none


integer, intent (in) :: i_start_at,i_final_at
integer, dimension(imm),intent(out)  :: local_d_n_neigh
integer, dimension(imm,imm_neigh), intent(out) :: local_d_kind_neigh
double precision,dimension(g2_dim,imm),intent(out) :: local_g2_out
double precision,dimension(g2_dim,imm, 0:imm_neigh,3),intent(out) :: local_g2_deriv_out
integer, optional, intent(in) :: iconf

real(double), dimension(:,:), allocatable :: xpnp
real(double), dimension(3) :: dxp_ji, ds
integer :: iw,iw1,iw2

integer :: p
real(double) :: fcut,r_ji,dfcut
real(double) :: factor_ia, factor_ja, g2_func

integer :: ia,ja,ia_n
logical :: small

if ((i_start_at==0) .and. (i_final_at==0)) then
  local_d_kind_neigh(:,:)=0
  local_d_n_neigh(:)=0
  local_g2_out(:,:)=0.d0
  local_g2_deriv_out(:,:,:,:)=0.d0
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


local_d_kind_neigh(:,:)=0
local_d_n_neigh(:)=0
local_g2_out(:,:)=0d0
local_g2_deriv_out(:,:,:,:)=0d0


!  if (allocated(m_typ)) deallocate(m_typ); allocate(m_typ(ntyp,ntyp))
!  call mconf(m_typ(:,:))

  if (i_start_at==1)  iw2=0
  if (i_start_at > 1) iw2=iwmax2(i_start_at-1)

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

    if (weighted) then
      factor_ja=config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))!/factor_weight_mass
    else
      factor_ja=1.d0
    endif



     ia_n=0


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
           factor_ia=config_real(iconf)%weight_per_type(config_real(iconf)%itype(ia))!/factor_weight_mass
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
          r_ji = dsqrt( Sum( dxp_ji(1:3)**2 ) )
        end if

        if (r_ji >= r_cut) cycle
        ia_n = ia_n + 1

        fcut=0.5d0*(cos(one_pi*r_ji/r_cut)+1.d0)*factor_ia
        dfcut=-0.5d0*one_pi*sin(one_pi*r_ji/r_cut)/r_cut*factor_ia
        !debug write(*,*) fcut, dfcut, factor_ia
        do p=1,g2_dim
              g2_func = dexp(-g2_eta(p)*(r_ji-g2_rs(p))**2)
              local_g2_out(p,ja) = local_g2_out(p, ja) + fcut * g2_func
              if (desc_forces) then
                 local_g2_deriv_out(p,ja, ia_n, 1:3) = local_g2_deriv_out(p,ja, ia_n, 1:3) + (-2d0*g2_eta(p)*(r_ji-g2_rs(p))*fcut + dfcut)  * &
                                                          g2_func*dxp_ji(1:3)/r_ji

                local_g2_deriv_out(p,ja, 0, 1:3) = local_g2_deriv_out(p,ja, 0, 1:3) - local_g2_deriv_out(p,ja, ia_n, 1:3)
              end if
        enddo
        local_d_kind_neigh(ja,ia_n)=ia
     end do  !ia
  local_g2_out(:,:) = local_g2_out(:,:) * factor_ja!/dble(ia_n)
  if (desc_forces) local_g2_deriv_out(:,:,:,:) = local_g2_deriv_out(:,:,:,:) * factor_ja!/dble(ia_n)
  local_d_n_neigh(ja)=ia_n
  end do  !ja

 deallocate(xpnp)

return
end subroutine compute_g2
end module
