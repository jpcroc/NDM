module compute_g3_mod
        implicit none
        contains
subroutine compute_g3(i_start_at,i_final_at,  local_d_n_neigh, local_d_kind_neigh,   local_g3_out,local_g3_deriv_out, iconf)

USE T_kind_param_m, ONLY:  double
use gen_com_m, ONLY: imm,A2cm,lperiod,bg,at,indi2
use tab_imm_m, ONLY : iwmax2,xp
use ml_in_ndm_module, ONLY: tconf,w3_rho,weighted, &
                            one_pi,r_cut,g3_eta,g3_zeta,g3_lambda,g3_dim, imm_neigh, factor_weight_mass, linvisible

use derived_types, only: config_real
use notperiod_mod
use cryst_to_cart_mod
implicit none

integer, intent (in) :: i_start_at,i_final_at
integer, dimension(imm),intent(out)  :: local_d_n_neigh
integer, dimension(imm,imm_neigh), intent(out) :: local_d_kind_neigh
double precision,dimension(g3_dim, imm),intent(out) :: local_g3_out
double precision,dimension(g3_dim, imm,0:imm_neigh, 3),intent(out) :: local_g3_deriv_out
integer, optional, intent(in) :: iconf

real(double), dimension(:,:),allocatable :: xpnp
real(double), dimension(3) :: dxp_ji,dxp_jk,dxp_ik, ds
integer :: ia,ia_n,ja,ka,ka_n, iw,iw1,iw2

integer :: iz,p
double precision :: fcut_ji,fcut_jk,fcut_ik
double precision :: r_ji,r_jk,r_ik
double precision :: r2_ji,r2_jk,r2_ik

double precision :: cos_jik,ang,rad, dang
double precision,dimension(3) :: dfcut_ji   ,dfcut_jk   ,dfcut_ik, cdxp_ik, cdxp_ji, cdxp_jk
double precision,dimension(3) :: dcos_jik_ji,dcos_jik_jk
double precision,dimension(3) :: drad_ji    ,    drad_jk,drad_ik
double precision :: factor_ia, factor_ja, factor_ka

logical :: small

if ((i_start_at==0) .and. (i_final_at==0)) then
  local_d_kind_neigh(:,:)=0
  local_d_n_neigh(:)=0
  local_g3_out(:,:)=0.d0
  local_g3_deriv_out(:,:,:,:)=0.d0
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

local_d_kind_neigh(:,:)=0
local_d_n_neigh(:)=0
local_g3_out(:,:)=0.d0
local_g3_deriv_out(:,:,:,:)=0.d0

!  if (allocated(t_typ)) deallocate(t_typ) ; allocate(t_typ(ntyp,ntyp,ntyp))
!  call tconf(t_typ(:,:,:))

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


    if (linvisible.and.weighted) then
        if( config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ja)) ) cycle
    end if


     if (weighted) then
        !factor_ja=massat(ityp(ja))
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
           dxp_ji = MatMul(at,ds)/A2cm
           r2_ji = Sum( dxp_ji(1:3)**2 )
           r_ji = dsqrt(r2_ji)
        end if
        if (r_ji >= r_cut) cycle
        ia_n = ia_n + 1
        !
        cdxp_ji(1:3) = dxp_ji(1:3)/r_ji
        fcut_ji = 0.5d0*(cos(one_pi*r_ji/r_cut)+1.d0)*factor_ia*factor_ja
        dfcut_ji(1:3)=-0.5d0*one_pi*sin(one_pi*r_ji/r_cut)*cdxp_ji(1:3)/r_cut*factor_ia*factor_ja
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
           if (r_ik >= r_cut) cycle
           if (r_ik <= 1d-20) cycle

           cdxp_ik(1:3) = dxp_ik(1:3)/r_ik
           cdxp_jk(1:3) = dxp_jk(1:3)/r_jk

           fcut_jk = 0.5d0*(cos(one_pi*r_jk/r_cut)+1d0)*factor_ja*factor_ka
           fcut_ik = 0.5d0*(cos(one_pi*r_ik/r_cut)+1d0)*factor_ia*factor_ka

           dfcut_jk(1:3)=-0.5d0*one_pi*sin(one_pi*r_jk/r_cut)*cdxp_jk(1:3)/r_cut*factor_ja*factor_ka
           dfcut_ik(1:3)=-0.5d0*one_pi*sin(one_pi*r_ik/r_cut)*cdxp_ik(1:3)/r_cut*factor_ia*factor_ka

           cos_jik=dot_product(dxp_ji(1:3)/r_ji,dxp_jk(1:3)/r_jk)
           dcos_jik_ji(1:3)=(cdxp_jk(1:3)-cos_jik*cdxp_ji(1:3))/r_ji
           dcos_jik_jk(1:3)=(cdxp_ji(1:3)-cos_jik*cdxp_jk(1:3))/r_jk

           do p=1, g3_dim
                    ang = (1.d0+g3_lambda(p)*cos_jik)**g3_zeta(p)
                    rad = dexp(-g3_eta(p)*(r2_ji+r2_jk+r2_ik))

                    dang  = g3_zeta(p)*(1.d0+g3_lambda(p)*cos_jik)**(g3_zeta(p)-1.d0)*g3_lambda(p)
                    drad_ji(1:3) = -2.d0*g3_eta(p)*rad*dxp_ji(1:3)
                    drad_jk(1:3) = -2.d0*g3_eta(p)*rad*dxp_jk(1:3)
                    drad_ik(1:3) = -2.d0*g3_eta(p)*rad*dxp_ik(1:3)

                    local_g3_out(p, ja) = local_g3_out(p, ja) + 0.5d0*2.d0**(1.d0-g3_zeta(p)) * ang * rad * fcut_ji * fcut_jk * fcut_ik
                    local_g3_deriv_out(p, ja, ia_n, 1:3) = local_g3_deriv_out(p, ja, ia_n,  1:3) + 0.5d0*2d0**(1d0-g3_zeta(p))*( &
                    (dang*rad*dcos_jik_ji(1:3) + ang*(drad_ji(1:3) - drad_ik(1:3)) )*fcut_ji*fcut_jk*fcut_ik + &
                       ang*rad*fcut_jk*(dfcut_ji(1:3)*fcut_ik - fcut_ji*dfcut_ik(1:3) )  )


                    local_g3_deriv_out(p, ja, ka_n, 1:3) = local_g3_deriv_out(p, ja, ka_n,  1:3) + 0.5d0*2d0**(1d0-g3_zeta(p))*( &
                    (dang*rad*dcos_jik_jk(1:3) + ang*(drad_jk(1:3) + drad_ik(1:3)) )*fcut_ji*fcut_jk*fcut_ik + &
                       ang*rad*fcut_ji*(dfcut_jk(1:3)*fcut_ik + fcut_jk*dfcut_ik(1:3) )  )

                    !local_g3_deriv_out(p, ja, 0, 1:3) = local_g3_deriv_out(p, ja, 0,  1:3) - 0.5d0*2d0**(1d0-g3_zeta(p))*( &
                    !(dang*rad*(dcos_jik_ji(1:3) + dcos_jik_jk) + ang*(drad_ji(1:3) + drad_jk(1:3)) )*fcut_ji*fcut_jk*fcut_ik + &
                    !   ang*rad*fcut_ik*(dfcut_ji(1:3)*fcut_jk + fcut_ji*dfcut_jk(1:3) )  )

           end do
        end do
        local_d_kind_neigh(ja,ia_n)=ia
     end do

  local_d_n_neigh(ja)=ia_n
  end do
  local_g3_deriv_out(:,:,0,:) = - SUM(local_g3_deriv_out(:,:,1:imm_neigh,:),dim=3)

 deallocate(xpnp)
return
end subroutine compute_g3
end module
