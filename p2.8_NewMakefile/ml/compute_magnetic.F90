module magnetic_parameters_rho
USE T_kind_param_m, ONLY:  double

real(double), dimension(:), allocatable  :: a_m, b_m
real(double), dimension(:), allocatable :: c_spline, k_spline

end module magnetic_parameters_rho

subroutine compute_magnetic_sld(i_start_at,i_final_at, local_d_n_neigh, local_d_kind_neigh,  iconf)

USE T_kind_param_m, ONLY:  double
use gen_com_m, ONLY: imm,A2cm,lperiod,bg,at,indi2
use tab_imm_m, ONLY : iwmax2,xp
use ml_in_ndm_module, ONLY: mconf,w2_rho,weighted, &
                            one_pi, imm_neigh, &
                            magnetic_alpha , magnetic_sld_j_dim, magnetic_sld_s2_dim, magnetic_sld_s4_dim, &
                            magnetic_rho_s2, magnetic_rho_s4, &
                            r_cut,  r_cut_magnetic, &
                            magnetic_sld_dim, factor_weight_mass, desc_forces, linvisible
use derived_types, only: config_real, config_desc
implicit none


integer, intent (in) :: i_start_at,i_final_at
integer, dimension(imm),intent(out)  :: local_d_n_neigh
integer, dimension(imm,imm_neigh), intent(out) :: local_d_kind_neigh
!double precision,dimension(magnetic_sld_dim,imm),intent(out) :: local_magnetic_sld_out
!double precision,dimension(magnetic_sld_dim,imm, 0:imm_neigh,3),intent(out) :: local_magnetic_sld_deriv_out
integer, optional, intent(in) :: iconf

real(double), dimension(:,:), allocatable :: xpnp
real(double), dimension(3) :: dxp_ji, ds
integer :: iw,iw1,iw2

integer :: p, p_desc
real(double) :: fcut,r_ji,dfcut, f_out, df_out
real(double) :: factor_ia, factor_ja
real(double) :: magnetic_sld_func, tij_mag, d_tij_mag, rho_ia, rho_ja, A_out, B_out, dA_out, dB_out

integer :: ia,ja,ia_n
logical :: small

if ((i_start_at==0) .and. (i_final_at==0)) then
  local_d_kind_neigh(:,:)=0
  local_d_n_neigh(:)=0
  config_desc(iconf)%energy(:,:)=0.d0
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
config_desc(iconf)%energy(:,:)=0.d0
config_desc(iconf)%force(:,:,:,:)=0.d0
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
      !debug write(*,'("rangml, debug neigh",4i5)') rangml, ja, iw1,  iw2
      !debug write(*,*) rangml, config_real(iconf)%filename
     else
      iw1=iw2+1
      iw2=iwmax2(ja)
      !write(*,*) 'debug neigh', iw1, iw2
     end if
     ia_n=0
     rho_ja=0.d0
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
           !factor_ja=massat(ityp(ja))
           !factor_ia=massat(ityp(ia))
           factor_ia=config_real(iconf)%mass_per_type(config_real(iconf)%itype(ia))/factor_weight_mass
           factor_ja=config_real(iconf)%mass_per_type(config_real(iconf)%itype(ja))/factor_weight_mass


        else
            factor_ia=1.d0
            factor_ja=1.d0
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

        if (r_ji > r_cut_magnetic) cycle
        ia_n = ia_n + 1

        fcut=factor_ia*DOT_PRODUCT(config_real(iconf)%atomic_spin(:,ia), config_real(iconf)%atomic_spin(:,ja))
        !dfcut=factor_ia**config_real(iconf)%atomic_spin(ia)
        !debug write(*,*) fcut, dfcut, factor_ia
        do p=1,magnetic_sld_j_dim
              call exchange_function(r_ji, magnetic_alpha(p), r_cut_magnetic, f_out, df_out)
              magnetic_sld_func = f_out*fcut
              config_desc(iconf)%energy(p,ja) = config_desc(iconf)%energy(p, ja) + magnetic_sld_func
              if (desc_forces) config_desc(iconf)%force(p,ja, ia_n, 1:3) = config_desc(iconf)%force(p,ja, ia_n, 1:3) + df_out*dxp_ji(1:3)/r_ji*fcut
              if (desc_forces) config_desc(iconf)%force(p,ja, 0, 1:3) = config_desc(iconf)%force(p,ja, 0, 1:3) - config_desc(iconf)%force(p,ja, ia_n, 1:3)
        enddo
        rho_ja = rho_ja + tij_mag(r_ji)**2*factor_ia
      local_d_kind_neigh(ja,ia_n)=ia
     end do  !ia
     do p=1 , magnetic_sld_s2_dim
        p_desc = p + magnetic_sld_j_dim
        call A_rho_magnetic(rho_ja,magnetic_rho_s2(p),A_out,dA_out)
        !write(*,*) rho_ja, magnetic_rho_s2(p), A_out, dA_out
        config_desc(iconf)%energy(p_desc,ja) = config_desc(iconf)%energy(p_desc, ja) + A_out*SUM(config_real(iconf)%atomic_spin(:,ja)**2)
        !config_desc(iconf)%force(p,ja, ia_n, 1:3) = config_desc(iconf)%force(p,ja, ia_n, 1:3) + df_out*dxp_ji(1:3)/r_ji*fcut
        !config_desc(iconf)%force(p,ja, 0, 1:3) = config_desc(iconf)%force(p,ja, 0, 1:3) - config_desc(iconf)%force(p,ja, ia_n, 1:3)
     end do
     do p=1 , magnetic_sld_s4_dim
        p_desc = magnetic_sld_s2_dim+magnetic_sld_j_dim+p
        call B_rho_magnetic(rho_ja,magnetic_rho_s4(p),B_out,dB_out)
        config_desc(iconf)%energy(p_desc,ja) = config_desc(iconf)%energy(p_desc, ja) + B_out*SUM(config_real(iconf)%atomic_spin(:,ja)**2)**2
        !config_desc(iconf)%force(p,ja, ia_n, 1:3) = config_desc(iconf)%force(p,ja, ia_n, 1:3) + df_out*dxp_ji(1:3)/r_ji*fcut
        !config_desc(iconf)%force(p,ja, 0, 1:3) = config_desc(iconf)%force(p,ja, 0, 1:3) - config_desc(iconf)%force(p,ja, ia_n, 1:3)
     end do


     config_desc(iconf)%energy(:,:) = config_desc(iconf)%energy(:,:) * factor_ja !/dble(ia_n)
     if (desc_forces) config_desc(iconf)%force(:,:,:,:) = config_desc(iconf)%force(:,:,:,:) * factor_ja!/dble(ia_n)
     local_d_n_neigh(ja)=ia_n
  end do  !ja

 deallocate(xpnp)

return
end subroutine compute_magnetic_sld


subroutine exchange_function(r, alpha, r_cut, f_out, df_out)
USE T_kind_param_m, ONLY:  double
real(double), intent(in)  :: r, alpha, r_cut
real(double), intent(out) ::  f_out, df_out
real(double) :: r_temp


unit=0.176d0
r_temp = 1.d0 - r / r_cut

f_out = unit*r_temp**alpha
df_out = - unit*alpha / r_cut * f_out / r_temp

return
end subroutine exchange_function


subroutine  init_magnetic_sld
use T_kind_param_m, only: double
use ml_in_ndm_module, only: magnetic_sld_dim, magnetic_sld_j_dim, magnetic_sld_s2_dim, magnetic_sld_s4_dim, &
                            magnetic_alpha, magnetic_rho_s2, magnetic_rho_s4

real(double) :: alpha_ini, alpha_fin


if (allocated(magnetic_alpha)) deallocate(magnetic_alpha) ; allocate(magnetic_alpha(magnetic_sld_j_dim))
if (allocated(magnetic_rho_s2)) deallocate(magnetic_rho_s2) ; allocate(magnetic_rho_s2(magnetic_sld_s2_dim))
if (allocated(magnetic_rho_s4)) deallocate(magnetic_rho_s4) ; allocate(magnetic_rho_s4(magnetic_sld_s4_dim))

!..................init J_ij part
alpha_ini=2.0
alpha_fin=7.0


if (magnetic_sld_j_dim==1) then
   magnetic_alpha(1) = (alpha_ini  + alpha_fin)/2.d0
else
  do i=1,magnetic_sld_j_dim
    magnetic_alpha(i) = alpha_ini + (alpha_fin - alpha_ini)* dble(i-1) / dble(magnetic_sld_j_dim-1)
  end do
end if

!...............init S2 part
alpha_ini=0.d0
alpha_fin=1.d0


if (magnetic_sld_s2_dim==1) then
   magnetic_rho_s2(1) = (alpha_ini  + alpha_fin)/2.d0
else
  do i=1,magnetic_sld_s2_dim
    magnetic_rho_s2(i) = alpha_ini + (alpha_fin - alpha_ini)* dble(i-1) / dble(magnetic_sld_s2_dim-1)
  end do
end if

!.................init S4 part
alpha_ini=0.d0
alpha_fin=1.d0


if (magnetic_sld_s4_dim==1) then
   magnetic_rho_s4(1) = (alpha_ini  + alpha_fin)/2.d0
else
  do i=1,magnetic_sld_s4_dim
    magnetic_rho_s4(i) = alpha_ini + (alpha_fin - alpha_ini)* dble(i-1) / dble(magnetic_sld_s4_dim-1)
  end do
end if



!...............final dimension of magnetic desriptor.
magnetic_sld_dim = magnetic_sld_j_dim + magnetic_sld_s2_dim  + magnetic_sld_s4_dim

return
end subroutine init_magnetic_sld


subroutine init_magnetic_parameters_rho

use magnetic_parameters_rho, only: a_m, b_m, c_spline, k_spline

if (allocated(a_m)) deallocate(a_m) ; allocate(a_m(3))
if (allocated(b_m)) deallocate(b_m) ; allocate(b_m(3))
if (allocated(c_spline)) deallocate(c_spline) ; allocate(c_spline(7))
if (allocated(k_spline)) deallocate(k_spline) ; allocate(k_spline(7))


! Landau A fitted parameters
a_m(1)=-2.3827723674043900D-01
a_m(2)=1.2945703172205700D-02
a_m(3)=-1.1518969922985000D-04

! Landau B fitted parameters
b_m(1)=1.06003150785869D-02
b_m(2)=1.6104913287021000D-03
b_m(3)=-4.3178188078544200D-05

! Hopping parameters:
!   r^t_n = c_spline
!   t_n   = k_spline
c_spline(1)=2.0000000000000000D+00
c_spline(2)=2.2000000000000000D+00
c_spline(3)=2.6000000000000000D+00
c_spline(4)=3.2000000000000000D+00
c_spline(5)=3.8000000000000000D+00
c_spline(6)=4.6000000000000000D+00
c_spline(7)=5.3000000000000000D+00

k_spline(1)=2.5999782982854347D+00
k_spline(2)=2.9319480072508499D+00
k_spline(3)=-2.8388905185188360D+00
k_spline(4)=-1.0267419494754382D-01
k_spline(5)=1.5484736035888333D-02
k_spline(6)=-7.2805743511785065D-02
k_spline(7)=-3.6343523861565924D-03



return
end subroutine init_magnetic_parameters_rho


function  tij_mag(x) result (value)
  use magnetic_parameters_rho, only: c_spline, k_spline
  implicit none
  real(kind(0.d0)), intent(in) :: x
  real(kind(0.d0)) :: value,step
  integer :: i
  !
  value=0.d0
  do i =1,size(c_spline,1)
    value = value + c_spline(i) *step(k_spline(i)-x)*(k_spline(i) - x)**3
  end do
  !
end function tij_mag



function  d_tij_mag(x) result (value)
  use magnetic_parameters_rho, only: c_spline, k_spline
  implicit none
  real(kind(0.d0)), intent(in) :: x
   real(kind(0.d0)) ::  value, step
  integer :: i
  !
  value=0.d0
  do i =1,size(c_spline,1)
    value = value - 2.d0*c_spline(i) *step(k_spline(i)-x)*(k_spline(i) - x)**2
  end do
  !
end function d_tij_mag


subroutine A_rho_magnetic(x, alpha, A_out, dA_out)
use T_kind_param_m, only : double
use magnetic_parameters_rho, only: a_m
real(double), intent(in) :: x, alpha
real(double), intent(out) :: A_out, dA_out

A_out = a_m(1) + 2.0*a_m(2) * x *(1.d0-alpha)+ alpha*2.0*a_m(3)*x**2
dA_out = A_out

end subroutine A_rho_magnetic


subroutine B_rho_magnetic(x, alpha, B_out, dB_out)
use T_kind_param_m, only : double
use magnetic_parameters_rho, only: b_m
real(double), intent(in) :: x, alpha
real(double), intent(out) :: B_out, dB_out

B_out = b_m(1) + 2.0*b_m(2) * x *(1.d0-alpha)+ alpha*2.0*b_m(3)*x**2
dB_out=B_out

end subroutine B_rho_magnetic





function  step(x) result (value)
  implicit none
  real(kind(0.d0)) :: x, value

  if (x>0.d0) then
    value=1.d0
  else
    value =0.d0
  end if

end function step
