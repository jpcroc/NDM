!subroutine compute_mtp(i_start_at,i_final_at,d_n_neigh, d_kind_neigh, local_mtp_out,local_mtp_deriv_out, iconf)
!LM
subroutine compute_mtp(i_start_at,i_final_at,d_n_neigh, d_kind_neigh, iconf)

use T_kind_param_m, ONLY:  double
use gen_com_m, ONLY: imm, A2cm,lperiod,bg,at,indi2
use tab_imm_m, ONLY : iwmax2,xp
use angular_functions
use ml_in_ndm_module, ONLY: rangml,pi,r_cut, imm_neigh, mconf,   &
                            w2_rho,weighted,mtp_dim, mtp_rad_order, &
                            mtp_poly_min, factor_weight_mass, desc_forces, linvisible
use derived_types, only : config_real, config_desc

implicit none
integer, intent (in) :: i_start_at,i_final_at
integer, dimension(imm),intent(out)  :: d_n_neigh
integer, dimension(imm,imm_neigh), intent(out) :: d_kind_neigh
!real(kind(0.d0)),dimension(mtp_dim,imm),intent(out) :: local_mtp_out
!real(kind(0.d0)),dimension(mtp_dim,imm, 0:imm_neigh, 3),intent(out) :: local_mtp_deriv_out
integer, optional :: iconf

logical :: small

real(kind(0.d0)),dimension(mtp_rad_order)                ::   tmp_mtp1
real(kind(0.d0)),dimension(mtp_rad_order, 3)             ::   tmp_mtp2
real(kind(0.d0)),dimension(mtp_rad_order, 3, 3)          ::   tmp_mtp3
real(kind(0.d0)),dimension(mtp_rad_order, 3, 3, 3)       ::   tmp_mtp4
real(kind(0.d0)),dimension(mtp_rad_order, 3, 3, 3, 3)    ::   tmp_mtp5
real(kind(0.d0)),dimension(mtp_rad_order, 3, 3, 3, 3, 3) ::   tmp_mtp6
real(kind(0.d0)),dimension(:, :, :),                  allocatable :: d_tmp_mtp1
real(kind(0.d0)),dimension(:, :, :, :),               allocatable :: d_tmp_mtp2
real(kind(0.d0)),dimension(:, :, :, :, :),            allocatable :: d_tmp_mtp3
real(kind(0.d0)),dimension(:, :, :, :, :, :),         allocatable :: d_tmp_mtp4
real(kind(0.d0)),dimension(:, :, :, :, :, :, :),      allocatable :: d_tmp_mtp5
real(kind(0.d0)),dimension(:, :, :, :, :, :, :, : ),  allocatable :: d_tmp_mtp6



real(kind(0.d0)),dimension(:, :),       allocatable :: vx_ji
real(kind(0.d0)),dimension(:),          allocatable :: vr_ji


real(kind(0.d0)),dimension(:, :),       allocatable :: tmpx_ji
real(kind(0.d0)),dimension(:),          allocatable :: tmpr_ji


real(kind(0.d0)),dimension(:,:,:),          allocatable :: f_vr_ji
real(kind(0.d0)),dimension(:,:,:),        allocatable :: d_f_vr_ji



real(double), dimension(:,:), allocatable :: xpnp
real(double), dimension(3) :: dxp_ji, ds
real(kind(0.d0)) :: tmp1, tmp2, tmp3, tmp4, tmp5
integer :: iw,iw1,iw2

integer :: ip, ipo, ip1, ip2, ip3, ix, iy, iz
integer :: b1,b2,b3,b4,b5,nu
integer :: ia, ja, ia_n, icount, icount_e, icount_f
double precision :: r2_ji,r_ji, a_ip1_ip2, b_ip1_ip2, c_ip1_ip2
double precision :: factor_ia, factor_ja, fcut, dfcut,  tmp_local, tmp_local_3d(3)


if (desc_forces) then
  if(allocated(d_tmp_mtp1)) deallocate(d_tmp_mtp1) ;    allocate(d_tmp_mtp1(mtp_rad_order, imm_neigh, 3))
  if(allocated(d_tmp_mtp2)) deallocate(d_tmp_mtp2) ;    allocate(d_tmp_mtp2(mtp_rad_order, 3, imm_neigh, 3))
  if(allocated(d_tmp_mtp3)) deallocate(d_tmp_mtp3) ;    allocate(d_tmp_mtp3(mtp_rad_order, 3, 3, imm_neigh, 3))
  if(allocated(d_tmp_mtp4)) deallocate(d_tmp_mtp4) ;    allocate(d_tmp_mtp4(mtp_rad_order, 3, 3, 3, imm_neigh, 3))
  if(allocated(d_tmp_mtp5)) deallocate(d_tmp_mtp5) ;    allocate(d_tmp_mtp5(mtp_rad_order, 3, 3, 3, 3, imm_neigh, 3))
  if(allocated(d_tmp_mtp6)) deallocate(d_tmp_mtp6) ;    allocate(d_tmp_mtp6(mtp_rad_order, 3, 3, 3, 3, 3, imm_neigh, 3))
end if



if ((i_start_at==0) .and. (i_final_at==0)) then
  d_n_neigh(:)=0
  d_kind_neigh(:,:)=0
  config_desc(iconf)%energy(:,:)=0.d0
!  local_mtp_deriv_out(:,:,:,:)=0.d0
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
config_desc(iconf)%energy(:,:)=0.d0
if (desc_forces) config_desc(iconf)%force(:,:,:,:)=0.d0

!if (allocated(m_typ)) deallocate(m_typ); allocate(m_typ(ntyp,ntyp))
!call mconf(m_typ(:,:))

if (i_start_at==1)  iw2=0
if (i_start_at > 1) iw2=iwmax2(i_start_at-1)

do ja=i_start_at,i_final_at
    if (allocated(tmpx_ji)) deallocate(tmpx_ji) ; allocate (tmpx_ji(3,imm_neigh))
    if (allocated(tmpr_ji)) deallocate(tmpr_ji) ; allocate (tmpr_ji(imm_neigh))
    !begin small box or not 1/


    if (linvisible.and.weighted) then
        if( config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ja)) ) cycle
    end if


    if (weighted) then
        factor_ja=config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))!/factor_weight_mass
    else
        factor_ja=1.d0
    endif



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
    tmp_mtp1(:)=0.d0
    if (desc_forces) d_tmp_mtp1(:,:,:)=0.d0

    tmp_mtp2(:,:)=0.d0
    if (desc_forces) d_tmp_mtp2(:,:,:,:)=0.d0

    tmp_mtp3(:,:,:)=0.d0
    if (desc_forces) d_tmp_mtp3(:,:,:,:,:)=0.d0

    tmp_mtp4(:,:,:,:)=0.d0
    if (desc_forces) d_tmp_mtp4(:,:,:,:,:,:)=0.d0

    tmp_mtp5(:,:,:,:,:)=0.d0
    if (desc_forces) d_tmp_mtp5(:,:,:,:,:,:,:)=0.d0


    tmp_mtp6(:,:,:,:,:,:)=0.d0
    if (desc_forces) d_tmp_mtp6(:,:,:,:,:,:,:,:)=0.d0


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
        tmpx_ji(:,ia_n) = dxp_ji(:)
        tmpr_ji(ia_n) = r_ji
        fcut=0.5d0*(cos(pi*r_ji/r_cut)+1d0)*factor_ia
        dfcut= -0.5d0*pi/r_cut*sin(pi*r_ji/r_cut)*factor_ia
        d_kind_neigh(ja,ia_n)=ia
    end do  !iw
    !ia_n is the max of neighbours ...
    d_n_neigh(ja)=ia_n
    if (ia_n > imm_neigh) then
      write(6,*) 'Fatal for atom ja for which the number of neighbours is higher than the admitted MAXIMUM', ja, imm_neigh
    end if


    if (allocated(vx_ji)) deallocate(vx_ji) ; allocate (vx_ji(3,ia_n))
    if (allocated(vr_ji)) deallocate(vr_ji) ; allocate (vr_ji(ia_n))

    vx_ji(1:3,1:ia_n) = tmpx_ji(1:3,1:ia_n)
    deallocate(tmpx_ji)
    vr_ji(1:ia_n) = tmpr_ji(1:ia_n)
    deallocate(tmpr_ji)

    if (allocated(f_vr_ji)) deallocate(f_vr_ji) ; allocate (f_vr_ji(mtp_rad_order,0:5,ia_n))
    if (allocated(d_f_vr_ji)) deallocate(d_f_vr_ji) ; allocate (d_f_vr_ji(mtp_rad_order,0:5,ia_n))


    ! initialize the function ...
    do ip = 1,mtp_rad_order
      ipo = -ip-mtp_poly_min+1
      do nu=0,5
        f_vr_ji(ip,nu, :) = vr_ji(:)**(ipo-nu)
        d_f_vr_ji(ip,nu, :) = dble(ipo)*vr_ji(:)**(ipo-2-nu)
     end do
    end do


    do ip = 1,mtp_rad_order
      ipo = -ip-mtp_poly_min+1

      ! zero order ...
      !------------------------------------------
      ! mtp(ipo, 0) ~ mtp(ipo)
      !------------------------------------------


      tmp_mtp1(ip) = SUM(f_vr_ji(ip,0,:))

      if (desc_forces) then
        do ia=1,ia_n
          d_tmp_mtp1(ip,ia,:) =  d_f_vr_ji(ip,0,ia)*vx_ji(:,ia)
        end do
      end if

      !second order ...
      !------------------------------------------!
      ! mtp(ip, 1) ~ mtp(ip, b1)                 !
      !------------------------------------------!
      do b1=1,3
          do ia=1,ia_n
             tmp_mtp2(ip,b1) = tmp_mtp2(ip,b1) + f_vr_ji(ip,1,ia)*vx_ji(b1,ia)
          end do
      end do

      if (desc_forces) then
        do b1=1,3
          do ia=1,ia_n
            do ix=1,3
              tmp1=0.d0
              if (ix==b1) tmp1 = f_vr_ji(ip,1,ia)
              d_tmp_mtp2(ip,b1,ia,ix)=d_f_vr_ji(ip,1,ia)*vx_ji(ix,ia)*vx_ji(b1,ia) + tmp1
            end do
          end do
        end do
      end if


      !third order ...
      !------------------------------------------!
      ! mtp(ip, 2) ~ mtp(ip, b1, b2)             !
      !------------------------------------------!
      do b1=1,3
        do b2=1,3
         do ia=1,ia_n
         tmp_mtp3(ip,b1,b2) =tmp_mtp3(ip,b1,b2) +  f_vr_ji(ip,2,ia)*vx_ji(b1,ia)*vx_ji(b2,ia)
         end do
        end do
      end do

      if (desc_forces) then
        do b1=1,3
          do b2=1,3
            do ia=1,ia_n
              do ix=1,3
                tmp1=0.d0
                tmp2=0.d0
                if (b1==ix) tmp1 =   f_vr_ji(ip,2,ia)*vx_ji(b2,ia)
                if (b2==ix) tmp2 =   f_vr_ji(ip,2,ia)*vx_ji(b1,ia)
                d_tmp_mtp3(ip,b1,b2,ia,ix)=d_f_vr_ji(ip,2,ia)*vx_ji(ix,ia)*vx_ji(b1,ia)*vx_ji(b2,ia) + tmp1 + tmp2
              end do
            end do
          end do
        end do
      end if



      !fourth order ...
      !------------------------------------------!
      ! mtp(ip, 3) ~ mtp(ip, b1, b2,b3)          !
      !------------------------------------------!
      do b1=1,3
        do b2=1,3
          do b3=1,3
            do ia=1,ia_n
              tmp_mtp4(ip,b1,b2,b3) =tmp_mtp4(ip,b1,b2,b3) +  f_vr_ji(ip,3,ia)*vx_ji(b1,ia)*vx_ji(b2,ia)*vx_ji(b3,ia)
            end do
          end do
        end do
      end do

      if (desc_forces) then
        do b1=1,3
          do b2=1,3
            do b3=1,3
              do ia=1,ia_n
                do ix=1,3
                  tmp1=0.d0
                  tmp2=0.d0
                  tmp3=0.d0
                  if (b1==ix) tmp1 =   f_vr_ji(ip,3,ia)*vx_ji(b2,ia)*vx_ji(b3,ia)
                  if (b2==ix) tmp2 =   f_vr_ji(ip,3,ia)*vx_ji(b1,ia)*vx_ji(b3,ia)
                  if (b3==ix) tmp3 =   f_vr_ji(ip,3,ia)*vx_ji(b1,ia)*vx_ji(b2,ia)
                  d_tmp_mtp4(ip,b1,b2,b3,ia,ix)=d_f_vr_ji(ip,3,ia)*vx_ji(ix,ia)*vx_ji(b1,ia)*vx_ji(b2,ia)*vx_ji(b3,ia) + tmp1 + tmp2 + tmp3
                end do
              end do
            end do
          end do
        end do
      end if



      !fifth order ...
      !------------------------------------------!
      ! mtp(ip, 4) ~ mtp(ip, b1, b2, b3, b4)     !
      !------------------------------------------!
      do b1=1,3
        do b2=1,3
          do b3=1,3
            do b4=1,3
              do ia=1,ia_n
               tmp_mtp5(ip,b1,b2,b3,b4) =tmp_mtp5(ip,b1,b2,b3,b4) +  f_vr_ji(ip,4,ia)*vx_ji(b1,ia)*vx_ji(b2,ia)*vx_ji(b3,ia)*vx_ji(b4,ia)
              end do
            end do
          end do
        end do
      end do

      if (desc_forces) then
        do b1=1,3
          do b2=1,3
            do b3=1,3
              do b4=1,3
                do ia=1,ia_n
                  do ix=1,3
                    tmp1=0.d0
                    tmp2=0.d0
                    tmp3=0.d0
                    tmp4=0.d0
                    if (b1==ix) tmp1 =   f_vr_ji(ip,4,ia)*vx_ji(b2,ia)*vx_ji(b3,ia)*vx_ji(b4,ia)
                    if (b2==ix) tmp2 =   f_vr_ji(ip,4,ia)*vx_ji(b1,ia)*vx_ji(b3,ia)*vx_ji(b4,ia)
                    if (b3==ix) tmp3 =   f_vr_ji(ip,4,ia)*vx_ji(b1,ia)*vx_ji(b2,ia)*vx_ji(b4,ia)
                    if (b4==ix) tmp4 =   f_vr_ji(ip,4,ia)*vx_ji(b1,ia)*vx_ji(b2,ia)*vx_ji(b3,ia)
                    d_tmp_mtp5(ip,b1,b2,b3,b4,ia,ix)=d_f_vr_ji(ip,4,ia)*vx_ji(ix,ia)*vx_ji(b1,ia)*vx_ji(b2,ia)*vx_ji(b3,ia)*vx_ji(b4,ia) + tmp1 + tmp2 + tmp3 + tmp4
                  end do
                end do
              end do
            end do
          end do
        end do
      end if


      !sixth order ...
      !------------------------------------------!
      ! mtp(ip, 5) ~ mtp(ip, b1, b2, b3, b4, b5) !
      !------------------------------------------!
      do b1=1,3
        do b2=1,3
          do b3=1,3
            do b4=1,3
              do b5=1,3
                do ia=1,ia_n
                tmp_mtp6(ip,b1,b2,b3,b4,b5) =tmp_mtp6(ip,b1,b2,b3,b4,b5) +  f_vr_ji(ip,5,ia)*vx_ji(b1,ia)*vx_ji(b2,ia)*vx_ji(b3,ia)*vx_ji(b4,ia)*vx_ji(b5,ia)
                end do
              end do
            end do
          end do
        end do
      end do

      if (desc_forces) then
        do b1=1,3
          do b2=1,3
            do b3=1,3
              do b4=1,3
                do b5=1,3
                  do ia=1,ia_n
                    do ix=1,3
                      tmp1=0.d0
                      tmp2=0.d0
                      tmp3=0.d0
                      tmp4=0.d0
                      tmp5=0.d0
                      if (b1==ix) tmp1 =   f_vr_ji(ip,5,ia)*vx_ji(b2,ia)*vx_ji(b3,ia)*vx_ji(b4,ia)*vx_ji(b5,ia)
                      if (b2==ix) tmp2 =   f_vr_ji(ip,5,ia)*vx_ji(b1,ia)*vx_ji(b3,ia)*vx_ji(b4,ia)*vx_ji(b5,ia)
                      if (b3==ix) tmp3 =   f_vr_ji(ip,5,ia)*vx_ji(b1,ia)*vx_ji(b2,ia)*vx_ji(b4,ia)*vx_ji(b5,ia)
                      if (b4==ix) tmp4 =   f_vr_ji(ip,5,ia)*vx_ji(b1,ia)*vx_ji(b2,ia)*vx_ji(b3,ia)*vx_ji(b5,ia)
                      if (b5==ix) tmp5 =   f_vr_ji(ip,5,ia)*vx_ji(b1,ia)*vx_ji(b2,ia)*vx_ji(b3,ia)*vx_ji(b4,ia)
                      d_tmp_mtp6(ip,b1,b2,b3,b4,b5,ia,ix)=d_f_vr_ji(ip,5,ia)*vx_ji(ix,ia)*vx_ji(b1,ia)*vx_ji(b2,ia)*vx_ji(b3,ia)*vx_ji(b4,ia)*vx_ji(b5,ia) &
                                                          + tmp1 + tmp2 + tmp3 + tmp4 + tmp5
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end if

    end do !...first ip.
    ! compute the descriptor  ... the reduction of the above tensors
    !01---------------------------------!
    ! alpha = ( \mu  )           !
    !-----------------------------------!
    config_desc(iconf)%energy(1:mtp_rad_order,ja) = tmp_mtp1(1:mtp_rad_order)
    if (desc_forces) config_desc(iconf)%force(1:mtp_rad_order,ja,1:ia_n,1:3) = d_tmp_mtp1(1:mtp_rad_order,1:ia_n,1:3)

    !02---------------------------------!
    ! alpha = ( \mu_1      1)           !
    !         (     1  \mu_2)           !
    !             -----                 !
    ! yields \mu1*\mu2 components       !
    !-----------------------------------!
    !----------------------------!
    config_desc(iconf)%energy(mtp_rad_order+1:mtp_rad_order + mtp_rad_order**2,ja) = RESHAPE(MATMUL(tmp_mtp2,TRANSPOSE(tmp_mtp2)), (/mtp_rad_order**2/))
    icount_e=mtp_rad_order + mtp_rad_order**2+1

    if (desc_forces) then
        icount_f = mtp_rad_order+1
        do ip1=1,mtp_rad_order
          do ip2=1,mtp_rad_order
             do b1=1,3
               do ia=1,ia_n
               do ix=1,3
               a_ip1_ip2 =  d_tmp_mtp2(ip1, b1, ia, ix) * tmp_mtp2(ip2, b1)
               b_ip1_ip2 =  tmp_mtp2(ip1, b1) * d_tmp_mtp2 (ip2, b1, ia, ix)
               config_desc(iconf)%force(icount_f,ja,ia,ix) = config_desc(iconf)%force(icount_f,ja,ia,ix)  + a_ip1_ip2  + b_ip1_ip2
               !config_desc(iconf)%force(icount,ja,ia,iy) =  a_ip1_ip2  + b_ip1_ip2
               end do
               end do
             end do
             icount_f = icount_f + 1
          end do
         end do
    end if

    !02a---------------------------------!
    ! alpha = ( \mu_1      2)           !
    !         (     2  \mu_2)           !
    !             -----                 !
    ! yields \mu1*\mu2 components       !
    !-----------------------------------!


    do ip1=1,mtp_rad_order
      do ip2=1,mtp_rad_order
        do b1=1,3
          do b2=1,3
            !config_desc(iconf)%energy(icount,ja)= config_desc(iconf)%energy(icount,ja) + tmp_mtp3(ip1,ix,iy)*tmp_mtp3(ip2,ix,iy)
            config_desc(iconf)%energy(icount_e,ja)= config_desc(iconf)%energy(icount_e,ja) + tmp_mtp3(ip1,b1,b2)*tmp_mtp3(ip2,b1,b2)
          end do
        end do
        icount_e=icount_e+1
      end do
    end do
    !
    if (desc_forces) then
        do ip1=1,mtp_rad_order
          do ip2=1,mtp_rad_order
             do b1=1,3
               do b2=1,3
                do ia=1,ia_n
                 do ix=1,3

                   a_ip1_ip2 =  d_tmp_mtp3(ip1, b1, b2, ia, ix) * tmp_mtp3(ip2, b1, b2)
                   b_ip1_ip2 =    tmp_mtp3(ip1, b1, b2) *       d_tmp_mtp3(ip2, b1, b2, ia, ix)
                   config_desc(iconf)%force(icount_f,ja,ia,ix) = config_desc(iconf)%force(icount_f,ja,ia,ix)  + a_ip1_ip2 + b_ip1_ip2
                   !config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2
                 end do
                end do
               end do
             end do
             icount_f = icount_f + 1
          end do
         end do
    end if


    !02b---------------------------------!
    ! alpha = ( \mu_1      3)           !
    !         (     3  \mu_2)           !
    !             -----                 !
    ! yields \mu1*\mu2 components       !
    !-----------------------------------!


    do ip1=1,mtp_rad_order
      do ip2=1,mtp_rad_order
        do b1=1,3
          do b2=1,3
            do b3=1,3
              !config_desc(iconf)%energy(icount,ja)= config_desc(iconf)%energy(icount,ja) + tmp_mtp3(ip1,ix,iy)*tmp_mtp3(ip2,ix,iy)
              config_desc(iconf)%energy(icount_e,ja)= config_desc(iconf)%energy(icount_e,ja) + tmp_mtp4(ip1,b1,b2,b3)*tmp_mtp4(ip2,b1,b2,b3)
            end do
          end do
        end do
        icount_e=icount_e+1
      end do
    end do
    !
    if (desc_forces) then
        do ip1=1,mtp_rad_order
          do ip2=1,mtp_rad_order
            do b1=1,3
              do b2=1,3
                do b3=1,3
                  do ia=1,ia_n
                    do ix=1,3

                      a_ip1_ip2 =  d_tmp_mtp4(ip1, b1, b2, b3, ia, ix) * tmp_mtp4(ip2, b1, b2,b3)
                      b_ip1_ip2 =    tmp_mtp4(ip1, b1, b2, b3) *       d_tmp_mtp4(ip2, b1, b2, b3, ia, ix)
                      config_desc(iconf)%force(icount_f,ja,ia,ix) = config_desc(iconf)%force(icount_f,ja,ia,ix)  + a_ip1_ip2 + b_ip1_ip2
                      !config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2
                    end do
                  end do
                end do
              end do
            end do
            icount_f = icount_f + 1
          end do
         end do
    end if


    !02c---------------------------------!
    ! alpha = ( \mu_1      4)           !
    !         (     4  \mu_2)           !
    !             -----                 !
    ! yields \mu1*\mu2 components       !
    !-----------------------------------!


    do ip1=1,mtp_rad_order
      do ip2=1,mtp_rad_order
        do b1=1,3
          do b2=1,3
            do b3=1,3
              do b4=1,3
                !config_desc(iconf)%energy(icount,ja)= config_desc(iconf)%energy(icount,ja) + tmp_mtp3(ip1,ix,iy)*tmp_mtp3(ip2,ix,iy)
                config_desc(iconf)%energy(icount_e,ja)= config_desc(iconf)%energy(icount_e,ja) + tmp_mtp5(ip1,b1,b2,b3,b4)*tmp_mtp5(ip2,b1,b2,b3,b4)
              end do
            end do
          end do
        end do
        icount_e=icount_e+1
      end do
    end do
    !
    if (desc_forces) then
        do ip1=1,mtp_rad_order
          do ip2=1,mtp_rad_order
            do b1=1,3
              do b2=1,3
                do b3=1,3
                  do b4=1,3
                    do ia=1,ia_n
                      do ix=1,3

                        a_ip1_ip2 =  d_tmp_mtp5(ip1, b1, b2, b3, b4, ia, ix) * tmp_mtp5(ip2, b1, b2, b3, b4)
                        b_ip1_ip2 =    tmp_mtp5(ip1, b1, b2, b3, b4) *       d_tmp_mtp5(ip2, b1, b2, b3, b4, ia, ix)
                        config_desc(iconf)%force(icount_f,ja,ia,ix) = config_desc(iconf)%force(icount_f,ja,ia,ix)  + a_ip1_ip2 + b_ip1_ip2
                        !config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2
                      end do
                    end do
                  end do
                end do
              end do
            end do
            icount_f = icount_f + 1
          end do
         end do
    end if


    !02d---------------------------------!
    ! alpha = ( \mu_1      5)           !
    !         (     5  \mu_2)           !
    !             -----                 !
    ! yields \mu1*\mu2 components       !
    !-----------------------------------!


    do ip1=1,mtp_rad_order
      do ip2=1,mtp_rad_order
        do b1=1,3
          do b2=1,3
            do b3=1,3
              do b4=1,3
                do b5=1,3
                  config_desc(iconf)%energy(icount_e,ja)= config_desc(iconf)%energy(icount_e,ja) + tmp_mtp6(ip1,b1,b2,b3,b4,b5)*tmp_mtp6(ip2,b1,b2,b3,b4,b5)
                end do
              end do
            end do
          end do
        end do
        icount_e=icount_e+1
      end do
    end do
    !
    if (desc_forces) then
        do ip1=1,mtp_rad_order
          do ip2=1,mtp_rad_order
            do b1=1,3
              do b2=1,3
                do b3=1,3
                  do b4=1,3
                    do b5=1,3
                      do ia=1,ia_n
                        do ix=1,3
                          a_ip1_ip2 =  d_tmp_mtp6(ip1, b1, b2, b3, b4, b5, ia, ix) * tmp_mtp6(ip2, b1, b2, b3, b4, b5)
                          b_ip1_ip2 =    tmp_mtp6(ip1, b1, b2, b3, b4, b5) *       d_tmp_mtp6(ip2, b1, b2, b3, b4, b5, ia, ix)
                          config_desc(iconf)%force(icount_f,ja,ia,ix) = config_desc(iconf)%force(icount_f,ja,ia,ix)  + a_ip1_ip2 + b_ip1_ip2
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
            icount_f = icount_f + 1
          end do
        end do
    end if

    !03---------------------------------!
    ! alpha = ( \mu_1      1      1)    !
    !         (     1  \mu_2      0)    !
    !         (     1      0  \mu_3)    !
    !             -----                 !
    ! yields \mu1*\mu2*\mu3 components  !
    !-----------------------------------!


    do ip1=1,mtp_rad_order
      do ip2=1,mtp_rad_order
        do ip3=1,mtp_rad_order
          do b1=1,3
            do b2=1,3
              config_desc(iconf)%energy(icount_e,ja)=  config_desc(iconf)%energy(icount_e,ja) +  tmp_mtp2(ip1,b1)*tmp_mtp3(ip2,b1,b2)*tmp_mtp2(ip3,b2)
            end do
          end do
          icount_e=icount_e+1
        end do
       end do
    end do


    if (desc_forces) then
        do ip1=1,mtp_rad_order
          do ip2=1,mtp_rad_order
            do ip3=1,mtp_rad_order
             do b1=1,3
               do b2=1,3
                do ia=1,ia_n
                 do ix=1,3
                   a_ip1_ip2=  tmp_mtp2(ip1,b1)*    d_tmp_mtp3(ip2,b1,b2,ia,ix)*tmp_mtp2(ip3,b2)
                   b_ip1_ip2=d_tmp_mtp2(ip1,b1,ia,ix)*tmp_mtp3(ip2,b1,b2)*      tmp_mtp2(ip3,b2)
                   c_ip1_ip2=  tmp_mtp2(ip1,b1)*      tmp_mtp3(ip2,b1,b2)*    d_tmp_mtp2(ip3,b2,ia,ix)
                   config_desc(iconf)%force(icount_f,ja,ia,ix) = config_desc(iconf)%force(icount_f,ja,ia,ix)  +  a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                  !config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2

                 end do
                end do
               end do
             end do
             icount_f = icount_f + 1
          end do
         end do
        end do
    end if


    !03a--------------------------------!
    ! alpha = ( \mu_1      1      2)    !
    !         (     1  \mu_2      0)    !
    !         (     2      0  \mu_3)    !
    !             -----                 !
    ! yields \mu1*\mu2*\mu3 components  !
    !-----------------------------------!


    do ip1=1,mtp_rad_order
      do ip2=1,mtp_rad_order
        do ip3=1,mtp_rad_order
          do b1=1,3
            do b2=1,3
              do b3=1,3
                config_desc(iconf)%energy(icount_e,ja)=  config_desc(iconf)%energy(icount_e,ja) +  tmp_mtp2(ip1,b1)*tmp_mtp4(ip2,b1,b2,b3)*tmp_mtp3(ip3,b2,b3)
              end do
            end do
          end do
          icount_e=icount_e+1
        end do
       end do
    end do


    if (desc_forces) then
        do ip1=1,mtp_rad_order
          do ip2=1,mtp_rad_order
            do ip3=1,mtp_rad_order
             do b1=1,3
               do b2=1,3
                do b3=1,3
                  do ia=1,ia_n
                    do ix=1,3
                      a_ip1_ip2=  tmp_mtp2(ip1,b1)      *d_tmp_mtp4(ip2,b1,b2,b3,ia,ix)  *  tmp_mtp3(ip3,b2,b3)
                      b_ip1_ip2=d_tmp_mtp2(ip1,b1,ia,ix)*  tmp_mtp4(ip2,b1,b2,b3)        *  tmp_mtp3(ip3,b2,b3)
                      c_ip1_ip2=  tmp_mtp2(ip1,b1)      *  tmp_mtp4(ip2,b1,b2,b3)        *d_tmp_mtp3(ip3,b2,b3,ia,ix)
                      config_desc(iconf)%force(icount_f,ja,ia,ix) = config_desc(iconf)%force(icount_f,ja,ia,ix)  +  a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                      !config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                    end do
                  end do
                end do
               end do
             end do
             icount_f = icount_f + 1
          end do
         end do
        end do
    end if


    !03b--------------------------------!
    ! alpha = ( \mu_1      1      3)    !
    !         (     1  \mu_2      0)    !
    !         (     3      0  \mu_3)    !
    !             -----                 !
    ! yields \mu1*\mu2*\mu3 components  !
    !-----------------------------------!


    do ip1=1,mtp_rad_order
      do ip2=1,mtp_rad_order
        do ip3=1,mtp_rad_order
          do b1=1,3
            do b2=1,3
              do b3=1,3
                do b4=1,3
                    config_desc(iconf)%energy(icount_e,ja)=  config_desc(iconf)%energy(icount_e,ja) +  tmp_mtp5(ip1,b1,b2,b3,b4)*tmp_mtp2(ip2,b1)*tmp_mtp4(ip3,b2,b3,b4)
                end do
              end do
            end do
          end do
          icount_e=icount_e+1
        end do
       end do
    end do


    if (desc_forces) then
        do ip1=1,mtp_rad_order
          do ip2=1,mtp_rad_order
            do ip3=1,mtp_rad_order
             do b1=1,3
               do b2=1,3
                do b3=1,3
                  do b4=1,3
                      do ia=1,ia_n
                        do ix=1,3
                          a_ip1_ip2=  tmp_mtp5(ip1,b1,b2,b3,b4)       *d_tmp_mtp2(ip2,b1,ia,ix)*tmp_mtp4(ip3,b2,b3,b4)
                          b_ip1_ip2=d_tmp_mtp5(ip1,b1,b2,b3,b4,ia,ix) *  tmp_mtp2(ip2,b1)      *tmp_mtp4(ip3,b2,b3,b4)
                          c_ip1_ip2=  tmp_mtp5(ip1,b1,b2,b3,b4)       *  tmp_mtp2(ip2,b1)    *d_tmp_mtp4(ip3,b2,b3,b4,ia,ix)
                          config_desc(iconf)%force(icount_f,ja,ia,ix) = config_desc(iconf)%force(icount_f,ja,ia,ix)  +  a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                          !config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                        end do
                      end do
                  end do
                end do
               end do
             end do
             icount_f = icount_f + 1
          end do
         end do
        end do
    end if


    !03c--------------------------------!
    ! alpha = ( \mu_1      1      4)    !
    !         (     1  \mu_2      0)    !
    !         (     4      0  \mu_3)    !
    !             -----                 !
    ! yields \mu1*\mu2*\mu3 components  !
    !-----------------------------------!


    do ip1=1,mtp_rad_order
      do ip2=1,mtp_rad_order
        do ip3=1,mtp_rad_order
          do b1=1,3
            do b2=1,3
              do b3=1,3
                do b4=1,3
                  do b5=1,3
                    config_desc(iconf)%energy(icount_e,ja)=  config_desc(iconf)%energy(icount_e,ja) +  tmp_mtp6(ip1,b1,b2,b3,b4,b5)*tmp_mtp2(ip2,b1)*tmp_mtp5(ip3,b2,b3,b4,b5)
                  end do
                end do
              end do
            end do
          end do
          icount_e=icount_e+1
        end do
       end do
    end do


    if (desc_forces) then
        do ip1=1,mtp_rad_order
          do ip2=1,mtp_rad_order
            do ip3=1,mtp_rad_order
             do b1=1,3
               do b2=1,3
                do b3=1,3
                  do b4=1,3
                    do b5=1,3
                      do ia=1,ia_n
                        do ix=1,3
                            config_desc(iconf)%force(icount_f,ja,ia,ix) = config_desc(iconf)%force(icount_f,ja,ia,ix) +  &
                             tmp_mtp6(ip1,b1,b2,b3,b4,b5)       *d_tmp_mtp2(ip2,b1,ia,ix)*tmp_mtp5(ip3,b2,b3,b4,b5) +    &
                           d_tmp_mtp6(ip1,b1,b2,b3,b4,b5,ia,ix) *  tmp_mtp2(ip2,b1)      *tmp_mtp5(ip3,b2,b3,b4,b5) +    &
                             tmp_mtp6(ip1,b1,b2,b3,b4,b5)       *  tmp_mtp2(ip2,b1)    *d_tmp_mtp5(ip3,b2,b3,b4,b5,ia,ix)
                          !config_desc(iconf)%force(icount_f,ja,ia,ix) = config_desc(iconf)%force(icount_f,ja,ia,ix)  +  a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                          !config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                        end do
                      end do
                    end do
                  end do
                end do
               end do
             end do
             icount_f = icount_f + 1
          end do
         end do
        end do
    end if





    if (desc_forces) then
      do ia=1,ia_n
       config_desc(iconf)%force(:,ja, 0,:) =  config_desc(iconf)%force(:,ja, 0,:) -  config_desc(iconf)%force(:,ja, ia,:)
      enddo  !ia from ia_m
    end if
    !NON NORMALIZED VERSION
    config_desc(iconf)%energy(:,ja) = config_desc(iconf)%energy(:,ja)*factor_ja!/dble(ia_n)
    !debug if (rangml==0) write(6,*) local_mtp_out(1,ja), ja
    if (desc_forces) config_desc(iconf)%force(:,ja, :,:) =  config_desc(iconf)%force(:,ja, :,:)*factor_ja!/dble(ia_n)
enddo !ja
 deallocate(xpnp)

return
end subroutine compute_mtp



subroutine  gen_dimension_for_mtp
use ml_in_ndm_module, ONLY: rangml, mtp_dim, mtp_poly_min, mtp_poly_max, mtp_rad_order
implicit none

if (mtp_poly_max <= mtp_poly_min ) then
  if (rangml==0) write(6,*) 'Error in setting MTP descriptor. mpt_poly_max sould be larger than mtp_poly_min'
  stop 'error in MPT in gen_dimension_for_mtp'
end if

mtp_rad_order = mtp_poly_max - mtp_poly_min + 1
   mtp_dim=mtp_rad_order + 5*mtp_rad_order**2 + 4*mtp_rad_order**3
!6  mtp_dim=mtp_rad_order + 5*mtp_rad_order**2 + 2*mtp_rad_order**3

if (rangml==0) write(6,*) 'ML: dimension of the descriptor space: ',mtp_dim

return
end subroutine gen_dimension_for_mtp


!real(kind=kind(1.d0)) function mtp_radial_tensor(type_mtp, r, ipo) result(func, d_func)
!integer :: type_mtp,ipo
!real(kind=kind(1.d0)) :: r, func, d_func
!
!  select  case (type_mtp)
!
!  case(1)
!    func= r**ipo
!    d_func= dble(ipo)*r**(ipo-2)
!  end select
!
!end function mtp_radial_tensor
