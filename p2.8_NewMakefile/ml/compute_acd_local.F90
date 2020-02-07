subroutine coord_acd_local(i)

  USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: imm,lperiod,bg,at,indi2,nvois
  use tab_imm_m, ONLY : xp,iwmax2
  use ml_in_ndm_module, ONLY: r_acd,at_acd,iwmax2_acd,indi2_acd

  implicit none

  integer,intent(in) :: i
  real(double), dimension(:,:), allocatable :: xpnp

  iwmax2_acd(:,i)=iwmax2(:)
  indi2_acd(1:nvois,i)=indi2(1:nvois)
  at_acd(:,:,i)=at
  ALLOCATE(xpnp(3,imm))
  if (lperiod) then
   xpnp(:,:)=xp(:,:)
  else
   call notperiod(xp,xpnp)
  end if
  call cryst_to_cart (imm, xpnp, bg, -1)
  r_acd(:,:,i)=xpnp

  deallocate(xpnp)

return
end subroutine coord_acd_local

subroutine compute_kernel_acd_local(ns_data,k_acd_out,distance_acd_out)

 USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: A2cm
  use tab_imm_m, ONLY : ityp
  use ml_in_ndm_module, ONLY: rangml,pi,r_cut,alpha_acd,kappa_acd,acd_fcut,ksi_ini,temp_ini,tau,mc_step,rotate,&
                              mconf,w2_rho,acd_weighted,massat,reject,sparsification_by_acd,max_data, &
                              data_im,data_natm,max_ntyp,r_acd,seed,at_acd,iwmax2_acd,indi2_acd
#if(PARAML)
!for curie      use mkl_service
  use mpi
  use mod_mpi_ml
#endif
  use set_limits

 implicit none

  integer,intent(in) :: ns_data
  double precision,dimension(ns_data,ns_data),intent(out) :: k_acd_out,distance_acd_out

  integer :: i_start_at,i_final_at
  integer :: i,k,j,ji,jk,it,iwj,iwj1,iwj2,iwk,iwk1,iwk2,ij,ik,count_mc,nstep,i_test,nk,ns_data_k
  double precision :: kappa,r,r1,r2,r3,r_ji,r1_ji,r2_ji,r3_ji,r_jk,r1_jk,r2_jk,r3_jk,fcut_ji,fcut_jk,fcut1_ji,fcut1_jk,fcut2_ji,fcut2_jk,fcut3_ji,fcut3_jk
  double precision,dimension(maxval(data_im)) :: norm_fcut_ji,norm_fcut_jk,norm_fcut1_ji,norm_fcut1_jk,norm_fcut2_ji,norm_fcut2_jk,norm_fcut3_ji,norm_fcut3_jk
  double precision,dimension(maxval(data_im)) :: local_norm_fcut_ji,local_norm_fcut1_ji,local_norm_fcut2_ji,local_norm_fcut3_ji
  double precision, dimension(3) :: dxp_ji,dxp1_ji,dxp2_ji,dxp3_ji,dxp_jk,dxp1_jk,dxp2_jk,dxp3_jk
  double precision :: factor
  double precision,dimension(3,maxval(data_im),ns_data) :: r_acd_new
  double precision,dimension(ns_data,ns_data) :: k_acd_tmp,local_k_acd_tmp
  double precision,dimension(max_ntyp,ns_data,ns_data) :: local_distance_acd_tmp1,local_distance_acd_tmp2,local_distance_acd_tmp3
  double precision,dimension(max_ntyp,ns_data,ns_data) :: distance_acd_tmp1,distance_acd_tmp2,distance_acd_tmp3
  double precision,dimension(3) :: angle,angle_old,rd_angle
  double precision :: rand,acc,gain,distance_old,k_old,temp,ksi,test_distance,test_distance_old,test_distance_tmp

  k_acd_out(:,:)=0d0
  k_acd_tmp(:,:)=0d0
  local_k_acd_tmp(:,:)=0d0
  local_distance_acd_tmp1(:,:,:)=0d0
  local_distance_acd_tmp2(:,:,:)=0d0
  local_distance_acd_tmp3(:,:,:)=0d0
  distance_acd_tmp1(:,:,:)=0d0
  distance_acd_tmp2(:,:,:)=0d0
  distance_acd_tmp3(:,:,:)=0d0
  distance_acd_out(:,:)=0d0

  kappa=(2d0*pi/alpha_acd)**(1.5)

if (sparsification_by_acd) then
   reject(1)=.false.
   ns_data_k=ns_data
else
   max_data=2
endif

do nk=1,max_data-1
  test_distance_old=0d0

  do i=1,ns_data
!     do k=1,i
     test_distance=0d0
     test_distance_tmp=1d0
     if (.not.sparsification_by_acd) ns_data_k=i
     do k=1,ns_data_k
        if (sparsification_by_acd.and.(k.eq.i)) cycle
        if (sparsification_by_acd.and.reject(k)) cycle
        if (sparsification_by_acd.and.(.not.reject(i))) cycle

        count_mc=0
        temp=temp_ini
        ksi=ksi_ini
!        call DLARNV( 3, seed, 3, angle )
        if (k.ne.i) then
            nstep=mc_step
        else
            nstep=0
        endif
        do it=0,nstep
           if (it==0) then
              angle=0d0
              r_acd_new(1:3,:,k)=r_acd(1:3,:,k)
           else
              angle_old=angle
              call DLARNV( 3, seed, 3, rd_angle )
              angle=angle+ksi*rd_angle
              do jk=1,data_im(k)
                 call rotate(angle,r_acd(1:3,jk,k),r_acd_new(1:3,jk,k))
              enddo
           endif

           local_k_acd_tmp(k,i)=0d0
           local_distance_acd_tmp1(:,k,i)=0d0
           local_distance_acd_tmp2(:,k,i)=0d0
           local_distance_acd_tmp3(:,k,i)=0d0
           do j=1,max_ntyp

              local_norm_fcut1_ji(:)=0d0
              local_norm_fcut2_ji(:)=0d0
              call set_limit_for_atoms(rangml,data_natm(j,i),i_start_at,i_final_at)
              if (i_start_at==1)  iwj2=0
              if (i_start_at > 1) iwj2=iwmax2_acd(i_start_at-1,i)

              do ji=i_start_at,i_final_at
                 iwj1=iwj2+1
                 iwj2=iwmax2_acd(ji,i)

                 iwk2=0
                 do jk=1,data_natm(j,i)
                    iwk1=iwk2+1
                    iwk2=iwmax2_acd(jk,i)

                    local_norm_fcut1_ji(ji)=0d0
                    do iwj=iwj1,iwj2
                       ij=indi2_acd(iwj,i)
                       if (ji==ij) cycle

                       dxp1_ji(1:3) = r_acd(1:3,ij,i) - r_acd(1:3,ji,i)
                       WHERE ( (dxp1_ji.GT.0.5d0).OR.(dxp1_ji.LT.-0.5d0) )
                         dxp1_ji(1:3) = dxp1_ji(1:3) - Dble(Nint(dxp1_ji(1:3)))
                       END WHERE
                       dxp1_ji = MatMul(at_acd(:,:,i),dxp1_ji)/A2cm
                       r1_ji = dsqrt(sum(dxp1_ji(1:3)**2))
                       if (r1_ji.gt.r_cut) cycle
                       if (acd_fcut) then
                          fcut1_ji=0.5d0*(cos(pi*r1_ji/r_cut)+1d0)
                       else
                          fcut1_ji=1d0
                       endif
                       local_norm_fcut1_ji(ji)=local_norm_fcut1_ji(ji)+fcut1_ji

                       norm_fcut1_jk(jk)=0d0
                       do iwk=iwk1,iwk2
                          ik=indi2_acd(iwk,i)
                          if (jk==ik) cycle

                          dxp1_jk(1:3) = r_acd(1:3,ik,i) - r_acd(1:3,jk,i)
                          WHERE ( (dxp1_jk.GT.0.5d0).OR.(dxp1_jk.LT.-0.5d0) )
                            dxp1_jk(1:3) = dxp1_jk(1:3) - Dble(Nint(dxp1_jk(1:3)))
                          END WHERE
                          dxp1_jk = MatMul(at_acd(:,:,i),dxp1_jk)/A2cm
                          r1_jk = dsqrt(sum(dxp1_jk(1:3)**2))
                          if (r1_jk.gt.r_cut) cycle
                          if (acd_fcut) then
                             fcut1_jk=0.5d0*(cos(pi*r1_jk/r_cut)+1d0)
                          else
                             fcut1_jk=1d0
                          endif
                          norm_fcut1_jk(jk)=norm_fcut1_jk(jk)+fcut1_jk

                          if (acd_weighted) then
                             factor=massat(ityp(j))/sum(massat)
                          else
                             factor=1.d0
                          endif
!                          r1 = dsqrt(sum((dxp1_jk(1:3)-dxp1_ji(1:3))**2))
                          r1 = r1_jk - r1_ji ! densité radiale
                          local_distance_acd_tmp1(j,k,i) = local_distance_acd_tmp1(j,k,i) + dexp(-0.5*alpha_acd*(r1**2)) * fcut1_ji * fcut1_jk * factor
                       enddo
                    enddo
                 enddo !jk

                 iwk2=0
                 do jk=1,data_natm(j,k)
                    iwk1=iwk2+1
                    iwk2=iwmax2_acd(jk,k)

                    local_norm_fcut2_ji(ji)=0d0
                    do iwj=iwj1,iwj2
                       ij=indi2_acd(iwj,i)
                       if (ji==ij) cycle

                       dxp2_ji(1:3) = r_acd(1:3,ij,i) - r_acd(1:3,ji,i)
                       WHERE ( (dxp2_ji.GT.0.5d0).OR.(dxp2_ji.LT.-0.5d0) )
                         dxp2_ji(1:3) = dxp2_ji(1:3) - Dble(Nint(dxp2_ji(1:3)))
                       END WHERE
                       dxp2_ji = MatMul(at_acd(:,:,i),dxp2_ji)/A2cm
                       r2_ji = dsqrt(sum(dxp2_ji(1:3)**2))
                       if (r2_ji.gt.r_cut) cycle
                       if (acd_fcut) then
                          fcut2_ji=0.5d0*(cos(pi*r2_ji/r_cut)+1d0)
                       else
                          fcut2_ji=1d0
                       endif
                       local_norm_fcut2_ji(ji)=local_norm_fcut2_ji(ji)+fcut2_ji

                       norm_fcut2_jk(jk)=0d0
                       do iwk=iwk1,iwk2
                          ik=indi2_acd(iwk,k)
                          if (jk==ik) cycle

                          dxp2_jk(1:3) = r_acd_new(1:3,ik,k) - r_acd_new(1:3,jk,k)
                          WHERE ( (dxp2_jk.GT.0.5d0).OR.(dxp2_jk.LT.-0.5d0) )
                            dxp2_jk(1:3) = dxp2_jk(1:3) - Dble(Nint(dxp2_jk(1:3)))
                          END WHERE
                          dxp2_jk = MatMul(at_acd(:,:,k),dxp2_jk)/A2cm
                          r2_jk = dsqrt(sum(dxp2_jk(1:3)**2))
                          if (r2_jk.gt.r_cut) cycle
                          if (acd_fcut) then
                             fcut2_jk=0.5d0*(cos(pi*r2_jk/r_cut)+1d0)
                          else
                             fcut2_jk=1d0
                          endif
                          norm_fcut2_jk(jk)=norm_fcut2_jk(jk)+fcut2_jk

                          if (acd_weighted) then
                             factor=massat(ityp(j))/sum(massat)
                          else
                             factor=1.d0
                          endif
!                          r2 = dsqrt(sum((dxp2_jk(1:3)-dxp2_ji(1:3))**2))
                          r2 = r2_jk - r2_ji ! densité radiale
                          local_distance_acd_tmp2(j,k,i) = local_distance_acd_tmp2(j,k,i) + dexp(-0.5*alpha_acd*(r2**2)) * fcut2_ji * fcut2_jk * factor
                       enddo
                    enddo
                 enddo !jk
#if (PARAML)
#else
                 norm_fcut1_ji(ji)=local_norm_fcut1_ji(ji)
                 norm_fcut2_ji(ji)=local_norm_fcut2_ji(ji)
#endif
              enddo !ji

#if (PARAML)
              call MPI_ALLREDUCE(local_norm_fcut1_ji,norm_fcut1_ji,data_natm(j,i),MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
              call MPI_ALLREDUCE(local_norm_fcut2_ji,norm_fcut2_ji,data_natm(j,i),MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
              if (((sum(norm_fcut1_ji)==0d0).or.(sum(norm_fcut1_jk)==0d0)).and.(rangml==0)) then
                 write(*,*)"sum(norm_fcut1_ji) or sum(norm_fcut1_jk) = 0"
                 write(*,*)"Take larger r_cut"
                 stop
              elseif (((sum(norm_fcut2_ji)==0d0).or.(sum(norm_fcut2_jk)==0d0)).and.(rangml==0)) then
                 write(*,*)"sum(norm_fcut2_ji) or sum(norm_fcut2_jk) = 0"
                 write(*,*)"Take larger r_cut"
                 stop
              endif
              local_distance_acd_tmp1(j,k,i) = local_distance_acd_tmp1(j,k,i) / (kappa*sum(norm_fcut1_ji)*sum(norm_fcut1_jk))
              local_distance_acd_tmp2(j,k,i) = local_distance_acd_tmp2(j,k,i) / (kappa*sum(norm_fcut2_ji)*sum(norm_fcut2_jk))
              call MPI_REDUCE(local_distance_acd_tmp1(j,k,i),distance_acd_tmp1(j,k,i),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,codeml)
              call MPI_REDUCE(local_distance_acd_tmp2(j,k,i),distance_acd_tmp2(j,k,i),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,codeml)
#else
              if (((sum(norm_fcut1_ji)==0d0).or.(sum(norm_fcut1_jk)==0d0)).and.(rangml==0)) then
                 write(*,*)"sum(norm_fcut1_ji) or sum(norm_fcut1_jk) = 0"
                 write(*,*)"Take larger r_cut"
                 stop
              elseif (((sum(norm_fcut2_ji)==0d0).or.(sum(norm_fcut2_jk)==0d0)).and.(rangml==0)) then
                 write(*,*)"sum(norm_fcut2_ji) or sum(norm_fcut2_jk) = 0"
                 write(*,*)"Take larger r_cut"
                 stop
              endif
              distance_acd_tmp1(j,k,i) = local_distance_acd_tmp1(j,k,i) / (kappa*sum(norm_fcut1_ji)*sum(norm_fcut1_jk))
              distance_acd_tmp2(j,k,i) = local_distance_acd_tmp2(j,k,i) / (kappa*sum(norm_fcut2_ji)*sum(norm_fcut2_jk))
#endif
              local_norm_fcut3_ji(:)=0d0
              call set_limit_for_atoms(rangml,data_natm(j,k),i_start_at,i_final_at)
              if (i_start_at==1)  iwj2=0
              if (i_start_at > 1) iwj2=iwmax2_acd(i_start_at-1,k)

              do ji=i_start_at,i_final_at
                 iwj1=iwj2+1
                 iwj2=iwmax2_acd(ji,k)

                 iwk2=0
                 do jk=1,data_natm(j,k)
                    iwk1=iwk2+1
                    iwk2=iwmax2_acd(jk,k)

                    local_norm_fcut3_ji(ji)=0d0
                    do iwj=iwj1,iwj2
                       ij=indi2_acd(iwj,k)
                       if (ji==ij) cycle

                       dxp3_ji(1:3) = r_acd_new(1:3,ij,k) - r_acd_new(1:3,ji,k)
                       WHERE ( (dxp3_ji.GT.0.5d0).OR.(dxp3_ji.LT.-0.5d0) )
                         dxp3_ji(1:3) = dxp3_ji(1:3) - Dble(Nint(dxp3_ji(1:3)))
                       END WHERE
                       dxp3_ji = MatMul(at_acd(:,:,k),dxp3_ji)/A2cm
                       r3_ji = dsqrt(sum(dxp3_ji(1:3)**2))
                       if (r3_ji.gt.r_cut) cycle
                       if (acd_fcut) then
                          fcut3_ji=0.5d0*(cos(pi*r3_ji/r_cut)+1d0)
                       else
                          fcut3_ji=1d0
                       endif
                       local_norm_fcut3_ji(ji)=local_norm_fcut3_ji(ji)+fcut3_ji

                       norm_fcut3_jk(jk)=0d0
                       do iwk=iwk1,iwk2
                          ik=indi2_acd(iwk,k)
                          if (jk==ik) cycle

                          dxp3_jk(1:3) = r_acd_new(1:3,ik,k) - r_acd_new(1:3,jk,k)
                          WHERE ( (dxp3_jk.GT.0.5d0).OR.(dxp3_jk.LT.-0.5d0) )
                            dxp3_jk(1:3) = dxp3_jk(1:3) - Dble(Nint(dxp3_jk(1:3)))
                          END WHERE
                          dxp3_jk = MatMul(at_acd(:,:,k),dxp3_jk)/A2cm
                          r3_jk = dsqrt(sum(dxp3_jk(1:3)**2))
                          if (r3_jk.gt.r_cut) cycle
                          if (acd_fcut) then
                             fcut3_jk=0.5d0*(cos(pi*r3_jk/r_cut)+1d0)
                          else
                             fcut3_jk=1d0
                          endif
                          norm_fcut3_jk(jk)=norm_fcut3_jk(jk)+fcut3_jk

                          if (acd_weighted) then
                             factor=massat(ityp(j))/sum(massat)
                          else
                             factor=1.d0
                          endif
!                          r3 = dsqrt(sum((dxp3_jk(1:3)-dxp3_ji(1:3))**2))
                          r3 = r3_jk - r3_ji ! densité radiale
                          local_distance_acd_tmp3(j,k,i) = local_distance_acd_tmp3(j,k,i) + dexp(-0.5*alpha_acd*(r3**2)) * fcut3_ji * fcut3_jk * factor
                       enddo
                    enddo

                 enddo !jk
#if (PARAML)
#else
                 norm_fcut3_ji(ji)=local_norm_fcut3_ji(ji)
#endif
              enddo !ji

#if (PARAML)
             call MPI_ALLREDUCE(local_norm_fcut3_ji,norm_fcut3_ji,data_natm(j,k),MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
             if (((sum(norm_fcut3_ji)==0d0).or.(sum(norm_fcut3_jk)==0d0)).and.(rangml==0)) then
                write(*,*)"sum(norm_fcut3_ji) or sum(norm_fcut3_jk) = 0"
                write(*,*)"Take larger r_cut"
                stop
             endif
             local_distance_acd_tmp3(j,k,i) = local_distance_acd_tmp3(j,k,i) / (kappa*sum(norm_fcut3_ji)*sum(norm_fcut3_jk))
             call MPI_REDUCE(local_distance_acd_tmp3(j,k,i),distance_acd_tmp3(j,k,i),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,codeml)
#else
             if (((sum(norm_fcut3_ji)==0d0).or.(sum(norm_fcut3_jk)==0d0)).and.(rangml==0)) then
                write(*,*)"sum(norm_fcut3_ji) or sum(norm_fcut3_jk) = 0"
                write(*,*)"Take larger r_cut"
                stop
             endif
             distance_acd_tmp3(j,k,i) = local_distance_acd_tmp3(j,k,i) / (kappa*sum(norm_fcut3_ji)*sum(norm_fcut3_jk))
#endif
           enddo !j

           local_norm_fcut_ji(:)=0d0
           call set_limit_for_atoms(rangml,data_im(i),i_start_at,i_final_at)
           if (i_start_at==1)  iwj2=0
           if (i_start_at > 1) iwj2=iwmax2_acd(i_start_at-1,i)

           do ji=i_start_at,i_final_at
              iwj1=iwj2+1
              iwj2=iwmax2_acd(ji,i)

              iwk2=0
              do jk=1,data_im(k)
                 iwk1=iwk2+1
                 iwk2=iwmax2_acd(jk,k)

                 local_norm_fcut_ji(ji)=0d0
                 do iwj=iwj1,iwj2
                    ij=indi2_acd(iwj,i)
                    if (ji==ij) cycle

                    dxp_ji(1:3) = r_acd(1:3,ij,i) - r_acd(1:3,ji,i)
                    WHERE ( (dxp_ji.GT.0.5d0).OR.(dxp_ji.LT.-0.5d0) )
                      dxp_ji(1:3) = dxp_ji(1:3) - Dble(Nint(dxp_ji(1:3)))
                    END WHERE
                    dxp_ji = MatMul(at_acd(:,:,i),dxp_ji)/A2cm
                    r_ji = dsqrt(sum(dxp_ji(1:3)**2))
                    if (r_ji.gt.r_cut) cycle
                    if (acd_fcut) then
                       fcut_ji=0.5d0*(cos(pi*r_ji/r_cut)+1d0)
                    else
                       fcut_ji=1d0
                    endif
                    local_norm_fcut_ji(ji)=local_norm_fcut_ji(ji)+fcut_ji

                    norm_fcut_jk(jk)=0d0
                    do iwk=iwk1,iwk2
                       ik=indi2_acd(iwk,k)
                       if (jk==ik) cycle

                       dxp_jk(1:3) = r_acd_new(1:3,ik,k) - r_acd_new(1:3,jk,k)
                       WHERE ( (dxp_jk.GT.0.5d0).OR.(dxp_jk.LT.-0.5d0) )
                         dxp_jk(1:3) = dxp_jk(1:3) - Dble(Nint(dxp_jk(1:3)))
                       END WHERE
                       dxp_jk = MatMul(at_acd(:,:,k),dxp_jk)/A2cm
                       r_jk = dsqrt(sum(dxp_jk(1:3)**2))
                       if (r_jk.gt.r_cut) cycle
                       if (acd_fcut) then
                          fcut_jk=0.5d0*(cos(pi*r_jk/r_cut)+1d0)
                       else
                          fcut_jk=1d0
                       endif
                       norm_fcut_jk(jk)=norm_fcut_jk(jk)+fcut_jk

                       if (acd_weighted) then
                          factor=w2_rho(massat(ityp(ji)),massat(ityp(jk)))
                       else
                          factor=1.d0
                       endif
!                       r = dsqrt(sum((dxp_jk(1:3)-dxp_ji(1:3))**2))
                       r = r_jk - r_ji ! densité radiale
                       local_k_acd_tmp(k,i) = local_k_acd_tmp(k,i) + dexp(-0.5*alpha_acd*(r**2)) * fcut_ji * fcut_jk * factor

                    enddo
                 enddo
              enddo !jk
#if (PARAML)
#else
              norm_fcut_ji(ji)=local_norm_fcut_ji(ji)
#endif
           enddo !ji
#if (PARAML)
           call MPI_ALLREDUCE(local_norm_fcut_ji,norm_fcut_ji,data_im(i),MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
           if (((sum(norm_fcut_ji)==0d0).or.(sum(norm_fcut_jk)==0d0)).and.(rangml==0)) then
              write(*,*)"sum(norm_fcut_ji) or sum(norm_fcut_jk) = 0"
              write(*,*)"Take larger r_cut"
              stop
           endif
           local_k_acd_tmp(k,i) = local_k_acd_tmp(k,i) / (kappa*sum(norm_fcut_ji)*sum(norm_fcut_jk))
           call MPI_REDUCE(local_k_acd_tmp(k,i),k_acd_tmp(k,i),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,codeml)
#else
              if (((sum(norm_fcut_ji)==0d0).or.(sum(norm_fcut_jk)==0d0)).and.(rangml==0)) then
                 write(*,*)"sum(norm_fcut_ji) or sum(norm_fcut_jk) = 0"
                 write(*,*)"Take larger r_cut"
                 stop
              endif
              k_acd_tmp(k,i) = local_k_acd_tmp(k,i) / (kappa*sum(norm_fcut_ji)*sum(norm_fcut_jk))
#endif

!           if (rangml==0) then

              distance_acd_out(k,i)=0d0
              do j=1,max_ntyp
                 distance_acd_out(k,i) = distance_acd_out(k,i) + (distance_acd_tmp1(j,k,i) - 2d0*distance_acd_tmp2(j,k,i) + distance_acd_tmp3(j,k,i))
              enddo
              distance_acd_out(k,i) = dsqrt(abs(distance_acd_out(k,i)))

              if (it==0) then
                 distance_old=distance_acd_out(k,i)
                 k_old=k_acd_tmp(k,i)
              else
                 gain=distance_acd_out(k,i)-distance_old
                 if (gain>0d0) then
                    acc=min(1d0,dexp(-gain/temp))
                    call DLARNV( 1, seed, 1, rand )
                    if (rand<acc) then
                       count_mc=count_mc+1
                    else
                       angle=angle_old
                       distance_acd_out(k,i)=distance_old
                       k_acd_tmp(k,i)=k_old
                    endif
                 else
                    distance_old=distance_acd_out(k,i)
                    k_old=k_acd_tmp(k,i)
                 endif
                 temp=tau*temp
                 ksi=tau*ksi
              endif

!              if (sparsification_by_acd.and.(distance_acd_out(k,i)<acd_threshold).and.(k.ne.i)) then
!                 reject(i)=.true.
!              endif
              if (distance_acd_out(k,i)==0d0) exit

!           endif

        enddo !it

        test_distance=test_distance_tmp*distance_acd_out(k,i)
        test_distance_tmp=test_distance

     enddo !k

     if (test_distance>test_distance_old) then
        i_test=i
        test_distance_old=test_distance
     endif

  enddo !i

  if (sparsification_by_acd) reject(i_test)=.false.

enddo ! nk

  if (rangml==0) then
     do i=1,ns_data
        do k=1,i
           k_acd_out(k,i) = ( k_acd_tmp(k,i) / dsqrt( k_acd_tmp(k,k) * k_acd_tmp(i,i) ) )**kappa_acd
        enddo
     enddo
  endif


 deallocate(iwmax2_acd,indi2_acd,at_acd)

return
end subroutine compute_kernel_acd_local
