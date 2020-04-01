module compute_acd_mod
        use notperiod_mod
        implicit none
        contains
subroutine compute_kernel_acd(ns_data,k_acd_out,distance_acd_out)

 USE T_kind_param_m, ONLY:  double
  use tab_imm_m, ONLY : ityp
  use ml_in_ndm_module, ONLY: rangml,pi,r_cut,alpha_acd,kappa_acd,ksi_ini,temp_ini,tau,mc_step,rotate,&
                              mconf,w2_rho,acd_weighted,acd_fcut,massat,reject,sparsification_by_acd,max_data, &
                              data_im,data_natm,max_ntyp,r_acd,seed
#ifdef PARAML 
!for curie      use mkl_service
  use mpi
  use mod_mpi_ml 
#endif
  use set_limits

 implicit none

  integer,intent(in) :: ns_data
  double precision,dimension(ns_data,ns_data),intent(out) :: k_acd_out,distance_acd_out

  integer :: i_start_at,i_final_at
  integer :: i,k,j,ji,jk,it,count_mc,nstep,i_test,nk,ns_data_k
  double precision :: kappa,r,r1,r2,r3,fcut_ji,fcut_jk,fcut1_ji,fcut1_jk,fcut2_ji,fcut2_jk,fcut3_ji,fcut3_jk
  double precision :: norm_fcut_ji,norm_fcut_jk,norm_fcut1_ji,norm_fcut1_jk,norm_fcut2_ji,norm_fcut2_jk,norm_fcut3_ji,norm_fcut3_jk
  double precision :: local_norm_fcut_ji,local_norm_fcut1_ji,local_norm_fcut2_ji,local_norm_fcut3_ji
  double precision, dimension(3) :: dxp,dxp1,dxp2,dxp3
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
   
              local_norm_fcut1_ji=0d0
              local_norm_fcut2_ji=0d0
   
              call set_limit_for_atoms(rangml,data_natm(j,i),i_start_at,i_final_at)
              do ji=i_start_at,i_final_at
                 if (dsqrt(sum(r_acd(1:3,ji,i)**2)).gt.r_cut) cycle
                 if (acd_fcut) then
                    fcut1_ji=0.5d0*(cos(pi*dsqrt(sum(r_acd(1:3,ji,i)**2))/r_cut)+1d0)
                    fcut2_ji=0.5d0*(cos(pi*dsqrt(sum(r_acd(1:3,ji,i)**2))/r_cut)+1d0)
                 else
                    fcut1_ji=1d0
                    fcut2_ji=1d0
                 endif
      
                 if (rangml==0) norm_fcut1_jk=0d0
                 do jk=1,data_natm(j,i)
   
                    if (dsqrt(sum(r_acd(1:3,jk,i)**2)).gt.r_cut) cycle
                    if (acd_fcut) then
                       fcut1_jk=0.5d0*(cos(pi*dsqrt(sum(r_acd(1:3,jk,k)**2))/r_cut)+1d0)
                    else
                       fcut1_jk=1d0
                    endif
   
                    if (acd_weighted) then
                       factor=massat(ityp(j))/sum(massat)
                    else
                       factor=1.d0
                    endif
                    dxp1(1:3) = r_acd(1:3,jk,i) - r_acd(1:3,ji,i)
                    r1 = dsqrt(sum(dxp1(1:3)**2))
                    local_distance_acd_tmp1(j,k,i) = local_distance_acd_tmp1(j,k,i) + dexp(-0.5*alpha_acd*(r1**2)) * fcut1_ji * fcut1_jk * factor
                    if (rangml==0) norm_fcut1_jk=norm_fcut1_jk+fcut1_jk
                 enddo !jk
      
                 if (rangml==0) norm_fcut2_jk=0d0
                 do jk=1,data_natm(j,k)
                    if (dsqrt(sum(r_acd_new(1:3,jk,k)**2)).gt.r_cut) cycle
                    if (acd_fcut) then
                       fcut2_jk=0.5d0*(cos(pi*dsqrt(sum(r_acd_new(1:3,jk,k)**2))/r_cut)+1d0)
                    else
                       fcut2_jk=1d0
                    endif
               
                    if (acd_weighted) then
                       factor=massat(ityp(j))/sum(massat)
                    else
                       factor=1.d0
                    endif
          
                    dxp2(1:3) = r_acd_new(1:3,jk,k) - r_acd(1:3,ji,i)
                    r2 = dsqrt(sum(dxp2(1:3)**2))
                    local_distance_acd_tmp2(j,k,i) = local_distance_acd_tmp2(j,k,i) + dexp(-0.5*alpha_acd*(r2**2)) * fcut2_ji * fcut2_jk * factor
                    if (rangml==0) norm_fcut2_jk=norm_fcut2_jk+fcut2_jk
   
                 enddo !jk
                 local_norm_fcut1_ji=local_norm_fcut1_ji+fcut1_ji
                 local_norm_fcut2_ji=local_norm_fcut2_ji+fcut2_ji
#ifdef PARAML
#else 
                 norm_fcut1_ji=local_norm_fcut1_ji
                 norm_fcut2_ji=local_norm_fcut2_ji
#endif
              enddo !ji           
      
#ifdef PARAML 
              call MPI_ALLREDUCE(local_norm_fcut1_ji,norm_fcut1_ji,1,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
              call MPI_ALLREDUCE(local_norm_fcut2_ji,norm_fcut2_ji,1,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
              call MPI_BCAST(norm_fcut1_jk,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,codeml)
              call MPI_BCAST(norm_fcut2_jk,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,codeml)
              if (((norm_fcut1_ji==0d0).or.(norm_fcut1_jk==0d0)).and.(rangml==0)) then
                 write(*,*)"norm_fcut1_ji or norm_fcut1_jk = 0"
                 write(*,*)"Take larger r_cut"
                 stop
              elseif (((norm_fcut2_ji==0d0).or.(norm_fcut2_jk==0d0)).and.(rangml==0)) then
                 write(*,*)"norm_fcut2_ji or norm_fcut2_jk = 0"
                 write(*,*)"Take larger r_cut"
                 stop
              endif
              local_distance_acd_tmp1(j,k,i) = local_distance_acd_tmp1(j,k,i) / (kappa*norm_fcut1_ji*norm_fcut1_jk)
              local_distance_acd_tmp2(j,k,i) = local_distance_acd_tmp2(j,k,i) / (kappa*norm_fcut2_ji*norm_fcut2_jk)
              call MPI_REDUCE(local_distance_acd_tmp1(j,k,i),distance_acd_tmp1(j,k,i),1,MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,codeml)
              call MPI_REDUCE(local_distance_acd_tmp2(j,k,i),distance_acd_tmp2(j,k,i),1,MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,codeml)
#else 
              if (((norm_fcut1_ji==0d0).or.(norm_fcut1_jk==0d0)).and.(rangml==0)) then
                 write(*,*)"norm_fcut1_ji or norm_fcut1_jk = 0"
                 write(*,*)"Take larger r_cut"
                 stop
              elseif (((norm_fcut2_ji==0d0).or.(norm_fcut2_jk==0d0)).and.(rangml==0)) then
                 write(*,*)"norm_fcut2_ji or norm_fcut2_jk = 0"
                 write(*,*)"Take larger r_cut"
                 stop
              endif
              local_distance_acd_tmp1(j,k,i) = local_distance_acd_tmp1(j,k,i) / (kappa*norm_fcut1_ji*norm_fcut1_jk)
              local_distance_acd_tmp2(j,k,i) = local_distance_acd_tmp2(j,k,i) / (kappa*norm_fcut2_ji*norm_fcut2_jk)
              distance_acd_tmp1(j,k,i) = local_distance_acd_tmp1(j,k,i)
              distance_acd_tmp2(j,k,i) = local_distance_acd_tmp2(j,k,i)
#endif
              local_norm_fcut3_ji=0d0
              call set_limit_for_atoms(rangml,data_natm(j,k),i_start_at,i_final_at)
              do ji=i_start_at,i_final_at
                 if (dsqrt(sum(r_acd_new(1:3,ji,k)**2)).gt.r_cut) cycle
                 if (acd_fcut) then
                    fcut3_ji=0.5d0*(cos(pi*dsqrt(sum(r_acd_new(1:3,ji,k)**2))/r_cut)+1d0)
                 else
                    fcut3_ji=1d0
                 endif
                 if (rangml==0) norm_fcut3_jk=0d0
                 do jk=1,data_natm(j,k)
   
                    if (dsqrt(sum(r_acd_new(1:3,jk,k)**2)).gt.r_cut) cycle
                    if (acd_fcut) then
                       fcut3_jk=0.5d0*(cos(pi*dsqrt(sum(r_acd_new(1:3,jk,k)**2))/r_cut)+1d0)
                    else
                       fcut3_jk=1d0
                    endif
   
                    if (acd_weighted) then
                       factor=massat(ityp(j))/sum(massat)
                    else
                       factor=1.d0
                    endif
                    dxp3(1:3) = r_acd_new(1:3,jk,k) - r_acd_new(1:3,ji,k)
                    r3 = dsqrt(sum(dxp3(1:3)**2))
      
                    local_distance_acd_tmp3(j,k,i) = local_distance_acd_tmp3(j,k,i) + dexp(-0.5*alpha_acd*(r3**2)) * fcut3_ji * fcut3_jk * factor
                    if (rangml==0) norm_fcut3_jk=norm_fcut3_jk+fcut3_jk
      
                 enddo !jk
                 local_norm_fcut3_ji=local_norm_fcut3_ji+fcut3_ji
#ifdef PARAML 
#else 
                 norm_fcut3_ji=local_norm_fcut3_ji
#endif
              enddo !ji           
      
#ifdef PARAML 
             call MPI_ALLREDUCE(local_norm_fcut3_ji,norm_fcut3_ji,1,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
             call MPI_BCAST(norm_fcut3_jk,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,codeml)
             if (((norm_fcut3_ji==0d0).or.(norm_fcut3_jk==0d0)).and.(rangml==0)) then
                write(*,*)"norm_fcut3_ji or norm_fcut3_jk = 0"
                write(*,*)"Take larger r_cut"
                stop
             endif
             local_distance_acd_tmp3(j,k,i) = local_distance_acd_tmp3(j,k,i) / (kappa*norm_fcut3_ji*norm_fcut3_jk)
             call MPI_REDUCE(local_distance_acd_tmp3(j,k,i),distance_acd_tmp3(j,k,i),1,MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,codeml)
#else 
             if (((norm_fcut3_ji==0d0).or.(norm_fcut3_jk==0d0)).and.(rangml==0)) then
                write(*,*)"norm_fcut3_ji or norm_fcut3_jk = 0"
                write(*,*)"Take larger r_cut"
                stop
             endif
             local_distance_acd_tmp3(j,k,i) = local_distance_acd_tmp3(j,k,i) / (kappa*norm_fcut3_ji*norm_fcut3_jk)
             distance_acd_tmp3(j,k,i) = local_distance_acd_tmp3(j,k,i)
#endif
           enddo !j
   
           local_norm_fcut_ji=0d0
           call set_limit_for_atoms(rangml,data_im(i),i_start_at,i_final_at)
           do ji=i_start_at,i_final_at
              if (dsqrt(sum(r_acd(1:3,ji,i)**2)).gt.r_cut) cycle
              if (acd_fcut) then
                 fcut_ji=0.5d0*(cos(pi*dsqrt(sum(r_acd(1:3,ji,i)**2))/r_cut)+1d0)
              else
                 fcut_ji=1d0
              endif
   
              if (rangml==0) norm_fcut_jk=0d0
              do jk=1,data_im(k)
   
                 if (dsqrt(sum(r_acd_new(1:3,jk,k)**2)).gt.r_cut) cycle
                 if (acd_fcut) then
                    fcut_jk=0.5d0*(cos(pi*dsqrt(sum(r_acd_new(1:3,jk,k)**2))/r_cut)+1d0)
                 else
                    fcut_jk=1d0
                 endif
            
                 if (acd_weighted) then
                    factor=w2_rho(massat(ityp(ji)),massat(ityp(jk)))
                 else
                    factor=1.d0
                 endif
       
                 dxp(1:3) = r_acd_new(1:3,jk,k) - r_acd(1:3,ji,i)
                 r = dsqrt(sum(dxp(1:3)**2))
                 local_k_acd_tmp(k,i) = local_k_acd_tmp(k,i) + dexp(-0.5*alpha_acd*(r**2)) * fcut_ji * fcut_jk * factor
                 if (rangml==0) norm_fcut_jk=norm_fcut_jk+fcut_jk

              enddo !jk
              local_norm_fcut_ji=local_norm_fcut_ji+fcut_ji
#ifdef PARAML 
#else 
              norm_fcut_ji=local_norm_fcut_ji
#endif
           enddo !ji           
#ifdef PARAML 
           call MPI_ALLREDUCE(local_norm_fcut_ji,norm_fcut_ji,1,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
           call MPI_BCAST(norm_fcut_jk,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,codeml)
           if (((norm_fcut_ji==0d0).or.(norm_fcut_jk==0d0)).and.(rangml==0)) then
              write(*,*)"norm_fcut_ji or norm_fcut_jk = 0"
              write(*,*)"Take larger r_cut"
              stop
           endif
           local_k_acd_tmp(k,i) = local_k_acd_tmp(k,i) / (norm_fcut_ji*norm_fcut_jk)
           call MPI_REDUCE(local_k_acd_tmp(k,i),k_acd_tmp(k,i),1,MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,codeml)
#else 
           if (((norm_fcut_ji==0d0).or.(norm_fcut_jk==0d0)).and.(rangml==0)) then
              write(*,*)"norm_fcut_ji or norm_fcut_jk = 0"
              write(*,*)"Take larger r_cut"
              stop
           endif
           k_acd_tmp(k,i) = local_k_acd_tmp(k,i) / (norm_fcut_ji*norm_fcut_jk)
#endif
   
!           if (rangml==0) then

              distance_acd_out(k,i)=0d0
              do j=1,max_ntyp
                 distance_acd_out(k,i)=distance_acd_out(k,i)+(distance_acd_tmp1(j,k,i)-2d0*distance_acd_tmp2(j,k,i)+distance_acd_tmp3(j,k,i))
              enddo
              distance_acd_out(k,i)=dsqrt(abs(distance_acd_out(k,i)))
  
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
  
!              if (sparsification_by_acd.and.(distance_acd_out(k,i)>acd_threshold).and.(k.ne.i)) then
!                 reject(i)=.false.
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
  
return
end subroutine compute_kernel_acd
end module
