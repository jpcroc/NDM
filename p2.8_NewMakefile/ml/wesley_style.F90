
subroutine coord_soap(i)

 USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: imm,lperiod,bg,at,indi2,nvois
  use tab_imm_m, ONLY : iwmax2,xp
  use ml_in_ndm_module, ONLY: r_soap,iwmax2_soap,indi2_soap,at_soap

  implicit none

  integer,intent(in) :: i
  real(double), dimension(:,:), allocatable :: xpnp

  iwmax2_soap(:,i)=iwmax2(:)
  indi2_soap(1:nvois,i)=indi2(1:nvois)
  at_soap(:,:,i)=at
  ALLOCATE(xpnp(3,imm))
  if (lperiod) then
   xpnp(:,:)=xp(:,:)
  else
   call notperiod(xp,xpnp)
  end if
  call cryst_to_cart (imm, xpnp, bg, -1)
  r_soap(:,:,i)=xpnp

return
end subroutine coord_soap

subroutine compute_kernel_soap(k_soap_out)

 USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: A2cm
  use tab_imm_m, ONLY : ityp
  use angular_functions, only : spherical_harm
  use ml_in_ndm_module, ONLY: rangml,ns_data,pi,r_cut,l_max,alpha_soap,kappa_acd,n_soap,&
                              mconf,w2_rho,weighted,massat, &
                              data_im,r_soap,iwmax2_soap,indi2_soap,at_soap,lsoap_fcut, cg_vector
#if(PARAML)
!for curie      use mkl_service
  use mpi
  use mod_mpi_ml
#endif
  use set_limits

 implicit none

  double precision,dimension(sum(data_im),sum(data_im)),intent(out) :: k_soap_out

  integer :: i_start_at,i_final_at
  integer :: i,k,ji,jk,l,m,mp,l1,l2,m1,m1p,iwj,iwj1,iwj2,iwk,iwk1,iwk2,ij,ik
  double precision :: r_ji,r_jk,fcut_ji,fcut_jk,arg_bessel
  double precision :: norm_fcut_ji,norm_fcut_jk
  double precision, dimension(3) :: dxp_ji,dxp_jk
  double precision,dimension(0:l_max) :: ms_bessel
  double complex,dimension(-l_max:l_max,-l_max:l_max,0:l_max) :: c
  double complex :: sum_bi
  double precision :: factor_ji,factor_jk
  double precision,dimension(sum(data_im),sum(data_im)) :: k_soap_tmp,local_k_soap_tmp


!  if (allocated(m_typ)) deallocate(m_typ)
!  allocate(m_typ(ntyp,ntyp))
!  call mconf(m_typ(:,:))

  k_soap_out(:,:)=0d0
  k_soap_tmp(:,:)=0d0
  local_k_soap_tmp(:,:)=0d0

  do i=1,ns_data
     do k=1,i

        call set_limit_for_atoms(rangml,data_im(i),i_start_at,i_final_at)
        if (i_start_at==1)  iwj2=0
        if (i_start_at > 1) iwj2=iwmax2_soap(i_start_at-1,i)

        do ji=i_start_at,i_final_at
           iwj1=iwj2+1
           iwj2=iwmax2_soap(ji,i)

           iwk2=0
           do jk=1,data_im(k)
              iwk1=iwk2+1
              iwk2=iwmax2_soap(jk,k)

              c(:,:,:)=(0d0,0d0)

              norm_fcut_ji=0d0
              do iwj=iwj1,iwj2
                 ij=indi2_soap(iwj,i)
                 if (ji==ij) cycle

                    dxp_ji(1:3) = r_soap(1:3,ij,i) - r_soap(1:3,ji,i)
                    WHERE ( (dxp_ji.GT.0.5d0).OR.(dxp_ji.LT.-0.5d0) )
                      dxp_ji(1:3) = dxp_ji(1:3) - Dble(Nint(dxp_ji(1:3)))
                    END WHERE

                    dxp_ji = MatMul(at_soap(:,:,i),dxp_ji)/A2cm
                    r_ji = dsqrt(sum(dxp_ji(1:3)**2))
                    if (r_ji.gt.r_cut) cycle
                    if (lsoap_fcut) then
                       fcut_ji=0.5d0*(cos(pi*r_ji/r_cut)+1d0)
                    else
                       fcut_ji=1d0
                    endif
                    norm_fcut_ji=norm_fcut_ji+fcut_ji

                 norm_fcut_jk=0d0
                 do iwk=iwk1,iwk2
                    ik=indi2_soap(iwk,k)
                    if (jk==ik) cycle

                    if (weighted) then
!                       ctyp_ji=1
!                       ctyp_jk=1
                       factor_ji=w2_rho(massat(ityp(ji)),massat(ityp(ij)))
                       factor_jk=w2_rho(massat(ityp(jk)),massat(ityp(ik)))
                    else
!                       ctyp_ji=m_typ(ityp(ji),ityp(ij))
!                       ctyp_jk=m_typ(ityp(jk),ityp(ik))
                       factor_ji=1.d0
                       factor_jk=1.d0
                    endif

                    dxp_jk(1:3) = r_soap(1:3,ik,k) - r_soap(1:3,jk,k)
                    WHERE ( (dxp_jk.GT.0.5d0).OR.(dxp_jk.LT.-0.5d0) )
                      dxp_jk(1:3) = dxp_jk(1:3) - Dble(Nint(dxp_jk(1:3)))
                    END WHERE

                    dxp_jk = MatMul(at_soap(:,:,k),dxp_jk)/A2cm
                    r_jk = dsqrt(sum(dxp_jk(1:3)**2))

                    if (r_jk.gt.r_cut) cycle
                    if (lsoap_fcut) then
                       fcut_jk=0.5d0*(cos(pi*r_jk/r_cut)+1d0)
                    else
                       fcut_jk=1d0
                    endif
                    norm_fcut_jk=norm_fcut_jk+fcut_jk

                    arg_bessel=alpha_soap*r_ji*r_jk
                    if (arg_bessel==0) then
                       ms_bessel(0)=1d0
                       ms_bessel(1:l_max)=0d0
                    else
                       do l=0,l_max
                          if (l==0) then
                             ms_bessel(0)=sinh(arg_bessel)/arg_bessel
                          elseif (l==1) then
                             ms_bessel(1)=cosh(arg_bessel)/arg_bessel - sinh(arg_bessel)/(arg_bessel)**2
                          else
                             ms_bessel(l)=ms_bessel(l-2)-(2*l-1)*ms_bessel(l-1)/arg_bessel
                          endif
                       enddo
                    endif

                    do l=0,l_max
                       do m=-l,l
                          do mp=-l,l
                             c(mp,m,l) = c(mp,m,l) + 4d0 * pi * dexp(-alpha_soap*0.5*(r_ji**2+r_jk**2)) * fcut_ji * fcut_jk * factor_ji * factor_jk *&
                                               ms_bessel(l) * spherical_harm(l,m,dxp_ji(1:3)) * conjg(spherical_harm(l,mp,dxp_jk(1:3)))
                          enddo
                       enddo
                    enddo

                 enddo !ik
              enddo !ij

              if (((norm_fcut_ji==0d0).or.(norm_fcut_jk==0d0)).and.(rangml==0)) then
                 write(*,*)"norm_fcut_ji or norm_fcut_jk = 0"
                 write(*,*)"Take larger r_cut"
                 stop
              endif
              c=c/(norm_fcut_ji*norm_fcut_jk)

              if (n_soap==2) then
                 do l=0,l_max
                    do m=-l,l
                       do mp=-l,l
                          local_k_soap_tmp(sum(data_im(0:k-1))+jk,sum(data_im(0:i-1))+ji) = local_k_soap_tmp(sum(data_im(0:k-1))+jk,sum(data_im(0:i-1))+ji) + &
                                                                                            conjg(c(mp,m,l)) * c(mp,m,l)
                       enddo
                    enddo
                 enddo
              elseif (n_soap==3) then
                 do l1=0,l_max
!                    do l2=0,l_max
                    l2=l1
                    do l=abs(l1-l2),min(l_max,l1+l2)
                       do m=-l,l
                          do mp=-l,l

                             sum_bi=(0d0,0d0)
                             do m1=max(-l1,m-l2),min(l1,m+l2)
                                do m1p=max(-l1,mp-l2),min(l1,mp+l2)
                                   sum_bi = sum_bi + cg_vector(l1,m1,l2,m-m1,l,m)*cg_vector(l1,m1p,l2,mp-m1p,l,mp)*c(m1p,m1,l1)* &
                                                     c(mp-m1p,m-m1,l2)
                                enddo
                             enddo

                             local_k_soap_tmp(sum(data_im(0:k-1))+jk,sum(data_im(0:i-1))+ji) = &
                             local_k_soap_tmp(sum(data_im(0:k-1))+jk,sum(data_im(0:i-1))+ji) + &
                             conjg(c(mp,m,l)) * sum_bi

                          enddo !mp
                       enddo !m
                    enddo !l
!                 enddo !l2
                 enddo !l1
              endif

           enddo !jk
        enddo !ji
#if (PARAML)
!       call MPI_ALLREDUCE(local_k_soap_tmp,k_soap_tmp,sum(data_im)**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,codeml)
       call MPI_ALLREDUCE(local_k_soap_tmp(sum(data_im(0:k-1))+1:sum(data_im(0:k)),sum(data_im(0:i-1))+1:sum(data_im(0:i))), &
                k_soap_tmp(sum(data_im(0:k-1))+1:sum(data_im(0:k)),sum(data_im(0:i-1))+1:sum(data_im(0:i))),data_im(i)*data_im(k),&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,codeml)
! call MPI_BARRIER(MPI_COMM_WORLD,codeml)
#else
       k_soap_tmp(sum(data_im(0:k-1))+1:sum(data_im(0:k)),sum(data_im(0:i-1))+1:sum(data_im(0:i))) = &
           local_k_soap_tmp(sum(data_im(0:k-1))+1:sum(data_im(0:k)),sum(data_im(0:i-1))+1:sum(data_im(0:i)))
#endif

     enddo !k
  enddo !i

  do i=1,ns_data
     do k=1,i
        do ji=1,data_im(i)
           do jk=1,data_im(k)
              k_soap_out(sum(data_im(0:k-1))+jk,sum(data_im(0:i-1))+ji)=( k_soap_tmp(sum(data_im(0:k-1))+jk,sum(data_im(0:i-1))+ji) / &
                         dsqrt( k_soap_tmp(sum(data_im(0:k-1))+jk,sum(data_im(0:k-1))+jk) * k_soap_tmp(sum(data_im(0:i-1))+ji,sum(data_im(0:i-1))+ji) ) )**kappa_acd
!              if ((k_soap_out(sum(data_im(0:k-1))+jk,sum(data_im(0:i-1))+ji).gt.1d0).and.(rangml==0)) then
!                 write(*,*) i,k,ji,jk,k_soap_out(sum(data_im(0:k-1))+jk,sum(data_im(0:i-1))+ji)
!                 write(*,*)k_soap_tmp(sum(data_im(0:k-1))+jk,sum(data_im(0:i-1))+ji), k_soap_tmp(sum(data_im(0:k-1))+jk,sum(data_im(0:k-1))+jk),k_soap_tmp(sum(data_im(0:i-1))+ji,sum(data_im(0:i-1))+ji)
!              endif
           enddo
        enddo
     enddo
  enddo


 deallocate(r_soap,iwmax2_soap,indi2_soap,at_soap)

return
end subroutine compute_kernel_soap



















subroutine  wesley_fill_desc(n_count)
 USE T_kind_param_m, ONLY:  double
 use gen_com_m
 use var_pot
 use tab_imm_m
 use ml_in_ndm_module
 use temporary_data_cov, ONLY:  yfunc, yfunc_nd, xdesc, dim_data, dim_valid, dim_train, dim_xdesc, dim_data_train, dim_data_test, &
                                i_final_cov, i_start_cov, i_local_cov, &
                                error_valid, error_test, &
                                yfunc_valid, yfunc_train, y_test, &
                                xdesc_train, xdesc_valid,xdesc_test,y_extra,eigen_values, &
                                dim_yfunc, yfunc_average
use extrapolation
use k_cross_validation
use set_limits
use def_kernels, ONLY : length_kse
use opt_marginal_likelihood

#if(PARAML)
!for curie      use mkl_service
      use mpi
      use mod_mpi_ml
#endif
implicit none
interface compute_descriptors
   subroutine compute_descriptors(xdesc_out,icount)
    implicit none
    double precision,dimension(:,:),allocatable :: xdesc_out
    integer, optional ::  icount
   end subroutine compute_descriptors
end interface compute_descriptors
integer, intent(in) :: n_count
integer :: i, i2_count, i_count, nb
double precision :: y_target
double precision,dimension(:,:),allocatable :: xdesc_i

i_count=0
    select case (target_type)
       case (target_energy)
          dim_data=n_count!ns_data
          if (allocated(yfunc)) deallocate(yfunc); allocate(yfunc(dim_data))
          if (allocated(xdesc)) deallocate(xdesc)
          i2_count=0
          do i=i_begin+1,i_begin+ns_data
             i2_count=i2_count+1
             if (reject(i2_count)) cycle
             i_count=i_count+1
             i_poscar=i
             call read_poscar(pref,i,y_target) ! Mise à jour de im,at,bg,xp,ntyp,ityp
            call convert_A2cm(1,.true.,.true.)
             yfunc(i_count)=y_target
             im_glob=im                        ! Mise à jour de im_glob pour divid
             call alloc_typ_ml()
             call Deallocatecel()
             natperc=-1                        ! Force le calcul de natperc dans divid
             call divid(1)
             call DynamicalAllocationCell()    ! Réallocation du pointeur last(natperc,:noxyz) pour caltabt
             call neigcel()
             call caltabt()
             call caltabi()
             if (lsoap.or.acd) then
                data_im(i_count)=im
                if (size_indi2<size(indi2)) then
                    size_indi2=size(indi2)
                endif
                if (max_ntyp<size(natm)) then
                    max_ntyp=size(natm)
                endif
             endif
             call allocate_ml()
             call compute_descriptors(xdesc_i)
             if (.not.allocated(xdesc)) allocate(xdesc(dim_xdesc,dim_data))
             xdesc(1:dim_xdesc,i_count)=xdesc_i(1:dim_xdesc,1)
             call deallocate_ml()
          enddo
       case (target_force)
          dim_data=0
          i2_count=0
          do i=i_begin+1,i_begin+ns_data
             i2_count=i2_count+1
             if (reject(i2_count)) cycle
             call read_poscar(pref,i,y_target)
             dim_data=dim_data+im
          enddo
          if (allocated(yfunc)) deallocate(yfunc)
          allocate(yfunc(dim_data))
          if (allocated(xdesc)) deallocate(xdesc)
          nb=0
          i2_count=0
          do i=i_begin+1,i_begin+ns_data
             i2_count=i2_count+1
             if (reject(i2_count)) cycle
             i_count=i_count+1
             i_poscar=i
             call read_poscar(pref,i,y_target) ! Mise à jour de im,at,bg,xp,ntyp,ityp
             call convert_A2cm(1,.true.,.true.)
             call recips(at(:,1), at(:,2), at(:,3), bg(:,1), bg(:,2), bg(:,3))
             yfunc(nb+1:nb+im)=fp(force_comp,1:im)  ! the fp forces comes from shared module
             im_glob=im                        ! Mise à jour de im_glob pour divid
             call alloc_typ_ml()
             call Deallocatecel()
             natperc=-1                        ! Force le calcul de natperc dans divid
             call divid(1)
             call DynamicalAllocationCell()    ! Réallocation du pointeur last(natperc,:noxyz) pour caltabt
             call neigcel()
             call caltabt()
             call caltabi()
             if (lsoap.or.acd) then
                data_im(i_count)=im
                if (size_indi2<size(indi2)) then
                    size_indi2=size(indi2)
                endif
                if (max_ntyp<size(natm)) then
                    max_ntyp=size(natm)
                endif
             endif
             call allocate_ml()
             call compute_descriptors(xdesc_i)
             if (.not.allocated(xdesc)) allocate(xdesc(dim_xdesc,dim_data))
             xdesc(1:dim_xdesc,nb+1:nb+im)=xdesc_i(1:dim_xdesc,1:im)
             nb=nb+im
             call deallocate_ml()
          enddo
    end select

return
end subroutine  wesley_fill_desc


subroutine  wesley_fill_soap(n_count)
 USE T_kind_param_m, ONLY:  double
 use gen_com_m
 use var_pot
 use tab_imm_m
 use ml_in_ndm_module
 use temporary_data_cov, ONLY:  yfunc, yfunc_nd, xdesc, dim_data, dim_valid, dim_train, dim_xdesc, dim_data_train, dim_data_test, &
                                i_final_cov, i_start_cov, i_local_cov, &
                                error_valid, error_test, &
                                yfunc_valid, yfunc_train, y_test, &
                                xdesc_train, xdesc_valid,xdesc_test,y_extra,eigen_values, &
                                dim_yfunc, yfunc_average
use extrapolation
use k_cross_validation
use set_limits
use def_kernels, ONLY : length_kse
use opt_marginal_likelihood

#if(PARAML)
!for curie      use mkl_service
      use mpi
      use mod_mpi_ml
#endif
implicit none

integer, intent(in) :: n_count
integer :: i_count, i2_count,i
integer :: dim_soap
double precision :: y_target

if (rangml==0) write(6,*) 'ML: enter SOAP ...'

 if (allocated(r_soap))      deallocate(r_soap);      allocate(r_soap(3,maxval(data_im),n_count))
 if (allocated(iwmax2_soap)) deallocate(iwmax2_soap); allocate(iwmax2_soap(maxval(data_im),n_count))
 if (allocated(indi2_soap))  deallocate(indi2_soap);  allocate(indi2_soap(size_indi2,n_count))
 if (allocated(at_soap))     deallocate(at_soap);     allocate(at_soap(3,3,n_count))

 i_count=0
 i2_count=0
 do i=i_begin+1,i_begin+ns_data
    i2_count=i2_count+1
    if (reject(i2_count)) cycle
    i_count=i_count+1
    call read_poscar(pref,i,y_target) ! Mise à jour de im,at,bg,xp,ntyp,ityp
    call convert_A2cm(1,.true.,.true.)
    im_glob=im                        ! Mise à jour de im_glob pour divid
    call alloc_typ_ml()
    call Deallocatecel()
    natperc=-1                        ! Force le calcul de natperc dans divid
    call divid(1)
    call DynamicalAllocationCell()    ! Réallocation du pointeur last(natperc,:noxyz) pour caltabt
    call neigcel()
    call caltabt()
    call caltabi()
    call coord_soap(i_count)
 enddo

 dim_soap=0
 do i=1,n_count!ns_data
    dim_soap=dim_soap+data_im(i)
 enddo
  dim_soap=sum(data_im)
 if (allocated(k_soap)) deallocate(k_soap)
 allocate(k_soap(dim_soap,dim_soap))
 if (rangml==0) write(6,*) 'ML: compute SOAP ...'
 call compute_kernel_soap(k_soap)
 if (rangml==0) write(6,*) 'ML: end compute SOAP ...'
! if (debug) then
!    if (rangml==0) then
!       do i=1,n_count!ns_data
!          do k=1,i
!             write(99,*)i,k
!             do j=1,data_im(i)
!                write(99,'(<data_im(i)>(G16.8))')k_soap(sum(data_im(0:k-1))+1:sum(data_im(0:k)),sum(data_im(0:i-1))+j)
!             enddo
!          enddo
!       enddo
!    endif
! endif
 deallocate(k_soap,data_im)

end subroutine  wesley_fill_soap


subroutine  wesley_fill_acd(n_count)
 USE T_kind_param_m, ONLY:  double
 use gen_com_m
 use var_pot
 use tab_imm_m
 use ml_in_ndm_module
 use temporary_data_cov, ONLY:  yfunc, yfunc_nd, xdesc, dim_data, dim_valid, dim_train, dim_xdesc, dim_data_train, dim_data_test, &
                                i_final_cov, i_start_cov, i_local_cov, &
                                error_valid, error_test, &
                                yfunc_valid, yfunc_train, y_test, &
                                xdesc_train, xdesc_valid,xdesc_test,y_extra,eigen_values, &
                                dim_yfunc, yfunc_average
use extrapolation
use k_cross_validation
use set_limits
use def_kernels, ONLY : length_kse
use opt_marginal_likelihood

#if(PARAML)
!for curie      use mkl_service
      use mpi
      use mod_mpi_ml
#endif
implicit none
interface compute_descriptors
   subroutine compute_descriptors(xdesc_out,icount)
    implicit none
    double complex,dimension(:,:),allocatable :: xdesc_out
    integer, optional ::  icount
   end subroutine compute_descriptors
end interface compute_descriptors
integer, intent(in) :: n_count
integer :: i, j, k, i_count, i2_count
double precision :: y_target
double precision,dimension(:,:),allocatable :: distance_desc

 if (rangml==0) write(6,*) 'ML: compute ACD ...'
 sparsification_by_acd=.false.
 if (allocated(distance_desc)) deallocate(distance_desc); allocate(distance_desc(n_count,n_count))
 distance_desc(:,:)=0d0
 select case (target_type)
 case (target_energy)
 do i=1,n_count!ns_data
    do j=1,i
       distance_desc(j,i)=dsqrt(sum(abs(xdesc(:,j)-xdesc(:,i))**2))
    enddo
 enddo
 case (target_force)
 do i=1,n_count!ns_data
    do j=1,i
       do k=1,dim_xdesc
          distance_desc(j,i)=dsqrt(sum(abs(xdesc(k,sum(data_im(0:j-1))+1:sum(data_im(0:j)))-xdesc(k,sum(data_im(0:i-1))+1:sum(data_im(0:i))))**2)/min(data_im(i),data_im(j)))
       enddo
    enddo
 enddo
 end select
 !if (debug) then
!    if (rangml==0) then
!       write(*,*)"distance_desc"
!       do i=1,n_count!ns_data
!          write(*,'(<ns_data>(G16.10))')distance_desc(:,i)
!       enddo
!    endif
 !endif

 if (allocated(r_acd)) deallocate(r_acd)
 allocate(r_acd(3,maxval(data_im),n_count))
  allocate(r_acd(3,maxval(data_im),ns_data))
 if (allocated(data_natm)) deallocate(data_natm)
 allocate(data_natm(max_ntyp,n_count))
  allocate(data_natm(max_ntyp,ns_data))
 i_count=0
 i2_count=0
 select case (target_type)
 case (target_energy)
    do i=i_begin+1,i_begin+ns_data
       i2_count=i2_count+1
       if (reject(i2_count)) cycle
       i_count=i_count+1
       call read_poscar(pref,i,y_target) ! Mise à jour de im,at,bg,xp,ntyp,ityp
       r_acd(:,:,i_count)=xp
       data_natm(:,i_count)=natm(:)
    enddo
 case (target_force)
    if (allocated(iwmax2_acd)) deallocate(iwmax2_acd)
    if (allocated(indi2_acd)) deallocate(indi2_acd)
    if (allocated(at_acd)) deallocate(at_acd)
    allocate(iwmax2_acd(maxval(data_im),n_count))
    allocate(indi2_acd(size_indi2,n_count))
    allocate(at_acd(3,3,n_count))
    do i=i_begin+1,i_begin+ns_data
       i2_count=i2_count+1
       if (reject(i2_count)) cycle
       i_count=i_count+1
       call read_poscar(pref,i,y_target) ! Mise à jour de im,at,bg,xp,ntyp,ityp
       data_natm(:,i_count)=natm(:)
       call convert_A2cm(1,.true.,.true.)
       im_glob=im                        ! Mise à jour de im_glob pour divid
       call alloc_typ_ml()
       call Deallocatecel()
       natperc=-1                        ! Force le calcul de natperc dans divid
       call divid(1)
       call DynamicalAllocationCell()    ! Réallocation du pointeur last(natperc,:noxyz) pour caltabt
       call neigcel()
       call caltabt()
       call caltabi()
       call coord_acd_local(i_count)
    enddo
 end select

 if (allocated(k_acd)) deallocate(k_acd)
 allocate(k_acd(n_count,n_count))
 if (allocated(distance_acd)) deallocate(distance_acd)
 allocate(distance_acd(n_count,n_count))
 select case (target_type)
 case (target_energy)
    call compute_kernel_acd(n_count,k_acd,distance_acd)
 case (target_force)
    call compute_kernel_acd_local(n_count,k_acd,distance_acd)
 end select
! if (debug) then
!    if (rangml==0) then
!       write(*,*)"kernel_acd"
!       do i=1,n_count!ns_data
!          write(*,'(<n_count>(G16.10))')k_acd(:,i)
!       enddo
!       write(*,*)"distance_acd"
!       do i=1,n_count!ns_data
!          write(*,'(<n_count>(G16.10))')distance_acd(:,i)
!       enddo
!    endif
 !endif

 if (rangml==0) then
    open(65, file="distance.dat", status="unknown")
    do i=1,n_count!ns_data
       do j=1,i
          write(65,'(2(G20.10,1x))')distance_desc(j,i),distance_acd(j,i)
       enddo
    enddo
    close(65)
 endif

 deallocate(distance_desc,k_acd,distance_acd,r_acd,data_im)
 if (rangml==0) write(6,*) 'ML: end compute ACD ...'

end subroutine  wesley_fill_acd
