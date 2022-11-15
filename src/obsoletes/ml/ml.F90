!!this code is copyrighted
module ml_main_mod
        use read_ml_file_mod
        use snap
        use build_subdata_mod
        use regularization_mod
        use cholesky_mod
        use toy_models
        use wesley_style_mod
        use sparsification_mod
        use main_compute_descriptors_mod
        implicit none 
        contains

subroutine ml
!-----------------------------------------------
USE T_kind_param_m, ONLY:  double
use gen_com_m
use var_pot
use tab_imm_m
use ml_in_ndm_module
!use derived_types, only : config_real
use temporary_data_cov, ONLY:  yfunc, yfunc_nd, xdesc, dim_valid, dim_train, dim_xdesc, &
                                i_final_cov, i_start_cov, &
                                error_valid, error_test, &
                                yfunc_valid, y_test, &
                                xdesc_train, xdesc_valid,xdesc_test,y_extra, &
                                dim_yfunc, yfunc_average, &
                                dim_data, dim_data_train, dim_data_test
use extrapolation
use k_cross_validation
use set_limits
use def_kernels, ONLY : length_kse
use opt_marginal_likelihood
use compute_descriptors_mod
!use snap, only: i_fit_snap, dim_ene_train_snap,dim_force_train_snap
#ifdef PARAML
!for curie      use mkl_service
      use mpi
      use mod_mpi_ml
      use gen_mpi
#endif
      implicit none


integer :: i, ik, i_count, n_count
real(kind=kind(1.d0)) :: r_coeff,r_mean
logical :: exist
!double precision,dimension(:,:),allocatable :: xdesc_i

call print_milady(rangml)

if (rangml==0) write(6,*)
if (rangml==0) write(6,*)
if (rangml==0) write(6,*)'************ DEBUT DE MACHINE LEARNING ****************'
if (rangml==0) write(6,*)
if (rangml==0) write(6,*)


call read_ml_file()
call prepare_factorial()
call periodic_table()
!if (debug) call print_message(0, 'ML: reading data')

if (toy_model) call set_toy_model
call init_descriptors


if (ml_type==ml_type_basis) then

   !  fit with regularization ...
   if (snap_fit_type==fit_home_regularization) then
      call regularization
   else
   !  regular fit ... one shot.
      call train_snap
      if (.not.train_only) call test_snap_routine
   end if

end if

if (ml_type==ml_type_descriptors) then
   call main_compute_descriptors
end if


if (ml_type==ml_type_krr) then

   ! we read how many files are in the repository given by the variable path in order to provide  nd_data
   if (build_subdata) then
     call buildsubdata(nd_data)
   endif
    !
   if (lsoap.or.acd) then
      size_indi2=0
      max_ntyp=0
      if (allocated(data_im)) deallocate(data_im) ; allocate(data_im(0:ns_data))
      data_im(0)=0
   endif
!<-----------sparsification>

   if (allocated(reject)) deallocate(reject); allocate(reject(ns_data))
   reject(:)=.false.

   if (sparsification) then
      !if (debug)  call print_message(0,'ML: computing sparsification ...')
      call try_sparsification(n_count)
   else
      n_count=ns_data
   endif
!<-----------end of sparsification>

    !if (debug)  call print_message(0,'ML: computing descriptors ...')
    i_count=0
    call wesley_fill_desc(n_count)
    if (lsoap)  call wesley_fill_soap(n_count)
    if (acd)  call wesley_fill_acd(n_count)
    !up to here we fill xdesc
    if (rescale) then
       do i=1,dim_data
          xdesc(:,i)= (xdesc(:,i)-minval(xdesc))*(s_max_r-s_min_r)/(maxval(xdesc)-minval(xdesc))+s_min_r + &
                      (0d0,1d0)*((xdesc(:,i)-minval(xdesc))*(s_max_i-s_min_i)/(maxval(xdesc)-minval(xdesc))+s_min_i)
       enddo
    endif

    if (rangml==0)  write(*,*) "dim_xdesc",dim_xdesc,"dim_data", dim_data,"sub_data",n_count!ns_data

    if (debug) then
      do i=1,dim_data
         if (rangml==0) write(91,'(2D20.8, 12D16.8)') yfunc(i), abs(xdesc(:,i))**2, xdesc(:,i)
      end do
    end if

    if (allocated(yfunc_average)) deallocate(yfunc_average); allocate(yfunc_average(dim_yfunc))


! this part is completely strage ........ir seems that the database is not recentred.
    allocate(yfunc_nd(dim_yfunc,dim_data))

    yfunc_nd(dim_yfunc,1:dim_data)=yfunc(1:dim_data)
    call recentrate_database(dim_data,dim_yfunc,yfunc_nd,yfunc_average)
    yfunc(:)=yfunc_nd(dim_yfunc,:)

!end of the part completely strange. That now is recentred ...




    if (debug) then
       if (rangml==0) write(6,*) 'ML: entering kcross  ...'
       if (rangml==0) write(6,*) 'ML: dim_data ',dim_data
    end if


    ! fix dim_data = dim_data_train + dim_data_test (nf * dim_data)
    call allocate_xdesc(n_frac, rangml)


    if (kcross) then

      call set_kcross(dim_data_train,rangml)

      if (search_hyp) then
         if (rangml==0) then
           inquire(file="res", exist=exist)
           if (exist) then
             open(63, file="res", status="old", position="append", action="write")
           else
             open(63, file="res", status="new", action="write")
           endif
         endif
      endif ! search_hyp

      if ((rangml==0).and.(kernel_type==1)) write(6,*)'lambda_krr',lambda_krr,'length_kse',length_kse
      if (rangml==0) then
          inquire(file="corr.dat", exist=exist)
          if (exist) then
             open(62, file="corr.dat", status="unknown", position="rewind", action="write")
          else
            open(62, file="corr.dat", status="new",position="append", action="write")
          endif
      endif

      r_mean=0d0
      do ik=1,n_kcross
        !the original data is splitted in train and validate ...
        !dim_train, dim_valid are the dimension for the ik process of the n_kcross
        call fill_xdesc_yfunc_kcross(ik,rangml)
        call set_limit_for_cov(rangml, dim_train,i_final_cov,i_start_cov)
        call build_kernel_matrix(xdesc_train,dim_xdesc,dim_train)
        !debug #ifdef PARAML
        !debug         if (allocated(eigen_values)) deallocate(eigen_values)
        !debug         allocate(eigen_values(dim_train))
        !debug         call diago_scalapack(dim_train,eigen_values)
        !debug #endif
        call solve_kernel_matrix(rangml)
        if (allocated(y_extra))  deallocate(y_extra); allocate(y_extra(dim_valid))
        call krr_extrapolation (dim_valid,xdesc_valid, y_extra,error_valid)
        call correlation_coef(yfunc_valid, y_extra, dim_valid, r_coeff)
        if(rangml==0) write(6,'(A,I2,1x,3(G16.8,1x))')'ML:  kcross ....:', ik,lambda_krr,length_kse, r_coeff
        do i=1,dim_valid
           if (rangml==0) write(62,'(I6,1x,3(G16.8,1x))') i, yfunc_valid(i), y_extra(i), error_valid(i)
        end do
        r_mean=r_mean+r_coeff
         !debug  write(*,*) rangml,im, dim_train, i_start_at, i_final_at,i_final_cov,i_start_cov
      enddo !ik and nk_cross
      r_mean=r_mean/n_kcross
      if (rangml==0) write(6,'(A,1x,G16.8)')'r_mean :',r_mean
      if (search_hyp.and.(rangml==0)) write(63,'(3(G16.8))')lambda_krr,length_kse,r_mean
      close(62)

    endif !kcross



    if (search_hyp) close(63)

    if (marginal_likelihood) then
       call min_marginal_likelihood
    end if


    if ((.not.marginal_likelihood).and.(.not.kcross)) then
      dim_valid=dim_data_train*0.2
      dim_train=dim_data_train-dim_valid

      call fill_xdesc_yfunc_oneshot (rangml)

      call set_limit_for_cov(rangml,dim_train,i_final_cov,i_start_cov)
      call set_unlimit_for_atoms(rangml,im,i_start_at, i_final_at)


      ! call set_neighbors_size_cov_ml (i_start_cov,i_final_cov,i_local_cov, dim_data)
      call build_kernel_matrix(xdesc_train,dim_xdesc, dim_train)
      call solve_kernel_matrix(rangml)

      if (allocated(y_extra))  deallocate(y_extra); allocate(y_extra(dim_valid))
      call krr_extrapolation (dim_valid,xdesc_valid, y_extra,error_valid)
      call correlation_coef(yfunc_valid, y_extra, dim_valid, r_coeff)
      if(rangml==0) write(6,'(A,3(G16.8,1x))')'ML:  oneshot ....:', lambda_krr,length_kse, r_coeff
      do i=1,dim_valid
        if (rangml==0) write(62,'(I6,1x,3(G16.8,1x))') i, yfunc_valid(i), y_extra(i), error_valid(i)
      end do


      if (krr_error) then
        allocate(error_valid(dim_data))
        call krr_extrapolation(dim_data, xdesc ,yfunc_valid,error_valid)
      end if
        if (rangml==0) write(*,*) 'Correlation coefficient for (both github) validation', r_coeff

    end if

    if (n_frac.gt.0d0) then

      if (allocated(y_test))  deallocate(y_test); allocate(y_test(dim_data_test))
      if (allocated(error_test)) deallocate(error_test); allocate(error_test(dim_data_test))
      y_test(:)=yfunc(dim_data_train+1:dim_data)
      call krr_extrapolation (dim_data_test,xdesc_test, y_test,error_test)
      call correlation_coef(yfunc(dim_data_train+1:dim_data), y_test, dim_data_test, r_coeff)
      if(rangml==0) write(6,*)'Correlation coefficient for prediction :', r_coeff
      do i=1,dim_data_test
         if (rangml==0) write(63,'(I4,3(G16.8))') i, yfunc(dim_data_train+i), y_test(i), error_test(i)
      end do

    endif

    if (rangml==0) then
      write(*,*) 'ML: Do not forget young Jedi: the database recentred was'
      write(*,*) 'ML: .... and also young Jedi: the results corrected should (add this value, you should):'
      write(*,*) 'ML: yfunc_average........:', yfunc_average(1:dim_yfunc)
    end if

end if ! here get out from ml_type_krr

if (rangml==0) then
    write(6,'("ML:   ----------------c u later alligator----------------  ")')
end if

if (rangml==0) write(6,*)
if (rangml==0) write(6,*)
if (rangml==0) write(6,*)'************  END  DE MACHINE LEARNING ****************'
if (rangml==0) write(6,*)
if (rangml==0) write(6,*)

#ifdef PARAML
    call MPI_BARRIER(MPI_COMM_WORLD, codeml)
    call MPI_FINALIZE(codeml)
#endif

!debug stop
end subroutine ml
end module
