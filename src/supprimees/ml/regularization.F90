module regularization_mod
        use build_subdata_mod
        use compute_descriptors_mod
        use snap
        use neighbours_mod 
        !use snap_weights_trains_mod
        implicit none
        contains
!$-------------------------------------------------------------
subroutine regularization
!$-------------------------------------------------------------
#ifdef PARAML
use mpi
use mod_mpi_ml
#endif
use ml_in_ndm_module, only: lambda_krr,rangml, vector_lambda_krr, regularization_name
use snap, only: w_params
use temporary_data_cov, only: dim_xdesc
implicit none
integer :: i_lambda, ntest
character(len=70) :: tmp_name
call prepare_train_dimensions
call prepare_test_dimensions
!in this subroutine the desc are computed on all procs but they are pack one one proc.
call compute_descriptors_all_database

!The final fitting is only one proc:  to use DGEMM threading facilities for descriptors with big dimensions
if (rangml==0) then
  write(6,'(">>>------------All descriptors are computed -------------------------------")')
end if
!due to memory allocation we do the job only one CPU rangml==0). All the others were used only for descripors calculations ...
!if (rangml==0) then
  call prepare_train_Amatrix
  do i_lambda = 1 , size(vector_lambda_krr)
    ntest=10000 + i_lambda
    write(tmp_name,'(i5)') ntest
    regularization_name=tmp_name(2:5)

    lambda_krr = vector_lambda_krr(i_lambda)
    if (rangml==0) then
        write(6,'(">>>---------------------------------------------------------------------------")')
        write(6,'("ML:   ", a4, "     RIDGE REGRESSION              lambda_krr   ", e20.8)') regularization_name, lambda_krr
        write(6,'("------------------------------------------------------------------------------")')
    end if
if (rangml==0 ) then
    call train_snap_get_parameters
end if
    call MPI_BCAST(w_params,1+dim_xdesc,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,codeml)
    call MPI_BARRIER(MPI_COMM_WORLD, codeml)
    call get_train_errors
    if (rangml==0) write(6,'("------------------------------------------------------------------------------")')
    call get_test_errors
    if (rangml==0) write(6,'("<<<--------------------------------------------------------------------------|")')
    if (rangml==0) write(6,'("   ")')
  end do


end subroutine regularization
!<-------------------------------------------------------------


!$-------------------------------------------------------------
subroutine prepare_train_Amatrix
!$-------------------------------------------------------------
use ml_in_ndm_module, only: iconf_data, allocate_ml, deallocate_ml, rangml, optimize_weights, snap_fit_type, snap_class_constraints, &
                            fit_lapack_qr_constraints, tmp_val_desc_max, lambda_krr
use temporary_data_cov, only: dim_data, dim_data_train, dim_data_test, dim_data_constraints
use derived_types, only: config_real
!use snap, only: i_fit_snap, dim_ene_train_snap, dim_force_train_snap, dim_stress_train_snap, &
!                dim_force_constraints, dim_stress_constraints
use module_snap_quadratic, only: i_e_fit_snap, i_f_fit_snap, i_s_fit_snap

implicit none
integer :: i


  !debug write (*,*) dim_data, dim_data_test,  dim_data_train
  call  train_allocate_snap()
  i_fit_snap=0
  i_e_fit_snap=0
  i_f_fit_snap=0
  i_s_fit_snap=0
  !computing the descriptors for the training part ...
  do i =1,iconf_data
     if (.not.(config_real(i)%train)) cycle
     call train_fill_Amat_with_energy(i, pack_opt=.false.)
     call train_fill_Amat_with_force (i, pack_opt= .false.)
     call train_fill_Amat_with_stress(i, pack_opt=.false.)
  end do

return
end subroutine prepare_train_Amatrix



!$-------------------------------------------------------------
subroutine prepare_train_dimensions
!$-------------------------------------------------------------
use ml_in_ndm_module, only: iconf_data, allocate_ml, deallocate_ml, rangml, optimize_weights, snap_fit_type, snap_class_constraints, &
                            fit_lapack_qr_constraints, tmp_val_desc_max, lambda_krr
use temporary_data_cov, only: dim_data, dim_data_train, dim_data_test, dim_data_constraints
use derived_types, only: config_real
use snap, only: i_fit_snap, dim_ene_train_snap, dim_force_train_snap, dim_stress_train_snap, &
                dim_force_constraints, dim_stress_constraints

implicit none
integer :: i

  call prepare_database()

  dim_data=0
  dim_data_constraints=0
  dim_data_test=0
  dim_data_train=0
  dim_ene_train_snap=0
  dim_force_train_snap=0
  dim_stress_train_snap=0
  !set-up the dimensions of train matrix Amat
  do i =1,iconf_data
    if (.not.(config_real(i)%train)) cycle
    call read_poscar_sasha(trim(config_real(i)%filename),i)
     !config_real(i)%has_force=.false.
    call train_get_fit_dimensions(i)
  end do

if (rangml == 0) then
   write(6,'("ML: the number of datapoints in train, n_train, n_E, n_F, n_S:   ", 4i8)') dim_data_train, dim_ene_train_snap, dim_force_train_snap, dim_stress_train_snap
   !write(6,'("ML: the number of datapoints in test,  n_test                :   ",  i8)') dim_data-dim_data_train
end if


return
end subroutine prepare_train_dimensions
!<-------------------------------------------------------------



!$-------------------------------------------------------------
subroutine get_train_errors
!$-------------------------------------------------------------
use ml_in_ndm_module, only: iconf_data, allocate_ml, deallocate_ml, rangml, optimize_weights, snap_fit_type, snap_class_constraints, &
                            fit_lapack_qr_constraints, tmp_val_desc_max, lambda_krr, regularization_name
use temporary_data_cov, only: dim_data, dim_data_train, dim_data_test, dim_data_constraints
use derived_types, only: config_real
use snap, only: i_fit_snap, dim_ene_train_snap, dim_force_train_snap, dim_stress_train_snap, &
                i_e_train_snap, i_f_train_snap, i_s_train_snap, &
                y_e_train_snap, y_f_train_snap, y_s_train_snap, &
                y_e_train_base, y_f_train_base, y_s_train_base, &
                y_e_p_a_train_base, y_e_p_a_train_snap, i_constraints_snap,  &
                dim_force_constraints, dim_stress_constraints
implicit none
integer :: i,einp, finp, sinp
real(kind(1.d0)) :: tmp_val
  !testing the training ...
  einp=51
  finp=52
  sinp=53

  if (dim_ene_train_snap > 0) then
     if (allocated(y_e_train_snap)) deallocate(y_e_train_snap) ; allocate(y_e_train_snap(dim_ene_train_snap))
     if (allocated(y_e_train_base)) deallocate(y_e_train_base) ; allocate(y_e_train_base(dim_ene_train_snap))
     if (allocated(y_e_p_a_train_snap)) deallocate(y_e_p_a_train_snap) ; allocate(y_e_p_a_train_snap(dim_ene_train_snap))
     if (allocated(y_e_p_a_train_base)) deallocate(y_e_p_a_train_base) ; allocate(y_e_p_a_train_base(dim_ene_train_snap))
     open (file="train_energy"//regularization_name//".out", unit=einp, action='write')
  end if


  if (dim_force_train_snap > 0) then
     if (allocated(y_f_train_snap)) deallocate(y_f_train_snap) ; allocate(y_f_train_snap(dim_force_train_snap))
     if (allocated(y_f_train_base)) deallocate(y_f_train_base) ; allocate(y_f_train_base(dim_force_train_snap))
     open (file="train_force"//regularization_name//".out", unit=finp, action='write')
  end if


  if (dim_stress_train_snap > 0) then
     if (allocated(y_s_train_snap)) deallocate(y_s_train_snap) ; allocate(y_s_train_snap(dim_stress_train_snap))
     if (allocated(y_s_train_base)) deallocate(y_s_train_base) ; allocate(y_s_train_base(dim_stress_train_snap))
     open (file="train_stress"//regularization_name//".out", unit=sinp, action='write')
  end if

  i_e_train_snap=0
  i_f_train_snap=0
  i_s_train_snap=0
  do i=1,dim_data_train
      call train_snap_compute_energy_force_stress(i,einp, finp, sinp)
  end do

  if (optimize_weights) then
     if (rangml==0) write(6,'("ML: the errors before optimization ...")')
  end if
  if (rangml==0) then
  write(6,'("ML:        ERROR TYPE                             R^2             D             RMSE          MAE    ")')
  end if
 call train_error_snap()


  !(r)evolutianary part ...
  if (optimize_weights) then
     call snap_optimize_weights
  end if

  if (optimize_weights) then
     if (rangml==0) write(6,'("ML: the errors after  optimization ...")')
     i_e_train_snap=0
     i_f_train_snap=0
     i_s_train_snap=0
     do i=1,dim_data_train
       call train_snap_compute_energy_force_stress(i,einp, finp, sinp)
     end do
     call train_error_snap()
  end if
return

end subroutine get_train_errors
!<-------------------------------------------------------------


!$-------------------------------------------------------------
subroutine compute_descriptors_all_database
!$-------------------------------------------------------------
use ml_in_ndm_module, only: rangml, iconf_data, allocate_ml, deallocate_ml
use compute_descriptors_mod
implicit none
integer :: i
double precision,dimension(:,:),allocatable :: xdesc_i
real(kind(0.d0)) :: tmp_val, tmp_val_desc_max

  tmp_val=-1.d0
  do i =1,iconf_data
     !debug if (rangml==0) write(6,*) i
     !if (.not.(config_real(i)%train)) cycle
     call test_if_config_is_small(i)
     call calc_neighbours(i)
     call allocate_ml
     call compute_descriptors(xdesc_i, i)

     !if (rangml==0) then !mpi_rangml ... 0milady ...
       call snap_pack_energy_descriptor(i)
       call snap_pack_force_descriptor(i)
       call snap_pack_stress_descriptor(i)
       call val_renormalize_descritors(i)
       if (tmp_val_desc_max >= tmp_val)  then
          tmp_val = tmp_val_desc_max
       end if
     !end if !mpi_rangml ... 0milady ...

     call train_deallocate_desc(i)
     call deallocate_ml
     !call deallocate_real_config(i)
  end do

  if (rangml==0) then
    write(6,*) 'ML: the max value of the descriptor is   :', tmp_val
  end if
end subroutine compute_descriptors_all_database
!<-------------------------------------------------------------

!$-------------------------------------------------------------
subroutine prepare_test_dimensions
!$-------------------------------------------------------------
! prepare dimension of the fit from configuration by directly reading the poscar configurations
use ml_in_ndm_module, only: rangml, debug, iconf_data, allocate_ml, deallocate_ml
use temporary_data_cov, only:  dim_data_test
use derived_types, only: config_real
use snap, only: dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap, &
                y_e_test_snap, y_f_test_snap, y_s_test_snap, &
                y_e_test_base, y_f_test_base, y_s_test_base, &
                y_e_p_a_test_snap, y_e_p_a_test_base, &
                ene_snap, fp_snap, stress_snap
use compute_descriptors_mod
implicit none
integer :: i, ix, einp, finp, sinp
double precision,dimension(:,:),allocatable :: xdesc_i
integer :: i_e_test_snap, i_f_test_snap, i_s_test_snap

  dim_data_test=0
  dim_ene_test_snap=0
  dim_force_test_snap=0
  dim_stress_test_snap=0
  do i=1,iconf_data
        if (config_real(i)%train) cycle
        call read_poscar_sasha(trim(config_real(i)%filename),i)
        call test_get_fit_dimensions(i)
        ! PAY ATTENTION OF THE FACT THAT IN THE PRESENT VERSION dim_data_train WILL BE WRONG
  end do

  if (rangml == 0) then
     write(6,'("ML: the number of datapoints in  test,  n_test, n_E, n_F, n_S:   ", 4i8)') dim_data_test, dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap
   !write(6,'("ML: the number of datapoints in test,  n_test                :   ",  i8)') dim_data-dim_data_train
  end if
  if (debug) then
    if (rangml==0) then
       write(6,'("ML: Dim for test dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap", 3i10 )') dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap
    end if
  end if

return
end subroutine prepare_test_dimensions
!<-------------------------------------------------------------



!$-------------------------------------------------------------
subroutine get_test_errors
!$-------------------------------------------------------------
use ml_in_ndm_module, only: rangml, debug, iconf_data, allocate_ml, deallocate_ml, regularization_name
use temporary_data_cov, only:  dim_data_test
use derived_types, only: config_real
!use snap, only: dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap, &
!                y_e_test_snap, y_f_test_snap, y_s_test_snap, &
!                y_e_test_base, y_f_test_base, y_s_test_base, &
!                y_e_p_a_test_snap, y_e_p_a_test_base, &
!                ene_snap, fp_snap, stress_snap
use compute_descriptors_mod
implicit none
integer :: i, ix, einp, finp, sinp
double precision,dimension(:,:),allocatable :: xdesc_i
integer :: i_e_test_snap, i_f_test_snap, i_s_test_snap

  !testing the training ...
  einp=51
  finp=52
  sinp=53
  if (dim_ene_test_snap > 0) then
     open (file="test_energy"//regularization_name//".out", unit=einp, action='write')
     if(allocated(y_e_test_snap)) deallocate(y_e_test_snap) ; allocate(y_e_test_snap(dim_ene_test_snap))
     if(allocated(y_e_test_base)) deallocate(y_e_test_base) ; allocate(y_e_test_base(dim_ene_test_snap))
     if(allocated(y_e_p_a_test_snap)) deallocate(y_e_p_a_test_snap) ; allocate(y_e_p_a_test_snap(dim_ene_test_snap))
     if(allocated(y_e_p_a_test_base)) deallocate(y_e_p_a_test_base) ; allocate(y_e_p_a_test_base(dim_ene_test_snap))
  end if

  if (dim_force_test_snap > 0) then
     open (file="test_force"//regularization_name//".out", unit=finp, action='write')
     if(allocated(y_f_test_snap)) deallocate(y_f_test_snap) ; allocate(y_f_test_snap(dim_force_test_snap))
     if(allocated(y_f_test_base)) deallocate(y_f_test_base) ; allocate(y_f_test_base(dim_force_test_snap))
  end if

  if (dim_stress_test_snap > 0) then
     open (file="test_stress"//regularization_name//".out", unit=sinp, action='write')
     if(allocated(y_s_test_snap)) deallocate(y_s_test_snap) ; allocate(y_s_test_snap(dim_stress_test_snap))
     if(allocated(y_s_test_base)) deallocate(y_s_test_base) ; allocate(y_s_test_base(dim_stress_test_snap))
  end if

  i_e_test_snap=0
  i_f_test_snap=0
  i_s_test_snap=0
  do i =1,iconf_data
    if (config_real(i)%train) cycle
    !call read_poscar_sasha(trim(config_real(i)%filename),i)
    call test_if_config_is_small(i)
    call calc_neighbours(i)
    call allocate_ml
    call md_allocate_snap_desc
    ! call compute_descriptors(xdesc_i, i)
     if (config_real(i)%has_energy) then
        i_e_test_snap = i_e_test_snap + 1
        call md_snap_compute_energy(i, pack_opt=.false.)
        y_e_test_snap(i_e_test_snap)=ene_snap
        y_e_test_base(i_e_test_snap)=config_real(i)%energy(2)
        y_e_p_a_test_snap(i_e_test_snap)=ene_snap/dble(config_real(i)%nat)
        y_e_p_a_test_base(i_e_test_snap)=config_real(i)%energy(2)/dble(config_real(i)%nat)
        write(einp,'(2e20.10," ", (a), "  ", (a))') ene_snap, y_e_test_base(i_e_test_snap),   config_real(i)%class, config_real(i)%filename
     end if
     if (config_real(i)%has_force) then
        call md_snap_compute_force(i, pack_opt=.false.)
        do ix=1,3
          y_f_test_snap(i_f_test_snap+(ix-1)*config_real(i)%nat+1:i_f_test_snap+ix*config_real(i)%nat) =fp_snap(ix,1:config_real(i)%nat)
          y_f_test_base(i_f_test_snap+(ix-1)*config_real(i)%nat+1:i_f_test_snap+ix*config_real(i)%nat) =config_real(i)%force(ix,1:config_real(i)%nat)
        end do
        do ix=1,3*config_real(i)%nat
           write(finp,'(2e20.10," ", (a), "  ", (a))') y_f_test_snap(i_f_test_snap+ix), y_f_test_base(i_f_test_snap+ix), config_real(i)%class,   config_real(i)%filename
        end do
        i_f_test_snap = i_f_test_snap + 3*config_real(i)%nat
     end if

     if (config_real(i)%has_stress) then
        call md_snap_compute_stress(i, pack_opt=.false.)
        do ix=1,3
          y_s_test_snap(i_s_test_snap+1:i_s_test_snap+6) =stress_snap(1:6)
          y_s_test_base(i_s_test_snap+1:i_s_test_snap+6) =config_real(i)%stress(1:6)
        end do
        do ix=1,6
           write(sinp,'(2e20.10," ", (a),"  ", (a))') y_s_test_snap(i_s_test_snap+ix), y_s_test_base(i_s_test_snap+ix),  config_real(i)%class, config_real(i)%filename
        end do
        i_s_test_snap = i_s_test_snap + 6
     end if

     !call train_deallocate_desc(i)
     call deallocate_ml
  end do

  close(einp)
  close(finp)
  close(sinp)
  !compute the test
  call test_error_snap

return
end subroutine get_test_errors
!<-------------------------------------------------------------
end module
