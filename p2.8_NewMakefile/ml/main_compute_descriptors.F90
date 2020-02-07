
!$-------------------------------------------------------------
subroutine main_compute_descriptors
!$-------------------------------------------------------------
use ml_in_ndm_module, only: iconf_data, allocate_ml, deallocate_ml, rangml, optimize_weights, snap_fit_type, snap_class_constraints, &
                            fit_lapack_qr_constraints, tmp_val_desc_max
use temporary_data_cov, only: dim_data, dim_data_train, dim_data_test, dim_data_constraints
use derived_types, only: config_real
use snap, only: i_fit_snap, dim_ene_train_snap, dim_force_train_snap, dim_stress_train_snap, &
                i_e_train_snap, i_f_train_snap, i_s_train_snap, &
                y_e_train_snap, y_f_train_snap, y_s_train_snap, &
                y_e_train_base, y_f_train_base, y_s_train_base, &
                y_e_p_a_train_base, y_e_p_a_train_snap, i_constraints_snap,  &
                dim_force_constraints, dim_stress_constraints

use descriptors_interface
use tab_imm_m
implicit none
integer :: i,einp, finp, sinp
double precision,dimension(:,:),allocatable :: xdesc_i
real(kind(1.d0)) :: tmp_val
  call prepare_database()

  dim_data=0
  dim_data_constraints=0
  dim_data_test=0
  dim_data_train=0
  dim_ene_train_snap=0
  dim_force_train_snap=0
  dim_stress_train_snap=0
  !set-up the dimensions of train matrix Amat
  !do i =1,iconf_data
  !  if (.not.(config_real(i)%train)) cycle
  !  call read_poscar_sasha(trim(config_real(i)%filename),i)
     !config_real(i)%has_force=.false.
    !call train_get_fit_dimensions(i)
  !end do

!if (rangml == 0) then
!   write(6,'("ML: the number of datapoints in train, n_train, n_E, n_F, n_S:   ", 4i8)') dim_data_train, dim_ene_train_snap, dim_force_train_snap, dim_stress_train_snap
   !write(6,'("ML: the number of datapoints in test,  n_test                :   ",  i8)') dim_data-dim_data_train
!end if

  !debug write (*,*) dim_data, dim_data_test,  dim_data_train
  i_fit_snap=0
  i_constraints_snap=0
  !computing the descriptors for the training part ...
  tmp_val=-1.d0
  do i =1,iconf_data
     if (.not.(config_real(i)%train)) cycle
     call read_poscar_sasha(trim(config_real(i)%filename),i)
     !config_real(i)%has_force=.false.
     call test_if_config_is_small(i)
     !config_real(i)%small=.true.
     call calc_neighbours(i)
     call allocate_ml
     call compute_descriptors(xdesc_i, i)
     call val_renormalize_descritors(i)
     if (tmp_val_desc_max >= tmp_val)  then
                 !debug if (rangml) write(6,*) tmp_val_desc_max
                 tmp_val = tmp_val_desc_max
      end if
      call train_deallocate_desc(i)
     call deallocate_ml
     call deallocate_real_config(i)
     call dealloc_all_tab_imm
     call DeallocateAll
  end do

  if (rangml==0) then
    write(6,*) 'ML: the max value of the descriptor is   :', tmp_val
  end if

return
end subroutine main_compute_descriptors
!<-------------------------------------------------------------
