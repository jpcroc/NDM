module snap_interface
   interface

      subroutine train_fill_Amat_with_energy (iconf, pack_opt)
         integer, intent(in) :: iconf
         logical, intent(in), optional :: pack_opt
      end subroutine


      subroutine train_fill_Amat_with_force (iconf, pack_opt)
         integer, intent(in) :: iconf
         logical, intent(in), optional :: pack_opt
      end subroutine


      subroutine train_fill_Amat_with_stress (iconf, pack_opt)
         integer, intent(in) :: iconf
         logical, intent(in), optional :: pack_opt
      end subroutine

      subroutine md_snap_compute_energy(iconf, pack_opt)
         integer, intent(in) :: iconf
         logical, intent(in), optional :: pack_opt
      end subroutine

      subroutine md_snap_compute_force(iconf, pack_opt)
         integer, intent(in) :: iconf
         logical, intent(in), optional :: pack_opt
      end subroutine

      subroutine md_snap_compute_stress(iconf, pack_opt)
         integer, intent(in) :: iconf
         logical, intent(in), optional :: pack_opt
      end subroutine


  end interface

end module snap_interface


!$-------------------------------------------------------------
module module_snap_quadratic
!$-------------------------------------------------------------
! Store the tools for snap fitting
! Emat matrix for descriptors matrix to regress E ~ Descriptors in quadratic approximation .
! Fmat matrix for descriptors matrix to regress F ~ Descriptors in quadratic approximation .
! Smat matrix for descriptors matrix to regress S ~ Descriptors in quadratic approximation .
! In the case of snap:
!       Emat    (1 + dim_xdesc , dim_ene_train_snap)
!       sub_yfunc_E_train   (dim_ene_train_snap,              1)
!       Fmat    ((1 + dim_xdesc)^2 , dim_force_train_snap)
!       sub_yfunc_F_train    (dim_force_train_snap,              1)
!       w_params(1 + dim_xdesc ,   1 + dim_xdesc)
!$-------------------------------------------------------------
!training Emat and ymat_E
real(kind=kind(1.d0)), dimension(:,:), allocatable :: Emat
real(kind=kind(1.d0)), dimension(:,:), allocatable :: ymat_E
real(kind=kind(1.d0)), dimension(:)  , allocatable :: sub_yfunc_E_train
integer,               dimension(:)  , allocatable :: map_Emat_index_to_fit

real(kind=kind(1.d0)), dimension(:,:), allocatable :: Fmat
real(kind=kind(1.d0)), dimension(:,:), allocatable :: ymat_F
real(kind=kind(1.d0)), dimension(:),   allocatable :: sub_yfunc_F_train
integer,               dimension(:)  , allocatable :: map_Fmat_index_to_fit


real(kind=kind(1.d0)), dimension(:,:), allocatable :: Smat
real(kind=kind(1.d0)), dimension(:,:), allocatable :: ymat_S
real(kind=kind(1.d0)), dimension(:),   allocatable :: sub_yfunc_S_train
integer,               dimension(:)  , allocatable :: map_Smat_index_to_fit


integer :: i_e_fit_snap, i_f_fit_snap, i_s_fit_snap
!parameters
!real(kind=kind(1.d0)), dimension(:,:), allocatable :: w_params

end module module_snap_quadratic
!<-------------------------------------------------------------


!$-------------------------------------------------------------
module snap
        use build_subdata_mod
        use notperiod_mod
        use math
        
!$-------------------------------------------------------------
! Store the tools for snap fitting
! A matrix from Amat^T x  w_params = ymat.
! In the case of snap:
!       Amat    (1 + dim_xdesc , dim_data_train)
!       ymat    (dim_data_train,              1)
!       w_params(1 + dim_xdesc ,              1)
!$-------------------------------------------------------------
!training Amat and ymat
real(kind=kind(1.d0)), dimension(:,:), allocatable :: Amat
real(kind=kind(1.d0)), dimension(:,:), allocatable :: ymat
!constraints Bmat and zmat
real(kind=kind(1.d0)), dimension(:,:), allocatable :: Bmat
real(kind=kind(1.d0)), dimension(:,:), allocatable :: zmat

real(kind=kind(1.d0)), dimension(:), allocatable   :: weights_snap
!parameters
real(kind=kind(1.d0)), dimension(:,:), allocatable :: w_params

real(kind(0.d0)) :: ene_snap
real(kind(0.d0)), dimension(:,:), allocatable  :: fp_snap
real(kind(0.d0)), dimension(6)   :: stress_snap

integer :: i_fit_snap, i_constraints_snap

integer :: dim_ene_train_snap, dim_force_train_snap, dim_stress_train_snap
integer :: dim_ene_constraints, dim_force_constraints, dim_stress_constraints

integer :: i_e_train_snap, i_f_train_snap, i_s_train_snap
real(kind(0.d0)), dimension(:), allocatable:: y_e_train_snap, y_f_train_snap, y_s_train_snap
real(kind(0.d0)), dimension(:), allocatable:: y_e_train_base, y_f_train_base, y_s_train_base
real(kind(0.d0)), dimension(:), allocatable:: y_e_p_a_train_base, y_e_p_a_train_snap

!testing
integer :: dim_ene_test_snap , dim_force_test_snap , dim_stress_test_snap
real(kind(0.d0)), dimension(:), allocatable:: y_e_test_snap, y_f_test_snap, y_s_test_snap
real(kind(0.d0)), dimension(:), allocatable:: y_e_test_base, y_f_test_base, y_s_test_base
real(kind(0.d0)), dimension(:), allocatable:: y_e_p_a_test_base, y_e_p_a_test_snap


real(kind(0.d0)) :: train_rmse_energy, train_rmse_force, train_rmse_stress, &
                    train_mae_energy, train_mae_force, train_mae_stress, &
                    test_rmse_energy, test_rmse_force, test_rmse_stress, &
                    test_mae_energy, test_mae_force, test_mae_stress

integer :: dim_weights_function, dim_weights_function_full
! reduced dimension of the objective function arguments ...
real(kind(0.d0)), dimension(:), allocatable:: tmp_weights, lower_weights, upper_weights
! the original weigths vector of dimension dim_data_train
real(kind=kind(1.d0)), dimension(:), allocatable   :: tmp_weights_snap
real(kind(0.d0)) :: Jfunc

type fit_snap_type
    integer :: iconf
    logical :: force=.false.
    logical :: energy=.false.
    logical :: stress=.false.
    real(kind(0.d0)) :: weight
end type fit_snap_type

type weights_type
    integer :: db_line
    logical :: has_force=.false.
    logical :: has_energy=.false.
    logical :: has_stress=.false.
end type weights_type


type db_type
    integer :: i_e
    integer :: i_f
    integer :: i_s
end type db_type


type (fit_snap_type), dimension(:), allocatable :: fit_snap, test_snap
type (weights_type),  dimension(:), allocatable :: map_weights_in_db
type (db_type),  dimension(:), allocatable :: map_db_in_weights

!<-------------------------------------------------------------








!/------------------------------------------------------------\
!                                                             !
!                           TRAIN PART                        !
!                                                             !
!\------------------------------------------------------------/

contains

!$-------------------------------------------------------------
subroutine train_snap
!$-------------------------------------------------------------
use ml_in_ndm_module, only: iconf_data, allocate_ml, deallocate_ml, rangml, optimize_weights, snap_fit_type, snap_class_constraints, &
                            fit_lapack_qr_constraints, tmp_val_desc_max, lambda_krr, &
                            snap_order, snap_linear, snap_quadratic
use temporary_data_cov, only: dim_data, dim_data_train, dim_data_test, dim_data_constraints
use derived_types, only: config_real
!use snap, only: i_fit_snap, dim_ene_train_snap, dim_force_train_snap, dim_stress_train_snap, &
!                i_e_train_snap, i_f_train_snap, i_s_train_snap, &
!                y_e_train_snap, y_f_train_snap, y_s_train_snap, &
!                y_e_train_base, y_f_train_base, y_s_train_base, &
!                y_e_p_a_train_base, y_e_p_a_train_snap, i_constraints_snap,  &
!                dim_force_constraints, dim_stress_constraints
use module_snap_quadratic, only: i_e_fit_snap, i_f_fit_snap, i_s_fit_snap

use compute_descriptors_mod
!use snap_interface
implicit none
integer :: i,einp, finp, sinp
double precision,dimension(:,:),allocatable :: xdesc_i
real(kind(1.d0)) :: tmp_val
real(kind(0.d0)), dimension(11) :: vector_lambda
integer :: i_lambda

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

  !set-up the constraints, or the dimension of Bmat
  if (snap_fit_type == fit_lapack_qr_constraints) then
    do i =1,iconf_data
      if (.not.(config_real(i)%train) .or. config_real(i)%class /= snap_class_constraints) cycle
      call read_poscar_sasha(trim(config_real(i)%filename),i)
      !config_real(i)%has_force=.false.
      !debug write(*,*) 'iii', config_real(i)%train, config_real(i)%class, snap_class_constraints
      call train_get_constraints_dimensions(i)
    end do
  end if
  !debug write (*,*) dim_data, dim_data_test,  dim_data_train
  call  train_allocate_snap()
  i_fit_snap=0
  i_e_fit_snap=0
  i_f_fit_snap=0
  i_s_fit_snap=0
  i_constraints_snap=0
  !computing the descriptors for the training part ...
  tmp_val=-1.d0
  do i =1,iconf_data
     if (.not.(config_real(i)%train)) cycle
     !config_real(i)%has_force=.false.
     call test_if_config_is_small(i)
     !config_real(i)%small=.true.
     call calc_neighbours(i)
     call allocate_ml
     call compute_descriptors(xdesc_i, i)
     !mpi_rangml 0
     !if (rangml==0) then
     select case(snap_order)
        case(snap_linear)
          call train_fill_Amat_with_energy(i)
          call train_fill_Amat_with_force(i)
          call train_fill_Amat_with_stress(i)
          if ( (snap_fit_type == fit_lapack_qr_constraints) .and. config_real(i)%class==snap_class_constraints) then
            call train_fill_Bmat_constraints(i)
          end if
        case (snap_quadratic)
          !call train_fill_Emat_with_energy(i)
          !call train_fill_Fmat_with_force(i)
          !call train_fill_Smat_with_stress(i)
     end select
     call val_renormalize_descritors(i)
     if (tmp_val_desc_max >= tmp_val)  then
                 !debug if (rangml) write(6,*) tmp_val_desc_max
                 tmp_val = tmp_val_desc_max
     end if
     !end if
     !mpi_rangml end if  !mpi_rangml
     call train_deallocate_desc(i)

     call deallocate_ml
  end do

!#ifdef PARAML
!    call MPI_BARRIER(MPI_COMM_WORLD,codeml)
!#endif

  if (rangml==0) then
    write(6,*) 'ML: the max value of the descriptor is   :', tmp_val
  end if

  if (snap_fit_type == fit_lapack_qr_constraints ) then
     if (i_constraints_snap /= dim_data_constraints) then
        write(6,*) 'FATAL: Serious problems in imposong the constraints', i_constraints_snap, dim_data_constraints
        stop 'The dimension of constraints are ill defined'
     end if
   if ((dim_force_constraints > 0 ) .or. (dim_stress_constraints > 0 )) then
       write(6,*) 'FATAL: the constraints of forces and stress is not, YET, implemented. If you want to continue ... put the corresponding flags from TRUE to FLASE'
       stop ' not force or stress constraints'
   end if
  end if

  if (rangml==0) then !rangml=>0 running only on the proc 0 for trainning ... train 0milady ...
    call train_snap_get_parameters

    !testing the training ...
    einp=51
    finp=52
    sinp=53

    if (dim_ene_train_snap > 0) then
      if (allocated(y_e_train_snap)) deallocate(y_e_train_snap) ; allocate(y_e_train_snap(dim_ene_train_snap))
      if (allocated(y_e_train_base)) deallocate(y_e_train_base) ; allocate(y_e_train_base(dim_ene_train_snap))
      if (allocated(y_e_p_a_train_snap)) deallocate(y_e_p_a_train_snap) ; allocate(y_e_p_a_train_snap(dim_ene_train_snap))
      if (allocated(y_e_p_a_train_base)) deallocate(y_e_p_a_train_base) ; allocate(y_e_p_a_train_base(dim_ene_train_snap))
      open (file="train_energy.out", unit=einp, action='write')
    end if


    if (dim_force_train_snap > 0) then
      if (allocated(y_f_train_snap)) deallocate(y_f_train_snap) ; allocate(y_f_train_snap(dim_force_train_snap))
      if (allocated(y_f_train_base)) deallocate(y_f_train_base) ; allocate(y_f_train_base(dim_force_train_snap))
      open (file="train_force.out", unit=finp, action='write')
    end if


    if (dim_stress_train_snap > 0) then
      if (allocated(y_s_train_snap)) deallocate(y_s_train_snap) ; allocate(y_s_train_snap(dim_stress_train_snap))
      if (allocated(y_s_train_base)) deallocate(y_s_train_base) ; allocate(y_s_train_base(dim_stress_train_snap))
      open (file="train_stress.out", unit=sinp, action='write')
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

!end do

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
end if ! rangml=>0 running on the proc 0 for trainning ... train 0milady ...
return

end subroutine train_snap
!<-------------------------------------------------------------




!$-------------------------------------------------------------
subroutine  train_get_fit_dimensions(iconf)
!increment the dimension of the fit depending on the configuration and if
!the energy forces or stress are included in the fit
! Input: iconf
!        config_real(iconf)%
!                          %has_stress has_force has_energy
! Ouptput:    incremented value of
!             dim_data, dim_data_train, dim_data_test
!$-------------------------------------------------------------
use temporary_data_cov, only: dim_data, dim_data_test, dim_data_train
use derived_types, only:config_real
!use snap, only: dim_ene_train_snap, dim_force_train_snap, dim_stress_train_snap
implicit none
integer, intent(in) :: iconf
  if (config_real(iconf)%train)  then
      if (config_real(iconf)%has_force)  then
           dim_data_train=dim_data_train+3*config_real(iconf)%nat
           dim_force_train_snap = dim_force_train_snap + 3*config_real(iconf)%nat
      end if
      if (config_real(iconf)%has_energy) then
           dim_data_train=dim_data_train+1
           dim_ene_train_snap = dim_ene_train_snap + 1
      end if
      if (config_real(iconf)%has_stress) then
           dim_data_train=dim_data_train+6
           dim_stress_train_snap = dim_stress_train_snap + 6
      end if

  else
      if (config_real(iconf)%has_force)  dim_data_test=dim_data_test+3*config_real(iconf)%nat
      if (config_real(iconf)%has_energy) dim_data_test=dim_data_test+1
      if (config_real(iconf)%has_stress) dim_data_test=dim_data_test+6
  end if
  dim_data =  dim_data_test + dim_data_train
return
end subroutine  train_get_fit_dimensions
!<-------------------------------------------------------------


!$-------------------------------------------------------------
subroutine  train_get_constraints_dimensions(iconf)
!increment the dimension of the constraints matrix Bmat and Bvec (from the constraints Bmax x X)
!the energy forces or stress are included in the fit
! Input: iconf
!        config_real(iconf)%
!                          %has_stress has_force has_energy
! Ouptput:    incremented value of
!             dim_data_constraints
!$-------------------------------------------------------------
use ml_in_ndm_module, only: rangml
use temporary_data_cov, only: dim_data_constraints
use derived_types, only:config_real
!use snap, only: dim_ene_constraints, dim_force_constraints, dim_stress_constraints
implicit none
integer, intent(in) :: iconf

  if (config_real(iconf)%train)  then
      if (config_real(iconf)%has_force)  then
           dim_data_constraints=dim_data_constraints+3*config_real(iconf)%nat
           dim_force_constraints = dim_force_constraints + 3*config_real(iconf)%nat
      end if
      if (config_real(iconf)%has_energy) then
           dim_data_constraints=dim_data_constraints+1
           dim_ene_constraints = dim_ene_constraints + 1
      end if
      if (config_real(iconf)%has_stress) then
           dim_data_constraints=dim_data_constraints+6
           dim_stress_constraints = dim_stress_constraints + 6
      end if
  end if

  if (dim_data_constraints == 0) then
   if (rangml==0) then
      write(6,*) 'WARNING: despite the fact that snap_fit_type is with constraints  we do not found any valuable constraint to impose:'
      write(6,*) 'WARNING: Possibles sources:'
      write(6,*) 'WARNING:   - the classes for which the constraints are selected have no train configutations'
      write(6,*) 'WARNING:   - the configurations that have been chosen have F F F trainning flags '
   end if
  end if

  if (dim_force_constraints > 0) then
   if (rangml == 0 ) then
     write(6,*) 'WARNING: I not recommend at all to impose constraints on forces. Unless you really know what you do: goooood luck!'
   end if
  end if

return
end subroutine  train_get_constraints_dimensions
!<-------------------------------------------------------------


!$-------------------------------------------------------------
subroutine train_allocate_snap
!$-------------------------------------------------------------
! Allocate the matrix needed for snap fitting
!$-------------------------------------------------------------
! used only for training
use ml_in_ndm_module, only: rangml, snap_fit_type, fit_lapack_qr_constraints, &
                            snap_order, snap_linear, snap_quadratic
use temporary_data_cov, only: dim_data_train, &
                              xdesc_train, yfunc_train, dim_xdesc, dim_data_constraints
!use snap, only: w_params, Amat, ymat, Bmat, zmat, fit_snap, weights_snap, &
!                dim_ene_train_snap, dim_force_train_snap, dim_stress_train_snap
use module_snap_quadratic, only : Emat, Fmat, Smat, sub_yfunc_E_train, sub_yfunc_F_train, sub_yfunc_S_train, &
                           map_Emat_index_to_fit, map_Fmat_index_to_fit, map_Smat_index_to_fit
implicit none

if ((dim_xdesc == 0 ) .or. (dim_data_train ==0)) then
    if (rangml==0) write(6,*) 'The problem is not well posed dim_desc, dim_data_train', dim_xdesc, dim_data_train
    stop 'zeros in allocate_snap'
end if


    if (allocated(xdesc_train))  deallocate(xdesc_train) ; allocate(xdesc_train(dim_xdesc, dim_data_train))
    if (allocated(fit_snap))     deallocate(fit_snap)    ; allocate(fit_snap(dim_data_train))
select case(snap_order)

  case(snap_linear)
    if (allocated(yfunc_train))  deallocate(yfunc_train) ; allocate(yfunc_train(dim_data_train))


    if (allocated(Amat))         deallocate(Amat)        ; allocate(Amat(1+ dim_xdesc, dim_data_train))
    if (allocated(ymat))         deallocate(ymat)        ; allocate(ymat(dim_data_train,1))

    if (snap_fit_type == fit_lapack_qr_constraints) then
      if (allocated(Bmat))         deallocate(Bmat)          ; allocate(Bmat(1+ dim_xdesc, dim_data_constraints))
      if (allocated(zmat))         deallocate(zmat)          ; allocate(zmat(dim_data_constraints,1))
    end if

    if (allocated(w_params))        deallocate(w_params)     ; allocate(w_params(1+dim_xdesc,1))
    if (allocated(weights_snap))    deallocate(weights_snap) ; allocate(weights_snap(dim_data_train))

  case(snap_quadratic)

    if (allocated(yfunc_train))  deallocate(yfunc_train) ; allocate(yfunc_train(dim_data_train))
    if (allocated(Emat))  deallocate(Emat) ; allocate(Emat(1+dim_xdesc, dim_data_train))
    if (allocated(Fmat))  deallocate(Fmat) ; allocate(Fmat((1+dim_xdesc)**2, dim_data_train))
    !if (allocated(Smat))  deallocate(Smat) ; allocate(Smat((1+xdesc_train)**2, dim_data_train))
    if (allocated(sub_yfunc_E_train))  deallocate(sub_yfunc_E_train) ; allocate(sub_yfunc_E_train(dim_ene_train_snap))
    if (allocated(sub_yfunc_F_train))  deallocate(sub_yfunc_F_train) ; allocate(sub_yfunc_F_train(dim_force_train_snap))
    if (allocated(sub_yfunc_S_train))  deallocate(sub_yfunc_S_train) ; allocate(sub_yfunc_S_train(dim_stress_train_snap))
    if (allocated(map_Emat_index_to_fit))     deallocate(map_Emat_index_to_fit)    ; allocate(map_Emat_index_to_fit(dim_ene_train_snap))
    if (allocated(map_Fmat_index_to_fit))     deallocate(map_Fmat_index_to_fit)    ; allocate(map_Fmat_index_to_fit(dim_ene_train_snap))
    if (allocated(map_Smat_index_to_fit))     deallocate(map_Smat_index_to_fit)    ; allocate(map_Smat_index_to_fit(dim_ene_train_snap))
end select
return
end subroutine train_allocate_snap
!<-------------------------------------------------------------

!$-------------------------------------------------------------
subroutine snap_pack_energy_descriptor (iconf)
!$-------------------------------------------------------------
use temporary_data_cov, only:  dim_xdesc
use derived_types, only: config_desc, config_real
implicit none
integer, intent(in) :: iconf

if (.not.(config_real(iconf)%has_energy)) return
if (allocated(config_desc(iconf)%pack_energy))     deallocate(config_desc(iconf)%pack_energy ); allocate(config_desc(iconf)%pack_energy(dim_xdesc))
config_desc(iconf)%pack_energy(1:dim_xdesc)=SUM(config_desc(iconf)%energy(1:dim_xdesc,1:config_real(iconf)%nat), dim=2)

end subroutine snap_pack_energy_descriptor
!<-------------------------------------------------------------

!$-------------------------------------------------------------
subroutine  train_fill_Emat_with_energy (iconf, pack_opt)
!$-------------------------------------------------------------
use temporary_data_cov, only: xdesc_train, yfunc_train,  dim_xdesc
use derived_types, only: config_desc, config_real
!use snap, only: Amat, ymat, i_fit_snap, fit_snap,weights_snap
use module_snap_quadratic, only: Emat, sub_yfunc_E_train, map_Emat_index_to_fit, i_e_fit_snap

implicit none
integer, intent(in) :: iconf
!integer :: i
!real(kind(0.d0)) :: tmp(1:dim_xdesc)
logical, intent(in), optional :: pack_opt
logical :: lpack_tmp

lpack_tmp=.true.
if (present(pack_opt)) lpack_tmp = pack_opt

if (.not.(config_real(iconf)%has_energy)) return
i_fit_snap = i_fit_snap + 1
i_e_fit_snap = i_e_fit_snap + 1
if (lpack_tmp) call snap_pack_energy_descriptor(iconf)
!xdesc_train(1:dim_xdesc,i_fit_snap) = SUM(config_desc(iconf)%energy(1:dim_xdesc,1:config_real(iconf)%nat), dim=2)
xdesc_train(1:dim_xdesc,i_fit_snap) = config_desc(iconf)%pack_energy(1:dim_xdesc)


yfunc_train(i_fit_snap) = config_real(iconf)%energy(2)
sub_yfunc_E_train(i_e_fit_snap) = config_real(iconf)%energy(2)
fit_snap(i_fit_snap)%iconf = iconf
fit_snap(i_fit_snap)%energy = .true.

Emat(1,i_e_fit_snap) = dble(config_real(iconf)%nat)
Emat(2:dim_xdesc+1,i_e_fit_snap)=xdesc_train(1:dim_xdesc, i_fit_snap)
map_Emat_index_to_fit(i_e_fit_snap) = i_fit_snap
ymat(i_fit_snap,1) = yfunc_train(i_fit_snap)
weights_snap(i_fit_snap) = config_real(iconf)%w_e

return
end subroutine  train_fill_Emat_with_energy
!<-------------------------------------------------------------

!$-------------------------------------------------------------
subroutine  train_fill_Amat_with_energy (iconf, pack_opt)
!$-------------------------------------------------------------
use temporary_data_cov, only: xdesc_train, yfunc_train,  dim_xdesc
use derived_types, only: config_desc, config_real
!use snap, only: Amat, ymat, i_fit_snap, fit_snap,weights_snap

implicit none
integer, intent(in) :: iconf
!integer :: i
!real(kind(0.d0)) :: tmp(1:dim_xdesc)
logical, intent(in), optional :: pack_opt
logical :: lpack_tmp

lpack_tmp=.true.
if (present(pack_opt)) lpack_tmp = pack_opt

if (.not.(config_real(iconf)%has_energy)) return
i_fit_snap = i_fit_snap + 1
if (lpack_tmp) call snap_pack_energy_descriptor(iconf)
!xdesc_train(1:dim_xdesc,i_fit_snap) = SUM(config_desc(iconf)%energy(1:dim_xdesc,1:config_real(iconf)%nat), dim=2)
xdesc_train(1:dim_xdesc,i_fit_snap) = config_desc(iconf)%pack_energy(1:dim_xdesc)


yfunc_train(i_fit_snap) = config_real(iconf)%energy(2)
fit_snap(i_fit_snap)%iconf = iconf
fit_snap(i_fit_snap)%energy = .true.

Amat(1,i_fit_snap) = dble(config_real(iconf)%nat)
Amat(2:dim_xdesc+1,i_fit_snap)=xdesc_train(1:dim_xdesc, i_fit_snap)
ymat(i_fit_snap,1) = yfunc_train(i_fit_snap)
weights_snap(i_fit_snap) = config_real(iconf)%w_e

return
end subroutine  train_fill_Amat_with_energy
!<-------------------------------------------------------------

!$-------------------------------------------------------------
subroutine  train_fill_Bmat_constraints (iconf)
!$-------------------------------------------------------------
use temporary_data_cov, only: dim_xdesc
use derived_types, only: config_desc, config_real
!use snap, only: Bmat, zmat, i_constraints_snap
implicit none
integer, intent(in) :: iconf
!integer :: i
!real(kind(0.d0)) :: tmp(1:dim_xdesc)

  if (.not.(config_real(iconf)%has_energy)) return
  i_constraints_snap = i_constraints_snap + 1
  Bmat(1,i_constraints_snap) = dble(config_real(iconf)%nat)
  Bmat(2:dim_xdesc+1,i_constraints_snap)= SUM(config_desc(iconf)%energy(1:dim_xdesc,1:config_real(iconf)%nat), dim=2)
  zmat(i_constraints_snap,1) = config_real(iconf)%energy(2)
return
end subroutine  train_fill_Bmat_constraints
!<-------------------------------------------------------------

!$--------------------------------------------------------------
subroutine snap_pre_pack_force_descriptor (iconf)
!$--------------------------------------------------------------
!--------------------------------------------------------------------!
! Pack the descritor for forces once the derivatives are computed.   !
! Input: The number of configuration: iconf                          !
! Input: The derivatives through the config_desc(iconf)%force vector !
! Output: the vector config_desc(iconf)%pack_force(dim_of_desc,3,imm)!
!--------------------------------------------------------------------!
use ml_in_ndm_module, only : i_start_at, i_final_at, imm_neigh, rangml
use gen_com_m, only: imm
use temporary_data_cov, only:  dim_xdesc
use derived_types, only: config_desc, config_real
implicit none
integer, intent(in) :: iconf
integer :: ik, ia, ib, inn1, inn2, i_ghost
integer, dimension(:), allocatable :: itemp
integer, dimension(:,:), allocatable :: kind_itemp, proc_itemp


if (.not.(config_real(iconf)%has_force)) return
config_desc(iconf)%n_neigh_ghost(:)=0
config_desc(iconf)%kind_neigh_ghost(:,:)=0
config_desc(iconf)%kind_neigh_proc(:,:)=0

if ((i_start_at==0) .and. (i_final_at==0)) then
  return
end if
!local temp
allocate(itemp(i_start_at:i_final_at), kind_itemp(i_start_at:i_final_at, imm_neigh), proc_itemp(i_start_at:i_final_at, imm_neigh))

!new version towards EAM like
 do ik=i_start_at, i_final_at ! derivative index \partial / \partial ik,alpha
    i_ghost=0
    do inn1 = 1, config_desc(iconf)%n_neigh(ik)  ! sum over all atoms ...
              ia = config_desc(iconf)%kind_neigh(ik,inn1)
              do inn2=1,config_desc(iconf)%n_neigh(ia)
                ib = config_desc(iconf)%kind_neigh(ia,inn2)
                if (ib /= ik) cycle
                if (config_real(iconf)%small) then
                   if (SUM( (config_real(iconf)%u_ij(ik,inn1,:) + config_real(iconf)%u_ij(ia,inn2,:))**2).lt.1.d-6) then
                      i_ghost = i_ghost + 1
                      kind_itemp(ik,i_ghost) = ia
                      !temp(1:dim_xdesc, 1:3) = temp(1:dim_xdesc, 1:3)  + config_desc(iconf)%force(1:dim_xdesc,ia,inn2,1:3 )
                  end if
                else
                      i_ghost = i_ghost + 1
                      kind_itemp(ik,i_ghost) = ia
                      proc_itemp(ik,i_ghost) = config_real(iconf)%proc_atom(ia)
                      !temp(1:dim_xdesc, 1:3) = temp(1:dim_xdesc, 1:3)  + config_desc(iconf)%force(1:dim_xdesc,ia,inn2,1:3 )
                end if
              end do
     end do
    itemp(ik) = i_ghost
    if (i_ghost > imm_neigh) then
      if (rangml == 0) then
        write(6,*) "ML WARNING ghost atoms biger than max value of imm_neigh. Probably the forces are WRONG. Increse imm_neighbours or ask the master!"
        write(6,*) "master_of_puppets1986@metallica.com"
      end if
    end if
end do
config_desc(iconf)%n_neigh_ghost(i_start_at:i_final_at) =  itemp(i_start_at:i_final_at)
config_desc(iconf)%kind_neigh_proc(i_start_at:i_final_at,:) = proc_itemp(i_start_at:i_final_at,:)
config_desc(iconf)%kind_neigh_ghost(i_start_at:i_final_at,:) = kind_itemp(i_start_at:i_final_at,:)

return

end subroutine snap_pre_pack_force_descriptor
!<-------------------------------------------------------------


!$--------------------------------------------------------------
subroutine snap_pack_force_descriptor (iconf)
!$--------------------------------------------------------------
#ifdef PARAML
use mpi
use mod_mpi_ml
#endif
use ml_in_ndm_module, only: descriptor_type,  &
                            descriptor_afs, afs_dim, &
                            descriptor_bispectrum_so4, bisso4_dim, descriptor_mtp, &
                            desc_forces
use temporary_data_cov, only : dim_xdesc
use derived_types, only: config_desc
use gen_com_m, only  : imm
implicit none
integer, intent(in) :: iconf
integer :: dim_reduce

if (descriptor_type==descriptor_afs) then
    if (allocated(config_desc(iconf)%pack_force))     deallocate(config_desc(iconf)%pack_force )     ; allocate(config_desc(iconf)%pack_force (dim_xdesc,3,imm))
    config_desc(iconf)%pack_force(:,:,:)=0.d0
    call para_snap_pack_force_descriptor(iconf)

#ifdef PARAML
    dim_reduce=dim_xdesc*imm*3
    if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_force ,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
#endif


else if (descriptor_type==descriptor_bispectrum_so4) then
    if (allocated(config_desc(iconf)%pack_force))     deallocate(config_desc(iconf)%pack_force )     ; allocate(config_desc(iconf)%pack_force (dim_xdesc,3,imm))
    config_desc(iconf)%pack_force(:,:,:)=0.d0
    call para_snap_pack_force_descriptor(iconf)

#ifdef PARAML
    dim_reduce=dim_xdesc*imm*3
    if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_force ,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
#endif


else if (descriptor_type==descriptor_mtp) then
    if (allocated(config_desc(iconf)%pack_force))     deallocate(config_desc(iconf)%pack_force )     ; allocate(config_desc(iconf)%pack_force (dim_xdesc,3,imm))
    config_desc(iconf)%pack_force(:,:,:)=0.d0
    call para_snap_pack_force_descriptor(iconf)

#ifdef PARAML
    dim_reduce=dim_xdesc*imm*3
    if (desc_forces) call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_force ,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
#endif

  else
    call serial_snap_pack_force_descriptor(iconf)
end if


return
end subroutine snap_pack_force_descriptor
!<-------------------------------------------------------------



!$--------------------------------------------------------------
subroutine serial_snap_pack_force_descriptor (iconf)
!$--------------------------------------------------------------
!--------------------------------------------------------------------!
! Pack the descritor for forces once the derivatives are computed.   !
! Input: The number of configuration: iconf                          !
! Input: The derivatives through the config_desc(iconf)%force vector !
! Output: the vector config_desc(iconf)%pack_force(dim_of_desc,3,imm)!
!--------------------------------------------------------------------!
use gen_com_m, only: imm
use temporary_data_cov, only:  dim_xdesc
use derived_types, only: config_desc, config_real
implicit none
integer, intent(in) :: iconf
integer :: ik, ia, ib, inn1, inn2
real(kind(0.d0)) :: temp(dim_xdesc,3)


if (.not.(config_real(iconf)%has_force)) return
if (allocated(config_desc(iconf)%pack_force))     deallocate(config_desc(iconf)%pack_force )     ; allocate(config_desc(iconf)%pack_force (dim_xdesc,3,imm))
!new version towards EAM like
 do ik=1, config_real(iconf)%nat ! derivative index \partial / \partial ik,alpha
    temp(:,:)=0.d0
    temp(1:dim_xdesc, 1:3) = temp(1:dim_xdesc, 1:3) + config_desc(iconf)%force(1:dim_xdesc,ik,0,1:3 )
     !icount1=0
    do inn1 = 1, config_desc(iconf)%n_neigh(ik)  ! sum over all atoms ...
              ia = config_desc(iconf)%kind_neigh(ik,inn1)
              !if (ia==ik) then ! for big boxes we are never here. Only for small ...
            !          temp(1:dim_xdesc, 1:3) = temp(1:dim_xdesc, 1:3) + config_desc(iconf)%force(1:dim_xdesc,ik,0,1:3 )
                      !debug stop 'for big never ever should be here'
             ! else
              do inn2=1,config_desc(iconf)%n_neigh(ia)
                ib = config_desc(iconf)%kind_neigh(ia,inn2)
                if (ib /= ik) cycle
                if (config_real(iconf)%small) then
                   if (SUM( (config_real(iconf)%u_ij(ik,inn1,:) + config_real(iconf)%u_ij(ia,inn2,:))**2).lt.1.d-6) then
                      temp(1:dim_xdesc, 1:3) = temp(1:dim_xdesc, 1:3)  + config_desc(iconf)%force(1:dim_xdesc,ia,inn2,1:3 )
                  end if
                else
                      temp(1:dim_xdesc, 1:3) = temp(1:dim_xdesc, 1:3)  + config_desc(iconf)%force(1:dim_xdesc,ia,inn2,1:3 )
                end if
                !icount1 = icount1 + 1
              end do
              !end if


     end do
     config_desc(iconf)%pack_force(1:dim_xdesc,1:3,ik)=-temp(1:dim_xdesc,1:3)
end do
!new version of forces pre-EAM like style

!temp2(:,:) = 0.d0
!new version towards EAM like
!d     icount2=0
!d     do ia = 1, config_real(iconf)%nat  ! sum over all atoms ...
!d              if (ia == ik) then
!d              temp2(1:dim_xdesc, 1:3) = temp2(1:dim_xdesc, 1:3) + config_desc(iconf)%force(1:dim_xdesc,ik,0,1:3 )
!d              else
!d                do inn2=1,config_desc(iconf)%n_neigh(ia)
!d                  ib = config_desc(iconf)%kind_neigh(ia,inn2)
!d                  if (ib == ik) then
!d                    temp2(1:dim_xdesc, 1:3) = temp2(1:dim_xdesc, 1:3)  + config_desc(iconf)%force(1:dim_xdesc,ia,inn2,1:3 )
!d                    icount2 = icount2 + 1
!d                  end if
!d                end do
!d              end if
!d
!d     end do
!d new version of forces pre-EAM like style

!d    write(45,'(a,2i8,3i4, 2i8, 3e21.14)') trim(config_real(iconf)%filename), ik, config_desc(iconf)%n_neigh(ik), &
!d                                           config_real(iconf)%nxCell, config_real(iconf)%nyCell, config_real(iconf)%nzCell,   icount1, icount2, &
!d                                           dsqrt(SUM((temp(:,1))**2)), dsqrt(SUM((temp(:,2))**2)), dsqrt(SUM((temp(:,3))**2))

end subroutine serial_snap_pack_force_descriptor
!<-------------------------------------------------------------

!$--------------------------------------------------------------
subroutine para_snap_pack_force_descriptor (iconf)
!$--------------------------------------------------------------
!--------------------------------------------------------------------!
! Pack the descritor for forces once the derivatives are computed.   !
! Input: The number of configuration: iconf                          !
! Input: The derivatives through the config_desc(iconf)%force vector !
! Output: the vector config_desc(iconf)%pack_force(dim_of_desc,3,imm)!
!--------------------------------------------------------------------!
use gen_com_m, only: imm
use temporary_data_cov, only:  dim_xdesc
use derived_types, only: config_desc, config_real
use ml_in_ndm_module, only : i_start_at, i_final_at
implicit none
integer, intent(in) :: iconf
integer :: ik, ia, ib, inn1, inn2
real(kind(0.d0)) :: temp(dim_xdesc,3)

if ((i_start_at==0).and.(i_final_at==0)) return
if (.not.(config_real(iconf)%has_force)) return
!if (allocated(config_desc(iconf)%pack_force))     deallocate(config_desc(iconf)%pack_force )     ; allocate(config_desc(iconf)%pack_force (dim_xdesc,3,imm))
!new version towards EAM like
 do ik=i_start_at, i_final_at ! derivative index \partial / \partial ik,alpha
    temp(:,:)=0.d0
    config_desc(iconf)%pack_force(1:dim_xdesc,1:3,ik) = config_desc(iconf)%pack_force(1:dim_xdesc,1:3,ik) - config_desc(iconf)%force(1:dim_xdesc,ik,0,1:3 )
     !icount1=0
    do inn1 = 1, config_desc(iconf)%n_neigh(ik)  ! sum over all atoms ...
              ia = config_desc(iconf)%kind_neigh(ik,inn1)
              config_desc(iconf)%pack_force(1:dim_xdesc,1:3,ia) = config_desc(iconf)%pack_force(1:dim_xdesc,1:3,ia) -  config_desc(iconf)%force(1:dim_xdesc,ik,inn1,1:3 )

              !do inn2=1,config_desc(iconf)%n_neigh(ia)
              !  ib = config_desc(iconf)%kind_neigh(ia,inn2)
              !  if (ib /= ik) cycle
              !  if (config_real(iconf)%small) then
              !     if (SUM( (config_real(iconf)%u_ij(ik,inn1,:) + config_real(iconf)%u_ij(ia,inn2,:))**2).lt.1.d-6) then
              !         temp(1:dim_xdesc, 1:3) = temp(1:dim_xdesc, 1:3)  + config_desc(iconf)%force(1:dim_xdesc,ia,inn2,1:3 )
              !    end if
              !  else
              !        temp(1:dim_xdesc, 1:3) = temp(1:dim_xdesc, 1:3)  + config_desc(iconf)%force(1:dim_xdesc,ia,inn2,1:3 )
              !  end if
              !end do


     end do
     !old -v config_desc(iconf)%pack_force(1:dim_xdesc,1:3,ik)=-temp(1:dim_xdesc,1:3)
end do
!new version of forces pre-EAM like style

end subroutine para_snap_pack_force_descriptor
!<-------------------------------------------------------------


!$-------------------------------------------------------------
subroutine  train_fill_Amat_with_force (iconf, pack_opt)
!$-------------------------------------------------------------
use temporary_data_cov, only: xdesc_train, yfunc_train,  dim_xdesc
use derived_types, only: config_desc, config_real
!use snap, only: Amat, ymat, i_fit_snap, fit_snap, weights_snap
implicit none
integer, intent(in) :: iconf
integer :: ik,  i_fit_snap_ini, i_fit_snap_end
logical, intent(in), optional :: pack_opt
logical :: lpack_tmp

lpack_tmp=.true.
if (present(pack_opt)) lpack_tmp = pack_opt


!integer :: icount1, icount2
  if (.not.(config_real(iconf)%has_force)) return
  i_fit_snap_ini = i_fit_snap
  i_fit_snap_end = i_fit_snap + 3*config_real(iconf)%nat

  if (lpack_tmp) call snap_pack_force_descriptor(iconf)

  do ik=1, config_real(iconf)%nat ! descriptor index
     xdesc_train(1:dim_xdesc,i_fit_snap+1:i_fit_snap+3)  = config_desc(iconf)%pack_force(1:dim_xdesc,1:3,ik)
     yfunc_train(i_fit_snap+1:i_fit_snap+3) = config_real(iconf)%force(1:3,ik)
     i_fit_snap = i_fit_snap + 3
  end do
 !  stop 'debug foirce Amat'
  fit_snap(i_fit_snap_ini+1:i_fit_snap_end)%force = .true.
  fit_snap(i_fit_snap_ini+1:i_fit_snap_end)%iconf=iconf
  Amat(1,i_fit_snap_ini+1:i_fit_snap_end) = 0.d0
  Amat(2:dim_xdesc+1,i_fit_snap_ini+1:i_fit_snap_end)=xdesc_train(1:dim_xdesc, i_fit_snap_ini+1:i_fit_snap_end)
  ymat(i_fit_snap_ini+1:i_fit_snap_end,1) = yfunc_train(i_fit_snap_ini+1:i_fit_snap_end)
  weights_snap(i_fit_snap_ini+1:i_fit_snap_end)=config_real(iconf)%w_f

!
return
end subroutine  train_fill_Amat_with_force
!<-------------------------------------------------------------


!$--------------------------------------------------------------
subroutine snap_pack_stress_descriptor (iconf)
!$--------------------------------------------------------------
#ifdef PARAML
use mpi
use mod_mpi_ml
#endif
use ml_in_ndm_module, only: descriptor_type, &
                            descriptor_afs, afs_dim, &
                            descriptor_bispectrum_so4, bisso4_dim, descriptor_mtp, &
                            desc_forces, rangml
use temporary_data_cov, only : dim_xdesc
use derived_types, only: config_desc
use gen_com_m, only  : imm
implicit none
integer, intent(in) :: iconf
integer :: dim_reduce

if (descriptor_type==descriptor_afs) then
    if (allocated(config_desc(iconf)%pack_stress))     deallocate(config_desc(iconf)%pack_stress )     ; allocate(config_desc(iconf)%pack_stress (dim_xdesc,6))
    config_desc(iconf)%pack_stress(:,:)=0.d0
    !write(*,*) size(config_desc(iconf)%pack_stress,1),rangml
    !write(*,*) size(config_desc(iconf)%pack_stress,2),rangml
    !call MPI_BARRIER(MPI_COMM_WORLD, codeml)
    call para_snap_pack_stress_descriptor(iconf)

#ifdef PARAML
    dim_reduce=dim_xdesc*6
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_stress ,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
#endif

else if (descriptor_type==descriptor_bispectrum_so4) then
    if (allocated(config_desc(iconf)%pack_stress))     deallocate(config_desc(iconf)%pack_stress )     ; allocate(config_desc(iconf)%pack_stress (dim_xdesc,6))
    config_desc(iconf)%pack_stress(:,:)=0.d0
    !write(*,*) size(config_desc(iconf)%pack_stress,1),rangml
    !write(*,*) size(config_desc(iconf)%pack_stress,2),rangml
    !call MPI_BARRIER(MPI_COMM_WORLD, codeml)
    call para_snap_pack_stress_descriptor(iconf)

#ifdef PARAML
    dim_reduce=dim_xdesc*6
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_stress ,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
#endif


else if (descriptor_type==descriptor_mtp) then
    if (allocated(config_desc(iconf)%pack_stress))     deallocate(config_desc(iconf)%pack_stress )     ; allocate(config_desc(iconf)%pack_stress (dim_xdesc,6))
    config_desc(iconf)%pack_stress(:,:)=0.d0
    call para_snap_pack_stress_descriptor(iconf)

#ifdef PARAML
    dim_reduce=dim_xdesc*6
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_desc(iconf)%pack_stress ,dim_reduce,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,codeml)
#endif


  else
    call serial_snap_pack_stress_descriptor(iconf)
end if


return
end subroutine snap_pack_stress_descriptor
!<-------------------------------------------------------------


!$-------------------------------------------------------------
subroutine serial_snap_pack_stress_descriptor(iconf)
!$-------------------------------------------------------------
!$ This subroutine has strong interaction with neighbours subroutines
!$ The one proposed by Milady and the one proposed by NDM depending if the cell box is small or large.
!$ small=.true. ->  MiLaDy
!$ small=.false. -> NDM

use gen_com_m, only : imm, lperiod, indi2, A2cm, bg, at
use tab_imm_m, only : xp, iwmax2
!MiLaDy_interaction
use ml_in_ndm_module, only: r_cut, sign_stress, sign_stress_big_box
use derived_types, only: config_desc, config_real
use temporary_data_cov, only: dim_xdesc
implicit none
integer, intent(in) :: iconf
integer :: ia, ic
real(kind(0.d0)) :: st_temp(1:dim_xdesc,6), rr, uu(3), ds(3), force_ij(dim_xdesc,3)
integer :: iw,iw1,iw2, ia_n

real(kind(0.d0)), dimension(:,:), allocatable :: xpnp
logical :: small

if (.not.(config_real(iconf)%has_stress)) return

if (allocated(config_desc(iconf)%pack_stress))     deallocate(config_desc(iconf)%pack_stress )  ; allocate(config_desc(iconf)%pack_stress(dim_xdesc,6))
small= config_real(iconf)%small
allocate(xpnp(3,imm))
if (lperiod) then
  xpnp(:,:)=xp(:,:)
else
call notperiod(xp,xpnp)
end if
!call cryst_to_cart (imm, xpnp, bg, -1)

 iw2=0
 st_temp(:,:)=0.d0
 do ia=1, config_real(iconf)%nat ! descriptor index
     !begin small box or not 1/
     if (small) then
      iw1=1
      iw2=config_real(iconf)%n_neigh(ia)
     else
      iw1=iw2+1
      iw2=iwmax2(ia)
      !write(*,*) 'debug neigh', iw1, iw2
     end if
     ia_n=0
     do iw = iw1, iw2
        if (small) then
           ic = config_real(iconf)%kind_neigh(ia,iw)
        else
           ic=indi2(iw)
           if (ic==ia) cycle
        end if

        if (small) then
          rr = config_real(iconf)%r_ij(ia,iw)
          uu (:) = config_real(iconf)%u_ij(ia,iw,:)
          force_ij (1:dim_xdesc,1:3) = config_desc(iconf)%force(1:dim_xdesc,ia,iw,1:3)
        else
          uu(1:3) = xpnp(1:3,ia) - xpnp(1:3,ic)
          ds=MatMul(uu,bg)
          WHERE ( (ds.GT.0.5d0).OR.(ds.LT.-0.5d0) )
            ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
          END WHERE
          uu(:) = MatMul(at(:,:),ds(:))/A2cm
          rr = dsqrt( Sum( uu(1:3)**2 ) )
          if (rr >= r_cut) cycle
          ia_n=ia_n+1
          force_ij(1:dim_xdesc,1:3) = sign_stress_big_box*config_desc(iconf)%force(1:dim_xdesc,ia,ia_n,1:3)
        end if
                !rr     = config_real(iconf)%r_ij(ia,ic)
                !uu(:)  = config_real(iconf)%u_ij(ia,ic,:)
                st_temp(1:dim_xdesc, 1) = st_temp(1:dim_xdesc, 1) + force_ij(1:dim_xdesc,1)*uu(1)/rr
                st_temp(1:dim_xdesc, 2) = st_temp(1:dim_xdesc, 2) + force_ij(1:dim_xdesc,2)*uu(2)/rr
                st_temp(1:dim_xdesc, 3) = st_temp(1:dim_xdesc, 3) + force_ij(1:dim_xdesc,3)*uu(3)/rr
                st_temp(1:dim_xdesc, 4) = st_temp(1:dim_xdesc, 4) + force_ij(1:dim_xdesc,2)*uu(3)/rr
                st_temp(1:dim_xdesc, 5) = st_temp(1:dim_xdesc, 5) + force_ij(1:dim_xdesc,1)*uu(3)/rr
                st_temp(1:dim_xdesc, 6) = st_temp(1:dim_xdesc, 6) + force_ij(1:dim_xdesc,1)*uu(2)/rr
     end do
  end do

config_desc(iconf)%pack_stress(1:dim_xdesc,1:6)=sign_stress*st_temp(1:dim_xdesc,1:6)*3.d0/(config_real(iconf)%volume)

end subroutine serial_snap_pack_stress_descriptor


!$-------------------------------------------------------------
subroutine  para_snap_pack_stress_descriptor(iconf)
!$-------------------------------------------------------------
!$ This subroutine has strong interaction with neighbours subroutines
!$ The one proposed by Milady and the one proposed by NDM depending if the cell box is small or large.
!$ small=.true. ->  MiLaDy
!$ small=.false. -> NDM

use gen_com_m, only : imm, lperiod, indi2, A2cm, bg, at
use tab_imm_m, only : xp, iwmax2
!MiLaDy_interaction
use ml_in_ndm_module, only: r_cut, sign_stress, sign_stress_big_box, i_start_at, i_final_at, rangml
use derived_types, only: config_desc, config_real
use temporary_data_cov, only: dim_xdesc
implicit none
integer, intent(in) :: iconf
integer :: ia, ic
real(kind(0.d0)) :: st_temp(1:dim_xdesc,6), rr, uu(3), ds(3), force_ij(dim_xdesc,3)
integer :: iw,iw1,iw2, ia_n

real(kind(0.d0)), dimension(:,:), allocatable :: xpnp
logical :: small

if ((i_start_at==0).and.(i_final_at==0)) return
if (.not.(config_real(iconf)%has_stress)) return

small= config_real(iconf)%small
allocate(xpnp(3,imm))
if (lperiod) then
  xpnp(:,:)=xp(:,:)
else
call notperiod(xp,xpnp)
end if
!call cryst_to_cart (imm, xpnp, bg, -1)


 iw2=0
 st_temp(:,:)=0.d0
 !debug write(*,*) i_start_at, i_final_at, rangml
 if (i_start_at==1)  iw2=0
 if (i_start_at > 1) iw2=iwmax2(i_start_at-1)
 do ia=i_start_at, i_final_at ! descriptor index
     !begin small box or not 1/
     if (small) then
      iw1=1
      iw2=config_real(iconf)%n_neigh(ia)
     else
      iw1=iw2+1
      iw2=iwmax2(ia)
      !write(*,*) 'debug neigh', iw1, iw2
     end if
     ia_n=0
     do iw = iw1, iw2
        if (small) then
           ic = config_real(iconf)%kind_neigh(ia,iw)
        else
           ic=indi2(iw)
           if (ic==ia) cycle
        end if

        if (small) then
          rr = config_real(iconf)%r_ij(ia,iw)
          uu (:) = config_real(iconf)%u_ij(ia,iw,:)
          force_ij (1:dim_xdesc,1:3) = config_desc(iconf)%force(1:dim_xdesc,ia,iw,1:3)
        else
          uu(1:3) = xpnp(1:3,ia) - xpnp(1:3,ic)
          ds=MatMul(uu,bg)
          WHERE ( (ds.GT.0.5d0).OR.(ds.LT.-0.5d0) )
            ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
          END WHERE
          uu(:) = MatMul(at(:,:),ds(:))/A2cm
          rr = dsqrt( Sum( uu(1:3)**2 ) )

          if (rr >= r_cut) cycle
          !write(*,*) rangml, rr, r_cut
          ia_n=ia_n+1
          force_ij(1:dim_xdesc,1:3) = sign_stress_big_box*config_desc(iconf)%force(1:dim_xdesc,ia,ia_n,1:3)
        end if
                !rr     = config_real(iconf)%r_ij(ia,ic)
                !uu(:)  = config_real(iconf)%u_ij(ia,ic,:)
                st_temp(1:dim_xdesc, 1) = st_temp(1:dim_xdesc, 1) + force_ij(1:dim_xdesc,1)*uu(1)/rr
                st_temp(1:dim_xdesc, 2) = st_temp(1:dim_xdesc, 2) + force_ij(1:dim_xdesc,2)*uu(2)/rr
                st_temp(1:dim_xdesc, 3) = st_temp(1:dim_xdesc, 3) + force_ij(1:dim_xdesc,3)*uu(3)/rr
                st_temp(1:dim_xdesc, 4) = st_temp(1:dim_xdesc, 4) + force_ij(1:dim_xdesc,2)*uu(3)/rr
                st_temp(1:dim_xdesc, 5) = st_temp(1:dim_xdesc, 5) + force_ij(1:dim_xdesc,1)*uu(3)/rr
                st_temp(1:dim_xdesc, 6) = st_temp(1:dim_xdesc, 6) + force_ij(1:dim_xdesc,1)*uu(2)/rr
     end do
  end do

config_desc(iconf)%pack_stress(1:dim_xdesc,1:6)=sign_stress*st_temp(1:dim_xdesc,1:6)*3.d0/(config_real(iconf)%volume)

end subroutine para_snap_pack_stress_descriptor


!$-------------------------------------------------------------
subroutine  train_fill_Amat_with_stress (iconf, pack_opt)
!$-------------------------------------------------------------

use temporary_data_cov, only: xdesc_train, yfunc_train,  dim_xdesc
use derived_types, only: config_desc, config_real
!use snap, only: Amat, ymat, i_fit_snap, fit_snap, weights_snap
implicit none
integer, intent(in) :: iconf
logical, intent(in), optional :: pack_opt
integer :: i_fit_snap_ini, i_fit_snap_end
logical :: lpack_tmp

lpack_tmp=.true.
if (present(pack_opt)) lpack_tmp = pack_opt

if (.not.(config_real(iconf)%has_stress)) return
i_fit_snap_ini = i_fit_snap
i_fit_snap_end = i_fit_snap + 6
!i_fit_snap = i_fit_snap + 3*config_real(iconf)%nat
if (lpack_tmp) call snap_pack_stress_descriptor(iconf)

!xdesc_train(1:dim_xdesc,i_fit_snap_ini+1:i_fit_snap_ini+6) = st_temp(1:dim_xdesc,1:6)/(dble(config_real(iconf)%nat)*config_real(iconf)%volume)
xdesc_train(1:dim_xdesc,i_fit_snap_ini+1:i_fit_snap_ini+6) = config_desc(iconf)%pack_stress(1:dim_xdesc,1:6)
yfunc_train(i_fit_snap_ini+1:i_fit_snap_ini+6) = config_real(iconf)%stress(1:6)
Amat(1,i_fit_snap_ini+1:i_fit_snap_ini+6) = 0.d0
Amat(2:dim_xdesc+1,i_fit_snap_ini+1:i_fit_snap_ini+6)=xdesc_train(1:dim_xdesc, i_fit_snap_ini+1:i_fit_snap_ini+6)
ymat(i_fit_snap_ini+1:i_fit_snap_end,1) = yfunc_train(i_fit_snap_ini+1:i_fit_snap_end)
weights_snap(i_fit_snap_ini+1:i_fit_snap_ini+6)=config_real(iconf)%w_s

fit_snap(i_fit_snap_ini+1:i_fit_snap_ini+6)%iconf=iconf
fit_snap(i_fit_snap_ini+1:i_fit_snap_ini+6)%stress = .true.

i_fit_snap=i_fit_snap_end
!
return
end subroutine  train_fill_Amat_with_stress
!<-------------------------------------------------------------




!$-------------------------------------------------------------
subroutine train_snap_get_parameters
!$-------------------------------------------------------------
use ml_in_ndm_module, only : rangml, char_desc,  snap_fit_type, &
                             fit_home_hb, fit_lapack_qr, fit_lapack_qr_constraints, &
                             fit_lapack_ortho, fit_home_regularization, fit_lapack_svd, lambda_krr, regularization_name, &
                             factor_weight_mass, r_cut, &
                             descriptor_type, descriptor_afs, descriptor_soap, descriptor_g2, descriptor_bispectrum_so4, &
                             lsoap_diag, lsoap_norm, lsoap_lnorm, alpha_soap, n_soap, l_max, alpha_soap, r_cut_width_soap, &
                             n_cheb, n_rbf, &
                             eta_max_g2, n_g2_eta, n_g2_rs, debug, &
                             j_max, lbso4_diag, weighted, periodic_table_element, fix_no_of_elements, fix_type_to_periodic
use temporary_data_cov, only: dim_data_train, dim_xdesc
!use snap, only: w_params, Amat, ymat, Bmat, zmat, i_fit_snap, weights_snap
implicit none
real(kind=kind(1.d0)), dimension(size(Amat,1), size(Amat,1)) :: phi, phi_inv
real(kind=kind(1.d0)), dimension(size(Amat,1), size(Amat,1)) :: phi_diag
real(kind=kind(1.d0)), dimension(size(Amat,1), size(ymat,2)) :: rhs
integer :: i
real(kind=kind(1.d0)), dimension(size(Amat,1), size(Amat,2)) :: Cmat
real(kind=kind(1.d0)), dimension(size(Amat,2), size(Amat,1)) :: Cmat_transpose
real(kind=kind(1.d0)), dimension(size(Amat,2), size(Amat,1)) :: Amat_w
real(kind=kind(1.d0)), dimension(size(ymat,1), size(ymat,2)) :: ymat_copy
real(kind=kind(1.d0)), dimension(:), allocatable :: work
real(kind=kind(1.d0)), dimension(:,:), allocatable :: y_svd
integer,               dimension(:), allocatable :: iwork
integer , dimension(:), allocatable :: JPVT
integer :: lwork, liwork, info, LDA, LDB, RANK, Mline, Ncolm, ii
real(kind=kind(1.d0)) :: RCOND
integer :: n1_lsoap, n2_lsoap, n3_lsoap

if (i_fit_snap /= dim_data_train) then
    if (rangml==0) write(6,*) 'The problem is not correct sized', i_fit_snap, dim_data_train
    stop 'incorrect sized problem train_snap_get_parameters'
end if


Cmat(:,:)=Amat(:,:)
! the solution is (A^T w A)^-1 A^T W ymat or (Amat w Amat^T)^-1 Amat w ymat
! Cmat = Amat x w
do i=1, size(Cmat,2)
    !if (optimize_weights) then
    !  call dscal (size(Cmat,1), tmp_weights_snap(i), Cmat(1,i),1 )
    !else
      !call dscal (size(Cmat,1), weights_snap(i), Cmat(1,i),1 )
      Cmat(:,i)=Amat(:,i)*weights_snap(i)
    !end if
end do

select case (snap_fit_type)


   case (fit_home_hb)
     !this_take_long_time_replace_with_dgemm: phi = matmul (Amat, transpose(Cmat))
     Cmat_transpose(:,:) = transpose(Cmat)
     call dreal_matmul (Amat, size(Amat,1), size(Amat,2), Cmat_transpose, size(Cmat_transpose,1), size(Cmat_transpose,2), &
                        phi, size(phi,1), size(phi,2))


      if (lambda_krr >= 0.d0) then
        do i =1,size(Amat,1)
          phi_diag(i,i) = phi(i,i)
          phi_diag(i,i) = 1.d0
        end do
        phi = phi + lambda_krr * phi_diag
      end if

     call dreal_inverse(size(phi,1), phi, phi_inv)
     rhs = matmul (Cmat, ymat)
     w_params(:,:)=matmul(phi_inv,rhs)



   case (fit_home_regularization)
     phi_diag(:,:) = 0.d0
     !this_take_long_time_replace_with_dgemm: phi = matmul (Amat, transpose(Cmat))
     Cmat_transpose(:,:) = transpose(Cmat)
     call dreal_matmul (Amat, size(Amat,1), size(Amat,2), Cmat_transpose, size(Cmat_transpose,1), size(Cmat_transpose,2), &
                        phi, size(phi,1), size(phi,2))
     do i =1,size(Amat,1)
       phi_diag(i,i) = phi(i,i)
       phi_diag(i,i) = 1.d0
     end do

     phi = phi + lambda_krr * phi_diag
     call dreal_inverse(size(phi,1), phi, phi_inv)
     rhs = matmul (Cmat, ymat)
     w_params(:,:)=matmul(phi_inv,rhs)


   case (fit_lapack_qr)
     lwork=-1
     if (allocated(work)) deallocate(work) ; allocate(work(1))
     LDA=size(Cmat,2)
     LDB=size(ymat,1)
     ymat_copy = ymat
     call dgels ('N', size(Cmat,2), size(Cmat,1), size(ymat,2), transpose(Cmat), &
                LDA, ymat_copy, LDB, work, lwork, info)
     lwork=int(work(1))+2
     if (allocated(work)) deallocate(work) ; allocate(work(lwork))
     call dgels ('N', size(Cmat,2), size(Cmat,1), size(ymat_copy,2), transpose(Cmat), &
                LDA, ymat_copy, LDB, work, lwork, info)

     if( info > 0 ) then
         WRITE(*,*)'The diagonal element ',INFO,' of the triangular '
         WRITE(*,*)'factor of A is zero, so that A does not have full '
         WRITE(*,*)'rank; the least squares solution could not be '
         WRITE(*,*)'computed.'
         STOP
     end if
     w_params(:,1) = ymat_copy(1:size(Cmat,1),1)

   !snap_fit_type=3 ortho decomposition, rank estimation ...
   case (fit_lapack_ortho)
     if (allocated(work)) deallocate(work) ; allocate(work(1))
     Amat_w=transpose(Cmat)
     Mline=size(Cmat,2)
     Ncolm=size(Cmat,1)
     if (allocated(JPVT)) deallocate(JPVT) ; allocate(JPVT(Ncolm))
     JPVT(:)=0
     LDA=max(Mline,1)
     LDB=max(Mline, max(Ncolm,1))
     ymat_copy = ymat
     !RCOND=1.d-20
     RCOND=-1
     lwork=-1
     call dgelsy (Mline, Ncolm, size(ymat,2), Amat_w, &
                LDA, ymat_copy, LDB, JPVT, RCOND, RANK, work, lwork, info)
     lwork=int(work(1))+2
     if (allocated(work)) deallocate(work) ; allocate(work(lwork))
     call dgelsy (Mline, Ncolm, size(ymat_copy,2), Amat_w, &
            LDA, ymat_copy, LDB, JPVT, RCOND, RANK, work, lwork, info)

     if (rangml==0) then
        write(6,'("ML: dgelsd inversion info lwork RANK ....:", i4, i6,i5)') &
              info, lwork, rank
     end if

     if( info > 0 ) then
         WRITE(*,*)'The diagonal element ',INFO,' of the triangular '
         WRITE(*,*)'factor of A is zero, so that A does not have full '
         WRITE(*,*)'rank; the least squares solution could not be '
         WRITE(*,*)'computed.'
         STOP
     end if
     w_params(:,1) = ymat_copy(1:Ncolm,1)

   !snap_fit_type=4 SVD, rank estimation
   case (fit_lapack_svd)
     if (allocated(work)) deallocate(work) ; allocate(work(1))
     if (allocated(iwork)) deallocate(iwork) ; allocate(iwork(1))

     if (lambda_krr < 0 ) then
        Amat_w=transpose(Cmat)
        Mline=size(Amat_w,1)
        Ncolm=size(Amat_w,2)



        if (allocated(JPVT)) deallocate(JPVT) ; allocate(JPVT(Ncolm))
        JPVT(:)=0
        LDA=max(Mline,1)
        LDB=max(Mline, max(Ncolm,1))

        ymat_copy = ymat
        !RCOND=1.d-20

        RCOND=-1
        lwork=-1
        call dgelsd (Mline, Ncolm, size(ymat,2), Amat_w, &
                LDA, ymat_copy, LDB, JPVT, RCOND, RANK, work, lwork, iwork, info)
        lwork=int(work(1))+2
        liwork=int(iwork(1))+2

        if (allocated(work)) deallocate(work) ; allocate(work(lwork))
        if (allocated(iwork)) deallocate(iwork) ; allocate(iwork(liwork))

        call dgelsd (Mline, Ncolm, size(ymat_copy,2), Amat_w, &
            LDA, ymat_copy, LDB, JPVT, RCOND, RANK, work, lwork, iwork, info)

        w_params(:,1) = ymat_copy(1:Ncolm,1)

     else

        phi_diag(:,:) = 0.d0
        !this_take_long_time_replace_with_dgemm: phi = matmul (Amat, transpose(Cmat))
        Cmat_transpose(:,:) = transpose(Cmat)
        call dreal_matmul (Amat, size(Amat,1), size(Amat,2), Cmat_transpose, size(Cmat_transpose,1), size(Cmat_transpose,2), &
                        phi, size(phi,1), size(phi,2))
        do i =1,size(Amat,1)
          phi_diag(i,i) = phi(i,i)
          phi_diag(i,i) = 1.d0
        end do

        phi = phi + lambda_krr * phi_diag


        Mline=size(phi,1)
        Ncolm=size(phi,2)

        if (allocated(JPVT)) deallocate(JPVT) ; allocate(JPVT(Ncolm))
        JPVT(:)=0
        LDA=max(Mline,1)
        LDB=max(Mline, max(Ncolm,1))
        if (allocated(y_svd)) deallocate(y_svd) ; allocate(y_svd(size(Cmat,1),size(ymat,2)))

      !debug if (rangml==0) then
      !debug   write(6,*) size(Cmat,1), size(Cmat,2)
      !debug write(6,*) size(ymat,1), size(ymat,2)
      !debug write(6,*) size(y_svd,1), size(y_svd,2)
      !debug
      !debug end if


        call dreal_matmul(Cmat, size(Cmat,1), size(Cmat,2), ymat, size(ymat,1), size(ymat,2), &
                                 y_svd, size(y_svd,1), size(y_svd,2))

        RCOND=-1
        lwork=-1
        call dgelsd (Mline, Ncolm, size(y_svd,2), phi, &
                LDA, y_svd, LDB, JPVT, RCOND, RANK, work, lwork, iwork, info)
        lwork=int(work(1))+2
        liwork=int(iwork(1))+2

        if (allocated(work)) deallocate(work) ; allocate(work(lwork))
        if (allocated(iwork)) deallocate(iwork) ; allocate(iwork(liwork))

        call dgelsd (Mline, Ncolm, size(y_svd,2), phi, &
            LDA, y_svd, LDB, JPVT, RCOND, RANK, work, lwork, iwork, info)

        w_params(:,1) = y_svd(1:Ncolm,1)

     end if
     if (debug) then
      if (rangml==0) then
          write(6,'("ML: dgelsd inversion info lwork liwork RANK ....:", i4, 2i6,i5)') &
              info, lwork, liwork, rank
      end if 
     end if
     if( info > 0 ) then
      if (rangml==0) then
          write(6,'("ML: dgelsd inversion info lwork liwork RANK ....:", i4, 2i6,i5)') &
              info, lwork, liwork, rank
          WRITE(*,*)'The diagonal element ',INFO,' of the triangular '
          WRITE(*,*)'factor of A is zero, so that A does not have full '
          WRITE(*,*)'rank; the least squares solution could not be '
          WRITE(*,*)'computed.'
          STOP
      end if
     end if

   case (fit_lapack_qr_constraints)
     !call dgglse(m, n, p, A = Amat^T, lda, B=Bmat^T, ldb, c=ymat, d=zmat, x, work, lwork, info)
     !integer m: The number of rows of the matrix A (m ≥ 0).
     !integer n: The number of columns of the matrices A and B (n ≥ 0).
     !integer p:  The number of rows of the matrix B (0 ≤ p ≤ n ≤ m+p).
     !integer lda: max(1,m)
     !integer ldb: max(1,p)
     !lwork= max(1,m+n+p), if lwork=-1 this is calculated automatically
     lwork = max(1, size(Cmat,2)+size(Cmat,1)+size(Bmat,2)+size(Bmat,1))
     if (allocated(work)) deallocate(work) ; allocate(work(lwork))
     if (size(Bmat,2) > size(Cmat,1)) then
       if (rangml==0) then
          write(6,*) 'The constraints problem is ill posed.'
          write(6,*) 'The number of constraints should be lower than the dimesion of the descriptor.'
          write(6,*) 'xdesc dimesion ', size(Cmat,1), 'number of constraints', size(Bmat,2)
       end if
     end if
     call dgglse(size(Cmat,2), size(Cmat,1), size(Bmat,2), transpose(Cmat), max(1,size(Cmat,2)), &
               transpose(Bmat), max(1,size(Bmat,2)), ymat(:,1), zmat(:,1), w_params(:,1), work, lwork, info)
     write(*,*) 'LSE + Constraints with info', info
!

end select


!>begin older version without weights ....
!rhs = matmul (Amat, ymat)
!phi = matmul (Amat, transpose(Amat))
!call dreal_inverse(size(phi,1), phi, phi_inv)
!w_params(:,:)=matmul(phi_inv,rhs)
!>end older version without weights ....
if (rangml==0) then
   if (fix_no_of_elements /= 1) then
      write(6,*) 'WARNING NOT YET IMPLEMENTED FOR fix_no_of_elements > 1', fix_no_of_elements
   end if

   if (snap_fit_type == fit_home_regularization) then
    open(file=char_desc//'_snap1_params.pot_'//regularization_name, unit=51, action='write')
    open(file='lammps_'//char_desc//'_snap1_params.pot_'//regularization_name, unit=52, action='write')
   else
    open(file=char_desc//'_snap1_params.pot', unit=51, action='write')
    open(file='lammps_'//char_desc//'_snap1_params.pot', unit=52, action='write')
   end if

   write(52,'("# nelements elements")')
   write(52,'("# descriptor type, nweights, ndetails, weight_factor")')
   write(52,'("# mass_1,...mass_nelements")')
   write(52,'("# detail_1,...detail_ndetails (e.g. no_of_radial_channels no_of_angular_channels)")')
   write(52,'("# r_cut_ij for i=1,N j=(i,N) where N= number of elements; i.e. r_11 r_12 r_22 for N=2")')
   write(52,'("# w_ij  for i=1,N j=(i,N) ( i.e. w_11 w_12 w_22 per line for N=2)")')
   !write(CHFMT,*)'(i3, 1x, ',int(fix_no_of_elements),'(a2, " ")')'
   write(52,*) fix_no_of_elements, (('  '//periodic_table_element(fix_type_to_periodic(ii))%symbol), ii=1,fix_no_of_elements)
   select case(descriptor_type)
   case (descriptor_g2)
    if (weighted) then
       write(52,'(i4, i6, i4, f20.10)') descriptor_type, size(w_params,1), 3, factor_weight_mass
    else
       write(52,'(i4, i6, i4, f20.10)') descriptor_type, size(w_params,1), 3, 1.0
    end if
    write(52,'(f20.4)') periodic_table_element(fix_type_to_periodic(fix_no_of_elements))%mass
    write(52,'(f20.6, 2i5)') eta_max_g2, n_g2_eta, n_g2_rs
   case(descriptor_afs)
    if (weighted) then
       write(52,'(i4, i6, i4, f20.10)') descriptor_type, size(w_params,1), 2, factor_weight_mass
    else
       write(52,'(i4, i6, i4, f20.10)') descriptor_type, size(w_params,1), 2, 1.0
    end if
    write(52,'(f20.4)') periodic_table_element(fix_type_to_periodic(fix_no_of_elements))%mass
    write(52,'(2i5)') n_rbf, n_cheb
   case (descriptor_soap)
    if (weighted) then
       write(52,'(i4, i6, i4, f20.10)') descriptor_type, size(w_params,1), 7, factor_weight_mass
    else
       write(52,'(i4, i6, i4, f20.10)') descriptor_type, size(w_params,1), 7, 1.0
    end if

    write(52,'(f20.4)') periodic_table_element(fix_type_to_periodic(fix_no_of_elements))%mass
    n1_lsoap=0 ; if (lsoap_diag)  n1_lsoap=1
    n2_lsoap=0 ; if (lsoap_norm)  n2_lsoap=1
    n3_lsoap=0 ; if (lsoap_lnorm) n3_lsoap=1
    write(52,'(2i4, f12.5, f12.5, 3i4 )' ) n_soap, l_max, alpha_soap, r_cut_width_soap,n1_lsoap, n2_lsoap, n3_lsoap
   case (descriptor_bispectrum_so4)
    if (weighted) then
       write(52,'(i4, i6, i4, f20.10)') descriptor_type, size(w_params,1), 2, factor_weight_mass
    else
       write(52,'(i4, i6, i4, f20.10)') descriptor_type, size(w_params,1), 2, 1.0
    end if

    write(52,'(f20.4)') periodic_table_element(fix_type_to_periodic(fix_no_of_elements))%mass
    n1_lsoap=0 ; if (lbso4_diag)  n1_lsoap=1
    write(52,'(f5.1,i4)') j_max, n1_lsoap
   end select
   write(52,'(f20.4)') r_cut*2.d0


   if (snap_fit_type == fit_home_regularization) then
    open(file=char_desc//'_snap1_params.pot_'//regularization_name, unit=51, action='write')
   else
    open(file=char_desc//'_snap1_params.pot', unit=51, action='write')
   end if
   write(51,*) dim_xdesc + 1
   do i=1,1+ dim_xdesc
      write(51,*) w_params(i,1)
      write(52,*) w_params(i,1)
   end do
   close(51)
   close(52)
end if ! rangml==0

return
end subroutine train_snap_get_parameters
!<-------------------------------------------------------------

!$-------------------------------------------------------------
subroutine train_snap_compute_energy_force_stress(idata, einp, finp, sinp)
!$-------------------------------------------------------------
!use snap, only : w_params, Amat, ene_snap, fit_snap, &
!                 i_e_train_snap, i_f_train_snap, i_s_train_snap, &
!                 y_e_train_base, y_f_train_base, y_s_train_base, &
!                 y_e_train_snap, y_f_train_snap, y_s_train_snap, &
!                 y_e_p_a_train_snap, y_e_p_a_train_base

use temporary_data_cov, only:  yfunc_train
use derived_types, only: config_real
implicit none
integer, intent(in) :: idata, einp, finp, sinp
real(kind(0.d0)) :: force_snap,stress_snap

if (fit_snap(idata)%energy) then
     i_e_train_snap=i_e_train_snap + 1
     ene_snap = dot_product(w_params(:,1), Amat(:,idata))
     write(einp,'(2e20.10," ", (a)," ",(a))') ene_snap, yfunc_train(idata),   config_real(fit_snap(idata)%iconf)%class, config_real(fit_snap(idata)%iconf)%filename
     !debug write(73,'(2e20.10," ", (a)," ",(a),"  ", 2e20.10)') ene_snap, yfunc_train(idata),   config_real(fit_snap(idata)%iconf)%class, &
     !debug      config_real(fit_snap(idata)%iconf)%filename, config_real(fit_snap(idata)%iconf)%w_e, weights_snap(idata)

     y_e_train_snap(i_e_train_snap)=ene_snap
     y_e_train_base(i_e_train_snap)=yfunc_train(idata)
     y_e_p_a_train_snap(i_e_train_snap)=ene_snap/dble(config_real(fit_snap(idata)%iconf)%nat)
     y_e_p_a_train_base(i_e_train_snap)=yfunc_train(idata)/dble(config_real(fit_snap(idata)%iconf)%nat)
end if


if (fit_snap(idata)%force) then
     i_f_train_snap=i_f_train_snap + 1
     force_snap = dot_product(w_params(:,1), Amat(:,idata))
     write(finp,'(2e20.10," ",(a)," ",(a) )') force_snap, yfunc_train(idata),  config_real(fit_snap(idata)%iconf)%class, config_real(fit_snap(idata)%iconf)%filename
     y_f_train_snap(i_f_train_snap)=force_snap
     y_f_train_base(i_f_train_snap)=yfunc_train(idata)
end if


if (fit_snap(idata)%stress) then
     i_s_train_snap=i_s_train_snap + 1
     stress_snap = dot_product(w_params(:,1), Amat(:,idata))
     write(sinp,'(2e20.10," ",(a), " ", (a))') stress_snap, yfunc_train(idata),    config_real(fit_snap(idata)%iconf)%class, config_real(fit_snap(idata)%iconf)%filename
     y_s_train_snap(i_s_train_snap)=stress_snap
     y_s_train_base(i_s_train_snap)=yfunc_train(idata)
end if


return
end subroutine train_snap_compute_energy_force_stress
!<-------------------------------------------------------------



!$-------------------------------------------------------------
subroutine train_error_snap
!$-------------------------------------------------------------
use ml_in_ndm_module, only: rangml
!use snap, only: dim_ene_train_snap, dim_force_train_snap, dim_stress_train_snap, &
!                y_e_train_snap, y_f_train_snap, y_s_train_snap, &
!                y_e_train_base, y_f_train_base, y_s_train_base, &
!                y_e_p_a_train_base, y_e_p_a_train_snap,  &
!                train_rmse_energy, train_rmse_force, train_rmse_stress, &
!                train_mae_energy, train_mae_force, train_mae_stress
implicit none

real(kind(0.d0)) :: corr, detr, rmse, mae


  if (dim_ene_train_snap > 0) then
      call correlation_coef  (y_e_train_base, y_e_train_snap, dim_ene_train_snap, corr)
      call determination_coef(y_e_train_base, y_e_train_snap, dim_ene_train_snap, detr)
      call rmse_mae          (y_e_train_base, y_e_train_snap, dim_ene_train_snap, rmse, mae)
      if (rangml==0) then
          write(6,'("ML: stat  Total_Energy            train:", 4f15.6)')  corr, detr, rmse, mae
      end if

      call correlation_coef  (y_e_p_a_train_base, y_e_p_a_train_snap, dim_ene_train_snap, corr)
      call determination_coef(y_e_p_a_train_base, y_e_p_a_train_snap, dim_ene_train_snap, detr)
      call rmse_mae          (y_e_p_a_train_base, y_e_p_a_train_snap, dim_ene_train_snap, rmse, mae)
      if (rangml==0) then
          write(6,'("ML: stat  Total_Energy_per_atom   train:", 4f15.6)')  corr, detr, rmse, mae
      end if
      train_rmse_energy = rmse
      train_mae_energy = mae
  end if

  if (dim_force_train_snap > 0) then

      call correlation_coef  (y_f_train_base, y_f_train_snap, dim_force_train_snap, corr)
      call determination_coef(y_f_train_base, y_f_train_snap, dim_force_train_snap, detr)
      call rmse_mae          (y_f_train_base, y_f_train_snap, dim_force_train_snap, rmse, mae)
      if (rangml==0) then
          write(6,'("ML: stat  Force_each_component    train:", 4f15.6)')  corr, detr, rmse, mae
      end if
      train_rmse_force = rmse
      train_mae_force = mae
  end if

  if (dim_stress_train_snap > 0) then

      call correlation_coef  (y_s_train_base, y_s_train_snap, dim_stress_train_snap, corr)
      call determination_coef(y_s_train_base, y_s_train_snap, dim_stress_train_snap, detr)
      call rmse_mae          (y_s_train_base, y_s_train_snap, dim_stress_train_snap, rmse, mae)
      if (rangml==0) then
         write(6,'("ML: stat  Stress_virial           train:", 4f15.6)')  corr, detr, rmse, mae
      end if
      train_rmse_stress = rmse
      train_mae_stress = mae
  end if

return
end subroutine train_error_snap
!<-------------------------------------------------------------




!/------------------------------------------------------------\
!                                                             !
!                           TEST PART                         !
!                                                             !
!\------------------------------------------------------------/

!$-------------------------------------------------------------
subroutine test_snap_routine
!$-------------------------------------------------------------
use ml_in_ndm_module, only: rangml, debug, iconf_data, allocate_ml, deallocate_ml
use temporary_data_cov, only:  dim_data_test
use derived_types, only: config_real
!use snap, only: dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap, &
!                y_e_test_snap, y_f_test_snap, y_s_test_snap, &
!                y_e_test_base, y_f_test_base, y_s_test_base, &
!                y_e_p_a_test_snap, y_e_p_a_test_base, &
!                ene_snap, fp_snap, stress_snap
use compute_descriptors_mod
!use snap_interface
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
  !testing the training ...
  einp=51
  finp=52
  sinp=53
  if (dim_ene_test_snap > 0) then
     open (file="test_energy.out", unit=einp, action='write')
     if(allocated(y_e_test_snap)) deallocate(y_e_test_snap) ; allocate(y_e_test_snap(dim_ene_test_snap))
     if(allocated(y_e_test_base)) deallocate(y_e_test_base) ; allocate(y_e_test_base(dim_ene_test_snap))
     if(allocated(y_e_p_a_test_snap)) deallocate(y_e_p_a_test_snap) ; allocate(y_e_p_a_test_snap(dim_ene_test_snap))
     if(allocated(y_e_p_a_test_base)) deallocate(y_e_p_a_test_base) ; allocate(y_e_p_a_test_base(dim_ene_test_snap))
  end if

  if (dim_force_test_snap > 0) then
     open (file="test_force.out", unit=finp, action='write')
     if(allocated(y_f_test_snap)) deallocate(y_f_test_snap) ; allocate(y_f_test_snap(dim_force_test_snap))
     if(allocated(y_f_test_base)) deallocate(y_f_test_base) ; allocate(y_f_test_base(dim_force_test_snap))
  end if

  if (dim_stress_test_snap > 0) then
     open (file="test_stress.out", unit=sinp, action='write')
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
    call compute_descriptors(xdesc_i, i)
    !mpi_rangml =0
    !mpi_rangml
    !if (rangml==0) then
     if (config_real(i)%has_energy) then
        i_e_test_snap = i_e_test_snap + 1
        call md_snap_compute_energy(i, pack_opt=.true.)
        y_e_test_snap(i_e_test_snap)=ene_snap
        !debug if (ene_snap == 0.d0) then
        !debug    write(6,*)  rangml, trim(config_real(i)%filename)
        !debug end if
        y_e_test_base(i_e_test_snap)=config_real(i)%energy(2)
        y_e_p_a_test_snap(i_e_test_snap)=ene_snap/dble(config_real(i)%nat)
        y_e_p_a_test_base(i_e_test_snap)=config_real(i)%energy(2)/dble(config_real(i)%nat)
        if (rangml==0) write(einp,'(2e20.10," ", (a), "  ", (a))') ene_snap, y_e_test_base(i_e_test_snap),   config_real(i)%class, config_real(i)%filename
     end if
     if (config_real(i)%has_force) then
        call md_snap_compute_force(i)
        do ix=1,3
          y_f_test_snap(i_f_test_snap+(ix-1)*config_real(i)%nat+1:i_f_test_snap+ix*config_real(i)%nat) =fp_snap(ix,1:config_real(i)%nat)
          y_f_test_base(i_f_test_snap+(ix-1)*config_real(i)%nat+1:i_f_test_snap+ix*config_real(i)%nat) =config_real(i)%force(ix,1:config_real(i)%nat)
        end do
        do ix=1,3*config_real(i)%nat
         if (rangml==0)   write(finp,'(2e20.10," ", (a), "  ", (a))') y_f_test_snap(i_f_test_snap+ix), y_f_test_base(i_f_test_snap+ix), config_real(i)%class,   config_real(i)%filename
        end do
        i_f_test_snap = i_f_test_snap + 3*config_real(i)%nat
     end if

     if (config_real(i)%has_stress) then
        call md_snap_compute_stress(i)
        do ix=1,3
          y_s_test_snap(i_s_test_snap+1:i_s_test_snap+6) =stress_snap(1:6)
          y_s_test_base(i_s_test_snap+1:i_s_test_snap+6) =config_real(i)%stress(1:6)
        end do
        do ix=1,6
          if(rangml==0)  write(sinp,'(2e20.10," ", (a),"  ", (a))') y_s_test_snap(i_s_test_snap+ix), y_s_test_base(i_s_test_snap+ix),  config_real(i)%class, config_real(i)%filename
        end do
        i_s_test_snap = i_s_test_snap + 6
     end if
     !end if !mpi_rangml
     call train_deallocate_desc(i)
     !mpi_rangml end if ! mpi_rangml ==0
     !end if
     call deallocate_ml
  end do

  close(einp)
  close(finp)
  close(sinp)
  !compute the test mpi_rangml==0
  if (rangml==0) call test_error_snap
  !call test_error_snap

return
end subroutine test_snap_routine
!<-------------------------------------------------------------

!$-------------------------------------------------------------
subroutine  test_get_fit_dimensions(iconf)
!$-------------------------------------------------------------
!increment the dimension of the fit depending on the configuration and if
!the energy forces or stress are included in the fit
! Input: iconf
!        config_real(iconf)%
!                          %has_stress has_force has_energy
! Ouptput:    incremented value of
!             dim_data, dim_data_train, dim_data_test
!$-------------------------------------------------------------
use temporary_data_cov, only: dim_data, dim_data_test, dim_data_train
use derived_types, only:config_real
!use snap, only: dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap
implicit none
integer, intent(in) :: iconf

  if (.not.(config_real(iconf)%train))  then

      if (config_real(iconf)%has_force)  then
           dim_data_test=dim_data_test+3*config_real(iconf)%nat
           dim_force_test_snap = dim_force_test_snap + 3*config_real(iconf)%nat
      end if
      if (config_real(iconf)%has_energy) then
           dim_data_test=dim_data_test+1
           dim_ene_test_snap = dim_ene_test_snap + 1
      end if
      if (config_real(iconf)%has_stress) then
           dim_data_test=dim_data_test+6
           dim_stress_test_snap = dim_stress_test_snap + 6
      end if

  else
      if (config_real(iconf)%has_force)  dim_data_train=dim_data_train+3*config_real(iconf)%nat
      if (config_real(iconf)%has_energy) dim_data_train=dim_data_train+1
      if (config_real(iconf)%has_stress) dim_data_train=dim_data_train+6
  end if
  dim_data =  dim_data_test + dim_data_train
return
end subroutine  test_get_fit_dimensions
!<-------------------------------------------------------------




!$-------------------------------------------------------------
subroutine test_error_snap
!$-------------------------------------------------------------
use ml_in_ndm_module, only: rangml
!use snap, only: dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap, &
!                y_e_test_snap, y_f_test_snap, y_s_test_snap, &
!                y_e_test_base, y_f_test_base, y_s_test_base, &
!                y_e_p_a_test_snap, y_e_p_a_test_base, &
!                test_rmse_energy, test_rmse_force, test_rmse_stress, &
!                test_mae_energy, test_mae_force, test_mae_stress
implicit none
real(kind(0.d0)) :: corr, detr, rmse, mae


  if (dim_ene_test_snap > 0) then
      call correlation_coef  (y_e_test_base, y_e_test_snap, dim_ene_test_snap, corr)
      call determination_coef(y_e_test_base, y_e_test_snap, dim_ene_test_snap, detr)
      call rmse_mae          (y_e_test_base, y_e_test_snap, dim_ene_test_snap, rmse, mae)
      if (rangml==0) then
          write(6,'("ML: stat  Total_Energy             test:", 4f15.6)')  corr, detr, rmse, mae
      end if

      call correlation_coef  (y_e_p_a_test_base, y_e_p_a_test_snap, dim_ene_test_snap, corr)
      call determination_coef(y_e_p_a_test_base, y_e_p_a_test_snap, dim_ene_test_snap, detr)
      call rmse_mae          (y_e_p_a_test_base, y_e_p_a_test_snap, dim_ene_test_snap, rmse, mae)
      if (rangml==0) then
          write(6,'("ML: stat  Total_Energy_per_atom    test:", 4f15.6)')  corr, detr, rmse, mae
      end if
      test_rmse_energy = rmse
      test_mae_energy = mae
  end if

  if (dim_force_test_snap > 0) then

      call correlation_coef  (y_f_test_base, y_f_test_snap, dim_force_test_snap, corr)
      call determination_coef(y_f_test_base, y_f_test_snap, dim_force_test_snap, detr)
      call rmse_mae          (y_f_test_base, y_f_test_snap, dim_force_test_snap, rmse, mae)
      if (rangml==0) then
         write(6,'("ML: stat  Force_each_component     test:", 4f15.6)')  corr, detr, rmse, mae
      end if
      test_rmse_force = rmse
      test_mae_force = mae

  end if

  if (dim_stress_test_snap > 0) then

      call correlation_coef  (y_s_test_base, y_s_test_snap, dim_stress_test_snap, corr)
      call determination_coef(y_s_test_base, y_s_test_snap, dim_stress_test_snap, detr)
      call rmse_mae          (y_s_test_base, y_s_test_snap, dim_stress_test_snap, rmse, mae)
      if (rangml==0) then
         write(6,'("ML: stat  Stress_virial            test:", 4f15.6)')  corr, detr, rmse, mae
      end if
      test_rmse_stress = rmse
      test_mae_stress = mae

  end if
  !debug write (*,*) dim_data, dim_data_test,  dim_data_train


return
end subroutine test_error_snap
!<-------------------------------------------------------------


!/------------------------------------------------------------\
!                                                             !
!                           MD   PART                         !
!                                                             !
!\------------------------------------------------------------/


!$-------------------------------------------------------------
subroutine md_allocate_snap_params
!$-------------------------------------------------------------
!MD deivers for SNAP
!used olny for md
use ml_in_ndm_module, only: rangml
use temporary_data_cov, only: dim_xdesc
!use snap, only: w_params
implicit none

if (dim_xdesc == 0) then
    if (rangml==0) write(6,'("ML: dim_xdesc Fatal Error")')
    stop "dim_xdesc is zero in md_allocate_snap"
end if

if (allocated(w_params))      deallocate(w_params)     ; allocate(w_params(1+dim_xdesc,1))

return
end subroutine md_allocate_snap_params
!<-------------------------------------------------------------

!$-------------------------------------------------------------
subroutine md_allocate_snap_desc
!$-------------------------------------------------------------
!MD deivers for SNAP
!used olny for md
use ml_in_ndm_module, only: rangml
use temporary_data_cov, only: xdesc_emd, xdesc_fmd, xdesc_smd, dim_xdesc
!use snap, only: fp_snap
use gen_com_m, only: im, imm
implicit none

if (dim_xdesc == 0) then
    if (rangml==0) write(6,'("ML: dim_xdesc Fatal Error")')
    stop "dim_xdesc is zero in md_allocate_snap"
end if


if (allocated(xdesc_emd))     deallocate(xdesc_emd)    ; allocate(xdesc_emd(1+dim_xdesc))
if (im==0) then
    if (rangml==0) write(6,'("ML: im Fatal Error")')
    stop "im is zero in md_allocate_snap"
end if

if (im /= imm) then
   if (rangml==0) write(6,*) 'im should be equal to imm', im, imm
end if


if (allocated(xdesc_fmd))     deallocate(xdesc_fmd)    ; allocate(xdesc_fmd(1+dim_xdesc,im,3))
if (allocated(xdesc_smd))     deallocate(xdesc_smd)    ; allocate(xdesc_smd(1+dim_xdesc,6))
if (allocated(fp_snap))      deallocate(fp_snap)       ; allocate(fp_snap(3,im))

return
end subroutine md_allocate_snap_desc
!<-------------------------------------------------------------

!$-------------------------------------------------------------
subroutine md_snap_compute_energy(iconf, pack_opt)
!$-------------------------------------------------------------
use temporary_data_cov, only: dim_xdesc, xdesc_emd
use derived_types, only: config_desc, config_real
!use snap, only : w_params, ene_snap
implicit none
integer , intent(in) :: iconf
logical, intent(in), optional :: pack_opt
logical :: lpack_tmp

lpack_tmp=.true.
if (present(pack_opt)) lpack_tmp = pack_opt

xdesc_emd(1) = dble(config_real(iconf)%nat)
if (lpack_tmp) call snap_pack_energy_descriptor(iconf)
!xdesc_emd(2:dim_xdesc+1) = SUM(config_desc(iconf)%energy(1:dim_xdesc,1:config_real(iconf)%nat), dim=2)
xdesc_emd(2:dim_xdesc+1) = config_desc(iconf)%pack_energy(1:dim_xdesc)

ene_snap = dot_product(w_params(:,1), xdesc_emd(:))
!write(*,*)  'test', config_desc(iconf)%energy(1:dim_xdesc,1)
return
end subroutine md_snap_compute_energy
!<-------------------------------------------------------------


!$-------------------------------------------------------------
subroutine md_snap_compute_force(iconf, pack_opt)
!$-------------------------------------------------------------
use temporary_data_cov, only: dim_xdesc, xdesc_fmd
use derived_types, only: config_desc, config_real
!use snap, only : w_params, fp_snap
implicit none
integer , intent(in) :: iconf
integer ::  ik, ix
logical, intent(in), optional :: pack_opt
logical :: lpack_tmp

lpack_tmp=.true.
if (present(pack_opt)) lpack_tmp = pack_opt

!write(*,*) 'sizes in md_snap_compute_force', size(xdesc_fmd,1), size(xdesc_fmd,2), size(xdesc_fmd,3)
xdesc_fmd(1,:,:)=0.d0
if (lpack_tmp) call snap_pack_force_descriptor(iconf)
do ik=1, config_real(iconf)%nat ! force index with respect that derivative ik
   xdesc_fmd(2:dim_xdesc+1,ik,1:3)  = config_desc(iconf)%pack_force(1:dim_xdesc,1:3,ik)
   !write(*,*)  xdesc_fmd(5,ik,1)
   do ix=1,3
      fp_snap(ix,ik)=dot_product(w_params(:,1), xdesc_fmd(:,ik,ix))
   end do
end do

return
end subroutine md_snap_compute_force
!<-------------------------------------------------------------


!$-------------------------------------------------------------
subroutine md_snap_compute_stress(iconf, pack_opt)
!$-------------------------------------------------------------
!$ This subroutine has strong interaction with neighbours subroutines
!$ The one proposed by Milady and the one proposed by NDM depending if the cell box is small or large.
!$ small=.true. ->  MiLaDy
!$ small=.false. -> NDM
!NDM_interaction
!MiLaDy_interaction
use temporary_data_cov, only: dim_xdesc, xdesc_smd
use derived_types, only: config_desc
!use snap, only : w_params, stress_snap
implicit none
integer , intent(in) :: iconf
integer :: ix
logical, intent(in), optional :: pack_opt
logical :: lpack_tmp

lpack_tmp=.true.
if (present(pack_opt)) lpack_tmp = pack_opt

if (lpack_tmp)   call snap_pack_stress_descriptor(iconf)
!write(*,*) "md_snap_compute_stress", iconf, config_real(iconf)%volume
xdesc_smd(1,1:6)=0.d0
xdesc_smd(2:dim_xdesc+1,1:6) = config_desc(iconf)%pack_stress(1:dim_xdesc,1:6)

do ix=1,6
   stress_snap(ix)=dot_product(w_params(:,1), xdesc_smd(:,ix))
end do

return
end subroutine md_snap_compute_stress
!<-------------------------------------------------------------



!/------------------------------------------------------------\
!                                                             !
!                 General/restart/input  PART                 !
!                                                             !
!\------------------------------------------------------------/




!$-------------------------------------------------------------
subroutine dreal_inverse (dim_mat,mat,mat_inv)
!$-------------------------------------------------------------
use ml_in_ndm_module, only : rangml
implicit none
integer, intent(in) :: dim_mat
real(kind=kind(0.d0)), dimension(dim_mat,dim_mat), intent(in)  :: mat
real(kind=kind(0.d0)), dimension(dim_mat,dim_mat), intent(out) :: mat_inv
! driver inversion
integer :: dim_work,dim_ipiv,info
real(kind=kind(0.d0)), dimension(:), allocatable :: work
integer, dimension(:),allocatable  :: ipiv
real(kind=kind(0.d0)), dimension(dim_mat, dim_mat) :: id_matrix, c_matrix
integer :: i
!  if (rangml==0) write(6,*) 'Trying brute inversion ...'
  id_matrix(:,:)=0.d0
  do i=1,dim_mat
     id_matrix(i,i)=1.d0
  end do

  mat_inv(:,:) = mat(:,:)
  dim_work=max(1, 2*dim_mat)
  dim_ipiv=max(1, dim_mat)
  info=0
  if (allocated(work)) deallocate(work) ; allocate (work(dim_work))
  if (allocated(ipiv)) deallocate(ipiv) ; allocate (ipiv(dim_ipiv))
  call dsytrf( 'L', dim_mat, mat_inv, dim_mat, ipiv, work, dim_work, info ) ! factorization
  !debug if (rangml==0) write(6,'("ML: dsytrf info for factorization...:", i6)') info
  if (info /= 0) then
     if (rangml==0) write(6,'("ML: WARNING problems factorization  dsytrf info is ...:", i6)') info
  end if
  !dim_work=max(1, 2*dim_mat)
  !if (allocated(work)) deallocate(work) ; allocate (work(dim_work))
  call dsytri( 'L', dim_mat, mat_inv, dim_mat, ipiv, work, info )  !inversion
  !debug if (rangml==0) write(6,'("ML: dsytri info for inversion.......:", i6)') info
  if (info /= 0) then
     if (rangml==0) write(6,'("ML: WARNING problems inversion dsytri info is ...:", i6)') info
  end if
  call dsymm('L','L',dim_mat, dim_mat, 1.d0, mat_inv, dim_mat, id_matrix,dim_mat, 0.d0,c_matrix,dim_mat)
  mat_inv=c_matrix
return
end subroutine dreal_inverse
!<-------------------------------------------------------------


subroutine dreal_matmul (A, nrow_A, ncol_A, B, nrow_B, ncol_B, C, nrow_C, ncol_C)
!--------------------------------------------------------------
! perform C = A*B using dgemm ....
!
!On input:
!           the matrix A  of size  nrow_A x ncol_A
!           the matrix B  of size  nrow_B x ncol_B
!
!On output
!           the matrix C  of size nrow_C, ncol_C
!Observation:
!            ncol_A = nrow_B
!            nrow_A = nrow_C
!            ncol_B = nrow_C
!-------------------------------------------------------------

implicit none
integer, intent(in) :: nrow_A, ncol_A, nrow_B, ncol_B, nrow_C, ncol_C
real(kind=kind(0.d0)), dimension(nrow_A,ncol_A), intent(in)  :: A
real(kind=kind(0.d0)), dimension(nrow_B,ncol_B), intent(in)  :: B
real(kind=kind(0.d0)), dimension(nrow_C,ncol_C), intent(out)  :: C

if ((ncol_A /= nrow_B) .or. ( ncol_B /= ncol_C .or. (nrow_A /= nrow_C) )) then
   write(6,*) "ML: multiplication wrong. stop"
   write(6,*) "nrow_A, ncol_A", nrow_A, ncol_A
   write(6,*) "nrow_B, ncol_B", nrow_B, ncol_B
   write(6,*) "nrow_C, ncol_C", nrow_C, ncol_C

   stop "dreal_matmul wrong dimensions"
end if

call  dgemm ( 'N', 'N', nrow_A, ncol_B, ncol_A, 1.d0, A, nrow_A, B, nrow_B, 0.d0, C, nrow_C )

!deallocate (A, B)
end subroutine dreal_matmul

!$-------------------------------------------------------------
subroutine read_parameters_snap
!$-------------------------------------------------------------
! This subroutine read the parameters file of the snap potential
use temporary_data_cov, only : dim_xdesc
use ml_in_ndm_module, only: char_desc, rangml
!use snap, only: w_params
implicit none
character(len=80) :: file_params
integer :: dim_pot,i
logical :: ok

file_params=trim(adjustl(char_desc//'_snap1_params.pot'))
inquire (file=file_params, exist=ok)
if (ok) then
  open (file=file_params, unit=51, action='read')
  read(51,*) dim_pot
  if (dim_pot /= (dim_xdesc+1))  then
    if (rangml==0) write (6,*) 'ML: Fatal the parameters file and the descrptors no have the same dimension'
    stop "read_parameters_snap dimension wrong"
  end if
  do i=1,dim_pot
     read(51,*) w_params(i,1)
  end do
  close(51)
else
  if (rangml==0) write(6,*) 'ML: Fatal Error: the potential file is not there. The present name: ', file_params
  stop "read_parameters_snap"
endif

return
end subroutine read_parameters_snap
!<-------------------------------------------------------------
end module snap
