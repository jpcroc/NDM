
!/------------------------------------------------------------\
!                                                             !
!                    WEIGHTS OPTIMIZATION                     !
!                                                             !
!\------------------------------------------------------------/

module data_type
  implicit none
  integer(kind=4), parameter :: IB=4, RP=8
  integer, save :: icall
  integer :: nunit
end module data_type


!$-------------------------------------------------------------
subroutine snap_optimize_weights
!$-------------------------------------------------------------
use ml_in_ndm_module, only: rangml, max_iter_optimize_weights, optimize_ga_population, optimize_weights_Le
use snap, only:  dim_weights_function, dim_weights_function_full, &
                 upper_weights, lower_weights
use data_type, only :icall, nunit
implicit none
real(kind(0.d0)) :: VTR,F_XC, CR_XC, F_CR
integer :: NP, itermax, strategy, refresh, iwrite, nfeval
real(kind(0.d0)), dimension(dim_weights_function_full) :: bestmem_XC
real(kind(0.d0)) :: bestval
integer, dimension(3), parameter:: method=(/0, 1, 0/)
integer :: ic
external function_to_min

     ! size the dimension of J(w_1, w_2 ....)
    !debug if (rangml==0) write²(6,*) 'enter get_dim'
     call get_dim_snap_objective_function
     ! From db put the values of w_1, w_2 ...
    !debug if (rangml==0) write(6,*) 'enter init '
     call init_weigths_from_db_in_snap_objective_function
     ! From w_1, w_2  fill the matrix of weights - the one, of dimension dim_data_train,
     ! Fill tmp_weights_snap the used by Amat in order to get the parameters


     !oldversion call fill_tmp_weigths_snap_with_w
     !oldversion weights_snap = tmp_weights_snap
     !oldversion call train_snap_get_parameters
     !oldversion call train_objective_function
     !oldversion write(*,*) 'Jfunc', Jfunc

     icall = 0
     !the expected fitness value to reach.
     VTR=1.d-04
     ! Population size.old was 40
     NP=optimize_ga_population
     ! The maximum number of iteration.
     itermax=max_iter_optimize_weights
     !Mutation scaling factor for real decision parameters.
     F_XC=0.8d0
     !Crossover factor for real decision parameters.
     CR_XC=0.8d0
     !The strategy of the mutation operations is used in HDE 1 to 6
     strategy=6
     !The unit specfier for writing to an external data file.
     iwrite=7
      !The intermediate output will be produced after "refresh"
      !iterations. No intermediate output will be produced if
      !"refresh < 1".
      refresh=200
      F_CR=0.8d0
      !method(1) = 0, Fixed mutation scaling factors (F_XC)
      !          = 1, Random mutation scaling factors F_XC=[0, 1]
      !          = 2, Random mutation scaling factors F_XC=[-1, 1]
      !method(2) = 1, Random combined factor (F_CR) used for strategy = 6
      !               in the mutation operation
      !          = other, fixed combined factor provided by the user
      !method(3) = 1, Saving results in a data file.
      !          = other, displaying results only.
     nunit=73
     open(unit=nunit, file='optimization_error.dat', status='unknown')
     !debug if (rangml==0) write(6,*) 'Enter DE_Fortran'
     call DE_Fortran90(function_to_min, dim_weights_function_full, lower_weights, upper_weights, &
                        VTR, NP, itermax, F_XC,CR_XC, strategy, refresh, iwrite, bestmem_XC, &
				bestval, nfeval, F_CR, method)
        if (rangml==0) write(6,'("ML: Number of function evaluation for genetic algo minimization ...:", i8)') nfeval
        if (rangml==0) write(6,'("ML: The best value of the objective function ......................:", e20.10)') bestval
        if (rangml==0) write(6,'("ML: The set of weigths ............................................:", i8)')
        do ic=1, dim_weights_function
          if (rangml==0) write(6,'(e20.10)') bestmem_XC(ic)
        end do
        if (rangml==0) write(6,'("ML: The best lambda_krr L2 found .....................................:", e20.10)') bestmem_XC(dim_weights_function_full)
        if (optimize_weights_Le) then
          if (rangml==0) write(6,'("ML: The best lambda_krr L1 found .....................................:", 2e20.10)') bestmem_XC(dim_weights_function_full-1)
        else
        end if

        call ga_best_evaluation (bestmem_XC)
return
end subroutine snap_optimize_weights
!<-------------------------------------------------------------

subroutine ga_best_evaluation (xval)
  use ml_in_ndm_module, only : lambda_krr, optimize_weights_Le
  use snap, only : tmp_weights, weights_snap, dim_weights_function, dim_weights_function_full
  implicit none
  real(kind(0.d0)), dimension(dim_weights_function_full), intent(in) :: xval

  tmp_weights(1:dim_weights_function) = xval(1:dim_weights_function)
  lambda_krr = xval (dim_weights_function_full)
  call fill_tmp_weigths_snap_with_w
  call train_snap_get_parameters

  return
end subroutine ga_best_evaluation


subroutine function_to_min(xval,objval)
use data_type, only : icall, nunit
use ml_in_ndm_module, only: rangml, debug, factor_energy_error, factor_force_error, factor_stress_error, lambda_krr, &
                            optimize_weights_L1, optimize_weights_L2, optimize_weights_Le, lambda_krr_fake, lambda_krr_2
use snap, only: dim_weights_function, dim_weights_function_full, tmp_weights, weights_snap, tmp_weights_snap
implicit none
real(kind(0.d0)), dimension(dim_weights_function+1), intent(in) :: xval
real(kind(0.d0)), intent(out) :: objval
real(kind(0.d0)) :: mae_energy, mae_force, mae_stress
     tmp_weights(1:dim_weights_function) = xval(1:dim_weights_function)
    !debug if (rangml==0) write(6,'("00000",i6,e20.10)') icall, lambda_krr
     call fill_tmp_weigths_snap_with_w
    !debug if (rangml==0) write(6,'("11111", i6, 3e20.10)') icall, lambda_krr,  minval(tmp_weights_snap(:)), maxval(tmp_weights_snap)

     weights_snap(:) = tmp_weights_snap(:)
     !debug if (rangml==0) write(6,*) icall, 'before train_snap'
     lambda_krr = lambda_krr_fake
     lambda_krr = xval (dim_weights_function_full)
     call train_snap_get_parameters
     lambda_krr = xval (dim_weights_function_full)
     if (optimize_weights_Le) lambda_krr_2 = xval (dim_weights_function_full-1)
     !debug if (rangml==0) write(6,'("22222", i6,3e20.10)') icall, lambda_krr,  minval(tmp_weights_snap(:)), maxval(tmp_weights_snap)
     !old call train_objective_function
     !old objval = Jfunc
     !debug if (rangml==0) write(6,*) icall, 'before train_minimize'
     !debug if (rangml==0) write(6,'("33333", 6e20.10)') factor_energy_error, mae_energy, factor_force_error, mae_force, factor_stress_error, mae_stress
     call train_minimize_error (mae_energy, mae_force,mae_stress)
     !debug if (rangml==0) write(6,'("44444", 6e20.10)') factor_energy_error, mae_energy, factor_force_error, mae_force, factor_stress_error, mae_stress
     icall = icall + 1
     objval= factor_energy_error*mae_energy + factor_force_error*mae_force + factor_stress_error*mae_stress
     if (lambda_krr >=0) then
        if (optimize_weights_L2) objval = objval + lambda_krr * sum(weights_snap(:)**2)
        if (optimize_weights_L1) objval = objval + lambda_krr * sum(dabs(weights_snap(:)))
        if (optimize_weights_Le) objval = objval + lambda_krr_2 * sum(dabs(weights_snap(:))) + lambda_krr * sum(weights_snap(:)**2)
       !else
       ! objval= factor_energy_error*mae_energy + factor_force_error*mae_force + factor_stress_error*mae_stress
     end if
     !debug if (optimize_weights_Le) then
     !debug       if (rangml==0) write(6,'("ML:   Jfunc optimization :",i7, 4e20.10)') icall, lambda_krr, lambda_krr_2, objval, sum(weights_snap(:)**2)
     !debug else
    !debug     if (rangml==0) write(6,'("ML:   Jfunc optimization :",i7, 3e20.10)') icall, lambda_krr,  objval, sum(weights_snap(:)**2)
    !debug  end if
     if (debug) then

      if (optimize_weights_Le) then
          if (rangml==0) write(6,'("ML:   Jfunc optimization :",i7, 3e20.10)') icall, lambda_krr, lambda_krr_2, objval, sum(weights_snap(:)**2)
      else
        if (rangml==0) write(6,'("ML:   Jfunc optimization :",i7, 3e20.10)') icall, lambda_krr, objval, sum(weights_snap(:)**2)
      end if

     end if
     write(nunit,'(i7,4e20.10)') icall, objval, mae_energy, mae_force, mae_stress
return
end subroutine function_to_min




!$-------------------------------------------------------------
subroutine get_dim_snap_objective_function
!$-------------------------------------------------------------
use ml_in_ndm_module, only : iconf_data, rangml, no_class_weights, optimize_weights_L1, optimize_weights_L2, optimize_weights_Le
use snap, only : dim_weights_function, dim_weights_function_full, map_weights_in_db, map_db_in_weights
use derived_types, only : db_model, config_real
implicit none
integer :: i_w,idb, jc, i, iclass
integer :: ienergy, iforce, istress
logical :: l_energy, l_force, l_stress, skip_line

l_energy=.false.
l_force=.false.
l_stress=.false.

i_w=0
do idb=1,size(db_model)
   !some classes are skipped
   skip_line=.false.
   do iclass=1,size(no_class_weights)
      if (db_model(idb)%class==no_class_weights(iclass)) then
        skip_line=.true.
        db_model(idb)%has_optimize_weights_energy=.false.
        db_model(idb)%has_optimize_weights_force =.false.
        db_model(idb)%has_optimize_weights_stress=.false.
        cycle
      end if
   end do
   if (skip_line) cycle
   !end class selection

   ienergy=0
   iforce=0
   istress=0
   !depending of fitting parameter T/F is in the db line the ienergy, iforce, istress is set to 1/0
   if (db_model(idb)%has_db_energy) then
         ienergy=1
   end if
   if (db_model(idb)%has_db_force) then
         iforce=1
   end if
   if (db_model(idb)%has_db_stress) then
         istress=1
   end if

   !if all the componenets are false the ienergy, iforce. istress are putted to zero whatever is the initial value

   do jc=1, iconf_data
      !selecting only trainning data ...
      if (.not.(config_real(jc)%train)) cycle
         ! in which line the jc conf is located ....
         if (config_real(jc)%db_line==idb) then
             l_energy=l_energy.or.config_real(jc)%has_energy
             l_force=l_force.or.config_real(jc)%has_force
             l_stress=l_stress.or.config_real(jc)%has_stress
         end if
   end do

   if (.not.l_energy) ienergy=0
   if (.not.l_force)  iforce=0
   if (.not.l_stress) istress=0
   i_w = i_w + ienergy + iforce + istress
   db_model(idb)%has_optimize_weights_energy=db_model(idb)%has_db_energy.and.l_energy
   db_model(idb)%has_optimize_weights_force =db_model(idb)%has_db_force.and.l_force
   db_model(idb)%has_optimize_weights_stress=db_model(idb)%has_db_stress.and.l_stress
end do

dim_weights_function=i_w

if (optimize_weights_L1) dim_weights_function_full=dim_weights_function+1
if (optimize_weights_L2) dim_weights_function_full=dim_weights_function+1
if (optimize_weights_Le) dim_weights_function_full=dim_weights_function+2

if (allocated(map_weights_in_db)) deallocate(map_weights_in_db) ; allocate(map_weights_in_db(dim_weights_function))
if (allocated(map_db_in_weights)) deallocate(map_db_in_weights) ; allocate(map_db_in_weights(size(db_model)))

!if (rangml==0) then
!    write(6,'("ML: the dimension of the objective function, to optimize, is: ", i7)') dim_weights_function
!end if


i_w=0
do i=1,size(db_model)

   !some classes are skipped
   skip_line=.false.
   do iclass=1,size(no_class_weights)
      if (db_model(i)%class==no_class_weights(iclass)) then
        skip_line=.true.
        cycle
      end if
   end do
   if (skip_line) cycle
   !end class selection


   if (db_model(i)%has_optimize_weights_energy) then
         i_w = i_w + 1
         map_weights_in_db(i_w)%db_line=i
         map_weights_in_db(i_w)%has_energy=.true.
         map_db_in_weights(i)%i_e=i_w
         !map_db_in_weights(i)%has_energy=.true.
   end if
   if (db_model(i)%has_optimize_weights_force) then
         i_w = i_w + 1
         map_weights_in_db(i_w)%db_line=i
         map_weights_in_db(i_w)%has_force=.true.
         map_db_in_weights(i)%i_f=i_w
   end if
   if (db_model(i)%has_optimize_weights_stress) then
         i_w = i_w + 1
         map_weights_in_db(i_w)%db_line=i
         map_weights_in_db(i_w)%has_stress=.true.
         map_db_in_weights(i)%i_s=i_w
   end if
end do

if (i_w /= dim_weights_function) then
   if(rangml==0)  write(6,'("There are some problems in the definition of the dimension of objective function, old and new :", 2i6)') i_w, dim_weights_function
   stop 'dimension problem in get_dim_snap_objective_function'
end if
return
end subroutine get_dim_snap_objective_function
!<-------------------------------------------------------------

!$-------------------------------------------------------------
subroutine init_weigths_from_db_in_snap_objective_function
!$-------------------------------------------------------------
use temporary_data_cov, only:  dim_data_train
use snap, only: dim_weights_function, dim_weights_function_full,  tmp_weights, upper_weights, lower_weights, &
                map_weights_in_db, tmp_weights_snap
use derived_types, only:  db_model
use ml_in_ndm_module, only : rangml, lambda_krr, min_lambda_krr, max_lambda_krr, optimize_weights_Le
implicit none
integer :: i
real(kind(0.d0)) :: tmp_ini, tmp_fin

if (allocated(tmp_weights_snap))     deallocate(tmp_weights_snap)    ; allocate(tmp_weights_snap(dim_data_train))


if (allocated(tmp_weights))   deallocate(tmp_weights)   ; allocate(tmp_weights(dim_weights_function_full))
if (allocated(lower_weights)) deallocate(lower_weights) ; allocate(lower_weights(dim_weights_function_full))
if (allocated(upper_weights)) deallocate(upper_weights) ; allocate(upper_weights(dim_weights_function_full))

do i=1,dim_weights_function
   if (map_weights_in_db(i)%has_energy)  then
        tmp_ini=db_model(map_weights_in_db(i)%db_line)%w_e
        tmp_fin=db_model(map_weights_in_db(i)%db_line)%w_e_end
    end if

   if (map_weights_in_db(i)%has_force)  then
        tmp_ini=db_model(map_weights_in_db(i)%db_line)%w_f
        tmp_fin=db_model(map_weights_in_db(i)%db_line)%w_f_end
    end if


   if (map_weights_in_db(i)%has_stress)  then
        tmp_ini=db_model(map_weights_in_db(i)%db_line)%w_s
        tmp_fin=db_model(map_weights_in_db(i)%db_line)%w_s_end
    end if

   !tmp_weights(i) = (tmp_ini + tmp_fin)/2.d0
   tmp_weights(i) = tmp_ini
   lower_weights(i)= tmp_ini
   upper_weights(i)= tmp_fin
end do

!tmp_weights(dim_weights_function+1)=(min_lambda_krr + max_lambda_krr)/2.d0
if (optimize_weights_Le) then
  tmp_weights(dim_weights_function_full-1)=min_lambda_krr
  lower_weights(dim_weights_function_full-1)=min_lambda_krr
  upper_weights(dim_weights_function_full-1)=max_lambda_krr


  tmp_weights(dim_weights_function_full)=min_lambda_krr
  lower_weights(dim_weights_function_full)=min_lambda_krr
  upper_weights(dim_weights_function_full)=max_lambda_krr

else
  tmp_weights(dim_weights_function_full)=min_lambda_krr
  lower_weights(dim_weights_function_full)=min_lambda_krr
  upper_weights(dim_weights_function_full)=max_lambda_krr
end if

if (rangml==0) then
   write(6,*) 'ML:-------------------Genetic Optimization----------------------'
   if (lambda_krr > 0.d0) then
     write(6,*) 'ML: The optimization on weigths is performed in the same time with L2 norm min and max values if lambda: ', min_lambda_krr, max_lambda_krr
   end if
   write(6,*) 'ML: The size of the function to be optimized: ', dim_weights_function_full
end if


return
end subroutine init_weigths_from_db_in_snap_objective_function
!<-------------------------------------------------------------


!$-------------------------------------------------------------
subroutine fill_tmp_weigths_snap_with_w
!$-------------------------------------------------------------
use ml_in_ndm_module, only: iconf_data, rangml!, no_class_weights
use temporary_data_cov, only : dim_data_train
use derived_types, only: config_real, db_model
use snap, only: map_db_in_weights, tmp_weights, tmp_weights_snap
implicit none
integer :: i, itmp, itmp_f
!integer ::  jc
!logical :: lskip

itmp=0
  do i =1,iconf_data
     if (.not.(config_real(i)%train)) cycle
    !here is the problem

     !lskip=.false.
     !do jc=1,size(no_class_weights)
      !    if (config_real(i)%class==no_class_weights(jc)) then
    !         lskip=.true.
    !         cycle
    !       end if
     !end do
     !if (lskip) cycle

     if (config_real(i)%has_energy) then
        itmp=itmp+1
          if (db_model(config_real(i)%db_line)%has_optimize_weights_energy) then
              tmp_weights_snap(itmp)=tmp_weights(map_db_in_weights(config_real(i)%db_line)%i_e)
            else
              tmp_weights_snap(itmp)=config_real(i)%w_e
          end if
     end if

     if (config_real(i)%has_force) then
        itmp_f=3*config_real(i)%nat

          if (db_model(config_real(i)%db_line)%has_optimize_weights_force) then
              tmp_weights_snap(itmp+1:itmp+itmp_f)=tmp_weights(map_db_in_weights(config_real(i)%db_line)%i_f)
            else
              tmp_weights_snap(itmp+1:itmp+itmp_f)=config_real(i)%w_f
          end if
        itmp=itmp+itmp_f
     end if

     if (config_real(i)%has_stress) then
          if (db_model(config_real(i)%db_line)%has_optimize_weights_stress) then
              tmp_weights_snap(itmp+1:itmp+6)=tmp_weights(map_db_in_weights(config_real(i)%db_line)%i_s)
            else
              tmp_weights_snap(itmp+1:itmp+6)=config_real(i)%w_s
          end if
         itmp=itmp+6
     end if
   end do

if (itmp /= dim_data_train) then
   if(rangml==0)  write(6,'("There are some problems in tmp_weights_snap vector, old and new dimension:", 2i6)') dim_data_train, itmp
   stop 'dimension problem in fill_tmp_weigths_snap_with_w'
end if
return
end subroutine fill_tmp_weigths_snap_with_w
!<-------------------------------------------------------------






!$-------------------------------------------------------------
subroutine train_minimize_error (mae_energy, mae_force, mae_stress)
!$-------------------------------------------------------------
use ml_in_ndm_module, only : no_class_weights
use snap, only : w_params, Amat, ene_snap, fit_snap!, &
                 !i_e_train_snap, i_f_train_snap, i_s_train_snap, &
                 !y_e_train_base, y_f_train_base, y_s_train_base, &
                 !y_e_train_snap, y_f_train_snap, y_s_train_snap, &
                 !y_e_p_a_train_snap, y_e_p_a_train_base, Jfunc, tmp_weights_snap

use temporary_data_cov, only:  yfunc_train, dim_data_train
use derived_types, only: config_real
implicit none
real(kind(0.d0)),intent(out)  :: mae_energy, mae_force, mae_stress
integer :: idata, i_w, jc, dim_y_force_error, dim_y_energy_error, dim_y_stress_error
real(kind(0.d0)) :: force_snap,stress_snap
!old things
!character(len=2), dimension(2) :: class_energy_error
!character(len=2), dimension(2) :: class_force_error
real(kind(0.d0)), dimension(:), allocatable :: y_energy_error_snap, y_energy_error_base, &
                                               y_force_error_snap , y_force_error_base , &
                                               y_stress_error_snap, y_stress_error_base
real(kind(0.d0)) :: rmse_force, rmse_energy, rmse_stress
!if (config_real(fit_snap(idata)%iconf)%class=="02") then
!write(*,*) 'objective_test', idata,  config_real(fit_snap(idata)%iconf)%class, config_real(fit_snap(idata)%iconf)%filename
!end if

!class_energy_error(1)='02'
!class_energy_error(2)='06'

!class_force_error(1)='03'
!class_force_error(2)='06'

i_w=0
do idata=1,dim_data_train
   if (fit_snap(idata)%energy) then

       do jc=1,size(no_class_weights)
          if (config_real(fit_snap(idata)%iconf)%class==no_class_weights(jc)) cycle
          i_w= i_w+1
       end do

  end if
end do

dim_y_energy_error=i_w


i_w=0
do idata=1,dim_data_train
   if (fit_snap(idata)%force) then

       do jc=1,size(no_class_weights)
          if (config_real(fit_snap(idata)%iconf)%class==no_class_weights(jc)) cycle
          i_w= i_w+1
       end do

  end if
end do

dim_y_force_error=i_w


i_w=0
do idata=1,dim_data_train
   if (fit_snap(idata)%stress) then

       do jc=1,size(no_class_weights)
          if (config_real(fit_snap(idata)%iconf)%class==no_class_weights(jc)) cycle
          i_w= i_w+1
       end do

  end if
end do

dim_y_stress_error=i_w



if (allocated (y_force_error_base))   deallocate (y_force_error_base) ; allocate(y_force_error_base(dim_y_force_error))
if (allocated (y_force_error_snap))   deallocate (y_force_error_snap) ; allocate(y_force_error_snap(dim_y_force_error))

if (allocated (y_energy_error_base))  deallocate (y_energy_error_base) ; allocate(y_energy_error_base(dim_y_energy_error))
if (allocated (y_energy_error_snap))  deallocate (y_energy_error_snap) ; allocate(y_energy_error_snap(dim_y_energy_error))


if (allocated (y_stress_error_base))  deallocate (y_stress_error_base) ; allocate(y_stress_error_base(dim_y_stress_error))
if (allocated (y_stress_error_snap))  deallocate (y_stress_error_snap) ; allocate(y_stress_error_snap(dim_y_stress_error))


i_w=0
do idata=1,dim_data_train
  if (fit_snap(idata)%energy) then
       do jc=1,size(no_class_weights)
          if (config_real(fit_snap(idata)%iconf)%class==no_class_weights(jc)) cycle
          ene_snap = dot_product(w_params(:,1), Amat(:,idata))
          i_w= i_w+1
          y_energy_error_snap(i_w) = ene_snap
          y_energy_error_base(i_w) = yfunc_train(idata)
       end do
  end if
end do

if (i_w>0) then
  call rmse_mae (y_energy_error_base, y_energy_error_snap, dim_y_energy_error, rmse_energy, mae_energy)
else
  rmse_force=0.d0
  mae_force=0.d0
end if

i_w=0
do idata=1,dim_data_train
  if (fit_snap(idata)%force) then
       do jc=1,size(no_class_weights)
          if (config_real(fit_snap(idata)%iconf)%class==no_class_weights(jc)) cycle
          force_snap = dot_product(w_params(:,1), Amat(:,idata))
          i_w= i_w+1
          y_force_error_snap(i_w) = force_snap
          y_force_error_base(i_w) = yfunc_train(idata)
       end do
  end if
end do

if (i_w>0) then
  call rmse_mae (y_force_error_base, y_force_error_snap, dim_y_force_error, rmse_force, mae_force)
else
  rmse_force=0.d0
  mae_force=0.d0
end if


i_w=0
do idata=1,dim_data_train
  if (fit_snap(idata)%stress) then
       do jc=1,size(no_class_weights)
          if (config_real(fit_snap(idata)%iconf)%class==no_class_weights(jc)) cycle
          stress_snap = dot_product(w_params(:,1), Amat(:,idata))
          i_w= i_w+1
          y_stress_error_snap(i_w) = force_snap
          y_stress_error_base(i_w) = yfunc_train(idata)
       end do
  end if
end do

if (i_w>0) then
  call rmse_mae (y_stress_error_base, y_stress_error_snap, dim_y_stress_error, rmse_stress, mae_stress)
else
  rmse_force=0.d0
  mae_force=0.d0
end if


return
end subroutine train_minimize_error
!<-------------------------------------------------------------






! Berkley implementation from Taiwan ....


subroutine DE_Fortran90(obj, Dim_XC, XCmin, XCmax, VTR, NP, itermax, F_XC, &
           CR_XC, strategy, refresh, iwrite, bestmem_XC, bestval, nfeval, &
		   F_CR, method)
!.......................................................................
!
! Differential Evolution for Optimal Control Problems
!
!.......................................................................
!  This Fortran 90 program translates from the original MATLAB
!  version of differential evolution (DE). This FORTRAN 90 code
!  has been tested on Compaq Visual Fortran v6.1.
!  Any users new to the DE are encouraged to read the article of Storn and Price.
!
!  Refences:
!  Storn, R., and Price, K.V., (1996). Minimizing the real function of the
!    ICEC'96 contest by differential evolution. IEEE conf. on Evolutionary
!    Comutation, 842-844.
!
!  This Fortran 90 program written by Dr. Feng-Sheng Wang
!  Department of Chemical Engineering, National Chung Cheng University,
!  Chia-Yi 621, Taiwan, e-mail: chmfsw@ccunix.ccu.edu.tw
!.........................................................................
!                obj : The user provided file for evlauting the objective function.
!                      subroutine obj(xc,fitness)
!                      where "xc" is the real decision parameter vector.(input)
!                            "fitness" is the fitness value.(output)
!             Dim_XC : Dimension of the real decision parameters.
!      XCmin(Dim_XC) : The lower bound of the real decision parameters.
!      XCmax(Dim_XC) : The upper bound of the real decision parameters.
!                VTR : The expected fitness value to reach.
!                 NP : Population size.
!            itermax : The maximum number of iteration.
!               F_XC : Mutation scaling factor for real decision parameters.
!              CR_XC : Crossover factor for real decision parameters.
!           strategy : The strategy of the mutation operations is used in HDE.
!            refresh : The intermediate output will be produced after "refresh"
!                      iterations. No intermediate output will be produced if
!                      "refresh < 1".
!             iwrite : The unit specfier for writing to an external data file.
! bestmen_XC(Dim_XC) : The best real decision parameters.
!              bestval : The best objective function.
!             nfeval : The number of function call.
!         method(1) = 0, Fixed mutation scaling factors (F_XC)
!                   = 1, Random mutation scaling factors F_XC=[0, 1]
!                   = 2, Random mutation scaling factors F_XC=[-1, 1]
!         method(2) = 1, Random combined factor (F_CR) used for strategy = 6
!                        in the mutation operation
!                   = other, fixed combined factor provided by the user
!         method(3) = 1, Saving results in a data file.
!                   = other, displaying results only.

	 use data_type, only : IB, RP
	 implicit none
	 integer(kind=IB), intent(in) :: NP, Dim_XC, itermax, strategy,   &
	                                 iwrite, refresh
     real(kind=RP), intent(in) :: VTR, CR_XC
	 real(kind=RP) :: F_XC, F_CR
	 real(kind=RP), dimension(Dim_XC), intent(in) :: XCmin, XCmax
     real(kind=RP), dimension(Dim_XC), intent(inout) :: bestmem_XC
	 real(kind=RP), intent(out) :: bestval
	 integer(kind=IB), intent(out) :: nfeval
     real(kind=RP), dimension(NP,Dim_XC) :: pop_XC, bm_XC, mui_XC, mpo_XC,   &
	                                        popold_XC, rand_XC, ui_XC
     integer(kind=IB) :: i, ibest, iter
     integer(kind=IB), dimension(NP) :: rot, a1, a2, a3, a4, a5, rt
     integer(kind=IB), dimension(4) :: ind
     real(kind=RP) :: tempval
	 real(kind=RP), dimension(NP) :: val
     real(kind=RP), dimension(Dim_XC) :: bestmemit_XC
     real(kind=RP), dimension(Dim_XC) :: rand_C1
	 integer(kind=IB), dimension(3), intent(in) :: method
	 external  obj
	 intrinsic max, min, random_number, mod, abs, any, all, maxloc
	 interface
        function randperm(num)
		   use data_type, only : IB
		   implicit none
	       integer(kind=IB), intent(in) :: num
           integer(kind=IB), dimension(num) :: randperm
        end function randperm
     end interface
 !!-----Initialize a population --------------------------------------------!!

        pop_XC=0.0_RP
        do i=1,NP
           call random_number(rand_C1)
           pop_XC(i,:)=XCmin+rand_C1*(XCmax-XCmin)
		end do

!!--------------------------------------------------------------------------!!

!!------Evaluate fitness functions and find the best member-----------------!!
     val=0.0_RP
     nfeval=0
     ibest=1
     call obj(pop_XC(1,:), val(1))
     bestval=val(1)
     nfeval=nfeval+1
     do i=2,NP
        call obj(pop_XC(i,:), val(i))
        nfeval=nfeval+1
        if (val(i) < bestval) then
           ibest=i
	       bestval=val(i)
        end if
     end do
     bestmemit_XC=pop_XC(ibest,:)
     bestmem_XC=bestmemit_XC
!!--------------------------------------------------------------------------!!

     bm_XC=0.0_RP
     rot=(/(i,i=0,NP-1)/)
     iter=1
!!------Perform evolutionary computation------------------------------------!!

     do while (iter <= itermax)
        popold_XC=pop_XC

!!------Mutation operation--------------------------------------------------!!
        ind=randperm(4)
        a1=randperm(NP)
        rt=mod(rot+ind(1),NP)
        a2=a1(rt+1)
        rt=mod(rot+ind(2),NP)
        a3=a2(rt+1)
        rt=mod(rot+ind(3),NP)
        a4=a3(rt+1)
        rt=mod(rot+ind(4),NP)
        a5=a4(rt+1)
        bm_XC=spread(bestmemit_XC, DIM=1, NCOPIES=NP)

!----- Generating a random sacling factor--------------------------------!
		select case (method(1))
		case (1)
		   call random_number(F_XC)
		case(2)
		   call random_number(F_XC)
           F_XC=2.0_RP*F_XC-1.0_RP
		end select

!---- select a mutation strategy-----------------------------------------!
		select case (strategy)
        case (1)
           ui_XC=bm_XC+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))

	    case default
           ui_XC=popold_XC(a3,:)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))

	    case (3)
           ui_XC=popold_XC+F_XC*(bm_XC-popold_XC+popold_XC(a1,:)-popold_XC(a2,:))

	    case (4)
           ui_XC=bm_XC+F_XC*(popold_XC(a1,:)-popold_XC(a2,:)+popold_XC(a3,:)-popold_XC(a4,:))

		case (5)
		   ui_XC=popold_XC(a5,:)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:)+popold_XC(a3,:) &
		         -popold_XC(a4,:))
        case (6) ! A linear crossover combination of bm_XC and popold_XC
           if (method(2) == 1) call random_number(F_CR)
		   ui_XC=popold_XC+F_CR*(bm_XC-popold_XC)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))

        end select
!!--------------------------------------------------------------------------!!
!!------Crossover operation-------------------------------------------------!!
        call random_number(rand_XC)
           mui_XC=0.0_RP
           mpo_XC=0.0_RP
        where (rand_XC < CR_XC)
           mui_XC=1.0_RP
!           mpo_XC=0.0_RP
        elsewhere
!           mui_XC=0.0_RP
           mpo_XC=1.0_RP
        end where

		ui_XC=popold_XC*mpo_XC+ui_XC*mui_XC
!!--------------------------------------------------------------------------!!
!!------Evaluate fitness functions and find the best member-----------------!!
        do i=1,NP
!!------Confine each of feasible individuals in the lower-upper bound-------!!
		   ui_XC(i,:)=max(min(ui_XC(i,:),XCmax),XCmin)
           call obj(ui_XC(i,:), tempval)
	       nfeval=nfeval+1
	       if (tempval < val(i)) then
	          pop_XC(i,:)=ui_XC(i,:)
	          val(i)=tempval
	          if (tempval < bestval) then
	             bestval=tempval
		         bestmem_XC=ui_XC(i,:)
              end if
           end if
        end do
        bestmemit_XC=bestmem_XC
	    if( (refresh > 0) .and. (mod(iter,refresh)==0)) then
   		     if (method(3)==1) write(unit=iwrite,FMT=203) iter
		     write(unit=*, FMT=203) iter
 		     do i=1,Dim_XC
		         if (method(3)==1) write(unit=iwrite, FMT=202) i, bestmem_XC(i)
				 write(*,FMT=202) i,bestmem_XC(i)
             end do
			 if (method(3)==1) write(unit=iwrite, FMT=201) bestval
			 write(unit=*, FMT=201) bestval
        end if
        iter=iter+1
        if ( bestval <= VTR .and. refresh > 0) then
           write(unit=iwrite, FMT=*) ' The best fitness is smaller than VTR'
		   write(unit=*, FMT=*) 'The best fitness is smaller than VTR'
           exit
        endif
	 end do
!!------end the evolutionary computation------------------------------!!
201 format(2x, 'bestval =', ES14.7, /)
202 format(5x, 'bestmem_XC(', I3, ') =', ES12.5)
203 format(2x, 'No. of iteration  =', I8)
end subroutine DE_Fortran90


function randperm(num)
  use data_type, only : IB, RP
  implicit none
  integer(kind=IB), intent(in) :: num
  integer(kind=IB) :: number, i, j, k
  integer(kind=IB), dimension(num) :: randperm
  real(kind=RP), dimension(num) :: rand2
  intrinsic random_number
  call random_number(rand2)
  do i=1,num
     number=1
     do j=1,num
        if (rand2(i) > rand2(j)) then
	       number=number+1
        end if
     end do
     do k=1,i-1
        if (rand2(i) <= rand2(k) .and. rand2(i) >= rand2(k)) then
	       number=number+1
        end if
     end do
     randperm(i)=number
  end do
  return
end function randperm
