subroutine read_ml_file
use T_kind_param_m, ONLY:  double
use gen_com_m, ONLY: lenfnam,fnam,angst,ev2erg,rangml

use  ml_in_ndm_module, ONLY: ml_type, one_pi, &
                             descriptor_type, target_type, force_comp, &
                             kernel_type,      &
                             nd_data, nd_fingerprint,n_frac,&
                             iread_ml, isave_ml, &
                             toy_model,krr_error,debug,rescale,s_max_r,s_min_r,s_max_i,s_min_i,&
                             search_hyp,&
                             n_g2_eta,n_g2_rs,n_g3_eta, n_g3_zeta,n_g3_lambda, &
                             n_rbf,n_cheb,l_max,j_max,r0,r_cut, &
                             kappa_acd,acd_threshold,  &
                             alpha_soap,n_soap,lsoap,lsoap_fcut, lsoap_diag,  &
                             sparsification,sparsification_by_acd,sparsification_by_entropy,sparsification_by_cur, max_data, &
                             alpha_acd,acd,acd_fcut,acd_weighted,temp_ini,ksi_ini,tau,mc_step, &
                             path,build_subdata,selection_type,pref,ns_data,i_begin,kelem,seed,write_desc,weighted,massat, &
                             marginal_likelihood, jj_max, db_file, db_path, strict_behler, descriptor_behler, &
                             descriptor_g2, descriptor_g3, descriptor_afs, descriptor_g2_afs, descriptor_behler, descriptor_soap, descriptor_pow_so3, &
                             descriptor_pow_so4, descriptor_bispectrum_so3, descriptor_bispectrum_so4, descriptor_g2_bispectrum_so4, descriptor_mtp, &
                             descriptor_magnetic_sld, descriptor_magnetic_sld_afs, &
                             optimize_weights, class_no_optimize_weights, &
                             optimize_weights_L1, optimize_weights_L2, optimize_weights_Le, &
                             optimize_ga_population, &
                             descriptor_g2_pow_so4, sign_stress, sign_stress_big_box, snap_fit_type, snap_class_constraints, &
                             fit_home_hb, fit_home_regularization, train_only, &
                             mtp_poly_min, mtp_poly_max, val_desc_max, lbso4_diag,  &
                             no_class_weights, max_iter_optimize_weights, factor_force_error, factor_energy_error, factor_stress_error, &
                             inv_r0, desc_forces, lbso3_diag, eta_max_g2, factor_weight_mass, lsoap_diag, lsoap_norm, lsoap_lnorm, nspecies_soap, r_cut_width_soap, &
                             lambda_krr, min_lambda_krr, max_lambda_krr, n_values_lambda_krr, vector_lambda_krr, ml_type_descriptors, &
                             magnetic_sld_j_dim, magnetic_sld_s2_dim, magnetic_sld_s4_dim, r_cut_magnetic, &
                             fix_no_of_elements, weight_per_element, chemical_elements, chemical_elements_invisible, snap_order
use k_cross_validation, ONLY :  kcross,n_kcross
use def_kernels, ONLY : length_kse, sigma_kse
use temporary_data_cov, ONLY : dim_data
implicit none

namelist /input_ml/ ml_type, &
                    descriptor_type, target_type, force_comp, &
                    nd_data,  dim_data, nd_fingerprint, &
                    kcross, n_kcross,n_frac,  &
                    iread_ml,isave_ml, kernel_type, &
                    marginal_likelihood, &
                    toy_model,krr_error,debug, &
                    search_hyp,&
                    length_kse,sigma_kse,rescale,s_max_r,s_min_r,s_max_i,s_min_i, &
                    n_g2_eta,n_g2_rs,n_g3_eta, n_g3_zeta,n_g3_lambda, &
                    n_rbf,n_cheb,l_max,j_max,r0,r_cut, &
                    kappa_acd,acd_threshold,   &
                    alpha_soap,n_soap,lsoap,lsoap_fcut, &
                    sparsification,sparsification_by_acd,sparsification_by_entropy,sparsification_by_cur, max_data, &
                    alpha_acd,acd,acd_fcut,acd_weighted,temp_ini,ksi_ini,tau,mc_step, &
                    path,build_subdata,selection_type,pref,ns_data,i_begin,kelem,seed, &
                    write_desc,weighted,massat, db_file, db_path, strict_behler, &
                    optimize_weights, class_no_optimize_weights, optimize_weights_L1, optimize_weights_L2, optimize_weights_Le, &
                    optimize_ga_population, &
                    factor_force_error, factor_energy_error, factor_stress_error, &
                    sign_stress, sign_stress_big_box, snap_fit_type, &
                    snap_class_constraints, train_only, mtp_poly_min, mtp_poly_max, val_desc_max, lbso4_diag, max_iter_optimize_weights, inv_r0, &
                    desc_forces, lbso3_diag, eta_max_g2, factor_weight_mass, nspecies_soap, r_cut_width_soap, lsoap_diag, lsoap_norm, lsoap_lnorm, &
                    lambda_krr, min_lambda_krr, max_lambda_krr, n_values_lambda_krr, &
                    magnetic_sld_j_dim, magnetic_sld_s2_dim, magnetic_sld_s4_dim, r_cut_magnetic, &
                    fix_no_of_elements, weight_per_element, chemical_elements, chemical_elements_invisible, snap_order

character(len=128) :: fnamtin
integer :: luml, none_class, i, nitems
!integer :: i

ml_type=1    ! 1 - KRR ; 0 SNAP_1 ; 2 - GP
iread_ml=0   ! 0 compute in the fly, 1 read from previous run
isave_ml=0   ! 0 do nothing, 1 write on the HDD and run, 2 MPI and threading ...
db_file="db_model.in"
db_path="DB/"
nd_data=0    ! the number of files in the repository
dim_data=100 ! the dimension of the database

desc_forces=.true.
lbso4_diag=.true.
strict_behler=.false.
 kcross=.false.
 marginal_likelihood=.false.
 sparsification=.false.
 sparsification_by_entropy=.false.
 sparsification_by_acd=.false.
 sparsification_by_cur=.false.
 n_frac=0d0
 n_kcross=0
 kernel_type=1
 krr_error=.false.
 toy_model=.true.
 search_hyp=.false.

 lambda_krr=-1.d0
 n_values_lambda_krr=21
 min_lambda_krr=1.d-10
 max_lambda_krr=1.d+10

 length_kse=1.d0
 sigma_kse=1.d0
 fnamtin = fnam(1:lenfnam)//'.ml'
 descriptor_type=1
 target_type=1
 force_comp=1
 snap_fit_type=fit_home_hb
 snap_class_constraints="02"
 sign_stress=1.d0
 sign_stress_big_box=1.d0
 train_only=.false.
 factor_weight_mass=50.d0
 snap_order=1
 optimize_ga_population=40
 optimize_weights=.false.
 optimize_weights_L2=.false.
 optimize_weights_L1=.false.
 optimize_weights_Le=.false.
 class_no_optimize_weights=" "
 max_iter_optimize_weights=40
 factor_force_error=1.d0
 factor_energy_error=1.d0
 factor_stress_error=1.d0


!g2
n_g2_eta=3
n_g2_rs=1
eta_max_g2=0.8d0

!g3
n_g3_eta=3
n_g3_zeta=2
n_g3_lambda=2

!mtp
mtp_poly_min=2
mtp_poly_max=5


!afs
 n_rbf=4          ! afs, pow_so3, bispectrum_so3
 n_cheb=5         ! afs

!magnetic
magnetic_sld_j_dim=10
magnetic_sld_s2_dim=10
magnetic_sld_s4_dim=10


l_max=4          ! pow_so3, bispectrum_so3, soap
lbso3_diag=.true. ! bispectrum_so3

j_max=1.5        ! pow_so4, bispectrum_so4
r0=20d0          ! pow_so4, bispectrum_so4
inv_r0=(1.d0 - 0.02d0/one_pi)
jj_max=int(2*j_max)
r_cut=5.d0          !active g2, g3, rbf, pow_so4, bispectrum_so4, soap

n_soap=2            ! radial part of soap
lsoap=.false.       ! soap
lsoap_fcut=.false.  ! acd
r_cut_width_soap=0.5d0
!debug pourposes
!alpha_soap=2d0      ! soap
!atom_sigma=0.5, alpha_soap = 0.5/atom_sigma**2
!atom_sigma_soap=0.5d0
!alpha_soap = 0.5d0/(atom_sigma_soap**2)
alpha_soap = 2.d0
lsoap_diag=.false.   !
lsoap_norm=.true.  !
lsoap_lnorm=.true.
nspecies_soap=1     ! 1 - is our style, 2 - is Gabor style

 kappa_acd=1          ! soap, acd
 acd=.false.      ! acd
 acd_threshold=1d0! acd
 alpha_acd=1d0    ! acd
 acd_fcut=.false. ! acd
 acd_weighted=.false. ! acd
 temp_ini=1d0
 ksi_ini=1d0
 tau=0.95d0
 mc_step=1000
 rescale=.false.
 s_max_r=1d1
 s_min_i=1d0
 s_max_r=1d1
 s_min_i=1d0
 r_cut_magnetic=5.3
 build_subdata=.false.
 selection_type=1 ! 1 - selects first "ns" elements of the database
                  ! 2 - selects last "ns" elements of the database
                  ! 3 - selects randomly "ns" subsets of "kelem" elements of the database
 pref="00"        ! prefixe of subdatabase name
 ns_data=100
 i_begin=0
 kelem=100
 seed=11          ! seed for random generator
 write_desc=.false. ! stores descriptor
 path="/home/wesley/database/"
val_desc_max = 1.d0



 weighted=.false. ! if true : weighted descriptors and number of type of descriptors n2_typ,n3_typ = 1
 massat=(/51d0,4d0/)  ! weight
fix_no_of_elements=1
chemical_elements="Fe"
chemical_elements_invisible=" "
weight_per_element="1.d0"


!debug  write(*,*) fnamtin
 luml=621
 open(unit=luml, file=fnamtin, status='unknown')
 read (luml, nml=input_ml)
 close (luml)

if (snap_fit_type == fit_home_regularization) call set_grid_lambda_krr

if ((snap_fit_type==fit_home_hb) .and.(optimize_weights)) then
   if (lambda_krr < 0.d0) then
      lambda_krr = (min_lambda_krr + max_lambda_krr)/2.d0
      if (rangml==0) write(6,*) 'ML: The initial lambda_krr is negative and set-up to the average value  between lambda_krr_min and lambda_krr_max', lambda_krr
   end if
   if ((min_lambda_krr<0.d0).or.(max_lambda_krr<0.d0)) then
     if (rangml==0) write(6,*) 'ml: the min_lambda_krr or max_lambda_krr is negative', min_lambda_krr, max_lambda_krr
     stop 'in read_ml_file one of the value of min_ or max_ lambda_krr is negative ...'
   end if

   if (optimize_weights_L1.and.optimize_weights_L2) then
     if (rangml==0) write(6,*) 'ML: optimize_weigths_L1 and optimize_weigths_L2 cannot be true simulanously ', optimize_weights_L1, optimize_weights_L2
     stop 'in read_ml_file both optimize_weigths_L1 and optimize_weigths_L2 are true'
   end if

   if (optimize_weights_L1.and.optimize_weights_Le) then
     if (rangml==0) write(6,*) 'ML: optimize_weigths_L1 and optimize_weigths_Le cannot be true simulanously ', optimize_weights_L1, optimize_weights_Le
     stop 'in read_ml_file both optimize_weigths_L1 and optimize_weigths_Le are true'
   end if
   if (optimize_weights_L2.and.optimize_weights_Le) then
     if (rangml==0) write(6,*) 'ML: optimize_weigths_L2 and optimize_weigths_Le cannot be true simulanously ', optimize_weights_L2, optimize_weights_Le
     stop 'in read_ml_file both optimize_weigths_L2 and optimize_weigths_Le are true'
   end if

   if (min_lambda_krr > max_lambda_krr) then

      if (rangml==0) write(6,*) 'ML: lambda_krr_min should be lower thzn  lambda_krr_max', min_lambda_krr, max_lambda_krr
      stop 'in read_ml_file the order of lambda_krr. Fatal error '

   end if
end if

call periodic_table
call fix_type_of_atoms





inv_r0=one_pi*inv_r0
jj_max = int(2*j_max)


if (.not.((snap_order==1).or.(snap_order==2))) then
  if (rangml==0) write(6,*) 'ML: not implemented this order of snap potetial snap_order is ', snap_order
  stop ' read_ml snap_order is out of range'
end if

if (.not.toy_model) call get_number_of_files(nd_data)


if ( .not.(   (snap_fit_type==0)  .or. &
              (snap_fit_type==1) .or. &
              (snap_fit_type==2) .or. &
              (snap_fit_type==10) .or. &
              (snap_fit_type==4) .or. &
              (snap_fit_type==3) )) then

    if (rangml==0) write(6,*) 'ML error: this type of fitting is not implemented'
    if (rangml==0) write(6,*) 'ML error: snap_fit_type should be 0, 1, 2, 3, 4 ot 10'
    stop 'in read_ml_file invalid snap_fit_type'
end if


 if (kcross.and.marginal_likelihood) then
   if (rangml==0) then
    write(6,*) 'Both kcross or marginal_likelihood cannot be true.'
    write(6,*) 'stop in read_ml_file'
   end if
   stop
 end if

if (ml_type /= ml_type_descriptors) then
    if (.not.(desc_forces)) then
     if (rangml==0) then
        write(6,'("ML: ---------------------------------WARNING------------------------------------")')
        write(6,'("ML: WARNING THE FORCES ARE DESACTIVATED. IF YOU WANT FORCES PLEASE PUT desc_forces=.true.")')
        write(6,'("ML: ---------------------------------WARNING------------------------------------")')
     end if
   end if
end if

 if (sparsification_by_entropy.or.sparsification_by_acd) sparsification=.true.

 if (strict_behler) then
   if (descriptor_type /= descriptor_behler ) then
      if (rangml==0) write(6,*)'If strict_behler is  T, the descriptor should be Behler i.e. descriptor_type=3 '
   end if
 end if
 if ((descriptor_type==descriptor_g2).and.((n_g2_eta<=0).or.(n_g2_rs<=0))) then
    if (rangml==0) write(6,*)"n_g2_eta and n_g2_rs should be larger than 0"
 elseif ((descriptor_type==descriptor_g3).and.((n_g3_eta<=0).or.(n_g3_zeta<=0).or.(n_g3_lambda<=0))) then
    if (rangml==0) write(6,*)"n_g3_eta,n_g3_zeta,n_g3_lambda should be larger than 0"
 elseif ((descriptor_type==descriptor_behler).and.((n_g2_eta<=0).or.(n_g2_rs<=0).or.(n_g3_zeta<=0).or.(n_g3_lambda<=0))) then
    if (rangml==0) write(6,*)"n_g2 g3_eta,n_g2_rs,n_g3_zeta,n_g3_lambda should be larger than 0"
 elseif ((descriptor_type==descriptor_afs).and.((n_rbf<=0).or.(n_cheb<0))) then
    if (rangml==0) then
       write(6,*)"n_rbf should be larger than 0"
       write(6,*)"n_cheb should be larger or at least equals 0"
    endif
 elseif ((descriptor_type==6).and.((n_rbf<=0).or.(l_max<0))) then
    if (rangml==0) then
       write(6,*)"n_rbf should be larger than 0"
       write(6,*)"l_max should be larger or at least equals 0"
    endif
 elseif ((descriptor_type==7).and.((n_rbf<=0).or.(l_max<0))) then
    if (rangml==0) then
       write(6,*)"n_rbf should be larger than 0"
       write(6,*)"l_max should be larger or at least equals 0"
    endif
 elseif (descriptor_type==descriptor_pow_so4) then
    if (rangml==0) then
         if (j_max<0) then
             write(6,'("ML: pow_so4 parametrization j_max should be larger or at least equals 0")')
             stop 'read_ml_file pow_so4 j_max negative'
         end if
         if (.not.(inv_r0 > 0.d0).or. (inv_r0 < 1.d0) ) then
            write(6,'("ML: pow_so4 parametrization, inv_r0 should be larger than 0 and lower than 1")')
             stop 'read_ml_file pwo_so4 inv_r0 beyound the limits'
         endif
    endif
 elseif ((descriptor_type==descriptor_bispectrum_so4)) then
    if (rangml==0) then
         if (j_max<0) then
             write(6,'("ML: bi_so4 parametrization j_max should be larger or at least equals 0")')
             stop 'read_ml_file bi_so4 j_max negative'
         end if
         if (.not.(inv_r0 > 0.d0).or. (inv_r0 < 1.d0) ) then
            write(6,'("ML: bi_so4 parametrization, inv_r0 should be larger than 0 and lower than 1")')
             stop 'read_ml_file bi_so4 inv_r0 beyound the limits'
         endif
    endif
 endif

 if (build_subdata) then
    if (selection_type==1) then
       if (rangml==0) write(6,*)"ns_data should be larger than 0 and lower or equal nd_data"
    elseif (selection_type==2) then
       if (rangml==0) write(6,*)"ns_data should be larger than 0 lower or equal nd_data"
    elseif (selection_type==3) then
       if (rangml==0) write(6,*)"ns_data,kelem should be larger than 0 lower or equal nd_data"
    else
       if (rangml==0) write(6,*)"selection_type should be 1,2 or 3"
       stop
    endif
 endif

if (.not.(desc_forces)) then
   if (rangml==0) write(6,'("ML:  WARNING THE DESCRIPTORS FOR FORCES ARE NOT COMPUTED i.e. desc_forces = .false. . You are adult and you do whatever you want ... even without THE FORCE")')
end if


if ((n_frac > 0.5) .or. (n_frac < 0.d0)) then

   if (rangml==0) write(6,*) 'ML: < read_ml_file > error n_frac should be in [0.0, 0.5]. Present value: ', n_frac
   stop

end if

if (train_only) then
    if (rangml==0) then
       write(6,*) 'ML: Only TRAINING. If you want to test also please put train_only = .false. '
    end if
end if

 if (kcross) then

   if (n_kcross==0) then
    if (rangml==0) then
     write(6,*) 'ML: < read_ml_file > ERROR :n_kcross should be bigger than 0. Between  5 and  10 is a good choice'
     write(6,*) 'ML: < read_ml_file > ERROT :but it can be any value. You figure out that the trainning time will be x n_kcross'
     stop
    end if
   end if

   if (nd_data==0) then
    if (rangml==0) then
     write(6,*) 'n_data - size if the database - should be bigger than 0'
     stop
    end if
   end if

end if    !kcross

if (rangml==0) then
        select case (ml_type)
           case (0)
                  write(6,*) 'ML: Machine Learning type with linear basis,  LML '
           case (1)
                  write(6,*) 'ML: Machine learning type with  Kernel Ridge Regression'
           case (2)
                  write(6,*) 'ML: Machine learning type with Gaussian Precesses  (not implemented, yet)'
                  stop
        end select
end if


if ((ml_type .eq. 1).or.(ml_type .eq. 2)) then
   select case (kernel_type)
     case (1)
       if (rangml==0) write(6,*) ' KERNEL SQUARED EXPONENTIAL'
     case (2)
       if (rangml==0) write(6,*) ' KERNEL ORSTEIN UHLENBECK'
     case (3)
       if (rangml==0) write(6,*) ' KERNEL MATERN CLASS'
     case (4)
       if (rangml==0) write(6,*) ' KERNEL SOAP'
     case default
       if (rangml==0) write(6,*) 'No implementation for kernel_type...', kernel_type
       stop
   end select
end if

if (toy_model) then
      if (rangml==0) write(6,*) ' TOY model - ML_perfect descriptors'
else
  select case (descriptor_type)
     case (descriptor_g2)
       if (rangml==0) write(6,*) ' NDM + ML_G2'
     case (descriptor_g3)
       if (rangml==0) write(6,*) ' NDM + ML_G3'
     case (descriptor_behler)
       if (rangml==0) write(6,*) ' NDM + ML_BEHLER'
     case (descriptor_afs)
       if (rangml==0) write(6,*) 'NDM + ML_AFS'
     case (descriptor_g2_afs)
       if (rangml==0) write(6,*) 'NDM + ML_G2_AFS'
     case (descriptor_soap)
       if (rangml==0) write(6,*) ' NDM + ML_SOAP'
     case (descriptor_pow_so3)
       if (rangml==0) write(6,*) ' NDM + ML_POW_SO3'
     case (descriptor_bispectrum_so3)
       if (rangml==0) write(6,*) ' NDM + ML_BISPECTRUM_SO3'
     case (descriptor_pow_so4)
       if (rangml==0) write(6,*) ' NDM + ML_POW_SO4'
     case (descriptor_g2_pow_so4)
       if (rangml==0) write(6,*) ' NDM + ML_G2_POW_SO4'
     case (descriptor_bispectrum_so4)
       if (rangml==0) write(6,*) ' NDM + ML_BISPECTRUM_SO4'
     case (descriptor_mtp)
       if (rangml==0) write(6,*) ' NDM + MTP'
     case (descriptor_g2_bispectrum_so4)
       if (rangml==0) write(6,*) ' NDM + ML_G2_BISPECTRUM_SO4'
     case (descriptor_magnetic_sld)
       if (rangml==0) write(6,*) ' NDM + ML_MAGNETIC_SLD'
     case (descriptor_magnetic_sld_afs)
       if (rangml==0) write(6,*) ' NDM + ML_MAGNETIC_SLD_AFS'
     case default
       if (rangml==0) write(6,*) 'No implementation for descriptor_type...', descriptor_type
       stop
  end select
end if

if (optimize_weights) then
    none_class=nitems(class_no_optimize_weights)
    if (rangml==0) then
      write(6,'("ML: Weights optimization will be performed")')
    end if
    if (allocated(no_class_weights)) deallocate(no_class_weights) ; allocate(no_class_weights(none_class))
    read(class_no_optimize_weights,*) no_class_weights

    if (debug) then
       if (rangml==0) then
         write(6,'("ML: in the following classes there no weights optimization ...")')
         do i=1,none_class
           write(6,*) i, 'class..', no_class_weights(i)
         end do
       end if
    end if
end if


return
end subroutine  read_ml_file



subroutine print_milady(rangloc)
implicit none
integer, intent(in) :: rangloc
character(len=1) :: quote,dquote


 quote=char(39)
dquote=char(34)

if (rangloc==0) then


!write(6,'("__/\\\\____________/\\\\________/\\\_____________________________/\\\\\\\\\\\\__________________         ")')
!write(6,'(" _\/\\\\\\________/\\\\\\_______\/\\\____________________________\/\\\////////\\\________________        ")')
!write(6,'("  _\/\\\//\\\____/\\\//\\\__/\\\_\/\\\____________________________\/\\\______\//\\\____/\\\__/\\\_       ")')
!write(6,'("   _\/\\\\///\\\/\\\/_\/\\\_\///__\/\\\______________/\\\\\\\\\____\/\\\_______\/\\\___\//\\\/\\\__      ")')
!write(6,'("    _\/\\\__\///\\\/___\/\\\__/\\\_\/\\\_____________\////////\\\___\/\\\_______\/\\\____\//\\\\\___     ")')
!write(6,'("     _\/\\\____\///_____\/\\\_\/\\\_\/\\\_______________/\\\\\\\\\\__\/\\\_______\/\\\_____\//\\\____    ")')
!write(6,'("      _\/\\\_____________\/\\\_\/\\\_\/\\\______________/\\\/////\\\__\/\\\_______/\\\___/\\_/\\\_____   ")')
!write(6,'("       _\/\\\_____________\/\\\_\/\\\_\/\\\\\\\\\\\\\\\_\//\\\\\\\\/\\_\/\\\\\\\\\\\\/___\//\\\\/______  ")')
!write(6,'("        _\///______________\///__\///__\///////////////___\////////\//__\////////////______\////________ ")')


write(6,'(",---.    ,---.-./`)   .---.       ____    ______        ____     __  ")')
write(6,'("|    \  /    \ .-.,",a,")  | ,_|     .",a,"  __ `.|    _ `",a,a,".    \   \   /  / ")')quote,quote,quote,quote
write(6,'("|  ,  \/  ,  / `-",a," \,-./  )    /   ",a,"  \  \ _ | ) _  \    \  _. /  ",a,"  ")')quote, quote, quote
write(6,'("|  |\_   /|  |`-",a,"`",a,"`\  ",a,"_ ",a,"`)  |___|  /  |( ",a,a,"_",a,"  ) |     _( )_ .",a,"   ")') quote, dquote, quote, quote, quote, quote, quote, quote
write(6,'("|  _( )_/ |  |.---.  > (_)  )     _.-`   | . (_) `. | ___(_ o _)",a,"    ")') quote
write(6,'("| (_ o _) |  ||   | (  .  .-",a,"  .",a,"   _    |(_    ._) ",a,"|   |(_,_)",a,"     ")')quote, quote,quote,quote
write(6,'("|  (_,_)  |  ||   |  `-",a,"`-",a,"|___|  _( )_  |  (_.\.",a," / |   `-",a,"  /      ")')quote, quote, quote, quote
write(6,'("|  |      |  ||   |   |        \ (_ o _) /       .",a,"   \      /       ")') quote
write(6,'("",a,"--",a,"      ",a,"--",a,a,"---",a,"   `--------`",a,".(_,_).",a,a,"-----",a,"`      `-..-",a,"        ")') quote, quote, quote, quote, quote, quote,&
quote,quote,quote,quote, quote
write(6,'("contributions:                                             ")')
write(6,'("A. M. Goryaeva, W. Unn-Toc, MCM                            ")')
write(6,'("copyright CEA ...                                          ")')
write(6,'("email:mihai-cosmin.marinica@cea.fr                         ")')


end if


end subroutine print_milady



subroutine set_grid_lambda_krr
use ml_in_ndm_module, only: min_lambda_krr, max_lambda_krr, vector_lambda_krr, n_values_lambda_krr, lambda_krr, rangml
implicit none
integer :: i
real(kind(0.d0)) :: delta_lambda

 if (n_values_lambda_krr .lt. 0) then
    if (rangml==0) then
        write(6,'("ML: fatal error, n_values_lambda_krr cannot be negative, currently is:", i6)') n_values_lambda_krr
        stop "in set_grid_lambda_krr wrong value for n_values_lambda_krr"
    end if
 end if


 if (n_values_lambda_krr == 1) then
    if ( (min_lambda_krr /= max_lambda_krr) ) then
      if (rangml==0) then
        write(6,'("ML: fatal error, when n_values_lambda_krr==1 min_lambda_krr and  max_lambda_krr should be equal :", 2e20.10)') &
          min_lambda_krr, max_lambda_krr
        stop "in set_grid_lambda_krr wrong value for n_values_lambda_krr and min_lambda_krr, max_lambda_krr"
      end if
    end if
 end if



 if ( (min_lambda_krr .lt. 0) .or. (max_lambda_krr .lt. 0 ) ) then
  n_values_lambda_krr =21
  if (allocated(vector_lambda_krr))  deallocate(vector_lambda_krr) ; allocate(vector_lambda_krr(n_values_lambda_krr))
  vector_lambda_krr=(/1.d-10, 1.d-09, 1.d-08, 1.d-07, &
                  1.d-06, 1.d-05, 1.d-04, 1.d-03, &
                  1.d-02, 1.d-01, 1.d+00, 1.d+01, &
                  1.d+02, 1.d+03, 1.d+04, 1.d+05, &
                  1.d+06, 1.d+07, 1.d+08, 1.d+09, &
                  1.d+10/)
  else if ( (min_lambda_krr == max_lambda_krr) ) then
    n_values_lambda_krr=1
    if (allocated(vector_lambda_krr))  deallocate(vector_lambda_krr) ; allocate(vector_lambda_krr(n_values_lambda_krr))
    vector_lambda_krr=(/ min_lambda_krr /)
    lambda_krr = min_lambda_krr
  else
    if (allocated(vector_lambda_krr))  deallocate(vector_lambda_krr) ; allocate(vector_lambda_krr(n_values_lambda_krr))
    delta_lambda= (log10(max_lambda_krr) - log10(min_lambda_krr)) / dble(n_values_lambda_krr-1)
    do i=1,n_values_lambda_krr
     vector_lambda_krr(i) = 10.d0**(log10(min_lambda_krr) + dble(i-1) * delta_lambda)
    end do
 end if

return
end subroutine set_grid_lambda_krr



subroutine fix_type_of_atoms
  use ml_in_ndm_module, only: rangml, chemical_elements, weight_per_element, fix_no_of_elements, periodic_table_element, size_periodic_table, &
                              fix_ch_elements, fix_mass_elements, fix_type_to_periodic, fix_Z_elements, fix_weighted_for_element, &
                              lfix_weight_auto, &
                              fix_no_of_elements_invisible, fix_ch_elements_invisible, fix_mass_elements_invisible, fix_type_to_periodic_invisible, fix_Z_elements_invisible, fix_weighted_for_element_invisble, &
                              weighted, chemical_elements_invisible, linvisible
  use math
  implicit none
  integer :: int_local, nitems2, i_p, icount
  character(len=1) :: quote,dquote
  !character (len=60) :: READCH

  quote=char(39)
  dquote=char(34)


  !write(*,*) 'int', int_local
  !write(*,*) 'int', chemical_elements


  int_local=nitems2(chemical_elements)
  if (int_local==0) then
    if (rangml==0) write(6,'("ML error: The list of chemical elemnts is empty. Please put the appropiate chemical_elements list")')
    if (rangml==0) write(6,'("ML error: An example chemical_elements=",a,"Fe W Os",a)') dquote, dquote
    stop "read_ml_file with chemical_elements"
  end if


  if (int_local/=fix_no_of_elements) then
    if (rangml==0) write(6,'("ML error: The list of chemical elemnts and the list of weigths for each elemets has not the same lengths")')
    if (rangml==0) write(6,'("ML error: in that case we expect something with ",i6," values")') fix_no_of_elements
    stop "read_ml_file with weight_per_elements"
  end if



  fix_no_of_elements=int_local
  if (allocated(fix_ch_elements)) deallocate(fix_ch_elements) ; allocate(fix_ch_elements(int_local))
  if (allocated(fix_Z_elements)) deallocate(fix_Z_elements)  ; allocate(fix_Z_elements(int_local))
  if (allocated(fix_mass_elements)) deallocate(fix_mass_elements)  ; allocate(fix_mass_elements(int_local))
  if (allocated(fix_type_to_periodic)) deallocate(fix_type_to_periodic)  ; allocate(fix_type_to_periodic(int_local))
  read(chemical_elements, *) fix_ch_elements(:)

  icount = 0
  do int_local=1,size(fix_ch_elements, dim=1)
    do i_p=1, size_periodic_table
        if (periodic_table_element(i_p)%symbol == fix_ch_elements(int_local)) then
                    icount = icount + 1
                    fix_Z_elements(int_local)=periodic_table_element(i_p)%Z
                    fix_mass_elements(int_local)=periodic_table_element(i_p)%mass
                    fix_type_to_periodic(int_local) = i_p
       end if
    end do
    if (rangml==0) write(6,'("-----------------the list of possible atoms in the database-----------------------")')
    if (rangml==0) write(6,'(" id element, element, Z, mass", 2i6, " ", a2, " ",  f7.1, f15.5)') int_local, fix_type_to_periodic(int_local), &
                   fix_ch_elements(int_local), fix_Z_elements(int_local), fix_mass_elements(int_local)
    if (rangml==0) write(6,'("----------------------------------------------------------------------------------")')
  end do

if (weighted) then
  int_local=nitems2(weight_per_element)
  lfix_weight_auto=.false.
  if (int_local==0) then
    lfix_weight_auto=.true.
  else
    if (int_local/=fix_no_of_elements) then
       if (rangml==0) write(6,'("ML error: The list of chemical elemnts and the list of weigths for each elemets has not the same lengths")')
       if (rangml==0) write(6,'("ML error: in that case we expect something with ",i6," values")') fix_no_of_elements
       stop "read_ml_file with weight_per_elements"
    end if
  end if
  if (allocated(fix_weighted_for_element)) deallocate(fix_weighted_for_element) ; allocate(fix_weighted_for_element(int_local))
  read(weight_per_element, *) fix_weighted_for_element(:)
end if

! ....... invisible chemical element set up ... I do it what
  int_local=nitems2(chemical_elements_invisible)
  linvisible=.false.
  if (int_local>0) linvisible=.true.
    if (linvisible) then
      fix_no_of_elements_invisible=int_local

      if (allocated(fix_ch_elements_invisible)) deallocate(fix_ch_elements_invisible) ; allocate(fix_ch_elements_invisible(int_local))
      if (allocated(fix_Z_elements_invisible)) deallocate(fix_Z_elements_invisible)  ; allocate(fix_Z_elements_invisible(int_local))
      if (allocated(fix_mass_elements_invisible)) deallocate(fix_mass_elements_invisible)  ; allocate(fix_mass_elements_invisible(int_local))
      if (allocated(fix_type_to_periodic_invisible)) deallocate(fix_type_to_periodic_invisible)  ; allocate(fix_type_to_periodic_invisible(int_local))


      icount = 0
      read(chemical_elements_invisible, *) fix_ch_elements_invisible(:)
      do int_local=1,size(fix_ch_elements_invisible, dim=1)
        do i_p=1, size_periodic_table
          if (periodic_table_element(i_p)%symbol == fix_ch_elements_invisible(int_local)) then
                    icount = icount + 1
                    fix_Z_elements_invisible(int_local)=periodic_table_element(i_p)%Z
                    fix_mass_elements_invisible(int_local)=periodic_table_element(i_p)%mass
                    fix_type_to_periodic_invisible(int_local) = i_p
          end if
         end do

      if (rangml==0) write(6,'("-----------------the list of invisible atoms -----------------------")')
      if (rangml==0) write(6,'(" id element, element, Z, mass", 2i6, " ", a2, " ",  f7.1, f15.5)') int_local, fix_type_to_periodic_invisible(int_local), &
                   fix_ch_elements_invisible(int_local), fix_Z_elements_invisible(int_local), fix_mass_elements_invisible(int_local)
      if (rangml==0) write(6,'("----------------------------------------------------------------------------------")')
      end do
    end if
!...... end of invisible elements



  return
end subroutine fix_type_of_atoms
