&input_ml
iread_energy = 1
!--------------------------------------------------------------------\
!                           ML Mode                                  |
!--------------------------------------------------------------------/

ml_type=0 ! 0            ! -1 - just compute descriptors ...
                     !  0 - functions SNAP etc,
                     !  1 - Kernel RR
                     !  2 - GAP, not yet implemented
                     !  3 - n - linear
                     ! -2 - analyse the data or dump kernel.
!-------------------------------------------------:

!--------------------------------------------------------------------\
!                           ML Model                                 |
!--------------------------------------------------------------------/

!---------------ML model--------------------------!
snap_order= 2 !1             ! 1 - linear , 2 quadratic , 11 - n-linear, 3 PolyC, 7 - kernel, 11 - nlinear
snap_type_quadratic= 1 !1    ! 1 default +LML precondition, 2 - full quadratic, 3 - bi-linear
!this is for n-linear
order_nlinear = 2 ! is the order of n-linear form. Active if snap_order=11. For the moment is implemented only the case 2.
!this is for PolyC
polyc_n_poly = 3    ! max Poly: 1(lieanr) to 3
polyc_n_hermite = 2 ! Max Hermite degree: 1(identity) to 4

!--------------------------------------------------------------------\
!                           Database parameters                      |
!--------------------------------------------------------------------/
!build_subdata=.false.
selection_type=3 ! 1 - selects first    "ns" elements of the database
                 ! 2 - selects last     "ns" elements of the database
                 ! 3 - selects randomly "ns" subsets of "kelem" elements of the database
seed=42          ! seed for random generator
db_file="db_model.in" !"db_model.in"
db_path="/ccc/scratch/cont002/dam/hellierd/DB/" !"DB/"

!-------------------------------------------------
!If you want to read the desgn matrix.
write_design_matrix = .false.
!-------------------------------------------------


!--------------------------------------------------------------------\
!                           ML Train                                 |
!--------------------------------------------------------------------/

! if only training:
train_only=.false.

!--------------train options----------------------!
!snap_fit_type
!0  - home made BH ;
!10 - home made with regularization ;
!1 - lapack QR with norm minimization ;
!2 - lapack + constraints QR based snap-class_constraint should be written  ;
!3 - lapack full ortho decomposition with rank estimation. Norm minimization.
!4 - lapack full SVD  with rank estimation. Norm minimization.
snap_fit_type=4 !4
lambda_krr=1.d-7

optimize_weights_db=.false.

!--------------------------------------------------------------------\
!                        Descriptors parameters                      |
!--------------------------------------------------------------------/

weighted=.false. ! if true : weighted descriptors and number of type of descriptors n2_typ,n3_typ = 1
weighted_3ch=.false.
weighted_auto=.false.
fix_no_of_elements=1 !1
chemical_elements=" Ge "

!--> cut-off of the descriptor ...
r_cut=4.9 !4.7d0

activate_k2b =.true. !.true.
sigma_2b=0.9d0 !0.3d0
delta_2b=5.5d0 !1.d0
np_radial_2b=50 !60
r_cut_in = 1.5d0 ! 1.5d0          !  where the many-body interaction is shut dow (in A)n, and a smooth link with ZBL is prepared. Below this distance, there is no fitting.
r_cut_width_in = 0.3d0 ! 0.3d0    ! controls the smoothness of the inner cutoff
type_fcut = 3 !2 or 3             ! the type of fcut function (there are several options). 

!---> The descriptord type
descriptor_type=9 !9   

!---> set up SO4 and bi-SO4
j_max=4.0 !4          ! pow_so4, bispectrum_so4
lbso4_diag=.true.  !lbso4_diag=.true.  ... only diagonal (Gabor like)
                    !lbso4_diag=.false. ... only diagonal (non diagonal -SNAP like)

zbl_potential = .true. !.true.    ! activate or not
r1_zbl = 1.0d0 !1.0d0             ! for r < r1_zbl, ONLY ZBL is used in the radial part of your ML
r2_zbl = 2.2d0 !2.2d0             ! internal parameter that controls the smoothness of ZBL

&end
