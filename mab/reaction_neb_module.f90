module reaction_neb_module
USE T_kind_param_m, ONLY:  double
USE gen_com_m, ONLY:imm
implicit none 
integer :: npoints
real(double) :: lambda_min,lambda_max, delta_lambda,delta_neb
real(double), dimension(:), allocatable :: llambda
real(double), dimension(:), allocatable :: coord_neb
real(double), dimension(:), allocatable :: force_defect,free_energy_force_defect


 type discrete_path
    integer :: imm
    real(double), dimension(3,3)   :: at       ! the box in the at(xyz,123) format
    real(double), dimension(3,3)   :: bg      ! the inverse of the box in the at(xyz,123) format
    real(double), dimension(:,:),pointer :: xp_neb !  the  positions in neb trajectories
    real(double), dimension(:,:),pointer :: dxp_neb !  the derivatives positions in neb trajectories

    real(double), dimension(:,:),allocatable :: fp_neb

    integer,dimension(:), pointer            :: ityp             !  type of each the atoms
    integer                                  :: ntyp             !  how many types
    integer, dimension(:),allocatable        :: natoms_type      !  how many atoms for each type  
    integer, dimension(:),allocatable        :: ielat,iwmax      !  profile, this information is used 
 end type discrete_path


 type coeff_spline
    real(double),dimension(:), pointer :: b,c,d
 end type coeff_spline


 type(discrete_path), dimension(:), allocatable :: neb_images,lambda_images
 type(coeff_spline),dimension(:,:),allocatable :: atoms_on_spline



end module reaction_neb_module

