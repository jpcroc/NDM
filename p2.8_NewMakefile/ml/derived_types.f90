module derived_types
use T_kind_param_m, ONLY:  double


type descriptor_cell

!dimension of the descriptor
integer :: dim_desc

!desc for energy for each atom
!energy(1:dim_desc,1:imm)
!energy(d,ia) is the d^th dimension of the descriptor on the i^th atom
real(double), dimension(:,:), allocatable :: energy

!force(1:dim_desc,1:imm,0:imm_neigh,1:3)
!force(d,ia,ja,1) the d component of force descriptor for ia^th atom derived
!on (ja ,x). If ja=0 then the derivative is with the atom itself.
real(double), dimension(:,:,:,:), allocatable :: force
!pack energy - the descritor for the energy of some configuration (1:dim_desc)
real(double), dimension(:), allocatable ::  pack_energy
!pack force - the descritor for force of some atom (1:dim_desc, 3, imm)
real(double), dimension(:,:,:), allocatable ::  pack_force
!pack stress from the derivatives ...(1:dim_desc,1:6)
real(double), dimension(:,:), allocatable ::  pack_stress

!n_neigh(ia) number of neighbours of the ia^th atom
integer, dimension(:), allocatable:: n_neigh
!number of ghost atoms for forces ...
integer, dimension(:), allocatable:: n_neigh_ghost
! type (species) of neighbours (ia, ith neighbours) and kind (if it is a replica for the atom 1, 2, 3  ...)
integer, dimension(:,:), allocatable :: type_neigh, kind_neigh
! the_same thing but gfor ghost atoms ... type (species) of neighbours (ia, ith neighbours) and kind (if it is a replica for the atom 1, 2, 3  ...)
integer, dimension(:,:), allocatable :: kind_neigh_ghost, kind_neigh_proc
end type descriptor_cell

type system_state
character(len=80) :: filename
!a system state
!in which class / klm belongs and the i_order of the file in orginal database
character(len=2) :: class
character(len=3) :: klm
character(len=6) :: cnumber
!from which line from database the file configuration comes ...
integer :: db_line
!the number of the configuration in that line
integer :: no_file_in_db_line

integer :: i_weights
! if this config is used for train or test
! train=T for training
! train=F for testing
logical :: train=.false.
logical :: selected=.false.
!the forces are reliable T/F
logical :: has_force=.false.
!the energy is reliable  T/F
logical :: has_energy=.false.
!the stress is reliable T/F
logical ::  has_stress=.false.
!the spins are reliable T/F
logical :: has_atomic_spin=.false.
!copy of weights from the databse line.
real(double) :: w_e, w_f, w_s, w_e_end, w_f_end, w_s_end
! if the box is small or big
logical :: small=.false.
!>>>begin small_box all these quantities are active for small boxes
! these values will be active in the case small=.true.
integer :: nxCell, nyCell, nzCell
! number of neighbours for each atom
integer, dimension(:),  allocatable :: n_neigh
! on which processor the atom is located ...
integer, dimension(:),  allocatable :: proc_atom
! type (species) of neighbours (ia, ith neighbours) and kind (if it is a replica for the atom 1, 2, 3  ...)
integer, dimension(:,:), allocatable :: type_neigh, kind_neigh
! distance between atoms (ia, ith neighbours)
real(kind=kind(1.d0)), dimension(:,:), allocatable :: r_ij
! distane projected on each direction  u(ia, ith neighbours, 1:3)= r(ia,1:3) - r(ith niegh, 1:3)
real(kind=kind(1.d0)), dimension(:,:,:), allocatable :: u_ij, u_per
!>>>end   small_box
character(len=3) :: state
!number of atoms
integer :: nat
!number of types of atoms
integer :: ntypes
! volume of the box
real(kind=kind(1.d0)) :: volume
!types of atoms
integer, dimension(:), allocatable :: itype
! link the type of the atom with the Z in the periodic table
real(double), dimension(:), allocatable :: fix_type_poscar_to_periodic
! direct cell(1:3, component_j), reciprocal bg_cell(1:3,component_j)
real(double), dimension(3,3)   :: cell, bg_cell
real(double), dimension(:,:), allocatable :: pos_cart,pos_crst
real(double), dimension(:,:), allocatable :: force
real(double), dimension(:,:), allocatable :: atomic_spin
real(double), dimension(3) :: energy
real(double), dimension(6)   :: stress
! this will be allocated of dimension(ntypes) and contains the atomic mass per type
real(double), dimension(:), allocatable :: mass_per_type
! Z per type
real(double), dimension(:), allocatable :: Z_per_type, weight_per_type
logical, dimension(:), allocatable :: invisible_per_type
integer :: spin
end type system_state

type db_line
!an object that englobe one line in the db
character(len=2) :: class
character(len=3) :: klm
integer :: no_total
integer :: no_selec
integer :: no_start
! The T/F on db line that correspond to the energy/force/stress optimization
logical :: has_db_force, has_db_energy, has_db_stress
real(double) :: w_e, w_f, w_s, w_e_end, w_f_end, w_s_end
! The T/F on db line that correspond to the energy/force/stress weights optimization
! not necessarly the same ... e.g. if all configurations on the line have  the flag F there is no optimization for the line
logical :: has_optimize_weights_energy, has_optimize_weights_force, has_optimize_weights_stress
end type db_line




type(descriptor_cell), dimension(:), allocatable :: config_desc, cofig_desc_copy
type(system_state)   , dimension(:), allocatable :: config_real, config_real_copy
type(db_line), dimension(:), allocatable :: db_model

end module derived_types
