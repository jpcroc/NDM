module derived_types_ph
USE T_kind_param_m, ONLY:  double



type system_state
character(len=80) :: filename !a system state

character(len=2) :: class !in which class / klm belongs and the i_order of the file in orginal database
character(len=3) :: klm
character(len=6) :: cnumber

integer :: db_line !from which line from database the file configuration comes ...

integer :: no_file_in_db_line !the number of the configuration in that line

integer :: i_weights

logical :: train=.false. ! if this config is USEd for train or test ! train=T for training ! train=F for testing

logical :: has_force=.false. !the forces are reliable T/F

logical :: has_energy=.false. !the energy is reliable  T/F

logical ::  has_stress=.false. !the stress is reliable T/F

logical :: has_atomic_spin=.false. !the spins are reliable T/F

real(double) :: w_e, w_f, w_s, w_e_end, w_f_end, w_s_end !copy of weights from the databse line.

logical :: small=.false. ! if the box is small or big

integer :: nxCell, nyCell, nzCell !>>>begin small_box all these quantities are active for small boxes ! these values will be active in the case small=.true.

integer, dimension(:),  allocatable :: n_neigh ! number of neighbours for each atom

integer, dimension(:),  allocatable :: proc_atom ! on which processor the atom is located ...

integer, dimension(:,:), allocatable :: type_neigh, kind_neigh, kind_neigh_big ! type (species) of neighbours (ia, ith neighbours) and kind (if it is a replica for the atom 1, 2, 3  ...)

real(kind=kind(1.d0)), dimension(:,:), allocatable :: r_ij ! distance between atoms (ia, ith neighbours)

real(kind=kind(1.d0)), dimension(:,:,:), allocatable :: u_ij, uperiod_ij ! distane projected on each direction  u(ia, ith neighbours, 1:3)= r(ia,1:3) - r(ith niegh, 1:3)
!>>>end   small_box
character(len=3) :: state

integer :: nat !number of atoms

integer :: ntypes !number of types of atoms

integer :: im, imm !for the NDM interface purposes

real(kind=kind(1.d0)) :: volume ! volume of the box

integer, dimension(:), allocatable :: itype !types of atoms

real(double), dimension(3,3)   :: cell, bg_cell ! direct cell(1:3, component_j), reciprocal bg_cell(1:3,component_j)
real(double), dimension(:,:), allocatable :: pos_cart,pos_crst
real(double), dimension(:,:), allocatable :: force
real(double), dimension(:,:), allocatable :: atomic_spin
real(double), dimension(3) :: energy
real(double), dimension(6)   :: stress
! this will be allocated of dimension(ntypes) and contains the atomic mass per type
real(double), dimension(:), allocatable :: mass_per_type
integer :: spin
real(double), dimension(:,:), allocatable :: Rperiodic
integer, dimension(:), allocatable :: ia_ini
end type system_state

type atom_type_ph
  integer :: nn1, nn2 ,nn3, ia
end type atom_type_ph


type(system_state)   , dimension(:), allocatable :: config_real
type(atom_type_ph), dimension(:), allocatable    :: atom_ph
integer, dimension(:,:,:,:), allocatable :: inv_atom_ph

end module derived_types_ph
