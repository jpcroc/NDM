module Tpara
#ifdef PARA
  use mpi
  integer, parameter :: NDM_MPI_REAL_DOUBLE = MPI_REAL8
  integer::MPI_COMM_space
  integer:: grp_world
#endif
  integer :: myidsp,nprocspace,nprocs 			! numero de process mis là pour être utilisé en sequentiel
  integer,dimension(MPI_STATUS_SIZE):: status  ! statut de la communication
  integer::ierr
  type para_space_config
     integer, allocatable :: res_cpu(:,:)   	!stocke le nombre de cellules de chaques decoupages pour le meilleur decoupage
     integer, allocatable :: proc_cell(:)          !proc_cell(i) : Numero du proc associe a la cellule i
     integer, allocatable :: proc_voisin(:)        ! liste des processeurs voisins du processeur courant
     integer, allocatable :: cell_frontiere(:,:)   ! (i,j) jeme cellule frontiere associee au ieme processeur voisin
     integer, allocatable :: nbr_cell_frontiere(:) ! nbre de cellules frontieres associees au ieme processeur voisin
     integer :: nbr_cell_ftm                       ! nbr de cellules fantomes du processeur courant
     integer, allocatable :: cell_ftm(:)           ! liste des cellules fantomes du processeur courant
     integer :: nbr_proc_voisin            ! nbre de processeurs voisins du processeur courant
  end type para_space_config

  type(para_space_config)::psc0
  
end module Tpara
