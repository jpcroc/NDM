module Tpara
#ifdef PARA
  use mpi
  integer, parameter :: NDM_MPI_REAL_DOUBLE = MPI_REAL8
  integer, parameter :: NDM_MPI_COMPLEX_DOUBLE = MPI_COMPLEX16
  integer::MPI_COMM_space
  integer:: grp_world
  integer,dimension(MPI_STATUS_SIZE):: status
#else
   integer:: status 
#endif
  integer :: myidsp,nprocspace,nprocs 			! numero de process mis là pour être utilisé en sequentiel

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
     integer :: cell_debx, cell_deby, cell_debz     !numero de la premiere cellule locale suivant x, y et z
     integer :: cell_finx, cell_finy, cell_finz     !numero de la derniere cellule locale  suivant x, y et z
     integer :: nb_cell_x, nb_cell_y, nb_cell_z     !nb de cel locales suivant x y z
  end type para_space_config



end module Tpara
