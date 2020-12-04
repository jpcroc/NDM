#ifdef PARA
module init_vois_mod
  use cellconfig,only:cell_config
  implicit none

  
        contains
subroutine init_voisinage (cellv)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  use gen_com_m,only:
  use tab_imm_m
!  use mod_para,only:MPI_COMM_space,
!  USE mpi
  use mod_para,only:MPI_COMM_space, nprocspace,myidsp,NDM_MPI_REAl_DOUBLE,proc_voisin,nbr_proc_voisin,nbr_cell_ftm,&
       &NBR_CELL_FRONTIERE,RES_CPU,CELL_FRONTIERE,cell_ftm

  implicit none
  type(cell_config),intent(in)::cellv
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: cell
  integer :: icell
  integer :: cell_vois
  integer :: est_present
  integer :: nb_cell_frontieres
  integer :: nb_frontieres
  integer :: nb_internes
  integer :: nb_fantomes_max
  integer :: i
  integer :: i_cell_ftm
  integer :: num_proc_vois

  ! intialisations preliminaires
  allocate(proc_voisin(min(nprocspace,26)))
  proc_voisin(:)=-1
  nbr_proc_voisin = 0
  nbr_cell_ftm  = 0
  allocate(nbr_cell_frontiere(min(nprocspace,26)))
  nbr_cell_frontiere(:) = 0
  ! Calcul du nombre de cellules frontieres
  nb_internes = max(res_cpu(myidsp,1)-2,0) * max(res_cpu(myidsp,2)-2,0) * max(res_cpu(myidsp,3)-2,0)
  nb_frontieres = res_cpu(myidsp,1)*res_cpu(myidsp,2)*res_cpu(myidsp,3) - nb_internes
  allocate(cell_frontiere(size(nbr_cell_frontiere,1),nb_frontieres))
  cell_frontiere(:,:) = 0
  nb_fantomes_max = (res_cpu(myidsp,1) + 2) * (res_cpu(myidsp,2) + 2) * (res_cpu(myidsp,3) + 2) 
  nb_fantomes_max = nb_fantomes_max - res_cpu(myidsp,1)*res_cpu(myidsp,2)*res_cpu(myidsp,3)
  allocate(cell_ftm(nb_fantomes_max))

  ! Calcul du nombre de cellules fantomes 

  ! boucle sur toutes les cellules
  do cell=1,cellv%noxyz

     if ( cellv%proc_cell(cell)==myidsp ) then
        ! si la cellule est locale

        do icell=1,26
           ! boucle sur les cellules voisines

           cell_vois = cellv%ncel(cell,icell)
           if (cell_vois/=0) then
              ! si la cellule voisine existe

              if ( cellv%proc_cell(cell_vois)/=myidsp ) then
                 ! Si la cellule voisine n'est pas locale

                 ! on stocke le processeur voisin si il n'est
                 ! pas deja connu
                 est_present=0
                 do i=1,nbr_proc_voisin
                    if ( proc_voisin(i)==cellv%proc_cell(cell_vois) ) then
                       est_present=1
                       num_proc_vois=i
                    endif
                 enddo
                 if (est_present==0) then
                    nbr_proc_voisin = nbr_proc_voisin + 1
                    proc_voisin(nbr_proc_voisin) = cellv%proc_cell(cell_vois)
                    num_proc_vois = nbr_proc_voisin
                 endif

                 ! on stocke la cellule courante comme cellule
                 ! a emettre si elle n'est pas deja connue
                 est_present=0
                 do i=1,nbr_cell_frontiere(num_proc_vois)
                    if (cell_frontiere(num_proc_vois,i)==cell) then
                       est_present=1
                    endif
                 enddo
                 if (est_present==0) then 
                    nbr_cell_frontiere(num_proc_vois) = nbr_cell_frontiere(num_proc_vois) + 1
                    cell_frontiere(num_proc_vois,nbr_cell_frontiere(num_proc_vois)) = cell
                 endif

              endif ! la cellule voisine n'est pas locale
           endif ! la cellule voisine existe
        enddo ! boucle sur les cellules voisines

     else  ! la cellule est non locale
        ! On regarde si elle est fantome

        do icell=1,26
           ! boucle sur les cellules voisines

           cell_vois = cellv%ncel(cell,icell)
           if (cell_vois/=0) then
              if ( cellv%proc_cell(cell_vois)==myidsp ) then

                 ! il s'agit bien d'une cellule fantome car une cellule voisine est
                 ! locale
                 ! on regarde si elle n'est pas deja presente et on l'ajoute
                 est_present=0
                 do i_cell_ftm=1,nbr_cell_ftm
                    if (cell_ftm(i_cell_ftm)==cell) est_present=1 
                 enddo
                 if (est_present==0) then
                    nbr_cell_ftm = nbr_cell_ftm + 1
                    cell_ftm(nbr_cell_ftm) = cell
                 endif

              endif ! la cellule voisine est locale

           endif ! si la cellule voisine existe

        enddo ! boucle sur les cellules voisines

     endif ! la cellule est locale

  enddo

end subroutine init_voisinage
end module init_vois_mod
#endif
