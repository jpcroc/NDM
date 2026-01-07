#ifdef PARA
module init_vois_mod
  use cellconfig,only:cell_config
   USE arret_ndm_mod,only:arret_ndm

  use Tpara,only:para_space_config
  implicit none

  
        contains
subroutine init_voisinage (cellv,psc,lwrite)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  use gen_com_m,only:rang
  use Tpara,only: nprocspace,myidsp,NDM_MPI_REAl_DOUBLE
!  use mod_para,only:nbr_cell_ftm,NBR_CELL_FRONTIERE,RES_CPU,CELL_FRONTIERE,cell_ftm

  implicit none
  type(cell_config),intent(in)::cellv
  type(para_space_config)::psc ! para_space_config
  logical,optional::lwrite
  logical::lwrt=.false.
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: cell
  integer :: icell
  integer :: cell_vois
  integer :: est_present

  integer :: nb_frontieres
  integer :: nb_internes
  integer :: nb_fantomes_max
  integer :: i
  integer :: i_cell_ftm
  integer :: num_proc_vois
  if (present(lwrite))lwrt=lwrite
  ! intialisations preliminaires
  if( allocated(psc%proc_voisin)) deallocate(psc%proc_voisin)
  allocate(psc%proc_voisin(min(nprocspace,cellv%ncelvmax)))
!  call cellv%print
  psc%proc_voisin(:)=-1
  psc%nbr_proc_voisin = 0
  psc%nbr_cell_ftm  = 0

  if (allocated(psc%nbr_cell_frontiere)) deallocate(psc%nbr_cell_frontiere)
  allocate(psc%nbr_cell_frontiere(min(nprocspace,cellv%ncelvmax)))

  psc%nbr_cell_frontiere(:) = 0
  ! Calcul du nombre de cellules frontieres
  nb_internes = max(psc%res_cpu(myidsp,1)-2,0) * max(psc%res_cpu(myidsp,2)-2,0) * max(psc%res_cpu(myidsp,3)-2,0)
  nb_frontieres = psc%res_cpu(myidsp,1)*psc%res_cpu(myidsp,2)*psc%res_cpu(myidsp,3) - nb_internes

  if (allocated(psc%cell_frontiere )) deallocate(psc%cell_frontiere )
  allocate(psc%cell_frontiere(size(psc%nbr_cell_frontiere,1),nb_frontieres))

  psc%cell_frontiere(:,:) = 0
  nb_fantomes_max = (psc%res_cpu(myidsp,1) + 2) * (psc%res_cpu(myidsp,2) + 2) * (psc%res_cpu(myidsp,3) + 2) 
  nb_fantomes_max = nb_fantomes_max - psc%res_cpu(myidsp,1)*psc%res_cpu(myidsp,2)*psc%res_cpu(myidsp,3)

  if( allocated(psc%cell_ftm )) deallocate(psc%cell_ftm )
  allocate(psc%cell_ftm(nb_fantomes_max))

  ! Calcul du nombre de cellules fantomes 

  ! boucle sur toutes les cellules

  do cell=1,cellv%noxyz

     if ( cellv%proc_cell(cell)==myidsp ) then
        ! si la cellule est locale

        do icell=1,cellv%ncelvois(cell)
           ! boucle sur les cellules voisines

           cell_vois = cellv%ncel(cell,icell)
           if (cell_vois/=0) then
              ! si la cellule voisine existe

              if ( cellv%proc_cell(cell_vois)/=myidsp ) then
                 ! Si la cellule voisine n'est pas locale

                 ! on stocke le processeur voisin si il n'est
                 ! pas deja connu
                 est_present=0
                 do i=1,psc%nbr_proc_voisin
                    if ( psc%proc_voisin(i)==cellv%proc_cell(cell_vois) ) then
                       est_present=1
                       num_proc_vois=i
                    endif
                 enddo
                 if (est_present==0) then
                    psc%nbr_proc_voisin = psc%nbr_proc_voisin + 1
                    psc%proc_voisin(psc%nbr_proc_voisin) = cellv%proc_cell(cell_vois)
                    num_proc_vois = psc%nbr_proc_voisin
                 endif

                 ! on stocke la cellule courante comme cellule
                 ! a emettre si elle n'est pas deja connue
                 est_present=0
                 do i=1,psc%nbr_cell_frontiere(num_proc_vois)
                    if (psc%cell_frontiere(num_proc_vois,i)==cell) then
                       est_present=1
                    endif
                 enddo
                 if (est_present==0) then 
                    psc%nbr_cell_frontiere(num_proc_vois) = psc%nbr_cell_frontiere(num_proc_vois) + 1
                    psc%cell_frontiere(num_proc_vois,psc%nbr_cell_frontiere(num_proc_vois)) = cell
                 endif

              endif ! la cellule voisine n'est pas locale
           endif ! la cellule voisine existe
        enddo ! boucle sur les cellules voisines

     else  ! la cellule est non locale
        ! On regarde si elle est fantome
        do icell=1,cellv%ncelvois(cell)
           ! boucle sur les cellules voisines

           cell_vois = cellv%ncel(cell,icell)
           if (cell_vois/=0) then
              if ( cellv%proc_cell(cell_vois)==myidsp ) then

                 ! il s'agit bien d'une cellule fantome car une cellule voisine est
                 ! locale
                 ! on regarde si elle n'est pas deja presente et on l'ajoute
                 est_present=0
                 do i_cell_ftm=1,psc%nbr_cell_ftm
                    if (psc%cell_ftm(i_cell_ftm)==cell) est_present=1 
                 enddo
                 if (est_present==0) then
                    psc%nbr_cell_ftm = psc%nbr_cell_ftm + 1
                    psc%cell_ftm(psc%nbr_cell_ftm) = cell
                 endif

              endif ! la cellule voisine est locale

           endif ! si la cellule voisine existe

        enddo ! boucle sur les cellules voisines

     endif ! la cellule est locale

  enddo
  block
    integer::idx(1:psc%nbr_proc_voisin),unw

!    unw=700+myidsp
!    write(unw,*)'PVOIS',psc%proc_voisin
!    write(unw,*)'nbrcf',psc%nbr_cell_frontiere
 !   write(unw,*)psc%cell_frontiere
!    write(unw,*)
    call sort(psc%proc_voisin(1:psc%nbr_proc_voisin),idx)
    psc%nbr_cell_frontiere(1:psc%nbr_proc_voisin)= psc%nbr_cell_frontiere(idx(1:psc%nbr_proc_voisin))
    psc%cell_frontiere(1:psc%nbr_proc_voisin,:)=    psc%cell_frontiere(idx(1:psc%nbr_proc_voisin),:)

!    write(unw,*)'PVOIS',psc%proc_voisin
!    write(unw,*)'nbrcf',psc%nbr_cell_frontiere
!    write(unw,*)psc%cell_frontiere

  end block

    
  if ((myidsp.le.2).or.(rang.ge.nprocspace-2))write(6,*)'rang rangspace ',rang, myidsp,' nbr procs voisins ', psc%nbr_proc_voisin

!call arret_ndm
!  call psc%print(unit=700+myidsp)

end subroutine init_voisinage

subroutine sort(Q,idx)
  implicit none
  integer::n
  integer, intent(inout) :: Q(:)
  integer, intent(out) :: idx(size(Q))

    integer :: i, j
    real    :: key
    integer :: keyidx

    ! Initialize index array
    do i = 1, size(Q)
        idx(i) = i
    end do

    ! Insertion sort
    do i = 2, size(Q)
        key    = Q(i)
        keyidx = idx(i)
        j = i - 1

        do while (j >= 1 )
           if (Q(j) <= key) exit
            Q(j+1)   = Q(j)
            idx(j+1) = idx(j)
            j = j - 1
        end do

        Q(j+1)   = key
        idx(j+1) = keyidx
    end do
  end subroutine sort

end module init_vois_mod
#endif
