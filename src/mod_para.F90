module mod_para
  USE arret_ndm_mod,only:arret_ndm
#ifdef PARA
  use Tpara,only: NDM_MPI_REAL_DOUBLE,MPI_COMM_space, myidsp,nprocspace,nprocs,ierr,status,para_space_config			! numero de process mis là pour être utilisé en sequentiel
  use mpi


#endif
  use T_kind_param_m, ONLY:  double 
  use gen_com_m ,only:l2t,rang
  USE atomconfig,only:atom_config,atom_config_d,atom_config_e,atom_config_arps
  USE cellconfig,only:cell_config
  USE boxconfig, only:box_config
  implicit none
  class(atom_config),pointer:: atmp
  type(cell_config),pointer::celmp

  !  integer :: myidsp,nprocspace,nprocs 			! numero de process mis là pour être utilisé en sequentiel
#ifdef PARA

  integer :: nb_var_int                              ! nbr de variables entieres a envoyer lors des echanges entre proc
  integer :: nb_var_dbl                              ! nbr de variables reelles a envoyer lors des echanges entre proc
  integer :: nb_var_lgc                              ! nbr de variables logical a envoyer lors des echanges entre proc
  integer, allocatable :: send_nb_val(:)             ! buffer d'envoi du nombre de valeurs envoyees
  integer, allocatable :: send_buff_int(:,:,:)       ! buffer d'envoi des variables entieres
  logical, allocatable :: send_buff_lgc(:,:,:)       ! buffer d'envoi des variables logical
  real(double), allocatable :: send_buff_dbl(:,:,:)  ! buffer d'envoi des variables reelles

  integer, allocatable :: recv_nb_val(:)             ! buffer de receptio du nombre de valeurs envoyees
  integer, allocatable :: recv_buff_int(:,:,:)       ! buffer de reception des variables entieres
  logical, allocatable :: recv_buff_lgc(:,:,:)       ! buffer de reception des variables logical
  real(double), allocatable :: recv_buff_dbl(:,:,:)  ! buffer de reception des variables reelles

  integer, allocatable :: send_rqst(:,:)              ! tableau pour stocker les requetes en envoi
  integer, allocatable :: recv_rqst(:,:)              ! tableau pour stocker les requetes en reception

!!$  real(double), allocatable,dimension(:,:)::xp,vp,fp,xpp,ax,glangv
!!$  real(double),allocatable::eat(:),sigat(:,:,:)
!!$  logical, allocatable::lgul(:)
!!$  integer,allocatable,dimension(:)::ityp,ielat,iwmax,num_at_glob,indi
!!$  logical ::llangevin,lprteat,lsigat,ltbv,lax
!!$  integer::im,imm,nvois
!!$  integer:: nox,noy,noz,natperc,noxyz
!!$  integer,allocatable::ncel(:,:),nato(:),atincel(:,:),deltadist(:,:,:),proc_cell(:)
!!$  real(double)::celsize(3)
!!$  logical::ltpcel
!!$  real(double),allocatable::sigc(:,:,:),tempc(:)
  integer::ic1,ic2,ic3,ival
  !---------------------------------------------------------!
  !               Routines spécifique à MPI                 !
  !---------------------------------------------------------!

contains

  !------------------------------------------------------------------------!
  ! Procedure pour la mise a jour des atomes frontieres et fantomes sur
  ! les processeurs. Elle prend en compte la nouvelle repartition dans 
  ! les cellules suite a l'appel a caltabt

  subroutine maj_atomes_frt_ftm(atcf,cellcf,boxcf,psc)

    USE T_kind_param_m, ONLY:  double

    implicit none
    type(cell_config),target::cellcf
    class(atom_config),target::atcf
    type(para_space_config)::psc
    class(box_config)::boxcf

    integer::ne

    atmp=> atcf
    celmp=>cellcf

    call transfert_atomes_fantomes(psc)
    !    call envoi_atomes_fantomes(psc) ! On envoit les atomes qui n'appartiennent plus au processeur courant (qui sont passés  dans des cellules fantomes) caltabt les a mis dans ces cellules fantomes alors qu'ils étaient locaux avant
    !    call reception_nouveaux_atomes(psc) ! On recoit les nouveaux atomes locaux (qui viennent des fantomes des procs voisins)
    call elimine_atomes_fantomes(psc) ! On retire les atomes qui ne sont plus locaux (qui ont été envoyés par envoi_atomes_fantomes)
    !    ne=4
    !    call finalisation_envoi_atomes(ne,psc)     ! Finalisation de l'envoi des atomes pour liberer les buffers d'envoi
    ! En ce point les atomes du proc local sont à jours
    call envoi_atomes_frontieres(psc)     ! On envoit les atomes frontieres aux processeurs voisins
    call reception_atomes_fantomes (psc)    ! On receptionne les nouveaux atomes fantomes
    call finalisation_envoi_atomes(ne,psc)     ! Finalisation de l'envoi des atomes pour liberer les buffers d'envoi

  end subroutine maj_atomes_frt_ftm

  !------------------------------------------------------------------------!
  ! Procedure pour la mise a jour des valeurs tabdensity des atomes 
  ! fantomes sur les processeurs.

  subroutine maj_atomes_frt_part(atcf,cellcf,boxcf,psc,caracT)
    implicit none
    type(cell_config),target::cellcf
    class(atom_config),target::atcf
    type(para_space_config)::psc
    type(box_config)::boxcf
    character(len=*),intent(in)::caracT

    integer::ne    ,nvi,nvl,nvr
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   
    ! v=vp,r=xpp
    ! u=eat,g=glangv;a=ax;s=sigat; m=mov(arps)

    atmp=> atcf
    celmp=>cellcf
    ne=0
    call envoi_atomes_frontieres_part(psc,caracT,nvi,nvl,nvr)     ! On envoit les atomes frontieres aux processeurs voisins
    call reception_atomes_fantomes_part (psc,caracT,nvi,nvl,nvr)    ! On receptionne les nouveaux atomes fantomes
    call finalisation_envoi_atomes(ne,psc,nvi,nvl,nvr)     ! Finalisation de l'envoi des atomes pour liberer les buffers d'envoi


  end subroutine maj_atomes_frt_part
  subroutine maj_tabdensity_ftm(tabdensity, imm, num_at_glob, psc, im)

    use T_kind_param_m, only : double
    implicit none

    type(para_space_config) :: psc
    integer, intent(in)     :: imm, im
    integer, intent(in)     :: num_at_glob(imm)
    real(double)            :: tabdensity(imm)

    call exchange_tabdensity_blocking(tabdensity, imm, num_at_glob, psc, im)

  end subroutine maj_tabdensity_ftm


  subroutine exchange_tabdensity_blocking(tabdensity, imm, num_at_glob, psc, im)

    use T_kind_param_m, only : double
    implicit none

    type(para_space_config) :: psc
    integer, intent(in)     :: imm, im
    integer, intent(in)     :: num_at_glob(imm)
    real(double)            :: tabdensity(imm)

    integer :: nproc_voisin, procv
    integer :: nb_at_send, nb_at_recv
    integer :: ncell_front, koo, n_at, i_at
    integer :: i, ind_loc, ftm_at

    integer, allocatable :: send_ids(:), recv_ids(:)
    real(double), allocatable :: send_val(:), recv_val(:)

    do nproc_voisin = 1, psc%nbr_proc_voisin

       procv = psc%proc_voisin(nproc_voisin)

       !--------------------------------------------------
       ! Build send buffers
       !--------------------------------------------------
       nb_at_send = 0
       do ncell_front = 1, psc%nbr_cell_frontiere(nproc_voisin)
          koo = psc%cell_frontiere(nproc_voisin, ncell_front)
          nb_at_send = nb_at_send + celmp%nato(koo)
       end do

       allocate(send_ids(nb_at_send))
       allocate(send_val(nb_at_send))

       i = 0
       do ncell_front = 1, psc%nbr_cell_frontiere(nproc_voisin)
          koo = psc%cell_frontiere(nproc_voisin, ncell_front)
          do n_at = 1, celmp%nato(koo)
             i_at = celmp%atincel(n_at, koo)
             i = i + 1
             send_ids(i) = atmp%num_at_glob(i_at)
             send_val(i) = tabdensity(i_at)
          end do
       end do

       !--------------------------------------------------
       ! Ordered blocking communication (deadlock-safe)
       !--------------------------------------------------
       if (myidsp < procv) then

          call MPI_SEND(nb_at_send, 1, MPI_INTEGER, procv, 3001, MPI_COMM_space, ierr)
          call MPI_SEND(send_ids, nb_at_send, MPI_INTEGER, procv, 3002, MPI_COMM_space, ierr)
          call MPI_SEND(send_val, nb_at_send, NDM_MPI_REAL_DOUBLE, procv, 3003, MPI_COMM_space, ierr)

          call MPI_RECV(nb_at_recv, 1, MPI_INTEGER, procv, 3001, MPI_COMM_space, status, ierr)
          allocate(recv_ids(nb_at_recv), recv_val(nb_at_recv))
          call MPI_RECV(recv_ids, nb_at_recv, MPI_INTEGER, procv, 3002, MPI_COMM_space, status, ierr)
          call MPI_RECV(recv_val, nb_at_recv, NDM_MPI_REAL_DOUBLE, procv, 3003, MPI_COMM_space, status, ierr)

       else

          call MPI_RECV(nb_at_recv, 1, MPI_INTEGER, procv, 3001, MPI_COMM_space, status, ierr)
          allocate(recv_ids(nb_at_recv), recv_val(nb_at_recv))
          call MPI_RECV(recv_ids, nb_at_recv, MPI_INTEGER, procv, 3002, MPI_COMM_space, status, ierr)
          call MPI_RECV(recv_val, nb_at_recv, NDM_MPI_REAL_DOUBLE, procv, 3003, MPI_COMM_space, status, ierr)

          call MPI_SEND(nb_at_send, 1, MPI_INTEGER, procv, 3001, MPI_COMM_space, ierr)
          call MPI_SEND(send_ids, nb_at_send, MPI_INTEGER, procv, 3002, MPI_COMM_space, ierr)
          call MPI_SEND(send_val, nb_at_send, NDM_MPI_REAL_DOUBLE, procv, 3003, MPI_COMM_space, ierr)

       end if

       !--------------------------------------------------
       ! Update ghost atoms
       !--------------------------------------------------
       do i = 1, nb_at_recv
          ind_loc = -1
          do ftm_at = im+1, imm
             if (atmp%num_at_glob(ftm_at) == recv_ids(i)) then
                ind_loc = ftm_at
                exit
             end if
          end do

          if (ind_loc < 0) then
             write(6,*) myidsp, 'Error: ghost atom not found',i,recv_ids(i),nproc_voisin,procv
             write(6,*)recv_ids
             call MPI_ABORT(MPI_COMM_space, 1, ierr)
          end if

          tabdensity(ind_loc) = recv_val(i)
       end do

       deallocate(send_ids, send_val, recv_ids, recv_val)

    end do

  end subroutine exchange_tabdensity_blocking






  ! Procedure pour la mise a jour des valeurs fp des atomes frontieres
  ! du processeur courant avec leurs contributions des processeurs voisins

  subroutine maj_fp_frt(psc,atcf,celcf) !appelée SEULEMENT dans force_tersoff_cel !

    USE T_kind_param_m, ONLY:  double

    implicit none
    integer::ne
    type(para_space_config)::psc
    class(atom_config),intent(inout)::atcf
    type(cell_config),intent(in)::celcf

    ! On envoit les atomes fantomes vers les processeurs voisins
    call envoi_fp_fantomes(psc,atcf,celcf)

    ! On receptionne les contributions des processeurs voisins
    call reception_fp_frontieres(psc,atcf,celcf)

    ! Finalisation de l'envoi pour liberer les buffers d'envoi (identique a l'envoi des atomes)
    ne=3
    call finalisation_envoi_atomes(ne,psc)

  end subroutine maj_fp_frt

  subroutine transfert_atomes_fantomes(psc) ! seulement maj_atomes_frt_ftm
    type(para_space_config)::psc
    integer :: nb_at, nb_at_max, nb_at_max_tot
    integer :: nproc_voisin
    integer :: ncell_ftm
    integer :: procv,cellf,n_at,i_at
  integer :: nb_var_int                              ! nbr de variables entieres a envoyer lors des echanges entre proc
  integer :: nb_var_dbl                              ! nbr de variables reelles a envoyer lors des echanges entre proc
  integer :: nb_var_lgc                              ! nbr de variables logical a envoyer lors des echanges entre proc
  integer :: send_nb_val,nb_at_send
  integer, allocatable :: send_buff_int(:,:)       ! buffer d'envoi des variables entieres
  logical, allocatable :: send_buff_lgc(:,:)       ! buffer d'envoi des variables logical
  real(double), allocatable :: send_buff_dbl(:,:)  ! buffer d'envoi des variables reelles

  integer :: recv_nb_val             ! buffer de receptio du nombre de valeurs envoyees
  integer, allocatable :: recv_buff_int(:,:)       ! buffer de reception des variables entieres
  logical, allocatable :: recv_buff_lgc(:,:)       ! buffer de reception des variables logical
  real(double), allocatable :: recv_buff_dbl(:,:)  ! buffer de reception des variables reelles



    nb_var_int = 4 !ityp,num_at_glob,ielat,proc_at
    nb_var_lgc=1 !lgul
    !    nb_var_int = 4 avec iwmax aucun intérêt
    !LPARAFULLSEND
    !    nb_var_dbl = 15
!!$   if (lsuivinonpbc) then
!!$    nb_var_dbl = 18
!!$   else 
!!$    nb_var_dbl = 9
!!$ end if
    !    nb_var_dbl = 9 avec ax


    nb_var_dbl = 6  !  xp,fp
    select type (atmp)
    class is (atom_config_d)
       nb_var_dbl=nb_var_dbl+3 !vp xpp
    end select
    select type (atmp)
    class is (atom_config_e)
       !       nb_var_dbl=nb_var_dbl+6 !vp xpp
       if (atmp%lxpp)        nb_var_dbl=nb_var_dbl+3 !eat
       if (atmp%lprteat)        nb_var_dbl=nb_var_dbl+1 !eat
       if (atmp%lsigat)        nb_var_dbl=nb_var_dbl+9 !eat
       if (atmp%llangevin)        nb_var_dbl=nb_var_dbl+3 !eat
       if (atmp%lax)        nb_var_dbl=nb_var_dbl+3 !eat
    end select
    select type (atmp)
    class is (atom_config_arps)
       nb_var_dbl=nb_var_dbl+3 !fpr
       if (allocated(atmp%rho)) nb_var_dbl=nb_var_dbl+1  ! ,rho
       nb_var_int = nb_var_int+1 !mov
       if (allocated(atmp%fpg)) nb_var_dbl=nb_var_dbl+3
    end select




    do nproc_voisin = 1, psc%nbr_proc_voisin

       procv = psc%proc_voisin(nproc_voisin)
       nb_at_send = 0
       send_nb_val=0
       ! Boucle sur les cellules fantomes
       do ncell_ftm=1,psc%nbr_cell_ftm
          cellf = psc%cell_ftm(ncell_ftm)
          ! Sommation des atomes de la cellule
          if (celmp%proc_cell(cellf)==procv) nb_at_send = nb_at_send + celmp%nato(cellf)
       enddo ! fin boucle sur les cellules

       !--------------------------------------------------
       ! Build send buffers
       !--------------------------------------------------

       allocate(send_buff_int(nb_var_int,nb_at_send))
       allocate(send_buff_lgc(nb_var_lgc,nb_at_send))
       allocate(send_buff_dbl(nb_var_dbl,nb_at_send))
       !    allocate(recv_nb_val(psc%nbr_proc_voisin))   ! nombre effectif d'atomes reçus du proc voisin

       ! On boucle sur les cellules fantomes associees a ce processeur voisin
       do ncell_ftm= 1, psc%nbr_cell_ftm
          cellf = psc%cell_ftm(ncell_ftm)
          if (celmp%proc_cell(cellf)==procv) then

             ! On boucle sur les atomes de cette cellule (ce ne sont que des anciens atomes locaux car la liste est construite par caltabt sur les atomes de i à im).
             do n_at= 1, celmp%nato(cellf)
                i_at = celmp%atincel(n_at,cellf)
                ! On complete le buffer
                send_nb_val = send_nb_val + 1

                send_buff_lgc(1,send_nb_val) = atmp%lgul(i_at)

                send_buff_int(1,send_nb_val) = atmp%ityp(i_at)
                send_buff_int(2,send_nb_val) = atmp%ielat(i_at)
                send_buff_int(3,send_nb_val) = atmp%num_at_glob(i_at)
                select type (atmp)
                class is (atom_config_arps)
                   send_buff_int(4,send_nb_val) = atmp%mov(i_at)
                end select

                send_buff_dbl(1,send_nb_val) = atmp%xp(1,i_at)
                send_buff_dbl(2,send_nb_val) = atmp%xp(2,i_at)
                send_buff_dbl(3,send_nb_val) = atmp%xp(3,i_at)
                send_buff_dbl(4,send_nb_val) = atmp%fp(1,i_at)
                send_buff_dbl(5,send_nb_val) = atmp%fp(2,i_at)
                send_buff_dbl(6,send_nb_val) = atmp%fp(3,i_at)
                select type (atmp)
                class is (atom_config_d)
                   send_buff_dbl(7,send_nb_val) = atmp%vp(1,i_at)
                   send_buff_dbl(8,send_nb_val) = atmp%vp(2,i_at)
                   send_buff_dbl(9,send_nb_val) = atmp%vp(3,i_at)
                end select
                ival=9
                select type (atmp)
                class is (atom_config_e)
                   if (atmp%lxpp) then
                      do ic3=1,3
                         ival =ival+1; send_buff_dbl(ival,send_nb_val) = atmp%xpp(ic3,i_at)
                      end do
                   end if

                   if (atmp%lprteat)then
                      ival =ival+1; send_buff_dbl(ival,send_nb_val) = atmp%eat(i_at)
                   end if
                   if (atmp%lsigat)then
                      do ic1=1,3
                         do ic2=1,3
                            ival =ival+1; send_buff_dbl(ival,send_nb_val) = atmp%sigat(ic1,ic2,i_at)
                         end do
                      end do
                   end if
                   if (atmp%llangevin) then
                      do ic3=1,3
                         ival =ival+1; send_buff_dbl(ival,send_nb_val) = atmp%glangv(ic3,i_at)
                      end do
                   end if
                   if (atmp%lax) then
                      do ic3=1,3
                         ival =ival+1; send_buff_dbl(ival,send_nb_val) = atmp%ax(ic3,i_at)
                      end do
                   end if
                end select
                select type (atmp)
                class is (atom_config_arps)
                   do ic3=1,3
                      ival =ival+1; send_buff_dbl(ival,send_nb_val) = atmp%fpr(ic3,i_at)
                   end do
                   if (allocated(atmp%rho)) then 
                      ival =ival+1; send_buff_dbl(ival,send_nb_val) = atmp%rho(i_at)
                   end if
                   if (allocated(atmp%fpg)) then
                      do ic3=1,3
                         ival =ival+1; send_buff_dbl(ival,send_nb_val) = atmp%fpg(ic3,i_at)
                      end do
                   end if
                end select


             enddo  ! fin de boucle sur les atomes
          endif
       enddo    ! fin de boucle sur les cellules fantomes
       !       write(6,*)'nouveaux ATOMOUTP',rang,myidsp, im,send_nb_val(nproc_voisin),procv


       if (send_nb_val.ne.nb_at_send) then
          write(6,*)'rang send_nb_val.ne.nb_at_send',myidsp,nproc_voisin, procv,send_nb_val,nb_at_send
          call arret_ndm
       end if



       !--------------------------------------------------
       ! Ordered blocking communication (deadlock-safe)
       !--------------------------------------------------
       if (myidsp < procv) then
          call MPI_SEND(send_nb_val,      1,   MPI_INTEGER,        &
               procv,3001,MPI_COMM_space,ierr)

          call MPI_SEND(send_buff_int,nb_var_int*send_nb_val,MPI_INTEGER,        &
               procv,3002,MPI_COMM_space,ierr)

          call MPI_SEND(send_buff_dbl,nb_var_dbl*send_nb_val,NDM_MPI_REAL_DOUBLE,&
               procv,3003,MPI_COMM_space,ierr)
          call MPI_SEND(send_buff_lgc,nb_var_lgc*send_nb_val,MPI_LOGICAL,        &
               procv,3004,MPI_COMM_space,ierr)

          call MPI_RECV(recv_nb_val, 1, MPI_INTEGER, procv, 3001, MPI_COMM_space, status, ierr)
          allocate(recv_buff_int(nb_var_int,recv_nb_val))
          allocate(recv_buff_dbl(nb_var_dbl,recv_nb_val))
          allocate(recv_buff_lgc(nb_var_lgc,recv_nb_val))
          call MPI_RECV(recv_buff_int(1,1), nb_var_int*recv_nb_val, MPI_INTEGER, procv, 3002, &
               MPI_COMM_space, status,  ierr)
          call MPI_RECV(recv_buff_dbl(1,1), nb_var_dbl*recv_nb_val, NDM_MPI_REAL_DOUBLE, procv, 3003, &
               MPI_COMM_space, status,  ierr)
          call MPI_RECV(recv_buff_lgc(1,1), nb_var_lgc*recv_nb_val, MPI_LOGICAL, procv, 3004, &
               MPI_COMM_space, status,  ierr)

       else

          call MPI_RECV(recv_nb_val, 1, MPI_INTEGER, procv, 3001, MPI_COMM_space, status, ierr)
          allocate(recv_buff_int(nb_var_int,recv_nb_val))
          allocate(recv_buff_dbl(nb_var_dbl,recv_nb_val))
          allocate(recv_buff_lgc(nb_var_lgc,recv_nb_val))
          call MPI_RECV(recv_buff_int(1,1), nb_var_int*recv_nb_val, MPI_INTEGER, procv, 3002, &
               MPI_COMM_space, status,  ierr)
          call MPI_RECV(recv_buff_dbl(1,1), nb_var_dbl*recv_nb_val, NDM_MPI_REAL_DOUBLE, procv, 3003, &
               MPI_COMM_space, status,  ierr)
          call MPI_RECV(recv_buff_lgc(1,1), nb_var_lgc*recv_nb_val, MPI_LOGICAL, procv, 3004, &
               MPI_COMM_space, status,  ierr)

          call MPI_SEND(send_nb_val,      1,   MPI_INTEGER,        &
               procv,3001,MPI_COMM_space,ierr)

          call MPI_SEND(send_buff_int,nb_var_int*send_nb_val,MPI_INTEGER,        &
               procv,3002,MPI_COMM_space,ierr)

          call MPI_SEND(send_buff_dbl,nb_var_dbl*send_nb_val,NDM_MPI_REAL_DOUBLE,&
               procv,3003,MPI_COMM_space,ierr)
          call MPI_SEND(send_buff_lgc,nb_var_lgc*send_nb_val,MPI_LOGICAL,        &
               procv,3004,MPI_COMM_space,ierr)

       end if

       do i_at = 1, recv_nb_val


          ! On ajoute un atome a la liste
          atmp%im = atmp%im + 1
          !          imd = imd + 1
          !          imf = imf + 1
          !          imana = imana + 1
          atmp%lgul(atmp%im)        = recv_buff_lgc(1,i_at)

          ! mise a jour des variables entieres
          atmp%ityp(atmp%im)        = recv_buff_int(1,i_at)
          atmp%ielat(atmp%im)       = recv_buff_int(2,i_at)
          !          iwmax(im)       = recv_buff_int(3,i_at)
          atmp%num_at_glob(atmp%im) = recv_buff_int(3,i_at)
          select type (atmp)
          class is (atom_config_arps)
             atmp%mov(atmp%im) = recv_buff_int(4,i_at)
          end select
          atmp%proc_at(atmp%im)=rang
          ! mise a jour des donnees de la cellule correspondante
          celmp%nato(atmp%ielat(atmp%im)) = celmp%nato(atmp%ielat(atmp%im)) + 1
          celmp%atincel(celmp%nato(atmp%ielat(atmp%im)),atmp%ielat(atmp%im)) = atmp%im

          ! mise a jour des variables reelles
          atmp%xp(1,atmp%im) = recv_buff_dbl(1,i_at) 
          atmp%xp(2,atmp%im) = recv_buff_dbl(2,i_at) 
          atmp%xp(3,atmp%im) = recv_buff_dbl(3,i_at) 
          atmp%fp(1,atmp%im) = recv_buff_dbl(4,i_at) 
          atmp%fp(2,atmp%im) = recv_buff_dbl(5,i_at) 
          atmp%fp(3,atmp%im) = recv_buff_dbl(6,i_at) 


          select type (atmp)
          class is (atom_config_d)
             atmp%vp(1,atmp%im) = recv_buff_dbl(7,i_at) 
             atmp%vp(2,atmp%im) = recv_buff_dbl(8,i_at) 
             atmp%vp(3,atmp%im) = recv_buff_dbl(9,i_at)
          end select
          ival=9
          select type (atmp)
          class is (atom_config_e)
!!$             atmp%vp(1,atmp%im) = recv_buff_dbl(7,i_at) 
!!$             atmp%vp(2,atmp%im) = recv_buff_dbl(8,i_at) 
!!$             atmp%vp(3,atmp%im) = recv_buff_dbl(9,i_at)
             if (atmp%lxpp) then
                do ic3=1,3
                   ival =ival+1;atmp%xpp(ic3,atmp%im)= recv_buff_dbl(ival,i_at)
                end do
             end if

             if (atmp%lprteat)then
                ival =ival+1
                atmp%eat(atmp%im)=recv_buff_dbl(ival,i_at)
             end if
             if (atmp%lsigat)then
                do ic1=1,3
                   do ic2=1,3
                      ival =ival+1
                      atmp%sigat(ic1,ic2,atmp%im)=recv_buff_dbl(ival,i_at)
                   end do
                end do
             end if
             if (atmp%llangevin) then
                do ic3=1,3
                   ival =ival+1;atmp%glangv(ic3,atmp%im)= recv_buff_dbl(ival,i_at)
                end do
             end if
             if (atmp%lax) then
                do ic3=1,3
                   ival =ival+1;atmp%ax(ic3,atmp%im)= recv_buff_dbl(ival,i_at)
                end do
             end if
          end select
          select type (atmp)
          class is (atom_config_arps)
             do ic3=1,3
                ival =ival+1; atmp%fpr(ic3,atmp%im)= recv_buff_dbl(ival,i_at)
             end do
             if (allocated(atmp%rho)) then
                ival =ival+1; atmp%rho(atmp%im)= recv_buff_dbl(ival,i_at)
             end if
             if (allocated(atmp%fpg)) then
                do ic3=1,3
                   ival =ival+1; atmp%fpg(ic3,atmp%im)= recv_buff_dbl(ival,i_at)
                end do
             end if

          end select

       enddo

       deallocate(send_buff_int)
       deallocate(send_buff_lgc)
       deallocate(send_buff_dbl)
       deallocate(recv_buff_int)
       deallocate(recv_buff_dbl)
       deallocate(recv_buff_lgc)



    end do



  end subroutine transfert_atomes_fantomes



  !------------------------------------------------------------------------!
  ! Procedure dont le but est l'envoi des atomes qui sont sorti du domaine
  ! courant pour etre pris en charge par leur nouveau processeur

  subroutine envoi_atomes_fantomes(psc) ! seulement maj_atomes_frt_ftm


    USE T_kind_param_m, ONLY:  double
    implicit none
    type(para_space_config)::psc
    integer :: nb_at, nb_at_max, nb_at_max_tot
    integer :: nproc_voisin
    integer :: ncell_ftm
    integer :: procv,cellf,n_at,i_at

    ! Premier passage a vide pour allouer les buffers au plus juste
    !    write(6,*)'envoi_atomes_fantomes',rang,nbr_proc_voisin
    nb_at_max=0
    ! Boucle sur les processeurs voisins
    do nproc_voisin=1,psc%nbr_proc_voisin
       procv = psc%proc_voisin(nproc_voisin)
       nb_at = 0
       ! Boucle sur les cellules fantomes
       do ncell_ftm=1,psc%nbr_cell_ftm
          cellf = psc%cell_ftm(ncell_ftm)
          ! Sommation des atomes de la cellule
          if (celmp%proc_cell(cellf)==procv) nb_at = nb_at + celmp%nato(cellf)
       enddo ! fin boucle sur les cellules
       ! calcul du max des atomes a envoyer
       nb_at_max = max(nb_at_max, nb_at)
    enddo ! fin boucle sur les processeurs voisins
    call MPI_ALLREDUCE(nb_at_max,nb_at_max_tot,1,MPI_INTEGER,MPI_MAX,MPI_COMM_space,ierr)
    nb_at_max=nb_at_max_tot
    ! petite manip pour eviter le cas nb_at_max=0

    nb_at_max = max(nb_at_max,1) ! nb max d'atomes dans les cellules fantomes autour du proc courant

    ! Allocation des buffers


    nb_var_int = 4 !ityp,num_at_glob,ielat,proc_at
    nb_var_lgc=1 !lgul
    !    nb_var_int = 4 avec iwmax aucun intérêt
    !LPARAFULLSEND
    !    nb_var_dbl = 15
!!$   if (lsuivinonpbc) then
!!$    nb_var_dbl = 18
!!$   else 
!!$    nb_var_dbl = 9
!!$ end if
    !    nb_var_dbl = 9 avec ax


    nb_var_dbl = 6  !  xp,fp
    select type (atmp)
    class is (atom_config_d)
       nb_var_dbl=nb_var_dbl+3 !vp xpp
    end select
    select type (atmp)
    class is (atom_config_e)
       !       nb_var_dbl=nb_var_dbl+6 !vp xpp
       if (atmp%lxpp)        nb_var_dbl=nb_var_dbl+3 !eat
       if (atmp%lprteat)        nb_var_dbl=nb_var_dbl+1 !eat
       if (atmp%lsigat)        nb_var_dbl=nb_var_dbl+9 !eat
       if (atmp%llangevin)        nb_var_dbl=nb_var_dbl+3 !eat
       if (atmp%lax)        nb_var_dbl=nb_var_dbl+3 !eat
    end select
    select type (atmp)
    class is (atom_config_arps)
       nb_var_dbl=nb_var_dbl+3 !fpr
       if (allocated(atmp%rho)) nb_var_dbl=nb_var_dbl+1  ! ,rho
       nb_var_int = nb_var_int+1 !mov
       if (allocated(atmp%fpg)) nb_var_dbl=nb_var_dbl+3
    end select

!!$  if ((llangevin.eqv..true.).or.(l2T.eqv..true.))then
!!$     nb_var_dbl = nb_var_dbl+3
!!$  end if

    !    write(6,*)'nb_var_dbl1',nb_var_dbl
    allocate(send_nb_val(psc%nbr_proc_voisin))
    allocate(send_buff_int(nb_var_int,nb_at_max,psc%nbr_proc_voisin))
    allocate(send_buff_lgc(nb_var_lgc,nb_at_max,psc%nbr_proc_voisin))
    allocate(send_buff_dbl(nb_var_dbl,nb_at_max,psc%nbr_proc_voisin))
    allocate(send_rqst(psc%nbr_proc_voisin,4))
    allocate(recv_nb_val(psc%nbr_proc_voisin))   ! nombre effectif d'atomes reçus du proc voisin
    allocate(recv_buff_int(nb_var_int,nb_at_max,psc%nbr_proc_voisin))
    allocate(recv_buff_dbl(nb_var_dbl,nb_at_max,psc%nbr_proc_voisin))
    allocate(recv_buff_lgc(nb_var_lgc,nb_at_max,psc%nbr_proc_voisin))
    allocate(recv_rqst(psc%nbr_proc_voisin,4))

    ! Preparation des receptions
    do nproc_voisin= 1, psc%nbr_proc_voisin
       procv = psc%proc_voisin(nproc_voisin)
       call MPI_IRECV(recv_nb_val(nproc_voisin), 1, MPI_INTEGER, procv, 1001, MPI_COMM_space, recv_rqst(nproc_voisin,1), ierr)
       call MPI_IRECV(recv_buff_int(1,1,nproc_voisin), nb_var_int*nb_at_max, MPI_INTEGER, procv, 1002, &
            MPI_COMM_space, recv_rqst(nproc_voisin,2), ierr)
       call MPI_IRECV(recv_buff_dbl(1,1,nproc_voisin), nb_var_dbl*nb_at_max, NDM_MPI_REAL_DOUBLE, procv, 1003, &
            MPI_COMM_space, recv_rqst(nproc_voisin,3), ierr)
       call MPI_IRECV(recv_buff_lgc(1,1,nproc_voisin), nb_var_lgc*nb_at_max, MPI_LOGICAL, procv, 1004, &
            MPI_COMM_space, recv_rqst(nproc_voisin,4), ierr)
    enddo
    ! On boucle sur les processeurs voisins
    do nproc_voisin= 1, psc%nbr_proc_voisin
       procv = psc%proc_voisin(nproc_voisin)

       send_nb_val(nproc_voisin) = 0

       ! On boucle sur les cellules fantomes associees a ce processeur voisin
       do ncell_ftm= 1, psc%nbr_cell_ftm
          cellf = psc%cell_ftm(ncell_ftm)
          if (celmp%proc_cell(cellf)==procv) then

             ! On boucle sur les atomes de cette cellule (ce ne sont que des anciens atomes locaux car la liste est construite par caltabt sur les atomes de i à im).
             do n_at= 1, celmp%nato(cellf)
                i_at = celmp%atincel(n_at,cellf)
                ! On complete le buffer
                send_nb_val(nproc_voisin) = send_nb_val(nproc_voisin) + 1

                send_buff_lgc(1,send_nb_val(nproc_voisin),nproc_voisin) = atmp%lgul(i_at)

                send_buff_int(1,send_nb_val(nproc_voisin),nproc_voisin) = atmp%ityp(i_at)
                send_buff_int(2,send_nb_val(nproc_voisin),nproc_voisin) = atmp%ielat(i_at)
                !                send_buff_int(3,send_nb_val(nproc_voisin),nproc_voisin) = iwmax(i_at)
                send_buff_int(3,send_nb_val(nproc_voisin),nproc_voisin) = atmp%num_at_glob(i_at)
                select type (atmp)
                class is (atom_config_arps)
                   send_buff_int(4,send_nb_val(nproc_voisin),nproc_voisin) = atmp%mov(i_at)
                end select

                send_buff_dbl(1,send_nb_val(nproc_voisin),nproc_voisin) = atmp%xp(1,i_at)
                send_buff_dbl(2,send_nb_val(nproc_voisin),nproc_voisin) = atmp%xp(2,i_at)
                send_buff_dbl(3,send_nb_val(nproc_voisin),nproc_voisin) = atmp%xp(3,i_at)
                send_buff_dbl(4,send_nb_val(nproc_voisin),nproc_voisin) = atmp%fp(1,i_at)
                send_buff_dbl(5,send_nb_val(nproc_voisin),nproc_voisin) = atmp%fp(2,i_at)
                send_buff_dbl(6,send_nb_val(nproc_voisin),nproc_voisin) = atmp%fp(3,i_at)
!!$                send_buff_dbl(4,send_nb_val(nproc_voisin),nproc_voisin) = ax(1,i_at)
!!$                send_buff_dbl(5,send_nb_val(nproc_voisin),nproc_voisin) = ax(2,i_at)
!!$                send_buff_dbl(6,send_nb_val(nproc_voisin),nproc_voisin) = ax(3,i_at)
                select type (atmp)
                class is (atom_config_d)
                   send_buff_dbl(7,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(1,i_at)
                   send_buff_dbl(8,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(2,i_at)
                   send_buff_dbl(9,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(3,i_at)
                end select
                ival=9
                select type (atmp)
                class is (atom_config_e)
!!$                   send_buff_dbl(7,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(1,i_at)
!!$                   send_buff_dbl(8,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(2,i_at)
!!$                   send_buff_dbl(9,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(3,i_at)
                   if (atmp%lxpp) then
                      do ic3=1,3
                         ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%xpp(ic3,i_at)
                      end do
                   end if

                   if (atmp%lprteat)then
                      ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%eat(i_at)
                   end if
                   if (atmp%lsigat)then
                      do ic1=1,3
                         do ic2=1,3
                            ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%sigat(ic1,ic2,i_at)
                         end do
                      end do
                   end if
                   if (atmp%llangevin) then
                      do ic3=1,3
                         ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%glangv(ic3,i_at)
                      end do
                   end if
                   if (atmp%lax) then
                      do ic3=1,3
                         ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%ax(ic3,i_at)
                      end do
                   end if
                end select
                select type (atmp)
                class is (atom_config_arps)
                   do ic3=1,3
                      ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%fpr(ic3,i_at)
                   end do
                   if (allocated(atmp%rho)) then 
                      ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%rho(i_at)
                   end if
                   if (allocated(atmp%fpg)) then
                      do ic3=1,3
                         ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%fpg(ic3,i_at)
                      end do
                   end if
                end select


             enddo  ! fin de boucle sur les atomes
          endif
       enddo    ! fin de boucle sur les cellules fantomes
       !       write(6,*)'nouveaux ATOMOUTP',rang,myidsp, im,send_nb_val(nproc_voisin),procv

       ! On envoit les buffers vers le processeur
       call MPI_ISSEND(send_nb_val(nproc_voisin),      1,   MPI_INTEGER,        &
            procv,1001,MPI_COMM_space,send_rqst(nproc_voisin,1),ierr)

       call MPI_ISSEND(send_buff_int(1,1,nproc_voisin),nb_var_int*nb_at_max,MPI_INTEGER,        &
            procv,1002,MPI_COMM_space,send_rqst(nproc_voisin,2),ierr)

       call MPI_ISSEND(send_buff_dbl(1,1,nproc_voisin),nb_var_dbl*nb_at_max,NDM_MPI_REAL_DOUBLE,&
            procv,1003,MPI_COMM_space,send_rqst(nproc_voisin,3),ierr)
       call MPI_ISSEND(send_buff_lgc(1,1,nproc_voisin),nb_var_lgc*nb_at_max,MPI_LOGICAL,        &
            procv,1004,MPI_COMM_space,send_rqst(nproc_voisin,4),ierr)

    enddo    ! fin de boucle sur les processeurs


  end subroutine envoi_atomes_fantomes


  !------------------------------------------------------------------------!
  ! Procedure dont le but est la reception des nouveaux atomes locaux

  subroutine reception_nouveaux_atomes(psc) !seulment MAJ

    USE T_kind_param_m, ONLY:  double

    implicit none
    type(para_space_config)::psc
    integer :: proc_source
    integer :: i_at
    integer :: nproc_voisin
    integer :: ind_recv

    ! On boucle sur les processeurs voisins
    do nproc_voisin=1,psc%nbr_proc_voisin

       call MPI_WAITANY(psc%nbr_proc_voisin,recv_rqst(:,1),ind_recv,status,ierr)
       proc_source = psc%proc_voisin(ind_recv)

       ! Reception des nouveaux atomes issus de ce processeur voisin

       call MPI_WAIT(recv_rqst(ind_recv,2), status, ierr)
       call MPI_WAIT(recv_rqst(ind_recv,3), status, ierr)
       call MPI_WAIT(recv_rqst(ind_recv,4), status, ierr)

       ! recopie des infos dans les tableaux locaux
       !       write(6,*)'nouveaux ATOMINTP',rang,myidsp, im,recv_nb_val(ind_recv),proc_source
       !       write(6,*)'indice nbv',ind_recv,recv_nb_val(ind_recv)
       do i_at = 1, recv_nb_val(ind_recv)


          ! On ajoute un atome a la liste
          atmp%im = atmp%im + 1
          !          imd = imd + 1
          !          imf = imf + 1
          !          imana = imana + 1
          atmp%lgul(atmp%im)        = recv_buff_lgc(1,i_at,ind_recv)

          ! mise a jour des variables entieres
          atmp%ityp(atmp%im)        = recv_buff_int(1,i_at,ind_recv)
          atmp%ielat(atmp%im)       = recv_buff_int(2,i_at,ind_recv)
          !          iwmax(im)       = recv_buff_int(3,i_at,ind_recv)
          atmp%num_at_glob(atmp%im) = recv_buff_int(3,i_at,ind_recv)
          select type (atmp)
          class is (atom_config_arps)
             atmp%mov(atmp%im) = recv_buff_int(4,i_at,ind_recv)
          end select
          atmp%proc_at(atmp%im)=rang
          ! mise a jour des donnees de la cellule correspondante
          celmp%nato(atmp%ielat(atmp%im)) = celmp%nato(atmp%ielat(atmp%im)) + 1
          celmp%atincel(celmp%nato(atmp%ielat(atmp%im)),atmp%ielat(atmp%im)) = atmp%im

          ! mise a jour des variables reelles
          atmp%xp(1,atmp%im) = recv_buff_dbl(1,i_at,ind_recv) 
          atmp%xp(2,atmp%im) = recv_buff_dbl(2,i_at,ind_recv) 
          atmp%xp(3,atmp%im) = recv_buff_dbl(3,i_at,ind_recv) 
          atmp%fp(1,atmp%im) = recv_buff_dbl(4,i_at,ind_recv) 
          atmp%fp(2,atmp%im) = recv_buff_dbl(5,i_at,ind_recv) 
          atmp%fp(3,atmp%im) = recv_buff_dbl(6,i_at,ind_recv) 


          select type (atmp)
          class is (atom_config_d)
             atmp%vp(1,atmp%im) = recv_buff_dbl(7,i_at,ind_recv) 
             atmp%vp(2,atmp%im) = recv_buff_dbl(8,i_at,ind_recv) 
             atmp%vp(3,atmp%im) = recv_buff_dbl(9,i_at,ind_recv)
          end select
          ival=9
          select type (atmp)
          class is (atom_config_e)
!!$             atmp%vp(1,atmp%im) = recv_buff_dbl(7,i_at,ind_recv) 
!!$             atmp%vp(2,atmp%im) = recv_buff_dbl(8,i_at,ind_recv) 
!!$             atmp%vp(3,atmp%im) = recv_buff_dbl(9,i_at,ind_recv)
             if (atmp%lxpp) then
                do ic3=1,3
                   ival =ival+1;atmp%xpp(ic3,atmp%im)= recv_buff_dbl(ival,i_at,ind_recv)
                end do
             end if

             if (atmp%lprteat)then
                ival =ival+1
                atmp%eat(atmp%im)=recv_buff_dbl(ival,i_at,ind_recv)
             end if
             if (atmp%lsigat)then
                do ic1=1,3
                   do ic2=1,3
                      ival =ival+1
                      atmp%sigat(ic1,ic2,atmp%im)=recv_buff_dbl(ival,i_at,ind_recv)
                   end do
                end do
             end if
             if (atmp%llangevin) then
                do ic3=1,3
                   ival =ival+1;atmp%glangv(ic3,atmp%im)= recv_buff_dbl(ival,i_at,ind_recv)
                end do
             end if
             if (atmp%lax) then
                do ic3=1,3
                   ival =ival+1;atmp%ax(ic3,atmp%im)= recv_buff_dbl(ival,i_at,ind_recv)
                end do
             end if
          end select
          select type (atmp)
          class is (atom_config_arps)
             do ic3=1,3
                ival =ival+1; atmp%fpr(ic3,atmp%im)= recv_buff_dbl(ival,i_at,ind_recv)
             end do
             if (allocated(atmp%rho)) then
                ival =ival+1; atmp%rho(atmp%im)= recv_buff_dbl(ival,i_at,ind_recv)
             end if
             if (allocated(atmp%fpg)) then
                do ic3=1,3
                   ival =ival+1; atmp%fpg(ic3,atmp%im)= recv_buff_dbl(ival,i_at,ind_recv)
                end do
             end if

          end select

       enddo
    enddo   ! fin de boucle sur les processeurs voisins

  end subroutine reception_nouveaux_atomes


  !------------------------------------------------------------------------!
  ! Procedure en charge de l'elimination des atomes locaux qui sont 
  ! maintenant fantomes (ils ont ete prealablement envoyes aux bons 
  ! processeurs par la routine envoi_atomes_fantomes)
  ! On elimine aussi les atomes deja fantomes pour faire la place
  ! a la mise a jour des atomes fantomes realisees dans 
  ! reception_atomes_fantomes

  subroutine elimine_atomes_fantomes(psc)

    USE T_kind_param_m, ONLY:  double


    implicit none
    type(para_space_config)::psc
    integer :: i_at
    integer :: i_new
    integer :: j_at
    integer :: koo
    integer :: nb_at_a_eliminer                 ! Nbre d'atomes a eliminer
    integer :: pt_at_elimine                    ! Pointeur sur le dernier atome elimine
    integer, allocatable :: at_a_eliminer(:)    ! Liste des atomes a eliminer


    allocate(at_a_eliminer(atmp%im))
    at_a_eliminer = 0

    ! initialisation de la liste des atomes a elatmp%iminer
    nb_at_a_eliminer = 0
    do i_at = 1, atmp%im
       koo = atmp%ielat(i_at)
       if (celmp%proc_cell(koo).ne.myidsp) then
          nb_at_a_eliminer = nb_at_a_eliminer + 1
          at_a_eliminer(nb_at_a_eliminer) = i_at
       endif
    enddo

    if (nb_at_a_eliminer.ne.0) then
       i_new = 1
       pt_at_elimine = 1
       ! On boucle sur tous les atomes
       do i_at = 1, atmp%im

          ! Si il s'agit d'un atome a eliminer
          if (i_at.eq.at_a_eliminer(pt_at_elimine)) then

             ! On met a jour les caracteristiques de la cellule correspondante
             koo = atmp%ielat(i_at)
             do j_at = 1, celmp%nato(koo)
                if ( celmp%atincel(j_at,koo).eq.i_at ) then
                   if (j_at.eq.celmp%nato(koo)) then
                      celmp%atincel(j_at,koo) = 0
                   else
                      celmp%atincel(j_at:celmp%nato(koo)-1,koo) = celmp%atincel(j_at+1:celmp%nato(koo),koo)
                   endif
                endif
             enddo
             celmp%nato(koo) = celmp%nato(koo) - 1

          endif

          ! Si les deux pointeurs ne sont pas au meme point, on deplace l'atome courant
          if ( i_new.ne.i_at ) then
             atmp%xp(:,i_new)  = atmp%xp(:,i_at)
             atmp%lgul(i_new)  = atmp%lgul(i_at)
             !LPARAFULLSEND

             !             ax(:,i_new)  = ax(:,i_at)
!!$	     if (lsuivinonpbc)  xpnonpbc(:,i_new)  = xpnonpbc(:,i_at)
!!$	     if (lsuivinonpbc)  tmpsuivi(:,i_new)  = tmpsuivi(:,i_at)
!!$	     if (lsuivinonpbc)  axnonpbc(:,i_new)  = axnonpbc(:,i_at)
             atmp%fp(:,i_new)  = atmp%fp(:,i_at)
             atmp%ityp(i_new)        = atmp%ityp(i_at)
             atmp%ielat(i_new)       = atmp%ielat(i_at)
             !             iwmax(i_new)       = iwmax(i_at)
             atmp%num_at_glob(i_new) = atmp%num_at_glob(i_at)

             select type (atmp)
             class is (atom_config_d)
                atmp%vp(:,i_new)  = atmp%vp(:,i_at)
             end select
             select type (atmp)
             class is (atom_config_e)
                if (atmp%lxpp)        atmp%xpp(:,i_new)  = atmp%xpp(:,i_at)
                if (atmp%lprteat)    atmp%eat(i_new)  = atmp%eat(i_at)
                if (atmp%lsigat)        atmp%sigat(:,:,i_new)  = atmp%sigat(:,:,i_at)
                if (atmp%llangevin)        atmp%glangv(:,i_new)  = atmp%glangv(:,i_at)
                if (atmp%lax)        atmp%ax(:,i_new)  = atmp%ax(:,i_at)
             end select
             select type (atmp)
             class is (atom_config_arps)
                atmp%fpr(:,i_new) = atmp%fpr(:,i_at)
                atmp%mov(i_new) = atmp%mov(i_at)
                if (allocated(atmp%rho)) atmp%rho(i_new)  = atmp%rho(i_at)
                if (allocated(atmp%fpg))atmp%fpg(:,i_new)  = atmp%fpg(:,i_at)
             end select




             ! On met aussi a jour le numero local de l'atome dans la liste de la cellule
             do j_at=1,celmp%nato(atmp%ielat(i_at))
                if (celmp%atincel(j_at,atmp%ielat(i_at)).eq.i_at) then
                   celmp%atincel(j_at,atmp%ielat(i_at))=i_new
                   exit
                endif
             enddo

          endif

          ! Mise a jour des pointeurs d'avancement
          if (i_at.eq.at_a_eliminer(pt_at_elimine)) then
             pt_at_elimine = pt_at_elimine + 1
          else
             i_new = i_new + 1
          endif

       enddo

       atmp%im = atmp%im - nb_at_a_eliminer
       !       imd = imd - nb_at_a_eliminer
       !       imf = imf - nb_at_a_eliminer
       !       imana = imana - nb_at_a_eliminer

    endif

    deallocate(at_a_eliminer)

    ! On verifie qu'il n'y a plus d'atomes a l'exterieur du domaine local
    do koo=1,celmp%noxyz
       if (celmp%proc_cell(koo).ne.myidsp .and. celmp%nato(koo).ne.0) print *,'ERREUR !!!',&
            myidsp,'possede encore',celmp%nato(koo),'at. dans la cellule',koo
    enddo

  end subroutine elimine_atomes_fantomes


  !------------------------------------------------------------------------!
  ! Procedure en charge de l'envoi des atomes frontieres du processeur 
  ! courant vers les processeurs voisins concernes

  subroutine envoi_atomes_frontieres(psc)

    USE T_kind_param_m, ONLY:  double

    implicit none
    type(para_space_config)::psc
    integer :: nproc_voisin
    integer :: ncell_front
    integer :: procv
    integer :: koo
    integer :: nb_at,nb_at_max,nb_at_max_tot
    integer :: n_at
    integer :: i_at

    ! Boucle a vide pour determiner au mieux la taille du buffer d'envoi 
    nb_at_max = 0
    do nproc_voisin = 1, psc%nbr_proc_voisin
       nb_at=0
       do ncell_front = 1, psc%nbr_cell_frontiere(nproc_voisin)
          nb_at = nb_at + celmp%nato(psc%cell_frontiere(nproc_voisin,ncell_front))
       enddo
       nb_at_max = max(nb_at_max, nb_at)
    enddo
    call MPI_ALLREDUCE(nb_at_max,nb_at_max_tot,1,MPI_INTEGER,MPI_MAX,MPI_COMM_space,ierr)
    nb_at_max=nb_at_max_tot
    ! petite manip pour eviter le cas nb_at_max=0
    nb_at_max = max(nb_at_max,1)

    ! Allocation des buffers
    nb_var_int = 4
    nb_var_lgc = 1
    !LPARAFULLSEND
    !    nb_var_dbl = 15
!!$   if  (lsuivinonpbc) then
!!$    nb_var_dbl = 18
!!$    else  
!!$    nb_var_dbl = 9
!!$   end if 

    nb_var_dbl = 6  !! xp,fp
    select type (atmp)
    class is (atom_config_d)
       nb_var_dbl=nb_var_dbl+3 !vp 
    end select
    select type (atmp)

    class is (atom_config_e)
       !       nb_var_dbl=nb_var_dbl+6 !vp xpp
       if (atmp%lxpp)        nb_var_dbl=nb_var_dbl+3 !eat
       if (atmp%lprteat)        nb_var_dbl=nb_var_dbl+1 !eat
       if (atmp%lsigat)        nb_var_dbl=nb_var_dbl+9 !eat
       if (atmp%llangevin)        nb_var_dbl=nb_var_dbl+3 !eat
       if (atmp%lax)        nb_var_dbl=nb_var_dbl+3 !eat
       !       write(6,*)'FLAGS', atmp%lprteat,atmp%lsigat,atmp%llangevin,atmp%lax
    end select
    select type (atmp)
    class is (atom_config_arps)
       nb_var_dbl=nb_var_dbl+3
       if (allocated(atmp%rho))        nb_var_dbl=nb_var_dbl+1
       if (allocated(atmp%fpg))        nb_var_dbl=nb_var_dbl+3
       !fprrho!vp xpp
       nb_var_int =     nb_var_int +1 !mov


    end select
    allocate(send_nb_val(psc%nbr_proc_voisin))
    allocate(send_buff_int(nb_var_int,nb_at_max,psc%nbr_proc_voisin))
    allocate(send_buff_lgc(nb_var_lgc,nb_at_max,psc%nbr_proc_voisin))
    allocate(send_buff_dbl(nb_var_dbl,nb_at_max,psc%nbr_proc_voisin))
    allocate(send_rqst(psc%nbr_proc_voisin,4))
    allocate(recv_nb_val(psc%nbr_proc_voisin))
    allocate(recv_buff_int(nb_var_int,nb_at_max,psc%nbr_proc_voisin))
    allocate(recv_buff_lgc(nb_var_lgc,nb_at_max,psc%nbr_proc_voisin))
    allocate(recv_buff_dbl(nb_var_dbl,nb_at_max,psc%nbr_proc_voisin))
    allocate(recv_rqst(psc%nbr_proc_voisin,4))

    ! Preparation des receptions
    do nproc_voisin = 1, psc%nbr_proc_voisin
       procv = psc%proc_voisin(nproc_voisin)
       call MPI_IRECV(recv_nb_val(nproc_voisin), 1, MPI_INTEGER, procv, 2001, MPI_COMM_space, recv_rqst(nproc_voisin,1), ierr)
       call MPI_IRECV(recv_buff_int(1,1,nproc_voisin), nb_var_int*nb_at_max, MPI_INTEGER, procv, 2002, &
            MPI_COMM_space, recv_rqst(nproc_voisin,2), ierr)
       call MPI_IRECV(recv_buff_dbl(1,1,nproc_voisin), nb_var_dbl*nb_at_max, NDM_MPI_REAL_DOUBLE, procv, 2003, &
            MPI_COMM_space, recv_rqst(nproc_voisin,3), ierr)
       call MPI_IRECV(recv_buff_lgc(1,1,nproc_voisin), nb_var_lgc*nb_at_max, MPI_LOGICAL, procv, 2004, &
            MPI_COMM_space, recv_rqst(nproc_voisin,4), ierr)

    enddo


    ! On boucle sur les processeurs voisins
    do nproc_voisin = 1, psc%nbr_proc_voisin
       procv = psc%proc_voisin(nproc_voisin)

       send_nb_val(nproc_voisin) = 0

       ! On boucle sur les cellules frontieres associees au processeur
       do ncell_front = 1, psc%nbr_cell_frontiere(nproc_voisin)
          koo = psc%cell_frontiere(nproc_voisin,ncell_front)

          ! On copie le contenu de la cellule dans le buffer d'envoi
          do n_at = 1, celmp%nato(koo)
             i_at = celmp%atincel(n_at,koo)

             ! On complete le buffer
             send_nb_val(nproc_voisin) = send_nb_val(nproc_voisin) + 1

             send_buff_lgc(1,send_nb_val(nproc_voisin),nproc_voisin) = atmp%lgul(i_at)

             send_buff_int(1,send_nb_val(nproc_voisin),nproc_voisin) = atmp%ityp(i_at)
             send_buff_int(2,send_nb_val(nproc_voisin),nproc_voisin) = atmp%ielat(i_at)
             !             send_buff_int(3,send_nb_val(nproc_voisin),nproc_voisin) = iwmax(i_at)
             send_buff_int(3,send_nb_val(nproc_voisin),nproc_voisin) = atmp%num_at_glob(i_at)
             select type (atmp)
             class is (atom_config_arps)
                send_buff_int(4,send_nb_val(nproc_voisin),nproc_voisin) = atmp%mov(i_at)
             end select

             send_buff_dbl(1,send_nb_val(nproc_voisin),nproc_voisin) = atmp%xp(1,i_at)
             send_buff_dbl(2,send_nb_val(nproc_voisin),nproc_voisin) = atmp%xp(2,i_at)
             send_buff_dbl(3,send_nb_val(nproc_voisin),nproc_voisin) = atmp%xp(3,i_at)

             send_buff_dbl(4,send_nb_val(nproc_voisin),nproc_voisin) = atmp%fp(1,i_at)
             send_buff_dbl(5,send_nb_val(nproc_voisin),nproc_voisin) = atmp%fp(2,i_at)
             send_buff_dbl(6,send_nb_val(nproc_voisin),nproc_voisin) = atmp%fp(3,i_at)

             select type (atmp)
             class is (atom_config_d)
                send_buff_dbl(7,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(1,i_at)
                send_buff_dbl(8,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(2,i_at)
                send_buff_dbl(9,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(3,i_at)
             end select
             ival=9
             select type (atmp)
             class is (atom_config_e)
!!$                send_buff_dbl(7,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(1,i_at)
!!$                send_buff_dbl(8,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(2,i_at)
!!$                send_buff_dbl(9,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(3,i_at)
                if (atmp%lxpp) then
                   do ic3=1,3
                      ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%xpp(ic3,i_at)
                   end do
                end if

                if (atmp%lprteat)then
                   ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%eat(i_at)
                end if
                if (atmp%lsigat)then
                   do ic1=1,3
                      do ic2=1,3
                         ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%sigat(ic1,ic2,i_at)
                      end do
                   end do
                end if
                if (atmp%llangevin) then
                   do ic3=1,3
                      ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%glangv(ic3,i_at)
                   end do
                end if
                if (atmp%lax) then
                   do ic3=1,3
                      ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%ax(ic3,i_at)
                   end do
                end if
             end select
             select type (atmp)
             class is (atom_config_arps)
                do ic3=1,3
                   ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%fpr(ic3,i_at)
                end do
                if (allocated(atmp%rho)) then
                   ival =ival+1; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%rho(i_at)
                end if
                if (allocated(atmp%fpg)) then
                   do ic3=1,3
                      ival =ival+1;; send_buff_dbl(ival,send_nb_val(nproc_voisin),nproc_voisin) = atmp%fpg(ic3,i_at)
                   end do
                end if

             end select
          enddo

       enddo   ! fin de boucle sur les cellules

       ! On envoit les buffers vers le processeur
       call MPI_ISEND(send_nb_val(nproc_voisin),      1,   MPI_INTEGER,        &
            procv,2001,MPI_COMM_space,send_rqst(nproc_voisin,1),ierr)

       call MPI_ISEND(send_buff_int(1,1,nproc_voisin),nb_var_int*nb_at_max,MPI_INTEGER,        &
            procv,2002,MPI_COMM_space,send_rqst(nproc_voisin,2),ierr)

       call MPI_ISEND(send_buff_dbl(1,1,nproc_voisin),nb_var_dbl*nb_at_max,NDM_MPI_REAL_DOUBLE,&
            procv,2003,MPI_COMM_space,send_rqst(nproc_voisin,3),ierr)
       call MPI_ISEND(send_buff_lgc(1,1,nproc_voisin),nb_var_lgc*nb_at_max,MPI_LOGICAL,        &
            procv,2004,MPI_COMM_space,send_rqst(nproc_voisin,4),ierr)


    enddo   ! fin de boucle sur les processeurs voisins


  end subroutine envoi_atomes_frontieres

  !------------------------------------------------------------------------!
  ! Procedure en charge d'attendre la fin des envois des atomes (frontieres
  ! ou fantomes) et la liberation des buffers d'envoi

  subroutine finalisation_envoi_atomes(ne,psc,nb_var_int,nb_var_lgc,nb_var_dbl)

    USE T_kind_param_m, ONLY:  double

    implicit none
    type(para_space_config)::psc
    integer,intent(in)::ne
    integer,optional::nb_var_int,nb_var_lgc,nb_var_dbl
    integer :: nproc_voisin

    !    integer, allocatable :: send_status(:,:,:)


    !    allocate(send_status(MPI_STATUS_SIZE,psc%nbr_proc_voisin,4))

    ! Attente de finalisation des envois 
    if(present(nb_var_int)) then
       do nproc_voisin= 1, psc%nbr_proc_voisin
          call MPI_Wait( send_rqst(nproc_voisin,1), status, ierr )
          call MPI_Wait( send_rqst(nproc_voisin,2), status, ierr )
          if (nb_var_dbl.gT.0) call MPI_Wait( send_rqst(nproc_voisin,3), status, ierr )
          if (nb_var_lgc.gT.0) call MPI_Wait( send_rqst(nproc_voisin,4), status, ierr )
       end do
    else

       do nproc_voisin= 1, psc%nbr_proc_voisin
          call MPI_Wait( send_rqst(nproc_voisin,1), status, ierr )
          call MPI_Wait( send_rqst(nproc_voisin,2), status, ierr )
          call MPI_Wait( send_rqst(nproc_voisin,3), status, ierr )
          if (ne==4) call MPI_Wait( send_rqst(nproc_voisin,4), status, ierr )
       enddo
    end if
    ! Liberation des buffers

    !    deallocate(send_status)
    if(present(nb_var_int)) then
       if (nb_var_dbl.gT.0) then
          deallocate(send_buff_dbl)
          deallocate(recv_buff_dbl)
       end if
       if (nb_var_lgc.gT.0)then
          deallocate(send_buff_lgc)
          deallocate(recv_buff_lgc)
       end if
       deallocate(send_rqst)
       deallocate(send_nb_val)
       deallocate(send_buff_int)
       deallocate(recv_rqst)
       deallocate(recv_nb_val)
       deallocate(recv_buff_int)
    else
       deallocate(send_rqst)
       deallocate(send_nb_val)
       deallocate(send_buff_int)
       if (ne==4)deallocate(send_buff_lgc)
       deallocate(send_buff_dbl)

       deallocate(recv_rqst)
       deallocate(recv_nb_val)
       deallocate(recv_buff_int)
       if (ne==4)    deallocate(recv_buff_lgc)
       deallocate(recv_buff_dbl)
    end if
  end subroutine finalisation_envoi_atomes


  !------------------------------------------------------------------------!
  ! Procedure en charge de la reception des nouveaux atomes fantomes en 
  ! provenance des processeurs voisins

  subroutine reception_atomes_fantomes(psc) !seulement MAJ

    USE T_kind_param_m, ONLY:  double

    implicit none
    type(para_space_config)::psc
    integer :: proc_source
    integer :: i_at
    integer :: nproc_voisin
    integer :: pt_at_ftm
    integer :: ind_recv

    integer :: nb_at_max,nb_at

    ! On place le pointeur de stockage des atomes fantomes a la suite des 
    ! atomes locaux
    pt_at_ftm = atmp%im

    ! On boucle sur les processeurs voisins
    do nproc_voisin=1,psc%nbr_proc_voisin

       call MPI_WAITANY(psc%nbr_proc_voisin,recv_rqst(:,1),ind_recv,status,ierr)
       proc_source = psc%proc_voisin(ind_recv)

       ! Reception des nouveaux atomes issus de ce processeur voisin
       call MPI_WAIT(recv_rqst(ind_recv,2), status, ierr)
       call MPI_WAIT(recv_rqst(ind_recv,3), status, ierr)
       call MPI_WAIT(recv_rqst(ind_recv,4), status, ierr)
       !       write(6,*)'indice nbv',ind_recv,recv_nb_val(ind_recv)
       ! recopie des infos dans les tableaux locaux au niveau des atomes fantomes
       do i_at = 1, recv_nb_val(ind_recv)

          ! On ajoute un atome fantome a la liste
          pt_at_ftm = pt_at_ftm + 1

          ! mise a jour des variables entieres
          atmp%lgul(pt_at_ftm)  = recv_buff_lgc(1,i_at,ind_recv)

          atmp%ityp(pt_at_ftm)  = recv_buff_int(1,i_at,ind_recv)
          atmp%ielat(pt_at_ftm) = recv_buff_int(2,i_at,ind_recv)
          !          iwmax(pt_at_ftm) = recv_buff_int(3,i_at,ind_recv)
          atmp%num_at_glob(pt_at_ftm) = recv_buff_int(3,i_at,ind_recv)
          select type (atmp)
          class is (atom_config_arps)
             atmp%mov(pt_at_ftm) = recv_buff_int(4,i_at,ind_recv)
          end select

          ! mise a jour des donnees de la cellule correspondante
          celmp%nato(atmp%ielat(pt_at_ftm)) = celmp%nato(atmp%ielat(pt_at_ftm)) + 1
          celmp%atincel(celmp%nato(atmp%ielat(pt_at_ftm)),atmp%ielat(pt_at_ftm)) = pt_at_ftm

          ! mise a jour des variables reelles
          atmp%xp(1,pt_at_ftm) = recv_buff_dbl(1,i_at,ind_recv) 
          atmp%xp(2,pt_at_ftm) = recv_buff_dbl(2,i_at,ind_recv) 
          atmp%xp(3,pt_at_ftm) = recv_buff_dbl(3,i_at,ind_recv) 
          !          ax(1,pt_at_ftm) = recv_buff_dbl(4,i_at,ind_recv)
          !          ax(2,pt_at_ftm) = recv_buff_dbl(5,i_at,ind_recv)
          !          ax(3,pt_at_ftm) = recv_buff_dbl(6,i_at,ind_recv)
          atmp%fp(1,pt_at_ftm) = recv_buff_dbl(4,i_at,ind_recv) 
          atmp%fp(2,pt_at_ftm) = recv_buff_dbl(5,i_at,ind_recv) 
          atmp%fp(3,pt_at_ftm) = recv_buff_dbl(6,i_at,ind_recv)

          select type (atmp)
          class is (atom_config_d)
             atmp%vp(1,pt_at_ftm) = recv_buff_dbl(7,i_at,ind_recv) 
             atmp%vp(2,pt_at_ftm) = recv_buff_dbl(8,i_at,ind_recv) 
             atmp%vp(3,pt_at_ftm) = recv_buff_dbl(9,i_at,ind_recv)
          end select
          ival=9
          select type (atmp)
          class is (atom_config_e)
!!$             atmp%vp(1,pt_at_ftm) = recv_buff_dbl(7,i_at,ind_recv) 
!!$             atmp%vp(2,pt_at_ftm) = recv_buff_dbl(8,i_at,ind_recv) 
!!$             atmp%vp(3,pt_at_ftm) = recv_buff_dbl(9,i_at,ind_recv)
             if (atmp%lxpp) then
                do ic3=1,3
                   ival =ival+1;atmp%xpp(ic3,pt_at_ftm)= recv_buff_dbl(ival,i_at,ind_recv)
                end do
             end if

             if (atmp%lprteat)then
                ival =ival+1
                atmp%eat(pt_at_ftm)=recv_buff_dbl(ival,i_at,ind_recv)
             end if
             if (atmp%lsigat)then
                do ic1=1,3
                   do ic2=1,3
                      ival =ival+1
                      atmp%sigat(ic1,ic2,pt_at_ftm)=recv_buff_dbl(ival,i_at,ind_recv)
                   end do
                end do
             end if
             if (atmp%llangevin) then
                do ic3=1,3
                   ival =ival+1;atmp%glangv(ic3,pt_at_ftm)= recv_buff_dbl(ival,i_at,ind_recv)
                end do
             end if
             if (atmp%lax) then
                do ic3=1,3
                   ival =ival+1;atmp%ax(ic3,pt_at_ftm)= recv_buff_dbl(ival,i_at,ind_recv)
                end do
             end if
          end select
          select type (atmp)
          class is (atom_config_arps)
             do ic3=1,3
                ival =ival+1;atmp%fpr(ic3,pt_at_ftm)= recv_buff_dbl(ival,i_at,ind_recv)
             end do
             if (allocated(atmp%rho)) then
                ival =ival+1;atmp%rho(pt_at_ftm)= recv_buff_dbl(ival,i_at,ind_recv)
             end if
             if (allocated(atmp%fpg)) then
                do ic3=1,3
                   ival =ival+1;atmp%fpg(ic3,pt_at_ftm)= recv_buff_dbl(ival,i_at,ind_recv)
                end do
             end if
          end select


       enddo

    enddo   ! fin de boucle sur les processeurs voisins
    atmp%imf=pt_at_ftm
  end subroutine reception_atomes_fantomes

  !------------------------------------------------------------------------!
  ! Procedure dont le but est l'envoi des valeurs de tabdensity pour les 
  ! atomes frontieres

  subroutine envoi_tabdensity_frontieres(tabdensity,imm,num_at_glob,psc)

    USE T_kind_param_m, ONLY:  double

    implicit none
    type(para_space_config)::psc
    integer::imm
    integer :: nproc_voisin
    integer :: ncell_front
    integer :: procv
    integer :: koo
    integer :: nb_at_max_tot,nb_at_max,nb_at
    integer :: n_at
    integer :: i_at
    real(double) :: tabdensity(imm)
    integer,intent(in)::num_at_glob(imm)

    ! Boucle a vide pour determiner au mieux la taille du buffer d'envoi 
    nb_at_max = 0
    do nproc_voisin = 1, psc%nbr_proc_voisin
       nb_at=0
       do ncell_front = 1, psc%nbr_cell_frontiere(nproc_voisin)
          nb_at = nb_at + celmp%nato(psc%cell_frontiere(nproc_voisin,ncell_front))
       enddo
       nb_at_max = max(nb_at_max, nb_at)
    enddo

    call MPI_ALLREDUCE(nb_at_max,nb_at_max_tot,1,MPI_INTEGER,MPI_MAX,MPI_COMM_space,ierr)
    nb_at_max=nb_at_max_tot
    ! petite manip pour eviter le cas nb_at_max=0
    nb_at_max = max(nb_at_max,1)

    ! Allocation des buffers
    allocate(send_nb_val(psc%nbr_proc_voisin))
    allocate(send_buff_int(1,nb_at_max,psc%nbr_proc_voisin))
    allocate(send_buff_dbl(1,nb_at_max,psc%nbr_proc_voisin))
    allocate(send_rqst(psc%nbr_proc_voisin,3))
    allocate(recv_nb_val(psc%nbr_proc_voisin))
    allocate(recv_buff_int(1,nb_at_max,psc%nbr_proc_voisin))
    allocate(recv_buff_dbl(1,nb_at_max,psc%nbr_proc_voisin))
    allocate(recv_rqst(psc%nbr_proc_voisin,3))

    ! Preparation des receptions
    do nproc_voisin = 1, psc%nbr_proc_voisin
       procv = psc%proc_voisin(nproc_voisin)
       call MPI_IRECV(recv_nb_val(nproc_voisin), 1, MPI_INTEGER, procv, 3001, MPI_COMM_space, recv_rqst(nproc_voisin,1), ierr)
       call MPI_IRECV(recv_buff_int(1,1,nproc_voisin), nb_at_max, MPI_INTEGER, procv, 3002, &
            MPI_COMM_space, recv_rqst(nproc_voisin,2), ierr)
       call MPI_IRECV(recv_buff_dbl(1,1,nproc_voisin), nb_at_max, NDM_MPI_REAL_DOUBLE, procv, 3003, &
            MPI_COMM_space, recv_rqst(nproc_voisin,3), ierr)
    enddo

    !    write(6,*)'psc%nbr_proc_voisin',rang,psc%nbr_proc_voisin
    !    write(6,*)'proc_voisin',rang,proc_voisin(1:psc%nbr_proc_voisin)
    !    write(6,*)'nbr_cell_frontiere',rang,nbr_cell_frontiere(1)
    !    write(6,*)'cell_frontiere',rang,cell_frontiere(1,1)

    ! On boucle sur les processeurs voisins
    do nproc_voisin = 1, psc%nbr_proc_voisin
       procv = psc%proc_voisin(nproc_voisin)

       send_nb_val(nproc_voisin) = 0

       ! On boucle sur les cellules frontieres associees au processeur
       do ncell_front = 1, psc%nbr_cell_frontiere(nproc_voisin)
          koo = psc%cell_frontiere(nproc_voisin,ncell_front)
          ! On copie le contenu de la cellule dans le buffer d'envoi
          do n_at = 1, celmp%nato(koo)
             i_at = celmp%atincel(n_at,koo)


             ! On complete le buffer
             send_nb_val(nproc_voisin) = send_nb_val(nproc_voisin) + 1

             send_buff_int(1,send_nb_val(nproc_voisin),nproc_voisin) = atmp%num_at_glob(i_at)

             send_buff_dbl(1,send_nb_val(nproc_voisin),nproc_voisin) = tabdensity(i_at)

          enddo

       enddo   ! fin de boucle sur les cellules
       !stop
       ! On envoit les buffers vers le processeur
       call MPI_ISEND(send_nb_val(nproc_voisin),      1,   MPI_INTEGER,        &
            procv,3001,MPI_COMM_space,send_rqst(nproc_voisin,1),ierr)

       call MPI_ISEND(send_buff_int(1,1,nproc_voisin),nb_at_max,MPI_INTEGER,        &
            procv,3002,MPI_COMM_space,send_rqst(nproc_voisin,2),ierr)

       call MPI_ISEND(send_buff_dbl(1,1,nproc_voisin),nb_at_max,NDM_MPI_REAL_DOUBLE,&
            procv,3003,MPI_COMM_space,send_rqst(nproc_voisin,3),ierr)

    enddo   ! fin de boucle sur les processeurs voisins


  end subroutine envoi_tabdensity_frontieres

  !------------------------------------------------------------------------!
  ! Procedure en charge de la reception des tabdensity des atomes fantomes
  ! en provenance des processeurs voisins

  subroutine reception_tabdensity_fantomes(tabdensity,imm,num_at_glob,psc,im)

    USE T_kind_param_m, ONLY:  double
    implicit none
    type(para_space_config)::psc
    integer::imm,im
    integer,intent(in)::num_at_glob(imm)
    integer :: proc_source
    integer :: i_at
    integer :: nproc_voisin
    integer :: pt_at_ftm
    integer :: ind_recv
    integer :: ind_loc,ftm_at

    real(double) :: tabdensity(imm)

    ! On place le pointeur de stockage des atomes fantomes a la suite des 
    ! atomes locaux
    pt_at_ftm = im

    ! On boucle sur les processeurs voisins
    do nproc_voisin=1,psc%nbr_proc_voisin

       call MPI_WAITANY(psc%nbr_proc_voisin,recv_rqst(:,1),ind_recv,status,ierr)
       proc_source = psc%proc_voisin(ind_recv)

       ! Reception des tabdensity issus de ce processeur voisin
       call MPI_WAIT(recv_rqst(ind_recv,2), status, ierr)
       call MPI_WAIT(recv_rqst(ind_recv,3), status, ierr)

       ! recopie des infos dans les tableaux locaux au niveau des atomes fantomes
       do i_at = 1, recv_nb_val(ind_recv)

          ! On boucle pour trouver l'indice local de l'atome fantome courant
          ind_loc=-1
          do ftm_at=im+1,imm
             if (atmp%num_at_glob(ftm_at)==recv_buff_int(1,i_at,ind_recv)) then
                ind_loc=ftm_at
                exit
             endif
          enddo
          if (ind_loc==-1) then
             print *,myidsp,'!!!Pb!!! Reception du proc',proc_source,'d''un atome fantome inexistant'
             call MPI_FINALIZE(ierr)
             call arret_ndm 
             !call arret_ndm
          endif

          ! On affecte a cet atome fantome la valeur de tabdensity recue

          tabdensity(ind_loc) = recv_buff_dbl(1,i_at,ind_recv)

       enddo

    enddo   ! fin de boucle sur les processeurs voisins

  end subroutine reception_tabdensity_fantomes


  !------------------------------------------------------------------------!
  ! Procedure dont le but est l'envoi des valeurs de fp pour les 
  ! atomes frontieres afin qu'elles soient sommees sur les processeurs
  ! possedant les atomes

  subroutine envoi_fp_fantomes(psc,atcf,celcf)

    USE T_kind_param_m, ONLY:  double

    implicit none
    class(atom_config),intent(inout)::atcf
    type(cell_config),intent(in)::celcf

    type(para_space_config)::psc
    integer :: nproc_voisin
    integer :: ncell_ftm
    integer :: procv
    integer :: koo
    integer :: nb_at,nb_at_max_tot
    integer :: n_at,nb_at_max
    integer :: i_at

    ! Boucle a vide pour determiner au mieux la taille du buffer d'envoi 
    nb_at_max = 0
    do nproc_voisin = 1, psc%nbr_proc_voisin
       nb_at=0
       do ncell_ftm = 1, psc%nbr_cell_ftm
          if (celmp%proc_cell(psc%cell_ftm(ncell_ftm)).eq.psc%proc_voisin(nproc_voisin)) then
             nb_at = nb_at + celcf%nato(psc%cell_ftm(ncell_ftm))
          endif
       enddo
       nb_at_max = max(nb_at_max, nb_at)
    enddo
    call MPI_ALLREDUCE(nb_at_max,nb_at_max_tot,1,MPI_INTEGER,MPI_MAX,MPI_COMM_space,ierr)
    nb_at_max=nb_at_max_tot
    ! petite manip pour eviter le cas nb_at_max=0
    nb_at_max = max(nb_at_max,1)

    ! Allocation des buffers
    allocate(send_nb_val(psc%nbr_proc_voisin))
    allocate(send_buff_int(1,nb_at_max,psc%nbr_proc_voisin))
    allocate(send_buff_dbl(3,nb_at_max,psc%nbr_proc_voisin))
    allocate(send_rqst(psc%nbr_proc_voisin,3))
    allocate(recv_nb_val(psc%nbr_proc_voisin))
    allocate(recv_buff_int(1,nb_at_max,psc%nbr_proc_voisin))
    allocate(recv_buff_dbl(3,nb_at_max,psc%nbr_proc_voisin))
    allocate(recv_rqst(psc%nbr_proc_voisin,3))

    ! Preparation des receptions
    do nproc_voisin = 1, psc%nbr_proc_voisin
       procv = psc%proc_voisin(nproc_voisin)
       call MPI_IRECV(recv_nb_val(nproc_voisin), 1, MPI_INTEGER, procv, 4001, MPI_COMM_space, recv_rqst(nproc_voisin,1), ierr)
       call MPI_IRECV(recv_buff_int(1,1,nproc_voisin), nb_at_max, MPI_INTEGER, procv, 4002, &
            MPI_COMM_space, recv_rqst(nproc_voisin,2), ierr)
       call MPI_IRECV(recv_buff_dbl(1,1,nproc_voisin), nb_at_max*3, NDM_MPI_REAL_DOUBLE, procv, 4003, &
            MPI_COMM_space, recv_rqst(nproc_voisin,3), ierr)
    enddo


    ! On boucle sur les processeurs voisins
    do nproc_voisin = 1, psc%nbr_proc_voisin
       procv = psc%proc_voisin(nproc_voisin)

       send_nb_val(nproc_voisin) = 0

       ! On boucle sur les cellules fantomes susceptibles d'appartenir au processeur
       do ncell_ftm = 1, psc%nbr_cell_ftm
          koo = psc%cell_ftm(ncell_ftm)

          ! appartient-elle au processeur voisin courant?
          if (celmp%proc_cell(koo).eq.procv) then

             ! On copie le contenu de la cellule dans le buffer d'envoi
             do n_at = 1, celcf%nato(koo)
                i_at = celcf%atincel(n_at,koo)

                ! On complete le buffer
                send_nb_val(nproc_voisin) = send_nb_val(nproc_voisin) + 1

                send_buff_int(1,send_nb_val(nproc_voisin),nproc_voisin) = atcf%num_at_glob(i_at)
                send_buff_dbl(1:3,send_nb_val(nproc_voisin),nproc_voisin) = atcf%fp(1:3,i_at)

             enddo
          endif

       enddo   ! fin de boucle sur les cellules

       ! On envoit les buffers vers le processeur
       call MPI_ISEND(send_nb_val(nproc_voisin),      1,   MPI_INTEGER,        &
            procv,4001,MPI_COMM_space,send_rqst(nproc_voisin,1),ierr)

       call MPI_ISEND(send_buff_int(1,1,nproc_voisin),nb_at_max,MPI_INTEGER,        &
            procv,4002,MPI_COMM_space,send_rqst(nproc_voisin,2),ierr)

       call MPI_ISEND(send_buff_dbl(1,1,nproc_voisin),nb_at_max*3,NDM_MPI_REAL_DOUBLE,&
            procv,4003,MPI_COMM_space,send_rqst(nproc_voisin,3),ierr)

    enddo   ! fin de boucle sur les processeurs voisins


  end subroutine envoi_fp_fantomes

  !------------------------------------------------------------------------!
  ! Procedure en charge de la reception des fp des atomes fantomes en
  ! provenance des processeurs voisins pour etre sommees en local

  subroutine reception_fp_frontieres(psc,atcf,celcf)

    USE T_kind_param_m, ONLY:  double
    implicit none
    class(atom_config),intent(inout)::atcf
    type(cell_config),intent(in)::celcf
    type(para_space_config)::psc

    integer :: proc_source
    integer :: i_at
    integer :: nproc_voisin
    integer :: ind_recv
    integer :: nb_at_max,koo
    integer :: ind_loc,i_at_loc,ind_glob


    ! On boucle sur les processeurs voisins
    do nproc_voisin=1,psc%nbr_proc_voisin

       call MPI_WAITANY(psc%nbr_proc_voisin,recv_rqst(:,1),ind_recv,status,ierr)
       proc_source = psc%proc_voisin(ind_recv)

       ! Reception des tabdensity issus de ce processeur voisin
       call MPI_WAIT(recv_rqst(ind_recv,2), status, ierr)
       call MPI_WAIT(recv_rqst(ind_recv,3), status, ierr)

       ! recopie des infos dans les tableaux locaux
       do i_at = 1, recv_nb_val(ind_recv)

          ind_glob = recv_buff_int(1,i_at,ind_recv)

          ! On boucle pour trouver l'indice local de l'atome fantome courant
          ind_loc=-1
          do i_at_loc=1,atcf%im
             if (atcf%num_at_glob(i_at_loc)==ind_glob) then
                ind_loc=i_at_loc
                exit
             endif
          enddo
          if (ind_loc==-1) then
             print *,myidsp,'!!!Pb!!! Reception du proc',proc_source,'d''un atome non local'
             call MPI_FINALIZE(ierr)
             call arret_ndm 
             call arret_ndm
             !call arret_ndm
          endif

          ! On ajoute a cet atome local la valeur de fp recue

          atcf%fp(1:3,ind_loc) = atcf%fp(1:3,ind_loc) + recv_buff_dbl(1:3,i_at,ind_recv)

       enddo

    enddo   ! fin de boucle sur les processeurs voisins

  end subroutine reception_fp_frontieres


  subroutine envoi_atomes_frontieres_part(psc,caracT,nb_var_int,nb_var_lgc,nb_var_dbl)

    USE T_kind_param_m, ONLY:  double

    implicit none
    type(para_space_config)::psc
    character(len=*),intent(in)::caracT
    integer::nb_var_int,nb_var_lgc,nb_var_dbl
    integer :: nproc_voisin
    integer :: ncell_front
    integer :: procv
    integer :: koo
    integer :: nb_at,nb_at_max,nb_at_max_tot
    integer :: n_at
    integer :: i_at
    integer::nvi,nvl,nvr
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   
    ! v=vp,r=xpp
    ! u=eat,g=glangv;a=ax;s=sigat; m=mov(arps)   

    ! Boucle a vide pour determiner au mieux la taille du buffer d'envoi 
    nb_at_max = 0
    do nproc_voisin = 1, psc%nbr_proc_voisin
       nb_at=0
       do ncell_front = 1, psc%nbr_cell_frontiere(nproc_voisin)
          nb_at = nb_at + celmp%nato(psc%cell_frontiere(nproc_voisin,ncell_front))
       enddo
       nb_at_max = max(nb_at_max, nb_at)
    enddo
    call MPI_ALLREDUCE(nb_at_max,nb_at_max_tot,1,MPI_INTEGER,MPI_MAX,MPI_COMM_space,ierr)
    nb_at_max=nb_at_max_tot
    ! petite manip pour eviter le cas nb_at_max=0
    nb_at_max = max(nb_at_max,1)

    ! Allocation des buffers
    nb_var_int = 0
    nb_var_lgc = 0
    nb_var_dbl = 0  

    !    if(scan('n',caract).ne.0)  then
    nb_var_int=nb_var_int+1
    !    end if
    if(scan('i',caracT).ne.0) then
       nb_var_int=nb_var_int+1
    end if
    if(scan('e',caracT).ne.0) then
       nb_var_int=nb_var_int+1
    end if
    if(scan('p',caracT).ne.0) then
       nb_var_int=nb_var_int+1
    end if
    if(scan('l',caracT).ne.0)  then
       nb_var_lgc=nb_var_lgc+1
    end if
    if(scan('x',caracT).ne.0)   then
       nb_var_dbl=nb_var_dbl+3
    end if
    if(scan('f',caracT).ne.0) then
       nb_var_dbl=nb_var_dbl+3
    end if


    select type (atmp)
    class is  (atom_config_d)
       if(scan('v',caracT).ne.0) then
          nb_var_dbl=nb_var_dbl+3
       end if
       if(scan('r',caracT).ne.0) then
          nb_var_dbl=nb_var_dbl+3
       end if
    end select
    select type (atmp)
    class is  (atom_config_e)
       if (atmp%lprteat)then
          if(scan('u',caracT).ne.0) then
             nb_var_dbl=nb_var_dbl+1
          end if
       end if
       if (atmp%llangevin)then
          if(scan('g',caracT).ne.0) then
             nb_var_dbl=nb_var_dbl+3
          end if
       end if
       if (atmp%lax)then
          if(scan('a',caracT).ne.0) then
             nb_var_dbl=nb_var_dbl+3
          end if
       end if
       if (atmp%lsigat)then
          if(scan('s',caracT).ne.0) then
             nb_var_dbl=nb_var_dbl+9
          end if
       end if
    end select
    select type (atmp)
    class is  (atom_config_arps)
       if(scan('m',caracT).ne.0) then
          nb_var_int=nb_var_int+1
       end if
    end select


    !    write(6,*)'nb_var_dbl2',nb_var_dbl
    allocate(send_nb_val(psc%nbr_proc_voisin))
    allocate(send_buff_int(nb_var_int,nb_at_max,psc%nbr_proc_voisin))
    if(nb_var_lgc.gt.0)    allocate(send_buff_lgc(nb_var_lgc,nb_at_max,psc%nbr_proc_voisin))
    if(nb_var_dbl.gt.0)    allocate(send_buff_dbl(nb_var_dbl,nb_at_max,psc%nbr_proc_voisin))
    allocate(send_rqst(psc%nbr_proc_voisin,4))
    allocate(recv_nb_val(psc%nbr_proc_voisin))
    allocate(recv_buff_int(nb_var_int,nb_at_max,psc%nbr_proc_voisin))
    if(nb_var_lgc.gt.0)    allocate(recv_buff_lgc(nb_var_lgc,nb_at_max,psc%nbr_proc_voisin))
    if(nb_var_dbl.gt.0)    allocate(recv_buff_dbl(nb_var_dbl,nb_at_max,psc%nbr_proc_voisin))
    allocate(recv_rqst(psc%nbr_proc_voisin,4))

    ! Preparation des receptions
    do nproc_voisin = 1, psc%nbr_proc_voisin
       procv = psc%proc_voisin(nproc_voisin)
       call MPI_IRECV(recv_nb_val(nproc_voisin), 1, MPI_INTEGER, procv, 2001, MPI_COMM_space, recv_rqst(nproc_voisin,1), ierr)
       call MPI_IRECV(recv_buff_int(1,1,nproc_voisin), nb_var_int*nb_at_max, MPI_INTEGER, procv, 2002, &
            MPI_COMM_space, recv_rqst(nproc_voisin,2), ierr)
       if(nb_var_dbl.gt.0) &
            &call MPI_IRECV(recv_buff_dbl(1,1,nproc_voisin), nb_var_dbl*nb_at_max, NDM_MPI_REAL_DOUBLE, procv, 2003, &
            MPI_COMM_space, recv_rqst(nproc_voisin,3), ierr)
       if(nb_var_lgc.gt.0)&
            &call MPI_IRECV(recv_buff_lgc(1,1,nproc_voisin), nb_var_lgc*nb_at_max, MPI_LOGICAL, procv, 2004, &
            MPI_COMM_space, recv_rqst(nproc_voisin,4), ierr)

    enddo


    ! On boucle sur les processeurs voisins
    do nproc_voisin = 1, psc%nbr_proc_voisin
       procv = psc%proc_voisin(nproc_voisin)

       send_nb_val(nproc_voisin) = 0

       ! On boucle sur les cellules frontieres associees au processeur
       do ncell_front = 1, psc%nbr_cell_frontiere(nproc_voisin)
          koo = psc%cell_frontiere(nproc_voisin,ncell_front)

          ! On copie le contenu de la cellule dans le buffer d'envoi
          do n_at = 1, celmp%nato(koo)
             i_at = celmp%atincel(n_at,koo)
             send_nb_val(nproc_voisin) = send_nb_val(nproc_voisin) + 1

             nvi=0;nvl=0;nvr=0

             !             if(scan('n',caracT).ne.0)  then
             nvi=nvi+1
             send_buff_int(nvi,send_nb_val(nproc_voisin),nproc_voisin) = atmp%num_at_glob(i_at)
             !             end if
             if(scan('i',caracT).ne.0) then
                nvi=nvi+1
                send_buff_int(nvi,send_nb_val(nproc_voisin),nproc_voisin) = atmp%ityp(i_at)
             end if
             if(scan('e',caracT).ne.0) then
                nvi=nvi+1
                send_buff_int(nvi,send_nb_val(nproc_voisin),nproc_voisin) = atmp%ielat(i_at)
             end if
             if(scan('p',caracT).ne.0) then
                nvi=nvi+1
                send_buff_int(nvi,send_nb_val(nproc_voisin),nproc_voisin) = atmp%proc_at(i_at)
             end if
             if(scan('l',caracT).ne.0)  then
                nvl=nvl+1
                send_buff_lgc(nvl,send_nb_val(nproc_voisin),nproc_voisin) = atmp%lgul(i_at)
             end if
             if(scan('x',caracT).ne.0)   then
                nvr=nvr+1
                send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%xp(1,i_at)
                nvr=nvr+1
                send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%xp(2,i_at)
                nvr=nvr+1
                send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%xp(3,i_at)
             end if
             if(scan('f',caracT).ne.0) then
                nvr=nvr+1
                send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%fp(1,i_at)
                nvr=nvr+1
                send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%fp(2,i_at)
                nvr=nvr+1
                send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%fp(3,i_at)
             end if

             select type (atmp)
             class is  (atom_config_d)
                if(scan('v',caracT).ne.0) then
                   nvr=nvr+1
                   send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(1,i_at)
                   nvr=nvr+1
                   send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(2,i_at)
                   nvr=nvr+1
                   send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%vp(3,i_at)
                end if
             end select


             select type (atmp)
             class is  (atom_config_e)
                if (atmp%lxpp) then
                   if(scan('r',caracT).ne.0) then
                      nvr=nvr+1
                      send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%xpp(1,i_at)
                      nvr=nvr+1
                      send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%xpp(2,i_at)
                      nvr=nvr+1
                      send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%xpp(3,i_at)

                   end if
                end if
                if (atmp%lprteat)then
                   if(scan('u',caracT).ne.0) then
                      nvr=nvr+1
                      send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%eat(i_at)
                   end if
                end if
                if (atmp%llangevin)then
                   if(scan('g',caracT).ne.0) then
                      do ic3=1,3
                         nvr=nvr+1; send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%glangv(ic3,i_at)
                      end do
                   end if
                end if
                if (atmp%lax)then
                   if(scan('a',caracT).ne.0) then
                      do ic3=1,3
                         nvr=nvr+1; send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%ax(ic3,i_at)
                      end do
                   end if
                end if
                if (atmp%lsigat)then
                   if(scan('s',caracT).ne.0) then
                      do ic1=1,3
                         do ic2=1,3
                            nvr=nvr+1; send_buff_dbl(nvr,send_nb_val(nproc_voisin),nproc_voisin) = atmp%sigat(ic1,ic2,i_at)
                         end do
                      end do

                   end if
                end if
             end select
             select type (atmp)
             class is  (atom_config_arps)
                if(scan('m',caracT).ne.0) then
                   nvi=nvi+1
                   send_buff_int(nvi,send_nb_val(nproc_voisin),nproc_voisin) = atmp%mov(i_at)
                end if
             end select

             if (nvi.ne.nb_var_int) then
                write(6,*)'pb val int part',rang,nvi,nb_var_int
                call arret_ndm
             end if
             if (nvl.ne.nb_var_lgc) then
                write(6,*)'pb val log part',rang,nvl,nb_var_lgc
                call arret_ndm
             end if
             if (nvr.ne.nb_var_dbl) then
                write(6,*)'pb val real part',rang,nvr,nb_var_dbl
                call arret_ndm
             end if

          enddo

       enddo   ! fin de boucle sur les cellules

       ! On envoit les buffers vers le processeur
       call MPI_ISEND(send_nb_val(nproc_voisin),      1,   MPI_INTEGER,        &
            procv,2001,MPI_COMM_space,send_rqst(nproc_voisin,1),ierr)


       call MPI_ISEND(send_buff_int(1,1,nproc_voisin),nb_var_int*nb_at_max,MPI_INTEGER,        &
            procv,2002,MPI_COMM_space,send_rqst(nproc_voisin,2),ierr)

       if (nb_var_dbl.ne.0)    &
            &call MPI_ISEND(send_buff_dbl(1,1,nproc_voisin),nb_var_dbl*nb_at_max,NDM_MPI_REAL_DOUBLE,&
            procv,2003,MPI_COMM_space,send_rqst(nproc_voisin,3),ierr)
       if (nb_var_lgc.ne.0)   &
            &call MPI_ISEND(send_buff_lgc(1,1,nproc_voisin),nb_var_lgc*nb_at_max,MPI_LOGICAL,        &
            procv,2004,MPI_COMM_space,send_rqst(nproc_voisin,4),ierr)


    enddo   ! fin de boucle sur les processeurs voisins


  end subroutine envoi_atomes_frontieres_part

  !------------------------------------------------------------------------!
  ! Procedure en charge d'attendre la fin des envois des atomes (frontieres
  ! ou fantomes) et la liberation des buffers d'envoi



  !------------------------------------------------------------------------!
  ! Procedure en charge de la reception des nouveaux atomes fantomes en 
  ! provenance des processeurs voisins

  subroutine reception_atomes_fantomes_part(psc,caracT,nb_var_int,nb_var_lgc,nb_var_dbl) !seulement MAJ

    USE T_kind_param_m, ONLY:  double

    implicit none
    type(para_space_config)::psc
    character(len=*),intent(in)::caracT
    integer::nb_var_int,nb_var_lgc,nb_var_dbl
    integer :: proc_source
    integer :: i_at
    integer :: nproc_voisin
    integer :: pt_at_ftm
    integer :: ind_recv,natgt

    integer :: nb_at_max,nb_at,ic,ic1,ic2,iloc,nvi,nvl,nvr
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   
    ! v=vp,r=xpp
    ! u=eat,g=glangv;a=ax;s=sigat; m=mov(arps)   

    ! On place le pointeur de stockage des atomes fantomes a la suite des 
    ! atomes locaux
    pt_at_ftm = atmp%im

!!$    nb_var_int = 0
!!$    nb_var_lgc = 0
!!$    nb_var_dbl = 0  
!!$
!!$    !    if(scan('n',carac).ne.0)  then 
!!$    'n'
!!$    nb_var_int=nb_var_int+1
!!$    !    end if
!!$    if(scan('i',carac).ne.0) then
!!$       nb_var_int=nb_var_int+1
!!$    end if
!!$    if(scan('e',carac).ne.0) then
!!$       nb_var_int=nb_var_int+1
!!$    end if
!!$    if(scan('p',carac).ne.0) then
!!$       nb_var_int=nb_var_int+1
!!$    end if
!!$    if(scan('l',carac).ne.0)  then
!!$       nvb_var_lgc=nb_var_lgc+1
!!$    end if
!!$    if(scan('x',carac).ne.0)   then
!!$       nb_var_dbl=nb_var_dbl+3
!!$    end if
!!$    if(scan('f',carac).ne.0) then
!!$       nb_var_dbl=nb_var_dbl+3
!!$    end if
!!$
!!$
!!$    select type (atmp)
!!$    class is  (atom_config_d)
!!$       if(scan('v',carac).ne.0) then
!!$          nb_var_dbl=nb_var_dbl+3
!!$       end if
!!$       if(scan('r',carac).ne.0) then
!!$          nb_var_dbl=nb_var_dbl+3
!!$       end if
!!$    end select
!!$    select type (atmp)
!!$    class is  (atom_config_e)
!!$       if (atmp%lprteat)then
!!$          if(scan('u',carac).ne.0) then
!!$             nb_var_dbl=nb_var_dbl+1
!!$          end if
!!$       end if
!!$       if (atmp%llangevin)then
!!$          if(scan('g',carac).ne.0) then
!!$             nb_var_dbl=nb_var_dbl+3
!!$          end if
!!$       end if
!!$       if (atmp%lax)then
!!$          if(scan('a',carac).ne.0) then
!!$             nb_var_dbl=nb_var_dbl+3
!!$          end if
!!$       end if
!!$       if (atmp%lsigat)then
!!$          if(scan('s',carac).ne.0) then
!!$             nb_var_dbl=nb_var_dbl+9
!!$          end if
!!$       end if
!!$    end select
!!$    select type (atmp)
!!$    class is  (atom_config_arps)
!!$       if(scan('m',carac).ne.0) then
!!$          nb_var_int=nb_var_int+1
!!$       end if
!!$    end select
!!$

    ! On boucle sur les processeurs voisins
    do nproc_voisin=1,psc%nbr_proc_voisin

       call MPI_WAITANY(psc%nbr_proc_voisin,recv_rqst(:,1),ind_recv,status,ierr)
       proc_source = psc%proc_voisin(ind_recv)

       ! Reception des nouveaux atomes issus de ce processeur voisin
       call MPI_WAIT(recv_rqst(ind_recv,2), status, ierr)
       if(nb_var_dbl.gt.0)        call MPI_WAIT(recv_rqst(ind_recv,3), status, ierr)
       if(nb_var_lgc.gt.0)        call MPI_WAIT(recv_rqst(ind_recv,4), status, ierr)
       !       write(6,*)'indice nbv',ind_recv,recv_nb_val(ind_recv)
       ! recopie des infos dans les tableaux locaux au niveau des atomes fantomes
       loopfant:       do i_at = 1, recv_nb_val(ind_recv)


          nvi=0;nvl=0;nvr=0

          do iloc=atmp%im+1,atmp%imm
             natgt=recv_buff_int(1,i_at,ind_recv)
             if (atmp%num_at_glob(iloc)==recv_buff_int(1,i_at,ind_recv)) then
                !                if (natgt==1) write(6,*)'RANG iloc TROUVE', rang,natgt,iloc
                nvi=1
                if(scan('i',caracT).ne.0) then
                   nvi=nvi+1
                   atmp%ityp(iloc)  = recv_buff_int(nvi,i_at,ind_recv)
                end if
                if(scan('e',caracT).ne.0) then
                   nvi=nvi+1
                   atmp%ielat(iloc)  = recv_buff_int(nvi,i_at,ind_recv)
                end if
                if(scan('p',caracT).ne.0) then
                   nvi=nvi+1
                   atmp%proc_at(iloc)  = recv_buff_int(nvi,i_at,ind_recv)
                end if
                if(scan('l',caracT).ne.0)  then
                   nvl=nvl+1
                   atmp%lgul(iloc)  = recv_buff_lgc(nvl,i_at,ind_recv)

                end if
                if(scan('x',caracT).ne.0)   then
                   nvr=nvr+1
                   atmp%xp(1,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                   nvr=nvr+1
                   atmp%xp(2,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                   nvr=nvr+1
                   atmp%xp(3,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 


                end if
                if(scan('f',caracT).ne.0) then
                   nvr=nvr+1
                   atmp%fp(1,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                   nvr=nvr+1
                   atmp%fp(2,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                   nvr=nvr+1
                   atmp%fp(3,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                end if


                select type (atmp)
                class is  (atom_config_d)
                   if(scan('v',caracT).ne.0) then
                      nvr=nvr+1
                      atmp%vp(1,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                      nvr=nvr+1
                      atmp%vp(2,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                      nvr=nvr+1
                      atmp%vp(3,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                   end if
                end select
                select type (atmp)
                class is  (atom_config_e)
                   if (atmp%lxpp) then 
                      if(scan('r',caracT).ne.0) then
                         nvr=nvr+1
                         atmp%xpp(1,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                         nvr=nvr+1
                         atmp%xpp(2,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                         nvr=nvr+1
                         atmp%xpp(3,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                      end if
                   end if

                   if (atmp%lprteat)then
                      if(scan('u',caracT).ne.0) then
                         nvr=nvr+1
                         atmp%eat(iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                      end if
                   end if
                   if (atmp%llangevin)then
                      if(scan('g',caracT).ne.0) then
                         nvr=nvr+1
                         atmp%glangv(1,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                         nvr=nvr+1
                         atmp%glangv(2,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                         nvr=nvr+1
                         atmp%glangv(3,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 

                      end if
                   end if
                   if (atmp%lax)then
                      if(scan('a',caracT).ne.0) then
                         nvr=nvr+1
                         atmp%ax(1,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                         nvr=nvr+1
                         atmp%ax(2,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 
                         nvr=nvr+1
                         atmp%ax(3,iloc) = recv_buff_dbl(nvr,i_at,ind_recv) 

                      end if
                   end if
                   if (atmp%lsigat)then
                      if(scan('s',caracT).ne.0) then
                         do ic1=1,3
                            do ic2=1,3
                               nvr=nvr+1;  atmp%sigat(ic1,ic2,iloc)= recv_buff_dbl(nvr,i_at,ind_recv) 
                            end do
                         end do

                      end if
                   end if
                end select
                select type (atmp)
                class is  (atom_config_arps)
                   if(scan('m',caracT).ne.0) then
                      nvi=nvi+1;atmp%mov(iloc)= recv_buff_dbl(nvr,i_at,ind_recv) 
                   end if
                end select
                cycle loopfant

             end if
          end do
          write(6,*)'RANG iloc non trouvé', natgt
       end do loopfant

    end do
  end subroutine reception_atomes_fantomes_part

#endif

end module mod_para

