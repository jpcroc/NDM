module mod_para
#ifdef PARA
  use Tpara,only: NDM_MPI_REAL_DOUBLE
  
#endif
  use T_kind_param_m, ONLY:  double 

  use gen_com_m ,only:l2t,rang
  USE atomconfig,only:atom_config,atom_config_d,atom_config_e,ndm2config,config2ndm
  USE cellconfig,only:cell_config,ndm2cellconfig,cellconfig2ndm

  !  use mpi
  implicit none
  integer :: myid,nprocspace,nprocs 			! numero de process mis là pour être utilisé en sequentiel
#ifdef PARA

  include 'mpif.h'

  integer::MPI_COMM_space

  ! Module de declaration des variables MPI pour le code NDM

  !Entiers :



  integer :: ierr 			! erreur MPI
  integer,dimension(MPI_STATUS_SIZE):: status  ! statut de la communication
  integer:: grp_world
  integer :: nbr_proc_voisin            ! nbre de processeurs voisins du processeur courant

  !Tableaux specifiques :

  !Tableaux liés au decoupage :

  integer :: nbr_cell_max		!plus grand nombre de cellules sur tous les processeurs
  integer :: nbr_atom_max		!plus grand nombre d'atomes sur tous les processeurs
  integer :: nbr_cell_max_vois		!plus grand nombre de cellules voisines à un processeur
  integer :: nbr_atom_max_vois		!plus grand nombre d'atomes sur toutes les cellules voisines à un processeur
  integer :: resultat(4) 		!stocke le resultat du meilleur decoupage

  integer, allocatable :: res_cpu(:,:)   	!stocke le nombre de cellules de chaques decoupages pour le meilleur decoupage
  !res_cpu est initialisé dans decoup3D à (0:nprocspace-1,3)
  integer, allocatable :: coord_min(:,:)	!stocke la "coordonnée" de la premiere cellule du découpage selon x,y,z
  !initialisée à (0:nprocspace-1,3) dans decoup3D
  integer, allocatable :: coord_max(:,:)	!stocke la "coordonnée" de la derniere cellule du découpage selon x,y,z
  !initialisée à (0:nprocspace-1,3) dans decoup3D

  integer, allocatable :: proc_cell(:)          !proc_cell(i) : Numero du proc associe a la cellule i
  integer, allocatable :: proc_voisin(:)        ! liste des processeurs voisins du processeur courant
  integer, allocatable :: cell_frontiere(:,:)   ! (i,j) jeme cellule frontiere associee au ieme processeur voisin
  integer, allocatable :: nbr_cell_frontiere(:) ! nbre de cellules frontieres associees au ieme processeur voisin

  integer :: nbr_cell_ftm                       ! nbr de cellules fantomes du processeur courant
  integer, allocatable :: cell_ftm(:)           ! liste des cellules fantomes du processeur courant

  integer :: nb_var_int                              ! nbr de variables entieres a envoyer lors des echanges entre proc
  integer :: nb_var_dbl                              ! nbr de variables reelles a envoyer lors des echanges entre proc
  integer, allocatable :: send_nb_val(:)             ! buffer d'envoi du nombre de valeurs envoyees
  integer, allocatable :: send_buff_int(:,:,:)       ! buffer d'envoi des variables entieres
  real(double), allocatable :: send_buff_dbl(:,:,:)  ! buffer d'envoi des variables reelles

  integer, allocatable :: recv_nb_val(:)             ! buffer d'envoi du nombre de valeurs envoyees
  integer, allocatable :: recv_buff_int(:,:,:)       ! buffer d'envoi des variables entieres
  real(double), allocatable :: recv_buff_dbl(:,:,:)  ! buffer d'envoi des variables reelles

  integer, allocatable :: send_rqst(:,:)              ! tableau pour stocker les requetes en envoi
  integer, allocatable :: recv_rqst(:,:)              ! tableau pour stocker les requetes en reception

  real(double) :: temps_deb

  ! Variables pour faire des mesures de temps dans le code

  real(double) :: temps_init_deb, temps_init
  real(double) :: temps_initspeed_deb, temps_initspeed
  real(double) :: temps_dmloop_deb, temps_dmloop
  real(double) :: temps_input_deb, temps_input
  real(double) :: temps_config_deb, temps_config
  real(double) :: temps_para,temps_debpara,temps_finpara

  real(double), allocatable,dimension(:,:)::xp,vp,fp,xpp,ax,glangv
  real(double),allocatable::eat(:),sigat(:,:,:)
  integer,allocatable,dimension(:)::ityp,ielat,iwmax,num_at_glob,indi
  logical ::llangevin,lprteat,lsigat,ltbv,lax
  integer::im,imm,nvois
  integer:: nox,noy,noz,natperc,noxyz
  integer,allocatable::ncel(:,:),nato(:),atincel(:,:),deltadist(:,:,:)
  real(double)::celsize(3)
  logical::ltpcel
  real(double),allocatable::sigc(:,:,:),tempc(:)
  !---------------------------------------------------------!
  !               Routines spécifique à MPI                 !
  !---------------------------------------------------------!

contains

  !------------------------------------------------------------------------!
  ! Procedure pour la mise a jour des atomes frontieres et fantomes sur
  ! les processeurs. Elle prend en compte la nouvelle repartition dans 
  ! les cellules suite a l'appel a caltabt

  subroutine maj_atomes_frt_ftm(atcf,cellcf)

    USE T_kind_param_m, ONLY:  double
!    use tab_imm_m

    implicit none
    type(cell_config)::cellcf
    class(atom_config)::atcf
    integer::i
    ltbv=atcf%ltabvois ;
    lsigat=.false.; lprteat=.false. ; llangevin=.false.;lax=.false.
    select type (atcf)
    type is (atom_config_e)
       lsigat=atcf%lsigat;lprteat=atcf%lprteat; llangevin=atcf%llangevin;lax=atcf%lax
    end select
    call config2ndm(atcf,im,imm,xp,fp,ityp,ielat,num_at_glob,ltbv,iwmax,indi,vp,xpp,eat,sigat,ax,ldeall=.true.)
    call cellconfig2ndm (cellcf,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize,ltpcel,sigc,tempc,proc_cell)
    ! On envoit les atomes qui n'appartiennent plus au processeur courant
    call envoi_atomes_fantomes
    ! On recoit les nouveaux atomes locaux
    call reception_nouveaux_atomes
    ! On retire les atomes qui ne sont plus locaux
    call elimine_atomes_fantomes
    ! Finalisation de l'envoi des atomes pour liberer les buffers d'envoi
    call finalisation_envoi_atomes
    ! On envoit les atomes frontieres aux processeurs voisins
    call envoi_atomes_frontieres
    ! On receptionne les nouveaux atomes fantomes
    call reception_atomes_fantomes
    ! Finalisation de l'envoi des atomes pour liberer les buffers d'envoi
    call finalisation_envoi_atomes
    xpp=0;fp=0
    call ndm2config (atcf,im,imm,xp,fp,ityp,ielat,num_at_glob,ltbv,iwmax,indi,nvois,vp,xpp,ldeall=.true.)
    call ndm2cellconfig(cellcf,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize,proc_cell=proc_cell)


end subroutine maj_atomes_frt_ftm

  !------------------------------------------------------------------------!
  ! Procedure pour la mise a jour des valeurs tabdensity des atomes 
  ! fantomes sur les processeurs.

  subroutine maj_tabdensity_ftm(tabdensity,imm,natR,num_at_glob) !appelée dans calfoeamcel

    USE T_kind_param_m, ONLY:  double
!    use tab_imm_m

    implicit none
    integer::imm
    integer,intent(in)::num_at_glob(imm)
    real(double) :: tabdensity(imm)
    integer::natr(:)
    nato=natr 

    ! On envoit les atomes frontieres aux processeurs voisins
    call envoi_tabdensity_frontieres(tabdensity,imm,num_at_glob)

    ! On receptionne les nouveaux atomes fantomes
    call reception_tabdensity_fantomes(tabdensity,imm,num_at_glob)

    ! Finalisation de l'envoi pour liberer les buffers d'envoi (identique a l'envoi des atomes)
    call finalisation_envoi_atomes

  end subroutine maj_tabdensity_ftm

  !------------------------------------------------------------------------!
  ! Procedure pour la mise a jour des valeurs fp des atomes frontieres
  ! du processeur courant avec leurs contributions des processeurs voisins
 
  subroutine maj_fp_frt !appelée SEULEMENT dans force_tersoff_cel !

    USE T_kind_param_m, ONLY:  double
!    use tab_imm_m

    implicit none

 
    ! On envoit les atomes fantomes vers les processeurs voisins
    call envoi_fp_fantomes

    ! On receptionne les contributions des processeurs voisins
    call reception_fp_frontieres

    ! Finalisation de l'envoi pour liberer les buffers d'envoi (identique a l'envoi des atomes)
    call finalisation_envoi_atomes

  end subroutine maj_fp_frt

  !------------------------------------------------------------------------!
  ! Procedure dont le but est l'envoi des atomes qui sont sorti du domaine
  ! courant pour etre pris en charge par leur nouveau processeur

  subroutine envoi_atomes_fantomes ! seulement maj_atomes_frt_ftm

    USE T_kind_param_m, ONLY:  double
!    use tab_imm_m

    implicit none

    integer :: nb_at, nb_at_max, nb_at_max_tot
    integer :: nproc_voisin
    integer :: ncell_ftm
    integer :: procv,cellf,n_at,i_at

    ! Premier passage a vide pour allouer les buffers au plus juste

    nb_at_max=0
    ! Boucle sur les processeurs voisins
    do nproc_voisin=1,nbr_proc_voisin
       procv = proc_voisin(nproc_voisin)

       nb_at = 0

       ! Boucle sur les cellules fantomes
       do ncell_ftm=1,nbr_cell_ftm
          cellf = cell_ftm(ncell_ftm)

          ! Sommation des atomes de la cellule
          if (proc_cell(cellf)==procv) nb_at = nb_at + nato(cellf)

       enddo ! fin boucle sur les cellules

       ! calcul du max des atomes a envoyer
       nb_at_max = max(nb_at_max, nb_at)

    enddo ! fin boucle sur les processeurs voisins
    call MPI_ALLREDUCE(nb_at_max,nb_at_max_tot,1,MPI_INTEGER,MPI_MAX,MPI_COMM_space,ierr)
    nb_at_max=nb_at_max_tot
    ! petite manip pour eviter le cas nb_at_max=0
    nb_at_max = max(nb_at_max,1)

    ! Allocation des buffers
    nb_var_int = 3
!    nb_var_int = 4 avec iwmax aucun intérêt
!LPARAFULLSEND
!    nb_var_dbl = 15
!!$   if (lsuivinonpbc) then
!!$    nb_var_dbl = 18
!!$   else 
!!$    nb_var_dbl = 9
!!$ end if
!    nb_var_dbl = 9 avec ax
    nb_var_dbl = 6
  if ((llangevin.eqv..true.).or.(l2T.eqv..true.))then
     nb_var_dbl = nb_var_dbl+3
  end if


    allocate(send_nb_val(nbr_proc_voisin))
    allocate(send_buff_int(nb_var_int,nb_at_max,nbr_proc_voisin))
    allocate(send_buff_dbl(nb_var_dbl,nb_at_max,nbr_proc_voisin))
    allocate(send_rqst(nbr_proc_voisin,3))
    allocate(recv_nb_val(nbr_proc_voisin))
    allocate(recv_buff_int(nb_var_int,nb_at_max,nbr_proc_voisin))
    allocate(recv_buff_dbl(nb_var_dbl,nb_at_max,nbr_proc_voisin))
    allocate(recv_rqst(nbr_proc_voisin,3))

    ! Preparation des receptions
    do nproc_voisin= 1, nbr_proc_voisin
       procv = proc_voisin(nproc_voisin)
       call MPI_IRECV(recv_nb_val(nproc_voisin), 1, MPI_INTEGER, procv, 1001, MPI_COMM_space, recv_rqst(nproc_voisin,1), ierr)
       call MPI_IRECV(recv_buff_int(1,1,nproc_voisin), nb_var_int*nb_at_max, MPI_INTEGER, procv, 1002, &
            MPI_COMM_space, recv_rqst(nproc_voisin,2), ierr)
       call MPI_IRECV(recv_buff_dbl(1,1,nproc_voisin), nb_var_dbl*nb_at_max, NDM_MPI_REAL_DOUBLE, procv, 1003, &
            MPI_COMM_space, recv_rqst(nproc_voisin,3), ierr)
    enddo


    ! On boucle sur les processeurs voisins
    do nproc_voisin= 1, nbr_proc_voisin
       procv = proc_voisin(nproc_voisin)

       send_nb_val(nproc_voisin) = 0

       ! On boucle sur les cellules fantomes associees a ce processeur voisin
       do ncell_ftm= 1, nbr_cell_ftm
          cellf = cell_ftm(ncell_ftm)
          if (proc_cell(cellf)==procv) then

             ! On boucle sur les atomes de cette cellule
             do n_at= 1, nato(cellf)
                i_at = atincel(n_at,cellf)

                ! On complete le buffer
                send_nb_val(nproc_voisin) = send_nb_val(nproc_voisin) + 1

                send_buff_int(1,send_nb_val(nproc_voisin),nproc_voisin) = ityp(i_at)
                send_buff_int(2,send_nb_val(nproc_voisin),nproc_voisin) = ielat(i_at)
!                send_buff_int(3,send_nb_val(nproc_voisin),nproc_voisin) = iwmax(i_at)
                send_buff_int(3,send_nb_val(nproc_voisin),nproc_voisin) = num_at_glob(i_at)

                send_buff_dbl(1,send_nb_val(nproc_voisin),nproc_voisin) = xp(1,i_at)
                send_buff_dbl(2,send_nb_val(nproc_voisin),nproc_voisin) = xp(2,i_at)
                send_buff_dbl(3,send_nb_val(nproc_voisin),nproc_voisin) = xp(3,i_at)
!!$                send_buff_dbl(4,send_nb_val(nproc_voisin),nproc_voisin) = ax(1,i_at)
!!$                send_buff_dbl(5,send_nb_val(nproc_voisin),nproc_voisin) = ax(2,i_at)
!!$                send_buff_dbl(6,send_nb_val(nproc_voisin),nproc_voisin) = ax(3,i_at)
                send_buff_dbl(4,send_nb_val(nproc_voisin),nproc_voisin) = vp(1,i_at)
                send_buff_dbl(5,send_nb_val(nproc_voisin),nproc_voisin) = vp(2,i_at)
                send_buff_dbl(6,send_nb_val(nproc_voisin),nproc_voisin) = vp(3,i_at)
                
!!$               if (lsuivinonpbc) then 
!!$		send_buff_dbl(10,send_nb_val(nproc_voisin),nproc_voisin) = xpnonpbc(1,i_at)
!!$                send_buff_dbl(11,send_nb_val(nproc_voisin),nproc_voisin) = xpnonpbc(2,i_at)
!!$                send_buff_dbl(12,send_nb_val(nproc_voisin),nproc_voisin) = xpnonpbc(3,i_at)
!!$
!!$		send_buff_dbl(13,send_nb_val(nproc_voisin),nproc_voisin) = tmpsuivi(1,i_at)
!!$                send_buff_dbl(14,send_nb_val(nproc_voisin),nproc_voisin) = tmpsuivi(2,i_at)
!!$                send_buff_dbl(15,send_nb_val(nproc_voisin),nproc_voisin) = tmpsuivi(3,i_at)
!!$
!!$		send_buff_dbl(16,send_nb_val(nproc_voisin),nproc_voisin) = axnonpbc(1,i_at)
!!$                send_buff_dbl(17,send_nb_val(nproc_voisin),nproc_voisin) = axnonpbc(2,i_at)
!!$                send_buff_dbl(18,send_nb_val(nproc_voisin),nproc_voisin) = axnonpbc(3,i_at)
!!$
!!$               end if
               if ((llangevin.eqv..true.).or.(l2T.eqv..true.))then
                  send_buff_dbl(nb_var_dbl-2,send_nb_val(nproc_voisin),nproc_voisin) = Glangv(1,i_at)
                  send_buff_dbl(nb_var_dbl-1,send_nb_val(nproc_voisin),nproc_voisin) = Glangv(2,i_at)
                  send_buff_dbl(nb_var_dbl,send_nb_val(nproc_voisin),nproc_voisin) = Glangv(3,i_at)
               end if

!LPARAFULLSEND
!                   send_buff_dbl(10,send_nb_val(nproc_voisin),nproc_voisin) = xpp(1,i_at)
!                   send_buff_dbl(11,send_nb_val(nproc_voisin),nproc_voisin) = xpp(2,i_at)
!                   send_buff_dbl(12,send_nb_val(nproc_voisin),nproc_voisin) = xpp(3,i_at)
!                   send_buff_dbl(13,send_nb_val(nproc_voisin),nproc_voisin) = fp(1,i_at)
!                   send_buff_dbl(14,send_nb_val(nproc_voisin),nproc_voisin) = fp(2,i_at)
!                   send_buff_dbl(15,send_nb_val(nproc_voisin),nproc_voisin) = fp(3,i_at)

             enddo  ! fin de boucle sur les atomes

          endif
       enddo    ! fin de boucle sur les cellules fantomes

       ! On envoit les buffers vers le processeur
       call MPI_ISSEND(send_nb_val(nproc_voisin),      1,   MPI_INTEGER,        &
            procv,1001,MPI_COMM_space,send_rqst(nproc_voisin,1),ierr)

       call MPI_ISSEND(send_buff_int(1,1,nproc_voisin),nb_var_int*nb_at_max,MPI_INTEGER,        &
            procv,1002,MPI_COMM_space,send_rqst(nproc_voisin,2),ierr)

       call MPI_ISSEND(send_buff_dbl(1,1,nproc_voisin),nb_var_dbl*nb_at_max,NDM_MPI_REAL_DOUBLE,&
            procv,1003,MPI_COMM_space,send_rqst(nproc_voisin,3),ierr)

    enddo    ! fin de boucle sur les processeurs


  end subroutine envoi_atomes_fantomes


  !------------------------------------------------------------------------!
  ! Procedure dont le but est la reception des nouveaux atomes locaux

  subroutine reception_nouveaux_atomes !seulment MAJ

    USE T_kind_param_m, ONLY:  double
!    use tab_imm_m

    implicit none

    integer :: nb_at_recv
    integer :: proc_source
    integer :: i_at
    integer :: nproc_voisin
    integer :: ind_recv

    ! On boucle sur les processeurs voisins
    do nproc_voisin=1,nbr_proc_voisin

       call MPI_WAITANY(nbr_proc_voisin,recv_rqst(:,1),ind_recv,status,ierr)
       proc_source = proc_voisin(ind_recv)

       ! Reception des nouveaux atomes issus de ce processeur voisin

       call MPI_WAIT(recv_rqst(ind_recv,2), status, ierr)
       call MPI_WAIT(recv_rqst(ind_recv,3), status, ierr)

       ! recopie des infos dans les tableaux locaux
       do i_at = 1, recv_nb_val(ind_recv)

          ! On ajoute un atome a la liste
          im = im + 1
!          imd = imd + 1
!          imf = imf + 1
!          imana = imana + 1


          ! mise a jour des variables entieres
          ityp(im)        = recv_buff_int(1,i_at,ind_recv)
          ielat(im)       = recv_buff_int(2,i_at,ind_recv)
!          iwmax(im)       = recv_buff_int(3,i_at,ind_recv)
          num_at_glob(im) = recv_buff_int(3,i_at,ind_recv)

          ! mise a jour des donnees de la cellule correspondante
          nato(ielat(im)) = nato(ielat(im)) + 1
          atincel(nato(ielat(im)),ielat(im)) = im

          ! mise a jour des variables reelles
          xp(1,im) = recv_buff_dbl(1,i_at,ind_recv) 
          xp(2,im) = recv_buff_dbl(2,i_at,ind_recv) 
          xp(3,im) = recv_buff_dbl(3,i_at,ind_recv) 
!          ax(1,im) = recv_buff_dbl(4,i_at,ind_recv)
!          ax(2,im) = recv_buff_dbl(5,i_at,ind_recv)
!          ax(3,im) = recv_buff_dbl(6,i_at,ind_recv)
          vp(1,im) = recv_buff_dbl(4,i_at,ind_recv) 
          vp(2,im) = recv_buff_dbl(5,i_at,ind_recv) 
          vp(3,im) = recv_buff_dbl(6,i_at,ind_recv) 
	  
!!$	  if (lsuivinonpbc) then
!!$           xpnonpbc(1,im) = recv_buff_dbl(10,i_at,ind_recv)
!!$           xpnonpbc(2,im) = recv_buff_dbl(11,i_at,ind_recv)
!!$           xpnonpbc(3,im) = recv_buff_dbl(12,i_at,ind_recv)
!!$ 
!!$           tmpsuivi(1,im) = recv_buff_dbl(13,i_at,ind_recv)
!!$           tmpsuivi(2,im) = recv_buff_dbl(14,i_at,ind_recv)
!!$           tmpsuivi(3,im) = recv_buff_dbl(15,i_at,ind_recv)
!!$ 
!!$           axnonpbc(1,im) = recv_buff_dbl(16,i_at,ind_recv)
!!$           axnonpbc(2,im) = recv_buff_dbl(17,i_at,ind_recv)
!!$           axnonpbc(3,im) = recv_buff_dbl(18,i_at,ind_recv)
!!$        end if
        if ((llangevin.eqv..true.).or.(l2T.eqv..true.))then
            Glangv(1,im)= recv_buff_dbl(nb_var_dbl-2,i_at,ind_recv)
            Glangv(2,im)= recv_buff_dbl(nb_var_dbl-1,i_at,ind_recv)
            Glangv(3,im)= recv_buff_dbl(nb_var_dbl,i_at,ind_recv)
        end if

!          if(lfrozen)free(im)=.true.
!LPARAFULLSEND
!          xpp(1,im) = recv_buff_dbl(10,i_at,ind_recv)
!          xpp(2,im) = recv_buff_dbl(11,i_at,ind_recv)
!          xpp(3,im) = recv_buff_dbl(12,i_at,ind_recv)
!          fp(1,im) = recv_buff_dbl(13,i_at,ind_recv)
!          fp(2,im) = recv_buff_dbl(14,i_at,ind_recv)
!          fp(3,im) = recv_buff_dbl(15,i_at,ind_recv)

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

  subroutine elimine_atomes_fantomes

    USE T_kind_param_m, ONLY:  double

!    use tab_imm_m

    implicit none

    integer :: i_at
    integer :: i_new
    integer :: j_at
    integer :: koo
    integer :: nb_at_a_eliminer                 ! Nbre d'atomes a eliminer
    integer :: pt_at_elimine                    ! Pointeur sur le dernier atome elimine
    integer, allocatable :: at_a_eliminer(:)    ! Liste des atomes a eliminer


    allocate(at_a_eliminer(im))
    at_a_eliminer = 0

    ! initialisation de la liste des atomes a eliminer
    nb_at_a_eliminer = 0
    do i_at = 1, im
       koo = ielat(i_at)
       if (proc_cell(koo).ne.myid) then
          nb_at_a_eliminer = nb_at_a_eliminer + 1
          at_a_eliminer(nb_at_a_eliminer) = i_at
       endif
    enddo

    if (nb_at_a_eliminer.ne.0) then
       i_new = 1
       pt_at_elimine = 1
       ! On boucle sur tous les atomes
       do i_at = 1, im

          ! Si il s'agit d'un atome a eliminer
          if (i_at.eq.at_a_eliminer(pt_at_elimine)) then

             ! On met a jour les caracteristiques de la cellule correspondante
             koo = ielat(i_at)
             do j_at = 1, nato(koo)
                if ( atincel(j_at,koo).eq.i_at ) then
                   if (j_at.eq.nato(koo)) then
                      atincel(j_at,koo) = 0
                   else
                      atincel(j_at:nato(koo)-1,koo) = atincel(j_at+1:nato(koo),koo)
                   endif
                endif
             enddo
             nato(koo) = nato(koo) - 1

          endif
          
          ! Si les deux pointeurs ne sont pas au meme point, on deplace l'atome courant
          if ( i_new.ne.i_at ) then
             xp(:,i_new)  = xp(:,i_at)
!LPARAFULLSEND
!             xpp(:,i_new) = xpp(:,i_at)
             vp(:,i_new)  = vp(:,i_at)
!             ax(:,i_new)  = ax(:,i_at)
!!$	     if (lsuivinonpbc)  xpnonpbc(:,i_new)  = xpnonpbc(:,i_at)
!!$	     if (lsuivinonpbc)  tmpsuivi(:,i_new)  = tmpsuivi(:,i_at)
!!$	     if (lsuivinonpbc)  axnonpbc(:,i_new)  = axnonpbc(:,i_at)
!             fp(:,i_new)  = fp(:,i_at)
!             if(lfrozen)free(i_new)=free(i_at)

             ityp(i_new)        = ityp(i_at)
             ielat(i_new)       = ielat(i_at)
!             iwmax(i_new)       = iwmax(i_at)
             num_at_glob(i_new) = num_at_glob(i_at)

             ! On met aussi a jour le numero local de l'atome dans la liste de la cellule
             do j_at=1,nato(ielat(i_at))
                if (atincel(j_at,ielat(i_at)).eq.i_at) then
                   atincel(j_at,ielat(i_at))=i_new
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

       im = im - nb_at_a_eliminer
!       imd = imd - nb_at_a_eliminer
!       imf = imf - nb_at_a_eliminer
!       imana = imana - nb_at_a_eliminer

    endif

    deallocate(at_a_eliminer)

    ! On verifie qu'il n'y a plus d'atomes a l'exterieur du domaine local
    do koo=1,noxyz
       if (proc_cell(koo).ne.myid .and. nato(koo).ne.0) print *,'ERREUR !!!',&
            myid,'possede encore',nato(koo),'at. dans la cellule',koo
    enddo

  end subroutine elimine_atomes_fantomes


  !------------------------------------------------------------------------!
  ! Procedure en charge de l'envoi des atomes frontieres du processeur 
  ! courant vers les processeurs voisins concernes

  subroutine envoi_atomes_frontieres

    USE T_kind_param_m, ONLY:  double
!    use tab_imm_m

    implicit none

    integer :: nproc_voisin
    integer :: ncell_front
    integer :: procv
    integer :: koo
    integer :: nb_at,nb_at_max,nb_at_max_tot
    integer :: n_at
    integer :: i_at

    ! Boucle a vide pour determiner au mieux la taille du buffer d'envoi 
    nb_at_max = 0
    do nproc_voisin = 1, nbr_proc_voisin
       nb_at=0
       do ncell_front = 1, nbr_cell_frontiere(nproc_voisin)
          nb_at = nb_at + nato(cell_frontiere(nproc_voisin,ncell_front))
       enddo
       nb_at_max = max(nb_at_max, nb_at)
    enddo
    call MPI_ALLREDUCE(nb_at_max,nb_at_max_tot,1,MPI_INTEGER,MPI_MAX,MPI_COMM_space,ierr)
    nb_at_max=nb_at_max_tot
    ! petite manip pour eviter le cas nb_at_max=0
    nb_at_max = max(nb_at_max,1)

    ! Allocation des buffers
    nb_var_int = 4
!LPARAFULLSEND
!    nb_var_dbl = 15
!!$   if  (lsuivinonpbc) then
!!$    nb_var_dbl = 18
!!$    else  
!!$    nb_var_dbl = 9
!!$   end if 
    nb_var_dbl = 9
    allocate(send_nb_val(nbr_proc_voisin))
    allocate(send_buff_int(nb_var_int,nb_at_max,nbr_proc_voisin))
    allocate(send_buff_dbl(nb_var_dbl,nb_at_max,nbr_proc_voisin))
    allocate(send_rqst(nbr_proc_voisin,3))
    allocate(recv_nb_val(nbr_proc_voisin))
    allocate(recv_buff_int(nb_var_int,nb_at_max,nbr_proc_voisin))
    allocate(recv_buff_dbl(nb_var_dbl,nb_at_max,nbr_proc_voisin))
    allocate(recv_rqst(nbr_proc_voisin,3))

    ! Preparation des receptions
    do nproc_voisin = 1, nbr_proc_voisin
       procv = proc_voisin(nproc_voisin)
       call MPI_IRECV(recv_nb_val(nproc_voisin), 1, MPI_INTEGER, procv, 2001, MPI_COMM_space, recv_rqst(nproc_voisin,1), ierr)
       call MPI_IRECV(recv_buff_int(1,1,nproc_voisin), nb_var_int*nb_at_max, MPI_INTEGER, procv, 2002, &
            MPI_COMM_space, recv_rqst(nproc_voisin,2), ierr)
       call MPI_IRECV(recv_buff_dbl(1,1,nproc_voisin), nb_var_dbl*nb_at_max, NDM_MPI_REAL_DOUBLE, procv, 2003, &
            MPI_COMM_space, recv_rqst(nproc_voisin,3), ierr)
    enddo


    ! On boucle sur les processeurs voisins
    do nproc_voisin = 1, nbr_proc_voisin
       procv = proc_voisin(nproc_voisin)

       send_nb_val(nproc_voisin) = 0

       ! On boucle sur les cellules frontieres associees au processeur
       do ncell_front = 1, nbr_cell_frontiere(nproc_voisin)
          koo = cell_frontiere(nproc_voisin,ncell_front)

          ! On copie le contenu de la cellule dans le buffer d'envoi
          do n_at = 1, nato(koo)
             i_at = atincel(n_at,koo)

             ! On complete le buffer
             send_nb_val(nproc_voisin) = send_nb_val(nproc_voisin) + 1

             send_buff_int(1,send_nb_val(nproc_voisin),nproc_voisin) = ityp(i_at)
             send_buff_int(2,send_nb_val(nproc_voisin),nproc_voisin) = ielat(i_at)
!             send_buff_int(3,send_nb_val(nproc_voisin),nproc_voisin) = iwmax(i_at)
             send_buff_int(3,send_nb_val(nproc_voisin),nproc_voisin) = num_at_glob(i_at)

             send_buff_dbl(1,send_nb_val(nproc_voisin),nproc_voisin) = xp(1,i_at)
             send_buff_dbl(2,send_nb_val(nproc_voisin),nproc_voisin) = xp(2,i_at)
             send_buff_dbl(3,send_nb_val(nproc_voisin),nproc_voisin) = xp(3,i_at)
!             send_buff_dbl(4,send_nb_val(nproc_voisin),nproc_voisin) = ax(1,i_at)
!             send_buff_dbl(5,send_nb_val(nproc_voisin),nproc_voisin) = ax(2,i_at)
!             send_buff_dbl(6,send_nb_val(nproc_voisin),nproc_voisin) = ax(3,i_at)
             send_buff_dbl(4,send_nb_val(nproc_voisin),nproc_voisin) = vp(1,i_at)
             send_buff_dbl(5,send_nb_val(nproc_voisin),nproc_voisin) = vp(2,i_at)
             send_buff_dbl(6,send_nb_val(nproc_voisin),nproc_voisin) = vp(3,i_at)
	     
!!$	      if (lsuivinonpbc) then
!!$              !
!!$	       send_buff_dbl(10,send_nb_val(nproc_voisin),nproc_voisin) = xpnonpbc(1,i_at)
!!$               send_buff_dbl(11,send_nb_val(nproc_voisin),nproc_voisin) = xpnonpbc(2,i_at)
!!$               send_buff_dbl(12,send_nb_val(nproc_voisin),nproc_voisin) = xpnonpbc(3,i_at)
!!$	      !
!!$	       send_buff_dbl(13,send_nb_val(nproc_voisin),nproc_voisin) = tmpsuivi(1,i_at)
!!$               send_buff_dbl(14,send_nb_val(nproc_voisin),nproc_voisin) = tmpsuivi(2,i_at)
!!$               send_buff_dbl(15,send_nb_val(nproc_voisin),nproc_voisin) = tmpsuivi(3,i_at)
!!$	      !
!!$	       send_buff_dbl(16,send_nb_val(nproc_voisin),nproc_voisin) = axnonpbc(1,i_at)
!!$               send_buff_dbl(17,send_nb_val(nproc_voisin),nproc_voisin) = axnonpbc(2,i_at)
!!$               send_buff_dbl(18,send_nb_val(nproc_voisin),nproc_voisin) = axnonpbc(3,i_at)
!!$	      !
!!$	      end if
             	     
!LPARAFULLSEND
!             send_buff_dbl(10,send_nb_val(nproc_voisin),nproc_voisin) = xpp(1,i_at)
!             send_buff_dbl(11,send_nb_val(nproc_voisin),nproc_voisin) = xpp(2,i_at)
!             send_buff_dbl(12,send_nb_val(nproc_voisin),nproc_voisin) = xpp(3,i_at)
!             send_buff_dbl(13,send_nb_val(nproc_voisin),nproc_voisin) = fp(1,i_at)
!             send_buff_dbl(14,send_nb_val(nproc_voisin),nproc_voisin) = fp(2,i_at)
!             send_buff_dbl(15,send_nb_val(nproc_voisin),nproc_voisin) = fp(3,i_at)

          enddo

       enddo   ! fin de boucle sur les cellules

       ! On envoit les buffers vers le processeur
       call MPI_ISEND(send_nb_val(nproc_voisin),      1,   MPI_INTEGER,        &
            procv,2001,MPI_COMM_space,send_rqst(nproc_voisin,1),ierr)

       call MPI_ISEND(send_buff_int(1,1,nproc_voisin),nb_var_int*nb_at_max,MPI_INTEGER,        &
            procv,2002,MPI_COMM_space,send_rqst(nproc_voisin,2),ierr)

       call MPI_ISEND(send_buff_dbl(1,1,nproc_voisin),nb_var_dbl*nb_at_max,NDM_MPI_REAL_DOUBLE,&
            procv,2003,MPI_COMM_space,send_rqst(nproc_voisin,3),ierr)

    enddo   ! fin de boucle sur les processeurs voisins


  end subroutine envoi_atomes_frontieres

  !------------------------------------------------------------------------!
  ! Procedure en charge d'attendre la fin des envois des atomes (frontieres
  ! ou fantomes) et la liberation des buffers d'envoi

  subroutine finalisation_envoi_atomes

    USE T_kind_param_m, ONLY:  double
!    use tab_imm_m

    implicit none

    integer :: nproc_voisin

    integer, allocatable :: send_status(:,:,:)


    allocate(send_status(MPI_STATUS_SIZE,nbr_proc_voisin,3))

    ! Attente de finalisation des envois 

    do nproc_voisin= 1, nbr_proc_voisin
       call MPI_Wait( send_rqst(nproc_voisin,1), send_status(1,nproc_voisin,1), ierr )
       call MPI_Wait( send_rqst(nproc_voisin,2), send_status(1,nproc_voisin,2), ierr )
       call MPI_Wait( send_rqst(nproc_voisin,3), send_status(1,nproc_voisin,3), ierr )
    enddo

    ! Liberation des buffers

    deallocate(send_status)
    deallocate(send_rqst)
    deallocate(send_nb_val)
    deallocate(send_buff_int)
    deallocate(send_buff_dbl)

    deallocate(recv_rqst)
    deallocate(recv_nb_val)
    deallocate(recv_buff_int)
    deallocate(recv_buff_dbl)

  end subroutine finalisation_envoi_atomes


  !------------------------------------------------------------------------!
  ! Procedure en charge de la reception des nouveaux atomes fantomes en 
  ! provenance des processeurs voisins

  subroutine reception_atomes_fantomes !seulement MAJ

    USE T_kind_param_m, ONLY:  double
!    use tab_imm_m

    implicit none

    integer :: nb_at_recv
    integer :: proc_source
    integer :: i_at
    integer :: nproc_voisin
    integer :: pt_at_ftm
    integer :: ind_recv

    integer :: nb_at_max,nb_at,koo

    ! On place le pointeur de stockage des atomes fantomes a la suite des 
    ! atomes locaux
    pt_at_ftm = im

    ! On boucle sur les processeurs voisins
    do nproc_voisin=1,nbr_proc_voisin

       call MPI_WAITANY(nbr_proc_voisin,recv_rqst(:,1),ind_recv,status,ierr)
       proc_source = proc_voisin(ind_recv)

       ! Reception des nouveaux atomes issus de ce processeur voisin
       call MPI_WAIT(recv_rqst(ind_recv,2), status, ierr)
       call MPI_WAIT(recv_rqst(ind_recv,3), status, ierr)

       ! recopie des infos dans les tableaux locaux au niveau des atomes fantomes
       do i_at = 1, recv_nb_val(ind_recv)

          ! On ajoute un atome fantome a la liste
          pt_at_ftm = pt_at_ftm + 1

          ! mise a jour des variables entieres
          ityp(pt_at_ftm)  = recv_buff_int(1,i_at,ind_recv)
          ielat(pt_at_ftm) = recv_buff_int(2,i_at,ind_recv)
!          iwmax(pt_at_ftm) = recv_buff_int(3,i_at,ind_recv)
          num_at_glob(pt_at_ftm) = recv_buff_int(3,i_at,ind_recv)

          ! mise a jour des donnees de la cellule correspondante
          nato(ielat(pt_at_ftm)) = nato(ielat(pt_at_ftm)) + 1
          atincel(nato(ielat(pt_at_ftm)),ielat(pt_at_ftm)) = pt_at_ftm

          ! mise a jour des variables reelles
          xp(1,pt_at_ftm) = recv_buff_dbl(1,i_at,ind_recv) 
          xp(2,pt_at_ftm) = recv_buff_dbl(2,i_at,ind_recv) 
          xp(3,pt_at_ftm) = recv_buff_dbl(3,i_at,ind_recv) 
!          ax(1,pt_at_ftm) = recv_buff_dbl(4,i_at,ind_recv)
!          ax(2,pt_at_ftm) = recv_buff_dbl(5,i_at,ind_recv)
!          ax(3,pt_at_ftm) = recv_buff_dbl(6,i_at,ind_recv)
          vp(1,pt_at_ftm) = recv_buff_dbl(4,i_at,ind_recv) 
          vp(2,pt_at_ftm) = recv_buff_dbl(5,i_at,ind_recv) 
          vp(3,pt_at_ftm) = recv_buff_dbl(6,i_at,ind_recv) 
!!$	  if (lsuivinonpbc) then
!!$          !
!!$	   xpnonpbc(1,pt_at_ftm) = recv_buff_dbl(10,i_at,ind_recv)
!!$           xpnonpbc(2,pt_at_ftm) = recv_buff_dbl(11,i_at,ind_recv)
!!$           xpnonpbc(3,pt_at_ftm) = recv_buff_dbl(12,i_at,ind_recv)
!!$	  !
!!$	   tmpsuivi(1,pt_at_ftm) = recv_buff_dbl(13,i_at,ind_recv)
!!$           tmpsuivi(2,pt_at_ftm) = recv_buff_dbl(14,i_at,ind_recv)
!!$           tmpsuivi(3,pt_at_ftm) = recv_buff_dbl(15,i_at,ind_recv)
!!$	  !
!!$	   axnonpbc(1,pt_at_ftm) = recv_buff_dbl(16,i_at,ind_recv)
!!$           axnonpbc(2,pt_at_ftm) = recv_buff_dbl(17,i_at,ind_recv)
!!$           axnonpbc(3,pt_at_ftm) = recv_buff_dbl(18,i_at,ind_recv)
!!$	  !
!!$	  end if
!!$	  
	  
!LPARAFULLSEND
!          xpp(1,pt_at_ftm) = recv_buff_dbl(10,i_at,ind_recv)
!          xpp(2,pt_at_ftm) = recv_buff_dbl(11,i_at,ind_recv)
!          xpp(3,pt_at_ftm) = recv_buff_dbl(12,i_at,ind_recv)
!          fp(1,pt_at_ftm) = recv_buff_dbl(13,i_at,ind_recv)
!          fp(2,pt_at_ftm) = recv_buff_dbl(14,i_at,ind_recv)
!          fp(3,pt_at_ftm) = recv_buff_dbl(15,i_at,ind_recv)

       enddo

    enddo   ! fin de boucle sur les processeurs voisins

  end subroutine reception_atomes_fantomes

  !------------------------------------------------------------------------!
  ! Procedure dont le but est l'envoi des valeurs de tabdensity pour les 
  ! atomes frontieres

  subroutine envoi_tabdensity_frontieres(tabdensity,imm,num_at_glob)

    USE T_kind_param_m, ONLY:  double
!    use tab_imm_m

    implicit none
    integer::imm
    integer :: nproc_voisin
    integer :: ncell_front
    integer :: procv
    integer :: koo
    integer :: nb_at,nb_at_max,nb_at_max_tot
    integer :: n_at
    integer :: i_at
    real(double) :: tabdensity(imm)
    integer,intent(in)::num_at_glob(imm)
    
    ! Boucle a vide pour determiner au mieux la taille du buffer d'envoi 
    nb_at_max = 0
    do nproc_voisin = 1, nbr_proc_voisin
       nb_at=0
       do ncell_front = 1, nbr_cell_frontiere(nproc_voisin)
          nb_at = nb_at + nato(cell_frontiere(nproc_voisin,ncell_front))
       enddo
       nb_at_max = max(nb_at_max, nb_at)
    enddo
    call MPI_ALLREDUCE(nb_at_max,nb_at_max_tot,1,MPI_INTEGER,MPI_MAX,MPI_COMM_space,ierr)
    nb_at_max=nb_at_max_tot
    ! petite manip pour eviter le cas nb_at_max=0
    nb_at_max = max(nb_at_max,1)

    ! Allocation des buffers
    allocate(send_nb_val(nbr_proc_voisin))
    allocate(send_buff_int(1,nb_at_max,nbr_proc_voisin))
    allocate(send_buff_dbl(1,nb_at_max,nbr_proc_voisin))
    allocate(send_rqst(nbr_proc_voisin,3))
    allocate(recv_nb_val(nbr_proc_voisin))
    allocate(recv_buff_int(1,nb_at_max,nbr_proc_voisin))
    allocate(recv_buff_dbl(1,nb_at_max,nbr_proc_voisin))
    allocate(recv_rqst(nbr_proc_voisin,3))

    ! Preparation des receptions
    do nproc_voisin = 1, nbr_proc_voisin
       procv = proc_voisin(nproc_voisin)
       call MPI_IRECV(recv_nb_val(nproc_voisin), 1, MPI_INTEGER, procv, 3001, MPI_COMM_space, recv_rqst(nproc_voisin,1), ierr)
       call MPI_IRECV(recv_buff_int(1,1,nproc_voisin), nb_at_max, MPI_INTEGER, procv, 3002, &
            MPI_COMM_space, recv_rqst(nproc_voisin,2), ierr)
       call MPI_IRECV(recv_buff_dbl(1,1,nproc_voisin), nb_at_max, NDM_MPI_REAL_DOUBLE, procv, 3003, &
            MPI_COMM_space, recv_rqst(nproc_voisin,3), ierr)
    enddo


    ! On boucle sur les processeurs voisins
    do nproc_voisin = 1, nbr_proc_voisin
       procv = proc_voisin(nproc_voisin)

       send_nb_val(nproc_voisin) = 0

       ! On boucle sur les cellules frontieres associees au processeur
       do ncell_front = 1, nbr_cell_frontiere(nproc_voisin)
          koo = cell_frontiere(nproc_voisin,ncell_front)

          ! On copie le contenu de la cellule dans le buffer d'envoi
          do n_at = 1, nato(koo)
             i_at = atincel(n_at,koo)

             ! On complete le buffer
             send_nb_val(nproc_voisin) = send_nb_val(nproc_voisin) + 1

             send_buff_int(1,send_nb_val(nproc_voisin),nproc_voisin) = num_at_glob(i_at)

             send_buff_dbl(1,send_nb_val(nproc_voisin),nproc_voisin) = tabdensity(i_at)

          enddo

       enddo   ! fin de boucle sur les cellules

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

  subroutine reception_tabdensity_fantomes(tabdensity,imm,num_at_glob)

    USE T_kind_param_m, ONLY:  double
!    use tab_imm_m
    implicit none
    integer::imm
    integer,intent(in)::num_at_glob(imm)
    integer :: nb_at_recv
    integer :: proc_source
    integer :: i_at
    integer :: nproc_voisin
    integer :: pt_at_ftm
    integer :: ind_recv
    integer :: nb_at_max,nb_at,koo
    integer :: ind_loc,ftm_at

    real(double) :: tabdensity(imm)
    real(double) :: temps_exe

    ! On place le pointeur de stockage des atomes fantomes a la suite des 
    ! atomes locaux
    pt_at_ftm = im

    ! On boucle sur les processeurs voisins
    do nproc_voisin=1,nbr_proc_voisin

       call MPI_WAITANY(nbr_proc_voisin,recv_rqst(:,1),ind_recv,status,ierr)
       proc_source = proc_voisin(ind_recv)

       ! Reception des tabdensity issus de ce processeur voisin
       call MPI_WAIT(recv_rqst(ind_recv,2), status, ierr)
       call MPI_WAIT(recv_rqst(ind_recv,3), status, ierr)

       ! recopie des infos dans les tableaux locaux au niveau des atomes fantomes
       do i_at = 1, recv_nb_val(ind_recv)

          ! On boucle pour trouver l'indice local de l'atome fantome courant
          ind_loc=-1
          do ftm_at=im+1,imm
             if (num_at_glob(ftm_at)==recv_buff_int(1,i_at,ind_recv)) then
                ind_loc=ftm_at
                exit
             endif
          enddo
          if (ind_loc==-1) then
             print *,myid,'!!!Pb!!! Reception du proc',proc_source,'d''un atome fantome inexistant'
             temps_exe = MPI_Wtime() - temps_deb
             if (myid==0) then
                print *, 'Temps d''execution : ', temps_exe
                print *, 'Temps d''init      : ', temps_init
                print *, 'Temps d''input     : ', temps_input
                print *, 'Temps de config   : ', temps_config
                print *, 'Temps d''initspeed : ', temps_initspeed
                print *, 'Temps para estime : ', temps_para
                print *, 'Temps dmloop : ', temps_dmloop
             endif
             call MPI_FINALIZE(ierr)
             stop 
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

  subroutine envoi_fp_fantomes

    USE T_kind_param_m, ONLY:  double
!    use tab_imm_m

    implicit none

    integer :: nproc_voisin
    integer :: ncell_ftm
    integer :: procv
    integer :: koo
    integer :: nb_at,nb_at_max,nb_at_max_tot
    integer :: n_at
    integer :: i_at

    ! Boucle a vide pour determiner au mieux la taille du buffer d'envoi 
    nb_at_max = 0
    do nproc_voisin = 1, nbr_proc_voisin
       nb_at=0
       do ncell_ftm = 1, nbr_cell_ftm
          if (proc_cell(cell_ftm(ncell_ftm)).eq.proc_voisin(nproc_voisin)) then
             nb_at = nb_at + nato(cell_ftm(ncell_ftm))
          endif
       enddo
       nb_at_max = max(nb_at_max, nb_at)
    enddo
    call MPI_ALLREDUCE(nb_at_max,nb_at_max_tot,1,MPI_INTEGER,MPI_MAX,MPI_COMM_space,ierr)
    nb_at_max=nb_at_max_tot
    ! petite manip pour eviter le cas nb_at_max=0
    nb_at_max = max(nb_at_max,1)

    ! Allocation des buffers
    allocate(send_nb_val(nbr_proc_voisin))
    allocate(send_buff_int(1,nb_at_max,nbr_proc_voisin))
    allocate(send_buff_dbl(3,nb_at_max,nbr_proc_voisin))
    allocate(send_rqst(nbr_proc_voisin,3))
    allocate(recv_nb_val(nbr_proc_voisin))
    allocate(recv_buff_int(1,nb_at_max,nbr_proc_voisin))
    allocate(recv_buff_dbl(3,nb_at_max,nbr_proc_voisin))
    allocate(recv_rqst(nbr_proc_voisin,3))

    ! Preparation des receptions
    do nproc_voisin = 1, nbr_proc_voisin
       procv = proc_voisin(nproc_voisin)
       call MPI_IRECV(recv_nb_val(nproc_voisin), 1, MPI_INTEGER, procv, 4001, MPI_COMM_space, recv_rqst(nproc_voisin,1), ierr)
       call MPI_IRECV(recv_buff_int(1,1,nproc_voisin), nb_at_max, MPI_INTEGER, procv, 4002, &
            MPI_COMM_space, recv_rqst(nproc_voisin,2), ierr)
       call MPI_IRECV(recv_buff_dbl(1,1,nproc_voisin), nb_at_max*3, NDM_MPI_REAL_DOUBLE, procv, 4003, &
            MPI_COMM_space, recv_rqst(nproc_voisin,3), ierr)
    enddo


    ! On boucle sur les processeurs voisins
    do nproc_voisin = 1, nbr_proc_voisin
       procv = proc_voisin(nproc_voisin)

       send_nb_val(nproc_voisin) = 0

       ! On boucle sur les cellules fantomes susceptibles d'appartenir au processeur
       do ncell_ftm = 1, nbr_cell_ftm
          koo = cell_ftm(ncell_ftm)

          ! appartient-elle au processeur voisin courant?
          if (proc_cell(koo).eq.procv) then

             ! On copie le contenu de la cellule dans le buffer d'envoi
             do n_at = 1, nato(koo)
                i_at = atincel(n_at,koo)
                
                ! On complete le buffer
                send_nb_val(nproc_voisin) = send_nb_val(nproc_voisin) + 1
                
                send_buff_int(1,send_nb_val(nproc_voisin),nproc_voisin) = num_at_glob(i_at)
                
                send_buff_dbl(1:3,send_nb_val(nproc_voisin),nproc_voisin) = fp(1:3,i_at)
                
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

  subroutine reception_fp_frontieres

    USE T_kind_param_m, ONLY:  double
!    use tab_imm_m
    implicit none

    integer :: nb_at_recv
    integer :: proc_source
    integer :: i_at
    integer :: nproc_voisin
    integer :: ind_recv
    integer :: nb_at_max,nb_at,koo
    integer :: ind_loc,i_at_loc,ind_glob
    real(double) :: temps_exe

 
    ! On boucle sur les processeurs voisins
    do nproc_voisin=1,nbr_proc_voisin

       call MPI_WAITANY(nbr_proc_voisin,recv_rqst(:,1),ind_recv,status,ierr)
       proc_source = proc_voisin(ind_recv)

       ! Reception des tabdensity issus de ce processeur voisin
       call MPI_WAIT(recv_rqst(ind_recv,2), status, ierr)
       call MPI_WAIT(recv_rqst(ind_recv,3), status, ierr)

       ! recopie des infos dans les tableaux locaux
       do i_at = 1, recv_nb_val(ind_recv)

          ind_glob = recv_buff_int(1,i_at,ind_recv)
 
          ! On boucle pour trouver l'indice local de l'atome fantome courant
          ind_loc=-1
          do i_at_loc=1,im
             if (num_at_glob(i_at_loc)==ind_glob) then
                ind_loc=i_at_loc
                exit
             endif
          enddo
          if (ind_loc==-1) then
             print *,myid,'!!!Pb!!! Reception du proc',proc_source,'d''un atome non local'
             temps_exe = MPI_Wtime() - temps_deb
             if (myid==0) then
                print *, 'Temps d''execution : ', temps_exe
                print *, 'Temps d''init      : ', temps_init
                print *, 'Temps d''input     : ', temps_input
                print *, 'Temps de config   : ', temps_config
                print *, 'Temps d''initspeed : ', temps_initspeed
                print *, 'Temps para estime : ', temps_para
                print *, 'Temps dmloop : ', temps_dmloop
             endif
             call MPI_FINALIZE(ierr)
             stop 
             stop
             !call arret_ndm
          endif

          ! On ajoute a cet atome local la valeur de fp recue
 
          fp(1:3,ind_loc) = fp(1:3,ind_loc) + recv_buff_dbl(1:3,i_at,ind_recv)
 
       enddo

    enddo   ! fin de boucle sur les processeurs voisins

  end subroutine reception_fp_frontieres

#endif
end module mod_para

