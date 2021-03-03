module endrunT_mod
  USE analyseT_mod,only:analyseT
  USE adf_mod,only:adf
  USE spebc_fin_mod,only:spebc_fin
  USE arret_ndm_mod,only:arret_ndm
  USE sauvegardeT_mod,only:sauvegardeT!,cin2gin
  USE rdf_mod,only:rdf
  USE rasmolT_mod,only:rasmolT
  USE gen_com_m, ONLY:itesauv,lprtfat,lwgin,angst,unitP,cunitP,erg2eV,itdes,&
       &iteanapos,iteangle,itecfg,iterasmol,itesigma,itetemp,ldesinteg,linstantfda,&
       &linstantrdf,lpkbar,lprteat,lprteattotm,lprtsigat,parallele,unitP,iterdf,&
       &lwgin,nstepdes, lposmoy,l2T,angst,dmtype,ibound,it,lenfnam,rang,timel,&
       &formatsauv,fnamcout,fnam,lspaceNDM
  use var_pot, only: eatref,eatref,eatref
  USE cellconfig,only:cell_config,caltabtC
  USE atomconfig,only:atom_config,atom_config_d,atom_config_e!, ndm2config, config2ndm
  use boxconfig,only: box_config!,ndm2boxconfig,boxconfig2ndm
  use posana,only:anapos

  implicit none
contains
  ! ****************************************************************
  subroutine endrunT(atdml,celndm,boxndm,latcomp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE tab_imm_m,only:posmoyx,ityp,xp,num_at_glob,fp,vp,iwmax
#if defined ML && defined PARAML
    USE time_measure
#endif
    USE posana,only:
    USE cfg_module,only:
    USE elec_cell, ONLY:  sauveelec
    !       version MPI du 07 f if (allocated(eatom)) eatom(:)=0

    ! ****************************************************************

    implicit none
#ifdef PARA
    include 'mpif.h'
    integer, dimension( MPI_STATUS_SIZE) :: statut2
#endif

    type(box_config)::boxndm
    class(atom_config)::atdml
    type(cell_config):: celndm
    logical,intent(in)::latcomp


    integer :: i,j, n, nAux_real
    CHARACTER(len=100) :: out_file
    REAL(kind(0.d0)), dimension(:,:), allocatable :: aux_real
    CHARACTER(len=20), dimension(:), allocatable :: aux_title
!!$#ifdef PARA
!!$    integer :: iproc
!!$    real(double), allocatable :: xp_loc(:,:),eatom_loc(:)
!!$    integer, allocatable      :: ityp_loc(:)
!!$    integer, allocatable      :: num_at_glob_loc(:)
!!$    integer :: im_loc
!!$    integer :: proc_source
!!$#endif
    !-----------------------------------------------
    !
    !
    !
    if (lPkbar) then
       unitP=1.0d-9
       cunitP='kbar'
    else
       unitP=1.0
       cunitP='d/cm2'
    endif

    ! Un dernier calcul des forces pour la route
    IF (iteTemp.GE.0) iteTemp=1
    IF (iteSigma.GE.0) iteSigma=1

    !flag_fin = .true. !*!
    if(ibound.ne.0) Call spebc_fin (.true.) !*!





 if (rang==0) then
#if defined ML && defined PARAML
    write (6, *) 'ML: neighbours  time',  temps_neigh
    write (6, *) 'ML: energy      time',  temps_energy
    write (6, *) 'ML: force       time',  temps_force
    write (6, *) 'ML: stress      time',  temps_stress
    write (6, *) 'ML: descriptors time',  temps_descripteurs
#endif

    write (6, *)
    write (6, *)

    write (6, *) '####### END OF RUN  ######## = ', it, '  time = ', timel
 endif

! if (lWgin.eqv..true.) call cin2gin
 IF (iteSauv.GE.0) then
    !    call boxndm%print
    !    call celndm%print
    !     call atdml%print
    formatsauv=3;fnamcout= fnam(1:lenfnam)//'.cout'
    call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp)     ! Modif E. Clouet: sauvegarde seulement si voulu
    if (l2T.and.rang==0) call sauveelec
 end IF

 if (.not.linstantrdf) then
    if (iterdf>=0) call rdf
 endif

 if (.not.linstantfda) then
    if (iteangle>=0) call adf
 endif
 if ((dmtype==2).or.(dmtype==3).or.(dmtype==30)) then
    it=0
 end if

 !  select type(atdml)
 !  type is (atom_config_d)
 !     call analyseT(atdml,celndm,boxndm)
 !  end select


 if (iterasmol.GE.0) call rasmolT (atdml,boxndm,it,latcomp=latcomp)
 if (.not.parallele.and.iteanapos>=0) call anapos (it)

 call arret_ndm


 stop
 return
end subroutine endrunT
end module endrunT_mod
