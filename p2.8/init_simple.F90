module init_simple_mod

  use init_pot_mod,only:init_pot,init_pot2
  USE setcell,only:setcellconf
  USE neigcel_mod,only: neigcel
!  USE dynalloccell
  USE initspeed_mod,only: initspeed
  USE caltabi_mod,only: caltabi
  USE prtplz_mod,only: prtplz
  USE atomconfig,only:atom_config,atom_config_d,atom_config_e
  USE cellconfig, only:cell_config,caltabtC,init_cel
  use boxconfig,only: box_config
  USE constrconf_mod, only :constrconf

#ifdef PARA
  USE init_vois_mod,only: init_voisinage
#endif
#ifdef LAMMPS_VERSION
  use lammps_util_mod
  use vars_lammps
#endif

  USE gen_com_m, ONLY:fnam,lenfnam,igen,lperiod,lrestart,rang,im_glob,tstep,two,usdh,&
       &lspacendm, posa, forca,latcomp
  use read_val,only:ltabvois
  USE var_pot, ONLY:ipotentiel
  implicit none

contains
  ! **************************************************************
  subroutine init_simple(atdml,celndm,boxndm,filename)



    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE arret_ndm_mod,only: arret_ndm
    use posana,only:anapos
    USE posana,only:
    USE defcdp, ONLY :itecdp
    USE elec_cell,ONLY: i2t,t_cpl, readelec
    USE eloss, ONLY : ibrake,ecelec,initeloss

#ifdef PARA
    use mpi
    USE mod_para,only:MPI_COMM_space,TEMPS_INPUT_DEB,TEMPS_INPUT,TEMPS_CONFIG_DEB,TEMPS_CONFIG,myidsp,&
         &NBR_PROC_VOISIN,TEMPS_INITSPEED_DEB,TEMPS_INITSPEED,nprocspace
#else
    use mod_para,only:nprocspace

#endif

    ! **************************************************************

    implicit none
    class(atom_config)::atdml
    type(cell_config),intent(out)::celndm
    type(box_config),intent(out)::boxndm
    character(*),optional::filename
    character*80::filenomIS
    integer :: i, lufilmpaf,itapp,j,lenfn2,ipath,ierr
    !-----------------------------------------------
    character*2::extension
    logical :: lrepart

    filenomIS=fnam(1:lenfnam)
    if (present(filename))filenomIS=filename
    call init_pot
    usdh = 1/(two*tstep)
!    if ((ipotentiel==-10).or.(ipotentiel==-11))then
       lrepart=.false. !TOUJOURS FALSE, repartition plus tard
!    else
!       lrepart=.true.
!    end if
    call constrconf(atdml,boxndm,celndm,lrepart,filenomIS)
    call init_pot2(boxndm)

#ifdef PARA
    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
       CALL MPI_BARRIER(MPI_COMM_space,ierr)
       call init_voisinage(celndm)

       if (rang==0)  write(6,*) 'NOMBRE DE CELLULES FRONTIERES ASSOCIEES A CHAQUE PROCESSEUR'
       write(6,*) 'Le proc ',myidsp,' a ',nbr_proc_voisin,' processeur voisin'
    end if
#endif
    !<---------end setting the cell diviion ----------------------
    if (.not.lrestart) then
       !    if (rang==0)     write(6,*)'>>>>>>>>>>>avant initspeed'
       select type(atdml)
          class is (atom_config_d)
          call initspeed(atdml,im_glob,boxndm)
       end select
    end if
    call caltabtC(celndm,atdml,lperiod,boxndm)
    if (ltabvois) then
       call caltabi(atdml,celndm,boxndm)
    end if

    return
  end subroutine init_simple

end module init_simple_mod
