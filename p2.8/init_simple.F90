module init_simple_mod

  use init_pot_mod,only:init_pot,init_pot2
  USE initspeed_mod,only: initspeed
  USE caltabi_mod,only: caltabi
  USE atomconfig,only:atom_config,atom_config_d,atom_config_e
  USE cellconfig, only:cell_config,caltabtC
  use boxconfig,only: box_config
  USE constrconf_mod, only :constrconf

#ifdef PARA
  USE init_vois_mod,only: init_voisinage
#endif
#ifdef LAMMPS_VERSION
  use lammps_util_mod
  use vars_lammps
#endif

  USE gen_com_m, ONLY:fnam,lenfnam,igen,lperiod,lrestart,rang,tstep,two,usdh,&
       &lspacendm
  use read_val,only:ltabvois
  USE var_pot, ONLY:ipotentiel
  use Tpara,only:para_space_config
  implicit none

contains
  ! **************************************************************
  subroutine init_simple(atdml,celndm,boxndm,filename,psc,linitpot)



    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE arret_ndm_mod,only: arret_ndm
    use posana,only:anapos
    USE posana,only:
    USE elec_cell,ONLY: i2t,t_cpl, readelec
    USE eloss, ONLY : ibrake,ecelec,initeloss

#ifdef PARA
    USE Tpara,only:COMM_space,myidsp,nprocspace

#else
    use Tpara,only:nprocspace

#endif

    ! **************************************************************

    implicit none
     type(para_space_config)::psc
     class(atom_config)::atdml
     type(cell_config),intent(out)::celndm
     class(box_config),intent(out)::boxndm
    character(*),optional::filename
    logical, optional::linitpot
    
    logical::linitpotW=.true.
    character*80::filenomIS
    integer :: i, lufilmpaf,itapp,j,lenfn2,ipath,ierr
    !-----------------------------------------------
    character*2::extension
    logical :: lrepart

    filenomIS=fnam(1:lenfnam)
    if (present(filename))filenomIS=filename
    if (present(linitpot))linitpotW=linitpot
    if (linitpotW)call init_pot
    usdh = 1/(two*tstep)
!    if ((ipotentiel==-10).or.(ipotentiel==-11))then
       lrepart=.false. !TOUJOURS FALSE, repartition plus tard
!    else
!       lrepart=.true.
       !    end if
       call constrconf(atdml,boxndm,celndm,lrepart,filenomIS,psc)
    call init_pot2(boxndm,atdml%imm)

#ifdef PARA
    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
       CALL comm_space%barrier
       call init_voisinage(celndm,psc)

       !if (rang==0)  write(6,*) 'NOMBRE DE CELLULES FRONTIERES ASSOCIEES A CHAQUE PROCESSEUR'
       !write(6,*) 'Le proc ',myidsp,' a ',psc%nbr_proc_voisin,' processeur voisin'
    end if
#endif
    !<---------end setting the cell diviion ----------------------
    if (.not.lrestart) then
       !    if (rang==0)     write(6,*)'>>>>>>>>>>>avant initspeed'
       select type(atdml)
          class is (atom_config_d)
             call initspeed(atdml,boxndm)
       end select
    end if
    if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       call caltabtC(celndm,atdml,lperiod,boxndm,psc=psc)
    else
       call caltabtC(celndm,atdml,lperiod,boxndm)
    end if
    if (ltabvois) then
       call caltabi(atdml,celndm,boxndm)
    end if

    return
  end subroutine init_simple

end module init_simple_mod
