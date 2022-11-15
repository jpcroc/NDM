module init_mod
   USE arret_ndm_mod,only:arret_ndm
  use init_pot_mod,only:init_pot,init_pot2
  USE transf_mod,only: transf
  USE initspeed_mod,only: initspeed
  USE sauvegardeT_mod,only: sauvegardeT!,cin2gin
  USE caltabi_mod,only: caltabi
  USE cdp_mod,only: initcdp
  USE initcasca_mod,only: initcasca
  USE deftimestep_mod,only: deftimestep
  USE rasmolT_mod,only: rasmolT
  USE neb_module,only: constrconfNEB,atneb,cellneb,boxneb
  USE atomconfig,only:atom_config,atom_config_d,atom_config_e
  USE cellconfig, only:cell_config,caltabtC
  use boxconfig,only: box_config
  USE constrconf_mod, only :constrconf
  USE arret_ndm_mod,only: arret_ndm

#ifdef PARA
  USE init_vois_mod,only: init_voisinage
#endif
#ifdef LAMMPS_VERSION
  use lammps_util_mod
  use vars_lammps
#endif
#ifdef ML
  use NDM_ML,only:init_config_ml
#endif
  

  
  use Tpara,only:para_space_config

  USE gen_com_m, ONLY:fnam,lenfnam,dmtype,fnamcout,igen,ilangevin,iteration,iteanapos,iterasmol,&
       &itetimestep,kinemean,lcasca,lperiod,lrestart,pmean,rang,timel,two,&
       &itmax,tmean,tstep,usdh,lspacendm,latcomp,l2T,lcdp
use read_val,only:ltabvois
USE var_pot, ONLY:ipotentiel
  implicit none

contains
  ! **************************************************************
  subroutine init(atdml,boxndm,celndm,psc)
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
    class(atom_config)::atdml
    type(cell_config),intent(out)::celndm
    class(box_config),intent(out)::boxndm
    type(para_space_config)::psc

    integer :: i, lufilmpaf,itapp,j,lenfn2,ipath,ierr,formatsauv
    !-----------------------------------------------
    character*2::extension
    logical :: lrepart

    tmean = 0.0
    pmean = 0.0
    timel = 0.0
    kinemean = 0.0
    lufilmpaf = 79

    !     write(6,*)'entree dans init.f'
    !potentiel BKS
#ifdef PARAPH
    rang=rangph
#endif

    call init_pot
    usdh = 1/(two*tstep)
    if (ibrake.gt.0) then
       call initeloss
    end if
    iteration=0
    !<---------setting the configuration by reading gin / cin file --------------

    if ((ipotentiel==-10).or.(ipotentiel==-11))then
       lrepart=.false.
    else
       lrepart=.true.
    end if
    call constrconf(atdml,boxndm,celndm,lrepart,psc=psc)
    if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       call caltabtC(celndm,atdml,lperiod,boxndm,psc=psc)
    else
       call caltabtC(celndm,atdml,lperiod,boxndm)
    end if
    call init_pot2(boxndm,atdml%imm)
#ifdef DECOUP
    ! Pas la peine d'aller plus loin dans l'initialisation
    return
#endif
#ifdef LAMMPS_VERSION
    if ((ipotentiel==-10).or.(ipotentiel==-11))then
       firsttime_lammps=.true.
       call init_lammps()
       if (rang==0) write(6,*)'postinitlammps'
    end if
#endif  
    !<---------setting the configuration by generation gin / cin file --------------

    select case (igen)
    case (-1)
       formatsauv = 2 ; fnamcout= fnam(1:lenfnam)//'.cout.'
       call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp=latcomp)
       if (rang==0) write (6, *) 'generation terminee'
       call arret_ndm
    case (2)
       !          call cin2gin
       call arret_ndm

    case (3)
       call transf(atdml)
       formatsauv = 2 ; fnamcout= fnam(1:lenfnam)//'.cout.'
       call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp=latcomp)
       if (rang==0) write (6, *) 'modification terminee'
       call arret_ndm
    case default
    end select
#ifdef PARA
    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
       CALL comm_space%BARRIER
       call init_voisinage(celndm,psc)

!       if (rang==0)  write(6,*) 'NOMBRE DE CELLULES FRONTIERES ASSOCIEES A CHAQUE PROCESSEUR'
!       write(6,*) 'Le proc ',myidsp,' a ',psc%nbr_proc_voisin,' processeur voisin'
    end if
#endif
    !<---------end setting the cell diviion ----------------------
    if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       call caltabtC(celndm,atdml,lperiod,boxndm,psc=psc)
    else
       call caltabtC(celndm,atdml,lperiod,boxndm)
    end if
    if (ltabvois) then
       call caltabi(atdml,celndm,boxndm)
    end if

#ifdef ML
    ! MiLaDy
    if(ipotentiel==20) then
       if (rang.eq.0) then
          write(6,*)
          write(6,*)' ML  ..... configuration MiLady '
          write(6,*)
       end if
       !This comes with MiLaDy Package
       call init_config_ml
    end if
#endif


    if (L2T.eqv..true.) then
       call readelec(celndm,boxndm)
       if (rang==0) write(6,*)'!*!*!*!*! 2T MD version =', i2t,'*!*!*!*!'
       dmtype=4
       ibrake=1
       ilangevin=1
       if((ecelec==0))then
          write(6,*) 'eccelec<>0  and l2T : STOP'
          call arret_ndm
       end if
       if ((i2T==0).and.(t_cpl.lt.0)) then
          write(6,*) 'i2T=0 t_cpl<0 and l2T : STOP'
          call arret_ndm
       end if
       if (celndm%nox.le.0 ) then
          write(6,*) 'nox noy noz MUST be defined in .din with 2T: STOP'
          call arret_ndm
       end if


    endif

    if (.not.lrestart) then
       !    if (rang==0)  

       ! input and initialization of 2T
       select type(atdml)
          class is (atom_config_d)
          call initspeed(atdml,boxndm)
       end select
       if (iterasmol>=0) then
          itapp=0
          call rasmolT (atdml,boxndm,itapp,latcomp=latcomp)
       end if
    end if
    !
    !end init the speed using Maxwell proba density-----------------
    select type(atdml)
       class is (atom_config_d)
       if ((itetimestep>0).and.(.not.lcasca)) call deftimestep(atdml,boxndm)
       if (lcasca) then
          select type(atdml)
          type is (atom_config_e)
             call initcasca(atdml,celndm,boxndm)
          class default
             write(6,*) 'lcasca and not atom_donfig_e ?'
             call arret_ndm
          end select

114       format(a3,1x,3(f10.4,1x),i5)
       end if
    end select

    if (itmax==0) stop
    if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       call caltabtC(celndm,atdml,lperiod,boxndm,psc=psc)
    else
       call caltabtC(celndm,atdml,lperiod,boxndm)
    end if
    if (ltabvois) then
       call caltabi(atdml,celndm,boxndm)
    end if
    if (dmtype==6) then
       call anapos (atdml,celndm,boxndm,iteration)
       call arret_ndm
    end if
    if (lcasca) then
       fnamcout = fnam(1:lenfnam)//'.0.cout'
       formatsauv=3
       call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp=latcomp)
    else
       fnamcout = fnam(1:lenfnam)//'.cout'
       select type(atdml)
       type is (atom_config)
          formatsauv=2
       class is (atom_config_d)
          formatsauv=3
       end select
       call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp=latcomp)
    end if
    
    if (itmax==0) call arret_ndm

    if (lcdp) call initcdp
    return
  end subroutine init

end module init_mod
