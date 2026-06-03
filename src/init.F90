module init_mod
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
  USE constrconf_mod, only :constrconf,lprt
  USE arret_ndm_mod,only: arret_ndm
  use vect_dist_mod,only:testclose
#ifdef PARA
  USE init_vois_mod,only: init_voisinage
#endif
#ifdef LAMMPS_VERSION
  use lammps_util_mod
  use vars_lammps
#endif
#ifdef ML
  use mld_interface_mod, only: mld_init_config
#endif
  USE montecarlo_mod, ONLY: distminat

  
  use Tpara,only:para_space_config

  USE gen_com_m, only:uwrt,lwrt,fnam,lenfnam,dmtype,fnamcout,igen,ilangevin,iteration,iteanapos,iterasmol,&
       &itetimestep,kinemean,lcasca,lperiod,lrestart,pmean,rang,timel,two,&
       &itmax,tmean,tstep,usdh,lspacendm,latcomp,l2T,lcdp,lspecialinit,lwgin,iattcl
  use read_val,only:ltabvois
  use specialinit_mod,only:specialinit
USE var_pot, ONLY:ipotentiel
  implicit none

contains
  ! **************************************************************
  subroutine init(atdml,boxndm,celndm,psc)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    use posana,only:anapos,  initanapos
    USE posana,only:
    USE elec_cell,ONLY: i2t,t_cpl, readelec
    USE eloss, ONLY : ibrake,ecelec,initeloss

#ifdef PARA
    USE Tpara,only:COMM_space,nprocspace
#else
    use Tpara,only:nprocspace

#endif

    ! **************************************************************

    implicit none
    class(atom_config)::atdml
    type(cell_config),intent(out)::celndm
    class(box_config),intent(out)::boxndm
    type(para_space_config)::psc

    integer :: lufilmpaf,itapp,formatsauv,i
    !-----------------------------------------------
    real(double)::xptclos(3)
    logical :: lrepart

    tmean = 0.0
    pmean = 0.0
    timel = 0.0
    kinemean = 0.0
    lufilmpaf = 79

    !     write(uwrt,*)'entree dans init.f'
    !potentiel BKS
#ifdef PARAPH
    rang=rangph
#endif

    call init_pot  ! contains calls to MLD
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
       call caltabtC(celndm,atdml,lperiod,boxndm,psc=psc,lchktrav=.true.)
    else
       call caltabtC(celndm,atdml,lperiod,boxndm,lchktrav=.false.)
    end if
!    call celndm%print
    call init_pot2(boxndm,atdml%imm)
#ifdef DECOUP
    ! Pas la peine d'aller plus loin dans l'initialisation
    return
#endif
#ifdef LAMMPS_VERSION
    if ((ipotentiel==-10).or.(ipotentiel==-11))then
       firsttime_lammps=.true.
       call init_lammps()
       if (rang==0) write(uwrt,*)'postinitlammps'
    end if
#endif  
#ifdef ML
    ! MiLaDy
    if(ipotentiel==20) then
      if (rang.eq.0) then
           write(uwrt,*)
           write(uwrt,*)' ML  ..... configuration MiLady '
           write(uwrt,*)
      end if
      !  !This comes with MiLaDy Package

      call mld_init_config(atdml) 

    !call init ! mld init
    end if
#endif
    !<---------setting the configuration by generation gin / cin file --------------

#ifdef PARA


    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
       CALL comm_space%BARRIER
       call init_voisinage(celndm,psc,lwrite=.true.)

!       if (rang==0)  write(uwrt,*) 'NOMBRE DE CELLULES FRONTIERES ASSOCIEES A CHAQUE PROCESSEUR'
!       write(uwrt,*) 'Le proc ',myidsp,' a ',psc%nbr_proc_voisin,' processeur voisin'
    end if
#endif
    !<---------end setting the cell diviion ----------------------
    if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       call caltabtC(celndm,atdml,lperiod,boxndm,psc=psc,lchktrav=.true.)
    else
       call caltabtC(celndm,atdml,lperiod,boxndm,lchktrav=.false.)
    end if
    if (ltabvois) then
       call caltabi(atdml,celndm,boxndm)
    end if


    if (L2T.eqv..true.) then
       call readelec(celndm,boxndm)
       if (rang==0) write(uwrt,*)'!*!*!*!*! 2T MD version =', i2t,'*!*!*!*!'
       dmtype=4
       if ((ibrake.ne.1).and.(ibrake.ne.3)) then
          if(rang==0) write(uwrt,*) 'ibrake1 or 3 for l2T'
          call arret_ndm
       end if
       ilangevin=1
       if((ecelec==0))then
          write(uwrt,*) 'eccelec<>0  and l2T : STOP'
          call arret_ndm
       end if
       if ((i2T==0).and.(t_cpl.lt.0)) then
          write(uwrt,*) 'i2T=0 t_cpl<0 and l2T : STOP'
          call arret_ndm
       end if
       if (celndm%nox(1).le.0 ) then
          write(uwrt,*) 'nox noy noz MUST be defined in .din with 2T: STOP'
          call arret_ndm
       end if


    endif

    if (.not.lrestart) then
       !    if (rang==0)  

       ! input and initialization of 2T
       select type(atdml)
          class is (atom_config_d)
          call initspeed(atdml,boxndm,lprt=lprt)
       end select
    end if
    select case (igen)
    case (-1)
       formatsauv = 2 ; fnamcout= fnam(1:lenfnam)//'.cout.'
       call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp=latcomp)
       call rasmolT (atdml,boxndm,-1,latcomp=latcomp,ivisumol=5)
       if (lwgin) call rasmolT (atdml,boxndm,latcomp=latcomp)
       if (rang==0) write (uwrt, *) 'generation terminee'
       call arret_ndm
    case (2)
       !          call cin2gin
       call arret_ndm

    case (3)
       call transf(atdml)
       formatsauv = 2 ; fnamcout= fnam(1:lenfnam)//'.cout.'
       call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp=latcomp)
       if (lwgin) call rasmolT (atdml,boxndm,-1,latcomp=latcomp,ivisumol=5)
       if (rang==0) write (uwrt, *) 'modification terminee'
       call arret_ndm
    case default
    end select

    !
    !end init the speed using Maxwell proba density-----------------
    select type(atdml)
    class is (atom_config_d)
       if ((itetimestep>0).and.(.not.lcasca)) call deftimestep(atdml,boxndm)
    end select
    if ((lspecialinit).and.(.not.(lrestart))) call specialinit(atdml,boxndm,celndm)
    if (.not.lrestart) then
       if (iterasmol>=0) then
          itapp=-1
          call rasmolT (atdml,boxndm,itapp,latcomp=latcomp)
          if (lwgin) call rasmolT (atdml,boxndm,itapp,latcomp=latcomp,ivisumol=5)
       end if
    end if

    if (itmax==0) stop
    if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       call caltabtC(celndm,atdml,lperiod,boxndm,psc=psc,lchktrav=.true.)
    else
       call caltabtC(celndm,atdml,lperiod,boxndm,lchktrav=.false.)
    end if
    if (ltabvois) then
       call caltabi(atdml,celndm,boxndm)
    end if
    if (dmtype==6) then
       call initanapos(atdml,celndm,boxndm)
       call anapos (atdml,celndm,boxndm,iteration)
       call arret_ndm
    end if
    if (lcasca) then
       fnamcout = fnam(1:lenfnam)//'.0.cout'
       formatsauv=5
       call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp=latcomp)
    else
       fnamcout = fnam(1:lenfnam)//'.cout'
       select type(atdml)
       type is (atom_config)
          formatsauv=2
       class is (atom_config_d)
          formatsauv=5
       end select
       call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp=latcomp)
    end if
    if (itmax==0) call arret_ndm

    if (lcdp) call initcdp

    if(iattcl.gt.0) then
       xptclos(1:3)=0.
       do i=1,atdml%im
          if (atdml%num_at_glob(i)==iattcl) then 
             xptclos(1:3)=atdml%xp(:,i)
          end if
       end do
       call comm_space%sum(xptclos)
       call testclose ( atdml,boxndm,celndm,lperiod,distminat,iattcl,xptclos)
    end if
    return
  end subroutine init

end module init_mod
