module NGC_mod
  USE arret_ndm_mod,only:arret_ndm
  USE atomconfig,only:atom_config,atom_config_d,atom_config_e
  USE cellconfig,only:cell_config
  USE boxconfig,only:box_config
  USE endrunT_mod,only:endrunT
  use Tpara,only:para_space_config
  USE recips_mod,only: calcvol
  USE Mat_utils_mod,only:  MatInv
  use constrconf_mod,only:repartition
  use WGC_mod,only:atcgcomp,atcgloc,boxcg,cellcgcomp,cellcgloc,F,ityprel,Nvar,R,pscCG,V,ncalls,betaguess,&
       &initsteep,back2ndm,final_tconv,nextsauv,nextmol,fpstop0,betaV,betaP,beta,gcpara,lchg,set_pointers_gc,&
       & unitgc,atcgmin,atcible,fpstopsig,betaV0,betaP0,ncgtry
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:itetemp2,imm_glob,dmtype,rang,iteration,mdcg_noise,iterasmol,lenfnam,fnam,&
         &angst,erg2ev,potist,lperiod,lspacendm,lprahman,dfpred,itesauv,unitP,fpstop,lcdp,fsumstop,sigstop
    USE var_pot, ONLY:ntyp
    use steepestdescent_mod, only: conjugategradient,adamrel
#ifdef PARA
    use paraconfig,only:para_config,initparapuresp
    USE parautils,only:initcomp,WORKER_TAG,tolstoi,STOP_TAG
    use Tpara,only:nprocspace,COMM_space
#else
    use Tpara,only:nprocspace
#endif
    USE cryst_to_cart_mod,only: cryst_to_cart


  implicit none


contains
  ! *************************************************************
  subroutine  NGC  (atcgin,celcgin,boxndm,psc)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------

    class(atom_config),target::atcgin
    type(cell_config),target::celcgin
    type(box_config),target::boxndm
    type(para_space_config)::psc
    logical ::lover
	integer::irel
    !-----------------------------------------------
    !
    !

#ifdef PARA


 !  integer :: im_loc



#endif
    open(unitgc, file=fnam(1:lenfnam)//'.GCout', form='formatted')

    pscCG=psc
    cellcgloc=>celcgin
    atcgloc=>atcgin
    boxcg=boxndm
    unitP=1d-9
    fpstop0=fpstop
    !    call arret_ndm

    lchg=.true.

    if (rang==0 )write(unitgc,*)'IN NGC',betaguess
!    if (rang==0 )write(6,*)'IN NGC',betaguess
    betaV=betaguess
    betaP=betaguess/3
    betaV0=betaguess
    betaP0=betaguess/3
    iteration=0

    call atcgin%deftype(atcgcomp)
    call atcgin%deftype(atcgmin)
    call atcgin%deftype(atcible)

#ifdef PARA
    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
       call atcgcomp%init(atcgin%im_glob,imm_glob,im_glob=atcgin%im_glob,rvois=atcgcomp%rvois)
       call initparapuresp(gcpara,rang,comm_space)
       call initcomp(atcgcomp,cellcgcomp,atcgin,celcgin,boxcg,gcpara,lperiod)
    else
       call initparapuresp(gcpara,rang,comm_space)
       atcgcomp=atcgin
       cellcgcomp=celcgin
    end if

#else
    atcgcomp=atcgin
    cellcgcomp=celcgin
#endif

    call set_pointers_gc ! initilisations des pointers pour tolstoi et calfo
#ifdef PARA
    if (gcpara%lmaster.neqv..true.) then
       call tolstoi (WORKER_TAG,gcpara) 
    else
#endif       

       NCALLS=0
       if (itesauv.gt.0)    nextsauv=itesauv
       if (iterasmol.gt.0)    nextsauv=iterasmol

       if (lprahman) then
          do irel=1,ncgtry
             write(6,*)
             write(unitgc,*)
             write(unitgc,*)
             ityprel=1
             fpstop=fpstop0
             beta=betaV
             call pilotcg(ityprel)
             call final_tconv(lover)
             write(unitgc,*)'loverout1',lover,irel 
             betaV=beta
             if (lover) exit
             write(6,*)
             write(unitgc,*)
             write(unitgc,*)
             ityprel=2
             beta=betaP
!             hold(:,:)=boxcg%at(:,:)


             call pilotcg(ityprel)
             call final_tconv(lover)
             betaP=beta 
             write(unitgc,*)'loverout2',lover,irel
             if (lover) exit
          end do
       else
          beta=betaV
          ityprel=1
          call pilotcg(ityprel)
          call final_tconv(lover)
       end if
       if (lover) then
          if (rang==0) then
             write(unitgc,*) 'MINIMUM REACHED after ',ncalls,' force calculations'
             write(unitgc,*) '*****************ENERGY erg eV ',Potist,Potist*erg2eV
             write(6,*) 'MINIMUM REACHED after ',ncalls,' force calculations'
             write(6,*) '*****************ENERGY erg eV ',Potist,Potist*erg2eV
          end if
       else
          if (rang==0) then
             write(unitgc,*) 'MINIMUM NOT REACHED !!!!!!!!!!!!!!!!!'
             write(unitgc,*) 'ENERGY erg eV ***',Potist,Potist*erg2eV
             write(6,*) 'MINIMUM NOT REACHED !!!!!!!!!!!!!!!!!'
             write(6,*) 'ENERGY erg eV ***',Potist,Potist*erg2eV
          end if
       end if

#ifdef PARA
       call tolstoi (STOP_TAG,gcpara) ! make servants return
    end if
#endif
    
    itesauv=0
    boxndm=boxcg
    return

  end subroutine NGC

  subroutine pilotCG (ityprel)
    integer,intent(in)::ityprel
    logical ::lover,lorig


    select case(ityprel)
    case(1)
       call initsteep
       select case (dmtype)
       case(32)
          write(unitgc,'(A,E15.8,A,2E15.8)')'******STEEPEST DESCENT*** beta init', beta, 'fpstop fsumstop ' ,fpstop,fsumstop
          write(6,'(A,E15.8,A,2E15.8)')'******STEEPEST DESCENT*** beta init', beta, 'fpstop fsumstop ' ,fpstop,fsumstop
          write(unitgc,*)' FPSTOP=', fpstop,fsumstop
             write(unitgc,'(2A)')' NCALLS      ENERGY (eV)       FORCTOT         FORMAX        ',&
                  &'   FORCESIGMA(eV/Ang)         SIGMA(kbar)   ***      energy gain eV]'
             call conjugategradient(Nvar,R,V,F,lover,lorig,betaV0) ! steepestdescent(Nvar,R,V,F,lover,betaV0)
       case(35)
          write(unitgc,'(A,E15.8,A,2E15.8)')'******ADAM RELAXATION*** beta init', beta, 'fpstop fsumstop ' ,fpstop,fsumstop
          write(6,'(A,E15.8,A,2E15.8)')'******ADAM RELAXATION*** beta init', beta, 'fpstop fsumstop ' ,fpstop,fsumstop
          write(unitgc,*)' FPSTOP=', fpstop,fsumstop
             write(unitgc,'(2A)')' NCALLS      ENERGY (eV)       FORCTOT         FORMAX        ',&
                  &'   FORCESIGMA(eV/Ang)         SIGMA(kbar)   ***      energy gain eV]'

             call Adamrel(Nvar,R,V,F,lover,betaV0) ! steepestdescent(Nvar,R,V,F,lover,betaV0)

       case(33,34)
          if (dmtype==34) then
             lorig=.false.
          else
             lorig=.true.
          end if
          write(unitgc,'(A,E20.8,A,2E20.8)')'******CG*** beta init', beta, 'fpstop fsumstop ' ,fpstop,fsumstop
!!$          if (lprahman) then
             write(unitgc,*)' FPSTOP=',fpstop,fsumstop
             write(unitgc,'(2A)')' NCALLS      ENERGY (eV)       FORCTOT         FORMAX        ',&
                  &'FORCESIGMA(eV/Ang)       SIGMA(kbar)   ***      energy gain eV]'
             write(6,'(A,E15.8,A,E12.6)')'******CONJUGATE GRADIENT *** beta init', beta
             write(6,*)' FPSTOP=',fpstop,fsumstop
             write(6,'(2A)')' NCALLS      ENERGY (eV)       FORCTOT         FORMAX        ',&
                  &'FORCESIGMA(eV/Ang)       SIGMA(kbar)   ***      energy gain eV]'
          call conjugategradient(Nvar,R,V,F,lover,lorig,betaV0)
       end select

    case(2)
       call initsteep
       select case (dmtype)
       case(32)
          write(unitgc,'(A,E15.8,A,2E12.6)')'******STEEPEST DESCENT  VARIABLE VOLUME  *** beta init', beta,&
               &' SIGSTOP=>FPSTOP=',sigstop,fpstopsig
             write(unitgc,'(2A)')' NCALLS      ENERGY (eV)       FORCTOT         FORMAX        ',&
                  &'   FORCESIGMA(eV/Ang)         SIGMA(kbar)   ***      energy gain eV]'
          write(6,'(A,E15.8,A,2E12.6)')'******STEEPEST DESCENT  VARIABLE VOLUME  *** beta init', beta,&
               &' SIGSTOP=>FPSTOP=',sigstop,fpstopsig
             write(6,'(2A)')' NCALLS      ENERGY (eV)       FORCTOT         FORMAX        ',&
                  &'    FORCESIGMA(eV/Ang)           SIGMA(kbar)   ***      energy gain eV]'
          !       call initsteep
          call conjugategradient(Nvar,R,V,F,lover,lorig,betaV0)

       case(33,34)
          if (dmtype==34) then
             lorig=.false.
          else
             lorig=.true.
          end if
          write(unitgc,'(A,E15.8,A,E12.6)')'******CONJUGATE GRADIENT VARIABLE VOLUME *** beta init', beta
          write(unitgc,*)' SIGSTOP=>FPSTOP=',sigstop, fpstopsig
             write(unitgc,'(2A)')' NCALLS      ENERGY (eV)       FORCTOT         FORMAX        ',&
                  &'FORCESIGMA(eV/Ang)       SIGMA(kbar)   ***      energy gain eV]'
          write(6,'(A,E15.8,A,E12.6)')'******CONJUGATE GRADIENT VARIABLE VOLUME *** beta init', beta
          write(6,*)' SIGSTOP=>FPSTOP=',sigstop, fpstopsig
             write(6,'(2A)')' NCALLS      ENERGY (eV)       FORCTOT         FORMAX        ',&
                  &'    FORCESIGMA(eV/Ang)         SIGMA(kbar)   ***      energy gain eV]'
!                 write(unitgc,*)'SIGSTOP==FPSTOP=',fpstop
             call conjugategradient(Nvar,R,V,F,lover,lorig,betaP0)
          case(35)
             
          write(unitgc,'(A,E15.8,A,2E15.8)')'******ADAM RELAXATION*** beta init', beta, 'fpstop fsumstop ' ,fpstop,fsumstop
          write(6,'(A,E15.8,A,2E15.8)')'******ADAM RELAXATION*** beta init', beta, 'fpstop fsumstop ' ,fpstop,fsumstop
          write(unitgc,*)' FPSTOP=', fpstop,fsumstop
             write(unitgc,'(2A)')' NCALLS      ENERGY (eV)       FORCTOT         FORMAX        ',&
                  &'   FORCESIGMA(eV/Ang)         SIGMA(kbar)   ***      energy gain eV]'

             call Adamrel(Nvar,R,V,F,lover,betaV0) ! steepestdescent(Nvar,R,V,F,lover,betaV0)

       end select
     
    end select
!    write(unitgc,*)
    call back2NDM(Nvar,R,V,F,lover)

    !       write(unitgc,*)'Post back2ndmT',lover


  end subroutine pilotCG

end module NGC_mod


