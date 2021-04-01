module NGC_mod
  USE endrun_mod,only: endrun
  USE analyse_mod,only: analyse
  USE atomconfig,only:atom_config
  USE cellconfig,only:cell_config
  USE boxconfig,only:box_config
  USE endrunT_mod,only:endrunT
  use Tpara,only:para_space_config
  USE recips_mod,only: calcvol
  USE Mat_utils_mod,only:  MatInv

  use WGC_mod,only:atcgcomp,atcgloc,boxcg,cellcgcomp,cellcgloc,F,ityprel,N,R,pscCG,V,ncalls,betaguess,&
       &initsteep,back2ndm,final_tconv,nextsauv,nextmol,fpstop0,betaV,betaP,beta,gcpara,lchg,set_pointers_gc
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:itetemp2,imm_glob,dmtype,rang,it,itmax,mdcg_noise,iterasmol,&
         &angst,erg2ev,potist,im_glob,lperiod,lspacendm,latcomp,lpr,dfpred,itesauv,unitP,fpstop
    USE var_pot, ONLY:ntyp
    use steepestdescent_mod, only: steepestdescent,conjugategradient
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

    type(atom_config),target::atcgin
    type(cell_config),target::celcgin
    type(box_config)::boxndm
    type(para_space_config)::psc
    logical ::lover

    !-----------------------------------------------
    !
    !

#ifdef PARA
    integer :: iproc
    integer, allocatable      :: num_at_glob_all(:)
    !  integer :: im_loc
    integer :: proc_source


#endif
    pscCG=psc
    cellcgloc=>celcgin
    atcgloc=>atcgin
    boxcg=boxndm
    unitP=1d-9
    fpstop0=fpstop
    !    stop

    lchg=.true.

    if (rang==0 )write(6,*)'IN NGC',betaguess
    betaV=betaguess
    betaP=betaguess/3
    it=0
#ifdef PARA
    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
       call atcgcomp%init(im_glob,imm_glob)
       call initparapuresp(gcpara,rang,comm_space,nprocspace)
       call initcomp(atcgcomp,cellcgcomp,atcgin,celcgin,boxcg,gcpara,lperiod)
    else
       call initparapuresp(gcpara,rang,comm_space,nprocspace)
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
       call tolstoi (WORKER_TAG,gcpara,'xft') 
    else
#endif       

       NCALLS=0
       if (itesauv.gt.0)    nextsauv=itesauv
       if (iterasmol.gt.0)    nextsauv=iterasmol

       if (lpr) then
          do it=1,10
             ityprel=1
             beta=betaV
             call pilotcg(ityprel)
             call final_tconv(lover)
             betaV=beta
             if (lover) exit
            ityprel=2
             beta=betaP
             call pilotcg(ityprel)
             call final_tconv(lover)
             betaP=beta
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
             write(6,*) 'MIMIMUM REACHED after ',ncalls,' force calculations'
             write(6,*) '*****************ENERGY erg eV ',Potist,Potist*erg2eV
          end if
       else
          if (rang==0) then
             write(6,*) 'MIMIMUM NOT REACHED !!!!!!!!!!!!!!!!!'
             write(6,*) 'ENERGY erg eV ***',Potist,Potist*erg2eV
          end if
       end if

#ifdef PARA
       call tolstoi (STOP_TAG,gcpara,'xft') ! make servants return
    end if
#endif
    
    latcomp=.true.
    itesauv=0
!    call atcgcomp%print(unit=100+rang)
    call endrunT(atcgcomp,cellcgcomp,boxcg,latcomp) 

    return

  end subroutine NGC

  subroutine pilotCG (ityprel)
    integer,intent(in)::ityprel
    logical ::lover,lorig
    
    select case(ityprel)
    case(2)


       call initsteep
       select case (dmtype)
       case(32)
          write(6,'(A,E15.8,A,E12.6)')'******STEEPEST DESCENT  VARIABLE VOLUME  *** beta init', beta,&
               &' SIGSTOP=>FPSTOP=',fpstop
          write(6,*)' NCALLS    ENERGY(erg)         ENERGY (eV)        FORCTOT/SIG      &
               &      SIGMA      ***          energy gain erg eV'
          !       call initsteep
!          write(6,*)'SIGSTOP==FPSTOP=',fpstop
          call steepestdescent(N,R,V,F,lover)

       case(33,34)
          if (dmtype==34) then
             lorig=.false.
          else
             lorig=.true.
          end if
          write(6,'(A,E15.8,A,E12.6)')'******CONJUGATE GRADIENT VARIABLE VOLUME *** beta init', beta,&
               &' SIGSTOP=>FPSTOP=',fpstop
!          write(6,'(A,E15.8)')'******CONJUGATE GRADIENT VARIABLE VOLUME *** beta init', beta
          write(6,*)' NCALLS    ENERGY(erg)         ENERGY (eV)        FORCTOT/SIG     &
               &    SIGMA        ***          energy gain erg eV'
!                 write(6,*)'SIGSTOP==FPSTOP=',fpstop
          call conjugategradient(N,R,V,F,lover,lorig)
       end select
     
    case(1)
       call initsteep

       select case (dmtype)
       case(32)
          write(6,'(A,E15.8)')'******STEEPEST DESCENT*** beta init', beta
          write(6,*)' NCALLS    ENERGY(erg)         ENERGY (eV)        FORCTOT      &
               &  FORMAX (eV/Ang)      ***          energy gain erg eV'
          !       call initsteep
          call steepestdescent(N,R,V,F,lover)

       case(33,34)
          if (dmtype==34) then
             lorig=.false.
          else
             lorig=.true.
          end if
          write(6,'(A,E15.8)')'******CONJUGATE GRADIENT*** beta init', beta
          write(6,*)' NCALLS    ENERGY(erg)         ENERGY (eV)        FORCTOT     &
               &   FORMAX (eV/Ang)      ***          energy gain erg eV'
          call conjugategradient(N,R,V,F,lover,lorig)
       end select
    end select
!    write(6,*)
    call back2NDM(N,R,V,F,lover)
    !       write(6,*)'Post back2ndmT',lover


  end subroutine pilotCG

end module NGC_mod


