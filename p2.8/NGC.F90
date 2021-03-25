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

  use WGC_mod


  implicit none


contains
  ! *************************************************************
  subroutine  NGC  (atcgin,celcgin,boxndm,psc)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:itetemp2,imm_glob,dmtype,rang,it,itmax,mdcg_noise,&
         &angst,erg2ev,potist,im_glob,lperiod,lspacendm,latcomp,lpr,dfpred
    USE var_pot, ONLY:ntyp
    use steepestdescent_mod, only: steepestdescent,conjugategradient
#ifdef PARA
    use paraconfig,only:para_config,initparapuresp
    USE parautils,only:initcomp
    use Tpara,only:nprocspace,COMM_space
#else
    use Tpara,only:nprocspace
#endif
    USE cryst_to_cart_mod,only: cryst_to_cart

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

    !    stop

!    write(6,*)'IN NGC'

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
    typrel=1
       NCALLS=0

    call initsteep
    select case (dmtype)
    case(32)
       write(6,*)'***STEEPEST DESCENT*** beta init', betaguess
       write(6,*)' NCALLS    ENERGY(erg)         ENERGY (eV)        FORCTOT      &
&  FORMAX (eV/Ang)      ***          energy gain erg eV'
!       call initsteep
       call steepestdescent(N,R,V,F,lover)

    case(33)
       write(6,*)'***CONJUGATE GRADIENT*** beta init', betaguess
       write(6,*)' NCALLS    ENERGY(erg)         ENERGY (eV)        FORCTOT     &
&   FORMAX (eV/Ang)      ***          energy gain erg eV'
       call conjugategradient(N,R,V,F,lover)
    end select

    call back2NDM(N,R,V,F,lover)
!       write(6,*)'Post back2ndmT',lover
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
    latcomp=.true.
    itesauv=0
    call endrunT(atcgcomp,cellcgcomp,boxcg,latcomp) 

    return

  end subroutine NGC


end module NGC_mod


