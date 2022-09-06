module prog_mod
  USE arret_ndm_mod,only:arret_ndm
  USE init_mod,only: init
  USE calfo_mod,only: calfo
  USE endrunT_mod,only: endrunT
  USE neb_mod,only: neb
!  USE dmloop_lpr_mod,only: dmloop_lpr
  USE gcII_mod,only: gcII
!  USE dmloop_vverlet_mod,only: dmloop_vverlet
!  USE dmloop_mod,only: dmloop
  USE analyseT_mod,only: analyseT
  USE arret_ndm_mod,only: arret_ndm
  USE controleT_mod,only: controleT
  USE neb_module,only:boxneb,init_neb0
  USE var_pot
  USE ForceMatrix_mod, only: calcFM, init_MPI_FM, pscFM,paraFM
  USE montecarlo_mod, only: montecarlo,atconf_n,cells_n,boxmcgc,init_mpi_mcgc,initNP1,pscgc,config_atom_n&
       &,config_atom_nplus1,config_cells_n,config_cells_nplus1,atconf_nplus1,nparapath,cells_nplus1,&
       &idirectionmcgc,initN,init_instyp,ins_typ,boxmcgc_p,boxmcgcpath,initmclpr
  USE init_simple_mod,only:init_simple
  USE boxconfig,only:box_config,box_config_lpr
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE cellconfig, only:cell_config
  USE gen_com_m, ONLY:potist,rang,sig,lspaceNDM,l2t,itmax,itloopmax,timemax,timeloopmax&
       &,lprteat,lsigat,dmtype,lax,llangevin,latcomp,imm_glob,lcdp,firsttime_lammps,lprahman
  
  use read_val,only:imm,ltabvois,rvois
  use NGC_mod,only:ngc
  use NDM_ML,only:init_config_ml
!!$#if defined ML || defined PARAML    
!!$  USE ml_main_mod,only: ml_main
!!$#endif
  USE Parrinello_Rahman,only:initlpr
#ifdef LAMMPS_VERSION
  use lammps_util_mod,only:init_lammps


#endif

  use cdp_mod,only:creadp
!  use one_calc_mod,only:one_calc
  use d_at_at_mod
  USE dmloop_pilot_mod,only:dmloop_pilot
  implicit none
contains
  subroutine prog
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:potist,rang,sig,lspaceNDM&
         &,lprteat,lsigat,dmtype,lax,llangevin,itetimestep

    use read_val,only:imm,ltabvois,rvois

#ifdef PARA
    USE Tpara,only:nprocspace,para_space_config
    USE mod_para,only:maj_atomes_frt_ftm
    USE neb_module,only:init_mpi_neb
#else
    USE Tpara,only:nprocspace,para_space_config
#endif
    USE neb_module,only:init_mpi_neb
    implicit none
    character :: extension*2
    integer::lenfn2,i,ko,im,nvois
    class(atom_config),pointer::atdml
    type(atom_config),target:: atdm
    type(atom_config_d),target:: atdmd
    type(atom_config_e),target:: atdme
    type(cell_config)::celndm
    class(box_config),pointer::boxndm
    type(box_config), target:: boxs
    type(box_config_LPR), target:: boxlpr
    type(para_space_config)::psc0
    real(double)::rv
    integer::ipp
    logical::linitpot
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e rs
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------

    ! Allocation des tableaux dimensionnes sur le nombre d'atomes
    !probablement inutile pour dmtype=9 ou 15
! choose actual data types for atmdl and boxndm depending on values read in readdm
    if (lPRahman) then
       boxndm=>boxlpr
    else
       boxndm=>boxs
    end if
    
    if ((lax).or.(lsigat).or.(lprteat).or.(llangevin).or.(l2t))then
       atdml=>atdme
       atdme%lax=lax
       atdme%lsigat=lsigat
       atdme%lprteat=lprteat
       atdme%llangevin=.false.
       if ((llangevin).or.(l2t)) atdme%llangevin=.true.
    elseif(itetimestep.gt.0) then
       atdml=>atdmd
    else
       select case(dmtype)
          case(30,32,34,33,19,35)
             atdml=>atdm
        case default
           atdml=>atdmd
        end select
    end if
    im=0 ; nvois=0
    atdml%imm_glob=imm
    imm_glob=imm
    atdml%ltabvois=ltabvois
       
    select case(dmtype)
    case default ! ALL EXCEPT 9 (NEB) OR 15 (MCGC) or 19 (ForceMatrix)


#ifdef PARA
       ! En parallle, on initialise le nombre maximum d'atomes d'un
       ! processus au nombre d'atomes locaux. Plus tard ce nombre sera
       ! complete par le nombre maximal d'atomes fantomes
       ! On suppose que la concentration max ne depasse pas 20%  de 
       ! la concentration moyenne
       if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
          imm      = min( imm_glob, int(1.2 * imm_glob / nprocspace) )
          if (rang==0) write(6,*)'IMM PARA = ',imm,imm_glob
       endif
#endif
       call init(atdml,boxndm,celndm,psc0)

#ifdef DECOUP
       ! Dans ce cas, pas la peine d'aller plus loin on peut terminer le programme
       return
#endif
#ifdef PARA
       if ((dmtype.ne.35).and.(dmtype.ne.30).and.(dmtype.ne.31).and.(dmtype.ne.32)&
            &.and.(dmtype.ne.34).and.(dmtype.ne.33))then
          if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
             call maj_atomes_frt_ftm(atdml,celndm,boxndm,psc0)
          end if
       end if
#endif

       select type(atdml)
       type is (atom_config)

          if (lcdp) then
             call creadp(atdml,celndm,boxndm,psc0)
          else
             itloopmax=itmax
             timeloopmax=timemax
!             write(6,*)'TIMELOOPMAX ITLOOPMAX',timeloopmax,itloopmax
             select case (dmtype) 
             case(30,31)
                call gcII (atdml,celndm,boxndm,psc0) ! ON PASSE LA VRAIE VARIABLE ET PAS LE POINTEUR !
             case(32,33,34,35)
                call NGC(atdml,celndm,boxndm,psc0)
                call endrunT(atdml,celndm,boxndm,latcomp)
             end select
          end if
       class is (atom_config_d)
          select case (dmtype) 
          case(5)
             write(6,*)'loopforcetest pas NDM2020' ; stop
          case(4,8,1,21,22,23,24,88)
             if (lcdp) then
                call creadp(atdml,celndm,boxndm,psc0)
             else
                itloopmax=itmax
                timeloopmax=timemax
                call dmloop_pilot(atdml,celndm,boxndm,psc0,linit=.true.)
             end if
          case(30,31)
             if (lcdp) then
                write(6,*)'noc cdp with old GC'
                call arret_ndm
             else
                call gcII (atdml%atom_config,celndm,boxndm,psc0) ! ON PASSE LA VRAIE VARIABLE ET PAS LE POINTEUR !
             endif
          case(32,33,34)
             if (lcdp) then
                call creadp(atdml,celndm,boxndm,psc0)
             else
                call NGC(atdml,celndm,boxndm,psc0)
             end if
          case(112)
             call d_at_at(atdml,celndm,boxndm)
          case(111)
             call arret_ndm
          case(11)
             if (rang==0) write (6, *) '***** PREMIERE ET UNIQUE ITERATION  ****'
             CALL CalFo(sig,potist,atdml,celndm,boxndm,psc=psc0) !(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
             call analyseT(atdml,celndm,boxndm,psc0)
             call controleT(atdml,celndm,boxndm,psc0)
             call endrunT(atdml,celndm,boxndm,latcomp)

#ifdef ART    
          case (12) 
             call art90
#endif

#ifdef SUNDAE    
          case (16) 
             call sundae
#endif

#ifdef MAB    
          case (17) 
             call mab
#endif


!#if defined ML || defined PARAML    ML est considéré comme un potentiel pas un DMTYPE, A CHANGER ???
!          case (18) 
!             call ml
!#endif
          case default
             write(6,*)'WTF dmtype',dmtype
          end select


       end select
    case(19) ! force matrix
       call init_mpi_FM
       if (ltabvois) then
          rv=rvois
       else
          rv=0
       end if
       call atdml%init(im,imm,ltabvois,nvois,rvois=rv)
       call init_simple(atdml,celndm,boxndm,psc=pscFM)
#ifdef LAMMPS_VERSION

    if ((ipotentiel==-10).or.(ipotentiel==-11))then
       firsttime_lammps=.true.
       call init_lammps()

    end if
#endif
#ifdef ML

    if (ipotentiel==20)    call init_config_ml


#endif

       call calcFM(atdml,celndm,boxndm)
       call arret_ndm
    case(9)
       !#ifdef PARA
       call init_mpi_neb
       !#endif
       call init_neb0 
       call neb  ! (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

    case(15)
       !#ifdef PARA
       call init_mpi_MCGC ! PARAPATH
       !#endif
#ifdef PARA
       ! En parallle, on initialise le nombre maximum d'atomes d'un
       ! processus au nombre d'atomes locaux. Plus tard ce nombre sera
       ! complete par le nombre maximal d'atomes fantomes
       ! On suppose que la concentration max ne depasse pas 20%  de 
       ! la concentration moyenne
       if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
          imm      = min( imm_glob, int(1.2 * imm_glob / nprocspace) )
          if (rang==0) write(6,*)'IMM PARA MCGC = ',imm,imm_glob
       endif
#endif
!!$ call MPI_FINALIZE(imm)
!!$ stop
       if (ltabvois) then
          rv=rvois
       else
          rv=0
       end if ! PARAPATH

       allocate (config_atom_n(nparapath))
       allocate (config_atom_nplus1(nparapath))
       allocate (config_cells_n(nparapath))
       allocate (config_cells_nplus1(nparapath))
       allocate(boxmcgcpath(nparapath))
       if (ins_typ==1) call init_instyp
       !       if (nparapath==1) then
       select type (boxndm)
       type is (box_config_lpr)
          boxmcgc=boxndm
       end select
       do ipp=1,nparapath
          boxmcgc_p=>boxmcgcpath(ipp)
          atconf_n=> config_atom_n(ipp)
          cells_n=>config_cells_n(ipp)
          atconf_nplus1=>config_atom_nplus1(ipp)
          cells_nplus1=>config_cells_nplus1(ipp)
          boxmcgc_p=boxmcgc
!          write(6,*)'IM',im,rang
          if (idirectionmcgc==0) then
             call atconf_n%init(im,imm_glob,ltabvois,nvois,rvois=rv)
          else
             call atconf_nplus1%init(im,imm_glob,ltabvois,nvois,rvois=rv)
          end if
          ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs

          if (ipp==1) then
             linitpot=.true.
          else
             linitpot=.false.
          end if
          if (idirectionmcgc==0) then
             call init_simple(atconf_n%atom_config_d,cells_n,boxmcgc_p,psc=pscgc,linitpot=linitpot) 
             call initNP1(ipp) ! initialise la configuration N+1
          else
             call init_simple(atconf_nplus1%atom_config_d,cells_nplus1,boxmcgc_p,psc=pscgc,linitpot=linitpot) 
             call initN(ipp) ! initialise la configuration N+1
          end if
       end do
       if (lprahman) then
          call initlpr(atconf_nplus1,cells_nplus1,boxmcgc_p,pscgc)
          call initMClpr(atconf_nplus1%im)
       end if
       !END PARAPATH
       atconf_n=> config_atom_n(1)
       cells_n=>config_cells_n(1)
       atconf_nplus1=>config_atom_nplus1(1)
       cells_nplus1=>config_cells_nplus1(1)

       call montecarlo

    end select
  end subroutine prog
end module prog_mod
