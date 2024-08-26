module prog_mod
  use arps_mod,only:dmloop_arps,initarps
  USE arret_ndm_mod,only:arret_ndm
  USE init_mod,only: init
  USE calfo_mod,only: calfo
  USE endrunT_mod,only: endrunT
  USE neb_mod,only: neb
  USE gcII_mod,only: gcII
  USE analyseT_mod,only: analyseT
  USE arret_ndm_mod,only: arret_ndm
  USE controleT_mod,only: controleT
  USE neb_module,only:boxneb,init_neb0
  USE var_pot
  USE ForceMatrix_mod, only: calcFM, init_MPI_FM, pscFM,paraFM
  USE montecarlo_mod, only: montecarlo,init_montecarlo,init_mpi_mcgc!!atconf_n,cells_n,boxmcgc,init_mpi_mcgc,initNP1,pscgc,config_atom_n&
  !       &,config_atom_nplus1,config_cells_n,config_cells_nplus1,atconf_nplus1,nparapath,cells_nplus1,&
  !       &idirectionmcgc,initN,init_instyp,ins_typ,boxmcgc_p,boxmcgcpath,paramcgc,seed!,initmclpr
  USE init_simple_mod,only:init_simple
  USE boxconfig,only:box_config,box_config_lpr
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e,atom_config_arps
  USE cellconfig, only:cell_config,cell_config_arps
  USE gen_com_m, ONLY:potist,rang,sig,lspaceNDM,l2t,itmax,itloopmax,timemax,timeloopmax,iseed,&
       &lprteat,lsigat,dmtype,lax,llangevin,latcomp,imm_glob,lcdp,firsttime_lammps,lprahman,lanaposart

  use read_val,only:imm,ltabvois,rvois
  use posana,only:initanapos
  use NGC_mod,only:ngc
  use NDM_ML,only:init_config_ml
#ifdef LAMMPS_VERSION
  use lammps_util_mod,only:init_lammps


#endif

  use cdp_mod,only:creadp
  !  use one_calc_mod,only:one_calc
  use d_at_at_mod
  USE dmloop_pilot_mod,only:dmloop_pilot
  use art_mod,only:art90
  use ndm2art2ndm,only:init_mpi_art



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
    integer::im,nvois
    class(atom_config),pointer::atdml
    type(atom_config),target:: atdm
    type(atom_config_arps),target:: atdmarps
    type(atom_config_d),target:: atdmd
    type(atom_config_e),target:: atdme
    class(cell_config),pointer::celndm
    class(box_config),pointer::boxndm
    type(box_config), target:: boxs
    type(box_config_LPR), target:: boxlpr
    type(cell_config),target::cellstd
    type(cell_config_arps),target::cellarps
    type(para_space_config)::psc0
    real(double)::rv
    integer::ipp

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
    select case(dmtype)
    case(41,42)
       celndm=>cellarps
    case default
       celndm=>cellstd
    end select
    if ((lPRahman).or.((dmtype == 15).or.(dmtype==151))) then
       boxndm=>boxlpr
    else
       boxndm=>boxs
    end if

    select case(dmtype)
    case(41,42)
       atdml=>atdmarps
       if ((lax).or.(lsigat).or.(lprteat).or.(llangevin)) then
          atdmarps%lax=lax
          atdmarps%lsigat=lsigat
          atdmarps%lprteat=lprteat
          if (llangevin) then
             atdmarps%llangevin=.true.
          else
             atdmarps%llangevin=.false.
          end if
       end if
    case default
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
          case(30,32,34,33,19,35,12)
             atdml=>atdm
          case default
             atdml=>atdmd
          end select
       end if
    end select
    im=0 ; nvois=0
    atdml%imm_glob=imm
    imm_glob=imm
    atdml%ltabvois=ltabvois
    atdml%rvois=rvois
    select case(dmtype)
    case default ! ALL EXCEPT 9 (NEB) OR 15 (MCGC) or 19 (ForceMatrix) or 12 ART


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

       select type (atdml)
       type is (atom_config_arps)
          select type(celndm)
          type is(cell_config_arps)
             call initarps(atdml,celndm)
          end select
       end select
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
       type is (atom_config) !no velocity

          if (lcdp) then ! special case defect creation
             call creadp(atdml,celndm,boxndm,psc0)
          else ! no change in atom number MD
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
       class is (atom_config_d) !velocities
          select case (dmtype) ! select from dmtype 
          case(5)
             write(6,*)'loopforcetest pas NDM2020' ; stop
          case(4,8,1,21,22,23,24,88,41,42) ! some form of MD, including quenchings
             if (lcdp) then ! special case defect creation
                call creadp(atdml,celndm,boxndm,psc0)
             else ! no change in atom number MD
                itloopmax=itmax
                timeloopmax=timemax

                call dmloop_pilot(atdml,celndm,boxndm,psc0,linit=.true.)
             end if
          case(30,31) ! old CG probably does not work anymore
             if (lcdp) then
                write(6,*)'noc cdp with old CG'
                call arret_ndm
             else
                call gcII (atdml%atom_config,celndm,boxndm,psc0) ! ON PASSE LA VRAIE VARIABLE ET PAS LE POINTEUR !
             endif
          case(32,33,34) 
             if (lcdp) then
                call creadp(atdml,celndm,boxndm,psc0)! special case defect creation
             else
                call NGC(atdml,celndm,boxndm,psc0) ! new CG
             end if
          case(112)
             call d_at_at(atdml,celndm,boxndm) ! simple calculations of intzeratomic distance...
          case(111)
             call arret_ndm
          case(11) ! one iteration
             if (rang==0) write (6, *) '***** PREMIERE ET UNIQUE ITERATION  ****'
             CALL CalFo(sig,potist,atdml,celndm,boxndm,psc=psc0) !(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
             call analyseT(atdml,celndm,boxndm,psc0)
             call controleT(atdml,celndm,boxndm,psc0)
             call endrunT(atdml,celndm,boxndm,latcomp)


!!$#ifdef SUNDAE    
!!$             case (16) 
!!$                call sundae
!!$#endif
!!$
!!$#ifdef MAB    
!!$             case (17) 
!!$                call mab
!!$#endif


          case default
             write(6,*)'WTF dmtype',dmtype
          end select
          !          end select
       end select
    case(19) ! force matrix  special case of case default ! ALL EXCEPT 9 (NEB) OR 15 (MCGC) or 19 (ForceMatrix) or 12 ART
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
    case(12) ! ART calculation special case of case default ! ALL EXCEPT 9 (NEB) OR 15 (MCGC) or 19 (ForceMatrix) or 12 ART
       call init_mpi_art
       if (ltabvois) then
          rv=rvois
       else
          rv=0
       end if
       call atdml%init(im,imm,ltabvois,nvois,rvois=rv) ! initialization of the complete structure (no spatial repartition)
       call init_simple(atdml,celndm,boxndm,psc=psc0)  ! in init_simple no spatial repartition
       if (lanaposart) call initanapos(atdml,celndm,boxndm)
       call art90(atdml,celndm,boxndm,psc0)

    case(9) !NEB calculation  special case of case default ! ALL EXCEPT 9 (NEB) OR 15 (MCGC) or 19 (ForceMatrix) or 12 ART
       !#ifdef PARA

       call init_mpi_neb
       !#endif
       call init_neb0
       call neb  ! (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

    case(15,151) ! Montecarlo (path or grand canonical)  special case of case default ! ALL EXCEPT 9 (NEB) OR 15 (MCGC) or 19 (ForceMatrix) or 12 ART
       !#ifdef PARA

       call init_mpi_MCGC ! PARAPATH
       !#endif
       if (ltabvois) then
          rv=rvois
       else
          rv=0
       end if ! PARAPATH
       call init_montecarlo(boxndm,rv)

       call montecarlo

    end select
  end subroutine prog
end module prog_mod
