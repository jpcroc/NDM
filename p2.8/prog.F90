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
       &idirectionmcgc,initN
  USE init_simple_mod,only:init_simple
  USE boxconfig,only:box_config,boxconfig2ndm,ndm2boxconfig
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE cellconfig, only:cell_config
  USE gen_com_m, ONLY:potist,rang,sig,lspaceNDM,l2t,itmax,itloopmax&
       &,lprteat,lsigat,dmtype,lax,llangevin,latcomp,imm_glob,lcdp
  
  use read_val,only:imm,ltabvois,rvois
  use NGC_mod,only:ngc
#if defined ML || defined PARAML    
  USE ml_main_mod,only: ml_main
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
    type(box_config)::boxndm
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
          case(30,32,34,33,19)
             atdml=>atdm
        case default
           atdml=>atdmd
        end select
    end if
    im=0 ; nvois=0
    atdml%imm_glob=imm
    imm_glob=imm

!!$    select type (atdml)
!!$    type is (atom_config)
!!$       write(6,*)'atomfig'
!!$    type is (atom_config_d)
!!$       write(6,*)'atomfigD'
!!$    type is (atom_config_e)
!!$       write(6,*)'atomfigE'
!!$!       write(6,*)'FLAGSFF', atdml%lprteat,atdml%lsigat,atdml%llangevin,atdml%lax
!!$    end select
       
    select case(dmtype)
    case default ! ALL EXCEPT 9 (NEB) OR 15 (MCGC) or 19 (ForceMatrix)


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
!       call atdml%print
#ifdef DECOUP
       ! Dans ce cas, pas la peine d'aller plus loin on peut terminer le programme
       return
#endif
!!$       select type (atdml)
!!$       type is (atom_config)
#ifdef PARA
       if ((dmtype.ne.30).and.(dmtype.ne.31).and.(dmtype.ne.32).and.(dmtype.ne.34).and.(dmtype.ne.33))then
          if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
             call maj_atomes_frt_ftm(atdml,celndm,boxndm,psc0)
          end if
       end if
#endif
!!$
!!$       select type (atdml)
!!$       type is (atom_config_d)
!!$          write(6,*)'typeDPROG111'
!!$          !    type is (atom_config)
!!$          !       write(6,*)'type0'
!!$       type is (atom_config_e)
!!$          write(6,*)'typeEPROG111',atdml%lprteat
!!$       end select

       select type(atdml)
       type is (atom_config)

          if (lcdp) then
             call creadp(atdml,celndm,boxndm,psc0)
          else
             itloopmax=itmax
             select case (dmtype) 
             case(30,31)
                call gcII (atdml,celndm,boxndm,psc0) ! ON PASSE LA VRAIE VARIABLE ET PAS LE POINTEUR !
             case(32,33,34)
                call NGC(atdml,celndm,boxndm,psc0)
                call endrunT(atdml,celndm,boxndm,.false.)
             end select
          end if
       class is (atom_config_d)
!!$          case default
!!$             write(6,*)'incohérence entre type(atom_config) et dmtype'
!!$             call arret_ndm
!!$          end select
!!$          class is (atom_config_d)
!!$          select case (dmtype) 
          select case (dmtype) 
          case(5)
             write(6,*)'loopforcetest pas NDM2020' ; stop
          case(4,10,8,1,21,22)
             if (lcdp) then
                call creadp(atdml,celndm,boxndm,psc0)
             else
!!$                select type (atdml)
!!$                type is (atom_config_d)
!!$                   write(6,*)'typeDPROG'
!!$                   !    type is (atom_config)
!!$                   !       write(6,*)'type0'
!!$                type is (atom_config_e)
!!$                   write(6,*)'typeEPROG',atdml%lprteat
!!$                end select
                itloopmax=itmax
                call dmloop_pilot(atdml,celndm,boxndm,psc0)
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

!!$          case(4,10)
!!$             call dmloop_vverlet (atdml,celndm,boxndm,psc0)
!!$          case(8)
!!$             call dmloop_lpr (atdml,celndm,boxndm,psc0)
!!$          case (1)
!!$             call dmloop (atdml,celndm,boxndm,psc0)
!!$          case (21)
!!$             call dmloop(atdml,celndm,boxndm,psc0)
!!$          case(22)
!!$             call dmloop_vverlet (atdml,celndm,boxndm,psc0)
!!$          case (3,30)
!!$             write(6,*)'incohérence entre type(atom_config_d) et dmtype=GC'
!!$             call arret_ndm
          case(112)
             call d_at_at(atdml,celndm,boxndm)
          case(111)
             call arret_ndm
!             if (rang==0) write (6, *) '***** PREMIERE ET UNIQUE ITERATION V2 ****'
!             CALL one_calc(atdml,celndm,boxndm,psc=psc0) 
          case(11)
             if (rang==0) write (6, *) '***** PREMIERE ET UNIQUE ITERATION  ****'
             CALL CalFo(sig,potist,atdml,celndm,boxndm,psc=psc0) !(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
             call analyseT(atdml,celndm,boxndm)
             call controleT(atdml,celndm,boxndm)
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

#if defined PHONDY || defined PARAPH    
          case (7) 
             call phondy
#endif

#if defined ML || defined PARAML    
          case (18) 
             call ml
#endif
             !          case (15)
             !             call montecarlo(atdml,celndm,boxndm)
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


       !       if (nparapath==1) then
       do ipp=1,nparapath
          atconf_n=> config_atom_n(ipp)
          cells_n=>config_cells_n(ipp)
          atconf_nplus1=>config_atom_nplus1(ipp)
          cells_nplus1=>config_cells_nplus1(ipp)
!          write(6,*)'IM',im,rang
          if (idirectionmcgc==0) then
             call atconf_n%init(im,imm_glob,ltabvois,nvois,rvois=rv)
          else
             call atconf_nplus1%init(im,imm_glob,ltabvois,nvois,rvois=rv)
          end if
          ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs
          boxmcgc=boxndm
          if (ipp==1) then
             linitpot=.true.
          else
             linitpot=.false.
          end if
          if (idirectionmcgc==0) then
             call init_simple(atconf_n,cells_n,boxmcgc,psc=pscgc,linitpot=linitpot) 
             call initNP1(ipp) ! initialise la configuration N+1
          else
             call init_simple(atconf_nplus1,cells_nplus1,boxmcgc,psc=pscgc,linitpot=linitpot) 
             call initN(ipp) ! initialise la configuration N+1
          end if
       end do
       !END PARAPATH
       atconf_n=> config_atom_n(1)
       cells_n=>config_cells_n(1)
       atconf_nplus1=>config_atom_nplus1(1)
       cells_nplus1=>config_cells_nplus1(1)
       
       call montecarlo

    end select
  end subroutine prog
end module prog_mod
