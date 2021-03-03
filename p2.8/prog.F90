module prog_mod
  USE init_mod,only: init
  USE calfo_mod,only: calfo
  USE analyse_mod,only: analyse
  USE controle_mod,only: controle
  USE endrunT_mod,only: endrunT
  USE neb_mod,only: neb
  USE dmloop_lpr_mod,only: dmloop_lpr
  USE loopforcetest_mod,only: loopforcetest
  USE gcII_mod,only: gcII
  USE dmloop_vverlet_mod,only: dmloop_vverlet
  USE dmloop_mod,only: dmloop
  USE analyseT_mod,only: analyseT
  USE controleT_mod,only: controleT
  USE neb_module,only:boxneb,init_neb0
  USE var_pot
  USE montecarlo_mod, only: montecarlo,config_atom_n,cells_n,boxmcgc,init_mpi_mcgc,initNP1,pscgc
  USE init_simple_mod,only:init_simple
  USE boxconfig,only:box_config,boxconfig2ndm,ndm2boxconfig
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e,ndm2config, config2ndm
  USE cellconfig, only:cell_config,ndm2cellconfig,cellconfig2ndm
  USE gen_com_m, ONLY:parallele,potist,rang,sig,lspaceNDM&
       &,lprteat,lsigat,imm_glob,dmtype,imm_glob,lax,llangevin,latcomp
  
  use read_val,only:imm,ltabvois,rvois

#if defined ML || defined PARAML    
  USE ml_main_mod,only: ml_main
#endif
  USE cellconfig, only:cell_config,ndm2cellconfig,cellconfig2ndm

  implicit none
contains
  subroutine prog
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:parallele,potist,rang,sig,lspaceNDM&
         &,lprteat,lsigat,imm_glob,dmtype,imm_glob,lax,llangevin

    use read_val,only:imm,ltabvois,rvois

    !    USE montecarlo_mod, ONLY: config_atom_n, cells_n

#ifdef PARA
!    use mpi
    USE Tpara,only:MPI_COMM_space,myidsp,nprocspace,para_space_config
    USE mod_para,only:maj_atomes_frt_ftm
    USE neb_module,only:init_mpi_neb
#else
    USE Tpara,only:nprocspace,para_space_config
#endif
    USE neb_module,only:init_mpi_neb
    implicit none
#ifdef PARA
 include 'mpif.h'
#endif
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
    if ((lax).or.(lsigat).or.(lprteat).or.(llangevin))then
       atdml=>atdme
    else
       if ((dmtype==3).or.(dmtype==30)) then
          atdml=>atdm
       else

          atdml=>atdmd
       end if
    end if
    im=0 ; nvois=0
    imm_glob=imm
    select case(dmtype)
    case default
       !if (dmtype.ne.9) then

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
       if (ltabvois) then
          rv=rvois
       else
          rv=0
       end if

       call atdml%init(im,imm,ltabvois,nvois,rvois=rv,lsigat=lsigat,lprteat=lprteat,llangevin=llangevin,lax=lax)
       ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs
       call init(atdml,boxndm,celndm,psc0)

#ifdef DECOUP
       ! Dans ce cas, pas la peine d'aller plus loin on peut terminer le programme
       return
#endif
       select type (atdml)
       type is (atom_config)
          select case (dmtype) 
          case(3,30)
             call gcII (atdm,celndm,boxndm,psc0) ! ON PASSE LA VRAIE VARIABLE ET PAS LE POINTEUR !
          case default
             write(6,*)'incohérence entre type(atom_config) et dmtype'
             stop
          end select
          class is (atom_config_d)
#ifdef PARA
          if ((dmtype.ne.3).and.(dmtype.ne.30))then
             if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
                call maj_atomes_frt_ftm(atdml,celndm,psc0)
             end if
          end if
#endif
          select case (dmtype) 
          case(5)
             write(6,*)'loopforcetest pas NDM2020' ; stop
             !       if (.not.parallele)    call loopforcetest (xp, xpp, vp, ax, fp, ielat, iwmax, ityp,num_at_glob)
          case(4,10)
             call dmloop_vverlet (atdml,celndm,boxndm,psc0)
          case(8)
             call dmloop_lpr (atdml,celndm,boxndm,psc0)
          case (1)
             if (.not.parallele)  call dmloop (atdml,celndm,boxndm,psc0)
          case (2)
             if (.not.parallele)  then
                call dmloop(atdml,celndm,boxndm,psc0)
             else
#ifdef PARA
                if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
                   if(rang==0) write (6,*)'DMTYPE 2 +PARA=DMLOOP_VVERLET_+OPTION'
                   call dmloop_vverlet (atdml,celndm,boxndm,psc0)
                end if
#endif
             endif
          case (3,30)
             write(6,*)'incohérence entre type(atom_config_d) et dmtype=GC'
             stop

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
          end select
       end select
    case(9)
       !#ifdef PARA
       call init_mpi_neb
       !#endif
       call init_neb0 
       call neb  ! (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
       
    case(15)
       !#ifdef PARA
       call init_mpi_MCGC
       !#endif
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
       if (ltabvois) then
          rv=rvois
       else
          rv=0
       end if
       call config_atom_n%init(im,imm,ltabvois,nvois,rvois=rv,lsigat=lsigat,lprteat=lprteat,llangevin=llangevin,lax=lax)
       ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs
       boxmcgc=boxndm

       call init_simple(config_atom_n,cells_n,boxmcgc,psc=pscgc) 
       call initNP1 ! initialise la configuration N+1 
       call montecarlo
       
    end select
  end subroutine prog
end module prog_mod
