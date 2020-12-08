module controleT_mod
  USE endrunT_mod,only: endrunT
  USE dynalloccell
  USE tempinst_mod,only: tempinst,andersenth
  USE jqbh_mod,only: jqbh
  USE period_mod,only: period
  USE caltabi_mod,only: caltabi
  USE heat_mod,only: heat
  USE creadp_mod,only: creadp
  USE deftimestep_mod,only: deftimestep
  USE atomconfig,only:atom_config,atom_config_d!,ndm2config,config2ndm
  USE cellconfig, only:cell_config!,ndm2cellconfig,cellconfig2ndm,caltabtC
  USE boxconfig,only:box_config,periodbox!,boxconfig2ndm,ndm2boxconfig
  USE mod_para,only:myidsp
#ifdef PARA
  use mpi
  USE mod_para,only:MPI_COMM_space,ierr,NDM_MPI_REAL_DOUBLE,nprocspace
#else
  USE mod_para,only:nprocspace
#endif

  implicit none
contains
  ! ***********************************************************
  !           sous-programme controle.f
  ! ***********************************************************

  subroutine controleT(atdml,celndm,boxndm)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:dmtype,unitP,unitE, timel,tempstop, sigtot,potist,maxtcel,tempstopcel,lpkbar,angst,leev,itmax,it,&
         &itetemp,fsumstop,fpstop,itetimestep,lprtrp,sigstop,temp,timemax,cunitE,cunitP,erg2eV!itab, ltabvois, itetabvois!deltaestop,epcou,fpstop,fsumstop,ibordcou,itab,itederive,iteheat,itetabvois,indi,&
    !         &landerscou,lastcool,lcdp,ljqbh,lprtrp,ltandersen,maxtcel,nbmoye,nuandersen,sigstop,tcooling,&
    !         &tempstop,tfroi,timemax,ttol,angst,bk,cunite,cunitp,dmtype,erg2ev,iko,im,imm,it,itdes,itetemp,itetimestep,&
    !         &itmax,ldesinteg,leev,lperiod,lpkbar,ltabvois,nstepdes,potist,sigtot,tcou,temp,text,tfcou,timel,tstep,unite,!&
    !         &unitp,zl,bg,nvois,nox,noy,noz,celsize

    USE var_pot, ONLY:
    USE suivinonpbc
    USE cryst_to_cart_mod,only: cryst_to_cart
    USE notperiod_mod,only: notperiod
    USE defcdp, ONLY :itecdp
    implicit none

    class(atom_config_d)::atdml
    type(cell_config):: celndm
    type(box_config)::boxndm

    integer :: nacou, i, ic, iti
    real(double) :: vv, a1, a2, a3, c1, c2, c3
    real(double), dimension(1,3) :: g1,aux
    real(double) :: ltc, ctime, tdev, tcool, epc1, epc2, epc3,masstot,massa,tclt
    real(double), dimension(1,3) :: xtr, cv
    real(double) :: fpmax,fpn,forctot,formax,fpsmax
    !    real(double),allocatable :: xpnp(:,:)
    real(double) :: potistmean,potistdif
    real(double),save :: potist1000
    real, allocatable,save :: potiststock(:)
    logical :: latcomp
#ifdef PARA
    real(double) :: tcou_glob
    integer      :: nacou_glob
    real(double) :: fpmax_glob,fpsmax_glob
    real(double) :: forctot_glob
    real(double) :: formax_glob
#endif
    save ltc
    !-----------------------------------------------
    !
    !
    !

    if (it>=itmax) then
       if (rang==0) write (6, *) '*******Derniere iteration **** '
       latcomp=.false.
       call endrunT(atdml,celndm,boxndm,latcomp)
       write (6, *) 'predeal '
       !       call DeallocateAll

       call arret_ndm

    endif

    if (timel>=timemax) then
       if (rang==0) write (6, *) '*******max time reached **** ',timel,timemax
       latcomp=.false.
       call endrunT(atdml,celndm,boxndm,latcomp)
       !       call DeallocateAll

       call arret_ndm

    endif
    if(lEev) then
       unitE=erg2eV
       cunitE='  eV'
    else
       unitE=1.0
       cunitE=' erg'
    end if
    if(lPkbar) then
       unitP=1.0d-9
       cunitP='kbar'
    else
       unitP=1.0
       cunitP='d/cm2'
    endif

    ! temperature is down enough ?
    if (itetemp>0) then
       if (mod(it,itetemp)==0) then
          if (temp<=tempstop) then
             if (rang==0)  write (6, *) 'temperature < tempstop '
             latcomp=.false.
             call endrunT(atdml,celndm,boxndm,latcomp)

             call arret_ndm

          endif
          if (tempstopcel.gt.0) then
             if (maxtcel<=tempstopcel) then
                write (6, *) 'temperature dans toutes les cels < tempstopcel '
                latcomp=.false.
                call endrunT(atdml,celndm,boxndm,latcomp)

                call arret_ndm
             end if
          endif
       endif
    endif


    ! change in time step ? itetimestep >0
    if (itetimestep>0) then
       if (mod(it,itetimestep)==0) call deftimestep
    endif


    select case (dmtype)

    case(2,10,8)
       if ((lprtrp.EQV..false.).and.((dmtype==10).or.(dmtype==8))) goto 123
       if ((fpstop>0.0).AND.(it.GE.1)) then

          fpmax=sqrt( MAXVAL( Sum(atdml%fp(1:3,1:atdml%im)**2,1) ) )
          fpSmax = MaxVal( Abs(atdml%fp(:,1:atdml%im)) )
#ifdef PARA
          if (nprocspace.gt.1) then
             call MPI_ALLREDUCE(fpmax,fpmax_glob,1,NDM_MPI_REAL_DOUBLE,MPI_MAX,MPI_COMM_space,ierr)
             fpmax=fpmax_glob
             call MPI_ALLREDUCE(fpSmax,fpSmax_glob,1,NDM_MPI_REAL_DOUBLE,MPI_MAX,MPI_COMM_space,ierr)
             fpSmax=fpSmax_glob
          end if
#endif

          fpn=fpSmax*erg2eV/angst
          !if (rang==0)     write(6,*)
          if (myidsp==0)      write(6,'("TR: force max, energy",i6,3E20.10)') it,fpn, potist*erg2eV
          if ( myidsp==0)     write (6, *) 'energie ',potist*erg2eV
!                write(6,'("TR: force max, energy",i6,2E20.10,I3)') it,fpn, potist*erg2eV,myidsp
!               write (6, *) 'energie ',potist*erg2eV,myidsp
          !if (rang==0)     write(6,'(a,2g20.12)')'force max cgs  ev/Ang ',fpmax, fpn
          if((myidsp==0).and.(sigstop.ge.0))write(6,*)'sigma max kbar', 1d-9*maxval(abs(sigtot))
          if (fpn.le.fpstop)then
             if (sigstop.ge.0) then
                if(maxval(abs(sigtot)).le.sigstop/1d-9) then
                   latcomp=.false.
                   call endrunT(atdml,celndm,boxndm,latcomp)
                endif
             else
                latcomp=.false.
                call endrunT(atdml,celndm,boxndm,latcomp)
             end if

             if (rang==0)      write(6,*)
          end if
       end if

       if ((fsumstop>0.0).AND.(it.GE.1)) then
          !        IF (lFrozen) THEN
          !fpmax=sqrt( Sum( SUM(fp(1:3,1:im)**2,1), Free(1:im) ) )
          !           fpmax=sqrt( SUM( fp(:,1:im)**2, .NOT.Frozen(:,1:im) ) )
          !        ELSE
          fpSmax=sqrt( SUM(atdml%fp(:,1:atdml%im)**2) )
          !        END IF
#ifdef PARA
          if (nprocspace.gt.1) then
             fpSmax=fpSmax**2
             fpmax_glob=0
             call MPI_ALLREDUCE(fpmax,fpmax_glob,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
             fpSmax=sqrt(fpmax_glob)
          end if
#endif


          fpn=fpmax*erg2eV/angst
          if (myidsp==0)      write(6,*)
          if (myidsp==0)      write(6,*)'  sqrt ( sum_f F_i^2 ):  cgs  ev/Ang ',fpmax, fpn
          if((myidsp==0).and.(sigstop.ge.0))write(6,*)'sigma max kbar',1d-9* maxval(abs(sigtot))
          if (fpn.le.fsumstop) then
             if (sigstop.ge.0) then
                if(maxval(abs(sigtot)).le.sigstop/1d-9) then
                   latcomp=.false.
                   call endrunT(atdml,celndm,boxndm,latcomp)
                end if
             else
                latcomp=.false.
                call endrunT(atdml,celndm,boxndm,latcomp)
             end if
             if (myidsp==0)      write(6,*)
          end if
       end if
       !     if (rang==0) write (6, '(I10,A,D21.12,A)') it,  '*Epot = ', potist*unitE, cunitE
123    continue
       
    case(3,30) ! Gradient conjugue sur coordonnee cartesiennes (3) ou reduites (30)


       if (it==1) then
          if (lEev.EQV..true.) then
             if (myidsp==0) write(6,*)'Resultats en eV, Ang'
          else
             if (myidsp==0) write(6,*)'Resultats en cgs'
          end if
          if (myidsp==0)      write(*,'(70("="))')
          if (myidsp==0)      write(*,'("CG:     ","iter",10(" "),"epsi",14(" "),"Fmax",14(" "), "Energy")')
          if (myidsp==0)      write(*,'(70("="))')
       end if


       !debug     write(*,*) 'DEBUG ALL IT IN CONTROLE',it

       IF (it.GE.1) THEN
          !        IF (lFrozen) THEN
          !forctot=sqrt( Sum( SUM(fp(1:3,1:im)**2,1), Free(1:im) ) )
          !formax=sqrt( MAXVAL( Sum(fp(1:3,1:im)**2,1), Free(1:im) ) )
          !           forctot = sqrt( SUM( fp(:,1:im)**2, .NOT.Frozen(:,1:im) ) )
          !           formax = MaxVal( Abs(fp(:,1:im)), .NOT.Frozen(:,1:im) )
          !        ELSE
          forctot=sqrt( SUM(atdml%fp(1:3,1:atdml%im)**2) )
          !formax=sqrt( MAXVAL( Sum(fp(1:3,1:im)**2,1) ) )
          formax = MaxVal( Abs(atdml%fp(:,1:atdml%im)) )
          !        END IF
#ifdef PARA
          if (nprocspace.gt.1) then
             fpmax_glob=0
             call MPI_ALLREDUCE(formax,fpmax_glob,1,NDM_MPI_REAL_DOUBLE,MPI_MAX,MPI_COMM_space,ierr)
             formax=fpmax_glob
             forctot=forctot**2
             fpmax_glob=0
             call MPI_ALLREDUCE(forctot,fpmax_glob,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
             forctot=sqrt(fpmax_glob)
          end if
#endif

          if (lEev.EQV..true.) then
             forctot = forctot*erg2eV/angst
             formax  = formax*erg2eV/angst
             if (myidsp==0) write(*,'("GC: ",i6,3E20.10)') it,forctot, formax, potist*erg2eV
             if (fpstop>0) then
                if (formax.le.fpstop) then
                   if (myidsp==0) write(6,*)'force par atome  max  ev/Ang ', formax
                   if (myidsp==0) write (6, *) 'energie ', potist*erg2eV
                   latcomp=.false.
                   call endrunT(atdml,celndm,boxndm,latcomp)
                end if
             end if
             if (fsumstop>0) then
                if (forctot.le.fsumstop) then
                   if (myidsp==0) write(6,*)'  sqrt ( sum_f F_i^2 ):   ev/Ang ', forctot
                   if (myidsp==0) write(6, *) 'energie ', potist*erg2eV
                   latcomp=.false.
                   call endrunT(atdml,celndm,boxndm,latcomp)
                end if
             end if

          else
             if (myidsp==0) write(*,'("GC: ",i6,3E20.10)') it,forctot, formax, potist
             if (fpstop>0) then
                if (formax.le.fpstop) then
                   if (myidsp==0) write(6,*)'force par atome  max cgs ',formax
                   if (myidsp==0) write (6, *) 'energie ', potist
                   latcomp=.false.
                   call endrunT(atdml,celndm,boxndm,latcomp)

                end if
             end if

             if (fsumstop>0) then
                if (forctot.le.fsumstop) then
                   if (myidsp==0) write(6,*)'  sqrt ( sum_f F_i^2 ):  cgs  ', forctot
                   if (myidsp==0) write (6, *) 'energie ', potist
                   latcomp=.false.
                   call endrunT(atdml,celndm,boxndm,latcomp)
                end if
             end if

          end if
       end if ! it .ge.1

    case default
    end select

!!$
!!$
!!$    if ((iteheat.gt.0).and.(mod(it,iteheat)==0))call heat(im,xp,vp,ityp)
!!$



    !
    return
  end subroutine controleT
   end module controleT_mod
