module controleT_mod
  USE arret_ndm_mod,only:arret_ndm
  USE endrunT_mod,only: endrunT
  USE deftimestep_mod,only: deftimestep
  USE atomconfig,only:atom_config,atom_config_d
  USE cellconfig, only:cell_config
  USE boxconfig,only:box_config
  USE analyseT_mod,only: analyseT
#ifdef PARA
  USE Tpara,only:COMM_space,nprocspace,myidsp
#else
  USE Tpara,only:nprocspace,myidsp
#endif
 USE arret_ndm_mod,only: arret_ndm
 use Tpara,only:para_space_config    
  implicit none
contains
  ! ***********************************************************
  !           sous-programme controle.f
  ! ***********************************************************

  subroutine controleT(atdml,celndm,boxndm,psc,lreturn)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:dmtype,unitP,unitE, timel,tempstop, sigtot,potist,maxtcel,tempstopcel,lpkbar,angst,leev,iteration,&
         &itetemp,fsumstop,fpstop,itetimestep,sigstop,temp,timemax,cunitE,cunitP,erg2eV, lspaceNDM,latcomp,rang,itesigma,&
         &itetemp2,ihbox0

    USE var_pot, ONLY:
    implicit none

    class(atom_config_d)::atdml
    type(cell_config):: celndm
    class(box_config)::boxndm
    type(para_space_config)::psc
    logical,optional::lreturn


    integer :: nacou, i, ic, iti,it1,it2,it3
    real(double) :: vv, a1, a2, a3, c1, c2, c3
    real(double), dimension(1,3) :: g1,aux
    real(double) :: ltc, ctime, tdev, tcool, epc1, epc2, epc3,masstot,massa,tclt
    real(double), dimension(1,3) :: xtr, cv
    real(double) :: fpmax,fpn,forctot,formax,fpsmax,sigtoth0(3,3)
    real(double) :: potistmean,potistdif
    real(double),save :: potist1000
    real, allocatable,save :: potiststock(:)
    save ltc
    !-----------------------------------------------
    !
    !
    sigtoth0(1:3,1:3)=sigtot(1:3,1:3)*ihbox0(1:3,1:3)    
    lreturn=.false.
    if (timel>=timemax) then
       if (rang==0) write (6, *) '*******max time reached **** ',timel,timemax
       if (present(lreturn)) then
          lreturn=.true.
          return
          else
             call endrunT(atdml,celndm,boxndm,latcomp)
             !       call DeallocateAll
             call arret_ndm
          end if

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
       if (mod(iteration,itetemp)==0) then
          if (temp<=tempstop) then
             if (rang==0)  write (6, *) 'temperature < tempstop '
             if (present(lreturn)) then
                lreturn=.true.
                return
             else
                call endrunT(atdml,celndm,boxndm,latcomp)
                !       call DeallocateAll
                call arret_ndm
             end if
          endif
          if (tempstopcel.gt.0) then
             if (maxtcel<=tempstopcel) then
                write (6, *) 'temperature dans toutes les cels < tempstopcel '
                if (present(lreturn)) then
                   lreturn=.true.
                   return
                else
                   call endrunT(atdml,celndm,boxndm,latcomp)
                   !       call DeallocateAll
                   call arret_ndm
                end if
             end if
          endif
       endif
    endif


    ! change in time step ? itetimestep >0
    if (itetimestep>0) then
       if (mod(iteration,itetimestep)==0) call deftimestep(atdml,boxndm)
    endif


    select case (dmtype)

    case(21,22,23,24)
       if ((fpstop>0.0).AND.(iteration.GE.1)) then

          fpmax=sqrt( MAXVAL( Sum(atdml%fp(1:3,1:atdml%im)**2,1) ) )
          fpSmax = MaxVal( Abs(atdml%fp(:,1:atdml%im)) )
#ifdef PARA
          if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
             call comm_space%max(fpmax)
             call comm_space%max(fpsmax)
          end if
#endif

          fpn=fpSmax*erg2eV/angst
          if (myidsp==0)      write(6,'("TR: force max, energy",i6,3E20.10)') iteration,fpn, potist*erg2eV
!          if ( myidsp==0)     write (6, *) 'energie ',potist*erg2eV
          if((myidsp==0).and.(sigstop.ge.0))write(6,*)'sigma max kbar', 1d-9*maxval(abs(sigtot)), 1d-9*maxval(abs(sigtoth0))
!!$                 write (unitgc, *)
!!$       write (unitgc, *) '************ STRESS in ', cunitP
!!$       do ic = 1, 3
!!$          write (unitgc, '(I1,3(A,I1),A,3G18.10)') ic,' sigma potentiel (1,', ic, ') (2,', ic, &
!!$               ') (3,', ic, ') =',sig(1:3,ic)*unitP
!!$          ppot = ppot+1.0/3.0*sig(ic,ic)
!!$       end do
!!$       write(unitgc,'(A,G18.10)')'PRESSURE',ppot*unitP

          if (fpn.le.fpstop)then
             select case(dmtype)
             case(22,24)
                if(maxval(abs(sigtoth0)).le.sigstop/1d-9) then
                   it1=itetemp;it2=itesigma;it3=itetemp2
                   itetemp=1;itesigma=1;itetemp2=1
                   call analyseT(atdml,celndm,boxndm,psc)
                   
                   itetemp=it1;itesigma=it2;itetemp2=it3
                   if (present(lreturn)) then
                      lreturn=.true.

                      return
                   else
                      call endrunT(atdml,celndm,boxndm,latcomp)
                      !       call DeallocateAll
                      call arret_ndm
                   end if
                endif
             case(21,23)
                itetemp=1;itesigma=1
                   if (present(lreturn)) then
                      lreturn=.true.

                      return
                   else
                      call endrunT(atdml,celndm,boxndm,latcomp)
                      !       call DeallocateAll
                      call arret_ndm
                   end if

                end select

             end if
          end if

       if ((fsumstop>0.0).AND.(iteration.GE.1)) then
          fpSmax=sqrt( SUM(atdml%fp(:,1:atdml%im)**2) )
#ifdef PARA
          if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
             fpSmax=fpSmax**2
             call comm_space%sum(fpsmax)
             fpsmax=sqrt(fpsmax)
          end if
#endif


          fpn=fpsmax*erg2eV/angst
          !          if (myidsp==0)      write(6,*)
          if (myidsp==0)      write(6,*)'TR:  sqrt ( sum_f F_i^2 ):  cgs  ev/Ang ', iteration,fpn, potist*erg2eV
          if((myidsp==0).and.(sigstop.ge.0))write(6,*)'sigma max kbar',1d-9* maxval(abs(sigtot)), 1d-9*maxval(abs(sigtoth0))
          if (fpn.le.fsumstop) then
             select case(dmtype)
             case(22,24)
                if(maxval(abs(sigtoth0)).le.sigstop/1d-9) then
                   itetemp=1;itesigma=1;itetemp2=1
                   call analyseT(atdml,celndm,boxndm,psc)
                   if (present(lreturn)) then
                      lreturn=.true.
                      return
                   else
                      call endrunT(atdml,celndm,boxndm,latcomp)
                      !       call DeallocateAll
                      call arret_ndm
                   end if

                end if
             case(21,23)
                if (present(lreturn)) then
                   lreturn=.true.
                   return
                else
                   call endrunT(atdml,celndm,boxndm,latcomp)
                   call arret_ndm
                end if
             end select
          end if
       end if
       !     if (rang==0) write (6, '(I10,A,D21.12,A)') it,  '*Epot = ', potist*unitE, cunitE
123    continue
       
    case(3,30) ! Gradient conjugue sur coordonnee cartesiennes (3) ou reduites (30)


       if (iteration==1) then
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

       IF (iteration.GE.1) THEN
          forctot=sqrt( SUM(atdml%fp(1:3,1:atdml%im)**2) )
          formax = MaxVal( Abs(atdml%fp(:,1:atdml%im)) )
#ifdef PARA
          if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
             call comm_space%sum(formax)
             forctot=forctot**2
             call comm_space%sum(forctot)
             forctot=sqrt(forctot)
          end if
#endif

          if (lEev.EQV..true.) then
             forctot = forctot*erg2eV/angst
             formax  = formax*erg2eV/angst
             if (myidsp==0) write(*,'("GC: ",i6,3E20.10)') iteration,forctot, formax, potist*erg2eV
             if (fpstop>0) then
                if (formax.le.fpstop) then
                   if (myidsp==0) write(6,*)'force par atome  max  ev/Ang ', formax
                   if (myidsp==0) write (6, *) 'energie ', potist*erg2eV
                   if (present(lreturn)) then
                      lreturn=.true.
                      return
                   else
                      call endrunT(atdml,celndm,boxndm,latcomp)
                      call arret_ndm
                   end if
                end if
             end if
             if (fsumstop>0) then
                if (forctot.le.fsumstop) then
                   if (myidsp==0) write(6,*)'  sqrt ( sum_f F_i^2 ):   ev/Ang ', forctot
                   if (myidsp==0) write(6, *) 'energie ', potist*erg2eV
                   if (present(lreturn)) then
                      lreturn=.true.
                      return
                   else
                      call endrunT(atdml,celndm,boxndm,latcomp)
                      !       call DeallocateAll
                      call arret_ndm
                   end if
                end if
             end if

          else
             if (myidsp==0) write(*,'("GC: ",i6,3E20.10)') iteration,forctot, formax, potist
             if (fpstop>0) then
                if (formax.le.fpstop) then
                   if (myidsp==0) write(6,*)'force par atome  max cgs ',formax
                   if (myidsp==0) write (6, *) 'energie ', potist
                   if (present(lreturn)) then
                      lreturn=.true.
                      return
                   else
                      call endrunT(atdml,celndm,boxndm,latcomp)
                      !       call DeallocateAll
                      call arret_ndm
                   end if

                end if
             end if

             if (fsumstop>0) then
                if (forctot.le.fsumstop) then
                   if (myidsp==0) write(6,*)'  sqrt ( sum_f F_i^2 ):  cgs  ', forctot
                   if (myidsp==0) write (6, *) 'energie ', potist
                   if (present(lreturn)) then
                      lreturn=.true.
                      return
                   else
                      call endrunT(atdml,celndm,boxndm,latcomp)
                      !       call DeallocateAll
                      call arret_ndm
                   end if
                end if
             end if

          end if
       end if ! it .ge.1

    case default
    end select
    !
    return
  end subroutine controleT
   end module controleT_mod
