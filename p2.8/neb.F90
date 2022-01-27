module neb_mod
  USE calfo_mod,only: calfo
  USE trempe_mod,only: trempe
  USE neb_controle_mod,only:neb_controle
  USE gen_com_m, ONLY:iteanaposneb,itesauvforce,itesauvposition,lfire,maxneb,neb_noise,nebrelaxation,cunitp,&
       &erg2ev,itesauv,lpkbar,nebtype,sig,unitp,potist,angst,itetabvois,rang,&
       &fnam,lenfnam,lfire,itesauv,itetabvois,iteanaposneb,maxneb,&
       &nebrelaxation,lperiod,lspacendm,latcomp

  use Tpara,only:para_space_config,ierr

  USE atomconfig,only:atom_config,atom_config_d
  USE cellconfig, only:cell_config
  USE boxconfig,only:box_config,periodbox
  use var_pot,only:rumax,ipotentiel
  use rasmolT_mod,only:rasmolT
  use sauvegardeT_mod,only:sauvegardet
  use neb_module,only:cellneb,atneb,sigpath,boxneb,npath,enepath,nebtype,enepathev,reaction_coord,&
       &lvzeroneb,dragtest,nebtest,force_neb,formax,init_neb,find_relax,bruit_neb,build_s_path_drag,&
       &force_projection,build_s_path_neb,force_projection_neb,paraneb,pscneb
  USE parautils,only:initloc,pointer_caltabt_calfo

#ifdef PARA
  use mpi ! utile pour barrière générale
  use Tpara,only: comm_space
  USE init_vois_mod,only: init_voisinage
#endif

  implicit none

contains
  subroutine neb 
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    !-----------------------------------------------
    USE posana
    USE FireModule
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------

    logical::lmaster
    logical::lcalcvois
    integer :: ineb,ii,it_neb_inter,ipath
    real(double)  :: a_local,forneb
    character::fnamcout*80,extension*9
    integer::formatsauv,i1
    ! Variables for Fire quench algorithm
    REAL(double), dimension(:), allocatable :: fire_dt, fire_alph
    INTEGER, dimension(:), allocatable :: fire_nstep,iter
    class(atom_config),pointer::atnebloc
    type(cell_config),pointer::cellnebloc
    type(cell_config),target:: cellcible ! ne sert qu'à faire pointer cellnebloc sur quelquechose
    type(atom_config_d),target::atcible
    integer::iun,i,ic
    logical::lchange

#ifdef PARA    
    real(double)::enertrf,sigpathtrf(3,3),rc_trf
    integer:: iproc,proc_source

    enepath(:)=0
    cellnebloc=>cellcible
    atnebloc=>atcible
#endif    

    latcomp=.true.! NEB=> latcomp=.true.
    allocate (iter(npath))
    do ii=1,npath
       if ((ii.ne.1).and.(ii.ne.npath))then
          cellneb(ii)=cellneb(1)
       end if
    end do
    !APRES    
    call init_neb(atneb(1)%im,atneb(1)%imm)
#ifdef PARA
    lmaster=paraneb%lmaster
    if ((paraneb%mpi_image%nproc.gt.1).and.(lspaceNDM.eqv..true.)) then
       call init_voisinage(cellneb(1),pscneb)
    end if


#else
    lmaster=.true.
#endif    
    if (neb_noise.eq.1) then
       if (rang==0) write(6,*)'NEB: We apply a random noise on the atoms '
       call bruit_neb(atneb(1)%im)

    end if
    !ENCORE AVANT

    !on connait atneb(ii),cellneb(ii) et boxneb
    if (nebrelaxation==1) then
       if (rang==0) write(6,*)'NEB: nebrelaxation == 1'
       if (rang==0) write(6,*)'NEB: relaxation NEB que pour les atomes' 
       call find_relax(atneb(1)%im)
    else
       do ipath=1,npath
          atneb(ipath)%lgul(:)=.true.
       end do
       !irelax(:)=1
       if (rang==0) write(6,*)'NEB: nebrelaxation /= 1'
       if (rang==0) write(6,*)'NEB: relaxation NEB pour TOUS les atomes'
    end if
    !????
    do ii=1,npath
       call periodbox (boxneb,atneb(ii)%atom_config_d)

    end do
    !AVANT
    if(lPkbar) then
       unitP=1.0d-9  ;     cunitP='kbar'
    else
       unitP=1.0
       cunitP='d/cm2'
    endif
    if (rang==0) then
       open(unit=55,file='image_col_relax.out')
       WRITE(55,'(a)')  '#  1: (i-1)/(nPath-1)'
       WRITE(55,'(a)')  '#  2: energy E(i) (eV)'
       WRITE(55,'(a)')  '#  3: energy difference E(i)-E(1) (eV)'
       WRITE(55,'(3a)') '#  4: stess s(1,1) (', cUnitP, ')'
       WRITE(55,'(a)')  '#  5:       s(2,2)'
       WRITE(55,'(a)')  '#  6:       s(3,3)'
       WRITE(55,'(a)')  '#  7:       s(2,3)'
       WRITE(55,'(a)')  '#  8:       s(1,3)'
       WRITE(55,'(a)')  '#  9:       s(1,2)'
       open(unit=56,file='react_col_relax.out')
       WRITE(56,'(a)')  '#  1: reaction coordinate z(i)'
       WRITE(56,'(a)')  '#  2: energy E(i) (eV)'
       WRITE(56,'(a)')  '#  3: energy difference E(i)-E(1) (eV)'
       WRITE(56,'(3a)') '#  4: stess s(1,1) (', cUnitP, ')'
       WRITE(56,'(a)')  '#  5:       s(2,2)'
       WRITE(56,'(a)')  '#  6:       s(3,3)'
       WRITE(56,'(a)')  '#  7:       s(2,3)'
       WRITE(56,'(a)')  '#  8:       s(1,3)'
       WRITE(56,'(a)')  '#  9:       s(1,2)'
    end if



    ! Initialization of fire quench algorithm
    IF (lFire) THEN
       ALLOCATE(fire_dt(1:npath))
       ALLOCATE(fire_nstep(1:npath))
       ALLOCATE(fire_alph(1:npath))
       do ii=1,npath
          CALL init_trempe_fire(fire_dt(ii), fire_nstep(ii), fire_alph(ii))
       END DO
    END IF

    iteration=1
#ifdef PARA
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif



    do i1=1,npath
#ifdef PARA

       !       if (i1==paraneb%image+2) then
       if ((i1==paraneb%image+2).or.((i1==1).and.(paraneb%image==0)).or.((i1==npath).and.(paraneb%image==paraneb%nimage-1))) then
          ii=i1
          if (i1==npath)ii=npath-1
          if (i1==npath-1)ii=npath
#else    
          ii=i1
#endif    
          call initloc(atneb(ii)%atom_config_d,cellneb(ii),atnebloc,cellnebloc,boxneb,paraneb,&
               &rumax,lperiod,psc=pscneb,lcalcvois=.true.) !initloc contient caltabtc sur atloc
          if (ipotentiel.lt.0) then
             lchange=.true.
          else
             lchange=.false.
          end if

          call pointer_caltabt_calfo(sig,potist,atneb(ii)%atom_config_d,cellneb(ii),boxneb,atnebloc,cellnebloc,paraneb,&
               &lperiod,lupdate=lchange,psc=pscneb)
          if (lmaster) then

             call neb_controle(ii,atneb(ii)%xp,atneb(ii)%fp,atneb(ii)%im)
!!!!!!!!!!!          broadcast de dragtest nebtest(ii) ou non ?
             enePATH(ii)=potist
             enePATHev(ii)=potist*erg2ev
             sigPATH(:,:,ii) = sig(:,:)      ! Contrainte

          endif
          !          on fait rien si pas master

#ifdef PARA
       endif
#endif       


    end do


#ifdef PARA

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

    if (paraneb%lmaster) then
       call paraneb%mpi_master%sum(enepath)
       call paraneb%mpi_master%sum(enepathev)
    end if
#endif

    select case (nebtype)
    case (1)
       if (rang==0) write(6,*) 'NEB: !!!!-------this is DRAG----------!!!!!!'
       iter=0

       if (lmaster) then
          call build_s_path_drag(atneb(1)%im,atneb(1)%imm)
       endif

       do ii=2,npath-1
#ifdef PARA
          if (ii==paraneb%image+2) then

             enepath(2:npath-1)=0; enepathev(2:npath-1)=0
             !             if (paraneb%image.ne.0) then
             !                enepath(1)=0;enepath(npath)=0;enepathev(1)=0;enepathev(npath)=0
             !             end if
             !             if ((paraneb%image.ne.0).and.(paraneb%image.ne.npath)) then
             !                enepath(1)=0;enepath(npath)=0;enepathev(1)=0;enepathev(npath)=0
             !             end if
#endif
             iteration=0; iter(ii)=0
             dragtest=0
             do while (dragtest==0)
                iteration=iteration+1

                if (lmaster)    call periodbox (boxneb,atneb(ii)%atom_config_d)
                if ((mod(iteration,itetabvois)==0).and.(atneb(ii)%ltabvois))then
                   lcalcvois=.true.
                else
                   lcalcvois=.false.
                end if
                call pointer_caltabt_calfo(sig,potist,atneb(ii)%atom_config_d,cellneb(ii),boxneb,&
                     &atnebloc,cellnebloc,paraneb,lperiod,lupdate=.true.,psc=pscneb,lcalcvois=lcalcvois)

#ifdef PARA

                if (paraneb%lmaster) then
#endif
                   !                   call atneb(ii)%print(unit=100+ii,natg1=1,natg2=2049)
                   !                   flush(100+ii)
                   call force_projection(ii,atneb(ii)%xp,  atneb(ii)%vp,  atneb(ii)%fp,  atneb(ii)%ityp,atneb(ii)%imm,atneb(ii)%im)
                   !                   call atneb(ii)%print(unit=300+ii,natg1=1,natg2=2049)
                   !                   flush(300+ii)
                   !                   call atneb(ii)%print(unit=400+ii)
                   !                   flush(400+ii)

                   IF (lFire) THEN
                      call trempe_fire(atneb(ii),fire_dt(ii), fire_nstep(ii), fire_alph(ii))
                   ELSE
                      call trempe(atneb(ii))
                   ENDIF

                   call neb_controle(ii,atneb(ii)%xp,atneb(ii)%fp,atneb(ii)%im)

#ifdef PARA

                end if
                call paraneb%mpi_image%bcast(0,dragtest)

#endif
             end do   ! end do for a while
             if (lmaster) then 
                enePATH(ii)=potist
                enePATHev(ii)=potist*erg2ev
                sigPATH(:,:,ii) = sig(:,:)      ! Contrainte
             end if
             iter(ii)=iteration
#ifdef PARA

          end if
#endif


       end do      !end ii,npath



#ifdef PARA
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       if(lmaster) then
          call paraneb%mpi_master%sum(enepath)
          call paraneb%mpi_master%sum(enepathev)
          call paraneb%mpi_master%sum(sigpath)
          call paraneb%mpi_master%sum(iter)
       end if
#endif
       if (rang==0) then
          do ii= 2,npath-1
             !             write(*,*) 'NEB: THIS IS THE DRAG IMAGE=====================', ii
             write(*,*) 'NEB: THIS IS THE DRAG IMAGE=====================', ii,iter(ii)
             write(*,*) 'NEB: THE ENERGY OF THIS IMAGE===================', enePATHev(ii)
             WRITE(6,'(3a)') 'NEB: stress tensor in Voigt notation (units: ', cunitP,' ):'
             WRITE(6,'(6g20.8)') unitP*sigPATH(1,1,ii), unitP*sigPATH(2,2,ii), unitP*sigPATH(3,3,ii), &
                  0.5*unitP*(sigPath(2,3,ii)+sigPath(3,2,ii)), &
                  0.5*unitP*(sigPath(1,3,ii)+sigPath(3,1,ii)), &
                  0.5*unitP*(sigPath(1,2,ii)+sigPath(2,1,ii))
          end do
       end if

    case(2)
       if (rang==0) write(6,*)'NEB: -------this is NEB--V2-------'
       if (rang==0) write(6,*)'NEB: The MAX steps in NEB        :',maxneb
       nebtest(:)=0
       do ineb=1,maxneb     ! Main loop for NEB
          if (lvzeroneb)then
             do ipath=2,npath-1
                atneb(ipath)%vp(:,:)=0
             end do
          end if
          iteration=0
          !   
          call build_s_path_neb(atneb(1)%im,atneb(1)%imm)
          ! 
          if (rang==0) write(6,*)'("NEB:===============================================")'
          if (rang==0) write(6,*)'("NEB:pas-neb image    force      force_NEB          energie      statut    energie/stable")'
          do ii=2,npath-1
#ifdef PARA
             if (ii==paraneb%image+2) then
                enepath(2:npath-1)=0; enepathev(2:npath-1)=0
#endif
                it_neb_inter=0  

                do while (it_neb_inter<=5)   ! drag-ize me that 5 steps while we keep NEB "attraction"
                   !
                   it_neb_inter=it_neb_inter+1
                   iteration=it_neb_inter

                   call periodbox (boxneb,atneb(ii)%atom_config_d)
                   if ((mod(iteration,itetabvois)==0).and.(atneb(ii)%ltabvois))then
                      lcalcvois=.true.
                   else
                      lcalcvois=.false.
                   end if
                   call pointer_caltabt_calfo(sig,potist,atneb(ii)%atom_config_d,cellneb(ii),&
                        &boxneb,atnebloc,cellnebloc,paraneb,lperiod,&
                        &lupdate=.true.,psc=pscneb,lcalcvois=lcalcvois)
                   if (lmaster) then
                      call force_projection_neb(ii,atneb(ii)%xp,  atneb(ii)%vp,  atneb(ii)%fp, atneb(ii)%ityp,&
                           &atneb(ii)%imm,atneb(ii)%im)

                      IF (lFire) THEN
                         ! CRC nettoyer ces appels !                   !
                         call trempe_fire(atneb(ii),fire_dt(ii), fire_nstep(ii), fire_alph(ii))
                      ELSE
                         call trempe(atneb(ii))
                      ENDIF
#ifdef PARA
                      nebtest(:)=0
#endif
                      call neb_controle(ii,atneb(ii)%xp,atneb(ii)%fp,atneb(ii)%im)
                   end if
                end do
                if (lmaster) then 
                   enePATH(ii)=potist
                   enePATHev(ii)=potist*erg2eV       
                   sigPATH(:,:,ii) = sig(:,:)      ! Contrainte

#ifdef PARA
                   if (paraneb%mpi_master%rank.lt.paraneb%nimage-1) then
                      call paraneb%mpi_master%send(enepath(ii),paraneb%mpi_master%rank+1,10001)
                   end if
                   if (paraneb%mpi_master%rank.gt.0) &
                        & call paraneb%mpi_master%recv(enepath(ii-1),paraneb%mpi_master%rank-1,10001)
                   if (paraneb%mpi_master%rank.lt.paraneb%nimage-1)   &
                        &call paraneb%mpi_master%send(atneb(ii)%xp(1:3,1:atneb(ii)%im),paraneb%mpi_master%rank+1,10002)
                   if (paraneb%mpi_master%rank.gt.0)  &
                        &call paraneb%mpi_master%recv(atneb(ii-1)%xp(1:3,1:atneb(ii)%im),paraneb%mpi_master%rank-1,10002)
                   if (paraneb%mpi_master%rank.gt.0) &
                        & call paraneb%mpi_master%send(enepath(ii),paraneb%mpi_master%rank-1,10003)
                   if (paraneb%mpi_master%rank.lt.paraneb%nimage-1)&
                        & call paraneb%mpi_master%recv(enepath(ii+1),paraneb%mpi_master%rank+1,10003)
                   if (paraneb%mpi_master%rank.gt.0)  &
                        &call paraneb%mpi_master%send(atneb(ii)%xp(1:3,1:atneb(ii)%im),paraneb%mpi_master%rank-1,10004)
                   if (paraneb%mpi_master%rank.lt.paraneb%nimage-1) &
                        &call paraneb%mpi_master%recv(atneb(ii+1)%xp(1:3,1:atneb(ii)%im),paraneb%mpi_master%rank+1,10004)
                   enepathev(:)=enepath(:)*erg2ev
                   !             call arret_ndm
#endif             

                   if (itesauv.GT.0) then
                      if (mod(ineb,itesauv)==0) then
                         !                   call into_path       (ii,2, xp, xpp, vp,  fp, ielat, iwmax, ityp,num_at_glob)
                         formatsauv = 2
                         write(extension,'(i9.9)') ii
                         fnamcout = fnam(1:lenfnam)//'.cout.'//extension
                         call sauvegardet(atneb(ii)%atom_config_d, cellneb(ii),boxneb,formatsauv,fnamcout,latcomp=latcomp)
                      end if
                   end if

                   !debug          print'("NEB: ",2i5,3E14.5,E20.10,i3)', ineb, ii,  formax,       &
                   !debug	            formaxperp, formaxparl, potist*erg2eV,nebtest(ii)
                   forneb = SQRT(MAXVAL(atneb(ii)%force_neb(1,:)**2 + atneb(ii)%force_neb(2,:)**2         &
                        + atneb(ii)%force_neb(3,:)**2))*erg2eV/angst 
                   write(6,'(A,2i5, 2g14.5,g20.10,i3,g15.8)')  'NEB: ', ineb, ii,  formax, forneb,       &
                        potist*erg2eV,nebtest(ii),potist*erg2eV-enepathev(1)
                   ! 
                endif
                !
#ifdef PARA
             end if
#endif

          end do      ! end ii,path

          !          write(6,*)'OUTloop',rang
#ifdef PARA
          if (lmaster) then
             call paraneb%mpi_master%sum(nebtest)

          end if

          call paraneb%mpi_image%bcast(0,nebtest)
#endif
          if (SUM(nebtest(2:npath-1))==(npath-2)) then
             if (rang==0) print'("NEB:=== bye,bye sweety is finished=================")'
             if (rang==0) print'("NEB:===============================================")'
             exit
          end if

       end do        !maxneb 

    case default

       if (rang==0) write(6,*) 'NO NEB-DYNAMICS FOR THIS  nebtype = ', nebtype
       if (rang==0) write(6,*) 'CHANGE nebtype AND TRY AGAIN.'
       if (rang==0) write(6,*) 'nebtype=1 for DRAG'
       if (rang==0) write(6,*) 'nebtype=2 for  NEB'
       if (rang==0) write(6,*) 'STOP in the neb.f90'
       call arret_ndm

    end select

    ! de_Initialization of fire quench algorithm
    IF (lFire) THEN
       DEALLOCATE(fire_dt)
       DEALLOCATE(fire_nstep)
       DEALLOCATE(fire_alph)
    END IF

    reaction_coord(1)=0
    reaction_coord(npath)=1
    a_local=SUM((atneb(npath)%xp(:,:)-atneb(1)%xp(:,:))**2)
    !    a_local=SUM((xp_n(:,:,npath)-xp_n(:,:,1))**2) 
    !


    do ii=2,npath-1
#ifdef PARA
       !       if ((ii==1).or.(ii==npath)) cycle
       if ((paraneb%lmaster).and.(ii==paraneb%image+2)) then

#endif

          formatsauv = 2
          write(extension,'(i9.9)') ii
          fnamcout = fnam(1:lenfnam)//'.cout.'//extension
          call sauvegardet(atneb(ii)%atom_config_d, cellneb(ii),boxneb,formatsauv,fnamcout,latcomp)
          call rasmolT(atneb(ii)%atom_config_d,boxneb,ii,latcomp=latcomp)
          if (iteanaposneb.gt.0) call anapos (atneb(ii),cellneb(ii),boxneb,ii)
          !	 
          reaction_coord(ii) = SUM((atneb(ii)%xp(:,:)-atneb(1)%xp(:,:))*(atneb(npath)%xp(:,:)-atneb(1)%xp(:,:)))/a_local
          !          reaction_coord(ii) = SUM((xp_n(:,:,ii)-xp_n(:,:,1))*(xp_n(:,:,npath)-xp_n(:,:,1)))/a_local 
          !
#ifdef PARA
       endif
#endif          

    end do

    if (lmaster) then 
#ifdef PARA


       if (paraneb%mpi_master%rank==0) then
          do iproc=1,paraneb%nimage-1
             ! Pour le processeur maitre il n'y a rien a faire
             ! reception des donnees des autres processeurs
             !          if (iproc.ne.0) then
             call  paraneb%mpi_master%probe(10001,sourceout=proc_source)
             call  paraneb%mpi_master%recv (enertrf,proc_source,10001)
             enepath(proc_source+2)=enertrf

             call  paraneb%mpi_master%probe(10002,sourceout=proc_source)
             call  paraneb%mpi_master%recv (sigpathtrf,proc_source,10002)
             sigpath(:,:,proc_source+2)=sigpathtrf(:,:)

             call  paraneb%mpi_master%probe(10005,sourceout=proc_source)
             call  paraneb%mpi_master%recv (rc_trf,proc_source,10005)
             reaction_coord(proc_source+2)=rc_trf
          end do
          call  paraneb%mpi_master%recv (enertrf,paraneb%nimage-1,10003)
          call  paraneb%mpi_master%recv (sigpathtrf,paraneb%nimage-1,10004)
          enepath(npath)=enertrf
          sigpath(:,:,npath)=sigpathtrf(:,:)
       else ! Les autres processeurs envoient leurs donnees locales
          call paraneb%mpi_master%send(enepath(paraneb%mpi_master%rank+2), 0,10001)
          call paraneb%mpi_master%send(sigpath(:,:,paraneb%mpi_master%rank+2), 0,10002)
          call paraneb%mpi_master%send(reaction_coord(paraneb%mpi_master%rank+2), 0,10005)
          if (paraneb%mpi_master%rank==paraneb%nimage-1)then
             call paraneb%mpi_master%send(enepath(npath), 0,10003)
             call paraneb%mpi_master%send(sigpath(:,:,npath), 0,10004)
          end if
       endif
       enepathev(:)=enepath(:)*erg2ev
#endif
       if (rang==0) WRITE(*,'(3a)') "NEB:--IMAGE-----REACT-COORD------ENERGY------ENERGY-ENERGY(1)&
            &-------------STRESS-sVoigt(1:6)-(units:-", cunitP, ")"
       do ii=1,npath
          !
          if (rang==0) write(*,'(i5,9(g20.8,1x))')  ii, reaction_coord(ii), enePATHev(ii),enePATHev(ii)-enePATHev(1), &
               unitP*sigPATH(1,1,ii), unitP*sigPATH(2,2,ii), unitP*sigPATH(3,3,ii), &
               0.5*unitP*(sigPath(2,3,ii)+sigPath(3,2,ii)), &
               0.5*unitP*(sigPath(1,3,ii)+sigPath(3,1,ii)), &
               0.5*unitP*(sigPath(1,2,ii)+sigPath(2,1,ii))
          if (rang==0)   write(55,'(9(g20.8,1x))') dble(ii-1)/dble(npath-1), enePATHev(ii),enePATHev(ii)-enePATHev(1), &
               unitP*sigPATH(1,1,ii), unitP*sigPATH(2,2,ii), unitP*sigPATH(3,3,ii), &
               0.5*unitP*(sigPath(2,3,ii)+sigPath(3,2,ii)), &
               0.5*unitP*(sigPath(1,3,ii)+sigPath(3,1,ii)), &
               0.5*unitP*(sigPath(1,2,ii)+sigPath(2,1,ii))
          if (rang==0)  write(56,'(9(g20.8,1x))') reaction_coord(ii), enePATHev(ii),enePATHev(ii)-enePATHev(1), &
               unitP*sigPATH(1,1,ii), unitP*sigPATH(2,2,ii), unitP*sigPATH(3,3,ii), &
               0.5*unitP*(sigPath(2,3,ii)+sigPath(3,2,ii)), &
               0.5*unitP*(sigPath(1,3,ii)+sigPath(3,1,ii)), &
               0.5*unitP*(sigPath(1,2,ii)+sigPath(2,1,ii))
          !       
       end do


       if (rang==0) print*,'-------------------------------------------'
       if (rang==0) write(*,'("          Min .. et .. Max")')
       if (rang==0) print '("E:",2g14.7)',minval(enePATHev), maxval(enePATHev)
       if (rang==0) write(*,'("          1 .. et .. NPATH")')
       if (rang==0) print '("B:",2g14.7)',enePATHev(1), enePATHev(npath)
       if (rang==0) write(*,'("les diffs")')
       if (rang==0) print*,'1-NPATH    :',enePATHev(1)-enePATHev(npath)
       if (rang==0) print*,'-------------------------------------------'
       if (rang==0) print*,'MAX-1      :',maxval(enePATHev)-enePATHev(1)
       if (rang==0) print*,'MAX-NPATH  :',maxval(enePATHev)-enePATHev(npath)
    end if

#ifdef PARA
    call mpi_barrier(mpi_comm_world,ierr)
#endif

    return

  end subroutine neb
end module neb_mod
