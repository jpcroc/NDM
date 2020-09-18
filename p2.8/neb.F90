module neb_mod
  USE calfo_mod,only: calfo
  USE analyse_mod,only: analyse
  USE trempe_mod,only: trempe
  USE neb_controle_mod,only:neb_controle
!  USE scalebox_mod,only: scalebox
  USE sauveforce_mod,only: sauveforce
  USE gen_com_m, ONLY:iteanaposneb,itesauvforce,itesauvposition,lfire,maxneb,neb_noise,nebrelaxation,cunitp,&
       &erg2ev,indi,itesauv,lpkbar,ltabvois,nebtype,potist,sig,unitp,potist,sigtot,angst,nvois,itetabvois,rang,&
       &celsize
  
  USE tab_imm_m,only: xp,xpp,vp,ityp,iwmax,fp,ielat,num_at_glob
  USE atomconfig,only:atom_config,atom_config_d,ndm2config,config2ndm
  USE cellconfig, only:cell_config,ndm2cellconfig,cellconfig2ndm,caltabtC
  use var_pot,only:coord
  use rasmol_mod,only:rasmol
  use calfoberend_mod,only:dynlangevin
  use period_mod,only:period
  use neb_module
  USE period_mod,only: period
  USE caltabi_mod,only: caltabi

#ifdef PARANEB
  USE paraneb_mod
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
    !integer  :: iph,i_dir_path
    !integer  :: ielat(imm)
    !integer  :: iwmax(imm)
    !integer  :: ityp(imm)
    !real(double)  :: xp(3,imm)
    !real(double)  :: xpp(3,imm)
    !real(double)  :: vp(3,imm)
    !real(double)  :: fp(3,imm)
    integer :: ineb,ii,it_neb_inter,ipath
    real(double)  :: a_local,forneb

    ! Variables for Fire quench algorithm
    REAL(double), dimension(:), allocatable :: fire_dt, fire_alph
    INTEGER, dimension(:), allocatable :: fire_nstep,iter
!    type(atom_config_d)::atdml
!    type(cell_config)::celndm
#ifdef PARANEB    
    real(double),allocatable:: enepathev_tot(:),enepath_tot(:),sigpath_tot(:,:,:),rc_tot(:)
    integer, allocatable:: nebtest_tot(:)
    real(double)::enertrf,sigpathtrf(3,3),rc_trf
    integer:: iproc,proc_source
    enepath(:)=0
    allocate(enepathev_tot(npath));    allocate(enepath_tot(npath)); allocate(sigpath_tot(3,3,npath))
    allocate(rc_tot(npath)); allocate (nebtest_tot(npath))
    enepathev_tot=0;enepath_tot=0
    if (npath.ne.(nprocs+2)) then
       write(6,*)'nrpocs<>npath-2 ; stop'
       stop
    end if
#endif    
     allocate (iter(npath))

    call ndm2cellconfig(cellneb(1),noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
    do ii=2,npath
       cellneb(ii)=cellneb(1)
    end do
!    call cellneb(npath)%print
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
    
    
    call init_neb 

    ! Initialization of fire quench algorithm
    IF (lFire) THEN
       ALLOCATE(fire_dt(1:npath))
       ALLOCATE(fire_nstep(1:npath))
       ALLOCATE(fire_alph(1:npath))
       do ii=1,npath
          CALL init_trempe_fire(fire_dt(ii), fire_nstep(ii), fire_alph(ii))
       END DO
    END IF

    if (nebrelaxation==1) then
       if (rang==0) write(6,*)'NEB: nebrelaxation == 1'
       if (rang==0) write(6,*)'NEB: relaxation NEB que pour les atomes' 
       call find_relax()
    else
       do ipath=1,npath
          atneb(ipath)%lgul(:)=.true.
       end do
       !irelax(:)=1
       if (rang==0) write(6,*)'NEB: nebrelaxation /= 1'
       if (rang==0) write(6,*)'NEB: relaxation NEB pour TOUS les atomes'
    end if

    if (neb_noise.eq.1) then
       if (rang==0) write(6,*)'NEB: We apply a random noise on the atoms '
       call bruit_neb
    end if
#ifdef PARANEB
    CALL MPI_BARRIER(MPI_COMM_space,ierr)
#endif
    !call calfo
    do ii=1,npath


       it=1
#ifdef PARANEB
       if ((ii==myid+2).or.((ii==1).and.(myid==0)).or.((ii==npath).and.(myid==nprocs-1))) then


#endif       
          !       call scalebox (xp, xpp, vp, ax, fp, ielat, iwmax, ityp,num_at_glob)
          if (lperiod)    call period (imm,atneb(ii)%xp,atneb(ii)%xpp)


          call caltabtC(cellneb(ii),atneb(ii),lperiod,bg)
          if (ltabvois.and.(dmtype==9).and.((it==1).or.(mod(it,itetabvois)==0)))&
               &call caltabi(atneb(ii)%atom_config,cellneb(ii))
          CALL CalFo(sig,potist,atneb(ii),cellneb(ii)) 


          call neb_controle(ii,atneb(ii)%xp,atneb(ii)%fp,atneb(ii)%im) 
          enePATH(ii)=potist
          enePATHev(ii)=potist*erg2ev
          sigPATH(:,:,ii) = sigtot(:,:)      ! Contrainte

#ifdef PARANEB
       endif
#endif       


    end do
#ifdef PARANEB
    write(6,*)
    CALL MPI_BARRIER(MPI_COMM_space,ierr)
    call MPI_ALLREDUCE(enepathev,enepathev_tot,npath,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
    call MPI_ALLREDUCE(enepath,enepath_tot,npath,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
    enepathev(:)=enepathev_tot ; enepath=enepath_tot
#endif
    if (rang==0) then
       do ii=1,npath
          write(*,'(i5,3(g20.8,1x))') ii, enePATHev(ii),enePATHev(ii)-enePATHev(1)
       end do
    end if
    !    CALL MPI_BARRIER(MPI_COMM_space,ierr)


    !stop
    select case (nebtype)


    case (1)
       if (rang==0) write(6,*) 'NEB: !!!!-------this is DRAG----------!!!!!!'
       call build_s_path_drag(ityp)

       do ii=2,npath-1
#ifdef PARANEB
          if (ii==myid+2) then
             enepath(2:npath-1)=0; enepathev(2:npath-1)=0
             if (myid.ne.0) then
                enepath(1)=0;enepath(npath)=0;enepathev(1)=0;enepathev(npath)=0
             end if
#endif
!             call into_path(ii,2,xp, xpp, vp,  fp, ielat, iwmax, ityp)

             it=0
             dragtest=0
             do while (dragtest==0)
                it=it+1
                if (lperiod)    call period (imm,atneb(ii)%xp,atneb(ii)%xpp)
!                call caltabt(im,xp,ielat)
!          call ndm2cellconfig(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
!                call ndm2config(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,&
!                     &iwmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
                call caltabtC(cellneb(ii),atneb(ii),lperiod,bg)
                if (ltabvois.and.(dmtype==9).and.((it==1).or.(mod(it,itetabvois)==0)))&
                     &call caltabi(atneb(ii)%atom_config,cellneb(ii))
                CALL CalFo(sig,potist,atneb(ii),cellneb(ii)) 
!                call config2ndm(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,xpp=xpp)
!    call cellconfig2ndm(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize) !inutile (calfo ne change pas celndm) mais laissé par sécurite
                call force_projection(ii,atneb(ii)%xp,  atneb(ii)%vp,  atneb(ii)%fp,  atneb(ii)%ityp)
                IF (lFire) THEN
                   call trempe_fire(atneb(ii)%xp, atneb(ii)%xpp, atneb(ii)%vp,  atneb(ii)%fp, atneb(ii)%ielat, &
                        &atneb(ii)%iwmax, atneb(ii)%ityp, &
                        fire_dt(ii), fire_nstep(ii), fire_alph(ii))
                ELSE
                   call trempe(atneb(ii)%xp, atneb(ii)%xpp, atneb(ii)%vp, atneb(ii)%fp, atneb(ii)%ielat, &
                        &atneb(ii)%iwmax, atneb(ii)%ityp)
                ENDIF
                !             call analyse  
                call neb_controle(ii,atneb(ii)%xp,atneb(ii)%fp,atneb(ii)%im)
                !call neb_controle(ii) 
             end do   ! end do for a while

             enePATH(ii)=potist
             enePATHev(ii)=potist*erg2ev
             sigPATH(:,:,ii) = sigtot(:,:)      ! Contrainte
!             call into_path(ii,1,xp, xpp, vp,  fp, ielat, iwmax, ityp)
#ifdef PARANEB

          end if
#endif
          iter(ii)=it

       end do      !end ii,npath
#ifdef PARANEB
       write(6,*)'rg ene',myid,enepathev
       CALL MPI_BARRIER(MPI_COMM_space,ierr)
       call MPI_ALLREDUCE(enepathev,enepathev_tot,npath,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
       call MPI_ALLREDUCE(enepath,enepath_tot,npath,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
       call MPI_ALLREDUCE(sigpath,sigpath_tot,9*npath,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
       call MPI_ALLREDUCE(iter,iter,npath,MPI_INTEGER,MPI_SUM,MPI_COMM_space,ierr)
       enepathev(:)=enepathev_tot ; enepath=enepath_tot;sigpath=sigpath_tot
#endif
       if (rang==0) then
          do ii= 2,npath-1
             write(6,*)
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
          it=0
          !   
          call build_s_path_neb(ityp)       
          ! 
          if (rang==0) print'("NEB:===============================================")'
          if (rang==0) print'("NEB:pas-neb image    force      force_NEB          energie      statut    energie/stable")'
          do ii=2,npath-1
#ifdef PARANEB
          if (ii==myid+2) then
             enepath(2:npath-1)=0; enepathev(2:npath-1)=0
!             if (myid.ne.0) then
!                enepath(1)=0;enepath(npath)=0;enepathev(1)=0;enepathev(npath)=0
             !             end if
#endif
             !
!             call into_path      (ii,2,xp, xpp, vp,  fp, ielat, iwmax, ityp)
             !
             it_neb_inter=0  
             !
             do while (it_neb_inter<=5)   ! drag-ize me that 5 steps while we keep NEB "attraction"
                !
                it_neb_inter=it_neb_inter+1
                it=it_neb_inter

                if (lperiod)    call period (imm,atneb(ii)%xp,atneb(ii)%xpp)
!                call caltabt(im,xp,ielat)
!                call ndm2cellconfig(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
!                call ndm2config(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,&
!                     &iwmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
                call caltabtC(cellneb(ii),atneb(ii),lperiod,bg)
                if (ltabvois.and.(dmtype==9).and.((it==1).or.(mod(it,itetabvois)==0)))&
                     &call caltabi(atneb(ii)%atom_config,cellneb(ii))
                CALL CalFo(sig,potist,atneb(ii),cellneb(ii)) 
!                call config2ndm(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,xpp=xpp)
!    call cellconfig2ndm(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize) !inutile (calfo ne change pas celndm) mais laissé par sécurite
!                write(6,*)'********PRE FPN ****************it',it,ii
!                call atneb(ii)%print(i1=2049,i2=2049)
                call force_projection_neb(ii,atneb(ii)%xp,  atneb(ii)%vp,  atneb(ii)%fp, atneb(ii)%ityp)
 !               write(6,*)'*********POST FPN***************it',it,ii
 !               call atneb(ii)%print(i1=2049,i2=2049)
                IF (lFire) THEN
! CRC nettoyer ces appels !                   !
                   call trempe_fire(atneb(ii)%xp, atneb(ii)%xpp,atneb(ii)%vp,  atneb(ii)%fp, atneb(ii)%ielat, &
                        &atneb(ii)%iwmax, atneb(ii)%ityp, &
                        fire_dt(ii), fire_nstep(ii), fire_alph(ii))
                ELSE
                   call trempe(atneb(ii)%xp, atneb(ii)%xpp, atneb(ii)%vp, atneb(ii)%fp, atneb(ii)%ielat,&
                        &atneb(ii)%iwmax, atneb(ii)%ityp)
                ENDIF
 !               write(6,*)'**********POST TRP ***********it',it,ii
 !               call atneb(ii)%print(i1=2049,i2=2049)
! CRC CA CRAINT !
!                call analyse
#ifdef PARANEB
                nebtest(:)=0
#endif
                call neb_controle(ii,atneb(ii)%xp,atneb(ii)%fp,atneb(ii)%im)
 !               write(6,*)'********POST NC*************it',it,ii
 !               call atneb(ii)%print(i1=2049,i2=2049)

                !call neb_controle(ii) 
                !
             end do
             enePATH(ii)=potist
             enePATHev(ii)=potist*erg2eV       
             sigPATH(:,:,ii) = sigtot(:,:)      ! Contrainte
!             call into_path       (ii,1, xp, xpp, vp,  fp, ielat, iwmax, ityp)

#ifdef PARANEB

             if (myid.lt.nprocs-1) call MPI_SEND(enepath(ii), 1,   NDM_MPI_REAL_DOUBLE,   myid+1,10001,MPI_COMM_space,ierr)
               if (myid.gt.0)  call MPI_RECV(enepath(ii-1),1,   NDM_MPI_REAL_DOUBLE,   myid-1,10001,MPI_COMM_space,status,ierr)

               if (myid.lt.nprocs-1)    call MPI_SEND(atneb(ii)%xp(1:3,1:im), 3*im,   NDM_MPI_REAL_DOUBLE,  myid+1,10002,MPI_COMM_space,ierr)
                if (myid.gt.0)                 call MPI_RECV(atneb(ii-1)%xp(1:3,1:im),3*im,   NDM_MPI_REAL_DOUBLE,   myid-1,10002,MPI_COMM_space,status,ierr)
              
                if (myid.gt.0) call MPI_SEND(enepath(ii), 1,   NDM_MPI_REAL_DOUBLE,   myid-1,10003,MPI_COMM_space,ierr)
              if (myid.lt.nprocs-1)    call MPI_RECV(enepath(ii+1),1,   NDM_MPI_REAL_DOUBLE,   myid+1,10003,MPI_COMM_space,status,ierr)

              if (myid.gt.0)  call MPI_SEND(atneb(ii)%xp(1:3,1:im), 3*im,   NDM_MPI_REAL_DOUBLE,  myid-1,10004,MPI_COMM_space,ierr)
               if (myid.lt.nprocs-1)   call MPI_RECV(atneb(ii+1)%xp(1:3,1:im),3*im,   NDM_MPI_REAL_DOUBLE,   myid+1,10004,MPI_COMM_space,status,ierr)
               
             enepathev(:)=enepath(:)*erg2ev
!             stop
#endif             
             ! Backup image ii
 !            write(6,*)'********PRE BACKUP*************it',it,ii
 !       call atneb(ii)%print(i1=2049,i2=2049)
        
             if (itesauv.GT.0) then
                if (mod(ineb,itesauv)==0) then
                   call into_path       (ii,2, xp, xpp, vp,  fp, ielat, iwmax, ityp,num_at_glob)
                   call sauveposition(ii)
                end if
             else if (itesauvposition.GT.0) then
                if (mod(ineb,itesauvposition)==0) then
                   call into_path       (ii,2, xp, xpp, vp,  fp, ielat, iwmax, ityp,num_at_glob)
                   call sauveposition ( ii)
                end if
             endif
             if (itesauvforce.GT.0) then
                if (mod(ineb,itesauvforce)==0) then
                   call into_path       (ii,2, xp, xpp, vp,  fp, ielat, iwmax, ityp,num_at_glob)
                   call sauveforce ( ii)
                end if
             endif

             !debug          print'("NEB: ",2i5,3E14.5,E20.10,i3)', ineb, ii,  formax,       &
             !debug	            formaxperp, formaxparl, potist*erg2eV,nebtest(ii)
             forneb = SQRT(MAXVAL(force_neb(1,:,ii)**2 + force_neb(2,:,ii)**2         &
                  + force_neb(3,:,ii)**2))*erg2eV/angst 
              print'("NEB: ",2i5, 2g14.5,g20.10,i3,g15.8)', ineb, ii,  formax, forneb,       &
                  potist*erg2eV,nebtest(ii),potist*erg2eV-enepathev(1)
             ! 

             !
#ifdef PARANEB
           end if
#endif
 !       write(6,*)'********PREENDO loop path*************it',it,ii
 !       call atneb(ii)%print(i1=2049,i2=2049)

        end do      ! end ii,path
 !       write(6,*)'********POST loop path*************it',it,ii
 !       call atneb(ii)%print(i1=2049,i2=2049)


#ifdef PARANEB
       !           write(6,'(A,10I4)')'rg AV nebtest',myid, nebtest(2:npath-1)
       CALL MPI_BARRIER(MPI_COMM_space,ierr)
       call MPI_ALLREDUCE(nebtest,nebtest_tot,npath,MPI_INTEGER,MPI_SUM,MPI_COMM_space,ierr)
       nebtest=nebtest_tot
!       write(6,'(A,10I4)')'rg AP nebtest',myid,nebtest(2:npath-1)

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
       stop

    end select

    ! Initialization of fire quench algorithm
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
#ifdef PARANEB
       if ((ii==1).or.(ii==npath)) cycle
       if (ii==myid+2) then
#endif
          call into_path      (ii,2,xp, xpp, vp,  fp, ielat, iwmax, ityp,num_at_glob)
          call sauveposition(ii)      
          call rasmol(ii)
          if (iteanaposneb.gt.0) call anapos(ii)
          !	 
          reaction_coord(ii) = SUM((atneb(ii)%xp(:,:)-atneb(1)%xp(:,:))*(atneb(npath)%xp(:,:)-atneb(1)%xp(:,:)))/a_local
!          reaction_coord(ii) = SUM((xp_n(:,:,ii)-xp_n(:,:,1))*(xp_n(:,:,npath)-xp_n(:,:,1)))/a_local 
          !
#ifdef PARANEB
       endif
#endif          

    end do

#ifdef PARANEB
    write(6,*)
    CALL MPI_BARRIER(MPI_COMM_space,ierr)

    if (myid==0) then
       do iproc=1,nprocs-1
          ! Pour le processeur maitre il n'y a rien a faire
          ! reception des donnees des autres processeurs
          !          if (iproc.ne.0) then
          call MPI_RECV(enertrf,               1,NDM_MPI_REAL_DOUBLE,      MPI_ANY_SOURCE, 10001, MPI_COMM_space, status, ierr)
          proc_source = status(MPI_SOURCE)
          enepath(proc_source+2)=enertrf
          call MPI_RECV(sigpathtrf,9,NDM_MPI_REAL_DOUBLE,      MPI_ANY_SOURCE, 10002, MPI_COMM_space, status, ierr)
          proc_source = status(MPI_SOURCE)
          sigpath(:,:,proc_source+2)=sigpathtrf(:,:)
          call MPI_RECV(rc_trf,1,NDM_MPI_REAL_DOUBLE,      MPI_ANY_SOURCE, 10005, MPI_COMM_space, status, ierr)
          proc_source = status(MPI_SOURCE)
          reaction_coord(proc_source+2)=rc_trf
          
          !          endif
       end do
       call MPI_RECV(enertrf,               1,NDM_MPI_REAL_DOUBLE,   nprocs-1, 10003, MPI_COMM_space, status, ierr)
       enepath(npath)=enertrf
       call MPI_RECV(sigpathtrf,9,NDM_MPI_REAL_DOUBLE,    nprocs-1, 10004, MPI_COMM_space, status, ierr)
       sigpath(:,:,npath)=sigpathtrf(:,:)


    else ! Les autres processeurs envoient leurs donnees locales
       call MPI_SEND(enepath(myid+2),               1,   NDM_MPI_REAL_DOUBLE,        0,10001,MPI_COMM_space,ierr)
       call MPI_SEND(sigpath(:,:,myid+2),               9,   NDM_MPI_REAL_DOUBLE,        0,10002,MPI_COMM_space,ierr)
       call MPI_SEND(reaction_coord(myid+2),               1,   NDM_MPI_REAL_DOUBLE,        0,10005,MPI_COMM_space,ierr)       
       if (myid==nprocs-1)then
          call MPI_SEND(enepath(npath),               1,   NDM_MPI_REAL_DOUBLE,        0,10003,MPI_COMM_space,ierr)
          call MPI_SEND(sigpath(:,:,npath),               9,   NDM_MPI_REAL_DOUBLE,        0,10004,MPI_COMM_space,ierr)
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



  end subroutine neb
end module neb_mod
