module work_cgII

  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY: im, imm,at, inv_angst, lperiod, rang,indi,itmax,leev,ltabvois,sig, &
                       it, itesauv, itesauvposition, itesauvforce,itmax, &
                       inv_angst, erg2ev, angst,fpstop,fsumstop,itetabvois, &
                       dmtype, potist,im_glob,nox,noy,noz,cell_finx,cell_finy,cell_finz,noxyz,nvois,&
                       &natperc,nato,ncel,atincel,deltadist,celsize,bg,mdcg_noise
  USE controle_mod,only: controle
  USE calfo_mod,only: calfo
  USE analyse_mod,only: analyse
  USE sauvegarde_mod,only: sauvegarde
  USE sauveposition_mod,only: sauveposition
  USE sauveforce_mod,only: sauveforce
  USE config_mod,only: config
  USE zero2all2zero_mod,only: zero2all,all2zero
  use period_mod,only:period
  USE endrun_mod,only: endrun
 USE dynalloccell,only:deallocateall
USE arret_ndm_mod,only: arret_ndm
USE caltabi_mod,only: caltabi
#ifdef PARA
USE mpi
use mod_para,only: status,ierr,nprocs,myid,NDM_MPI_REAl_DOUBLE,maj_atomes_frt_ftm
#endif
  USE tab_imm_m,only : xp, fp,num_at_glob,ax,vp,xpp,ityp,ielat,iwmax,bruitmd
  USE atomconfig,only : atom_config_d,ndm2config, config2ndm
  USE cellconfig, only:cell_config,ndm2cellconfig,cellconfig2ndm,caltabtC

  
  implicit none


contains

  subroutine FUNCT(N,X,F,G,NCALLS,                      &
       xp_local,  fp_local,   ityp_local,ims)

    real(double),intent(in):: X(N)
    real(double),intent(out)::G(N),F
    integer,intent(in) ::N,NCALLS,ims
    integer ::i
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer,intent(inout)  :: ityp_local(ims)
    real(double),intent(out)  :: xp_local(3,ims)
    real(double),intent(out)  :: fp_local(3,ims)
    !-----------------------------------------------
    !-----------------------------------------------
    integer::iproc,proc_source,cellx,celly,cellz
    integer::nag(imm)
    real(double)::aux,auy,auz
    character :: extension*2
    integer::lenfn2,ko
    real(double) :: fpmax,fpn,forctot,formax,fpmax_glob

    type(atom_config_d)::atcg
    type(cell_config):: celcg
    
!    write(6,*)'entree funct', it,ncalls
    it=NCALLS-1

!    IF (3*ims.NE.N) THEN
!       WRITE(0,'(a,i0)') "3*imm = ", 3*ims
!       WRITE(0,'(a,i0)') "N     = ", N
!       STOP "< work_cgII >"
!    END IF
    xp_local=0
    fp_local=0
    if (rang==0) then
       IF (dmtype.EQ.30) THEN ! Variables = reduced coordinates
          do i=1,im_glob
             xp_local(1:3,i) = MatMul( at, X(3*i-2:3*i) )
          end do
       ELSE ! Variables = cartesian coordinates (in A)
          do i=1,im_glob
             xp_local(1:3,i)=X(3*i-2:3*i)*inv_angst
          end do
       END IF
    end if
!#ifdef PARA
!    write(extension,'(i2.2)') rang
!    lenfn2=2   
!    open(unit=606, file='xp_local.'//extension(1:lenfn2)//'.csv', form='formatted', &
!             status='unknown')
!    do i=1,im_glob
!       !       write(6,*)rang,i,xp_all(:,i)
!       write(606,'(2I2,I6,3G18.9)') rang,ncalls,i,xp_local(:,i)
!    end do
!    call mpi_barrier(MPI_COMM_WORLD,ierr)
!#else
!     open(unit=606, file='xp_local.SEQ.csv', form='formatted', &
!             status='unknown')
!    do i=1,im_glob
!       !       write(6,*)rang,i,xp_all(:,i)
!       write(606,'(2I2,I6,3G18.9)') rang,ncalls,i,xp_local(:,i)
!    end do
!#endif
    
!    write(6,*)'taillesP', size(xp_local),size(xp)
    xp=0
    call zero2all(xp_local,xp,ityp_local)
!      write(6,*)'CALL W',it
#ifndef PARA
    call period (imm,xp)
#endif    

    call  ndm2cellconfig(celcg,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
    call ndm2config(atcg,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,&
         &iwmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
    call caltabtC(celcg,atcg,lperiod,bg)
    if (ltabvois.and.mod(it,itetabvois)==0) then

       call caltabi(atcg%atom_config,celcg)

    end if
       call config2ndm(atcg,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,xpp=xpp)
       call cellconfig2ndm(celcg,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize) !sans doute inutile
    
!    call period
    !    do ko=1,noxyz
    !       write(6,*)'0rg cel nat',rang, ko,nato(ko)
    !    end do
#ifdef PARA
    call maj_atomes_frt_ftm
#endif

    !back to internal units and JP world.......................................


    if (lperiod)          call period  (imm,xp)
!    call controle
     if (it==1) then
        if (lEev.EQV..true.) then 
           if (rang==0) write(6,*)'Resultats en eV, Ang'
        else
           if (rang==0) write(6,*)'Resultats en cgs'
        end if
        if (rang==0)      write(*,'(70("="))')
        if (rang==0)      write(*,'("CG:     ","iter",10(" "),"epsi",14(" "),"Fmax",14(" "), "Energy")')
        if (rang==0)      write(*,'(70("="))')
     end if


     !debug     write(*,*) 'DEBUG ALL IT IN CONTROLE',it

    


!do ko=1,noxyz
!    write(6,*)'1rg cel nat',rang, ko,nato(ko)
    !end do
!     write(6,*)'im_glob', im_glob
!    open(unit=607, file='xpG.csv', form='formatted', &
!             status='unknown')
!    do i=1,im_glob
!      !       write(6,*)rang,i,xp_all(:,i)
!       write(607,'(I6,3G22.13)') i,xp(:,i)
!    end do

      call ndm2cellconfig(celcg,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
  call ndm2config(atcg,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,&
       &iwmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
  CALL CalFo(sig,potist,atcg,celcg)
!  write(6,*)'dml potist ',potist,atdml%potist
    call config2ndm(atcg,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,xpp=xpp)
    call cellconfig2ndm(celcg,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize) !inutile (calfo ne change pas celndm) mais laissé par sécurite
    
!    call calfo
!    open(unit=606, file='fpG.csv', form='formatted', &
!             status='unknown')
!    do i=1,im_glob
!      !       write(6,*)rang,i,xp_all(:,i)
!       write(606,'(I6,3G22.13)') i,fp(:,i)
!    end do

!#ifdef PARA
!    lenfn2=2
!    write(extension,'(i2.2)') rang
!    open(unit=625, file='FORC.'//extension(1:lenfn2)//'.csv', form='formatted', &
!         status='unknown')
!#else    
!    open(unit=625, file='FORC.SEQ.csv', form='formatted', &
!         status='unknown')
!#endif
!    do i=1,im
!       write(625,'(2I2,2I6,3G18.9)') rang, ncalls,i,num_at_glob(i),fp(:,i)
!    end do
#ifdef PARA
    call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif    
!   stop

    !do ko=1,noxyz
!    write(6,*)'2rg cel 1nat',rang, ko,nato(ko)
! end do
     IF (it.GE.1) THEN
!        IF (lFrozen) THEN
           !forctot=sqrt( Sum( SUM(fp(1:3,1:im)**2,1), Free(1:im) ) )
           !formax=sqrt( MAXVAL( Sum(fp(1:3,1:im)**2,1), Free(1:im) ) )
!           forctot = sqrt( SUM( fp(:,1:im)**2, .NOT.Frozcalfoen(:,1:im) ) )
!           formax = MaxVal( Abs(fp(:,1:im)), .NOT.Frozen(:,1:im) ) 
!        ELSE
           forctot=sqrt( SUM(fp(1:3,1:im)**2) )
           !formax=sqrt( MAXVAL( Sum(fp(1:3,1:im)**2,1) ) )
           formax = MaxVal( Abs(fp(:,1:im)) )
!        END IF
!        write(6,*)'FF ', forctot,formax
#ifdef PARA
        call MPI_ALLREDUCE(formax,fpmax_glob,1,NDM_MPI_REAL_DOUBLE,MPI_MAX,MPI_COMM_WORLD,ierr)
        formax=fpmax_glob
        forctot=forctot**2
        call MPI_ALLREDUCE(forctot,fpmax_glob,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        forctot=sqrt(fpmax_glob)
#endif

        if (lEev.EQV..true.) then 
           forctot = forctot*erg2eV/angst
           formax  = formax*erg2eV/angst
           if (rang==0) write(*,'("GC: ",i6,3E20.10)') it,forctot, formax, potist*erg2eV
!           if (rang==0) write(*,'("GC: ",3E20.10)') forctot, formax, potist*erg2eV
           if (fpstop>0) then   
              if (formax.le.fpstop) then
                 if (rang==0) write(6,*)'force par atome  max  ev/Ang ', formax
                 if (rang==0) write (6, *) 'energie ', potist*erg2eV
                 if (it.le.1) xp(:,:)=ax(:,:)
                 call endrun
              end if
           end if
           if (fsumstop>0) then   
              if (forctot.le.fsumstop) then
                 if (rang==0) write(6,*)'  sqrt ( sum_f F_i^2 ):   ev/Ang ', forctot
                 if (rang==0) write(6, *) 'energie ', potist*erg2eV
                 if (it.le.1) xp(:,:)=ax(:,:)
                 call endrun
              end if
           end if

        else
           if (rang==0) write(*,'("GC: ",i6,3E20.10)') it,forctot, formax, potist
           if (fpstop>0) then   
              if (formax.le.fpstop) then
                 if (rang==0) write(6,*)'force par atome  max cgs ',formax
                 if (rang==0) write (6, *) 'energie ', potist
                 if (it.le.1) xp(:,:)=ax(:,:)
                 call endrun

              end if
           end if

           if (fsumstop>0) then   
              if (forctot.le.fsumstop) then
                 if (rang==0) write(6,*)'  sqrt ( sum_f F_i^2 ):  cgs  ', forctot
                 if (rang==0) write (6, *) 'energie ', potist
                 if (it.le.1) xp(:,:)=ax(:,:)
                 call endrun
              end if
           end if

        end if
     end if ! it .ge.1

    
  if (it>=itmax) then
     if (rang==0) write (6, *) '*******Derniere iteration **** '
     call endrun
     call DeallocateAll

     call arret_ndm

  endif






     call all2zero(fp_local,fp)


    call analyse 


    if (it.ne.0) then
!       if (rang==0) then
          !           write(6,*)'work_cg_II analyse -> sauvegarde',it
          if (itesauv.GT.0) then
             if (mod(it,itesauv)==0) call sauvegarde 
          endif

          !            write(6,*)'work_cg_II analyse -> sauveposition',it
          if (itesauvposition.GT.0) then
             if (mod(it,itesauvposition)==0) call sauveposition ( it)
          endif
          if (itesauvforce.GT.0) then
             if (mod(it,itesauvforce)==0) call sauveforce ( it)
          endif
          !            write(6,*)'work_cg_II sauvposition -> control',it
 !      endif                                   ! fin rang=0
    end if
    !go to into eV, ang and GC world............................................      

    if (rang==0)then
       F=potist*erg2eV


       IF (dmtype.EQ.30) THEN ! Variables = reduced coordinates
          do i=1,im_glob
             G(3*i-2:3*i)=-MatMul(fp_local(:,i), at)*erg2eV
          end do
       ELSE ! Variables = cartesian coordinates (in A)
          do i=1,im_glob
             G(3*i-2:3*i)=-fp_local(1:3,i)*erg2eV/angst
          end do
       END IF
       G(3*im_glob+1:N)=0.d0
    end if

    return
  end subroutine FUNCT



end module work_cgII
