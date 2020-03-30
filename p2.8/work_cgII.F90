module work_cgII

  USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: im, imm,at, inv_angst, lperiod, rang, &
                       it, itesauv, itesauvposition, itesauvforce, &
                       inv_angst, erg2ev, angst, &
                       dmtype, potist,im_glob,nox,noy,noz,cell_finx,cell_finy,cell_finz,noxyz
  use controle_mod
  use calfo_mod
  use analyse_mod
  use sauvegarde_mod
  use sauveposition_mod
  use sauveforce_mod
  use caltabt_mod
  use config_mod
  use zero2all2zero_mod
#ifdef PARA
  use mod_para
#endif
  use tab_imm_m, only : xp, fp,num_at_glob
  
  implicit none


contains

  subroutine FUNCT(N,X,F,G,NCALLS,                      &
       xp_local,  fp_local,   ityp_local,ims)
    use tab_imm_m, only : xp, fp,num_at_glob
    double precision X(N),G(N),F
    integer i,N,NCALLS,ims
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: ityp_local(ims)
    real(double)  :: xp_local(3,ims)
    real(double)  :: fp_local(3,ims)
    !-----------------------------------------------
    integer::iproc,proc_source,cellx,celly,cellz
    integer::nag(imm)
    real(double)::aux,auy,auz
    character :: extension*2
    integer::lenfn2,ko
    real(double) :: fpmax,fpn,forctot,formax,fpmax_glob

    IF (3*ims.NE.N) THEN
       WRITE(0,'(a,i0)') "3*imm = ", 3*ims
       WRITE(0,'(a,i0)') "N     = ", N
       STOP "< work_cgII >"
    END IF
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


call zero2all(xp_local,xp,ityp_local)
    
!#ifdef PARA
!    !    if (rang==0) then !envoi xp_local à tous!
!    call MPI_BCAST(xp_local,3*ims,NDM_MPI_REAL_DOUBLE,0,MPI_COMM_WORLD,ierr)
!    call MPI_BCAST(ityp_local,ims,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!    call MPI_BCAST(num_at_glob,ims,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    write(6,*)rang, 'POST BCAST',im_glob
!    write(extension,'(i2.2)') rang
!    lenfn2=2
!    open(unit=606, file='poB.'//extension(1:lenfn2)//'.csv', form='formatted', &
!             status='unknown')
!    do i=1,im_glob
!       !       write(6,*)rang,i,xp_all(:,i)
!       write(606,'(I2,I6,3G18.9)') rang,i,xp_local(:,i)
!    end do
!    call mpi_barrier(MPI_COMM_WORLD,ierr)


!    write(6,*)'rg cell',rang,cell_finx,cell_finy,cell_finz
!    i=0
!    call cryst_to_cart (im_glob, xp_local, bg, -1) !cart vers cryst
!    do it=1,im_glob
!       ! On teste si c'est un atome local pour le prendre
       ! en compte ou le retirer
!       aux = xp_local(1,it)*nox
!       auy = xp_local(2,it)*noy
!       auz = xp_local(3,it)*noz
!       cellx = int(aux)+1
!       celly = int(auy)+1
!       cellz = int(auz)+1
!       write(6,'(A,I2,I5,6G18.9,3I3)')'at cell',rang,it,xp_local(1,it),xp_local(2,it),xp_local(3,it),aux,auy,auz,cellx,celly,cellz
!       call coord_to_cellcoord(xp_local(1,i),xp_local(2,i),xp_local(3,i),cellx,celly,cellz)
!       if (cellx<cell_debx.or.cellx>cell_finx .or. &
!            celly<cell_deby.or.celly>cell_finy .or. &
!            cellz<cell_debz.or.cellz>cell_finz) then
!          ! l'atome n'est pas local, on l'elimine du processeur courant
!       else
!          i=i+1
!          xp(:,i)=xp_local(:,it)
!          ityp(i)=ityp_local(it)
!          num_at_glob(i)=it
!       endif
!    end do
!    call cryst_to_cart (im_glob, xp_local, at, 1)!cryst vers cart
!    call cryst_to_cart (im_glob, xp, at, 1)  !cryst vers cart
!    im=i
!!    lenfn2=2
!!    write(6,*)rang, 'POST TRI',im_glob,im
!    write(extension,'(i2.2)') rang
!    open(unit=618, file='poscont.'//extension(1:lenfn2)//'.csv', form='formatted', &
!         status='unknown')
!     do i=1,im
!       write(618,'(I2,2I6,3G18.9)') rang,i,num_at_glob(i),xp(:,i)
!    end do
!    
!    call mpi_barrier(MPI_COMM_WORLD,ierr)
!    stop
!
!    
!   call caltabt
!    do ko=1,noxyz
!       write(6,*)'0rg cel nat',rang, ko,nato(ko)
!    end do
!    
!    call maj_atomes_frt_ftm
!#else
!
!    xp(:,:) =  xp_local(:,:)
!#endif
    !back to internal units and JP world.......................................

    it=NCALLS-1

    if (lperiod)          call period 
!    call controle
    

!    write(*,*) 'inside FUNCT debug_incg1', xp_local(1,1), fp_local(1,1)

!do ko=1,noxyz
!    write(6,*)'1rg cel nat',rang, ko,nato(ko)
!end do
    call calfo

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
!       write(625,'(2I6,3G18.9)') i,num_at_glob(i),fp(:,i)
!    end do
!#ifdef PARA
!    call mpi_barrier(MPI_COMM_WORLD,ierr)
!#endif    
!   stop

    !do ko=1,noxyz
!    write(6,*)'2rg cel 1nat',rang, ko,nato(ko)
! end do

    
  if (it>=itmax) then
     if (rang==0) write (6, *) '*******Derniere iteration **** '
     call endrun
     call DeallocateAll

     call arret_ndm

  endif

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

     IF (it.GE.1) THEN
        IF (lFrozen) THEN
           !forctot=sqrt( Sum( SUM(fp(1:3,1:im)**2,1), Free(1:im) ) )
           !formax=sqrt( MAXVAL( Sum(fp(1:3,1:im)**2,1), Free(1:im) ) )
           forctot = sqrt( SUM( fp(:,1:im)**2, .NOT.Frozen(:,1:im) ) )
           formax = MaxVal( Abs(fp(:,1:im)), .NOT.Frozen(:,1:im) ) 
        ELSE
           forctot=sqrt( SUM(fp(1:3,1:im)**2) )
           !formax=sqrt( MAXVAL( Sum(fp(1:3,1:im)**2,1) ) )
           formax = MaxVal( Abs(fp(:,1:im)) )
        END IF
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




!do ko=1,noxyz
!   write(6,*)'3rg cel nat',rang, ko,nato(ko)
!end do
!    lenfn2=2
!    write(6,*)rang, 'Fpcf',im_glob,im
!    write(extension,'(i2.2)') rang
!#ifdef PARA    
!    open(unit=625, file='FORC.'//extension(1:lenfn2)//'.csv', form='formatted', &
!         status='unknown')
!
!#else    
!    open(unit=625, file='FORC.SEQ.csv', form='formatted', &
!         status='unknown')
!#endif
!    do i=1,im
!       write(625,'(I2,2I6,3G18.9)') rang,i,num_at_glob(i),fp(:,i)
!    end do
!#ifdef PARA
!    call mpi_barrier(MPI_COMM_WORLD,ierr)
!#endif    
!    stop
    
!    lenfn2=2
!    write(6,*)rang, 'POST calf',im_glob,im
!    write(extension,'(i2.2)') rang
!    open(unit=608, file='poscalf2.'//extension(1:lenfn2)//'.csv', form='formatted', &
!         status='unknown')
!    do i=1,im_glob
!       write(608,'(I2,I6,3G18.9)') rang,i,xp_local(:,i)
!    end do
!    call mpi_barrier(MPI_COMM_WORLD,ierr)
!    stop

     call all2zero(fp_local,fp)
     
!#ifdef PARA
!    if (rang==0) then
!       do iproc=0,nprocs-1
!          ! Pour le processeur maitre il n'y a rien a faire
!          ! reception des donnees des autres processeurs
!          if (iproc.ne.0) then
!             call MPI_RECV(im,               1,    MPI_INTEGER,      MPI_ANY_SOURCE, 10001, MPI_COMM_WORLD, status, ierr)
!             proc_source = status(MPI_SOURCE)
!             call MPI_RECV(fp(1:3,1:im),     3*im, NDM_MPI_REAL_DOUBLE, proc_source, 10002, MPI_COMM_WORLD, status, ierr)
!!             call MPI_RECV(num_at_glob(1:im),im,   MPI_INTEGER,         proc_source, 10004, MPI_COMM_WORLD, status, ierr)
!             call MPI_RECV(nag(1:im),     im, MPI_INTEGER, proc_source, 10003, MPI_COMM_WORLD, status, ierr)
!             !             call MPI_RECV(num_at_glob(1:im),im,   MPI_INTEGER,         proc_source, 10004, MPI_COMM_WORLD, status, ierr)
!          else
!             nag(1:im)=num_at_glob(1:im)
!          end if
!          do i=1,im
!             num_at_glob(i)=nag(i)
!             fp_local(:,num_at_glob(i))=fp(:,i)
!          end do
!!          do i=1,10
!!             write(6,'(I2,2I6,3G18.9)')iproc,i, num_at_glob(i),fp(:,i)
!!          end do
!       end do
!    else
!       call MPI_SEND(im,               1,   MPI_INTEGER,        0,10001,MPI_COMM_WORLD,ierr)
!       call MPI_SEND(fp(1:3,1:im),     3*im,NDM_MPI_REAL_DOUBLE,0,10002,MPI_COMM_WORLD,ierr)
!       call MPI_SEND(num_at_glob(1:im),     im, MPI_INTEGER, 0, 10003, MPI_COMM_WORLD, status, ierr)
!       !       call MPI_SEND(ityp(1:im),       im,  MPI_INTEGER,        0,10003,MPI_COMM_WORLD,ierr)
!       !       call MPI_SEND(num_at_glob(1:im),im,  MPI_INTEGER,        0,10004,MPI_COMM_WORLD,ierr)
!       !      xp_all=0; fp_all=0;ityp_all=0
!       !      X=0
!    end if
!#else    
!    fp_local(:,:) =  fp(:,:)
!#endif

!    lenfn2=2
!    write(6,*)rang, 'F1',im_glob,im
!    write(extension,'(i2.2)') rang
!#ifdef PARA    
!    open(unit=620, file='FO1.'//extension(1:lenfn2)//'.csv', form='formatted', &
!         status='unknown')
!    open(unit=621, file='FOT.'//extension(1:lenfn2)//'.csv', form='formatted', &
!         status='unknown')
!#else    
!    open(unit=620, file='FO1.SEQ.csv', form='formatted', &
!         status='unknown')
!    open(unit=621, file='FOT.SEQ.csv', form='formatted', &
!         status='unknown')
!#endif
!    open(unit=620, file='Fo1.'//extension(1:lenfn2)//'.csv', form='formatted', &
!         status='unknown')
!    do i=1,im
!       write(620,'(2I6,3G18.9)') i,num_at_glob(i),fp(:,i)
!    end do
!        open(unit=621, file='FoT.'//extension(1:lenfn2)//'.csv', form='formatted', &
!         status='unknown')
!    do i=1,im_glob
!       write(621,'(I6,3G18.9)') i,fp_local(:,i)
!   end do
!#ifdef PARA    
!    call mpi_barrier(MPI_COMM_WORLD,ierr)
!#endif    
!    stop


    !write(*,*) 'inside FUNCT debug_incg2', xp_local(1,1), fp_local(1,1)
    call analyse 
    !write(*,*) 'inside FUNCT debug_incg3', xp_local(1,1), fp_local(1,1)


    if (it.ne.0) then
       if (rang==0) then
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
       endif                                   ! fin rang=0
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
