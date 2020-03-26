module work_cgII

  USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: im, imm,at, inv_angst, lperiod, rang, &
                       it, itesauv, itesauvposition, itesauvforce, &
                       inv_angst, erg2ev, angst, &
                       dmtype, potist,im_glob
  use controle_mod
  use calfo_mod
  use analyse_mod
  use sauvegarde_mod
  use sauveposition_mod
  use sauveforce_mod
  use caltabt_mod
  implicit none


contains

  subroutine FUNCT(N,X,F,G,NCALLS,                      &
       xp_local,  fp_local,   ityp_local,ims)
    use tab_imm_m, only : xp, fp,num_at_glob
    real(double),intent(in):: X(N)
    real(double),intent(out)::G(N),F
    integer,intent(in) ::N,NCALLS,ims
    integer ::i
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer,intent(in)  :: ityp_local(ims)
    real(double),intent(out)  :: xp_local(3,ims)
    real(double),intent(out)  :: fp_local(3,ims)
    !-----------------------------------------------

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

#ifdef PARA
    !    if (rang==0) then !envoi xp_local à tous
    call MPI_BCAST(xp_local,3*ims,MPI_REAL_DOUBLE,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ityp_local,ims,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

    i=0
    do it=1,im_glob

       ! On teste si c'est un atome local pour le prendre
       ! en compte ou le retirer
       call coord_to_cellcoord(xp(1,i),xp(2,i),xp(3,i),cellx,celly,cellz)
       if (cellx<cell_debx.or.cellx>cell_finx .or. &
            celly<cell_deby.or.celly>cell_finy .or. &
            cellz<cell_debz.or.cellz>cell_finz) then
          ! l'atome n'est pas local, on l'elimine du processeur courant
       else
          i=i+1
          xp(:,i)=xp_local(:,it)
          ityp(i)=ityp_local(it)
          num_at_glob(i)=it
       endif


    end do
    call caltabt
    call maj_atomes_frt_ftm

#else

    xp(:,:) =  xp_local(:,:)
#endif
    !back to internal units and JP world.......................................

    it=NCALLS-1

    if (lperiod)          call period 

    call controle 
    !write(*,*) 'inside FUNCT debug_incg1', xp_local(1,1), fp_local(1,1)

    call calfo

#ifdef PARA

    if (rang==0) then
       do iproc=0,nprocs-1
          ! Pour le processeur maitre il n'y a rien a faire
          ! reception des donnees des autres processeurs
          if (iproc.ne.0) then
             call MPI_RECV(im,               1,    MPI_INTEGER,      MPI_ANY_SOURCE, 10001, MPI_COMM_WORLD, status, ierr)
             proc_source = status(MPI_SOURCE)
             call MPI_RECV(fp(1:3,1:im),     3*im, NDM_MPI_REAL_DOUBLE, proc_source, 10002, MPI_COMM_WORLD, status, ierr)
             !             call MPI_RECV(num_at_glob(1:im),im,   MPI_INTEGER,         proc_source, 10004, MPI_COMM_WORLD, status, ierr)
          endif
          do i=1,im
             fp_local(:,num_at_glob(i))=fp(:,i)
          end do
       end do

    else
       call MPI_SEND(im,               1,   MPI_INTEGER,        0,10001,MPI_COMM_WORLD,ierr)
       call MPI_SEND(fp(1:3,1:im),     3*im,NDM_MPI_REAL_DOUBLE,0,10002,MPI_COMM_WORLD,ierr)
       !       call MPI_SEND(ityp(1:im),       im,  MPI_INTEGER,        0,10003,MPI_COMM_WORLD,ierr)
       !       call MPI_SEND(num_at_glob(1:im),im,  MPI_INTEGER,        0,10004,MPI_COMM_WORLD,ierr)
       !      xp_all=0; fp_all=0;ityp_all=0
       !      X=0
    end if


#else    
    fp_local(:,:) =  fp(:,:)
#endif    
    call analyse 


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
