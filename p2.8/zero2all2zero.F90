module zero2all2zero_mod
  use T_kind_param_m, ONLY:  double
#ifdef PARA
  use mod_para
#endif
  use gen_com_m
  use tab_imm_m
  use cryst_to_cart_mod
contains
  subroutine zero2all(vectall,vectp,itypall)
    real(double),intent(in)::vectall(3,imm_glob)
    integer,intent(in)::itypall(imm_glob)
    real(double),intent(out)::vectp(3,imm)
    !    integer,intent(out)::itypp(im_glob)

    integer::i,it,cellx,celly,cellz
    real(double)::aux,auy,auz


    !    real(double)

!    write(6,*)'tailles', size(vectall),size(vectp)
#ifdef PARA
    real(double)::vectin(3,imm_glob)
    vectin=vectall
    call cryst_to_cart (im_glob, vectin, bg, -1) !cart vers cryst
    !    if (rang==0) then !envoi vectall à tous
    call MPI_BCAST(vectin,3*imm_glob,NDM_MPI_REAL_DOUBLE,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(itypall,imm_glob,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    !    call MPI_BCAST(num_at_glob,ims,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)



     do it=1,im_glob
         do ic=1,3
	   xpici=vectin(ic,it)
	   if ( (xpici < 0.d0 ).OR.( xpici >= 1.d0 ) ) then
              if ( (xpici > -low_limit).and.(xpici<0.d0) ) then
                 vectin(ic,it)=zero
              else
                 cpp  = Dble(Floor(vectin(ic,it)))
                 vectin (ic,it) = xpici     - cpp
              end if
	   end if
         end do
      end do
!      write(6,'(A,7I4)')"rang cellxyz",rang,cell_debx,cell_finx,cell_deby,cell_finy,cell_debz,cell_finz
      i=0
    do it=1,im_glob
       
       ! On teste si c'est un atome local pour le prendre
       ! en compte ou le retirer
       aux = vectin(1,it)*nox
       auy = vectin(2,it)*noy
       auz = vectin(3,it)*noz
       cellx = int(aux)+1
       celly = int(auy)+1
       cellz = int(auz)+1
!       write(6,'(A,I2,I5,6G18.9,3I3)')'at cell',rang,it,vectin(1,it),vectin(2,it),vectin(3,it),aux,auy,auz,cellx,celly,cellz


       !       call coord_to_cellcoord(vectin(1,i),vectin(2,i),vectin(3,i),cellx,celly,cellz)
       if (cellx<cell_debx.or.cellx>cell_finx .or. &
            celly<cell_deby.or.celly>cell_finy .or. &
            cellz<cell_debz.or.cellz>cell_finz) then
          ! l'atome n'est pas local, on l'elimine du processeur courant
       else
          i=i+1
          vectp(:,i)=vectin(:,it)
          ityp(i)=itypall(it)
          num_at_glob(i)=it
       endif


    end do
!    call cryst_to_cart (im_glob, vectall, at, 1)!cryst vers cart
    im=i
    call cryst_to_cart (im, vectp, at, 1)  !cryst vers cart
        !    lenfn2=2
    !    write(6,*)rang, 'POST TRI',im_glob,im
    !    write(extension,'(i2.2)') rang
    !    open(unit=618, file='poscont.'//extension(1:lenfn2)//'.csv', form='formatted', &
    !         status='unknown')
    !     do i=1,im
    !       write(618,'(I2,2I6,3G18.9)') rang,i,num_at_glob(i),vectp(:,i)
    !    end do
    !    
    !    call mpi_barrier(MPI_COMM_WORLD,ierr)
    !    stop

!    write(6,*)"rg im",rang,im


#else

    vectp(:,:) =  vectall(:,:)
#endif

  end subroutine zero2all

  subroutine all2zero(vectall,vectp)
    real(double),intent(out)::vectall(3,imm_glob)
!    integer,intent(out)::itypall(imm_glob)
    real(double),intent(in)::vectp(3,imm)

    integer::nag(imm)
    integer::iproc, proc_source,i

#ifdef PARA

    if (rang==0) then
       do iproc=0,nprocs-1
          ! Pour le processeur maitre il n'y a rien a faire
          ! reception des donnees des autres processeurs
          if (iproc.ne.0) then
             call MPI_RECV(im,               1,    MPI_INTEGER,      MPI_ANY_SOURCE, 10001, MPI_COMM_WORLD, status, ierr)
             proc_source = status(MPI_SOURCE)
             call MPI_RECV(vectp(1:3,1:im),     3*im, NDM_MPI_REAL_DOUBLE, proc_source, 10002, MPI_COMM_WORLD, status, ierr)
             !             call MPI_RECV(num_at_glob(1:im),im,   MPI_INTEGER,         proc_source, 10004, MPI_COMM_WORLD, status, ierr)
             call MPI_RECV(nag(1:im),     im, MPI_INTEGER, proc_source, 10003, MPI_COMM_WORLD, status, ierr)
             !             call MPI_RECV(num_at_glob(1:im),im,   MPI_INTEGER,         proc_source, 10004, MPI_COMM_WORLD, status, ierr)
          else
             nag(1:im)=num_at_glob(1:im)
          end if
          do i=1,im
             num_at_glob(i)=nag(i)
             vectall(:,num_at_glob(i))=vectp(:,i)
          end do
          !          do i=1,10
          !             write(6,'(I2,2I6,3G18.9)')iproc,i, num_at_glob(i),vectp(:,i)
          !          end do
       end do

    else
       call MPI_SEND(im,               1,   MPI_INTEGER,        0,10001,MPI_COMM_WORLD,ierr)
       call MPI_SEND(vectp(1:3,1:im),     3*im,NDM_MPI_REAL_DOUBLE,0,10002,MPI_COMM_WORLD,ierr)
       call MPI_SEND(num_at_glob(1:im),     im, MPI_INTEGER, 0, 10003, MPI_COMM_WORLD, status, ierr)
       !       call MPI_SEND(ityp(1:im),       im,  MPI_INTEGER,        0,10003,MPI_COMM_WORLD,ierr)
       !       call MPI_SEND(num_at_glob(1:im),im,  MPI_INTEGER,        0,10004,MPI_COMM_WORLD,ierr)
       !      xp_all=0; vectp_all=0;ityp_all=0
       !      X=0
    end if


#else    
    vectall(:,:) =  vectp(:,:)
#endif
  end subroutine all2zero

end module zero2all2zero_mod
