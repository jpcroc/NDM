module sauvegarde_mod

    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:at,bg,im,im_glob,imm,rang,at,fnamcout,formatsauv,im_glob,it,itesauvinter,&
         &pmean,rang,timel,tmean,tstep,fnam,lenfnam,lcasca,imm_glob
    USE tab_imm_m
    USE elec_cell, ONLY : sauveelec
    USE cryst_to_cart_mod,only: cryst_to_cart
#ifdef PARA
    USE mod_para
#endif



    implicit none

contains
  ! ********************************************************************
  subroutine sauvegarde
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    !         version du 29 novembre 2000
    ! ********************************************************************


    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: lucout, formatsauvmod,i
    character :: extension*9
#ifdef PARA
    integer,dimension(:),allocatable     :: ibuffer
    real(double), dimension(:,:),allocatable   :: buffer
    integer,      dimension(0:nprocs-1)   :: im_loc
    integer,      dimension(0:nprocs-1)   :: pt_im
    integer :: next_pt
    integer :: i_proc
    integer :: proc_source
    integer :: im_temp

    allocate (buffer(3,imm_glob))
    allocate (ibuffer(imm_glob))

#endif
    !-----------------------------------------------
    !      include 'pot.com'
    !      include 'potdyn.com'
    !      include 'dyn.com'
    !    open du fichier .cout

    formatsauvmod = mod(formatsauv,2)

    if (rang==0) then
       if(itesauvinter.gt.0) then
          if (mod(it,itesauvinter).eq.0) then
             write(extension,'(i9.9)') it
             fnamcout = fnam(1:lenfnam)//'.cout.'//extension
          else
             fnamcout = fnam(1:lenfnam)//'.cout'
          endif
       else
          fnamcout = fnam(1:lenfnam)//'.cout'
       end if
       if((it==0).and.lcasca)       fnamcout = fnam(1:lenfnam)//'.0.cout'

       lucout = 87
       write (6, *) ' sauvegarde it=', it, fnamcout
       open(unit=lucout, file=fnamcout, form='unformatted', status='unknown')
       write (lucout) formatsauv
       write (lucout) at
       write (lucout) im_glob

#ifdef PARA
       im_loc(0)=im
       ibuffer=0
       ibuffer(1:im)  = ityp(1:im)
       buffer=0
       buffer(:,1:im) = xp(:,1:im)
       pt_im(0)=1
       next_pt = pt_im(0) + im_loc(0)
       do i_proc=1,nprocs-1
          call MPI_RECV(im_temp,1, MPI_INTEGER, MPI_ANY_SOURCE, 11001, MPI_COMM_WORLD, status, ierr)
          proc_source = status(MPI_SOURCE)
          im_loc(proc_source)=im_temp
          pt_im(proc_source)=next_pt
          next_pt = pt_im(proc_source) + im_loc(proc_source)
          call MPI_RECV(ibuffer(pt_im(proc_source):pt_im(proc_source)+im_temp-1),    im_loc(proc_source),   &
               MPI_INTEGER,         proc_source, 11002, MPI_COMM_WORLD, status, ierr)
          call MPI_RECV(buffer(1:3,pt_im(proc_source):pt_im(proc_source)+im_temp-1),3*im_loc(proc_source), &
               NDM_MPI_REAL_DOUBLE, proc_source, 11003, MPI_COMM_WORLD, status, ierr)
       enddo
       write (lucout) ibuffer  ! Ecriture ityp
       write (lucout) buffer   ! Ecriture xp


       ibuffer(1:im) = num_at_glob(1:im)
       do i_proc=1,nprocs-1
          call MPI_RECV(ibuffer(pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),im_loc(i_proc), &
               MPI_INTEGER, i_proc, 11004, MPI_COMM_WORLD, status, ierr)
       enddo
       write (lucout) ibuffer   ! Ecriture num_at_glob

       if (formatsauvmod==1) then
          buffer(:,1:im) = xpp(:,1:im)
          do i_proc=1,nprocs-1
             call MPI_RECV(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),3*im_loc(i_proc), &
                  NDM_MPI_REAL_DOUBLE, i_proc, 11005, MPI_COMM_WORLD, status, ierr)
          enddo
          write (lucout) buffer   ! Ecriture xpp

          buffer(:,1:im) = vp(:,1:im)
          do i_proc=1,nprocs-1
             call MPI_RECV(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),3*im_loc(i_proc), &
                  NDM_MPI_REAL_DOUBLE, i_proc, 11006, MPI_COMM_WORLD, status, ierr)
          enddo
          write (lucout) buffer   ! Ecriture vp

          buffer(:,1:im) = ax(:,1:im)
          do i_proc=1,nprocs-1
             call MPI_RECV(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),3*im_loc(i_proc), &
                  NDM_MPI_REAL_DOUBLE, i_proc, 11007, MPI_COMM_WORLD, status, ierr)
          enddo
          write (lucout) buffer   ! Ecriture ax

          write (lucout) tstep
          write (lucout) tmean, pmean, it, timel
       endif
#else
       !
       ! Partie sequentielle de la sauvegarde :
       !
       !     do i=1,im
       !        write(1004,*)i,num_at_glob(i),xp(1,i)
       !     enddo
       write (lucout) ityp
       write (lucout) xp
       write (lucout) num_at_glob
       if (formatsauvmod==1) then
          write (lucout) xpp
          write (lucout) vp
          write (lucout) ax
          write (lucout) tstep
          write (lucout) tmean, pmean, it, timel
       endif
#endif

       close(unit=lucout)

       if (l2T)call sauveelec

    else ! rang different de 0 :
#ifdef PARA
       call MPI_SEND(im,          1,   MPI_INTEGER,        0,11001,MPI_COMM_WORLD,ierr)
       call MPI_SEND(ityp(1:im),  im,  MPI_INTEGER,        0,11002,MPI_COMM_WORLD,ierr)
       call MPI_SEND(xp(1:3,1:im),3*im,NDM_MPI_REAL_DOUBLE,0,11003,MPI_COMM_WORLD,ierr)
       call MPI_SEND(num_at_glob(1:im), im,    MPI_INTEGER,0,11004,MPI_COMM_WORLD,ierr)
       if (formatsauvmod==1) then
          call MPI_SEND(xpp(1:3,1:im),3*im,NDM_MPI_REAL_DOUBLE,0,11005,MPI_COMM_WORLD,ierr)
          call MPI_SEND(vp(1:3,1:im), 3*im,NDM_MPI_REAL_DOUBLE,0,11006,MPI_COMM_WORLD,ierr)
          call MPI_SEND(ax(1:3,1:im), 3*im,NDM_MPI_REAL_DOUBLE,0,11007,MPI_COMM_WORLD,ierr)
       endif
#endif
    endif
#ifdef PARA
    deallocate (buffer)
    deallocate (ibuffer)
#endif

    !  if (rang==0)write(6,*)'fin sauvegarde'
    return
  end subroutine sauvegarde

  subroutine cin2gin
    integer :: lugout,i
    character :: extension*9
    character :: fnamgout*80
#ifdef PARA
    integer,dimension(:),allocatable     :: ibuffer
    real(double), dimension(:,:),allocatable   :: buffer
    integer,      dimension(0:nprocs-1)   :: im_loc
    integer,      dimension(0:nprocs-1)   :: pt_im
    integer :: next_pt
    integer :: i_proc
    integer :: proc_source
    integer :: im_temp

    allocate (buffer(3,imm_glob))
    allocate (ibuffer(imm_glob))

#endif

    if (rang==0) then
       fnamgout = fnam(1:lenfnam)//'.newgin'



       lugout = 877
       write (6, *) ' ecriture newgin', fnamgout
       open(unit=lugout, file=fnamgout, form='formatted', status='unknown')
       write(lugout,*)'1 1 1 '
       at=at*1d8
       write(lugout,*)at(1,1),at(2,1),at(3,1) !a
       write(lugout,*)at(1,2),at(2,2),at(3,2) !b
       write(lugout,*)at(1,3),at(2,3),at(3,3) !c
       at=at*1d-8
      write(lugout,*)im_glob


#ifdef PARA
       im_loc(0)=im
       ibuffer=0
       ibuffer(1:im)  = ityp(1:im)
       buffer=0
       buffer(:,1:im) = xp(:,1:im)
       pt_im(0)=1
       next_pt = pt_im(0) + im_loc(0)
       do i_proc=1,nprocs-1
          call MPI_RECV(im_temp,1, MPI_INTEGER, MPI_ANY_SOURCE, 11001, MPI_COMM_WORLD, status, ierr)
          proc_source = status(MPI_SOURCE)
          im_loc(proc_source)=im_temp
          pt_im(proc_source)=next_pt
          next_pt = pt_im(proc_source) + im_loc(proc_source)
          call MPI_RECV(ibuffer(pt_im(proc_source):pt_im(proc_source)+im_temp-1),    im_loc(proc_source),   &
               MPI_INTEGER,         proc_source, 11002, MPI_COMM_WORLD, status, ierr)
          call MPI_RECV(buffer(1:3,pt_im(proc_source):pt_im(proc_source)+im_temp-1),3*im_loc(proc_source), &
               NDM_MPI_REAL_DOUBLE, proc_source, 11003, MPI_COMM_WORLD, status, ierr)
       enddo
       call cryst_to_cart(im_glob,buffer,bg,-1) !cart vers cryst

       do i=1,im_glob
          write(lugout,'(3F18.11,I6)')buffer(1,i),buffer(2,i),buffer(3,i),ibuffer(i) !coordonnes reduites des
       enddo
#else
       !
       ! Partie sequentielle de la sauvegarde :
       call cryst_to_cart(imm,xp,bg,-1) !cart vers cryst
       do i=1,im
!          write(6,*)xp(1,i),xp(2,i),xp(3,i),ityp(i) !coordonnes reduites des
          write(lugout,'(3F21.11,I6)')xp(1,i),xp(2,i),xp(3,i),ityp(i) !coordonnes reduites des
       enddo
       call cryst_to_cart(imm,xp,at,1) !cart vers cryst
#endif

       close(unit=lugout)


    else ! rang different de 0 :
#ifdef PARA
       call MPI_SEND(im,          1,   MPI_INTEGER,        0,11001,MPI_COMM_WORLD,ierr)
       call MPI_SEND(ityp(1:im),  im,  MPI_INTEGER,        0,11002,MPI_COMM_WORLD,ierr)
       call MPI_SEND(xp(1:3,1:im),3*im,NDM_MPI_REAL_DOUBLE,0,11003,MPI_COMM_WORLD,ierr)
#endif
    endif
#ifdef PARA
    deallocate (buffer)
    deallocate (ibuffer)
#endif

    !  if (rang==0)write(6,*)'fin sauvegarde'
    return
  end subroutine cin2gin



    
end module sauvegarde_mod
