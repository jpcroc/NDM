module sauvegardeT_mod

    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:rang,formatsauv,im_glob,it,itesauvinter,&
         &pmean,rang,timel,tmean,tstep,fnam,lenfnam,lcasca,imm_glob,l2T

    USE elec_cell, ONLY : sauveelec
    USE cryst_to_cart_mod,only: cryst_to_cart
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE cellconfig, only:cell_config
  USE boxconfig,only:box_config
#ifdef PARA

    USE mpi
    USE mod_para,only:MPI_COMM_space,status,ierr,nprocspace,myid,NDM_MPI_REAl_DOUBLE
#else
USE mod_para,only:nprocspace    ,myid
         
#endif



    implicit none

contains
  ! ********************************************************************
  subroutine sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout)
    !-----------------------------------------------
    !   M o d u l e s



    implicit none
    type(box_config)::boxndm
    class(atom_config)::atdml
    type(cell_config):: celndm

    integer :: lucout, formatsauvmod,i,formatsauv,im
    character :: extension*9
    logical :: lwax
    character::fnamcout*80
#ifdef PARA
    integer,dimension(:),allocatable     :: ibuffer
    real(double), dimension(:,:),allocatable   :: buffer
    integer,      dimension(0:nprocspace-1)   :: im_loc
    integer,      dimension(0:nprocspace-1)   :: pt_im
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
    im =atdml%im
 !      write (6, *) ' sauvegarde it=', it, rang,fnamcout
    if (myid==0) then

       lucout = 87
!       write (6, *) ' sauvegarde it=', it, rang,fnamcout
       open(unit=lucout, file=fnamcout, form='unformatted', status='unknown')
       write (lucout) formatsauv
       write (lucout) boxndm%at
       write (lucout) im_glob
#ifdef PARA
       if (nprocspace.gt.1) then
          im_loc(0)=im
          ibuffer=0
          ibuffer(1:im)  = atdml%ityp(1:im)
          buffer=0
          buffer(:,1:im) = atdml%xp(:,1:im)
          pt_im(0)=1
          next_pt = pt_im(0) + im_loc(0)
          
          do i_proc=1,nprocspace-1
             call MPI_RECV(im_temp,1, MPI_INTEGER, MPI_ANY_SOURCE, 11001, MPI_COMM_space, status, ierr)
             
             proc_source = status(MPI_SOURCE)
             im_loc(proc_source)=im_temp
             pt_im(proc_source)=next_pt
             next_pt = pt_im(proc_source) + im_loc(proc_source)
             call MPI_RECV(ibuffer(pt_im(proc_source):pt_im(proc_source)+im_temp-1),    im_loc(proc_source),   &
                  MPI_INTEGER,         proc_source, 11002, MPI_COMM_space, status, ierr)
             call MPI_RECV(buffer(1:3,pt_im(proc_source):pt_im(proc_source)+im_temp-1),3*im_loc(proc_source), &
                  NDM_MPI_REAL_DOUBLE, proc_source, 11003, MPI_COMM_space, status, ierr)
          enddo
          write (lucout) ibuffer  ! Ecriture ityp
          write (lucout) buffer   ! Ecriture xp
          
          ibuffer(1:im) = atdml%num_at_glob(1:im)
          do i_proc=1,nprocspace-1
             call MPI_RECV(ibuffer(pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),im_loc(i_proc), &
                  MPI_INTEGER, i_proc, 11004, MPI_COMM_space, status, ierr)
          enddo
          write (lucout) ibuffer   ! Ecriture num_at_glob
          
          if (formatsauvmod==1) then
             lwax=.false.
             select type(atdml)
             type is (atom_config_d)
                buffer(:,1:im) = atdml%xpp(:,1:im)
                do i_proc=1,nprocspace-1
                   call MPI_RECV(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),3*im_loc(i_proc), &
                        NDM_MPI_REAL_DOUBLE, i_proc, 11005, MPI_COMM_space, status, ierr)
                enddo
                write (lucout) buffer   ! Ecriture xpp
                
                buffer(:,1:im) = atdml%vp(:,1:im)
                do i_proc=1,nprocspace-1
                   call MPI_RECV(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),3*im_loc(i_proc), &
                        NDM_MPI_REAL_DOUBLE, i_proc, 11006, MPI_COMM_space, status, ierr)
                enddo
                write (lucout) buffer   ! Ecriture vp
             type is (atom_config_e)
                buffer(:,1:im) = atdml%xpp(:,1:im)
                do i_proc=1,nprocspace-1
                   call MPI_RECV(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),3*im_loc(i_proc), &
                        NDM_MPI_REAL_DOUBLE, i_proc, 11005, MPI_COMM_space, status, ierr)
                enddo
                write (lucout) buffer   ! Ecriture xpp
                
                buffer(:,1:im) = atdml%vp(:,1:im)
                do i_proc=1,nprocspace-1
                   call MPI_RECV(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),3*im_loc(i_proc), &
                        NDM_MPI_REAL_DOUBLE, i_proc, 11006, MPI_COMM_space, status, ierr)
                enddo
                write (lucout) buffer   ! Ecriture vp
                if (atdml%lax)then
                   buffer(:,1:im) = atdml%ax(:,1:im)
                   do i_proc=1,nprocspace-1
                      call MPI_RECV(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),3*im_loc(i_proc), &
                           NDM_MPI_REAL_DOUBLE, i_proc, 11007, MPI_COMM_space, status, ierr)
                   enddo
                   write (lucout) buffer   ! Ecriture ax
                end if
             end select
             if (.not.lwax)then
                buffer(:,1:im) = atdml%xp(:,1:im)
                do i_proc=1,nprocspace-1
                   call MPI_RECV(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),3*im_loc(i_proc), &
                        NDM_MPI_REAL_DOUBLE, i_proc, 11008, MPI_COMM_space, status, ierr)
                enddo
                write (lucout) buffer   ! Ecriture xpp
                
             end if
             write (lucout) tstep
             write (lucout) tmean, pmean, it, timel
          endif
       else
                 write (lucout) atdml%ityp
       write (lucout) atdml%xp
       write (lucout) atdml%num_at_glob
       if (formatsauvmod==1) then
          lwax=.false.
          select type (atdml)
          type is (atom_config_d)
             write (lucout) atdml%xpp
             write (lucout) atdml%vp
          type is (atom_config_e)
             write (lucout) atdml%xpp
             write (lucout) atdml%vp
             if (atdml%lax)then
                write (lucout) atdml%ax
                lwax=.true.
             end if
          end select
          if (.not.lwax)write (lucout) atdml%xp
          write (lucout) tstep
          write (lucout) tmean, pmean, it, timel
       endif
    end if
#else
       !
       ! Partie sequentielle de la sauvegarde :
       !
       !     do i=1,im
       !        write(1004,*)i,num_at_glob(i),xp(1,i)
       !     enddo
       write (lucout) atdml%ityp
       write (lucout) atdml%xp
       write (lucout) atdml%num_at_glob
       if (formatsauvmod==1) then
          lwax=.false.
          select type (atdml)
          type is (atom_config_d)
             write (lucout) atdml%xpp
             write (lucout) atdml%vp
          type is (atom_config_e)
             write (lucout) atdml%xpp
             write (lucout) atdml%vp
             if (atdml%lax)then
                write (lucout) atdml%ax
                lwax=.true.
             end if
          end select
          if (.not.lwax)write (lucout) atdml%xp
          write (lucout) tstep
          write (lucout) tmean, pmean, it, timel
       endif

#endif

       close(unit=lucout)

       if (l2T)call sauveelec

    else ! rang different de 0 :
#ifdef PARA
          if (nprocspace.gt.1) then
             
             call MPI_SEND(im,          1,   MPI_INTEGER,        0,11001,MPI_COMM_space,ierr)
             call MPI_SEND(atdml%ityp(1:im),  im,  MPI_INTEGER,        0,11002,MPI_COMM_space,ierr)
             call MPI_SEND(atdml%xp(1:3,1:im),3*im,NDM_MPI_REAL_DOUBLE,0,11003,MPI_COMM_space,ierr)
             call MPI_SEND(atdml%num_at_glob(1:im), im,    MPI_INTEGER,0,11004,MPI_COMM_space,ierr)
             
             if (formatsauvmod==1) then
                lwax=.false.
                select type (atdml)
                type is (atom_config_d)
                   call MPI_SEND(atdml%xpp(1:3,1:im),3*im,NDM_MPI_REAL_DOUBLE,0,11005,MPI_COMM_space,ierr)
                   call MPI_SEND(atdml%vp(1:3,1:im), 3*im,NDM_MPI_REAL_DOUBLE,0,11006,MPI_COMM_space,ierr)
                type is (atom_config_e)
                   if (atdml%lax)then
                      lwax=.true.
                      call MPI_SEND(atdml%ax(1:3,1:im), 3*im,NDM_MPI_REAL_DOUBLE,0,11007,MPI_COMM_space,ierr)
                   end if
                end select
                if (.not.lwax)call MPI_SEND(atdml%xp(1:3,1:im), 3*im,NDM_MPI_REAL_DOUBLE,0,11008,MPI_COMM_space,ierr)
             endif
          end if
#endif
       endif
#ifdef PARA
              if (nprocspace.gt.1) then

                 deallocate (buffer)
                 deallocate (ibuffer)
              end if
#endif

!    write (6, *) ' OUT sauvegarde it=', it,rang, myid,fnamcout
    return
  end subroutine sauvegardeT

  subroutine cin2gin
    USE gen_com_m, ONLY:at,bg,im,im_glob,imm,rang,at,fnamcout,formatsauv,im_glob,it,itesauvinter,&
         &pmean,rang,timel,tmean,tstep,fnam,lenfnam,lcasca,imm_glob,l2T
    USE tab_imm_m,only:xp,vp,fp,xpp,ax,ityp,num_at_glob
    integer :: lugout,i
    character :: extension*9
    character :: fnamgout*80
#ifdef PARA
    integer,dimension(:),allocatable     :: ibuffer
    real(double), dimension(:,:),allocatable   :: buffer
    integer,      dimension(0:nprocspace-1)   :: im_loc
    integer,      dimension(0:nprocspace-1)   :: pt_im
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
       do i_proc=1,nprocspace-1
          call MPI_RECV(im_temp,1, MPI_INTEGER, MPI_ANY_SOURCE, 11001, MPI_COMM_space, status, ierr)
          proc_source = status(MPI_SOURCE)
          im_loc(proc_source)=im_temp
          pt_im(proc_source)=next_pt
          next_pt = pt_im(proc_source) + im_loc(proc_source)
          call MPI_RECV(ibuffer(pt_im(proc_source):pt_im(proc_source)+im_temp-1),    im_loc(proc_source),   &
               MPI_INTEGER,         proc_source, 11002, MPI_COMM_space, status, ierr)
          call MPI_RECV(buffer(1:3,pt_im(proc_source):pt_im(proc_source)+im_temp-1),3*im_loc(proc_source), &
               NDM_MPI_REAL_DOUBLE, proc_source, 11003, MPI_COMM_space, status, ierr)
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
       call MPI_SEND(im,          1,   MPI_INTEGER,        0,11001,MPI_COMM_space,ierr)

       call MPI_SEND(ityp(1:im),  im,  MPI_INTEGER,        0,11002,MPI_COMM_space,ierr)
       call MPI_SEND(xp(1:3,1:im),3*im,NDM_MPI_REAL_DOUBLE,0,11003,MPI_COMM_space,ierr)
#endif
    endif
#ifdef PARA
    deallocate (buffer)
    deallocate (ibuffer)
#endif

    return
  end subroutine cin2gin



    
end module sauvegardeT_mod
