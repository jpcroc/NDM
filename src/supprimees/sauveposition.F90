module sauveposition_mod
  USE gen_com_m, ONLY:fnamcoutxp,im_glob,lcasca,rang,fnam,lenfnam,imm_glob
  USE temp_com,only:at,im

  implicit none
contains 
  ! ********************************************************************
  subroutine sauveposition(itapp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE tab_imm_m,only:xp,ityp,num_at_glob
    USE arret_ndm_mod,only: arret_ndm
#ifdef PARA
  use mpi
  USE Tpara,only:MPI_COMM_space,ierr,NDM_MPI_REAL_DOUBLE,status,nprocs

#endif

    !         version du 04 octobre 2000
    ! ********************************************************************


    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer :: itapp ! iteration apparente

    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: lucoutxp, formatsauvT,rgloc
    character :: extension*9
    logical::lcrcin
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
    rgloc=rang
#else
    rgloc=0
#endif    
    !-----------------------------------------------
    !


    ! conversion entier-->alphanumerique par transfert du nombre
    ! de l'iteration vers fichier tampon relu sous format caractere.

    lucoutxp = 97+rang
    if(rgloc==0) then
       if (itapp==0) then
          inquire (file=fnam(1:lenfnam)//'.crcin',EXIST=lcrcin)
          if (lcasca.EQV..true.)lcrcin=.false.
          if (lcrcin.EQV..false.)  open(lucoutxp, file=fnam(1:lenfnam)//'.crcin', form='unformatted', &
               status='unknown')
       else

          write(extension,'(i9.9)') itapp

          ! -------------------------------------------------------------
          ! -------------------------------------------------------------

          !            open(lucoutxp, file=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.cout', form='unformatted', &
          !                 status='unknown')


          !    ouverture d'un fichier .coutposition.(iteration) pour sauvegarde
          !    des positions toutes les itesauvposition iterations. Mode binaire
          !    pour l'ecriture.
          write (6, *)
          fnamcoutxp = fnam(1:lenfnam)//'.coutposition.'//extension




          write (6, *) 'sauvegarde positions it=', itapp, fnamcoutxp
          open(unit=lucoutxp, file=fnamcoutxp, form='unformatted', status='unknown')

       end if  !itapp  = 0 



       formatsauvT=2
       write (lucoutxp) formatsauvT
       write (lucoutxp) at
       write (lucoutxp) im_glob

#ifdef PARA
       im_loc(0)=im
       ibuffer=0
       ibuffer(1:im)  = ityp(1:im)
       buffer=0
       buffer(:,1:im) = xp(:,1:im)
       pt_im(0)=1
       next_pt = pt_im(0) + im_loc(0)
       do i_proc=1,nprocs-1
          call MPI_RECV(im_temp,1, MPI_INTEGER, MPI_ANY_SOURCE, 12001, MPI_COMM_space, status, ierr)
          proc_source = status(MPI_SOURCE)
          im_loc(proc_source)=im_temp
          pt_im(proc_source)=next_pt
          next_pt = pt_im(proc_source) + im_loc(proc_source)
          call MPI_RECV(ibuffer(pt_im(proc_source):pt_im(proc_source)+im_temp-1),    im_loc(proc_source),   &
               MPI_INTEGER,         proc_source, 12002, MPI_COMM_space, status, ierr)
          call MPI_RECV(buffer(1:3,pt_im(proc_source):pt_im(proc_source)+im_temp-1),3*im_loc(proc_source), &
               NDM_MPI_REAL_DOUBLE, proc_source, 12003, MPI_COMM_space, status, ierr)
       enddo

       !debug if (rang==0) write(*,*) 'Preparation ecriture ityp'
       write (lucoutxp) ibuffer  ! Ecriture ityp
       !debug     if (rang==0) write(*,*) 'Preparation ecriture ityp'
       write (lucoutxp) buffer   ! Ecriture xp
       !debug     if (rang==0) write(*,*) 'Ecriture xp'


       ibuffer(1:im) = num_at_glob(1:im)
       do i_proc=1,nprocs-1
          call MPI_RECV(ibuffer(pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),im_loc(i_proc), &
               MPI_INTEGER, i_proc, 12004, MPI_COMM_space, status, ierr)
       enddo
       write (lucoutxp) ibuffer   ! Ecriture num_at_glob

#else
       write (lucoutxp) ityp
       write (lucoutxp) xp
       write (lucoutxp) num_at_glob
#endif


    else ! rang.ne.0
       ! Pour etre en conformite avec la partie rang = 0 on sort si itapp n'est pas bon
       if (itapp>=99999999) then
          write (6, *) 'probleme de format dans sauveposition.F90'
          call arret_ndm
       endif

#ifdef PARA
       call MPI_SEND(im,          1,   MPI_INTEGER,        0,12001,MPI_COMM_space,ierr)
       call MPI_SEND(ityp(1:im),  im,  MPI_INTEGER,        0,12002,MPI_COMM_space,ierr)
       call MPI_SEND(xp(1:3,1:im),3*im,NDM_MPI_REAL_DOUBLE,0,12003,MPI_COMM_space,ierr)
       call MPI_SEND(num_at_glob(1:im), im,    MPI_INTEGER,0,12004,MPI_COMM_space,ierr)
#endif

    end if
    !100  format(i1)
200 format(i2)
300 format(i3)
400 format(i4)
500 format(i5)
600 format(i6)
700 format(i7)
800 format(i8)
1000 format(a5)
101 format(a1)
201 format(a2)
301 format(a3)
401 format(a4)
501 format(a5)
601 format(a6)
701 format(a7)
801 format(a8)

    if(rgloc==0)      close(lucoutxp)
    return
  end subroutine sauveposition
end module sauveposition_mod
