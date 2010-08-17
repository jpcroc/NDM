module suivinonpbc
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------

  !         version du 09 juin 2010, 11h23 AM, last change by MCM
  ! ********************************************************************

  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
#if(PARA)
  use mod_mpi
#endif

  
  



CONTAINS
  
subroutine init_suivinonpbc
implicit none
integer ::i
!
 DO i=1,imd
  tmpsuivi(1:3,i) = zero
  xpnonpbc(1:3,i)=axnonpbc(1:3,i)
 END DO
!
end subroutine  init_suivinonpbc

subroutine reset_suivinonpbc
implicit none
integer :: i 
!
 DO i=1,imd  
 xpnonpbc(1:3,i)= xpnonpbc(1:3,i) + tmpsuivi(1:3,i)
 tmpsuivi(1:3,i) = zero 
 END DO
!
end subroutine 
 
subroutine sauvepositionnonpbc(itapp)

  implicit none
  integer :: itapp ! iteration apparente

  integer ::  formatsauvT,lenfn2,lucoutnonpbcxp
  character :: extension*9

#if(PARA)
  integer,dimension(:),pointer     :: ibuffer
  real(double), dimension(:,:),pointer   :: buffer
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
  !

  ! conversion entier-->alphanumerique par transfert du nombre
  ! de l'iteration vers fichier tampon relu sous format caractere.

  lucoutnonpbcxp = 111
  if(rang==0) then
     open(unit=17, file='tampon', form='formatted', status='unknown')
        if (itapp<=9) then
           write (17, '(I1)') itapp
           rewind 17
           read (17, 101) extension
           write(6,*)extension
        end if
        if (itapp<=99.and.itapp>9) then
           write (17, 200) itapp
           rewind 17
           read (17, 201) extension
           write(6,*)extension
        end if
        if (itapp<=999.and.itapp>99) then
           write (17, 300) itapp
           rewind 17
           read (17, 301) extension
        end if
        if (itapp<=9999.and.itapp>999) then
           write (17, 400) itapp
           rewind 17
           read (17, 401) extension
        end if
        if (itapp<=99999.and.itapp>9999) then
           write (17, 500) itapp
           rewind 17
           read (17, 501) extension
        end if
        if (itapp<=999999.and.itapp>99999)  then
           write (17, 600) itapp
           rewind 17
           read (17, 601) extension
        end if
        if (itapp<=9999999.and.itapp>999999)  then
           write (17, 700) itapp
           rewind 17
           read (17, 701) extension
        end if
        if (itapp<=99999999.and.itapp>9999999)  then
           write (17, 800) itapp
           rewind 17
           read (17, 801) extension
        end if
        if (itapp>=99999999) then
           write (6, *) 'probleme de format dans sauveposition'
           call arret_ndm
        endif
        lenfn2=index(extension,' ')-1
        ! -------------------------------------------------------------
        ! -------------------------------------------------------------

        !            open(lucoutxp, file=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.cout', form='unformatted', &
        !                 status='unknown')


        !    ouverture d'un fichier .coutposition.(iteration) pour sauvegarde
        !    des positions toutes les itesauvposition iterations. Mode binaire
          fnamcoutnonpbcxp = fnam(1:lenfnam)//'.coutpositiononpbc.'//extension
	
	
        write (6, *) 'sauvegarde < non pbc >  positions it=', itapp, fnamcoutxp
        open(unit=lucoutnonpbcxp, file=fnamcoutnonpbcxp, form='unformatted', status='unknown')


	 formatsauvT=2
	 write (lucoutnonpbcxp) formatsauvT
         write (lucoutnonpbcxp) at
         write (lucoutnonpbcxp) im_glob

#if(PARA)
     im_loc(0)=im
     ibuffer=0
     ibuffer(1:im)  = ityp(1:im)
     buffer=0
     buffer(:,1:im) = xpnonpbc(:,1:im)
     pt_im(0)=1
     next_pt = pt_im(0) + im_loc(0)
     do i_proc=1,nprocs-1
        call MPI_RECV(im_temp,1, MPI_INTEGER, MPI_ANY_SOURCE, 12021, MPI_COMM_WORLD, status, ierr)
        proc_source = status(MPI_SOURCE)
        im_loc(proc_source)=im_temp
        pt_im(proc_source)=next_pt
        next_pt = pt_im(proc_source) + im_loc(proc_source)
        call MPI_RECV(ibuffer(pt_im(proc_source):pt_im(proc_source)+im_temp-1),    im_loc(proc_source),   &
             MPI_INTEGER,         proc_source, 12022, MPI_COMM_WORLD, status, ierr)
        call MPI_RECV(buffer(1:3,pt_im(proc_source):pt_im(proc_source)+im_temp-1),3*im_loc(proc_source), &
             NDM_MPI_REAL_DOUBLE, proc_source, 12023, MPI_COMM_WORLD, status, ierr)
     enddo
      write (lucoutnonpbcxp) ibuffer  ! Ecriture ityp
      write (lucoutnonpbcxp) buffer   ! Ecriture xpnonpbc

     ibuffer(1:im) = num_at_glob(1:im)
     do i_proc=1,nprocs-1
        call MPI_RECV(ibuffer(pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),im_loc(i_proc), &
             MPI_INTEGER, i_proc, 12024, MPI_COMM_WORLD, status, ierr)
     enddo
      
      write (lucoutnonpbcxp) ibuffer  ! Ecriture ityp

#else
     write (lucoutnonpbcxp) ityp
     write (lucoutnonpbcxp) xpnonpbc
     write (lucoutnonpbcxp) num_at_glob
#endif


  else ! rang.ne.0
     ! Pour etre en conformite avec la partie rang = 0 on sort si itapp n'est pas bon
     if (itapp>=99999999) then
        write (6, *) 'probleme de format dans sauvepositionnonpbc.F90'
        call arret_ndm
    endif

#if(PARA)
     call MPI_SEND(im,          1,   MPI_INTEGER,              0,12021,MPI_COMM_WORLD,ierr)
     call MPI_SEND(ityp(1:im),  im,  MPI_INTEGER,              0,12022,MPI_COMM_WORLD,ierr)
     call MPI_SEND(xpnonpbc(1:3,1:im),3*im,NDM_MPI_REAL_DOUBLE,0,12023,MPI_COMM_WORLD,ierr)
     call MPI_SEND(num_at_glob(1:im), im,    MPI_INTEGER,      0,12024,MPI_COMM_WORLD,ierr)
     deallocate (buffer)
     deallocate(ibuffer)
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

  if(rang==0)      close(17)
  if(rang==0)      close(lucoutnonpbcxp)
  return
end subroutine sauvepositionnonpbc

end module suivinonpbc
