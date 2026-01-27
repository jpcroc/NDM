module sauvegardeT_mod
   USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:rang,iteration,itesauvinter,lspaceNDM,&
       &pmean,timel,tmean,tstep,fnam,lenfnam,lcasca,l2T

  USE elec_cell, ONLY : sauveelec
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE cellconfig, only:cell_config
  USE boxconfig,only:box_config
#ifdef PARA
  USE Tpara,only:COMM_space,nprocspace,myidsp
#else
  USE Tpara,only:nprocspace    ,myidsp

#endif



  implicit none

contains

  subroutine sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp)
    implicit none
    class(box_config)::boxndm
    class(atom_config)::atdml
    type(cell_config):: celndm
    character::fnamcout*80
    logical, intent(in):: latcomp ! true= pas besoinde rapatrier atdml, false= il faut rapatrier atdml sur les masters
    integer :: formatsauv

    call sauvegardeT_para(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp)
    !call sauvegardeT_originale(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp)

  end subroutine sauvegardeT


  ! ********************************************************************
  subroutine sauvegardeT_originale(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp)
    !-----------------------------------------------
    !   M o d u l e s

    !latcomp= en PARA latcomp=.true.=> atmol est une cofiguration complète/latcomp=false=>atmol est distributé sur comm_space

    implicit none
    class(box_config)::boxndm
    class(atom_config)::atdml
    type(cell_config):: celndm
    character::fnamcout*80
    logical, intent(in):: latcomp ! true= pas besoinde rapatrier atdml, false= il faut rapatrier atdml sur les masters

    integer :: lucout, formatsauvmod,formatsauv,im!,formatsauvw
    logical :: lwax

#ifdef PARA
    integer,dimension(:),allocatable     :: ibuffer
    real(double), dimension(:,:),allocatable   :: buffer
    integer,      dimension(0:nprocspace-1)   :: im_loc
    integer,      dimension(0:nprocspace-1)   :: pt_im
    integer :: next_pt
    integer :: i_proc
    integer :: proc_source
    integer :: im_temp


#endif
!    if (rang==0) write(6,*)'insauvegarde',iteration,timel
    formatsauvmod = mod(formatsauv,2)
    im =atdml%im
    if (atdml%im_glob==0) then
       write(6,*)'sauvegarde imglob=0 stop'
       call arret_ndm
    end if

    lucout = 87
!!$    select type(atdml)
!!$       class is (atom_config_e)
!!$          if (atdml%lax)then
!!$             if (formatsauv>=5) formatsauvw=7
!!$          else
!!$             formatsauvw=formatsauv
!!$          end if
!!$       class default
!!$          formatsauvw=formatsauv
!!$       end select
#ifdef PARA
!       write(6,*)'SPDBG1 ',rang,formatsauv
    if (.not.latcomp) then 
       if (myidsp==0) then


          open(unit=lucout, file=fnamcout, form='unformatted', status='unknown')
          write (lucout) formatsauv
          write (lucout) boxndm%at
          write (lucout) atdml%im_glob
       end if
       if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
          allocate (buffer(3,atdml%imm_glob))
          allocate (ibuffer(atdml%imm_glob))
       else
          allocate (buffer(3,atdml%imm))
          allocate (ibuffer(atdml%imm))
       end if
!       write(6,*)'SPDBG2 ',rang,formatsauv,formatsauvmod
       if (myidsp==0) then
          im_loc(0)=im
          ibuffer=0
          ibuffer(1:im)  = atdml%ityp(1:im)
          buffer=0
          buffer(:,1:im) = atdml%xp(:,1:im)
          pt_im(0)=1
          next_pt = pt_im(0) + im_loc(0)
          if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
             do i_proc=1,nprocspace-1
                !call comm_space%probe(11001,sourceout=proc_source)
                proc_source=i_proc
                call comm_space%recv (im_temp,proc_source,11001)
                im_loc(proc_source)=im_temp
                pt_im(proc_source)=next_pt
                next_pt = pt_im(proc_source) + im_loc(proc_source)
                call comm_space%recv(ibuffer(pt_im(proc_source):pt_im(proc_source)+im_temp-1),proc_source,11002)
                call comm_space%recv(buffer(1:3,pt_im(proc_source):pt_im(proc_source)+im_temp-1),proc_source,11003)
             enddo
          end if
          write (lucout) ibuffer  ! Ecriture ityp
          write (lucout) buffer   ! Ecriture xp
          ibuffer(1:im) = atdml%num_at_glob(1:im)
          if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
             do i_proc=1,nprocspace-1
                call comm_space%recv(ibuffer(pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11004)
             enddo
          end if
          write (lucout) ibuffer   ! Ecriture num_at_glob
          if (formatsauvmod==1) then
             lwax=.false.
             select type(atdml)
             class is (atom_config_d)
                buffer(:,1:im) = atdml%vp(:,1:im)
                if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
                   do i_proc=1,nprocspace-1
                      call comm_space%recv(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11016)
                   enddo
                end if
                write (lucout) buffer   ! Ecriture vp
             class is (atom_config_e)
                buffer(:,1:im) = atdml%vp(:,1:im)
                if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
                   do i_proc=1,nprocspace-1
                      call comm_space%recv(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11006)
                   enddo
                end if
                write (lucout) buffer   ! Ecriture vp

             end select
          end if
          write (lucout) tstep
          write (lucout) tmean, pmean, iteration, timel
          if (l2T)call sauveelec
!       write(6,*)'SPDBG3 ',rang
    else ! myidsp different de 0 :
          if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
             call comm_space%send(im,0,11001)
             call comm_space%send(atdml%ityp(1:im),0,11002)
             call comm_space%send(atdml%xp(1:3,1:im),0,11003)
             call comm_space%send(atdml%num_at_glob(1:im),0,11004)
             if (formatsauvmod==1) then
                lwax=.false.
                select type (atdml)
                class is (atom_config_d)
                   call comm_space%send(atdml%vp(1:3,1:im),0,11016)
                class is (atom_config_e)
                   call comm_space%send(atdml%vp(1:3,1:im),0,11006)
!                   if ((atdml%lxpp).and.(formatsauv==3)) call comm_space%send(atdml%xpp(1:3,1:im),0,11015)
!!$                   if (formatsauvw==7)then
!!$                      lwax=.true.
!!$                      call comm_space%send(atdml%ax(1:3,1:im),0,11007)
!!$                   end if
                end select
!!$                write(6,*)'SPDBG251 ',rang,lwax
!!$                if (.not.lwax)call comm_space%send(atdml%xp(1:3,1:im),0,11008)
!                write(6,*)'SPDBG25 ',rang
             endif
          end if
       endif
    else !latcomp
       lucout = 87
       open(unit=lucout, file=fnamcout, form='unformatted', status='unknown')
       write (lucout) formatsauv
       write (lucout) boxndm%at
       write (lucout) atdml%im
       write (lucout) atdml%ityp
       write (lucout) atdml%xp
       write (lucout) atdml%num_at_glob
       if (formatsauvmod==1) then
          lwax=.false.
          select type (atdml)
          class is (atom_config_d)
             write (lucout) atdml%vp
          end select
          select type (atdml)
          class is (atom_config_e)

!             write (lucout) atdml%vp
!             if (atdml%lxpp)             write (lucout) atdml%xpp
!!$             if (formatsauvw==7)then
!!$                write (lucout) atdml%ax
!!$                lwax=.true.
!!$             end if
          end select
          if (.not.lwax)write (lucout) atdml%xp
          write (lucout) tstep
          write (lucout) tmean, pmean, iteration, timel
       endif

    endif

    if (allocated(buffer)) then
       deallocate (buffer)
       deallocate (ibuffer)
    end if
#else
    ! sauvegarde SEQ
    !
    ! Partie sequentielle de la sauvegarde :
    !
    !     do i=1,im
    !        write(1004,*)i,num_at_glob(i),xp(1,i)
    !     enddo
    lucout = 87
    open(unit=lucout, file=fnamcout, form='unformatted', status='unknown')
    write (lucout) formatsauv
    write (lucout) boxndm%at
    write (lucout) atdml%im
    write (lucout) atdml%ityp
    write (lucout) atdml%xp
    write (lucout) atdml%num_at_glob
    if (formatsauvmod==1) then
       lwax=.false.
       select type (atdml)
       class is (atom_config_d)
          write (lucout) atdml%vp
       end select
       select type (atdml)
       class is (atom_config_e)
!          write (lucout) atdml%xpp
!          write (lucout) atdml%vp
!!$          if (atdml%lxpp)          write (lucout) atdml%xpp
!!$          if (formatsauvw==7)then
!!$             write (lucout) atdml%ax
!!$             lwax=.true.
!!$          end if
       end select
!       if (.not.lwax)write (lucout) atdml%xp
       write (lucout) tstep
       write (lucout) tmean, pmean, iteration, timel
    endif

#endif
    lucout=87
    close(unit=lucout)
    if (l2T)call sauveelec
    return
  end subroutine sauvegardeT_originale




  ! ********************************************************************
  subroutine sauvegardeT_para(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp)
    !-----------------------------------------------
    !   Version parallèle (MPI-IO) de sauvegardeT
    !-----------------------------------------------
    
    ! formatsauv=2 => positions
    ! formatsauv=3 => positions + vitesses + xpp      + tstep, tmean, pmean, iteration, timel
    ! formatsauv=4 => positions 
    ! formatsauv=5 => positions + vitesses            + tstep, tmean, pmean, iteration, timel


#ifdef PARA
    use Tpara_io
    use Tpara, only: NDM_MPI_REAL_DOUBLE
#endif
    

    implicit none
    class(box_config)::boxndm
    class(atom_config)::atdml
    type(cell_config):: celndm
    character::fnamcout*80
    logical, intent(in):: latcomp ! true= pas besoinde rapatrier atdml, false= il faut rapatrier atdml sur les masters

    integer :: lucout, formatsauvmod,formatsauv,im
    integer :: mpi_size_double, mpi_size_int, i
    logical :: lwax
    integer, dimension(:), allocatable :: all_im

#ifdef PARA
    integer(KIND=MPI_OFFSET_KIND) :: offset, para_offset
#endif

    formatsauvmod = mod(formatsauv,2)
    im =atdml%im
    if (atdml%im_glob==0) then
       write(6,*)'sauvegarde imglob=0 stop'
       call arret_ndm
    end if

#ifdef PARA
    ! ******************* sauvegarde PARA *********************

    mpi_size_double = type_size(NDM_MPI_REAL_DOUBLE) ! mpi_size_double = double sinon erreurs
    mpi_size_int = type_size(MPI_INTEGER)

    
    ! Calcul de para_offset
    allocate(all_im(nprocspace))
    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
      ! Envoie et réception du nombre d'atomes de chaques procs
      call mpic_allgather_i(comm_space, atdml%im, all_im)
    else
      all_im=atdml%im ! un seul proc ou lspaceNDM .false. ?
    end if

    para_offset=0
    do i=1, myidsp
      para_offset = para_offset + all_im(i)
    end do
    deallocate(all_im)


    call mpic_file_open(comm_space, fnamcout ,lucout)
    offset = 0


    ! Entête
    if (myidsp==0) then
      call file_write_at(lucout, offset, formatsauv)      ! Ecriture formatsauv
      offset = offset + mpi_size_int
      call file_write_at(lucout, offset, boxndm%at)       ! Ecriture boxndm%at
      offset = offset + mpi_size_double*size(boxndm%at)
      call file_write_at(lucout, offset, atdml%im_glob)   ! Ecriture im_glob
      offset = offset + mpi_size_int
    end if


    ! Corps
    call file_write_at_all(lucout, offset + para_offset*3*mpi_size_double, atdml%xp(1:3,1:atdml%im))   ! Ecriture xp
    offset = offset + mpi_size_double*3*atdml%im_glob
    call file_write_at_all(lucout, offset + para_offset*mpi_size_int, atdml%ityp(1:atdml%im))   ! Ecriture ityp
    offset = offset + mpi_size_int*atdml%im_glob
    call file_write_at_all(lucout, offset + para_offset*mpi_size_int, atdml%num_at_glob(1:atdml%im))   ! Ecriture num_at_glob
    offset = offset + mpi_size_int*atdml%im_glob


    if (formatsauvmod==1) then
      if (formatsauv==3) then
         select type (atdml)
         class is (atom_config_e)
            call file_write_at_all(lucout, offset + para_offset*3*mpi_size_double, atdml%xpp(1:3,1:atdml%im))   ! Ecriture xpp
            offset = offset + mpi_size_double*3*atdml%im_glob
         class default
            write(6,*) "sauvegarde demandée avec xpp, mais atom_config n'a pas xpp, stop"
         call arret_ndm
         end select
      end if
      select type (atdml)
      class is (atom_config_d) ! atom_config_e extends atom_config_d, donc on entre ici aussi avec atom_config_e
        call file_write_at_all(lucout, offset + para_offset*3*mpi_size_double, atdml%vp(1:3,1:atdml%im))   ! Ecriture vp
        offset = offset + mpi_size_double*3*atdml%im_glob
      class default
        write(6,*) "sauvegarde demandée avec vp, mais atom_config n'a pas vp, stop"
        call arret_ndm
      end select

      if (myidsp==0) then
        call file_write_at(lucout, offset, tstep)      ! Ecriture tstep
        offset = offset + mpi_size_double
        call file_write_at(lucout, offset, tmean)      ! Ecriture tmean
        offset = offset + mpi_size_double
        call file_write_at(lucout, offset, pmean)      ! Ecriture pmean
        offset = offset + mpi_size_double
        call file_write_at(lucout, offset, iteration)      ! Ecriture iteration
        offset = offset + mpi_size_int
        call file_write_at(lucout, offset, timel)      ! Ecriture timel
        offset = offset + mpi_size_double

      end if
    endif

    call file_close(lucout)
    if (myidsp==0) then
      if (l2T) call sauveelec ! A faire absolument sur le proc rang=0
    end if

#else

    ! ******************* sauvegarde SEQ *********************
    lucout = 87
    open(unit=lucout, file=fnamcout, form='unformatted', access='stream', status='unknown')
    ! Entête
    write (lucout) formatsauv
    write (lucout) boxndm%at
    write (lucout) atdml%im
    ! Corps
    write (lucout) atdml%xp(1:3,1:atdml%im)
    write (lucout) atdml%ityp(1:atdml%im)
    write (lucout) atdml%num_at_glob(1:atdml%im)
    if (formatsauvmod==1) then
      if (formatsauvmod==3) then
         select type (atdml)
         class is (atom_config_e)
            write (lucout) atdml%xpp(1:3,1:atdml%im)
         class default
            write(6,*) "sauvegarde demandée avec xpp, mais atom_config n'a pas xpp, stop"
         call arret_ndm
         end select
      end if
      select type (atdml)
      class is (atom_config_d) ! atom_config_e extends atom_config_d, donc on entre ici aussi avec atom_config_e
        write (lucout) atdml%vp(1:3,1:atdml%im)
      class default
        write(6,*) "sauvegarde demandée avec vp, mais atom_config n'a pas vp, stop"
        call arret_ndm
      end select
      write (lucout) tstep                            ! Potentiellement à l'extérieur du if (formatsauvmod==1)
      write (lucout) tmean, pmean, iteration, timel   ! Potentiellement à l'extérieur du if (formatsauvmod==1)
    endif
    close(unit=lucout)
    if (l2T)call sauveelec ! A faire absolument sur le proc rang=0
#endif

  end subroutine sauvegardeT_para

end module sauvegardeT_mod
