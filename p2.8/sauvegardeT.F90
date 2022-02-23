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
  ! ********************************************************************
  subroutine sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp)
    !-----------------------------------------------
    !   M o d u l e s

    !latcomp= en PARA latcomp=.true.=> atmol est une cofiguration complète/latcomp=false=>atmol est distributé sur comm_space

    implicit none
    type(box_config)::boxndm
    class(atom_config)::atdml
    type(cell_config):: celndm
    character::fnamcout*80
    logical, intent(in):: latcomp ! true= pas besoinde rapatrier atdml, false= il faut rapatrier atdml sur les masters

    integer :: lucout, formatsauvmod,i,formatsauv,im
    character :: extension*9
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
    if (rang==0) write(6,*)'insauvegarde',iteration
    formatsauvmod = mod(formatsauv,2)
    im =atdml%im
    if (atdml%im_glob==0) then
       write(6,*)'sauvegarde imglob=0 stop'
       call arret_ndm
    end if

    lucout = 87
#ifdef PARA
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
                call comm_space%probe(11001,sourceout=proc_source)
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
             type is (atom_config_d)
                buffer(:,1:im) = atdml%xpp(:,1:im)
                if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
                   do i_proc=1,nprocspace-1
                      call comm_space%recv(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11015)
                   enddo
                end if
                write (lucout) buffer   ! Ecriture xpp
                buffer(:,1:im) = atdml%vp(:,1:im)
                if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
                   do i_proc=1,nprocspace-1
                      call comm_space%recv(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11016)
                   enddo
                end if
                write (lucout) buffer   ! Ecriture vp
             type is (atom_config_e)
                buffer(:,1:im) = atdml%xpp(:,1:im)
                if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
                   do i_proc=1,nprocspace-1
                      call comm_space%recv(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11005)
                   enddo
                end if
                write (lucout) buffer   ! Ecriture xpp

                buffer(:,1:im) = atdml%vp(:,1:im)
                if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
                   do i_proc=1,nprocspace-1
                      call comm_space%recv(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11006)
                   enddo
                end if
                write (lucout) buffer   ! Ecriture vp
                if (atdml%lax)then
                   lwax=.true.
                   buffer(:,1:im) = atdml%ax(:,1:im)
                   if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
                      do i_proc=1,nprocspace-1
                         call comm_space%recv(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11007)
                      enddo
                   end if
                   write (lucout) buffer   ! Ecriture ax
                end if
             end select
             if (.not.lwax)then
                buffer(:,1:im) = atdml%xp(:,1:im)
                if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
                   do i_proc=1,nprocspace-1
                      call comm_space%recv(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11008)
                   enddo
                end if
                write (lucout) buffer   ! Ecriture xpp

             end if
          end if
          write (lucout) tstep
          write (lucout) tmean, pmean, iteration, timel



          if (l2T)call sauveelec

       else ! myidsp different de 0 :

          if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
             call comm_space%send(im,0,11001)
             call comm_space%send(atdml%ityp(1:im),0,11002)
             call comm_space%send(atdml%xp(1:3,1:im),0,11003)
             call comm_space%send(atdml%num_at_glob(1:im),0,11004)

             if (formatsauvmod==1) then
                lwax=.false.
                select type (atdml)
                type is (atom_config_d)
                   call comm_space%send(atdml%xpp(1:3,1:im),0,11015)
                   call comm_space%send(atdml%vp(1:3,1:im),0,11016)
                type is (atom_config_e)
                   call comm_space%send(atdml%xpp(1:3,1:im),0,11005)
                   call comm_space%send(atdml%vp(1:3,1:im),0,11006)
                   if (atdml%lax)then
                      lwax=.true.
                      call comm_space%send(atdml%ax(1:3,1:im),0,11007)
                   end if
                end select
                if (.not.lwax)call comm_space%send(atdml%xp(1:3,1:im),0,11008)
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
       write (lucout) tmean, pmean, iteration, timel
    endif

#endif
    lucout=87
    close(unit=lucout)
    if (l2T)call sauveelec

    return
  end subroutine sauvegardeT




end module sauvegardeT_mod
