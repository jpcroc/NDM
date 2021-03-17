module sauvegardeT_mod

    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:rang,formatsauv,im_glob,it,itesauvinter,lspaceNDM,&
         &pmean,timel,tmean,tstep,fnam,lenfnam,lcasca,imm_glob,l2T

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
  subroutine sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp,lw0)
    !-----------------------------------------------
    !   M o d u l e s

    !latcomp= en PARA latcomp=.true.=> atmol est une cofiguration complète/latcomp=false=>atmol est distributé sur comm_space
    !lw0= .true. seul le proc 0 écrit la configuration

    implicit none
    type(box_config)::boxndm
    class(atom_config)::atdml
    type(cell_config):: celndm
    character::fnamcout*80
    logical, intent(in):: latcomp ! true= pas besoinde rapatrier atdml, false= il faut rapatrier atdml sur les masters
    logical, optional,intent(in):: lw0 ! seul le rang=0 écrit (implique latcomp=.true.)
    
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
    logical :: latcompin=.false.
    logical :: lw0in=.false.


#endif
    
    formatsauvmod = mod(formatsauv,2)
    im =atdml%im
    if (myidsp==0) then

       lucout = 87
       open(unit=lucout, file=fnamcout, form='unformatted', status='unknown')
       write (lucout) formatsauv
       write (lucout) boxndm%at
       write (lucout) im_glob
    end if
#ifdef PARA
    latcompin=latcomp

    if (present (lw0))lw0in=lw0
    if (lw0in) then
       if(latcompin.eqv..false.) then
          write(6,*)'comment sauvegarder seulement rang0 si latcomp=.false. ?'
          stop
       end if
       if (rang==0) then
          latcompin=.true.
       else
          latcompin=.false. !dans la suite latcompin intègre lw0 et rang=0 (NB on est dans ce if dans le cas lw0in=T)
       end if
    end if
       
    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.).and.(latcompin.eqv..false.)) then
       allocate (buffer(3,imm_glob))
       allocate (ibuffer(imm_glob))
    end if
    
    if (myidsp==0) then
       if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.).and.(latcompin.eqv..false.)) then
          im_loc(0)=im
          ibuffer=0
          ibuffer(1:im)  = atdml%ityp(1:im)
          buffer=0
          buffer(:,1:im) = atdml%xp(:,1:im)
          pt_im(0)=1
          next_pt = pt_im(0) + im_loc(0)

          do i_proc=1,nprocspace-1
             call comm_space%probe(11001,sourceout=proc_source)
             call comm_space%recv (im_temp,proc_source,11001)
             im_loc(proc_source)=im_temp
             pt_im(proc_source)=next_pt
             next_pt = pt_im(proc_source) + im_loc(proc_source)
             call comm_space%recv(ibuffer(pt_im(proc_source):pt_im(proc_source)+im_temp-1),proc_source,11002)
             call comm_space%recv(buffer(1:3,pt_im(proc_source):pt_im(proc_source)+im_temp-1),proc_source,11003)
          enddo
          write (lucout) ibuffer  ! Ecriture ityp
          write (lucout) buffer   ! Ecriture xp

          ibuffer(1:im) = atdml%num_at_glob(1:im)
          do i_proc=1,nprocspace-1
             call comm_space%recv(ibuffer(pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11004)
          enddo
          write (lucout) ibuffer   ! Ecriture num_at_glob

          if (formatsauvmod==1) then
             lwax=.false.
             select type(atdml)
             type is (atom_config_d)
                buffer(:,1:im) = atdml%xpp(:,1:im)
                do i_proc=1,nprocspace-1
                   call comm_space%recv(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11005)
                enddo
                write (lucout) buffer   ! Ecriture xpp

                buffer(:,1:im) = atdml%vp(:,1:im)
                do i_proc=1,nprocspace-1
                   call comm_space%recv(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11006)
                enddo
                write (lucout) buffer   ! Ecriture vp
             type is (atom_config_e)
                buffer(:,1:im) = atdml%xpp(:,1:im)
                do i_proc=1,nprocspace-1
                   call comm_space%recv(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11005)
                enddo
                write (lucout) buffer   ! Ecriture xpp

                buffer(:,1:im) = atdml%vp(:,1:im)
                do i_proc=1,nprocspace-1
                   call comm_space%recv(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11006)
                enddo
                write (lucout) buffer   ! Ecriture vp
                if (atdml%lax)then
                   lwax=.true.
                   buffer(:,1:im) = atdml%ax(:,1:im)
                   do i_proc=1,nprocspace-1
                      call comm_space%recv(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11007)
                   enddo
                   write (lucout) buffer   ! Ecriture ax
                end if
             end select
             if (.not.lwax)then
                buffer(:,1:im) = atdml%xp(:,1:im)
                do i_proc=1,nprocspace-1
                      call comm_space%recv(buffer(1:3,pt_im(i_proc):pt_im(i_proc)+im_loc(i_proc)-1),i_proc,11008)
                enddo
                write (lucout) buffer   ! Ecriture xpp

             end if
             write (lucout) tstep
             write (lucout) tmean, pmean, it, timel
          endif
       else
          if ((lw0in.eqv..false.).or.(lw0in.eqv..true.).and.(rang==0)) then
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
          endif
       end if

       close(unit=lucout)

       if (l2T)call sauveelec

    else ! myidsp different de 0 :

       if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.).and.(latcompin.eqv..false.)) then
          call comm_space%send(im,0,11001)
          call comm_space%send(atdml%ityp(1:im),0,11002)
          call comm_space%send(atdml%xp(1:3,1:im),0,11003)
          call comm_space%send(atdml%num_at_glob(1:im),0,11004)
          
          if (formatsauvmod==1) then
             lwax=.false.
             select type (atdml)
             type is (atom_config_d)
                call comm_space%send(atdml%xpp(1:3,1:im),0,11005)
                call comm_space%send(atdml%vp(1:3,1:im),0,11006)
             type is (atom_config_e)
                if (atdml%lax)then
                   lwax=.true.
                   call comm_space%send(atdml%ax(1:3,1:im),0,11007)
                end if
             end select
             if (.not.lwax)call comm_space%send(atdml%xp(1:3,1:im),0,11008)
          endif
       end if
    endif

    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.).and.(latcompin.eqv..false.)) then
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
    close(unit=lucout)
#endif
    if (l2T)call sauveelec

    return
  end subroutine sauvegardeT



    
end module sauvegardeT_mod
