module paraconfig
   USE arret_ndm_mod,only:arret_ndm

#ifdef PARA
  use mpi
  use Tpara,only:mpi_communicator
#else
  use Tpara,only:mpi_communicator
#endif
  use T_kind_param_m, ONLY:  double


  implicit none
#ifdef PARA

#endif
  type :: para_config
     type(mpi_communicator)::mpi_image,mpi_master,mpi_orig
     integer:: image ! n° de l'image
     integer::nimage ! = %mpi_master%nproc
     logical :: lmaster ! .true. si master
   contains
     procedure, pass::print


  end type para_config

#ifdef PARA

  
!!$  interface distribue
!!$     module procedure distribuereal, distribueinteger,distribuedouble
!!$  end interface distribue

#endif
contains
  subroutine commconstr(div)
    type(para_config)::div
#ifdef PARA
    !  integer:: grp_world,grp_masters,comm_masters,npm,rgm,imasters
    integer, allocatable::rgmasters(:)
    integer::ierr,ip

    integer::ipi,ipt,npi,npim,imasters,npm,reste,npr
    integer::clef, couleur,nimage,img,image
    integer,allocatable::npimg(:),GL(:),CL(:),ipimg(:,:)
    nimage=div%nimage
    if (div%nimage.gt.div%mpi_orig%nproc) then
       write(6,*)' division para impossible nimage > nprocs'
       call arret_ndm
       call MPI_FINALIZE(ierr)
    end if
    reste=mod(div%mpi_orig%nproc,div%nimage)
    allocate(rgmasters(div%nimage))
    allocate(ipimg(nimage,div%mpi_orig%nproc))
    allocate(npimg(nimage))
    allocate(GL(0:nimage-1));allocate(CL(0:nimage-1))
    if (reste==0) then
       npim=div%mpi_orig%nproc/div%nimage            ! nb de proc par image
       div%image=int(div%mpi_orig%rank/npim) !quel image pour le proc
       div%mpi_image%rank=mod(div%mpi_orig%rank,npim) ! quel rang dans le comm de l'image
       div%mpi_image%nproc=npim
       div%lmaster=.false.     
       if (div%mpi_image%rank==0) then
          div%lmaster=.true.
       end if
       ipi=0
       imasters=0
       npimg(:)=0
       ipimg(:,:)=-1
       do ip=0,div%mpi_orig%nproc-1
          image=1+int(ip/npim)
          npimg(image)=npimg(image)+1
          ipimg(image,npimg(image))=ip
          ipt=mod(ip,npim)
          if (ipt==0) then
             imasters=imasters+1
             rgmasters(imasters)=ip
          end if
       end do
       couleur=div%image
       clef=div%mpi_orig%rank

    else
     
       npim=div%mpi_orig%nproc/div%nimage
       npr=npim*(div%nimage-1)


       if (div%mpi_orig%rank.lt.npr) then
          div%mpi_image%nproc=npim
          div%mpi_image%rank=mod(div%mpi_orig%rank,npim)
          div%image=int(div%mpi_orig%rank/npim) !quel image pour le proc
       else
          div%image=div%nimage
          div%mpi_image%nproc=npim+reste
          div%mpi_image%rank=div%mpi_orig%rank-npr
       end if
       div%lmaster=.false.     
       if (div%mpi_image%rank==0) then
          div%lmaster=.true.
       end if
       ipi=0
       imasters=0

       do ip=0,div%mpi_orig%nproc-1
          if (ip.lt.npr) then
             image=1+int(ip/nimage)
          else
             image=nimage
          end if
          npimg(image)=npimg(image)+1
          ipimg(image,npimg(image))=ip
          ipt=mod(ip,npim)
          if (ipt==0) then
             imasters=imasters+1
             rgmasters(imasters)=ip
          end if
          
          ipt=mod(ip,npim)
          if (ipt==0) then
             imasters=imasters+1
             rgmasters(imasters)=ip
          end if
       end do
       couleur=div%image
       clef=div%mpi_orig%rank
    end if
Cl=0;GL=0
    do img=1,nimage
!       if (div%mpi_orig%rank==0) write(6,*)img,ipimg(img,1:npimg(img))
       call MPI_GROUP_INCL(div%mpi_orig%group,npimg(img),ipimg(img,1:npimg(img)),GL(img-1),ierr)
       call MPI_COMM_CREATE(div%mpi_orig%comm,GL(img-1),CL(img-1),ierr)
    end do
!    write(6,*)'TOTO', div%mpi_orig%rank,CL
!    write(6,*)'TATA', div%mpi_orig%rank,GL
!    write(6,*)'IMG',div%image,CL(div%image)
    call MPI_COMM_DUP(CL(div%image),div%mpi_image%comm,ierr)

    
!!$    call MPI_COMM_SPLIT (div%mpi_orig%comm,couleur,clef,div%mpi_image%comm,ierr)
!!$    call MPI_COMM_SIZE( div%mpi_image%comm, npi, ierr )
!!$    call MPI_COMM_RANK(div%mpi_image%comm, div%mpi_image%rank,ierr)
    if (div%mpi_orig%rank==0)then
!       write(6,*)'*************MPI DIVISION**************'
!       write(6,*)'orig_rank rank_in_image LMASTER Image'
    end if
!    write(6,*) div%mpi_orig%rank,div%mpi_image%rank,div%lmaster,div%image

    call MPI_BARRIER(div%mpi_orig%comm,ierr)
    
    call MPI_GROUP_INCL(div%mpi_orig%group,div%nimage,rgmasters,div%mpi_master%group,ierr)
    call MPI_COMM_CREATE(div%mpi_orig%comm,div%mpi_master%group,div%mpi_master%comm,ierr)

    if (div%lmaster) then
       call MPI_COMM_RANK(div%mpi_master%comm, div%mpi_master%rank,ierr)
       call MPI_COMM_SIZE( div%mpi_master%comm, npm, ierr )
       if (npm.ne.div%nimage) then
          write(6,*)'NPMPB',npm,div%nimage
          call arret_ndm
       end if
       div%mpi_master%nproc=npm
!      write(6,*)'Ranks among masters',div%mpi_orig%rank, div%mpi_master%rank
    end if
    call MPI_BARRIER(div%mpi_orig%comm,ierr)
    if (div%mpi_orig%rank==0)then
!       write(6,*)'*************MPI DIVISION**************'
    end if
    call MPI_BARRIER(div%mpi_orig%comm,ierr)
    div%nimage=div%nimage
#endif

!!$        call mpi_finalize(ierr)
!!$        call arret_ndm
    return

  end subroutine commconstr


  subroutine initparapuresp(div,rg,mpicsp)
    type(para_config)::div
    integer,intent(in)::rg!,nps
    type(mpi_communicator),intent(in)::mpicsp
!    integer,intent(in)::mpicsp
#ifdef PARA
    div%image=0
    div%nimage=1
    
    div%mpi_image%rank=rg
    div%mpi_image%comm=mpicsp%comm
    div%mpi_image%group=mpicsp%group
    div%mpi_image%nproc=mpicsp%nproc

    div%mpi_orig%comm=mpicsp%comm
    div%mpi_orig%rank=rg
    div%mpi_orig%nproc=mpicsp%nproc
    div%mpi_orig%group=mpicsp%group
    
    if (rg==0)then
       div%lmaster=.true.
       div%mpi_master%rank=0
       div%mpi_master%comm=-1
       div%mpi_master%nproc=1
    else
       div%lmaster=.false.
    end if
#else
    div%image=0
    div%nimage=1
    div%mpi_image%rank=0
    if (rg==0)then
       div%lmaster=.true.
       div%mpi_master%rank=0
       div%mpi_master%comm=-1
       div%mpi_master%nproc=1
    else
       div%lmaster=.false.
    end if
    
#endif
  end subroutine initparapuresp
!!$#ifdef PARA
!!$  subroutine distribuereal(div,x,nel,xrecv)
!!$    integer::nel
!!$    type(para_config)::div
!!$    real::x(0:div%nimage-1,nel)
!!$    real,allocatable::xrecv(:)
!!$    integer::img,ierr,isp
!!$    integer, dimension( MPI_STATUS_SIZE) :: statut
!!$
!!$    if (.not.(allocated(xrecv)))allocate(xrecv(1:nel))
!!$    !     write(6,*)div
!!$    !    if (div%lmaster) then
!!$    !       write(6,*)'PRE',div%lmaster,div%mpi_orig%rank,div%mpi_image%rank,div%mpi_master%rank
!!$    !    else
!!$    !       write(6,*)'PRE',div%lmaster,div%mpi_orig%rank,div%mpi_image%rank
!!$    !    end if
!!$    if (div%mpi_orig%rank==0) then
!!$       xrecv(:)=x(0,:)
!!$       write(6,*)'rg0 ',xrecv
!!$       do img=1,div%nimage-1
!!$          !           write(6,*)'send 0->',img,x(img,:)
!!$          call MPI_SEND (x(img,1:nel),nel, MPI_REAL ,img,100, div%mpi_master%comm,ierr)
!!$       end do
!!$       do isp=1,div%mpi_image%nproc-1
!!$          call MPI_SEND(x(0,1:nel),nel,MPI_REAL,isp,102,div%mpi_image%comm,ierr)
!!$       end do
!!$    else
!!$       if (div%lmaster) then
!!$          call MPI_RECV(xrecv(1:nel),nel,MPI_REAL,0,100,div%mpi_master%comm,statut,ierr)
!!$          !           write(6,*)'recv->',xrecv(:),div%mpi_master%rank
!!$          do isp=1,div%mpi_image%nproc-1
!!$             call MPI_SEND(xrecv(1:nel),nel,MPI_REAL,isp,101,div%mpi_image%comm,ierr)
!!$          end do
!!$       else
!!$          if (div%image==0) then 
!!$             call MPI_RECV(xrecv(1:nel),nel,MPI_REAL,0,102,div%mpi_image%comm,statut,ierr)
!!$          else
!!$             call MPI_RECV(xrecv(1:nel),nel,MPI_REAL,0,101,div%mpi_image%comm,statut,ierr)
!!$          end if
!!$       end if
!!$    end if
!!$
!!$  end subroutine distribuereal
!!$
!!$  subroutine distribuedouble(div,x,nel,xrecv)
!!$    integer::nel
!!$    type(para_config)::div
!!$    real::x(0:div%nimage-1,nel)
!!$    real(double),allocatable::xrecv(:)
!!$    integer::img,ierr,isp
!!$    integer, dimension( MPI_STATUS_SIZE) :: statut
!!$
!!$    if (.not.(allocated(xrecv)))allocate(xrecv(1:nel))
!!$    !     write(6,*)div
!!$    !    if (div%lmaster) then
!!$    !       write(6,*)'PRE',div%lmaster,div%mpi_orig%rank,div%mpi_image%rank,div%mpi_master%rank
!!$    !    else
!!$    !       write(6,*)'PRE',div%lmaster,div%mpi_orig%rank,div%mpi_image%rank
!!$    !    end if
!!$    if (div%mpi_orig%rank==0) then
!!$       xrecv(:)=x(0,:)
!!$       write(6,*)'rg0 ',xrecv
!!$       do img=1,div%nimage-1
!!$          !           write(6,*)'send 0->',img,x(img,:)
!!$          call MPI_SEND (x(img,1:nel),nel, NDM_MPI_REAL_DOUBLE ,img,100, div%mpi_master%comm,ierr)
!!$       end do
!!$       do isp=1,div%mpi_image%nproc-1
!!$          call MPI_SEND(x(0,1:nel),nel,NDM_MPI_REAL_DOUBLE,isp,102,div%mpi_image%comm,ierr)
!!$       end do
!!$    else
!!$       if (div%lmaster) then
!!$          call MPI_RECV(xrecv(1:nel),nel,NDM_MPI_REAL_DOUBLE,0,100,div%mpi_master%comm,statut,ierr)
!!$          !           write(6,*)'recv->',xrecv(:),div%mpi_master%rank
!!$          do isp=1,div%mpi_image%nproc-1
!!$             call MPI_SEND(xrecv(1:nel),nel,NDM_MPI_REAL_DOUBLE,isp,101,div%mpi_image%comm,ierr)
!!$          end do
!!$       else
!!$          if (div%image==0) then 
!!$             call MPI_RECV(xrecv(1:nel),nel,NDM_MPI_REAL_DOUBLE,0,102,div%mpi_image%comm,statut,ierr)
!!$          else
!!$             call MPI_RECV(xrecv(1:nel),nel,NDM_MPI_REAL_DOUBLE,0,101,div%mpi_image%comm,statut,ierr)
!!$          end if
!!$       end if
!!$    end if
!!$  end subroutine distribuedouble
!!$  subroutine distribueinteger(div,x,nel,xrecv)
!!$    integer::nel
!!$    type(para_config)::div
!!$    integer::x(0:div%nimage-1,nel)
!!$    integer,allocatable::xrecv(:)
!!$    integer::img,ierr,isp
!!$    integer, dimension( MPI_STATUS_SIZE) :: statut
!!$
!!$    if (.not.(allocated(xrecv)))allocate(xrecv(1:nel))
!!$    !     write(6,*)div
!!$    !    if (div%lmaster) then
!!$    !       write(6,*)'PRE',div%lmaster,div%mpi_orig%rank,div%mpi_image%rank,div%mpi_master%rank
!!$    !    else
!!$    !       write(6,*)'PRE',div%lmaster,div%mpi_orig%rank,div%mpi_image%rank
!!$    !    end if
!!$    if (div%mpi_orig%rank==0) then
!!$       xrecv(:)=x(0,:)
!!$       write(6,*)'rg0 ',xrecv
!!$       do img=1,div%nimage-1
!!$          !           write(6,*)'send 0->',img,x(img,:)
!!$          call MPI_SEND (x(img,1:nel),nel, MPI_INTEGER ,img,200, div%mpi_master%comm,ierr)
!!$       end do
!!$       do isp=1,div%mpi_image%nproc-1
!!$          call MPI_SEND(x(0,1:nel),nel,MPI_INTEGER,isp,202,div%mpi_image%comm,ierr)
!!$       end do
!!$    else
!!$       if (div%lmaster) then
!!$          call MPI_RECV(xrecv(1:nel),nel,MPI_INTEGER,0,200,div%mpi_master%comm,statut,ierr)
!!$          !           write(6,*)'recv->',xrecv(:),div%mpi_master%rank
!!$          do isp=1,div%mpi_image%nproc-1
!!$             call MPI_SEND(xrecv(1:nel),nel,MPI_INTEGER,isp,201,div%mpi_image%comm,ierr)
!!$          end do
!!$       else
!!$          if (div%image==0) then 
!!$             call MPI_RECV(xrecv(1:nel),nel,MPI_INTEGER,0,202,div%mpi_image%comm,statut,ierr)
!!$          else
!!$             call MPI_RECV(xrecv(1:nel),nel,MPI_INTEGER,0,201,div%mpi_image%comm,statut,ierr)
!!$          end if
!!$       end if
!!$    end if
!!$  end subroutine distribueinteger
!!$
!!$
!!$#endif

  subroutine print(paraprt,rang)
    class(para_config),intent(in)::paraprt
    integer,intent(in)::rang
!    write(6,*)'in print paraprt'
    write(rang+100,*)'rang rank comm group nproc'
    write(rang+100,*)'ORIG',rang,paraprt%mpi_orig%rank,paraprt%mpi_orig%comm,paraprt%mpi_orig%group,paraprt%mpi_orig%nproc
    write(rang+100,*)'NIMAGE',rang,paraprt%nimage
    write(rang+100,*)'IMAGEnum',rang,paraprt%image
    write(rang+100,*)'LMASTER',rang,paraprt%lmaster
    write(rang+100,*)'IMAGEcom',rang,paraprt%mpi_image%rank,paraprt%mpi_image%comm,paraprt%mpi_image%group,paraprt%mpi_image%nproc
    if (paraprt%lmaster)write(rang+100,*)'MASTER',rang,paraprt%mpi_master%rank,paraprt%mpi_master%comm,&
         &paraprt%mpi_master%group,paraprt%mpi_master%nproc
    flush(rang+100)
  end subroutine print
    


    
end module paraconfig



