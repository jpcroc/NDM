
module paraconfig


#ifdef PARA
  use Tpara,only:NDM_MPI_REAL_DOUBLE

#endif
  use T_kind_param_m, ONLY:  double


  implicit none
#ifdef PARA
  include 'mpif.h'
#endif
  type para_config
     integer :: rgim ! rang du proc dans l'image
     integer:: image ! n° de l'image
     integer:: npim ! nb de procs de l'image
     integer:: grp_master !
!     integer:: npmas ! nb de procs masters=nb d'image
     integer::nimage
     integer:: rgmas ! rang parmi les masters
     integer:: comm_image ! communicateur associé à l'image (peuvent être égaux pour toutes les imags, ça n'est pas le problème).
     integer:: comm_master ! communicateur associé aux master, défini si master
     logical :: lmaster ! .true. si master
     integer:: grp_orig ! groupe du comm à diviser
     integer:: comm_orig ! communicatuer à diviser
     integer:: rang_orig ! rang dans le comm à diviser
     integer::np_orig ! pombre de procs dans le comm à dvisier
   contains
     procedure, pass::print


  end type para_config

#ifdef PARA

  
  interface distribue
     module procedure distribuereal, distribueinteger,distribuedouble
  end interface distribue

#endif
contains
    
  subroutine commconstr(div)
    type(para_config)::div
#ifdef PARA
    !  integer:: grp_world,grp_masters,comm_masters,npm,rgm,imasters
    integer, allocatable::rgmasters(:)
    integer::ierr,ip

    integer::imt,imp,ipi,ipt,npi,npim,imasters,npm,reste,npr
    integer::clef, couleur,nimage,img,image
    integer,allocatable::npimg(:),GL(:),CL(:),procim(:),ipimg(:,:)
    nimage=div%nimage
    if (div%nimage.gt.div%np_orig) then
       write(6,*)' division para impossible nimage > nprocs'
       stop
       call MPI_FINALIZE(ierr)
    end if
    reste=mod(div%np_orig,div%nimage)
    allocate(rgmasters(div%nimage))
    allocate(ipimg(nimage,div%np_orig))
    allocate(npimg(nimage))
    allocate(GL(0:nimage-1));allocate(CL(0:nimage-1))
    if (reste==0) then
       npim=div%np_orig/div%nimage            ! nb de proc par image
       div%image=int(div%rang_orig/npim) !quel image pour le proc
       div%rgim=mod(div%rang_orig,npim) ! quel rang dans le comm de l'image
       div%npim=npim
       div%lmaster=.false.     
       if (div%rgim==0) then
          div%lmaster=.true.
       end if
       ipi=0
       imasters=0
       npimg(:)=0
       ipimg(:,:)=-1
       do ip=0,div%np_orig-1
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
       clef=div%rang_orig

    else
     
       npim=div%np_orig/div%nimage
       npr=npim*(div%nimage-1)


       if (div%rang_orig.lt.npr) then
          div%npim=npim
          div%rgim=mod(div%rang_orig,npim)
          div%image=int(div%rang_orig/npim) !quel image pour le proc
       else
          div%image=div%nimage
          div%npim=npim+reste
          div%rgim=div%rang_orig-npr
       end if
       div%lmaster=.false.     
       if (div%rgim==0) then
          div%lmaster=.true.
       end if
       ipi=0
       imasters=0

       do ip=0,div%np_orig-1
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
       clef=div%rang_orig
    end if
Cl=0;GL=0
    do img=1,nimage
       if (div%rang_orig==0) write(6,*)img,ipimg(img,1:npimg(img))
       call MPI_GROUP_INCL(div%grp_orig,npimg(img),ipimg(img,1:npimg(img)),GL(img-1),ierr)
       call MPI_COMM_CREATE(div%comm_orig,GL(img-1),CL(img-1),ierr)
    end do
!    write(6,*)'TOTO', div%rang_orig,CL
!    write(6,*)'TATA', div%rang_orig,GL
!    write(6,*)'IMG',div%image,CL(div%image)
    call MPI_COMM_DUP(CL(div%image),div%comm_image,ierr)

    
!!$    call MPI_COMM_SPLIT (div%comm_orig,couleur,clef,div%comm_image,ierr)
!!$    call MPI_COMM_SIZE( div%comm_image, npi, ierr )
!!$    call MPI_COMM_RANK(div%comm_image, div%rgim,ierr)
    write(6,*)'RGI', div%rang_orig,div%rgim,div%lmaster,div%image,div%comm_image,CL(div%image)

    call MPI_BARRIER(div%comm_orig)
    
    call MPI_GROUP_INCL(div%grp_orig,div%nimage,rgmasters,div%grp_master,ierr)
    call MPI_COMM_CREATE(div%comm_orig,div%grp_master,div%comm_master,ierr)

    if (div%lmaster) then
       call MPI_COMM_RANK(div%comm_master, div%rgmas,ierr)
       call MPI_COMM_SIZE( div%comm_master, npm, ierr )
       if (npm.ne.div%nimage) then
          write(6,*)'NPMPB',npm,div%nimage
          stop
       end if
      write(6,*)'RGM',div%rang_orig, div%rgmas
    end if

    div%nimage=div%nimage
#endif
!        call MPI_BARRIER(div%comm_orig)
!        call mpi_finalize(ierr)
!        stop
    return

  end subroutine commconstr


  subroutine initparapuresp(div,rg,mpicsp,nps)
    type(para_config)::div
    integer,intent(in)::rg,nps
    integer,intent(in)::mpicsp
#ifdef PARA
    div%image=0
    div%nimage=1
    div%rgim=rg
    div%comm_image=mpicsp
    div%comm_orig=mpicsp
    if (rg==0)then
       div%lmaster=.true.
       div%rgmas=0
       div%comm_master=-1
    else
       div%lmaster=.false.
    end if
    div%rang_orig=rg
    div%npim=nps
#endif
  end subroutine initparapuresp
#ifdef PARA
  subroutine distribuereal(div,x,nel,xrecv)
    integer::nel
    type(para_config)::div
    real::x(0:div%nimage-1,nel)
    real,allocatable::xrecv(:)
    integer::img,ierr,isp
    integer, dimension( MPI_STATUS_SIZE) :: statut

    if (.not.(allocated(xrecv)))allocate(xrecv(1:nel))
    !     write(6,*)div
    !    if (div%lmaster) then
    !       write(6,*)'PRE',div%lmaster,div%rang_orig,div%rgim,div%rgmas
    !    else
    !       write(6,*)'PRE',div%lmaster,div%rang_orig,div%rgim
    !    end if
    if (div%rang_orig==0) then
       xrecv(:)=x(0,:)
       write(6,*)'rg0 ',xrecv
       do img=1,div%nimage-1
          !           write(6,*)'send 0->',img,x(img,:)
          call MPI_SEND (x(img,1:nel),nel, MPI_REAL ,img,100, div%comm_master,ierr)
       end do
       do isp=1,div%npim-1
          call MPI_SEND(x(0,1:nel),nel,MPI_REAL,isp,102,div%comm_image,ierr)
       end do
    else
       if (div%lmaster) then
          call MPI_RECV(xrecv(1:nel),nel,MPI_REAL,0,100,div%comm_master,statut,ierr)
          !           write(6,*)'recv->',xrecv(:),div%rgmas
          do isp=1,div%npim-1
             call MPI_SEND(xrecv(1:nel),nel,MPI_REAL,isp,101,div%comm_image,ierr)
          end do
       else
          if (div%image==0) then 
             call MPI_RECV(xrecv(1:nel),nel,MPI_REAL,0,102,div%comm_image,statut,ierr)
          else
             call MPI_RECV(xrecv(1:nel),nel,MPI_REAL,0,101,div%comm_image,statut,ierr)
          end if
       end if
    end if

  end subroutine distribuereal

  subroutine distribuedouble(div,x,nel,xrecv)
    integer::nel
    type(para_config)::div
    real::x(0:div%nimage-1,nel)
    real(double),allocatable::xrecv(:)
    integer::img,ierr,isp
    integer, dimension( MPI_STATUS_SIZE) :: statut

    if (.not.(allocated(xrecv)))allocate(xrecv(1:nel))
    !     write(6,*)div
    !    if (div%lmaster) then
    !       write(6,*)'PRE',div%lmaster,div%rang_orig,div%rgim,div%rgmas
    !    else
    !       write(6,*)'PRE',div%lmaster,div%rang_orig,div%rgim
    !    end if
    if (div%rang_orig==0) then
       xrecv(:)=x(0,:)
       write(6,*)'rg0 ',xrecv
       do img=1,div%nimage-1
          !           write(6,*)'send 0->',img,x(img,:)
          call MPI_SEND (x(img,1:nel),nel, NDM_MPI_REAL_DOUBLE ,img,100, div%comm_master,ierr)
       end do
       do isp=1,div%npim-1
          call MPI_SEND(x(0,1:nel),nel,NDM_MPI_REAL_DOUBLE,isp,102,div%comm_image,ierr)
       end do
    else
       if (div%lmaster) then
          call MPI_RECV(xrecv(1:nel),nel,NDM_MPI_REAL_DOUBLE,0,100,div%comm_master,statut,ierr)
          !           write(6,*)'recv->',xrecv(:),div%rgmas
          do isp=1,div%npim-1
             call MPI_SEND(xrecv(1:nel),nel,NDM_MPI_REAL_DOUBLE,isp,101,div%comm_image,ierr)
          end do
       else
          if (div%image==0) then 
             call MPI_RECV(xrecv(1:nel),nel,NDM_MPI_REAL_DOUBLE,0,102,div%comm_image,statut,ierr)
          else
             call MPI_RECV(xrecv(1:nel),nel,NDM_MPI_REAL_DOUBLE,0,101,div%comm_image,statut,ierr)
          end if
       end if
    end if
  end subroutine distribuedouble
  subroutine distribueinteger(div,x,nel,xrecv)
    integer::nel
    type(para_config)::div
    integer::x(0:div%nimage-1,nel)
    integer,allocatable::xrecv(:)
    integer::img,ierr,isp
    integer, dimension( MPI_STATUS_SIZE) :: statut

    if (.not.(allocated(xrecv)))allocate(xrecv(1:nel))
    !     write(6,*)div
    !    if (div%lmaster) then
    !       write(6,*)'PRE',div%lmaster,div%rang_orig,div%rgim,div%rgmas
    !    else
    !       write(6,*)'PRE',div%lmaster,div%rang_orig,div%rgim
    !    end if
    if (div%rang_orig==0) then
       xrecv(:)=x(0,:)
       write(6,*)'rg0 ',xrecv
       do img=1,div%nimage-1
          !           write(6,*)'send 0->',img,x(img,:)
          call MPI_SEND (x(img,1:nel),nel, MPI_INTEGER ,img,200, div%comm_master,ierr)
       end do
       do isp=1,div%npim-1
          call MPI_SEND(x(0,1:nel),nel,MPI_INTEGER,isp,202,div%comm_image,ierr)
       end do
    else
       if (div%lmaster) then
          call MPI_RECV(xrecv(1:nel),nel,MPI_INTEGER,0,200,div%comm_master,statut,ierr)
          !           write(6,*)'recv->',xrecv(:),div%rgmas
          do isp=1,div%npim-1
             call MPI_SEND(xrecv(1:nel),nel,MPI_INTEGER,isp,201,div%comm_image,ierr)
          end do
       else
          if (div%image==0) then 
             call MPI_RECV(xrecv(1:nel),nel,MPI_INTEGER,0,202,div%comm_image,statut,ierr)
          else
             call MPI_RECV(xrecv(1:nel),nel,MPI_INTEGER,0,201,div%comm_image,statut,ierr)
          end if
       end if
    end if
  end subroutine distribueinteger


#endif

  subroutine print(paraprt,rang)
    class(para_config),intent(in)::paraprt
    integer,intent(in)::rang
    write(6,*)'PARAPRT',rang,paraprt%rang_orig
    write(6,*)'NIMAGE',rang,paraprt%nimage
    write(6,*)'IMAGE',rang,paraprt%image
    write(6,*)'NPIMAGE',rang,paraprt%npim
    write(6,*)'RGIMAGE',rang,paraprt%rgim
    write(6,*)'LMASTER',rang,paraprt%lmaster
    if (paraprt%lmaster)write(6,*)'RGMASTER',rang,paraprt%rgmas
  end subroutine print
    


    
end module paraconfig



