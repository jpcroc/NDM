module atomconfig
  USE T_kind_param_m
  implicit none
  integer:: incr=20 ! incrément des tailles de tableau 
  type atom_config ! type minimal des configurations atomiques. Tous les composants seront toujours alloué (im_glob seulement si PARA)
     integer::im=0
     real(double),allocatable:: xp(:,:)
     real(double),allocatable::fp(:,:)
     integer,allocatable::ityp(:)
     integer,allocatable::ielat(:)
     logical :: ltabvois
     integer,allocatable:: iwmax(:)
     integer,allocatable:: indi(:)
     logical, allocatable:: lgul(:)
     integer,allocatable::num_at_glob(:)
   contains
     procedure, pass::init=>init_atom_config
     procedure, pass::copy_atom=>copy_atom_config
     procedure, pass::dealloc=>dealloc_atom_config
     procedure, pass::copy_config
     procedure, pass::print
     procedure, pass::pack
     procedure, pass::fab
     procedure, pass::add2conf
     procedure, pass::extend
!#ifdef PARA
!     procedure, pass::send2proc=>s2p_atom
!     procedure, pass::send2all=>s2a_atom
!     procedure, pass::recv=>rcv_atom
!
!#endif     
  end type atom_config

  type, extends (atom_config):: atom_config_d ! type dynamique des configurations atomiques(+vp/+xpp). vp et xpp seront toujours allouées
     real(double),allocatable ::xpp(:,:)
     real(double),allocatable::vp(:,:)
   contains
     procedure, pass::copy_atom=>copy_atom_d
     procedure, pass::dealloc=>dealloc_atom_config_d
!#ifdef PARA
!     procedure, pass::send2proc=>s2p_atom_d
!     procedure, pass::send2all=>s2a_atom_d
!     procedure, pass::recv=>rcv_atom_d
!#endif     
  end type atom_config_d
  
  type, extends (atom_config_d):: atom_config_e ! type étendu des configurations atomiques avec quantités optionelles Ces quantités seront allouées en fonction dss logical
     logical::lprteat
     real(double),allocatable ::eat(:)
     logical::lsigat
     real(double),allocatable ::sigat(:,:,:)
   contains
     procedure, pass::copy_atom=>copy_atom_e
     procedure, pass::dealloc=>dealloc_atom_config_e
!#ifdef PARA
!     procedure, pass::send2proc=>s2p_atom_e
!     procedure, pass::send2all=>s2a_atom_e
!     procedure, pass::recv=>rcv_atom_e
     
!#endif     
  end type atom_config_e
  
contains
  !initialisations

  subroutine init_atom_config(atconf,imin,ltabvois,lsigat,lprteat)
    class(atom_config),intent(out)::atconf
    integer,intent(in):: imin
    logical,optional, intent(in)::ltabvois,lsigat,lprteat
    
    
    logical ::ltbv=.false.
!    write(6,*)'init_e'

!        write(6,*)'init'
    if (atconf%im.ne.0) then
       write(6,*)'atconf already allocated ; cannot be directly initiated ;  first deallocate'
       stop
    end if
    atconf%im=imin
    
    allocate(atconf%xp(3,imin));allocate(atconf%fp(3,imin));allocate(atconf%ityp(imin))
    allocate(atconf%ielat(imin));allocate(atconf%lgul(imin));allocate(atconf%num_at_glob(imin))
    atconf%ityp=0;atconf%xp=0;atconf%fp=0;atconf%ielat=0; atconf%lgul=.false.;atconf%num_at_glob=0
    if(present(ltabvois)) ltbv=ltabvois
    if(ltbv)then
       atconf%ltabvois=.true.
       allocate(atconf%iwmax(atconf%im)); atconf%iwmax=0
    end if
    select type (atconf)
    type is (atom_config_d)
!       write(6,*)'init_d'
       allocate(atconf%vp(3,imin));allocate(atconf%xpp(3,imin))
       atconf%vp=0;atconf%xpp=0
    type is (atom_config_e)
       allocate(atconf%vp(3,imin));allocate(atconf%xpp(3,imin))
       atconf%vp=0;atconf%xpp=0

       if(present(lprteat)) atconf%lprteat=lprteat
       if(present(lsigat)) atconf%lsigat=lsigat
       write(6,*)'init_e',atconf%lprteat,atconf%lsigat
       if(atconf%lprteat)then
          allocate(atconf%eat(imin))
          atconf%eat=0
       end if
       if(atconf%lsigat)then
          allocate(atconf%sigat(3,3,imin));atconf%sigat=0
       end if
    end select
  end subroutine init_atom_config
  
  !copie d'un élément
  
  subroutine copy_atom_config(atsource,i,atcible,j,lextend)
    class(atom_config), intent(in)::atsource
    integer,intent(in):: i
    class(atom_config), intent(inout)::atcible
    integer,intent(in):: j
    logical , optional, intent(in) :: lextend
    logical ::let=.false.

    integer::iw,nvj,idwi,idwj
    
    if (present(lextend))let= lextend
    if (j.gt.atcible%im) then
       if (let) then
          call atcible%extend(incr)
       else
          write(6,*)'copy of an atom element is not possible , target size too small' 
          stop
       end if
    end if
    
    atcible%xp(:,j)=atsource%xp(:,i)
    atcible%fp(:,j)=atsource%fp(:,i)
    atcible%ityp(j)=atsource%ityp(i)
    atcible%num_at_glob(j)=atsource%num_at_glob(i)
    atcible%ielat(j)=atsource%ielat(i)
    atcible%lgul(j)=atsource%lgul(i)
    if ((atcible%ltabvois).and.(atsource%ltabvois))then
       nvj=atsource%iwmax(i)-atsource%iwmax(i-1)
       atcible%iwmax(j)=nvj+atsource%iwmax(j-1)
       do iw=1,nvj
          idwj=iw+atcible%iwmax(j-1); idwi=iw+atsource%iwmax(i-1)
          atcible%indi(idwj)=atsource%indi(idwi)
       end do
    end if

    
  end subroutine copy_atom_config

  subroutine extend (atcf,iadd)
    class(atom_config),intent(inout)::atcf
    integer,intent(in)::iadd
    
    type(atom_config):: attemp
    type(atom_config_d):: attemp_d
    type(atom_config_e):: attemp_e
    integer::imcn,nvav
    
    imcn=atcf%im+iadd

    select type(atcf)
    type is (atom_config)
       call attemp%init(imcn,atcf%ltabvois)
    type is (atom_config_d)
       call attemp_d%init(imcn,atcf%ltabvois)
    type is (atom_config_e)
       call attemp_e%init(imcn,atcf%ltabvois,atcf%lsigat,atcf%lprteat)
    end select
    if (atcf%ltabvois)then
       nvav=imcn*int(1.1*(float(atcf%iwmax(atcf%im))/atcf%im+1))
       allocate(attemp%iwmax(nvav))
    end if
    call atcf%copy_config(attemp,lrescl=.false.)
    call attemp%copy_config(atcf,lrescl=.true.)
       

  end subroutine extend
    
  
  
  subroutine copy_atom_d(atsource,i,atcible,j,lextend)
    class(atom_config_d), intent(in)::atsource
    integer,intent(in):: i
    class(atom_config), intent(inout)::atcible
    integer,intent(in):: j
    logical, optional,intent(in):: lextend
    logical:: let=.false.
    if (present(lextend))let=lextend
!    call atsource%atom_config%copy_atom(i,atcible%atom_config,j,let)
    call copy_atom_config(atsource,i,atcible,j,let)
    select type(atcible)
    class is (atom_config_d)
       select type (atsource)
       class is (atom_config_d)
          atcible%vp(:,j)=atsource%vp(:,i)
          atcible%xpp(:,j)=atsource%xpp(:,i)
       end select
    end select
  end subroutine copy_atom_d
  
  subroutine copy_atom_e(atsource,i,atcible,j,lextend)
    class(atom_config_e), intent(in)::atsource
    integer,intent(in):: i
    class(atom_config), intent(inout)::atcible
    integer,intent(in):: j
    logical, optional,intent(in):: lextend
    logical:: let=.false.
    if (present(lextend))let=lextend
    call copy_atom_d(atsource,i,atcible,j,let)
    select type(atcible)
    class is (atom_config_e)
       select type (atsource)
       type is (atom_config_e)
          if ((atcible%lprteat).and.(atsource%lprteat)) atcible%eat(j)=atcible%eat(i)
          if ((atcible%lsigat).and.(atsource%lsigat)) atcible%sigat(:,:,j)=atcible%sigat(:,:,i)
       end select
    end select
  end subroutine copy_atom_e
  
  
  ! copie d'une config entière vers config de base
  subroutine copy_config (atsource,atcible,lrescl)
    class(atom_config),intent(in)::atsource
    class(atom_config)::atcible
    logical,intent(in)::lrescl

    logical :: lstop
    if (lrescl) then
       if (atcible%im.ne.atsource%im) then
          call atcible%dealloc
          call atcible%init(atsource%im,atcible%ltabvois)
          if ((atcible%ltabvois).and.(atsource%ltabvois)) then
             allocate (atcible%indi(size(atsource%indi)))
          end if
       end if
    else
       lstop=.false.
       if ((atcible%im.lt.atsource%im)) lstop=.true.
       if (atcible%ltabvois) then
          if (size(atcible%indi).lt.size(atsource%indi)) lstop=.true.
       end if
       if (lstop)then
          write(6,*)'(atcible%im < atsource%im  ou size(atcible%indi)<size(atsource%indi) )et lrescl = false ; pas possible'
          stop
       end if
    end if
    atcible%xp(:,1:atsource%im)=atsource%xp(:,1:atsource%im)
    atcible%fp(:,1:atsource%im)=atsource%fp(:,1:atsource%im)
    atcible%ielat(1:atsource%im)=atsource%ielat(1:atsource%im)
    atcible%lgul(1:atsource%im)=atsource%lgul(1:atsource%im)
    atcible%ityp(1:atsource%im)=atsource%ityp(1:atsource%im)
    atcible%num_at_glob(1:atsource%im)=atsource%num_at_glob(1:atsource%im)
    atcible%ltabvois=atsource%ltabvois
    if ((atsource%ltabvois).and.(atcible%ltabvois)) then
       atcible%iwmax(1:atsource%im)=atsource%iwmax(1:atsource%im)
       atcible%indi(1:atsource%iwmax(atsource%im))=atsource%indi(1:atsource%iwmax(atsource%im))
    end if
! atsource et atcible sont au moins_d
    select type (atsource)
    class is (atom_config_d)
       select type (atcible)
       class is (atom_config_d)
          atcible%vp(:,1:atsource%im)=atsource%xp(:,1:atsource%im)
          atcible%xpp(:,1:atsource%im)=atsource%fp(:,1:atsource%im)
       end select
    end select
! atsource et atcible sont _e    
    select type (atsource)
    class is (atom_config_e)
       select type (atcible)
       class is (atom_config_e)
          if((atcible%lprteat).and.(atsource%lprteat))atcible%eat(1:atsource%im)=atsource%eat(1:atsource%im)
          if((atcible%lsigat).and.(atsource%lsigat))atcible%sigat(:,:,1:atsource%im)=atsource%sigat(:,:,1:atsource%im)
       end select
    end select
  end subroutine copy_config
  
  subroutine dealloc_atom_config(atconf)
    class(atom_config), intent(inout)::atconf
    atconf%im=0
    if(allocated(atconf%xp))then
       deallocate(atconf%xp);deallocate(atconf%fp);deallocate(atconf%ielat)
       deallocate(atconf%ityp); deallocate(atconf%lgul); deallocate(atconf%num_at_glob)
    end if
    if (allocated (atconf%iwmax))deallocate (atconf%iwmax)
    if (allocated (atconf%indi))deallocate (atconf%indi)
    return
  end subroutine dealloc_atom_config
  
  subroutine dealloc_atom_config_d(atconf)
    class(atom_config_d), intent(inout)::atconf
    call atconf%atom_config%dealloc
    if(allocated(atconf%vp))then
       deallocate(atconf%vp); deallocate(atconf%xpp)
    end if
  end subroutine dealloc_atom_config_d
  subroutine dealloc_atom_config_e(atconf)
    class(atom_config_e), intent(inout)::atconf
    call atconf%atom_config_d%dealloc
    if(allocated(atconf%eat))deallocate(atconf%eat)
    if(allocated(atconf%sigat))deallocate(atconf%sigat)
    
  end subroutine dealloc_atom_config_e
  
#ifdef PARA  
  
  
!  subroutine s2p_atom(atconf_trf,pcible)
!    class(atom_config),intent(in)::atconf_trf
!    integer,intent(in):: pcible
    !MPI_SEND de xp
    !MPI_SEND de fp
    !MPI_SEND de ityp
    !MPI_SEND de ielat
    !MPI_SEND de num_at_glob
    
!  end subroutine s2p_atom
  
 ! subroutine s2p_atom_d(atomes_trf,pcible)
 !   class(atomes_types),intent(in)::atomes_trf
 !   integer,intent(in):: pcible
    
 !   call atomes_trf%atom_config%send2proc(pcible)
    ! MPI_SEND  de vp vers pcible
    ! MPI_SEND  de xpp vers pcible
    
 ! end subroutine s2p_atom_d
  
  
 ! subroutine s2p_atom_e(atomes_trf,pcible)
 !   class(atomes_types),intent(in)::atomes_trf
 !   integer,intent(in):: pcible
 !   call atomes_trf%atom_config_e%send2proc(pcible)
 !   if (allocated (ax)) then
       ! MPI_SEND  de ax vers pcible
 !   endif
    !    etc...
    
    
 ! end subroutine s2p_atom_e
#endif
  subroutine pack(at2pack)
    class(atom_config),intent(inout):: at2pack
    type(atom_config)::at
    type(atom_config_d)::atd
    type(atom_config_e)::ate
    integer::i,imn
       if (at2pack%ityp(at2pack%im).ne.0) return
    do i=at2pack%im-1,1,-1
       if (at2pack%ityp(i).ne.0)then
          imn=i
          exit
       end if
    end do
    select type (at2pack)
    type is (atom_config)
       call at%init(imn,at2pack%ltabvois)
       do i=1,imn
          call at2pack%copy_atom(i,at,i)
       end do
       call at2pack%dealloc
       call at%copy_config(at2pack,lrescl=.true.)
       
    type is(atom_config_d)
       call atd%init(imn,at2pack%ltabvois)
       do i=1,imn
          call at2pack%copy_atom(i,atd,i)
       end do
       call at2pack%dealloc
       call atd%copy_config(at2pack,lrescl=.true.)

    type is(atom_config_e)
       call ate%init(imn,at2pack%ltabvois,at2pack%lsigat,at2pack%lprteat)
       do i=1,imn
          call at2pack%copy_atom(i,ate,i)
       end do
          call at2pack%dealloc
       call ate%copy_config(at2pack,lrescl=.true.)
    end select

  end subroutine pack
    
  subroutine fab (atsource,atcible)
    class(atom_config),intent(in)::atsource
    class(atom_config),intent(out)::atcible

    integer::i2,imtrf,i
    call atcible%dealloc 

    imtrf=COUNT(atsource%lgul)
    select type (atsource)
    type is (atom_config_e)
       call atcible%init(imtrf,atsource%ltabvois,lsigat=atsource%lsigat,lprteat=atsource%lprteat)
    class is (atom_config)
       call atcible%init(imtrf,atsource%ltabvois)
    end select
    i2=0
    do i=1,atsource%im
       if(atsource%lgul(i)) then
          i2=i2+1
          call atsource%copy_atom(i,atcible,i2,lextend=.false.)
       end if
    end do
    if (i2.ne.imtrf) then
       write(6,*)'WTF ?'
       stop
    end if
  end subroutine fab
  
  
  subroutine add2conf (atsource,atcible,ldealloc)
    class(atom_config),intent(inout)::atsource
    class(atom_config),intent(inout)::atcible
    logical,optional,intent(in)::ldealloc


    
    logical :: ldal
    
    integer::i2,imtrf,i,imdeb,imnew
!    call atcible%dealloc 

    !    imtrf=COUNT(atsource%lgul)
    if(present(ldealloc)) ldal=ldealloc
    imtrf=atsource%im
    imdeb=atcible%im
    imnew=imtrf+imdeb

    if(atcible%im==0) then
       select type (atsource)
       type is (atom_config_e)
          call atcible%init(imtrf,atsource%ltabvois,lsigat=atsource%lsigat,lprteat=atsource%lprteat)
       class is (atom_config)
          call atcible%init(imtrf,atsource%ltabvois)
       end select
    else
       call atcible%extend(imnew)
    end if

    
    i2=0
    do i=1,atsource%im
       i2=imdeb+i
       call atsource%copy_atom(i,atcible,i2,lextend=.false.)
    end do
    if (i2.ne.atcible%im) then
       write(6,*)'WTF ?'
       stop
    end if
    if (ldal) call atsource%dealloc
  end subroutine add2conf
  
  
    
  subroutine print(atprt)
    class(atom_config), intent(in)::atprt
    type(atom_config_d):: td
    type(atom_config_e):: te
    integer::i,ic,im
    write(6,*)
    write(6,*)'in print'
     im=atprt%im
    
    write(6,*)'im = ',atprt%im
    write(6,*)'ltabvois ', atprt%ltabvois
    if (allocated(atprt%xp)) then
       do i=1,im
          write(6,*)'%xp= ', atprt%xp(:,i)
       end do
       do i=1,im
          write(6,*)'%fp= ', atprt%fp(:,i)
       end do
       do i=1,im
          write(6,*)'%ityp= ', atprt%ityp(i)
       end do
       do i=1,im
          write(6,*)'%ielat= ', atprt%ielat(i)
       end do
!       if (extends_type_of(atprt,td)) then
!          write(6,*)'prt_d'
!          do i=1,im
!             write(6,*)'%vp= ', atprt%vp(:,i)
!          end do
!          do i=1,im
!             write(6,*)'%xpp= ', atprt%xpp(:,i)
!          end do
!       end if
!       if (same_type_as(atprt,te))then
!                    if (atprt%lsigat) then
!             do i=1,im
!                write(6,*)'%sigat= ', atprt%sigat(:,:,i)
!             end do
!          end if
!          if (atprt%lprteat) then
!             do i=1,im
!                write(6,*)'%eat= ', atprt%eat(i)
!             end do
!          end if
!       end if
       select type (atprt)
       class is (atom_config_d)
          write(6,*)'prt_d'
          do i=1,im
             write(6,*)'%vp= ', i,atprt%vp(:,i)
          end do
          do i=1,im
             write(6,*)'%xpp= ', i,atprt%xpp(:,i)
          end do
       class is (atom_config_e)
          write(6,*)'prt_e'
          do i=1,im
             write(6,*)'%vp= ', i,atprt%vp(:,i)
          end do
          do i=1,im
             write(6,*)'%xpp= ', i,atprt%xpp(:,i)
          end do
          
          
          if (atprt%lsigat) then
             do i=1,im
                write(6,*)'%sigat= ',i, atprt%sigat(:,:,i)
             end do
          end if
          if (atprt%lprteat) then
             do i=1,im
                write(6,*)'%eat= ', i,atprt%eat(i)
             end do
          end if
       end select


       if (atprt%ltabvois) then
          do i=1,im
             write(6,*)'%iwmax= ', i,atprt%iwmax(i)
          end do
          
          if (allocated(atprt%indi))then
             do i=1,size(atprt%indi)
                write(6,*)'vois', i,atprt%indi(i)
             end do
          end if
       end if
    end if
    write(6,*)'out print'
    write(6,*)
  end subroutine print
  !MANQUE SIG AU MINIMUM
  subroutine ndm2config (atndm,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax,indi,nvois,vp,xpp,eat,sigat)
    class(atom_config)::atndm
    integer,intent(in)::im,imm
    real(double),intent(in),dimension(3,imm):: xp,fp
    integer,intent(in),dimension(imm)::ityp,ielat
    integer,optional,intent(in),dimension(imm)::num_at_glob
    logical, optional,intent(in)::ltabvois
    integer,optional,intent(in)::nvois
    integer, optional,intent(in) ::iwmax(imm)
    integer,optional,intent(in),allocatable:: indi(:)
    real(double),optional,intent(in),dimension(3,imm):: vp,xpp
    real(double),optional,intent(in):: eat(imm),sigat(3,3,imm)

    integer::is
    logical ::lprteat=.false.,lsigat=.false.,ltbv=.false.


    if (present(eat))lprteat=.true.; if(present(sigat))lsigat=.true.;  if(present(ltabvois))ltbv=ltabvois
!    write(6,*)'in ndm2conf'
    select type(atndm)
!    type is (atom_config)
    class is (atom_config)
       call atndm%init(im,ltabvois)
    type is (atom_config_e)
       call atndm%init(im,ltabvois,lsigat,lprteat)
    end select

    atndm%xp(:,1:im)=xp(:,1:im)
    atndm%fp(:,1:im)=fp(:,1:im)
    atndm%ityp(1:im)=ityp(1:im)
    if (present(num_at_glob))atndm%num_at_glob(1:im)=num_at_glob(1:im)
    atndm%ielat(1:im)=ielat(1:im)
    if (ltbv) then
!       write(6,*)'sizes ', size (iwmax),size(atndm%iwmax)
       atndm%ltabvois=.true.
       atndm%iwmax(1:im)=iwmax(1:im)
       allocate(atndm%indi(nvois))
       atndm%indi(1:nvois)=indi(1:nvois)
    end if

    select type(atndm)
    type is (atom_config_d)
       if (present(vp))atndm%vp(:,1:im)=vp(:,1:im)
       if(present(xpp))atndm%xpp(:,1:im)=xpp(:,1:im)

    type is (atom_config_e)
       if (present(sigat).and.(atndm%lsigat))then
        !  allocate(sigat(3,3,imm))
          atndm%sigat(:,:,1:atndm%im)=sigat(:,:,1:atndm%im)
       end if
       if (present(eat).and.(atndm%lprteat))then
!          allocate(eat(imm))
          atndm%eat(1:atndm%im)=eat(1:atndm%im)
       end if
    end select

  end subroutine ndm2config

  subroutine config2ndm (atndm,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax,indi,vp,xpp,eat,sigat)
    class(atom_config),intent(in)::atndm
    integer,intent(inout)::im
    integer,intent(in)::imm
    real(double),intent(inout),allocatable:: xp(:,:),fp(:,:)

    integer,intent(inout),dimension(:),allocatable::ityp,num_at_glob,ielat
    logical, intent(out)::ltabvois
    integer, optional,intent(inout),allocatable ::iwmax(:)
    integer,optional,intent(inout),allocatable:: indi(:)
    
    real(double),optional,intent(inout),allocatable:: vp(:,:),xpp(:,:)
    
    real(double),optional,intent(inout),allocatable::eat(:),sigat(:,:,:)
!    integer,intent(in)::nvois
    integer::is
!    write(6,*)'in conf2ndm'
    im=atndm%im
!    allocate(xp(3,imm));allocate(fp(3,imm));allocate(vp(3,imm));allocate(xpp(3,imm))
!    allocate(ityp(imm));allocate(ielat(imm));allocate(num_at_glob(imm))
    xp(:,1:im)=atndm%xp(:,1:im)

    fp(:,1:im)=atndm%fp(:,1:im)

    ityp(1:im)=atndm%ityp(1:im)
    num_at_glob(1:im)=atndm%num_at_glob(1:im)
    ielat(1:im)=atndm%ielat(1:im)
    if (atndm%ltabvois) then
!       write(6,*)'sizes ', size (iwmax),size(atndm%iwmax)
!       allocate(iwmax(imm))
       iwmax(1:im)=atndm%iwmax(1:im)
       is =size(atndm%indi)
!       write(6,*)'IS',is
!       allocate(indi(is))
       indi(1:is)=atndm%indi(1:is)
    end if
    select type(atndm)
    type is (atom_config_d)
       if(present(vp))vp(:,1:im)=atndm%vp(:,1:im)
       if(present(xpp))xpp(:,1:im)=atndm%xpp(:,1:im)
    type is (atom_config_e)
       if (present(sigat).and.(atndm%lsigat))then
!          allocate(sigat(3,3,imm))
          sigat(:,:,1:atndm%im)=atndm%sigat(:,:,1:atndm%im)
       end if
       if (present(eat).and.(atndm%lprteat))then
!          allocate(eat(imm))
          eat(1:atndm%im)=atndm%eat(1:atndm%im)
       end if
    end select
  end subroutine config2ndm


  
end module atomconfig


