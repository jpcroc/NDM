module deftypes
  USE T_kind_param_m
  implicit none
  integer:: incr=20 ! incrément des tailles de tableau 
  type atom_config ! type minimal des configurations atomiques. Tous les composants seront toujours alloué (im_glob seulement si PARA)
     integer::im=0
     real::potist
     real(double),allocatable:: xp(:,:)
     real(double),allocatable::fp(:,:)
     integer,allocatable::ityp(:)
     integer,allocatable::ielat(:)
     logical :: ltabvois
     integer,allocatable:: iwmax(:)
     integer,allocatable:: indi(:)
     logical, allocatable:: lgul(:)
#ifdef PARA
     integer,allocatable::num_at_glob(:)
#endif
   contains
     procedure, pass::init=>init_atom_config
     procedure, pass::copy_atom=>copy_atom_config
     procedure, pass::dealloc=>dealloc_atom_config
     procedure, pass::copy_config
     procedure, pass::print
     procedure, pass::pack
     procedure, pass::fab
#ifdef PARA
     procedure, pass::send2proc=>s2p_atom
     procedure, pass::send2all=>s2a_atom
     procedure, pass::recv=>rcv_atom

#endif     
  end type atom_config

  type, extends (atom_config):: atom_config_d ! type dynamique des configurations atomiques(+vp/+xpp). vp et xpp seront toujours allouées
     real(double),allocatable ::xpp(:,:)
     real(double),allocatable::vp(:,:)
   contains
     procedure, pass::copy_atom=>copy_atom_d
     procedure, pass::dealloc=>dealloc_atom_config_d
#ifdef PARA
     procedure, pass::send2proc=>s2p_atom_d
     procedure, pass::send2all=>s2a_atom_d
     procedure, pass::recv=>rcv_atom_d
#endif     
  end type atom_config_d
  
  type, extends (atom_config_d):: atom_config_e ! type étendu des configurations atomiques avec quantités optionelles Ces quantités seront allouées en fonction dss logical
     logical::lprteat
     real(double),allocatable ::eat(:)
     logical::lsigat
     real(double),allocatable ::sigat(:,:,:)
   contains
     procedure, pass::copy_atom=>copy_atom_e
     procedure, pass::dealloc=>dealloc_atom_config_e
#ifdef PARA
     procedure, pass::send2proc=>s2p_atom_e
     procedure, pass::send2all=>s2a_atom_e
     procedure, pass::recv=>rcv_atom_e
     
#endif     
  end type atom_config_e
  
contains
  !initialisations

  subroutine init_atom_config(atconf,imin,ltabvois,lsigat,lprteat)
    class(atom_config),intent(out)::atconf
    integer,intent(in):: imin
    logical,optional, intent(in)::ltabvois,lsigat,lprteat
    
    
    logical ::ltbv=.false.
    write(6,*)'init_e'

        write(6,*)'init'
    if (atconf%im.ne.0) then
       write(6,*)'atconf already allocated ; cannot be directly initiated ;  first deallocate'
       stop
    end if
    atconf%im=imin
    
    allocate(atconf%xp(3,imin));allocate(atconf%fp(3,imin));allocate(atconf%ityp(imin))
    allocate(atconf%ielat(imin));allocate(atconf%lgul(imin))
    atconf%ityp=0;atconf%xp=0;atconf%fp=0;atconf%ielat=0; atconf%lgul=.false.
#ifdef PARA    
    allocate(atconf%num_at_glob(im)) ; atconf%num_at_glob=0
#endif
    if(present(ltabvois)) ltbv=ltabvois
    if(ltbv)then
       atconf%ltabvois=.true.
       allocate(atconf%iwmax(atconf%im)); atconf%iwmax=0
    end if
    select type (atconf)
    type is (atom_config_d)
       write(6,*)'init_d'
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
    type(atom_config):: attemp
    integer:: imcn,nvav
    
    integer::iw,nvj,idwi,idwj
    
    if (present(lextend))let= lextend
    if (j.gt.atcible%im) then
       if (let) then
          imcn=atcible%im+incr
          call attemp%init(imcn,atcible%ltabvois)
          if (atcible%ltabvois)then
             nvav=imcn*int(1.1*(float(atcible%iwmax(atcible%im))/atcible%im+1))
             allocate(attemp%iwmax(nvav))
             !         write(6,*)'copy of an atom element is not possible , target size too small' 
             !         stop
          end if
          call atcible%copy_config(attemp,lrescl=.false.)
          call attemp%copy_config(atcible,lrescl=.true.)
       else
          write(6,*)'copy of an atom element is not possible , target size too small' 
          stop
       end if
    end if
    
    atcible%xp(:,j)=atsource%xp(:,i)
    atcible%fp(:,j)=atsource%fp(:,i)
    atcible%ityp(j)=atsource%ityp(i)
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
    atcible%ltabvois=atsource%ltabvois
    if ((atsource%ltabvois).and.(atcible%ltabvois)) then
       atcible%iwmax(1:atsource%im)=atsource%iwmax(1:atsource%im)
       atcible%indi(1:atsource%iwmax(atsource%im))=atsource%indi(1:atsource%iwmax(atsource%im))
    end if
! atsource et atcible sont moins _d
    select type (atsource)
    class is (atom_config_d)
       select type (atcible)
       class is (atom_config_d)
          atcible%vp(:,1:atsource%im)=atsource%xp(:,1:atsource%im)
          atcible%xpp(:,1:atsource%im)=atsource%fp(:,1:atsource%im)
       end select
    end select
! atsource et atcible sont moins _e    
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
       deallocate(atconf%ityp); deallocate(atconf%lgul)
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
  
  
  subroutine s2p_atom(atconf_trf,pcible)
    class(atom_config),intent(in)::atconf_trf
    integer,intent(in):: pcible
    !MPI_SEND de xp
    !MPI_SEND de fp
    !MPI_SEND de ityp
    !MPI_SEND de ielat
    !MPI_SEND de num_at_glob
    
  end subroutine s2p_atom
  
  subroutine s2p_atom_d(atomes_trf,pcible)
    class(atomes_types),intent(in)::atomes_trf
    integer,intent(in):: pcible
    
    call atomes_trf%atom_config%send2proc(pcible)
    ! MPI_SEND  de vp vers pcible
    ! MPI_SEND  de xpp vers pcible
    
  end subroutine s2p_atom_d
  
  
  subroutine s2p_atom_e(atomes_trf,pcible)
    class(atomes_types),intent(in)::atomes_trf
    integer,intent(in):: pcible
    call atomes_trf%atom_config_e%send2proc(pcible)
    if (allocated (ax)) then
       ! MPI_SEND  de ax vers pcible
    endif
    !    etc...
    
    
  end subroutine s2p_atom_e
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
  subroutine ndm2config (atndm,im,imm,potist,xp,fp,vp,xpp,ityp,ielat,ltabvois,iwmax,indi)
    type(atom_config_d)::atndm
    integer,intent(in)::im,imm
    real,intent(in)::potist
    real(double),intent(in),dimension(3,imm):: xp,vp,fp,xpp
    integer,intent(in),dimension(imm)::ityp,ielat
    logical, intent(in)::ltabvois
    integer, optional,intent(in) ::iwmax(imm)
    integer,optional,intent(in):: indi(:)
    integer::is
    call atndm%init(im,ltabvois)
    atndm%xp(:,1:im)=xp(:,1:im)
    atndm%vp(:,1:im)=vp(:,1:im)
    atndm%fp(:,1:im)=fp(:,1:im)
    atndm%xpp(:,1:im)=xpp(:,1:im)
    atndm%ityp(1:im)=ityp(1:im)
    atndm%ielat(1:im)=ielat(1:im)
    if (ltabvois) then
       atndm%iwmax(1:im)=iwmax(1:im)
       is=size(indi)
       allocate(atndm%indi(is))
       atndm%indi(1:is)=indi(1:is)
    end if
    atndm%potist=potist
  end subroutine ndm2config

  subroutine config2ndm (atndm,im,imm,potist,xp,fp,vp,xpp,ityp,ielat,ltabvois,iwmax,indi)
    type(atom_config_d),intent(in)::atndm
    integer,intent(out),allocatable::im
    integer, intent(in) ::imm
    real,intent(out)::potist
    real(double),intent(out),dimension(:,:),allocatable:: xp,vp,fp,xpp
    integer,intent(out),dimension(:),allocatable::ityp,ielat
    logical, intent(out)::ltabvois
    integer, optional,intent(out),allocatable ::iwmax(:)
    integer,optional,intent(out),allocatable:: indi(:)
    
    integer::is
    im=atndm%im
    potist=atndm%potist
    allocate(xp(3,imm));allocate(fp(3,imm));allocate(vp(3,imm));allocate(xpp(3,imm))
    allocate(ityp(imm));allocate(ielat(imm))
    xp(:,1:im)=atndm%xp(:,1:im)
    vp(:,1:im)=atndm%vp(:,1:im)
    fp(:,1:im)=atndm%fp(:,1:im)
    xpp(:,1:im)=atndm%xpp(:,1:im)
    ityp(1:im)=atndm%ityp(1:im)
    ielat(1:im)=atndm%ielat(1:im)
    if (atndm%ltabvois) then
       iwmax(1:im)=atndm%iwmax(1:im)
       is =size(atndm%indi)
       allocate(indi(is))
       indi(1:is)=atndm%indi(1:is)
    end if
  end subroutine config2ndm


  
end module deftypes


program test_generic
  use deftypes
  implicit none
  
  type(atom_config) :: a1,a2,a
  type(atom_config_d) :: d
  type(atom_config_e) :: e
  
  logical :: l
  integer:: i,im,imtrf,i2
  im=2
  l=.false.
  call a1%init(4,l)
  a1%xp=1
  a1%fp=2
  a1%ielat=3
  a1%ityp=6
  call a1%print
  call a1%copy_config(a2,.true.)
  write(6,*)
  write(6,*)'a2'
  call a2%print

  call d%init(im,l)
 
  call a1%copy_config(d,.true.)
  write(6,*)
  write(6,*)'d'
  call d%print



  call e%init(im,l,lprteat=.true.)
  
  call a1%copy_config(e,.true.)
  write(6,*)
  write(6,*)'e'
  call e%print

  write(6,*)'*******************************************'

  call e%dealloc
  write(6,*)
  write(6,*)'e2'  
  call e%init(2,lprteat=.true.,lsigat=.true.)
  e%xp=100; e%eat=101; e%sigat=102

  call e%print
  call e%dealloc
!  call e%init(3,lprteat=.true.)
  a1%ityp=12
  call a1%copy_atom(4,e,8,lextend=.true.)

  write(6,*)'e'  
  call e%print
  e%xp(1,1)=-40
  call e%copy_config(a1,lrescl=.true.)
  write(6,*)'a1'  
  call a1%print
  write(6,*)'a1pack'  
  call a1%pack
  call a1%print
  write(6,*)'TOTO LA'
  call a1%dealloc
  call a1%init(4,l)
  a1%xp=1
  a1%fp=2
  do i=1,4
     a1%ielat(i)=mod(i,2)
  end do
  a1%ityp=6
  write(6,*)'a1';call a1%print
  where (a1%ielat==1)
     a1%lgul=.true.
  elsewhere
     a1%lgul=.false.
  end where
  call a1%fab(a2)
!  call a2%dealloc
!  imtrf=COUNT(a1%lgul)
!  call a2%init(imtrf,a1%ltabvois)
!  i2= findloc(a1%lgul,.true.)
!  write(6,*)i2
!  i2=0
!  do i=1,a1%im
!     if(a1%lgul(i)) then
!        i2=i2+1
!        call a1%copy_atom(i,a2,i2,lextend=.false.)
!     end if
!  end do
!  if (i2.ne.imtrf) then
!     write(6,*)'WTF ?'
!     stop
!  end if
  call a2%print
end program test_generic
