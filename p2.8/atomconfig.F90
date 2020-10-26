module atomconfig
  USE T_kind_param_m
  implicit none
  integer:: incr=20 ! incrément des tailles de tableau 
  type atom_config ! type minimal des configurations atomiques. Tous les composants seront toujours alloué (im_glob seulement si PARA)
     integer::im=0,imm=0
     integer(long)::icaltabt 
     real(double),allocatable:: xp(:,:)
     real(double),allocatable::fp(:,:)
     integer,allocatable::ityp(:)
     integer,allocatable::ielat(:)
     logical :: ltabvois
     integer,allocatable:: iwmax(:) ! indice du dernier voisin de i
     integer:: nvois ! taille du tableau des voisins si ltabvois=.true.
     integer,allocatable:: indi(:) ! indice de tous les voisins
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
  
  type, extends (atom_config_d):: atom_config_e ! type étendu des configurations atomiques avec quantités optionelles Ces quantités seront allouées en fonction dss logical
     logical::lprteat
     real(double),allocatable ::eat(:)
     logical::lsigat
     real(double),allocatable ::sigat(:,:,:)
     logical::lLangevin
     real(double),allocatable ::Glangv(:,:)
     logical::lax
     real(double),allocatable ::ax(:,:)
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
  
  subroutine init_atom_config(atconf,imin,immin,ltabvois,nvois,lsigat,lprteat,lLangevin,lax)
    class(atom_config),intent(out)::atconf
    integer,intent(in):: imin
    logical,optional, intent(in)::ltabvois,lsigat,lprteat,lLangevin,lax
    integer, optional::nvois,immin

    logical ::ltbv
    ltbv=.false.
    if (atconf%im.ne.0) then
       write(6,*)'atconf already allocated ; cannot be directly initiated ;  first deallocate'
       stop
    end if
    atconf%im=imin       
    if (present(immin))then
       atconf%imm=immin
    else
       atconf%imm=imin
    end if

    allocate(atconf%xp(3,atconf%imm));allocate(atconf%fp(3,atconf%imm));allocate(atconf%ityp(atconf%imm))
    allocate(atconf%ielat(atconf%imm));allocate(atconf%lgul(atconf%imm));allocate(atconf%num_at_glob(atconf%imm))
    atconf%ityp=0;atconf%xp=0;atconf%fp=0;atconf%ielat=0; atconf%lgul=.false.;atconf%num_at_glob=0
    if(present(ltabvois)) ltbv=ltabvois
    if(ltbv)then
       atconf%ltabvois=.true.
       allocate(atconf%iwmax(atconf%imm)); atconf%iwmax=0
       if (nvois.ne.0) then
          atconf%nvois=nvois
          allocate(atconf%indi(nvois))
       end if
    end if
    select type (atconf)
    type is (atom_config_d)
       !       write(6,*)'init_d'
       allocate(atconf%vp(3,atconf%imm));allocate(atconf%xpp(3,atconf%imm))
       atconf%vp=0;atconf%xpp=0
    type is (atom_config_e)
       allocate(atconf%vp(3,atconf%imm));allocate(atconf%xpp(3,atconf%imm))
       atconf%vp=0;atconf%xpp=0
       atconf%lprteat=.false.
       atconf%lsigat=.false.
       atconf%lLangevin=.false.
       atconf%lax=.false.
       if(present(lprteat)) atconf%lprteat=lprteat
       if(present(lsigat)) atconf%lsigat=lsigat
       if(present(lLangevin)) atconf%lLangevin=lLangevin
       if(present(lax)) atconf%lax=lax
       !       write(6,*)'init_e',atconf%lprteat,atconf%lsigat
       if(atconf%lprteat)then
          allocate(atconf%eat(atconf%imm))
          atconf%eat=0
       end if
       if(atconf%lsigat)then
          allocate(atconf%sigat(3,3,atconf%imm));atconf%sigat=0
       end if
       if(atconf%llangevin)then
          allocate(atconf%Glangv(3,atconf%imm));atconf%Glangv=0
       end if
       if(atconf%lax)then
          allocate(atconf%ax(3,atconf%imm));atconf%ax=0
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
    logical ::let
    integer::iw,nvj,idwi,idwj,imm_min
    let=.false.
    if (present(lextend))let= lextend
    if (j.gt.atcible%imm) then
       if (let) then
          imm_min=j-atcible%imm
          call atcible%extend(imm_min)
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

  subroutine extend(atcf,iadd)
    class(atom_config),intent(inout)::atcf
    integer,intent(in)::iadd

    type(atom_config):: attemp
    type(atom_config_d):: attemp_d
    type(atom_config_e):: attemp_e
    integer::imcn,nvois,immcn
    immcn=atcf%imm+iadd
    imcn=atcf%im+iadd
    if (atcf%ltabvois)then
       nvois=int(1.1*imcn/atcf%im)*size(atcf%indi)
    else
       nvois=0
    end if
    select type(atcf)
    type is (atom_config)
       call attemp%init(imcn,immcn,atcf%ltabvois,nvois)
    type is (atom_config_d)
       call attemp_d%init(imcn,immcn,atcf%ltabvois,nvois)
    type is (atom_config_e)
       call attemp_e%init(imcn,immcn,atcf%ltabvois,nvois,atcf%lsigat,atcf%lprteat, atcf%llangevin,atcf%lax)
    end select
    call atcf%copy_config(attemp,lrescl=.false.)
    call attemp%copy_config(atcf,lrescl=.true.)


  end subroutine extend



  subroutine copy_atom_d(atsource,i,atcible,j,lextend)
    class(atom_config_d), intent(in)::atsource
    integer,intent(in):: i
    class(atom_config), intent(inout)::atcible
    integer,intent(in):: j
    logical, optional,intent(in):: lextend
    logical::let
    let=.false.
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
    logical:: let
    let=.false.
    if (present(lextend))let=lextend
    call copy_atom_d(atsource,i,atcible,j,let)
    select type(atcible)
       class is (atom_config_e)
       select type (atsource)
       type is (atom_config_e)
          if ((atcible%lprteat).and.(atsource%lprteat)) atcible%eat(j)=atcible%eat(i)
          if ((atcible%lsigat).and.(atsource%lsigat)) atcible%sigat(:,:,j)=atcible%sigat(:,:,i)
          if ((atcible%llangevin).and.(atsource%llangevin)) atcible%glangv(:,j)=atcible%glangv(:,i)
          if ((atcible%lax).and.(atsource%lax)) atcible%ax(:,j)=atcible%ax(:,i)
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
       if (atcible%imm.ne.atsource%imm) then
          call atcible%dealloc
          call atcible%init(atsource%im,atsource%imm,atcible%ltabvois,size(atsource%indi))
          atcible%icaltabt=atsource%icaltabt    
          !          if ((atcible%ltabvois).and.(atsource%ltabvois)) then
          !             allocate (atcible%indi(size(atsource%indi)))
          !          end if
       end if
    else
       lstop=.false.
       if ((atcible%im.lt.atsource%im).or.(atcible%imm.lt.atsource%imm)) lstop=.true.
       if (atcible%ltabvois) then
          if (size(atcible%indi).lt.size(atsource%indi)) lstop=.true.
       end if
       if (lstop)then
          write(6,*)'(atcible%im < atsource%im  ou size(atcible%indi)<size(atsource%indi) )et lrescl = false ; pas possible'
          stop
       end if
    end if
    !    atcible%im=atsource%im si rescale a deje été fait. Si pas rescale n'a pas a être fait. 
    !    atcible%imm=atsource%imm
    !    atcible%icaltabt=atsource%icaltabt    
    atcible%xp(:,1:atsource%imm)=atsource%xp(:,1:atsource%imm)
    atcible%fp(:,1:atsource%imm)=atsource%fp(:,1:atsource%imm)
    atcible%ielat(1:atsource%imm)=atsource%ielat(1:atsource%imm)
    atcible%lgul(1:atsource%imm)=atsource%lgul(1:atsource%imm)
    atcible%ityp(1:atsource%imm)=atsource%ityp(1:atsource%imm)
    atcible%num_at_glob(1:atsource%imm)=atsource%num_at_glob(1:atsource%imm)
    atcible%ltabvois=atsource%ltabvois
    if ((atsource%ltabvois).and.(atcible%ltabvois)) then
       atcible%iwmax(1:atsource%imm)=atsource%iwmax(1:atsource%imm)
       atcible%indi(1:atsource%iwmax(atsource%imm))=atsource%indi(1:atsource%iwmax(atsource%imm))
    end if
    ! atsource et atcible sont au moins_d
    select type (atsource)
       class is (atom_config_d)
       select type (atcible)
          class is (atom_config_d)
          atcible%vp(:,1:atsource%imm)=atsource%vp(:,1:atsource%imm)
          atcible%xpp(:,1:atsource%imm)=atsource%xp(:,1:atsource%imm)
       end select
    end select
    ! atsource et atcible sont _e    
    select type (atsource)
       class is (atom_config_e)
       select type (atcible)
          class is (atom_config_e)
          if((atcible%lprteat).and.(atsource%lprteat))atcible%eat(1:atsource%imm)=atsource%eat(1:atsource%imm)
          if((atcible%lsigat).and.(atsource%lsigat))atcible%sigat(:,:,1:atsource%imm)=atsource%sigat(:,:,1:atsource%imm)
          if((atcible%llangevin).and.(atsource%llangevin))atcible%glangv(:,1:atsource%imm)=atsource%glangv(:,1:atsource%imm)
          if((atcible%lax).and.(atsource%lax))atcible%ax(:,1:atsource%imm)=atsource%ax(:,1:atsource%imm)
       end select
    end select
  end subroutine copy_config

  subroutine dealloc_atom_config(atconf)
    class(atom_config), intent(inout)::atconf
    atconf%im=0 ; atconf%imm=0
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
    if(allocated(atconf%glangv))deallocate(atconf%glangv)
    if(allocated(atconf%ax))deallocate(atconf%ax)

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
  subroutine pack(at2pack,imm_in)
    class(atom_config),intent(inout):: at2pack
    integer,optional, intent(in):: imm_in
    type(atom_config)::at
    type(atom_config_d)::atd
    type(atom_config_e)::ate
    integer::i,imn,immn

    if (at2pack%ityp(at2pack%imm).ne.0) return
    do i=at2pack%imm-1,1,-1
       if (at2pack%ityp(i).ne.0)then
          imn=i
          exit
       end if
    end do
    if (present(imm_in))then
       if (imm_in.ge.imn) then
          immn=imm_in
       else
          write(6,*)'imm_in< imn ; stop'
          stop
       end if
    else
       immn=imn
    end if
    select type (at2pack)
    type is (atom_config)
       call at%init(imn,immn,at2pack%ltabvois,size(at2pack%indi))
       do i=1,immn
          call at2pack%copy_atom(i,at,i)
       end do
       call at2pack%dealloc
       call at%copy_config(at2pack,lrescl=.true.)

    type is(atom_config_d)
       call atd%init(imn,immn,at2pack%ltabvois,size(at2pack%indi))
       do i=1,immn
          call at2pack%copy_atom(i,atd,i)
       end do
       call at2pack%dealloc
       call atd%copy_config(at2pack,lrescl=.true.)

    type is(atom_config_e)
       call ate%init(imn,immn,at2pack%ltabvois,size(at2pack%indi),at2pack%lsigat,at2pack%lprteat,llangevin=at2pack%llangevin,&
            &lax=at2pack%lax)
       do i=1,immn
          call at2pack%copy_atom(i,ate,i)
       end do
       call at2pack%dealloc
       call ate%copy_config(at2pack,lrescl=.true.)
    end select

  end subroutine pack

  subroutine fab (atsource,atcible) ! construit atsource à partir de lgul de atcible , ecrase atcible
    class(atom_config),intent(in)::atsource
    class(atom_config),intent(out)::atcible

    integer::i2,imtrf,i
    call atcible%dealloc 

    imtrf=COUNT(atsource%lgul(1:atsource%im))

    select type (atsource)
    type is (atom_config_e)
       call atcible%init(imtrf,imtrf,atsource%ltabvois,size(atsource%indi),lsigat=atsource%lsigat,lprteat=atsource%lprteat,&
            &llangevin=atsource%llangevin,lax=atsource%lax)
       class is (atom_config)
       call atcible%init(imtrf,imtrf,atsource%ltabvois,size(atsource%indi))
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


  subroutine add2conf (atsource,atcible,lextend,ldealloc)
    class(atom_config),intent(inout)::atsource
    class(atom_config),intent(inout)::atcible
    logical,optional,intent(in)::ldealloc,lextend



    logical :: ldal ! par defaut pas de destruction de atsource
    logical :: lext ! par defaut on etend atcible si besoin

    integer::i2,i,imcib,imnew,immcib,imsrc,immsrc,immnew
    type(atom_config):: atcor
    type(atom_config_d):: atcor_d
    type(atom_config_e):: atcor_e
    ldal=.false.
    lext=.true.

    !    call atcible%dealloc 

    !    imtrf=COUNT(atsource%lgul)
    if(present(ldealloc)) ldal=ldealloc
    if(present(lextend)) lext=lextend

    immsrc=atsource%imm
    imsrc=atsource%im
    imcib=atcible%im
    immcib=atcible%imm
    imnew=imsrc+imcib
    immnew=immsrc+immcib

    if(atcible%imm==0) then
       select type (atsource)
       type is (atom_config_e)
          call atcible%init(imsrc,immsrc,atsource%ltabvois,size(atsource%indi),lsigat=atsource%lsigat,&
               &lprteat=atsource%lprteat,llangevin=atsource%llangevin,lax=atsource%lax)
          class is (atom_config)
          call atcible%init(imsrc,immsrc,atsource%ltabvois,size(atsource%indi))
       end select
    else
       select type(atcible)
       type is (atom_config)
          call atcor%init(imcib,immcib,atcible%ltabvois,size(atcible%indi))
          call atcible%copy_config(atcor,.false.)
       type is (atom_config_d)
          call atcor_d%init(imcib,immcib,atcible%ltabvois,size(atcible%indi))
          call atcible%copy_config(atcor_d,.false.)
       type is (atom_config_e)
          call atcor_e%init(imcib,immcib,atcible%ltabvois,size(atcible%indi),atcible%lsigat,atcible%lprteat,&
               &llangevin=atcible%llangevin,lax=atcible%lax)
          call atcible%copy_config(atcor_e,.false.)
       end select
       if (lext) then
          if (atcible%imm.lt.immnew)   call atcible%extend(immnew)
       else
          if (atcible%imm.lt.imnew)  then
             write(6,*)'addition de atsource a atcible pas possible'
             stop
          end if
       end if
       atcible%im=imsrc+imcib
       atcible%imm=immsrc+immcib
    end if


    do i=1,atsource%im
       i2=imcib+i
       call atsource%copy_atom(i,atcible,i2,lextend=.false.)
    end do

    select type(atcible)
    type is (atom_config)
       do i=atcor%im+1,atcor%imm
          i2=imcib+imsrc+i
          call atcor%copy_atom(i,atcible,i2,lextend=.false.)
       end do
    type is (atom_config_d)
       do i=atcor_d%im+1,atcor_d%imm
          i2=imcib+imsrc+i
          call atcor_d%copy_atom(i,atcible,i2,lextend=.false.)
       end do

    type is (atom_config_e)
       do i=atcor_e%im+1,atcor_e%imm
          i2=imcib+imsrc+i
          call atcor_e%copy_atom(i,atcible,i2,lextend=.false.)
       end do

    end select

    do i=atsource%im+1,atsource%imm
       i2=immcib+imsrc+i
       call atsource%copy_atom(i,atcible,i2,lextend=.false.)
    end do

    !    if (i2.ne.atcible%imm) then
    !       write(6,*)'WTF ?'
    !       stop
    !    end if
    if (ldal) call atsource%dealloc
  end subroutine add2conf



  subroutine print(atprt,i1,i2,iwr)
    class(atom_config), intent(in)::atprt
    integer,optional,intent(in)::i1,i2,iwr
    !    type(atom_config_d):: td
    !    type(atom_config_e):: te
    integer::i,im,ifin,ideb,ist,ifn,iw
    !    write(6,*)
    write(6,*)'in print'
    im=atprt%im
    iw=1
    if (present(iwr))iw=iwr

    write(6,*)'im = ',atprt%im
    write(6,*)'imm = ',atprt%imm
    write(6,*)'icaltabt = ',atprt%icaltabt
    write(6,*)'ltabvois ', atprt%ltabvois

    if (present(i1))then
       ideb=i1
    else
       ideb=1
    endif
    if (present(i2)) then
       ifin=i2
    else
       if (present(i1))then
          ifin=i1
       else
          ifin=im
       endif
    end if
    if (allocated(atprt%xp)) then
       do i=ideb,im
          write(6,*)'%xp= ', i,atprt%xp(:,i)
       end do
       do i=ideb,im
          write(6,*)'%ityp= ', i,atprt%ityp(i)
       end do
       do i=ideb,im
          write(6,*)'%num_at_glob= ', i,atprt%num_at_glob(i)
       end do
       if (iw==0) return
       do i=ideb,im
          write(6,*)'%fp= ', i,atprt%fp(:,i)
       end do
       do i=ideb,im
          write(6,*)'%ielat= ', i,atprt%ielat(i)
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
          do i=ideb,im
             write(6,*)'%vp= ', i,atprt%vp(:,i)
          end do
          do i=ideb,im
             write(6,*)'%xpp= ', i,atprt%xpp(:,i)
          end do
          class is (atom_config_e)
          write(6,*)'prt_e'
          do i=ideb,im
             write(6,*)'%vp= ', i,atprt%vp(:,i)
          end do
          do i=ideb,im
             write(6,*)'%xpp= ', i,atprt%xpp(:,i)
          end do


          if (atprt%lsigat) then
             do i=ideb,im
                write(6,*)'%sigat= ',i, atprt%sigat(:,:,i)
             end do
          end if
          if (atprt%lprteat) then
             do i=ideb,im
                write(6,*)'%eat= ', i,atprt%eat(i)
             end do
          end if
       end select


       if (atprt%ltabvois) then
          do i=ideb,im
             write(6,*)'%iwmax= ', i,atprt%iwmax(i)
          end do

          if (allocated(atprt%indi))then
             if (ideb==1) then
                ist=1
             else
                ist=atprt%iwmax(ideb-1)+1
             end if
             if (ifin==im) then
                ifn=size(atprt%indi)
             else
                ifn=atprt%iwmax(ifin)
             end if

             do i=ist,ifn,100
                write(6,*)'indi', i,atprt%indi(i)
             end do
          end if
       end if
    end if
    write(6,*)'out print'
    write(6,*)
  end subroutine print
!MANQUE SIG AU MINIMUM

  subroutine ndm2config (atndm,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax,indi,nvois,vp,xpp,&
       &lprteatR,eat,lsigatR,sigat,llangevinR,glangv,laxR,ax,ldeall)
    class(atom_config)::atndm
    integer,intent(in)::im,imm
    real(double),intent(inout),allocatable,dimension(:,:):: xp,fp
    integer,intent(inout),dimension(:),allocatable::ityp,ielat
    integer,optional,intent(inout),allocatable,dimension(:)::num_at_glob
    logical, optional,intent(in)::ltabvois
    integer,optional,intent(in)::nvois
    integer, optional,intent(inout),allocatable ::iwmax(:)
    integer,optional,intent(inout),allocatable:: indi(:)
    real(double),optional,intent(inout),dimension(:,:),allocatable:: vp,xpp
    real(double),optional,intent(inout),allocatable:: eat(:),sigat(:,:,:),glangv(:,:),ax(:,:)

    logical,intent(in),optional ::lprteatR,lsigatR,llangevinR,laxR
    logical ::lprteat,lsigat ,ltbv,llangevin,lax
    logical,optional,intent(in)::ldeall
    
    logical::ldealloc

    ldealloc=.false.
    if (present(ldeall))ldealloc=ldeall

    lprteat=.false.;lsigat=.false.;ltbv=.false.;llangevin=.false.;lax=.false.

!    write(6,*)'INNDM2CONF',allocated(xp),allocated(fp),allocated(vp),allocated(xpp),&
!         &allocated(ityp),allocated(ielat),allocated(num_at_glob)
    if (present(lprteatR))lprteat=lprteatR; if(present(lsigatR))lsigat=lsigatR;  if(present(ltabvois))ltbv=ltabvois
    if (present(llangevinR))llangevin=llangevinR; if(present(laxR))lax=laxR

    select type(atndm)
       !    type is (atom_config)
       class is (atom_config)
       call atndm%init(im,imm,ltabvois,nvois)
    type is (atom_config_e)
       call atndm%init(im,imm,ltabvois,nvois,lsigat,lprteat,llangevin,lax)
    end select

    atndm%icaltabt=0
    atndm%xp(:,1:imm)=xp(:,1:imm)
    atndm%fp(:,1:imm)=fp(:,1:imm)
    atndm%ityp(1:imm)=ityp(1:imm)
    atndm%ielat(1:im)=ielat(1:im)
    if (ldealloc) deallocate(xp,fp,ityp,ielat)
    if (present(num_at_glob))then
       atndm%num_at_glob(1:imm)=num_at_glob(1:imm)
       if (ldealloc) deallocate(num_at_glob)
    end if

    if (ltbv) then
       !       write(6,*)'sizes ', size (iwmax),size(atndm%iwmax)
       atndm%ltabvois=.true.
       atndm%iwmax(1:imm)=iwmax(1:imm)
       !       allocate(atndm%indi(nvois))
       atndm%indi(1:nvois)=indi(1:nvois)
       if (ldealloc) deallocate(iwmax,indi)
    end if

    select type(atndm)
    type is (atom_config_d)
       if (present(vp))then
          atndm%vp(:,1:imm)=vp(:,1:imm)
          if (ldealloc) deallocate(vp)
       end if
       if(present(xpp))then
          atndm%xpp(:,1:imm)=xpp(:,1:imm)
          if (ldealloc) deallocate(xpp)
       end if
    type is (atom_config_e)
       if (present(sigat).and.(atndm%lsigat))then
          !  allocate(sigat(3,3,imm))
          atndm%sigat(:,:,1:atndm%imm)=sigat(:,:,1:atndm%imm)
          if (ldealloc) deallocate(sigat)
       end if
       if (present(eat).and.(atndm%lprteat))then
          !          allocate(eat(imm))
          atndm%eat(1:atndm%imm)=eat(1:atndm%imm)
          if (ldealloc) deallocate(eat)
       end if
       if (present(glangv).and.(atndm%llangevin))then
          atndm%glangv(1:3,1:atndm%imm)=glangv(1:3,1:atndm%imm)
          if (ldealloc) deallocate(glangv)
       end if
       if (present(ax).and.(atndm%lax))then
          atndm%ax(1:3,1:atndm%imm)=ax(1:3,1:atndm%imm)
          if (ldealloc) deallocate(ax)
       end if
    end select

  end subroutine ndm2config

  subroutine config2ndm (atndm,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax,indi,vp,xpp,eat,sigat,ax,ldeall)
    class(atom_config),intent(inout)::atndm
    integer,intent(out)::im
    integer,intent(out)::imm
    real(double),intent(out),allocatable:: xp(:,:),fp(:,:)

    integer,intent(out),dimension(:),allocatable::ityp,num_at_glob,ielat
    logical, intent(out)::ltabvois
    integer, optional,intent(out),allocatable ::iwmax(:)
    integer,optional,intent(out),allocatable:: indi(:)

    real(double),optional,intent(inout),allocatable:: vp(:,:),xpp(:,:)

    real(double),optional,intent(inout),allocatable::eat(:),sigat(:,:,:),ax(:,:)
    logical,optional,intent(in)::ldeall
    
    logical::ldealloc
    integer::is
    ldealloc=.false.
    if (present(ldeall))ldealloc=ldeall
    
    !    integer,intent(in)::nvois

    imm=atndm%imm
    im=atndm%im
!    write(6,*)'INCONF2NDM',allocated(xp),allocated(fp),allocated(vp),allocated(xpp),&
!         &allocated(ityp),allocated(ielat),allocated(num_at_glob)
    if (.not.(allocated(xp)))then
       allocate(xp(3,imm));allocate(fp(3,imm));allocate(vp(3,imm));allocate(xpp(3,imm))
       allocate(ityp(imm));allocate(ielat(imm));allocate(num_at_glob(imm))
    end if
    xp(:,1:imm)=atndm%xp(:,1:imm)
    fp(:,1:imm)=atndm%fp(:,1:imm)
    ityp(1:imm)=atndm%ityp(1:imm)
    num_at_glob(1:imm)=atndm%num_at_glob(1:imm)
    ielat(1:imm)=atndm%ielat(1:imm)
    ltabvois=atndm%ltabvois
    if (atndm%ltabvois) then
       !       write(6,*)'sizes ', size (iwmax),size(atndm%iwmax)
!      allocate(iwmax(imm))
       if (.not.allocated(iwmax))allocate(iwmax(imm))
       iwmax(1:imm)=atndm%iwmax(1:imm)
       is =size(atndm%indi)
       !       write(6,*)'IS',is
      ! allocate(indi(is))
       if (.not.allocated(indi))allocate(indi(is))
       indi(1:is)=atndm%indi(1:is)
    end if
    select type(atndm)
    type is (atom_config_d)
       if(present(vp))vp(:,1:imm)=atndm%vp(:,1:imm)
       if(present(xpp))xpp(:,1:imm)=atndm%xpp(:,1:imm)
    type is (atom_config_e)
       if (present(sigat).and.(atndm%lsigat))then
          if (.not.(allocated(sigat)) )  allocate(sigat(3,3,imm))
          sigat(:,:,1:atndm%imm)=atndm%sigat(:,:,1:atndm%imm)
       end if
       if (present(eat).and.(atndm%lprteat))then
          if(.not.(allocated(eat)))   allocate(eat(imm))
          eat=0
          eat(1:atndm%im)=atndm%eat(1:atndm%im)
       end if
       if (present(ax).and.(atndm%lax))then
          if(.not.(allocated(ax)))   allocate(ax(3,imm))
          ax=0
          ax(:,1:atndm%im)=atndm%ax(:,1:atndm%im)
       end if
       if (ldealloc)call atndm%dealloc

    end select
  end subroutine config2ndm



end module atomconfig


