module atomconfig
  USE T_kind_param_m,only:double,long
  USE Mat_utils_mod,only: fillbuffer3D,fillbuffer1D,fillbuffer9D
#ifdef PARA
  use mpi
  USE Tpara,only:NDM_MPI_REAL_DOUBLE,ierr
  use gen_com_m,only:rang


#endif
  use paraconfig,only:para_config  

  implicit none
#ifdef PARA
  integer,dimension(MPI_STATUS_SIZE):: status  ! statut de la communication
#endif

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
     real(double)::rvois ! rayon des voisins
     integer,allocatable:: indi(:) ! indice de tous les voisins
     logical, allocatable:: lgul(:)
     integer,allocatable::num_at_glob(:)
#ifdef PARA
     integer,allocatable::proc_at(:) ! tableau de taille im_glob total indiquant le numéro du proc qui gère l'atome
#endif     
   contains
     procedure, pass::init=>init_atom_config
     procedure, pass::copy_atom=>copy_atom_b
     procedure, pass::dealloc=>dealloc_atom_config
     procedure, pass::vers_master=>vers_master_atom !(atcfloc,atcfcomp,div)
     procedure, pass::master2loc=>master2loc_atom !(atcfcomp,atcfloc,div)
     procedure, pass::copy_config
     procedure, pass::print
     procedure, pass::pack
     procedure, pass::fab
     procedure, pass::add2conf
     procedure, pass::extend

     procedure, pass::send2proc=>s2p_atom
     procedure, pass::send2all=>s2a_atom
     procedure, pass::recv=>rcv_atom
     procedure, pass::zero=>zero_atom
     !
  end type atom_config

  type, extends (atom_config):: atom_config_d ! type dynamique des configurations atomiques(+vp/+xpp). vp et xpp seront toujours allouées
     real(double),allocatable ::xpp(:,:)
     real(double),allocatable::vp(:,:)
   contains
     procedure, pass::copy_atom=>copy_atom_d
     procedure, pass::dealloc=>dealloc_atom_config_d
     procedure, pass::switch_atom
     !#ifdef PARA
     procedure, pass::send2proc=>s2p_atom_d
     procedure, pass::send2all=>s2a_atom_d
     procedure, pass::recv=>rcv_atom_d
     procedure, pass::zero=>zero_atom_d
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
     procedure, pass::send2proc=>s2p_atom_e
     procedure, pass::send2all=>s2a_atom_e
     procedure, pass::recv=>rcv_atom_e
     procedure, pass::zero=>zero_atom_e
   
     !#endif     
  end type atom_config_e
  
contains
  !initialisations
  
  subroutine init_atom_config(atconf,imin,immin,ltabvois,nvois,rvois,lsigat,lprteat,lLangevin,lax,lreallocate)
    class(atom_config),intent(inout)::atconf
    integer,intent(in):: imin
    logical,optional, intent(in)::ltabvois,lsigat,lprteat,lLangevin,lax,lreallocate
    integer, optional::nvois,immin
    real(double),optional::rvois
    integer::nv
    logical ::ltbv,lrealloc
    real(double)::rv
    nv=0 ;  if(present(nvois))nv=nvois
    rv=0;  if(present(rvois))rv=rvois
    lrealloc=.false.
    if (present(lreallocate))lrealloc=lreallocate
    ltbv=.false.
    atconf%im=imin       
    if (present(immin))then
       atconf%imm=immin
    else
       atconf%imm=imin
    end if
    if(lrealloc) then
       if (allocated(atconf%xp))then
          deallocate(atconf%xp);deallocate(atconf%fp);deallocate(atconf%ityp)
          deallocate(atconf%ielat);deallocate(atconf%lgul);deallocate(atconf%num_at_glob)
#ifdef PARA
          deallocate(atconf%proc_at)
#endif    
       end if
    end if
    if (.not.allocated(atconf%xp))then
       allocate(atconf%xp(3,atconf%imm));allocate(atconf%fp(3,atconf%imm));allocate(atconf%ityp(atconf%imm))
       allocate(atconf%ielat(atconf%imm));allocate(atconf%lgul(atconf%imm));allocate(atconf%num_at_glob(atconf%imm))
#ifdef PARA
       allocate(atconf%proc_at(atconf%imm))
#endif    
       
    end if
    atconf%ityp=0;atconf%xp=0;atconf%fp=0;atconf%ielat=0; atconf%lgul=.false.;atconf%num_at_glob=0
#ifdef PARA
    atconf%proc_at=-1
#endif    
    if(present(ltabvois)) then
       ltbv=ltabvois
    endif
    if(ltbv)then
       atconf%ltabvois=.true.
       atconf%rvois=rv
       if ((lrealloc).and.(allocated(atconf%iwmax)))deallocate(atconf%iwmax)
       if (.not.(allocated(atconf%iwmax)))then
          allocate(atconf%iwmax(atconf%imm)); atconf%iwmax=0
       endif
       if(.not.allocated(atconf%indi))then
          atconf%nvois=nv
          if (nv.ne.0)then
             allocate(atconf%indi(nv))
             atconf%indi=0
          end if
       end if
    else
       atconf%ltabvois=.false.
    end if

    select type (atconf)
    class is (atom_config_d)
       if ((lrealloc).and.(allocated(atconf%vp)))then
          deallocate(atconf%vp); deallocate(atconf%xpp)
       end if
       if (.not.allocated(atconf%vp))then
          allocate(atconf%vp(3,atconf%imm));allocate(atconf%xpp(3,atconf%imm))
       end if
       atconf%vp=0;atconf%xpp=0
    type is (atom_config_e)
!       if ((lrealloc).and.(allocated(atconf%vp)))then
!          deallocate(atconf%vp); deallocate(atconf%xpp)
!       end if
!       if (.not.allocated(atconf%vp))then
!          allocate(atconf%vp(3,atconf%imm));allocate(atconf%xpp(3,atconf%imm))
!       end if
!       atconf%vp=0;atconf%xpp=0
       atconf%lprteat=.false.
       atconf%lsigat=.false.
       atconf%lLangevin=.false.
       atconf%lax=.false.
       if(present(lprteat)) atconf%lprteat=lprteat
       if(present(lsigat)) atconf%lsigat=lsigat
       if(present(lLangevin)) atconf%lLangevin=lLangevin
       if(present(lax)) atconf%lax=lax
       if(atconf%lprteat)then
          if ((lrealloc).and.(allocated(atconf%eat)))deallocate(atconf%eat)
          if (.not.allocated(atconf%eat))allocate(atconf%eat(atconf%imm))
          atconf%eat=0
       end if
       if(atconf%lsigat)then
          if ((lrealloc).and.(allocated(atconf%sigat)))deallocate(atconf%sigat)
          if (.not.allocated(atconf%sigat)) allocate(atconf%sigat(3,3,atconf%imm))
          atconf%sigat=0
       end if
       if(atconf%llangevin)then
          if ((lrealloc).and.(allocated(atconf%Glangv)))deallocate(atconf%glangv)
          if (.not.allocated(atconf%Glangv))allocate(atconf%Glangv(3,atconf%imm))
          atconf%Glangv=0
       end if
       if(atconf%lax)then
          if ((lrealloc).and.(allocated(atconf%ax)))deallocate(atconf%ax)
          if (.not.allocated(atconf%ax))allocate(atconf%ax(3,atconf%imm))
          atconf%ax=0
       end if


    end select
    atconf%icaltabt=0
  end subroutine init_atom_config

  !copie d'un élément
  
  subroutine copy_atom_b(atsource,i,atcible,j,lextend)
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
#ifdef PARA
    atcible%proc_at(j)=atsource%proc_at(i)
#endif    
    atcible%ielat(j)=atsource%ielat(i)
    atcible%lgul(j)=atsource%lgul(i)
    if ((atcible%ltabvois).and.(atsource%ltabvois))then
       if (i.gt.1) then
          nvj=atsource%iwmax(i)-atsource%iwmax(i-1)
       else
          nvj=atsource%iwmax(i)
       end if
       if (j.gt.1) then
          atcible%iwmax(j)=nvj+atsource%iwmax(j-1)
       else
          atcible%iwmax(j)=nvj
       end if
       do iw=1,nvj
          if (i.gt.1) then
             idwi=iw+atsource%iwmax(i-1)
          else
             idwi=iw
          end if
          if (j.gt.1) then
             idwj=iw+atcible%iwmax(j-1)
          else
             idwj=iw
          end if
             
             atcible%indi(idwj)=atsource%indi(idwi)
       end do
    end if

  end subroutine copy_atom_b

  subroutine extend(atcf,iadd)
    class(atom_config),intent(inout)::atcf
    integer,intent(in)::iadd

    type(atom_config):: attemp
    type(atom_config_d):: attemp_d
    type(atom_config_e):: attemp_e
    integer::imcn,nvois,immcn
    real(double)::rvois
    immcn=atcf%imm+iadd
    imcn=atcf%im+iadd
    if (atcf%ltabvois)then
       nvois=int(1.1*imcn/atcf%im)*size(atcf%indi)
       rvois=atcf%rvois
    else
       nvois=0
       rvois=0
    end if
    select type(atcf)
    type is (atom_config)
       call attemp%init(imcn,immcn,atcf%ltabvois,nvois,rvois)
    type is (atom_config_d)
       call attemp_d%init(imcn,immcn,atcf%ltabvois,nvois,rvois)
    type is (atom_config_e)
       call attemp_e%init(imcn,immcn,atcf%ltabvois,nvois,rvois,atcf%lsigat,atcf%lprteat, atcf%llangevin,atcf%lax)
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
    call copy_atom_b(atsource,i,atcible,j,let)
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

  subroutine s2p_atom (atcf, rgcib,comm,caracT)
    class(atom_config):: atcf
    integer,intent(in)::rgcib,comm
    integer:: size1,size3,sizeV,ierr
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    
!x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   
    size1=atcf%imm;size3=3*size1
#ifdef PARA
    if (.not.present(caracT)) then
       carac='xfniewdlp'
    else
       carac=caracT
    end if
    if(scan('n',carac).ne.0)    call MPI_SEND(atcf%num_at_glob, size1, MPI_INTEGER, rgcib,104,comm,ierr)
    if(scan('i',carac).ne.0)    call MPI_SEND(atcf%ityp, size1, MPI_INTEGER, rgcib,105,comm,ierr)
    if(scan('e',carac).ne.0)    call MPI_SEND(atcf%ielat, size1, MPI_INTEGER, rgcib,106,comm,ierr)
    if(scan('p',carac).ne.0)    call MPI_SEND(atcf%proc_at, size1, MPI_INTEGER, rgcib,102,comm,ierr)
    if(scan('l',carac).ne.0)    call MPI_SEND(atcf%lgul, size1, MPI_LOGICAL, rgcib,103,comm,ierr)
    if(scan('x',carac).ne.0)   then
       call MPI_SEND(atcf%xp, size3, NDM_MPI_REAL_DOUBLE, rgcib,100,comm,ierr)
    end if
    if(scan('f',carac).ne.0)    call MPI_SEND(atcf%fp, size3, NDM_MPI_REAL_DOUBLE, rgcib,101,comm,ierr)
    if (atcf%ltabvois) then
       sizeV=size(atcf%indi)
       if(scan('w',carac).ne.0)    call MPI_SEND(atcf%iwmax, size1, MPI_INTEGER, rgcib,107,comm,ierr)
       if(scan('d',carac).ne.0)    call MPI_SEND(atcf%indi, sizeV, MPI_INTEGER, rgcib,108,comm,ierr)
    end if
#endif
  end subroutine s2p_atom

  subroutine s2p_atom_d (atcf, rgcib,comm,caracT)
    class(atom_config_d):: atcf
    integer,intent(in)::rgcib,comm
    integer:: size1,size3,sizeV,ierr
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at
!voir au dessus + v=vp,r=xpp
    

    size1=atcf%imm;size3=3*size1
#ifdef PARA
    if (.not.present(caracT)) then
       carac='xfniewdlpvr'
    else
       carac=caracT
    end if
    call s2p_atom(atcf,rgcib,comm,carac)
    if(scan('v',carac).ne.0)    call MPI_SEND(atcf%vp, size3, NDM_MPI_REAL_DOUBLE, rgcib,109,comm,ierr)
    if(scan('r',carac).ne.0)    call MPI_SEND(atcf%xpp, size3, NDM_MPI_REAL_DOUBLE, rgcib,110,comm,ierr)
#endif
  end subroutine s2p_atom_d

    subroutine s2p_atom_e (atcf, rgcib,comm,caracT)
    class(atom_config_e):: atcf
    integer,intent(in)::rgcib,comm
    integer:: size1,size3,sizeV,ierr
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at
!voir au dessus + v=vp,r=xpp
!voir au dessus + u=eat,g=glangv;a=ax;s=sigat
    

    size1=atcf%imm;size3=3*size1
#ifdef PARA
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugas'
    else
       carac=caracT
    end if
    call s2p_atom_d(atcf,rgcib,comm,carac)
    if (atcf%lprteat)then
    if(scan('u',carac).ne.0) call MPI_SEND(atcf%eat, size1, NDM_MPI_REAL_DOUBLE, rgcib,111,comm,ierr)       
    end if
    if (atcf%llangevin)then
    if(scan('g',carac).ne.0)  call MPI_SEND(atcf%glangv, size3, NDM_MPI_REAL_DOUBLE, rgcib,112,comm,ierr)     
    end if
    if (atcf%lax)then
     if(scan('a',carac).ne.0)call MPI_SEND(atcf%ax, size3, NDM_MPI_REAL_DOUBLE, rgcib,113,comm,ierr)           
    end if
    if (atcf%lsigat)then
     if(scan('s',carac).ne.0) call MPI_SEND(atcf%sigat, 3*size3, NDM_MPI_REAL_DOUBLE, rgcib,114,comm,ierr)          
    end if
#endif
  end subroutine s2p_atom_e


  subroutine zero_atom (atcf)
    class(atom_config):: atcf

    atcf%num_at_glob=0
    atcf%ityp=0
    atcf%ielat=0
#ifdef PARA
    atcf%proc_at=-1
#endif
    atcf%lgul=.false.
    atcf%xp=0.
    atcf%fp=0.
    
    if (atcf%ltabvois) then
       atcf%iwmax=0
       atcf%indi=0
    end if

  end subroutine zero_atom
  subroutine zero_atom_d (atcf)
    class(atom_config_d):: atcf

    call zero_atom(atcf)
    atcf%vp=0.
    atcf%xpp=0.
  end subroutine zero_atom_d

  subroutine zero_atom_e (atcf)
    class(atom_config_e):: atcf

    call zero_atom_d(atcf)

    if (atcf%lprteat)then

       atcf%eat=0
    end if
    if (atcf%llangevin)then
       atcf%glangv=0
    end if
    if (atcf%lax)then
       atcf%ax=0
    end if
    if (atcf%lsigat)then
       atcf%sigat=0
    end if


  end subroutine zero_atom_e
  


  subroutine rcv_atom (atcf, rgem,comm,caracT)
    class(atom_config):: atcf
    integer,intent(in)::rgem,comm
    integer:: size1,size3,sizeV,ierr
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    
!x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   
    size1=atcf%imm;size3=3*size1
#ifdef PARA
    if (.not.present(caracT)) then
       carac='xfniewdlp'
    else
       carac=caracT
    end if
    if(scan('n',carac).ne.0)    call MPI_RECV(atcf%num_at_glob, size1, MPI_INTEGER, rgem,104,comm,status,ierr)
    if(scan('i',carac).ne.0)    call MPI_RECV(atcf%ityp, size1, MPI_INTEGER, rgem,105,comm,status,ierr)
    if(scan('e',carac).ne.0)    call MPI_RECV(atcf%ielat, size1, MPI_INTEGER, rgem,106,comm,status,ierr)
    if(scan('p',carac).ne.0)    call MPI_RECV(atcf%proc_at, size1, MPI_INTEGER, rgem,102,comm,status,ierr)
    if(scan('l',carac).ne.0)    call MPI_RECV(atcf%lgul, size1, MPI_LOGICAL, rgem,103,comm,status,ierr)
    if(scan('x',carac).ne.0)    then
       call MPI_RECV(atcf%xp, size3, NDM_MPI_REAL_DOUBLE, rgem,100,comm,status,ierr)
    end if
    if(scan('f',carac).ne.0)    call MPI_RECV(atcf%fp, size3, NDM_MPI_REAL_DOUBLE, rgem,101,comm,status,ierr)
    if (atcf%ltabvois) then
       sizeV=size(atcf%indi)
       if(scan('w',carac).ne.0)    call MPI_RECV(atcf%iwmax, size1, MPI_INTEGER, rgem,107,comm,status,ierr)
       if(scan('d',carac).ne.0)    call MPI_RECV(atcf%indi, sizeV, MPI_INTEGER, rgem,108,comm,status,ierr)
    end if
#endif
  end subroutine rcv_atom

  subroutine rcv_atom_d (atcf, rgem,comm,caracT)
    class(atom_config_d):: atcf
    integer,intent(in)::rgem,comm
    integer:: size1,size3,sizeV,ierr
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at
!voir au dessus + v=vp,r=xpp
    

    size1=atcf%imm;size3=3*size1
#ifdef PARA
    if (.not.present(caracT)) then
       carac='xfniewdlpvr'
    else
       carac=caracT
    end if
    call rcv_atom(atcf,rgem,comm,carac)
    if(scan('v',carac).ne.0)    call MPI_RECV(atcf%vp, size3, NDM_MPI_REAL_DOUBLE, rgem,109,comm,status,ierr)
    if(scan('r',carac).ne.0)    call MPI_RECV(atcf%xpp, size3, NDM_MPI_REAL_DOUBLE, rgem,110,comm,status,ierr)
#endif
  end subroutine rcv_atom_d

    subroutine rcv_atom_e (atcf, rgem,comm,caracT)
    class(atom_config_e):: atcf
    integer,intent(in)::rgem,comm
    integer:: size1,size3,sizeV,ierr
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at
!voir au dessus + v=vp,r=xpp
!voir au dessus + u=eat,g=glangv;a=ax;s=sigat
    

    size1=atcf%imm;size3=3*size1
#ifdef PARA
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugas'
    else
       carac=caracT
    end if
    call rcv_atom_d(atcf,rgem,comm,carac)
    if (atcf%lprteat)then
    if(scan('u',carac).ne.0) call MPI_RECV(atcf%eat, size1, NDM_MPI_REAL_DOUBLE, rgem,111,comm,status,ierr)       
    end if
    if (atcf%llangevin)then
    if(scan('g',carac).ne.0)  call MPI_RECV(atcf%glangv, size3, NDM_MPI_REAL_DOUBLE, rgem,112,comm,status,ierr)     
    end if
    if (atcf%lax)then
     if(scan('a',carac).ne.0)call MPI_RECV(atcf%ax, size3, NDM_MPI_REAL_DOUBLE, rgem,113,comm,status,ierr)           
    end if
    if (atcf%lsigat)then
     if(scan('s',carac).ne.0) call MPI_RECV(atcf%sigat, 3*size3, NDM_MPI_REAL_DOUBLE, rgem,114,comm,status,ierr)          
    end if
#endif
  end subroutine rcv_atom_e

  
  subroutine s2a_atom(atcf,rgemet,comm,caracT)
    class(atom_config)::atcf
    integer,intent(in)::rgemet,comm
    integer:: size1,size3,sizeV
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at
    
    size1=atcf%imm;size3=3*size1
#ifdef PARA
    if (.not.present(caracT)) then
       carac='xfniewdlp'
    else
       carac=caracT
    end if
    if(scan('n',carac).ne.0) call MPI_BCAST(atcf%num_at_glob, size1,MPI_INTEGER, rgemet,comm,ierr)
    if(scan('i',carac).ne.0) call MPI_BCAST(atcf%ityp, size1,MPI_INTEGER, rgemet,comm,ierr)
    if(scan('e',carac).ne.0)call MPI_BCAST(atcf%ielat, size1,MPI_INTEGER, rgemet,comm,ierr)
    if(scan('p',carac).ne.0)call MPI_BCAST(atcf%proc_at, size1,MPI_INTEGER, rgemet,comm,ierr)
    if(scan('l',carac).ne.0)call MPI_BCAST(atcf%lgul, size1,MPI_LOGICAl, rgemet,comm,ierr)
    if(scan('x',carac).ne.0)then
       call MPI_BCAST(atcf%xp, size3, NDM_MPI_REAL_DOUBLE, rgemet,comm,ierr)
    end if
    if(scan('f',carac).ne.0)call MPI_BCAST(atcf%fp, size3, NDM_MPI_REAL_DOUBLE, rgemet,comm,ierr)
    if (atcf%ltabvois) then
       sizeV=size(atcf%indi)
       if(scan('w',carac).ne.0)call MPI_BCAST(atcf%iwmax, size1,MPI_INTEGER, rgemet,comm,ierr)
       if(scan('d',carac).ne.0)call MPI_BCAST(atcf%indi, sizeV,MPI_INTEGER, rgemet,comm,ierr)
    end if
#endif
  end subroutine s2a_atom
  
  subroutine s2a_atom_d(atcf,rgemet,comm,caracT)
    class(atom_config_d)::atcf
    integer,intent(in)::rgemet,comm
    integer:: size1,size3,sizeV
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at
!voir au dessus + v=vp,r=xpp
    if (.not.present(caracT)) then
       carac='xfniewdlpvr'
    else
       carac=caracT
    end if
    
    call s2a_atom(atcf,rgemet,comm,carac)
    size1=atcf%imm;size3=3*size1
#ifdef PARA
    if(scan('r',carac).ne.0)    call MPI_BCAST(atcf%xpp, size3, NDM_MPI_REAL_DOUBLE, rgemet,comm,ierr)
    if(scan('v',carac).ne.0)    call MPI_BCAST(atcf%vp, size3, NDM_MPI_REAL_DOUBLE, rgemet,comm,ierr)
#endif
  end subroutine s2a_atom_d

  
  subroutine s2a_atom_e(atcf,rgemet,comm,caracT)
    class(atom_config_e)::atcf
    integer,intent(in)::rgemet,comm
    integer:: size1,size3,sizeV
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at
!voir au dessus + v=vp,r=xpp
!voir au dessus + u=eat,g=glangv;a=ax;s=sigat
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugas'
    else
       carac=caracT
    end if

    call s2a_atom_d(atcf,rgemet,comm,carac)
    size1=atcf%imm;size3=3*size1
#ifdef PARA
    if (atcf%lprteat)then
    if(scan('u',carac).ne.0)       call MPI_BCAST(atcf%eat, size1, NDM_MPI_REAL_DOUBLE, rgemet,comm,ierr)
    end if
    if (atcf%llangevin)then
    if(scan('g',carac).ne.0)       call MPI_BCAST(atcf%Glangv, size3, NDM_MPI_REAL_DOUBLE, rgemet,comm,ierr)
    end if
    if (atcf%lax)then
     if(scan('a',carac).ne.0)      call MPI_BCAST(atcf%ax, size3, NDM_MPI_REAL_DOUBLE, rgemet,comm,ierr)
    end if
    if (atcf%lsigat)then
     if(scan('s',carac).ne.0)      call MPI_BCAST(atcf%sigat, 3*size3, NDM_MPI_REAL_DOUBLE, rgemet,comm,ierr)
    end if
#endif
  end subroutine s2a_atom_e

    
  ! copie d'une config entière vers config de base
  subroutine copy_config (atsource,atcible,lrescl)
    class(atom_config),intent(in)::atsource
    class(atom_config)::atcible
    logical,intent(in)::lrescl
    integer::nvois
    real(double)::rvois
    logical :: lstop
    if (lrescl) then
!       if (atcible%imm.ne.atsource%imm) then
       call atcible%dealloc
       if (atsource%ltabvois)then
          nvois=atsource%nvois; rvois=atsource%rvois
       else
          nvois=0;rvois=0
       end if
          call atcible%init(atsource%im,atsource%imm,atcible%ltabvois,nvois,rvois)
          atcible%icaltabt=atsource%icaltabt    

    else
       lstop=.false.
       if ((atcible%im.lt.atsource%im).or.(atcible%imm.lt.atsource%imm)) lstop=.true.

       if (atcible%ltabvois) then
          if (size(atcible%indi).lt.size(atsource%indi)) lstop=.true.
          if (atcible%nvois.lt.atsource%nvois)lstop=.true.
       end if
       if (lstop)then
          write(6,*)'(atcible%im < atsource%im  ou size(atcible%indi)<size(atsource%indi) )et lrescl = false ; pas possible'
          stop
       end if
    end if
    atcible%xp(1:3,1:atsource%imm)=atsource%xp(1:3,1:atsource%imm)
    atcible%fp(:,1:atsource%imm)=atsource%fp(:,1:atsource%imm)
    atcible%ielat(1:atsource%imm)=atsource%ielat(1:atsource%imm)
    atcible%lgul(1:atsource%imm)=atsource%lgul(1:atsource%imm)
    atcible%ityp(1:atsource%imm)=atsource%ityp(1:atsource%imm)
    atcible%num_at_glob(1:atsource%imm)=atsource%num_at_glob(1:atsource%imm)
#ifdef PARA
    atcible%proc_at(1:atsource%imm)=atsource%proc_at(1:atsource%imm)
#endif    
    
    if ((atsource%ltabvois).and.(atcible%ltabvois)) then
       atcible%iwmax(1:atsource%imm)=atsource%iwmax(1:atsource%imm)
       atcible%indi(1:atsource%nvois)=atsource%indi(1:atsource%nvois)
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
#ifdef PARA
       deallocate(atconf%proc_at)
#endif    
       
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



  !on inverse deux atomes dans la configuration
  subroutine switch_atom(atsource,ind_switch_1, ind_switch_2)
    class(atom_config_d)::atsource
    type(atom_config_d):: intermediaire
    integer :: ind_switch_1, ind_switch_2,nag1,nag2
    logical :: lex = .true.

    call intermediaire%init(imin=2)
    nag1=atsource%num_at_glob(ind_switch_1)
    nag2=atsource%num_at_glob(ind_switch_2)
    call atsource%copy_atom(ind_switch_1,intermediaire,1)
    call atsource%copy_atom(ind_switch_2,intermediaire,2)
    call intermediaire%copy_atom(2,atsource,ind_switch_1)
    call intermediaire%copy_atom(1,atsource,ind_switch_2)
    atsource%num_at_glob(ind_switch_1)=nag1
    atsource%num_at_glob(ind_switch_2)=nag2
    
!!$        select type (atsource)
!!$    type is (atom_config_d)
!!$     call atsource%copy_config(intermediaire, lex)
!!$     atsource%xp(:,ind_switch_1)=intermediaire%xp(:,2)
!!$     atsource%fp(:,ind_switch_1)=intermediaire%fp(:,2)
!!$     atsource%ielat(ind_switch_1)=intermediaire%ielat(ind_switch_2)
!!$     atsource%lgul(ind_switch_1)=intermediaire%lgul(ind_switch_2)
!!$     atsource%ityp(ind_switch_1)=intermediaire%ityp(ind_switch_2)
!!$     atsource%num_at_glob(ind_switch_1)=intermediaire%num_at_glob(ind_switch_2)
!!$     ! faut il changer des trucs concernant ltabvois ? iwmax ? ou indi ?
!!$     atsource%vp(:,ind_switch_1)=intermediaire%vp(:,ind_switch_2)
!!$     atsource%xpp(:,ind_switch_1)=intermediaire%xpp(:,ind_switch_2)
!!$
!!$     atsource%xp(:,ind_switch_2)=intermediaire%xp(:,ind_switch_1)
!!$     atsource%fp(:,ind_switch_2)=intermediaire%fp(:,ind_switch_1)
!!$     atsource%ielat(ind_switch_2)=intermediaire%ielat(ind_switch_1)
!!$     atsource%lgul(ind_switch_2)=intermediaire%lgul(ind_switch_1)
!!$     atsource%ityp(ind_switch_2)=intermediaire%ityp(ind_switch_1)
!!$     atsource%num_at_glob(ind_switch_2)=intermediaire%num_at_glob(ind_switch_1)
!!$     ! faut il changer des trucs concernant ltabvois ? iwmax ? ou indi ?
!!$     atsource%vp(:,ind_switch_2)=intermediaire%vp(:,ind_switch_1)
!!$     atsource%xpp(:,ind_switch_2)=intermediaire%xpp(:,ind_switch_1)
!!$    end select
  end subroutine




  
#ifdef PARA  
  
#endif
  subroutine pack(at2pack,imm_in)
    class(atom_config),intent(inout):: at2pack
    integer,optional, intent(in):: imm_in
    type(atom_config)::at
    type(atom_config_d)::atd
    type(atom_config_e)::ate
    integer::i,imn,immn
    real(double)::rvois
    integer::nvois

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
    if (at2pack%ltabvois)then
       nvois=at2pack%nvois; rvois=at2pack%rvois
    else
       nvois=0;rvois=0
    end if

    select type (at2pack)
    type is (atom_config)
       call at%init(imn,immn,at2pack%ltabvois,nvois,rvois)
       do i=1,immn
          call at2pack%copy_atom(i,at,i)
       end do
       call at2pack%dealloc
       call at%copy_config(at2pack,lrescl=.true.)

    type is(atom_config_d)
       call atd%init(imn,immn,at2pack%ltabvois,nvois,rvois)
       do i=1,immn
          call at2pack%copy_atom(i,atd,i)
       end do
       call at2pack%dealloc
       call atd%copy_config(at2pack,lrescl=.true.)

    type is(atom_config_e)
       call ate%init(imn,immn,at2pack%ltabvois,nvois,rvois,at2pack%lsigat,at2pack%lprteat,llangevin=at2pack%llangevin,&
            &lax=at2pack%lax)
       do i=1,immn
          call at2pack%copy_atom(i,ate,i)
       end do
       call at2pack%dealloc
       call ate%copy_config(at2pack,lrescl=.true.)
    end select

  end subroutine pack

  subroutine fab (atsource,atcible,lrescl) ! construit atsource à partir de lgul de atcible , ecrase atcible
    class(atom_config),intent(in)::atsource
    class(atom_config),intent(out)::atcible
    logical, optional::lrescl
    logical::lrescale=.true.
    integer::i2,imtrf,i
    integer::nvois
    real(double)::rvois

    if (present(lrescl))lrescale=lrescl
    if (lrescale) then
       call atcible%dealloc 
       if (atsource%ltabvois)then
          nvois=atsource%nvois; rvois=atsource%rvois
          
       else
          nvois=0;rvois=0
       end if
       imtrf=COUNT(atsource%lgul(1:atsource%im))

       select type (atsource)
       type is (atom_config_e)
          call atcible%init(imtrf,imtrf,atsource%ltabvois,nvois,rvois,lsigat=atsource%lsigat,lprteat=atsource%lprteat,&
               &llangevin=atsource%llangevin,lax=atsource%lax)
       class is (atom_config)
          call atcible%init(imtrf,imtrf,atsource%ltabvois,nvois,rvois)
       end select
    end if
    call atcible%zero
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
    integer::nvois
    real(double)::rvois
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
    if (atsource%ltabvois)then
       nvois=int(1.1*imcib/atsource%im)*size(atsource%indi)
       rvois=atsource%rvois
    else
       nvois=0
       rvois=0
    end if

    if(atcible%imm==0) then
       select type (atsource)
       type is (atom_config_e)
          call atcible%init(imsrc,immsrc,atsource%ltabvois,nvois,rvois,lsigat=atsource%lsigat,&
               &lprteat=atsource%lprteat,llangevin=atsource%llangevin,lax=atsource%lax)
          class is (atom_config)
          call atcible%init(imsrc,immsrc,atsource%ltabvois,nvois,rvois)
       end select
    else
       select type(atcible)
       type is (atom_config)
          call atcor%init(imcib,immcib,atcible%ltabvois,nvois,rvois)
          call atcible%copy_config(atcor,.false.)
       type is (atom_config_d)
          call atcor_d%init(imcib,immcib,atcible%ltabvois,nvois,rvois)
          call atcible%copy_config(atcor_d,.false.)
       type is (atom_config_e)
          call atcor_e%init(imcib,immcib,atcible%ltabvois,nvois,rvois,atcible%lsigat,atcible%lprteat,&
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
    if (ldal) call atsource%dealloc
  end subroutine add2conf



  subroutine print(atin,i1,i2,iwr,unit,natg1,natg2,caracT)
    class(atom_config), intent(in)::atin
    integer,optional::i1,i2,iwr,unit,natg1,natg2
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    integer::i,im,ifin,ideb,ist,ifn,iw,natpr,ig,iprt,unitw
    class (atom_config),allocatable::atprt
    unitw=6
    if (.not.present(caracT)) then
       carac='xfniewdlpvfrugas'
    else
       carac=caracT
    end if
    if (present(unit))unitw=unit

    iw=1
       select type (atin)
       type is (atom_config)
          allocate(atom_config::atprt)
       type is (atom_config_d)
          allocate(atom_config_d::atprt)
       type is (atom_config_e)
          allocate(atom_config_e::atprt)
       end select

    if (present(iwr))iw=iwr
    if (present(natg1)) then
       if (present(natg2)) then
          natpr=natg2-natg1+1
          if (natpr.lT.1) stop
       else
          natg2=natg1
       end if
       iprt=0
       do ig=natg1,natg2
          do i=1,atin%im
             if(atin%num_at_glob(i)==ig) then
                iprt=iprt+1
                if (iprt.gt.natpr) then
                   write(6,*)'num_at_glob multiples ?'
                   stop
                end if
             end if
          end do
       end do
       natpr=iprt 
       call atprt%init(natpr,rvois=0.d0)
      im=atprt%im
       iprt=0
       do ig=natg1,natg2
          do i=1,atin%im
             if(atin%num_at_glob(i)==ig) then
                iprt=iprt+1
                if (iprt.gt.natpr) then
                   write(6,*)'num_at_glob multiples ?'
                   stop
                end if
                call atin%copy_atom(i,atprt,iprt)
             end if
          end do
       end do
       im=atprt%im
       iprt=0
       do ig=natg1,natg2
          do i=1,atin%im
             if(atin%num_at_glob(i)==ig) then
                iprt=iprt+1
                if (iprt.gt.natpr) then
                   write(6,*)'num_at_glob multiples ?'
                   stop
                end if
                call atin%copy_atom(i,atprt,iprt)
             end if
          end do
       end do
       im=atprt%im
    else
       call atin%copy_config(atprt,lrescl=.true.)
       im=atin%im
    end if
          

    
    write(unitw,*)'im = ',atprt%im
    write(unitw,*)'imm = ',atprt%imm
    write(unitw,*)'icaltabt = ',atprt%icaltabt
    write(unitw,*)'ltabvois ', atprt%ltabvois
!    write(6,*)
    ideb=1
    ifin=atprt%im
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
    if(scan('x',carac).ne.0)then
       do i=ideb,im
          write(unitw,*)'%xp= ', i,atprt%num_at_glob(i),atprt%xp(:,i)
       end do
    end if
    if(scan('i',carac).ne.0)then
       do i=ideb,im
          write(unitw,*)'%ityp= ', i,atprt%num_at_glob(i),atprt%ityp(i)
       end do
    end if
    if(scan('n',carac).ne.0)then
       do i=ideb,im
          write(unitw,*)'%num_at_glob= ', i,atprt%num_at_glob(i)
       end do
    end if
#ifdef PARA
    if(scan('p',carac).ne.0)then
       do i=ideb,im
          write(unitw,*)'%proc_at= ', i,atprt%num_at_glob(i),atprt%proc_at(i)
       end do
    end if
#endif    
       
       if (iw==0) return
       if(scan('f',carac).ne.0)then
          do i=ideb,im
             write(unitw,*)'%fp= ', i,atprt%num_at_glob(i),atprt%fp(:,i)
          end do
       end if
       if(scan('e',carac).ne.0)then
          do i=ideb,im
             write(unitw,*)'%ielat= ', i,atprt%num_at_glob(i),atprt%ielat(i)
          end do
       end if
       select type (atprt)
          class is (atom_config_d)
             write(unitw,*)'prt_d'
             if(scan('v',carac).ne.0)then
                do i=ideb,im
                   write(unitw,*)'%vp= ', i,atprt%num_at_glob(i),atprt%vp(:,i)
                end do
             end if
             if(scan('r',carac).ne.0)then
                do i=ideb,im
                   write(unitw,*)'%xpp= ', i,atprt%num_at_glob(i),atprt%xpp(:,i)
                end do
             end if
          class is (atom_config_e)
             write(unitw,*)'prt_e'
             if(scan('v',carac).ne.0)then
                do i=ideb,im
                   write(unitw,*)'%vp= ', i,atprt%num_at_glob(i),atprt%vp(:,i)
                end do
             end if
             if(scan('r',carac).ne.0)then
                do i=ideb,im
                   write(unitw,*)'%xpp= ', i,atprt%num_at_glob(i),atprt%xpp(:,i)
                end do
             end if


             if (atprt%lsigat) then
                if(scan('g',carac).ne.0)then
                   do i=ideb,im
                      write(unitw,*)'%sigat= ',i,atprt%num_at_glob(i), atprt%sigat(:,:,i)
                   end do
                end if
             end if
          if (atprt%lprteat) then
             if(scan('u',carac).ne.0)then
             do i=ideb,im
                write(unitw,*)'%eat= ', i,atprt%num_at_glob(i),atprt%eat(i)
             end do
          end if
          end if
       end select


       if (atprt%ltabvois) then
          do i=ideb,im
             write(unitw,*)'%iwmax= ', i,atprt%num_at_glob(i),atprt%iwmax(i)
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
                write(unitw,*)'indi', i,atprt%indi(i)
             end do
          end if
       end if
    end if
    flush(unitw)
  end subroutine print
!MANQUE SIG AU MINIMUM

  subroutine ndm2config (atndm,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax,indi,nvois,vp,xpp,&
       &lprteatR,eat,lsigatR,sigat,llangevinR,glangv,laxR,ax,ldeall,lgul)
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
    logical,allocatable,optional::lgul(:)
    logical ::lprteat,lsigat ,ltbv,llangevin,lax
    logical,optional,intent(in)::ldeall

    real(double)::rvois
    
    logical::ldealloc

    ldealloc=.false.
    if (present(ldeall))ldealloc=ldeall

    lprteat=.false.;lsigat=.false.;ltbv=.false.;llangevin=.false.;lax=.false.
    if (present(lprteatR))lprteat=lprteatR; if(present(lsigatR))lsigat=lsigatR;  if(present(ltabvois))ltbv=ltabvois
    if (present(llangevinR))llangevin=llangevinR; if(present(laxR))lax=laxR

    select type(atndm)
       !    type is (atom_config)
       class is (atom_config)
       call atndm%init(im,imm,ltabvois,nvois,rvois=0.d0)
    type is (atom_config_e)
       call atndm%init(im,imm,ltabvois,nvois,rvois=0.d0,lsigat=lsigat,lprteat=lprteat,llangevin=llangevin,lax=lax)
    end select

    atndm%icaltabt=0
    atndm%xp(:,1:imm)=xp(:,1:imm)
    atndm%fp(:,1:imm)=fp(:,1:imm)
    atndm%ityp(1:imm)=ityp(1:imm)
    atndm%ielat(1:im)=ielat(1:im)
    if (present(lgul))atndm%lgul(1:imm)=lgul(1:imm)
    if (ldealloc) deallocate(xp,fp,ityp,ielat)
    if (present(lgul))deallocate(lgul)
    if (present(num_at_glob))then
       atndm%num_at_glob(1:imm)=num_at_glob(1:imm)
       if (ldealloc) deallocate(num_at_glob)
    end if

    if (ltbv) then
       atndm%ltabvois=.true.
       atndm%iwmax(1:imm)=iwmax(1:imm)
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
!    type is (atom_config_d)
       if (present(vp))then
!          atndm%vp(:,1:imm)=vp(:,1:imm)
          if ((ldealloc).and.allocated(vp)) deallocate(vp)
       end if
       if(present(xpp))then
!          atndm%xpp(:,1:imm)=xpp(:,1:imm)
          if  ((ldealloc).and.allocated(xpp)) deallocate(xpp)
       end if

  end subroutine ndm2config

  subroutine config2ndm (atndm,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax,indi,vp,xpp,eat,sigat,ax,ldeall,lgul)
    class(atom_config),intent(inout)::atndm
    integer,intent(out)::im
    integer,intent(out)::imm
    real(double),intent(out),allocatable:: xp(:,:),fp(:,:)
    logical,allocatable,optional::lgul(:)
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
    imm=atndm%imm
    im=atndm%im
    if (.not.(allocated(xp)))then
       allocate(xp(3,imm));allocate(fp(3,imm));allocate(vp(3,imm));allocate(xpp(3,imm))
       allocate(ityp(imm));allocate(ielat(imm));allocate(num_at_glob(imm))
       if (present(lgul))allocate(lgul(imm))
    end if
    vp=0;xpp=0
    xp(:,1:imm)=atndm%xp(:,1:imm)
    fp(:,1:imm)=atndm%fp(:,1:imm)
    ityp(1:imm)=atndm%ityp(1:imm)
    if (present(lgul))lgul(1:imm)=atndm%lgul(1:imm)
    num_at_glob(1:imm)=atndm%num_at_glob(1:imm)
    ielat(1:imm)=atndm%ielat(1:imm)
    ltabvois=atndm%ltabvois
    if (atndm%ltabvois) then
       if (.not.allocated(iwmax))allocate(iwmax(imm))
       iwmax(1:imm)=atndm%iwmax(1:imm)
       is =size(atndm%indi)
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


  
  subroutine vers_master_atom(atcfloc,atcfcomp,div)
    class(atom_config),intent(in)::atcfloc
    class(atom_config)::atcfcomp
    type(para_config)::div
#ifdef PARA
    integer::idmaster,idloc,icomm ! proc master,proclocal ,communicateur


    real(double),allocatable::buffer(:,:),buffer9(:,:,:),buffer1(:)
    integer,allocatable:: ibuffer(:)
    integer,allocatable::nag(:)
    logical,allocatable::lbuffer(:)
    integer::imloc,imloc3,iproc,imrecv,icomp,imrecv3,imrecv9,imloc9
    integer::imcomp,imtot,proc_source,npim,iloc
    idmaster=0
    idloc=div%mpi_image%rank
    npim=div%mpi_image%nproc
    icomm=div%mpi_image%comm
    if (idloc==idmaster) then
       !       allocate(buffer(3,atloc%im));allocate(ibuffer(atloc%im));allocate(lbuffer(atloc%im))
       imtot=0
!       ifin=0
       do iproc=0,npim-1
!          write(6,*)'IPROC',iproc
          ! Pour le processeur maitre il n'y a rien a faire
          ! reception des donnees des autres processeurs
          if (iproc==idmaster) then
             imrecv=atcfloc%im
             imtot=imtot+imrecv
             allocate(nag(imrecv))
             nag(1:imrecv)=atcfloc%num_at_glob(1:imrecv)
             do iloc=1,imrecv
!                write(6,*)'L2M',iproc,iloc,nag(iloc)
                atcfcomp%num_at_glob(nag(iloc))=nag(iloc)
                atcfcomp%xp(1:3,nag(iloc))=atcfloc%xp(1:3,iloc)
                atcfcomp%fp(1:3,nag(iloc))=atcfloc%fp(1:3,iloc)
                atcfcomp%ityp(nag(iloc))=atcfloc%ityp(iloc)
                if( allocated(atcfloc%lgul))  atcfcomp%lgul(nag(iloc))=atcfloc%lgul(iloc)
                atcfcomp%proc_at(nag(iloc))=idmaster 
                select type(atcfloc)
                class is (atom_config_d)
                   select type (atcfcomp)
                   class is (atom_config_d)
                      atcfcomp%vp(1:3,nag(iloc))=atcfloc%vp(1:3,iloc)
                      atcfcomp%xpp(1:3,nag(iloc))=atcfloc%xpp(1:3,iloc)
                   end select
                type is (atom_config_e)
                   select type (atcfcomp)
                   class is (atom_config_e)
                      if((atcfcomp%lsigat).and.(atcfloc%lsigat))then
                         atcfcomp%sigat(1:3,1:3,nag(iloc))=atcfloc%sigat(1:3,1:3,iloc)
                      endif
                      if((atcfcomp%lprteat).and.(atcfloc%lprteat))then
                         atcfcomp%eat(nag(iloc))=atcfloc%eat(iloc)
                      endif
                      if((atcfcomp%llangevin).and.(atcfloc%llangevin))then
                         atcfcomp%glangv(:,nag(iloc))=atcfloc%glangv(:,iloc)
                      endif
                   end select
                end select
             end do
             deallocate(nag)

          else
             if (allocated(buffer)) then
                deallocate(buffer);deallocate(ibuffer);deallocate(lbuffer);deallocate(buffer1);deallocate(buffer9)
             end if
             call MPI_RECV(imrecv,1, MPI_INTEGER,  MPI_ANY_SOURCE, 10001, icomm, status, ierr)
             imtot=imtot+imrecv
             proc_source = status(MPI_SOURCE)
             imrecv3=3*imrecv; imrecv9=3*imrecv3
             allocate(nag(imrecv))
             allocate(buffer(3,imrecv))
             allocate(buffer9(3,3,imrecv))
             allocate(buffer1(imrecv))
             allocate(ibuffer(imrecv))
             allocate(lbuffer(imrecv))

             call MPI_RECV(ibuffer(1:imrecv),imrecv, MPI_INTEGER, proc_source, 10004, icomm, status, ierr)
             nag(1:imrecv)=ibuffer(1:imrecv)
             do iloc=1,imrecv
                atcfcomp%num_at_glob(nag(iloc))=nag(iloc)
             end do

             call MPI_RECV(buffer(1:3,1:imrecv),imrecv3, NDM_MPI_REAL_DOUBLE, proc_source, 10002, icomm, status, ierr)
             do iloc=1,imrecv
                atcfcomp%xp(1:3,nag(iloc))=buffer(1:3,iloc)
             end do
                
             call MPI_RECV(buffer(1:3,1:imrecv),imrecv3, NDM_MPI_REAL_DOUBLE, proc_source, 10003, icomm, status, ierr)
             do iloc=1,imrecv
                atcfcomp%fp(1:3,nag(iloc))=buffer(1:3,iloc)
             end do
             call MPI_RECV(ibuffer(1:imrecv),imrecv, MPI_INTEGER, proc_source, 10005, icomm, status, ierr)
             do iloc=1,imrecv
                atcfcomp%ityp(nag(iloc))=ibuffer(iloc)
             end do
             if (allocated(atcfloc%lgul))then
                call MPI_RECV(lbuffer(1:imrecv),imrecv, MPI_LOGICAL, proc_source, 10006, icomm, status, ierr)
                do iloc=1,imrecv
                   atcfcomp%lgul(nag(iloc))=lbuffer(iloc)
                end do
             end if
             do iloc=1,imrecv
                atcfcomp%proc_at(nag(iloc))=proc_source !  ibuffer(1:imrecv)
             end do
             select type(atcfloc)
             class is (atom_config_d)
                select type(atcfcomp)
                class is (atom_config_d)
                   call MPI_RECV(buffer(1:3,1:imrecv),imrecv3, NDM_MPI_REAL_DOUBLE, proc_source, 10007, icomm, status, ierr)
                   do iloc=1,imrecv
                      atcfcomp%vp(1:3,nag(iloc))=buffer(1:3,iloc)
                   end do
                   call MPI_RECV(buffer(1:3,1:imrecv),imrecv3, NDM_MPI_REAL_DOUBLE, proc_source, 10008, icomm, status, ierr)
                   do iloc=1,imrecv
                      atcfcomp%xpp(1:3,nag(iloc))=buffer(1:3,iloc)
                   end do
                end select
             type is (atom_config_e)
                select type(atcfcomp)
                class is (atom_config_e)
                   if((atcfcomp%lsigat).and.(atcfloc%lsigat))then
                      call MPI_RECV(buffer9(1:3,1:3,1:imrecv),imrecv9, NDM_MPI_REAL_DOUBLE, proc_source, 10009, icomm, status, ierr)
                      do iloc=1,imrecv
                         atcfcomp%sigat(1:3,1:3,nag(iloc))=buffer9(1:3,1:3,iloc)
                      end do
                   endif
                   if((atcfcomp%lprteat).and.(atcfloc%lprteat))then
                      call MPI_RECV(buffer1(1:imrecv),imrecv, NDM_MPI_REAL_DOUBLE, proc_source, 10010, icomm, status, ierr)
                      do iloc=1,imrecv
                         atcfcomp%eat(nag(iloc))=buffer1(iloc)
                      end do
                   endif
                   if((atcfcomp%llangevin).and.(atcfloc%llangevin))then
                      call MPI_RECV(buffer(1:3,1:imrecv),imrecv3, NDM_MPI_REAL_DOUBLE, proc_source, 10011, icomm, status, ierr)
                      do iloc=1,imrecv
                         atcfcomp%glangv(1:3,nag(iloc))=buffer(1:3,iloc)
                      end do
                   endif
                end select
             end select
             deallocate(nag)
          end if
       end do


       if (imtot.ne.atcfcomp%im) then
          write(6,*)'atomes perdus ?',idloc, imtot,atcfcomp%im,div%mpi_orig%rank
          call MPI_finalize(ierr)
          stop
       end if
       do icomp=1,imtot
          if (atcfcomp%num_at_glob(icomp).ne.icomp) then
             write(6,*)'atomes mal rangés L2M ?',idloc, icomp,atcfcomp%num_at_glob(icomp)
          end if
       end do
    else
       imloc=atcfloc%im;imloc3=3*imloc; imloc9=3*imloc3
       call MPI_SEND(imloc, 1,   MPI_INTEGER,idmaster ,10001,icomm,ierr)
       call MPI_SEND(atcfloc%num_at_glob(1:imloc), imloc, MPI_INTEGER, idmaster, 10004, icomm, ierr)
       call MPI_SEND(atcfloc%xp(1:3,1:imloc), imloc3, NDM_MPI_REAL_DOUBLE,idmaster ,10002,icomm,ierr)
       call MPI_SEND(atcfloc%fp(1:3,1:imloc), imloc3, NDM_MPI_REAL_DOUBLE,idmaster ,10003,icomm,ierr)
       call MPI_SEND(atcfloc%ityp(1:imloc), imloc, MPI_INTEGER, idmaster, 10005, icomm, ierr)
       if (allocated(atcfloc%lgul))     &
            &call MPI_SEND(atcfloc%lgul(1:imloc), imloc, MPI_LOGICAL, idmaster, 10006, icomm, ierr)
       select type(atcfloc)
       class is (atom_config_d)
          select type(atcfcomp)
          class is (atom_config_d)
             call MPI_SEND(atcfloc%vp(1:3,1:imloc), imloc3, NDM_MPI_REAL_DOUBLE,idmaster ,10007,icomm,ierr)
             call MPI_SEND(atcfloc%xpp(1:3,1:imloc), imloc3, NDM_MPI_REAL_DOUBLE,idmaster ,10008,icomm,ierr)
          end select
       type is (atom_config_e)
          select type(atcfcomp)
          class is (atom_config_e)
             if((atcfcomp%lsigat).and.(atcfloc%lsigat))then
                call MPI_SEND(atcfloc%sigat(1:3,1:3,1:imloc), imloc9, NDM_MPI_REAL_DOUBLE,idmaster ,10009,icomm,ierr)
             endif
             if((atcfcomp%lprteat).and.(atcfloc%lprteat))then
                call MPI_SEND(atcfloc%eat(1:imloc), imloc, NDM_MPI_REAL_DOUBLE,idmaster ,10010,icomm,ierr)
             endif
             if((atcfcomp%llangevin).and.(atcfloc%llangevin))then
                call MPI_SEND(atcfloc%glangv(1:3,1:imloc), imloc3, NDM_MPI_REAL_DOUBLE,idmaster ,10011,icomm,ierr)
             endif
          end select
       end select
       
    end if
    call MPI_barrier(icomm,ierr)
#endif       

    return
  end subroutine vers_master_atom



  
  subroutine master2loc_atom(atcfcomp,atcfloc,div)

    class(atom_config)::atcfloc
    class(atom_config)::atcfcomp
    type(para_config)::div
#ifdef PARA
    integer::idmaster,idloc,icomm,npim ! proc master,proclocal ,communicateur
    real(double),allocatable::buffer(:,:),buffer9(:,:,:),buffer1(:)
    integer,allocatable:: ibuffer(:)
    logical,allocatable::lbuffer(:)
    integer::imloc,imloc3,iproc,imrecv,ideb,ifin,imrecv3,imrecv9,icomp
    integer::imcomp,imtot,proc_source,ic,ns,iloc,i,iu
    logical,allocatable::mask(:)
    idmaster=0
    idloc=div%mpi_image%rank
    npim=div%mpi_image%nproc
    icomm=div%mpi_image%comm
    imcomp=atcfcomp%im

    
    if (idloc==idmaster) then
    allocate(mask(atcfcomp%im))
    do icomp=1,atcfcomp%im
          if (atcfcomp%num_at_glob(icomp).ne.icomp) then
             write(6,*)'atomes mal rangés M2L?',icomp,atcfcomp%num_at_glob(icomp)
          end if
       end do
       
       imtot=0
       do iproc=0,npim-1
          mask(:)=(atcfcomp%proc_at(1:imcomp)==iproc)
          ns=count(atcfcomp%proc_at(1:imcomp)==iproc)
          if (iproc==idmaster) then

             iloc=0
             do i=1,imcomp
                if (mask(i))then

                   iloc=iloc+1
                   atcfloc%xp(:,iloc)=atcfcomp%xp(:,i)
                   atcfloc%fp(:,iloc)=atcfcomp%fp(:,i)
                   atcfloc%num_at_glob(iloc)=atcfcomp%num_at_glob(i)
                   atcfloc%ityp(iloc)=atcfcomp%ityp(i)
                   atcfloc%lgul(iloc)=atcfcomp%lgul(i)
                   atcfloc%proc_at(iloc)=idmaster
                   select type(atcfloc)
                   class is (atom_config_d)
                      select type (atcfcomp)
                      class is (atom_config_d)
                         atcfloc%vp(:,iloc)=atcfcomp%vp(:,i)
                         atcfloc%xpp(:,iloc)=atcfcomp%xpp(:,i)
                      end select
                   type is (atom_config_e)
                      select type (atcfcomp)
                      class is (atom_config_e)
                         if((atcfcomp%lsigat).and.(atcfloc%lsigat))then
                            atcfloc%sigat(:,:,iloc)=atcfcomp%sigat(:,:,i)
                         endif
                         if((atcfcomp%lprteat).and.(atcfloc%lprteat))then
                            atcfloc%eat(iloc)=atcfcomp%eat(i)
                         endif
                         if((atcfcomp%llangevin).and.(atcfloc%llangevin))then
                            atcfloc%glangv(:,iloc)=atcfcomp%glangv(:,i)
                         endif
                      end select
                   end select
                end if
             end do
             atcfloc%im=ns
          else
             if (allocated(buffer)) then
                deallocate(buffer);deallocate(ibuffer);deallocate(lbuffer)
             end if
             allocate(buffer(3,ns));allocate(ibuffer(ns));allocate(lbuffer(ns))

             call MPI_SEND(ns, 1,   MPI_INTEGER,iproc,20001,icomm,ierr)
             call fillbuffer3D(buffer,atcfcomp%xp,mask)
             call MPI_SEND(buffer, 3*ns, NDM_MPI_REAL_DOUBLE,iproc ,20002,icomm,ierr)
             call fillbuffer3D(buffer,atcfcomp%fp,mask)
             call MPI_SEND(buffer, 3*ns, NDM_MPI_REAL_DOUBLE,iproc ,20003,icomm,ierr)
             call fillbuffer1D(ibuffer,atcfcomp%num_at_glob,mask)
             call MPI_SEND(ibuffer, ns, MPI_INTEGER,iproc ,20004,icomm,ierr)
             call fillbuffer1D(ibuffer,atcfcomp%ityp,mask)
             call MPI_SEND(ibuffer, ns, MPI_INTEGER,iproc ,20005,icomm,ierr)                          
             if (allocated(atcfloc%lgul)) then
                call fillbuffer1D(lbuffer,atcfcomp%lgul,mask)
                call MPI_SEND(lbuffer, ns, MPI_LOGICAL,iproc ,20006,icomm,ierr)                          
             end if
             select type(atcfloc)
             class is (atom_config_d)
                select type (atcfcomp)
                class is (atom_config_d)
                   call fillbuffer3D(buffer,atcfcomp%vp,mask)
                   call MPI_SEND(buffer, 3*ns, NDM_MPI_REAL_DOUBLE,iproc ,20007,icomm,ierr)
                   call fillbuffer3D(buffer,atcfcomp%xpp,mask)
                   call MPI_SEND(buffer, 3*ns, NDM_MPI_REAL_DOUBLE,iproc ,20008,icomm,ierr)
                end select
             type is (atom_config_e)
                select type (atcfcomp)
                class is (atom_config_e)
                   if((atcfcomp%lsigat).and.(atcfloc%lsigat))then
                      allocate(buffer9(3,3,ns))
                      call fillbuffer9D(buffer9,atcfcomp%sigat,mask)
                      call MPI_SEND(buffer9, 9*ns, NDM_MPI_REAL_DOUBLE,iproc ,20009,icomm,ierr)
                   endif
                   if((atcfcomp%lprteat).and.(atcfloc%lprteat))then
                      allocate(buffer1(ns))
                      call fillbuffer1D(buffer1,atcfcomp%eat,mask)
                      call MPI_SEND(buffer1, ns, NDM_MPI_REAL_DOUBLE,iproc ,20010,icomm,ierr)
                   endif
                   if((atcfcomp%llangevin).and.(atcfloc%llangevin))then
                      call fillbuffer3D(buffer,atcfcomp%glangv,mask)
                      call MPI_SEND(buffer, 3*ns, NDM_MPI_REAL_DOUBLE,iproc ,20011,icomm,ierr)
                   endif
                end select
             end select
          end if
       end do
       
    else
       call MPI_RECV(imrecv,1, MPI_INTEGER, idmaster, 20001, icomm, status, ierr)
       atcfloc%im=imrecv;imrecv3=3*imrecv;imrecv9=9*imrecv
       allocate(buffer(3,imrecv));allocate(ibuffer(imrecv));allocate(lbuffer(imrecv))
       
       call MPI_RECV(buffer,imrecv3, NDM_MPI_REAL_DOUBLE, idmaster, 20002, icomm, status, ierr)
       atcfloc%xp(1:3,1:imrecv)=buffer(1:3,1:imrecv)
       call MPI_RECV(buffer,imrecv3, NDM_MPI_REAL_DOUBLE, idmaster, 20003, icomm, status, ierr)
       atcfloc%fp(1:3,1:imrecv)=buffer(1:3,1:imrecv)
       
       call MPI_RECV(ibuffer,imrecv, MPI_INTEGER, idmaster, 20004, icomm, status, ierr)
       atcfloc%num_at_glob(1:imrecv)=ibuffer(1:imrecv)
       call MPI_RECV(ibuffer,imrecv, MPI_INTEGER, idmaster, 20005, icomm, status, ierr)
       atcfloc%ityp(1:imrecv)=ibuffer(1:imrecv)
       if (allocated(atcfloc%lgul) )then
          call MPI_RECV(lbuffer,imrecv, MPI_LOGICAL, idmaster, 20006, icomm, status, ierr)
          atcfloc%lgul(1:imrecv)=lbuffer(1:imrecv)
       end if
       atcfloc%proc_at(1:imrecv)=idloc
             select type(atcfloc)
             class is (atom_config_d)
                select type (atcfcomp)
                class is (atom_config_d)
                   call MPI_RECV(buffer,imrecv3, NDM_MPI_REAL_DOUBLE, idmaster, 20007, icomm, status, ierr)
                   atcfloc%vp(1:3,1:imrecv)=buffer(1:3,1:imrecv)
                   call MPI_RECV(buffer,imrecv3, NDM_MPI_REAL_DOUBLE, idmaster, 20008, icomm, status, ierr)
                   atcfloc%xpp(1:3,1:imrecv)=buffer(1:3,1:imrecv)
                end select
             type is (atom_config_e)
                select type (atcfcomp)
                class is (atom_config_e)
                   if((atcfcomp%lsigat).and.(atcfloc%lsigat))then
                   call MPI_RECV(buffer9,imrecv9, NDM_MPI_REAL_DOUBLE, idmaster, 20009, icomm, status, ierr)
                   atcfloc%sigat(1:3,1:3,1:imrecv)=buffer9(1:3,1:3,1:imrecv)
                   endif
                   if((atcfcomp%lprteat).and.(atcfloc%lprteat))then
                   call MPI_RECV(buffer1,imrecv, NDM_MPI_REAL_DOUBLE, idmaster, 20010, icomm, status, ierr)
                   atcfloc%eat(1:imrecv)=buffer1(1:imrecv)
                   endif
                   if((atcfcomp%llangevin).and.(atcfloc%llangevin))then
                      call MPI_RECV(buffer,imrecv3, NDM_MPI_REAL_DOUBLE, idmaster, 20011, icomm, status, ierr)
                      atcfloc%glangv(1:3,1:imrecv)=buffer(1:3,1:imrecv)

                   endif
                end select
             end select
    end if
    
    return
#endif
 
  end subroutine master2loc_atom
 

    
end module atomconfig


