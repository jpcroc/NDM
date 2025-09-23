module atomconfig
   USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m,only:double,long
  USE Mat_utils_mod,only: fillbuffer3D,fillbuffer1D,fillbuffer9D
  use gen_com_m,only:rang,lspacendm
#ifdef PARA
  use mpi
  USE Tpara,only:NDM_MPI_REAL_DOUBLE,ierr,mpi_communicator,endmpi,nprocspace



#endif
  use paraconfig,only:para_config
  use Tpara,only:mpi_communicator,endmpi

  implicit none
#ifdef PARA
  integer,dimension(MPI_STATUS_SIZE):: status  ! statut de la communication
#endif

  type atom_config ! type minimal des configurations atomiques. Tous les composants seront toujours alloué (im_glob seulement si PARA)
     integer::im=0,imm=0,im_glob,imm_glob ! im nb d'atomes sur le proc, taille des tableaux sur le proc, im_glob nombre d'atomes en tout,imm_glob, taille des tableaux complets
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
     real(double),allocatable:: distance(:,:) ! indice de tous les distances
     logical, allocatable:: lgul(:)
     integer,allocatable::num_at_glob(:)
     logical :: lglock
#ifdef PARA
     integer,allocatable::proc_at(:) ! tableau de taille im_glob total indiquant le numéro du proc qui gère l'atome
     integer::imf=0 ! indice du dernier atome fantome (atomes fantomes entre im+1 et imf
#endif     
   contains
     procedure, pass::init=>init_atom_config
     procedure, pass::lgcheck
     procedure, pass::Eegal
     procedure, pass::copy_atom=>copy_atom_b
     procedure, pass::dealloc=>dealloc_atom_config
     procedure, pass::vers_master=>vers_master_atom !(atcfloc,atcfcomp,div)
     procedure, pass::master2loc=>master2loc_atom !(atcfcomp,atcfloc,div)
     procedure, pass::copy_config
     procedure, pass::print
     procedure, pass::pack
     procedure, pass::fab
     procedure, pass::sort
     procedure, pass::backto
     procedure, pass::add2conf
     procedure, pass::extend
     procedure, pass::deftype
     procedure, pass::send2all=>s2a_atom
     procedure, pass::send2proc=>s2p_atom
     procedure, pass::recv=>rcv_atom
     procedure, pass::s1at2p_rm
     procedure, pass::recv1at
     procedure, pass::zero=>zero_atom
     procedure, pass::switch_atom
     procedure, pass::print_type

     !
  end type atom_config

  type, extends (atom_config):: atom_config_d ! type dynamique des configurations atomiques(+vp). vp  seront toujours allouées
     real(double),allocatable::vp(:,:)
   contains
     procedure, pass::copy_atom=>copy_atom_d
     procedure, pass::dealloc=>dealloc_atom_config_d
     !#ifdef PARA
     procedure, pass::zero=>zero_atom_d
     !#endif     
  end type atom_config_d
  
  type, extends (atom_config_d):: atom_config_e ! type étendu des configurations atomiques avec quantités optionelles Ces quantités seront allouées en fonction dss logical
     real(double),allocatable ::xpp(:,:)
     logical ::lxpp
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
     procedure, pass::zero=>zero_atom_e
   
     !#endif     
  end type atom_config_e
  type, extends (atom_config_e)::atom_config_arps
     real(double),allocatable::rho(:),fpr(:,:),fpg(:,:)
     integer, allocatable::mov(:)
     contains
       procedure, pass::copy_atom=>copy_atom_arps
       procedure, pass::dealloc=>dealloc_atom_config_arps
       !#ifdef PARA
!       procedure, pass::zero=>zero_atom_arps
  end type atom_config_arps

  private ::buffersizes
contains
  !initialisations
  subroutine print_type(atcf,mess)
    class(atom_config),intent(in)::atcf
    character(len=*),optional::mess
    if (present(mess)) then
       write(6,*)'atomfigPRINT ',mess
    else
       write(6,*)'atomfigPRINT'
    end if
    write(6,*)'im imm im_glob imm_glob',atcf%im,atcf%imm,atcf%im_glob,atcf%imm_glob
    select type (atcf)
    type is (atom_config)
       write(6,*)'atomfig'
    type is (atom_config_d)
       write(6,*)'atomfigD'
    class is (atom_config_e)
       write(6,*)'atomfigE'
    type is (atom_config_arps)
       write(6,*)'atomfigARPS'
    end select
    select type (atcf)
    class is (atom_config_e)
       write(6,*)'FLAGSFF xpp eat sigat langevin ax'
       write(6,*)'FLAGSFF', atcf%lxpp,atcf%lprteat,atcf%lsigat,atcf%llangevin,atcf%lax
    end select
    write(6,*)'TYPE PRECISE ? SI NON extension'

  end subroutine print_type

!**************  INIT   *********
  
  subroutine init_atom_config(atconf,imin,immin,ltabvois,nvois,rvois,lreallocate,im_glob,imm_glob)
    class(atom_config),intent(inout)::atconf
    integer,intent(in):: imin
    logical,optional, intent(in)::ltabvois,lreallocate
    integer, optional::nvois,immin,im_glob,imm_glob
    real(double),optional::rvois
    integer::nv
    logical ::ltbv,lrealloc
    real(double)::rv
    
    lrealloc=.false.
    if (present(lreallocate))lrealloc=lreallocate
    ltbv=.false.
    if(present(ltabvois)) then
       ltbv=ltabvois
    end if
    nv=0 ;  if(present(nvois))nv=nvois
    rv=0;  if(present(rvois))rv=rvois

    if (ltbv) then
       if(.not.(present(rvois)))then
          write(6,*)'rvois must be set in initialization of atcf when ltabvois =True'
          call arret_ndm
       end if
       if (rv==0) then 
          write(6,*)'rvois must be set to non zero in initialization of atcf when ltabvois =True'
          call arret_ndm
       end if
    else
       if(present(rvois)) then
          if (rvois.ne.0)then
             write(6,*)'rvois must NOT be set in initialization of atcf when ltabvois =False'
             call arret_ndm
          end if
       end if
    end if
 


    atconf%im=imin
    if (present(im_glob))then
 !      write(6,*)'PRESENT imglob',im_glob
       atconf%im_glob=im_glob
!    else
!       atconf%im_glob=0
    end if
    if (present(imm_glob))then
!       write(6,*)'PRESENT imMglob',imm_glob
       atconf%imm_glob=imm_glob
!    else
!       atconf%imm_glob=0
    end if
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
    atconf%ityp=0;atconf%xp=0;atconf%fp=0;atconf%ielat=0; atconf%lgul=.false.;
    atconf%num_at_glob=0
#ifdef PARA
    atconf%proc_at=-1
#endif
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
       atconf%rvois=rv ! Ajouté pour éviter des erreurs d'initilaisations
    end if

    select type (atconf)
    class is (atom_config_d)
       if ((lrealloc).and.(allocated(atconf%vp)))then
          deallocate(atconf%vp)
       end if
       if (.not.allocated(atconf%vp))then
          allocate(atconf%vp(3,atconf%imm))
       end if
       atconf%vp=0
    end select
    select type (atconf)
    class is (atom_config_e)
       if (atconf%lxpp) then 
          if ((lrealloc).and.(allocated(atconf%xpp)))then
             deallocate(atconf%xpp)
          end if
          if (.not.allocated(atconf%xpp))then
             allocate(atconf%xpp(3,atconf%imm))
          end if
          atconf%xpp=0
       end if
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
    select type (atconf)
    class is (atom_config_arps)
       if ((lrealloc).and.(allocated(atconf%mov)))then
          deallocate(atconf%mov)
       end if
       if (.not.allocated(atconf%mov))then
          allocate(atconf%mov(atconf%imm))
       end if
       atconf%mov=0
       if ((lrealloc).and.(allocated(atconf%fpr)))then
          deallocate(atconf%fpr)
       end if
       if (.not.allocated(atconf%fpr))then
          allocate(atconf%fpr(3,atconf%imm))
       end if
       atconf%fpr=0
       if ((lrealloc).and.(allocated(atconf%fpg)))then
          deallocate(atconf%fpg)
       end if
!       if (.not.allocated(atconf%rho))then
!          allocate(atconf%rho(atconf%imm))
!       end if
    end select

    atconf%icaltabt=0
    atconf%lglock=.false.
  end subroutine init_atom_config

  !copie d'un élément
  
  subroutine copy_atom_b(atsource,i,atcible,j,lextend,caracT)
    class(atom_config), intent(in)::atsource
    integer,intent(in):: i
    class(atom_config), intent(inout)::atcible
    integer,intent(in):: j
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    logical , optional, intent(in) :: lextend
    logical ::let
    integer::iw,nvj,idwi,idwj,imm_min

    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   
    ! v=vp,r=xpp
    ! u=eat,g=glangv;a=ax;s=sigat; m=mov(arps)   
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugasm'
    else
       carac=caracT
    end if
    
    let=.false.
    if (present(lextend))let= lextend
    if (j.gt.atcible%imm) then
       if (let) then
          imm_min=j-atcible%imm
          call atcible%extend(imm_min)
       else
          write(6,*)'PB COPY', rang,j,atcible%imm,i,atsource%imm
          write(6,*)'copy of an atom element is not possible , target size too small' ,rang
          call arret_ndm
       end if
    end if

    if(scan('x',carac).ne.0)    atcible%xp(:,j)=atsource%xp(:,i)
    if(scan('f',carac).ne.0)    atcible%fp(:,j)=atsource%fp(:,i)
    if(scan('i',carac).ne.0)    atcible%ityp(j)=atsource%ityp(i)
    if(scan('n',carac).ne.0)    atcible%num_at_glob(j)=atsource%num_at_glob(i)
#ifdef PARA
    if(scan('p',carac).ne.0)    atcible%proc_at(j)=atsource%proc_at(i)
#endif    
    if(scan('e',carac).ne.0)    atcible%ielat(j)=atsource%ielat(i)
    if(scan('l',carac).ne.0)    atcible%lgul(j)=atsource%lgul(i)
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
    class (atom_config),allocatable::attemp
!!$    type(atom_config):: attemp_b
!!$    type(atom_config_d):: attemp_d
!!$    type(atom_config_e):: attemp_e
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
    call atcf%deftype(attemp)
    call attemp%init(imcn,immcn,atcf%ltabvois,nvois,rvois)
    call atcf%copy_config(attemp,lrescl=.false.)
    call attemp%copy_config(atcf,lrescl=.true.)


  end subroutine extend



  subroutine copy_atom_d(atsource,i,atcible,j,lextend,caracT)
    class(atom_config_d), intent(in)::atsource
    integer,intent(in):: i
    class(atom_config), intent(inout)::atcible
    integer,intent(in):: j
    logical, optional,intent(in):: lextend
    logical::let
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   
    ! v=vp,r=xpp
    ! u=eat,g=glangv;a=ax;s=sigat; m=mov(arps)   
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugasm'
    else
       carac=caracT
    end if
    
    let=.false.
    if (present(lextend))let=lextend
    call copy_atom_b(atsource,i,atcible,j,let,carac)
    select type(atcible)
       class is (atom_config_d)
       select type (atsource)
          class is (atom_config_d)
             if(scan('v',carac).ne.0)             atcible%vp(:,j)=atsource%vp(:,i)
       end select
    end select
  end subroutine copy_atom_d

  subroutine copy_atom_e(atsource,i,atcible,j,lextend,caracT)
    class(atom_config_e), intent(in)::atsource
    integer,intent(in):: i
    class(atom_config), intent(inout)::atcible
    integer,intent(in):: j
    logical, optional,intent(in):: lextend
    logical:: let
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   
    ! v=vp,r=xpp
    ! u=eat,g=glangv;a=ax;s=sigat; m=mov(arps)   
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugasm'
    else
       carac=caracT
    end if

    let=.false.
    if (present(lextend))let=lextend
    call copy_atom_d(atsource,i,atcible,j,let)
    select type(atcible)
       class is (atom_config_e)
       select type (atsource)
       class is (atom_config_e)
          if ((atcible%lxpp).and.(atsource%lxpp).and.(scan('r',carac).ne.0))        atcible%xpp(:,j)=atsource%xpp(:,i)
          if ((atcible%lprteat).and.(atsource%lprteat).and.(scan('u',carac).ne.0)) atcible%eat(j)=atsource%eat(i)
          if ((atcible%lsigat).and.(atsource%lsigat).and.(scan('s',carac).ne.0)) atcible%sigat(:,:,j)=atsource%sigat(:,:,i)
          if ((atcible%llangevin).and.(atsource%llangevin).and.(scan('g',carac).ne.0)) atcible%glangv(:,j)=atsource%glangv(:,i)
          if ((atcible%lax).and.(atsource%lax).and.(scan('a',carac).ne.0)) atcible%ax(:,j)=atsource%ax(:,i)
       end select
    end select
  end subroutine copy_atom_e
  subroutine copy_atom_arps(atsource,i,atcible,j,lextend,caracT)
    class(atom_config_arps), intent(in)::atsource
    integer,intent(in):: i
    class(atom_config), intent(inout)::atcible
    integer,intent(in):: j
    logical, optional,intent(in):: lextend
    logical:: let
        character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   
    ! v=vp,r=xpp
    ! u=eat,g=glangv;a=ax;s=sigat; m=mov(arps)   
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugasm'
    else
       carac=caracT
    end if

    let=.false.
    if (present(lextend))let=lextend
    call copy_atom_e(atsource,i,atcible,j,let,carac)
    if(scan('m',carac).ne.0) then 
       select type(atcible)
       class is (atom_config_arps)
          select type (atsource)
          class is (atom_config_arps)
             if(allocated(atcible%rho))atcible%rho(j)=atsource%rho(i)
             atcible%fpr(:,j)=atsource%fpr(:,i)
             if (allocated(atcible%fpg))atcible%fpg(:,j)=atsource%fpg(:,i)
             atcible%mov(j)=atsource%mov(i)
          end select
       end select
    end if
  end subroutine copy_atom_arps
  
  subroutine recv1at(atcf, rgem,mpic,mypp,caracT)
    class(atom_config):: atcf
    type(mpi_communicator),intent(in)::mpic
    integer,intent(in)::rgem,mypp
    integer:: size1,size3,size9
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    integer::nvi,nvr,sizeI,sizeR,nvl,sizel,ivi,ivr,ivl,csi,csr,csl
    integer, dimension (0:26):: Ipos1,Iposf,Rpos1,Rposf,Lpos1,Lposf
    integer,allocatable:: ibuffer(:)
    logical,allocatable::lbuffer(:)!,mask(:)
    real(double),allocatable::rbuffer(:)
    integer:: cst(3),nmask
    nmask=1
!    allocate(mask(atcf%imm))
!    mask(:)=.false.; mask(iat)=.true.

    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   
    ! v=vp,r=xpp
    ! u=eat,g=glangv;a=ax;s=sigat    
    size1=1;size3=3*size1;size9=3*size3
#ifdef PARA
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugasm'
    else
       carac=caracT
    end if

    nvi=0; sizeI=0
    Ipos1(:)=0;Iposf(:)=0
    nvR=0; sizeR=0
    Rpos1(:)=0;Rposf(:)=0
    nvl=0; sizel=0
    lpos1(:)=0;lposf(:)=0
    ivi=0;ivl=0;ivr=0

    call buffersizes(atcf,sizeI,sizeR,sizel,IposF,Rposf,Lposf,carac,nmask)
    allocate(ibuffer(sizeI));     allocate(Lbuffer(sizeL));     allocate(Rbuffer(sizeR)); 
    call mpic%recv (cst,rgem,112)
    if (cst(1).ne.sizeI) then
       write(6,*)'erreur CST1A ',sizeI,cst(1)
       call arret_ndm
    end if
    if (cst(2).ne.sizel) then
       write(6,*)'erreur CST2A ',sizel,cst(2)
       call arret_ndm
    end if
    if (cst(3).ne.sizeR) then
       write(6,*)'erreur CST3A ',sizer,cst(3)
       call arret_ndm
   end if
   
    if (cst(1).ne.0)call mpic%recv(ibuffer,rgem,314)
    if (cst(2).ne.0)call mpic%recv(lbuffer,rgem,315)
    if (cst(3).ne.0)call mpic%recv(Rbuffer,rgem,316)

    call addatim (atcf,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf,atcf%imm)
    write(6,*)'PATa',    atcf%proc_at(atcf%im)
    atcf%proc_at(atcf%im)=mypp
    write(6,*)'PATb',    atcf%proc_at(atcf%im)
    return

#endif

  end subroutine recv1at

  subroutine s1at2p_rm (atcf, rgcib,iat,mpic,caracT)
    class(atom_config):: atcf
    type(mpi_communicator),intent(in)::mpic
    integer,intent(in)::rgcib,iat! target proc and atom number
    integer:: size1,size3,size9
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    integer::nvi,nvr,sizeI,sizeR,nvl,sizel,csi,csr,csl
    integer, dimension (0:26):: Iposf,Rposf,Lposf
    integer,allocatable:: ibuffer(:)
    logical,allocatable::lbuffer(:)
    real(double),allocatable::rbuffer(:)
    integer:: cst(3),nat
    logical, allocatable::mask(:)
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   
    ! v=vp,r=xpp
    ! u=eat,g=glangv;a=ax;s=sigat; m=mov(arps)   
    size1=1;size3=3*size1; size9=3*size3
    allocate(mask(atcf%imm))
    mask(:)=.false.; mask(iat)=.true.
#ifdef PARA
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugasm'
    else
       carac=caracT
    end if
   write(6,*)'IAT',rang,iat
    nvi=0
    sizeI=0
    Iposf(:)=0
    nvR=0; sizeR=0
    Rposf(:)=0
    nvl=0
    sizel=0
    lposf(:)=0
    nat=1
    call buffersizes(atcf,sizeI,sizeR,sizel,IposF,Rposf,Lposf,carac,nat)
    allocate(ibuffer(sizeI));     allocate(Lbuffer(sizeL));     allocate(Rbuffer(sizeR));
    call buildbuff(atcf,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf,maskR=mask,nmask=nat)
    
    csT(1)=csi
    csT(2)=csl
    csT(3)=csR
    call mpic%send (cst,rgcib,112)
    if (csi.ne.0)call mpic%send(ibuffer,rgcib,314)
    if (csl.ne.0)call mpic%send(lbuffer,rgcib,315)
    if (csR.ne.0)call mpic%send(Rbuffer,rgcib,316)
    call atcf%switch_atom(iat,atcf%im)
    
    atcf%im=atcf%im-1


#endif
  end subroutine s1at2p_rm
  
  subroutine s2p_atom (atcf, rgcib,mpic,caracT)
    class(atom_config):: atcf
    type(mpi_communicator),intent(in)::mpic
    integer,intent(in)::rgcib
    integer:: size1,size3,size9
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    integer::nvi,nvr,sizeI,sizeR,nvl,sizel,csi,csr,csl
    integer, dimension (0:26):: Iposf,Rposf,Lposf
    integer,allocatable:: ibuffer(:)
    logical,allocatable::lbuffer(:)
    real(double),allocatable::rbuffer(:)
    integer:: cst(3)
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   
    ! v=vp,r=xpp
    ! u=eat,g=glangv;a=ax;s=sigat; m=mov(arps)   
    size1=atcf%imm;size3=3*size1; size9=3*size3
#ifdef PARA
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugasm'
    else
       carac=caracT
    end if
   
    nvi=0
    sizeI=0
    Iposf(:)=0
    nvR=0; sizeR=0
    Rposf(:)=0
    nvl=0
    sizel=0
    lposf(:)=0

    call buffersizes(atcf,sizeI,sizeR,sizel,IposF,Rposf,Lposf,carac)
    allocate(ibuffer(sizeI));     allocate(Lbuffer(sizeL));     allocate(Rbuffer(sizeR)); 
    call buildbuff(atcf,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf)
    
    csT(1)=csi
    csT(2)=csl
    csT(3)=csR
    call mpic%send (cst,rgcib,112)
    if (csi.ne.0)call mpic%send(ibuffer,rgcib,114)
    if (csl.ne.0)call mpic%send(lbuffer,rgcib,115)
    if (csR.ne.0)call mpic%send(Rbuffer,rgcib,116)

#endif
  end subroutine s2p_atom



  subroutine rcv_atom (atcf, rgem,mpic,caracT)
    class(atom_config):: atcf
    type(mpi_communicator),intent(in)::mpic
    integer,intent(in)::rgem
    integer:: size1,size3,size9
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    integer::nvi,nvr,sizeI,sizeR,nvl,sizel,ivi,ivr,ivl,csi,csr,csl
    integer, dimension (0:26):: Ipos1,Iposf,Rpos1,Rposf,Lpos1,Lposf
    integer,allocatable:: ibuffer(:)
    logical,allocatable::lbuffer(:)
    real(double),allocatable::rbuffer(:)
    integer:: cst(3)
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   
    ! v=vp,r=xpp
    ! u=eat,g=glangv;a=ax;s=sigat    
    size1=atcf%imm;size3=3*size1;size9=3*size3
#ifdef PARA
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugasm'
    else
       carac=caracT
    end if

    nvi=0; sizeI=0
    Ipos1(:)=0;Iposf(:)=0
    nvR=0; sizeR=0
    Rpos1(:)=0;Rposf(:)=0
    nvl=0; sizel=0
    lpos1(:)=0;lposf(:)=0
    ivi=0;ivl=0;ivr=0

    call buffersizes(atcf,sizeI,sizeR,sizel,IposF,Rposf,Lposf,carac)
    allocate(ibuffer(sizeI));     allocate(Lbuffer(sizeL));     allocate(Rbuffer(sizeR)); 
    call mpic%recv (cst,rgem,112)
    if (cst(1).ne.sizeI) then
       write(6,*)'erreur CST1 ',sizeI,cst(1)
       call arret_ndm
    end if
    if (cst(2).ne.sizel) then
       write(6,*)'erreur CST2 ',sizel,cst(2)
       call arret_ndm
    end if
    if (cst(3).ne.sizeR) then
       write(6,*)'erreur CST3 ',sizer,cst(3)
       call arret_ndm
   end if
   
    if (cst(1).ne.0)call mpic%recv(ibuffer,rgem,114)
    if (cst(2).ne.0)call mpic%recv(lbuffer,rgem,115)
    if (cst(3).ne.0)call mpic%recv(Rbuffer,rgem,116)

    call copybuff (atcf,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf,atcf%imm)

    return

#endif
  end subroutine rcv_atom



  subroutine s2a_atom(atcf,rgemet,mpic,caracT)
    type(mpi_communicator),intent(in)::mpic
    class(atom_config)::atcf
    integer,intent(in)::rgemet
    integer:: size1,size3,size9
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at
    ! v=vp,r=xpp
    ! u=eat,g=glangv;a=ax;s=sigat    
    integer::nvi,nvr,sizeI,sizeR,nvl,sizel,csi,csr,csl
    integer, dimension (0:26):: Iposf,Rposf,Lposf
    integer,allocatable:: ibuffer(:)
    logical,allocatable::lbuffer(:)
    real(double),allocatable::rbuffer(:)
    integer:: cst(3)


    size1=atcf%imm;size3=3*size1; size9=3*size3
#ifdef PARA
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugasm'
    else
       carac=caracT
    end if
   
    nvi=0
    sizeI=0
    Iposf(:)=0
    nvR=0; sizeR=0
    Rposf(:)=0
    nvl=0
    sizel=0
    lposf(:)=0

    call buffersizes(atcf,sizeI,sizeR,sizel,IposF,Rposf,Lposf,carac)
    allocate(ibuffer(sizeI));     allocate(Lbuffer(sizeL));     allocate(Rbuffer(sizeR)); 
    call buildbuff(atcf,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf)
    
    csT(1)=csi
    csT(2)=csl
    csT(3)=csR
    call mpic%bcast (rgemet,cst)
    if (csi.ne.0)call mpic%bcast(rgemet,ibuffer)
    if (csl.ne.0)call mpic%bcast(rgemet,lbuffer)
    if (csR.ne.0)call mpic%bcast(rgemet,Rbuffer)

    call copybuff (atcf,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf,atcf%imm)

#endif
  end subroutine s2a_atom
  

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

  end subroutine zero_atom_d

  subroutine zero_atom_e (atcf)
    class(atom_config_e):: atcf

    call zero_atom_d(atcf)

    if (atcf%lxpp)    atcf%xpp=0.
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
  
  subroutine Eegal(atsource,atcible)
    class (atom_config),intent(in)::atsource
    class(atom_config)::atcible
    select type (atcible)
    class is (atom_config_e)
       select type (atsource)
       class is (atom_config_e)
          atcible%lax=atsource%lax
          atcible%lprteat=atsource%lprteat
          atcible%llangevin=atsource%llangevin
          atcible%lsigat=atsource%lsigat
          atcible%lxpp=atsource%lxpp
       end select
    end select
    atcible%ltabvois=atsource%ltabvois
  end subroutine Eegal
  
  ! copie d'une config entière vers config de base
  subroutine copy_config (atsource,atcible,lrescl)
    class(atom_config),intent(in)::atsource
    class(atom_config)::atcible
    logical,intent(in)::lrescl
    integer::nvois
    real(double)::rvois
    logical :: lstop
    call atsource%Eegal(atcible)
    if (lrescl) then
!       if (atcible%imm.ne.atsource%imm) then
       call atcible%dealloc
       if (atsource%ltabvois)then
          nvois=atsource%nvois; rvois=atsource%rvois
       else
          nvois=0;rvois=0
       end if
       call atcible%init(atsource%im,atsource%imm,atcible%ltabvois,nvois,rvois,&
            &im_glob=atsource%im_glob,imm_glob=atsource%imm_glob)
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
          call arret_ndm
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
    
    if ((atsource%ltabvois).and.(atcible%ltabvois).and.(atsource%nvois.ne.0)) then
       atcible%iwmax(1:atsource%imm)=atsource%iwmax(1:atsource%imm)
       atcible%indi(1:atsource%nvois)=atsource%indi(1:atsource%nvois)
    end if
    ! atsource et atcible sont au moins_d
    select type (atsource)
       class is (atom_config_d)
       select type (atcible)
       class is (atom_config_d)
!          write(6,*)'cible',atcible%imm,atcible%im,size(atcible%vp)
!          write(6,*)'source',atsource%imm,atsource%im,size(atsource%vp)
          atcible%vp(:,1:atsource%imm)=atsource%vp(:,1:atsource%imm)
       end select
    end select
    ! atsource et atcible sont _e    
    select type (atsource)
       class is (atom_config_e)
       select type (atcible)
       class is (atom_config_e)
          if((atcible%lxpp).and.(atsource%lxpp))atcible%xpp(1:3,1:atsource%imm)=atsource%xpp(1:3,1:atsource%imm)
          if((atcible%lprteat).and.(atsource%lprteat))atcible%eat(1:atsource%imm)=atsource%eat(1:atsource%imm)
          if((atcible%lsigat).and.(atsource%lsigat))atcible%sigat(:,:,1:atsource%imm)=atsource%sigat(:,:,1:atsource%imm)
          if((atcible%llangevin).and.(atsource%llangevin))atcible%glangv(:,1:atsource%imm)=atsource%glangv(:,1:atsource%imm)
          if((atcible%lax).and.(atsource%lax))atcible%ax(:,1:atsource%imm)=atsource%ax(:,1:atsource%imm)
       end select
    end select
    select type (atsource)
       class is (atom_config_arps)
       select type (atcible)
          class is (atom_config_arps)
             atcible%mov(1:atsource%imm)=atsource%mov(1:atsource%imm)
            if(allocated(atcible%rho))atcible%rho(1:atsource%imm)=atsource%rho(1:atsource%imm)
             atcible%fpr(:,1:atsource%imm)=atsource%fpr(:,1:atsource%imm)
             if(allocated(atcible%fpg))             atcible%fpg(:,1:atsource%imm)=atsource%fpg(:,1:atsource%imm)
       end select
    end select

    atcible%im_glob=atsource%im_glob
    atcible%imm_glob=atsource%imm_glob
  end subroutine copy_config

  subroutine dealloc_atom_config(atconf)
    class(atom_config), intent(inout)::atconf
    atconf%im=0 ; atconf%imm=0 ; atconf%im_glob=0; atconf%imm_glob=0
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
       deallocate(atconf%vp)
    end if
  end subroutine dealloc_atom_config_d
  subroutine dealloc_atom_config_e(atconf)
    class(atom_config_e), intent(inout)::atconf
    call atconf%atom_config_d%dealloc
    if(allocated(atconf%eat))deallocate(atconf%eat)
    if(allocated(atconf%xpp))deallocate(atconf%xpp)
    if(allocated(atconf%sigat))deallocate(atconf%sigat)
    if(allocated(atconf%glangv))deallocate(atconf%glangv)
    if(allocated(atconf%ax))deallocate(atconf%ax)

  end subroutine dealloc_atom_config_e
  subroutine dealloc_atom_config_arps(atconf)
    class(atom_config_arps), intent(inout)::atconf
    call atconf%atom_config_e%dealloc
    if(allocated(atconf%rho))deallocate(atconf%rho)
    if(allocated(atconf%mov))deallocate(atconf%mov)
    if(allocated(atconf%fpr))deallocate(atconf%fpr)
    if(allocated(atconf%fpg))deallocate(atconf%fpg)

  end subroutine dealloc_atom_config_arps

  subroutine deftype(atsource,atcible)  ! initialize atcible to the type of atsource, inluding the values of lax, lpreeat, etc.
    class(atom_config),intent(in)::atsource
    class(atom_config),allocatable::atcible

       select type (atsource)
       class is (atom_config)
          allocate(atom_config::atcible)
       class is (atom_config_d)
          allocate(atom_config_d::atcible)
       class is (atom_config_e)
          allocate(atom_config_e::atcible)
       class is (atom_config_arps)
          allocate(atom_config_arps::atcible)
       end select
       call atsource%Eegal(atcible)
     end subroutine deftype
     
  !on inverse deux atomes dans la configuration
  subroutine switch_atom(atsource,ind_switch_1, ind_switch_2)
    class(atom_config)::atsource
    class(atom_config),allocatable:: intermediaire
    integer :: ind_switch_1, ind_switch_2,nag1,nag2
    call atsource%deftype(intermediaire)
    call intermediaire%init(imin=2)
    nag1=atsource%num_at_glob(ind_switch_1)
    nag2=atsource%num_at_glob(ind_switch_2)
    call atsource%copy_atom(ind_switch_1,intermediaire,1)
    call atsource%copy_atom(ind_switch_2,intermediaire,2)
    call intermediaire%copy_atom(2,atsource,ind_switch_1)
    call intermediaire%copy_atom(1,atsource,ind_switch_2)
    atsource%num_at_glob(ind_switch_1)=nag2
    atsource%num_at_glob(ind_switch_2)=nag1
  end subroutine




  
#ifdef PARA  
  
#endif
  subroutine pack(at2pack,imm_in)
    class(atom_config),intent(inout):: at2pack
    integer,optional, intent(in):: imm_in
    class(atom_config),allocatable::at
!!$    type(atom_config)::at
!!$    type(atom_config_d)::atd
!!$    type(atom_config_e)::ate
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
          call arret_ndm
       end if
    else
       immn=imn
    end if
    if (at2pack%ltabvois)then
       nvois=at2pack%nvois; rvois=at2pack%rvois
    else
       nvois=0;rvois=0
    end if

    call at2pack%deftype(at)
    call at%init(imn,immn,at2pack%ltabvois,nvois,rvois,im_glob=at2pack%im_glob,imm_glob=at2pack%imm_glob)
    do i=1,immn
       call at2pack%copy_atom(i,at,i)
    end do
    call at2pack%dealloc
    call at%copy_config(at2pack,lrescl=.true.)
    

  end subroutine pack

  subroutine fab(atsource,atcible,lback,lrescl,commsp) ! construit atsource à partir de lgul de atcible , ecrase atcible
    class(atom_config),intent(in)::atsource
    class(atom_config),intent(out)::atcible
    logical::lback
    logical, optional::lrescl
    type(mpi_communicator),optional::commsp
    logical::lrescale=.true.
    integer::i2,imtrf,i,immtrf
    integer::nvois
    real(double)::rvois

    if (present(lrescl))lrescale=lrescl
!    write(6,*)'FAB lrescale',lrescale
    if (lrescale) then
       call atcible%dealloc 
       if (atsource%ltabvois)then
          nvois=atsource%nvois; rvois=atsource%rvois
          
       else
          nvois=0;rvois=0
       end if
       imtrf=COUNT(atsource%lgul(1:atsource%im))
       if (imtrf==0) then
          immtrf=1
       else
          immtrf=imtrf
       end if
!       write(6,*)'FAB',imtrf
       call atsource%Eegal(atcible)
       call atcible%init(imtrf,immtrf,atsource%ltabvois,nvois,rvois)
    end if

    call atcible%zero
    i2=0
    do i=1,atsource%im
       if(atsource%lgul(i)) then
          i2=i2+1
          call atsource%copy_atom(i,atcible,i2,lextend=.false.)
          if (lback) atcible%num_at_glob(i2)=i
       end if
    end do
    if (i2.ne.imtrf) then
       write(6,*)'WTF ?'
       call arret_ndm
    end if
       atcible%im_glob=atcible%im
#ifdef PARA
    if (present(commsp)) then

!       atcible%imm_glob=atcible%imm
       call commsp%sum(atcible%im_glob)

    end if
#endif
    

  end subroutine fab



  subroutine backto(atfab,atback)
    class(atom_config),intent(in)::atfab
    class(atom_config)::atback
    integer::i2,imtrf,i
    do i2=1,atfab%im
       i=atfab%num_at_glob(i2)
       imtrf=atback%num_at_glob(i)
       call atfab%copy_atom(i2,atback,i,lextend=.false.)
       atback%num_at_glob(i)=imtrf
    end do
  end subroutine backto

  
  subroutine sort (atsource,atcible) ! construit atsource à partir de lgul de atcible , ecrase atcible
    class(atom_config),intent(in)::atsource
    class(atom_config),intent(out)::atcible
    integer::i2,imtrf,i,immtrf,i3
    integer::nvois
    real(double)::rvois
    logical::lback=.true.
    call atcible%dealloc 
    if (atsource%ltabvois)then
       nvois=atsource%nvois; rvois=atsource%rvois
    else
       nvois=0;rvois=0
    end if
    immtrf=atsource%im
    imtrf=COUNT(atsource%lgul(1:atsource%im))
    !    call atsource%Eegal(atcible)
    write(6,*)'SORT',imtrf,immtrf
    call atcible%init(imtrf,immtrf,atsource%ltabvois,nvois,rvois)
    call atcible%zero
    i2=0;i3=atcible%im
    do i=1,atsource%im
       if(atsource%lgul(i)) then
          i2=i2+1
          call atsource%copy_atom(i,atcible,i2,lextend=.false.)
          if (lback) atcible%num_at_glob(i2)=i
       else
          i3=i3+1
          call atsource%copy_atom(i,atcible,i3,lextend=.false.)
          if (lback) atcible%num_at_glob(i3)=i
       end if
    end do
    if (i2.ne.imtrf) then
       write(6,*)'SORT WTF2 ?'
       call arret_ndm
    end if
    if (i3.ne.immtrf) then
       write(6,*)'SORT WTF3 ?'
       call arret_ndm
    end if

  end subroutine sort

    
  subroutine add2conf (atsource,atcible,lextend,ldealloc)
    class(atom_config),intent(inout)::atsource
    class(atom_config),intent(inout)::atcible
    logical,optional,intent(in)::ldealloc,lextend



    logical :: ldal ! par defaut pas de destruction de atsource
    logical :: lext ! par defaut on etend atcible si besoin

    integer::i2,i,imcib,imnew,immcib,imsrc,immsrc,immnew
    integer::nvois
    real(double)::rvois
    ldal=.false.
    lext=.true.


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
       call atcible%init(imsrc,immsrc,atcible%ltabvois,nvois,rvois)
    else
       if (lext) then
          if (atcible%imm.lt.immnew)   call atcible%extend(immnew)
       else
          if (atcible%imm.lt.imnew)  then
             write(6,*)'addition de atsource a atcible pas possible'
             call arret_ndm
          end if
       end if
       atcible%im=imsrc+imcib
       atcible%imm=immsrc+immcib
    end if


    do i=1,atsource%im
       i2=imcib+i
       call atsource%copy_atom(i,atcible,i,lextend=.false.)
     end do

    if (ldal) call atsource%dealloc
  end subroutine add2conf



  subroutine print(atin,i1,i2,unit,natg1,natg2,caracT,mess,iter)
    class(atom_config), intent(in)::atin
    integer,optional::i1,i2,unit,natg1,natg2,iter
    integer::i1l,i2l,natg1l,natg2l
    character(len=*),optional,intent(in)::caracT
    character(len=*),optional::mess
    character(len=26)::carac
    integer::i,ig,iprt,unitw,imp
    unitw=6
    if (.not.present(caracT)) then
       carac='xfniewdlpvfrugasm'
    else
       carac=caracT
    end if
    if (present(unit))unitw=unit
    
    if (present(mess)) then
       write(unitw,*)'atomPRINT ',mess,rang
    end if
    if (present(iter)) then
       write(unitw,*)'ITERATION ',iter
    end if

    
#ifdef PARA
    if ((lspacendm).and.(nprocspace.gt.1)) then
       if (atin%imf.ne.0) then
          imp=atin%imf
       else
          imp=atin%im
       end if
       
    else
       imp=atin%im
    end if
!    write(6,*)'AAAAA',rang,atin%im,atin%imf
#else
    imp=atin%im
#endif

    write(unitw,*)'im = ',atin%im
    
    write(unitw,*)'imm = ',atin%imm
    write(unitw,*)'im_glob = ',atin%im_glob
    write(unitw,*)'imm_glob = ',atin%imm_glob
#ifdef PARA
    write(unitw,*)'imf = ',atin%imf
#endif
    
    write(unitw,*)'icaltabt = ',atin%icaltabt
    write(unitw,*)'ltabvois ', atin%ltabvois

    natg1l=1
    natg2l=maxval(atin%num_at_glob)
    if(present(natg1))natg1l=natg1
    if(present(natg2))natg2l=natg2
    i1l=1
    i2l=imp
    if(present(i1))i1l=i1
    if(present(i2))i2l=i2
    write(unitw,*)'borders',i1l,i2l,natg1l,natg2l
    do i=1,imp
       
       if((i.lt.i1l).or.(i.gt.i2l)) cycle
       if((atin%num_at_glob(i).lt.natg1l).or.(atin%num_at_glob(i).gt.natg2l)) cycle
       write(unitw,*)
       if(scan('x',carac).ne.0)then
          write(unitw,'(A,2i9,3E15.7)')'%xp= ', i,atin%num_at_glob(i),atin%xp(:,i)
       end if
       if(scan('i',carac).ne.0)then
          write(unitw,'(A,2i9,I3)')'%ityp= ', i,atin%num_at_glob(i),atin%ityp(i)
    end if
    if(scan('n',carac).ne.0)then
          write(unitw,*)'%num_at_glob= ', i,atin%num_at_glob(i)
    end if
#ifdef PARA
    if(scan('p',carac).ne.0)then
          write(unitw,*)'%proc_at= ', i,atin%num_at_glob(i),atin%proc_at(i)
    end if
#endif    
       
       if(scan('f',carac).ne.0)then
             write(unitw,'(A,2i9,3E15.7)')'%fp= ', i,atin%num_at_glob(i),atin%fp(:,i)
       end if
       if(scan('e',carac).ne.0)then
             write(unitw,*)'%ielat= ', i,atin%num_at_glob(i),atin%ielat(i)
       end if
       if(scan('l',carac).ne.0)then
             write(unitw,*)'%lgul= ', i,atin%num_at_glob(i),atin%lgul(i)
       end if

       select type (atin)
          class is (atom_config_d)
             write(unitw,*)'prt_d'
             if(scan('v',carac).ne.0)then
                   write(unitw,'(A,2i9,3E15.7)')'%vp= ', i,atin%num_at_glob(i),atin%vp(:,i)
             end if
          class is (atom_config_e)
             write(unitw,*)'prt_e'
             if(scan('v',carac).ne.0)then
                   write(unitw,'(A,2i9,3E15.7)')'%vp= ', i,atin%num_at_glob(i),atin%vp(:,i)
             end if
             if (atin%lxpp) then
                if(scan('r',carac).ne.0)then
                   write(unitw,'(A,2i9,3E15.7)')'%xpp= ', i,atin%num_at_glob(i),atin%xpp(:,i)
                end if
                
             end if

             if (atin%lsigat) then
                if(scan('g',carac).ne.0)then
                      write(unitw,'(A,2i9,9E15.7)')'%sigat= ',i,atin%num_at_glob(i), atin%sigat(:,:,i)
                end if
             end if
          if (atin%lprteat) then
             if(scan('u',carac).ne.0)then
                   write(unitw,*)'%eat= ', i,atin%num_at_glob(i),atin%eat(i)
             end if
          end if
          if (atin%lax) then
             if(scan('a',carac).ne.0)then
                   write(unitw,'(A,2i9,3E15.7)')'%ax= ', i,atin%num_at_glob(i),atin%ax(:,i)
             end if
          end if
       end select

       select type (atin)
          class is (atom_config_arps)
             write(unitw,*)'prt_arps'
             if(scan('g',carac).ne.0)then
                   write(unitw,'(A,2i9,I4)')'%mov= ', i,atin%num_at_glob(i),atin%mov(i)
                   write(unitw,'(A,2i9,3E15.7)')'%fpr= ', i,atin%num_at_glob(i),atin%fpr(:,i)
                   if (allocated(atin%fpg))write(unitw,'(A,2i9,3E15.7)')'%fpg= ', i,atin%num_at_glob(i),atin%fpg(:,i)
                   if (allocated(atin%rho))write(unitw,'(A,2i9,1E15.7)')'%rho= ', i,atin%num_at_glob(i),atin%rho(i)
             end if
          end select
       end do
       if (atin%ltabvois) then
          do i=1,atin%im
             write(unitw,*)'%iwmax= ', i,atin%num_at_glob(i),atin%iwmax(i)
          end do
       end if
       write(unitw,*)
    flush(unitw)
  end subroutine print
!MANQUE SIG AU MINIMUM
!!$
!!$  subroutine print(atin,i1,i2,unit,natg1,natg2,caracT)
!!$    class(atom_config), intent(in)::atin
!!$    integer,optional::i1,i2,unit,natg1,natg2
!!$    character(len=*),optional,intent(in)::caracT
!!$    character(len=26)::carac
!!$    integer::i,im,ifin,ideb,ist,ifn,natpr,ig,iprt,unitw,imp
!!$    class (atom_config),allocatable::atprt
!!$    unitw=6
!!$    if (.not.present(caracT)) then
!!$       carac='xfniewdlpvfrugasm'
!!$    else
!!$       carac=caracT
!!$    end if
!!$    if (present(unit))unitw=unit
!!$#ifdef PARA
!!$    if ((lspacendm).and.(nprocspace.gt.1)) then
!!$       imp=atin%imf
!!$    else
!!$       imp=atin%im
!!$    end if
!!$    write(6,*)'AAAAA',rang,atin%im,atin%imf
!!$#else
!!$    imp=atin%im
!!$#endif
!!$    call atin%deftype(atprt)
!!$    if (present(natg1)) then
!!$       if (present(natg2)) then
!!$          natpr=natg2-natg1+1
!!$          if (natpr.lT.1) stop
!!$       else
!!$          natg2=natg1
!!$       end if
!!$       iprt=0
!!$       do ig=natg1,natg2
!!$          do i=1,imp
!!$             if(atin%num_at_glob(i)==ig) then
!!$                iprt=iprt+1
!!$                if (iprt.gt.natpr) then
!!$                   write(6,*)'num_at_glob multiples ?'
!!$                   call arret_ndm
!!$                end if
!!$             end if
!!$          end do
!!$       end do
!!$       natpr=iprt 
!!$       call atprt%init(natpr,rvois=0.d0)
!!$!      im=atprt%im
!!$       iprt=0
!!$       do ig=natg1,natg2
!!$          do i=1,imp!atin%im
!!$             if(atin%num_at_glob(i)==ig) then
!!$                iprt=iprt+1
!!$                if (iprt.gt.natpr) then
!!$                   write(6,*)'num_at_glob multiples ?'
!!$                   call arret_ndm
!!$                end if
!!$                call atin%copy_atom(i,atprt,iprt)
!!$             end if
!!$          end do
!!$       end do
!!$!       im=atprt%im
!!$       iprt=0
!!$       do ig=natg1,natg2
!!$          do i=1,imp!atin%im
!!$             if(atin%num_at_glob(i)==ig) then
!!$                iprt=iprt+1
!!$                if (iprt.gt.natpr) then
!!$                   write(6,*)'num_at_glob multiples ?'
!!$                   call arret_ndm
!!$                end if
!!$                call atin%copy_atom(i,atprt,iprt)
!!$             end if
!!$          end do
!!$       end do
!!$!       im=atprt%im
!!$    else
!!$
!!$       natpr=atin%imf
!!$       call atprt%init(natpr,rvois=0.d0)
!!$       do i=1,atin%imf
!!$          call atin%copy_atom(i,atprt,i)
!!$       end do
!!$
!!$!       call atin%copy_config(atprt,lrescl=.true.)
!!$!       im=atin%im
!!$
!!$
!!$
!!$    end if
!!$          
!!$    atprt%im_glob=atin%im_glob
!!$    atprt%imm_glob=atin%imm_glob
!!$    
!!$    write(unitw,*)'im = ',atprt%im
!!$    
!!$    write(unitw,*)'imm = ',atin%imm
!!$    write(unitw,*)'im_glob = ',atprt%im_glob
!!$    write(unitw,*)'imm_glob = ',atprt%imm_glob
!!$#ifdef PARA
!!$    write(unitw,*)'imf = ',atin%imf
!!$#endif
!!$    
!!$    write(unitw,*)'icaltabt = ',atprt%icaltabt
!!$    write(unitw,*)'ltabvois ', atprt%ltabvois
!!$!    write(6,*)
!!$    ideb=1
!!$    ifin=atprt%im
!!$    if (present(i1))then
!!$       ideb=i1
!!$    else
!!$       ideb=1
!!$    endif
!!$    if (present(i2)) then
!!$       ifin=i2
!!$    else
!!$       if (present(i1))then
!!$          ifin=i1
!!$       else
!!$          ifin=atprt%im
!!$       endif
!!$    end if
!!$    if (allocated(atprt%xp)) then
!!$    if(scan('x',carac).ne.0)then
!!$       do i=ideb,ifin
!!$          write(unitw,'(A,2i9,3E15.7)')'%xp= ', i,atprt%num_at_glob(i),atprt%xp(:,i)
!!$       end do
!!$    end if
!!$    if(scan('i',carac).ne.0)then
!!$       do i=ideb,ifin
!!$          write(unitw,'(A,2i9,I3)')'%ityp= ', i,atprt%num_at_glob(i),atprt%ityp(i)
!!$       end do
!!$    end if
!!$    if(scan('n',carac).ne.0)then
!!$       do i=ideb,ifin
!!$          write(unitw,*)'%num_at_glob= ', i,atprt%num_at_glob(i)
!!$       end do
!!$    end if
!!$#ifdef PARA
!!$    if(scan('p',carac).ne.0)then
!!$       do i=ideb,ifin
!!$          write(unitw,*)'%proc_at= ', i,atprt%num_at_glob(i),atprt%proc_at(i)
!!$       end do
!!$    end if
!!$#endif    
!!$       
!!$       if(scan('f',carac).ne.0)then
!!$          do i=ideb,ifin
!!$             write(unitw,'(A,2i9,3E15.7)')'%fp= ', i,atprt%num_at_glob(i),atprt%fp(:,i)
!!$          end do
!!$       end if
!!$       if(scan('e',carac).ne.0)then
!!$          do i=ideb,ifin
!!$             write(unitw,*)'%ielat= ', i,atprt%num_at_glob(i),atprt%ielat(i)
!!$          end do
!!$       end if
!!$       select type (atprt)
!!$          class is (atom_config_d)
!!$             write(unitw,*)'prt_d'
!!$             if(scan('v',carac).ne.0)then
!!$                do i=ideb,ifin
!!$                   write(unitw,'(A,2i9,3E15.7)')'%vp= ', i,atprt%num_at_glob(i),atprt%vp(:,i)
!!$                end do
!!$             end if
!!$             if(scan('r',carac).ne.0)then
!!$                do i=ideb,ifin
!!$                   write(unitw,'(A,2i9,3E15.7)')'%xpp= ', i,atprt%num_at_glob(i),atprt%xpp(:,i)
!!$                end do
!!$             end if
!!$          class is (atom_config_e)
!!$             write(unitw,*)'prt_e'
!!$             if(scan('v',carac).ne.0)then
!!$                do i=ideb,ifin
!!$                   write(unitw,'(A,2i9,3E15.7)')'%vp= ', i,atprt%num_at_glob(i),atprt%vp(:,i)
!!$                end do
!!$             end if
!!$             if(scan('r',carac).ne.0)then
!!$                do i=ideb,ifin
!!$                   write(unitw,'(A,2i9,3E15.7)')'%xpp= ', i,atprt%num_at_glob(i),atprt%xpp(:,i)
!!$                end do
!!$             end if
!!$
!!$
!!$             if (atprt%lsigat) then
!!$                if(scan('g',carac).ne.0)then
!!$                   do i=ideb,ifin
!!$                      write(unitw,'(A,2i9,9E15.7)')'%sigat= ',i,atprt%num_at_glob(i), atprt%sigat(:,:,i)
!!$                   end do
!!$                end if
!!$             end if
!!$          if (atprt%lprteat) then
!!$             if(scan('u',carac).ne.0)then
!!$                do i=ideb,ifin
!!$                   write(unitw,*)'%eat= ', i,atprt%num_at_glob(i),atprt%eat(i)
!!$                end do
!!$             end if
!!$          end if
!!$          if (atprt%lax) then
!!$             if(scan('a',carac).ne.0)then
!!$                do i=ideb,ifin
!!$                   write(unitw,'(A,2i9,3E15.7)')'%ax= ', i,atprt%num_at_glob(i),atprt%ax(:,i)
!!$                end do
!!$             end if
!!$          end if
!!$       end select
!!$
!!$       select type (atprt)
!!$          class is (atom_config_arps)
!!$             write(unitw,*)'prt_arps'
!!$             if(scan('g',carac).ne.0)then
!!$                do i=ideb,ifin
!!$                   write(unitw,'(A,2i9,I4)')'%mov= ', i,atprt%num_at_glob(i),atprt%mov(i)
!!$                end do
!!$                do i=ideb,ifin
!!$                   write(unitw,'(A,2i9,3E15.7)')'%fpr= ', i,atprt%num_at_glob(i),atprt%fpr(:,i)
!!$                end do
!!$                if (allocated(atprt
!!$             end if
!!$             if(scan('r',carac).ne.0)then
!!$                do i=ideb,ifin
!!$                   write(unitw,'(A,2i9,3E15.7)')'%xpp= ', i,atprt%num_at_glob(i),atprt%xpp(:,i)
!!$                end do
!!$             end if
!!$
!!$       if (atprt%ltabvois) then
!!$          do i=ideb,imp
!!$             write(unitw,*)'%iwmax= ', i,atprt%num_at_glob(i),atprt%iwmax(i)
!!$          end do
!!$
!!$          if (allocated(atprt%indi))then
!!$             if (ideb==1) then
!!$                ist=1
!!$             else
!!$                ist=atprt%iwmax(ideb-1)+1
!!$             end if
!!$             if (ifin==im) then
!!$                ifn=size(atprt%indi)
!!$             else
!!$                ifn=atprt%iwmax(ifin)
!!$             end if
!!$
!!$             do i=ist,ifn,100
!!$                write(unitw,*)'indi', i,atprt%indi(i)
!!$             end do
!!$          end if
!!$       end if
!!$    end if
!!$    flush(unitw)
!!$  end subroutine print
!!$!MANQUE SIG AU MINIMUM
!!$

  
  subroutine vers_master_atom(atcfloc,atcfcomp,div,caracT)
    class(atom_config),intent(in)::atcfloc
    class(atom_config)::atcfcomp
    type(para_config)::div
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac

#ifdef PARA
    integer::idmaster,idloc ! proc master,proclocal 

    integer::sizeI,sizeR,sizel,csi,csr,csl
    integer, dimension (0:26):: Iposf,Rposf,Lposf
    integer,allocatable:: ibuffer(:)
    logical,allocatable::lbuffer(:)
    real(double),allocatable::rbuffer(:)
    integer:: cst(3)
    type(mpi_communicator)::mpic

    integer::natgM
    integer,allocatable::inag(:)

    integer,allocatable::nag(:)
    integer::imloc,imloc3,iproc,imrecv,icomp,imloc9
    integer::imtot,proc_source,npim,iloc
    logical,allocatable::mask(:)
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugasm'
    else
       carac=caracT//'np'
    end if
    idmaster=0
    idloc=div%mpi_image%rank
    npim=div%mpi_image%nproc
    mpic=div%mpi_image

    natgM=maxval(atcfloc%num_at_glob(1:atcfloc%im))
    call mpic%max(natgM)
    allocate(inag(natgM))
    inag=0
    if (idloc==idmaster) then
       imtot=0
       do iproc=0,npim-1
          if (iproc==idmaster) then
             imrecv=atcfloc%im
             imtot=imtot+imrecv
          else
             if (allocated(ibuffer)) then
                deallocate(rbuffer);deallocate(ibuffer);deallocate(lbuffer)
             end if
             call  mpic%probe(11011,sourceout=proc_source)
             call mpic%RECV(imrecv,proc_source, 11011)
             imtot=imtot+imrecv
          end if
       end do
       if (imtot.ne.atcfcomp%im) then
          write(6,*)'atomes perdus 1?',idloc, imtot,atcfcomp%im,div%mpi_orig%rank
          call MPI_finalize(ierr)
          call arret_ndm
       end if

       if (all(atcfcomp%num_at_glob(1:atcfcomp%im)==0)) then  ! This is anew atcfcomp with undefined atcfcomp :inag points to -1 to show that
          inag(1:natgM)=-1
       else if (any(atcfcomp%num_at_glob(1:atcfcomp%im)==0)) then !This is bulsshit (neither new nor pre-existing) smells like inconsistency
          write(6,*)'VERS MASTER au moins un NAG nul'
          call arret_ndm
       else ! This a return to an existing atcfcomp which has its own num_at_glob numbering
          do icomp=1,atcfcomp%im
             !          write(6,*)'atcfcomp',icomp,atcfcomp%num_at_glob(icomp)
             inag(atcfcomp%num_at_glob(icomp))=icomp
          end do
       end if
    else
       imloc=atcfloc%im
       call mpic%SEND(imloc,idmaster,11011)
    end if
       

    
    if (idloc==idmaster) then
       icomp=0
       imtot=0
       do iproc=0,npim-1

          ! Pour le processeur maitre il n'y a rien a faire
          ! reception des donnees des autres processeurs
          if (iproc==idmaster) then
             imrecv=atcfloc%im
             imtot=imtot+imrecv
             allocate(nag(imrecv))
             nag(1:imrecv)=atcfloc%num_at_glob(1:imrecv)
             do iloc=1,imrecv
                if (inag(nag(iloc))==-1) then ! num_at_glob pas defini pour atcfcomp
                   icomp=icomp+1
                   inag(nag(iloc))=icomp
                   atcfcomp%num_at_glob(icomp)=nag(iloc)
                else
                   icomp=inag(nag(iloc))
                end if
!                write(6,*)'L2M',iproc,iloc,nag(iloc)
                if (nag(iloc).ne.atcfcomp%num_at_glob(icomp))then
                   write(6,*)'erreur NATG',iloc,icomp,nag(iloc),atcfcomp%num_at_glob(icomp)
                   call MPI_finalize(ierr)
                   call arret_ndm
                end if
                if(scan('e',carac).ne.0)   atcfcomp%ielat(icomp)=atcfloc%ielat(iloc)
                if(scan('x',carac).ne.0)   atcfcomp%xp(1:3,icomp)=atcfloc%xp(1:3,iloc)
                if(scan('f',carac).ne.0)   atcfcomp%fp(1:3,icomp)=atcfloc%fp(1:3,iloc)
                if(scan('i',carac).ne.0)   atcfcomp%ityp(icomp)=atcfloc%ityp(iloc)
                if( allocated(atcfloc%lgul)) then
                   if(scan('l',carac).ne.0)  atcfcomp%lgul(icomp)=atcfloc%lgul(iloc)
                end if
                 if(scan('p',carac).ne.0) atcfcomp%proc_at(icomp)=idmaster 
                 if (atcfcomp%ltabvois) then 
                    if(scan('w',carac).ne.0)atcfcomp%iwmax(icomp)=atcfloc%iwmax(iloc)
                 end if
                 select type(atcfloc)
                 class is (atom_config_d)
                    select type (atcfcomp)
                    class is (atom_config_d)
                       if(scan('v',carac).ne.0) atcfcomp%vp(1:3,icomp)=atcfloc%vp(1:3,iloc)
                    end select
                 end select
                 select type(atcfloc)
                 class is (atom_config_e)
                    select type (atcfcomp)
                    class is (atom_config_e)
                       if((atcfcomp%lxpp).and.(atcfloc%lxpp))then
                          if(scan('r',carac).ne.0) atcfcomp%xpp(:,icomp)=atcfloc%xpp(:,iloc)
                       endif
                       if((atcfcomp%lsigat).and.(atcfloc%lsigat))then
                          if(scan('s',carac).ne.0) atcfcomp%sigat(1:3,1:3,icomp)=atcfloc%sigat(1:3,1:3,iloc)
                       endif
                       if((atcfcomp%lprteat).and.(atcfloc%lprteat))then
                          if(scan('u',carac).ne.0) atcfcomp%eat(icomp)=atcfloc%eat(iloc)
                       endif
                       if((atcfcomp%lax).and.(atcfloc%lax))then
                          if(scan('a',carac).ne.0) atcfcomp%ax(:,icomp)=atcfloc%ax(:,iloc)
                       endif
                       if((atcfcomp%llangevin).and.(atcfloc%llangevin))then
                          if(scan('g',carac).ne.0) atcfcomp%glangv(:,icomp)=atcfloc%glangv(:,iloc)
                       endif
                    end select
                 end select
                 select type(atcfloc)
                 class is (atom_config_arps)
                    select type (atcfcomp)
                    class is (atom_config_arps)
                       if(scan('m',carac).ne.0)  atcfcomp%mov(icomp)=atcfloc%mov(iloc)
                    end select
                 end select

             end do
             deallocate(nag)

          else
             if (allocated(ibuffer)) then
                deallocate(rbuffer);deallocate(ibuffer);deallocate(lbuffer)
             end if
             call  mpic%probe(11001,sourceout=proc_source)
             call mpic%RECV(imrecv,proc_source, 11001)
             call mpic%recv (cst,proc_source,312)
             imtot=imtot+imrecv
             allocate(nag(imrecv))
             call mpic%RECV(nag(1:imrecv), proc_source,10004)

             call buffersizes (atcfloc,sizeI,sizeR,sizel,IposF,Rposf,Lposf,carac,imrecv)
             if (cst(1).ne.sizeI) then
                write(6,*)'erreur CST1B ',sizeI,cst(1)
                call arret_ndm
             end if
             if (cst(2).ne.sizel) then
                write(6,*)'erreur CST2B ',sizel,cst(2)
                call arret_ndm
             end if
             if (cst(3).ne.sizeR) then
                write(6,*)'erreur CST3B ',sizer,cst(3)
                call arret_ndm
             end if
             allocate(Rbuffer(sizeR));allocate(ibuffer(sizeI));allocate(lbuffer(sizeL))
             if (cst(1).ne.0)call mpic%recv(ibuffer,proc_source,314)
             if (cst(2).ne.0)call mpic%recv(lbuffer,proc_source,315)
             if (cst(3).ne.0)call mpic%recv(Rbuffer,proc_source,316)
             do iloc=1,imrecv
                if (inag(nag(iloc))==-1) then ! num_at_glob pas defini pour atcfcomp
                   icomp=icomp+1
                   inag(nag(iloc))=icomp
                   atcfcomp%num_at_glob(icomp)=nag(iloc)
                else
                   icomp=inag(nag(iloc))
                end if
                if (nag(iloc).ne.atcfcomp%num_at_glob(icomp))then
                   write(6,*)'erreur NATG2',iloc,icomp,nag(iloc),atcfcomp%num_at_glob(icomp)
                   call MPI_finalize(ierr)
                   call arret_ndm
                end if
             end do
             call distribnag(nag,inag,imrecv,atcfcomp,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf)
             do iloc=1,imrecv
                atcfcomp%proc_at(inag(nag(iloc)))=proc_source !  ibuffer(1:imrecv)
             end do

             deallocate(nag)
             

          end if
       end do


       if (imtot.ne.atcfcomp%im) then
          write(6,*)'atomes perdus 3?',idloc, imtot,atcfcomp%im,div%mpi_orig%rank
          call MPI_finalize(ierr)
          call arret_ndm
       end if
    else
       imloc=atcfloc%im;imloc3=3*imloc; imloc9=3*imloc3
       allocate(mask(atcfloc%imm))
       mask=.false.
       mask(1:atcfloc%im)=.true.
       call buffersizes (atcfloc,sizeI,sizeR,sizel,IposF,Rposf,Lposf,carac,imloc) ! ns car on transfère seulement ns atomes vers atcfloc
       allocate(ibuffer(sizeI));     allocate(Lbuffer(sizeL));     allocate(Rbuffer(sizeR));
       call buildbuff(atcfloc,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf,mask,imloc)
       csT(1)=csi
       csT(2)=csl
       csT(3)=csR
       call mpic%SEND(imloc,idmaster,11001)
       call mpic%send (cst,idmaster,312)
       allocate(nag(imloc))
       nag(1:imloc)=atcfloc%num_at_glob(1:imloc)
       call mpic%SEND(nag, idmaster, 10004)
       deallocate(nag)
       if (csi.ne.0)call mpic%send(ibuffer,idmaster,314)
       if (csl.ne.0)call mpic%send(lbuffer,idmaster,315)
       if (csR.ne.0)call mpic%send(Rbuffer,idmaster,316)
    end if
    call mpic%barrier
#endif       
    return
  end subroutine vers_master_atom



  
  subroutine master2loc_atom(atcfcomp,atcfloc,div,caracT)

    class(atom_config)::atcfloc
    class(atom_config)::atcfcomp
    type(para_config)::div
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at
    ! v=vp,r=xpp
    ! u=eat,g=glangv;a=ax;s=sigat
#ifdef PARA
    integer::idmaster,idloc,npim ! proc master,proclocal ,communicateur
    integer::iproc,imrecv
    integer::imcomp,imtot,ns,iloc,i
    logical,allocatable::mask(:)

    integer::sizeI,sizeR,sizel,csi,csr,csl
    integer, dimension (0:26):: Iposf,Rposf,Lposf
    integer,allocatable:: ibuffer(:)
    logical,allocatable::lbuffer(:)
    real(double),allocatable::rbuffer(:)
    integer:: cst(3)
    type(mpi_communicator)::mpic

    idmaster=0
    idloc=div%mpi_image%rank
    npim=div%mpi_image%nproc
    mpic=div%mpi_image
    imcomp=atcfcomp%im


    if (.not.present(caracT)) then
       carac='xfniewdlpvrugasm'
    else
       carac=caracT//'np'
       
    end if



    if (idloc==idmaster) then
       allocate(mask(atcfcomp%imm))

       imtot=0
       do iproc=0,npim-1
          mask(:)=.false.
          mask(1:imcomp)=(atcfcomp%proc_at(1:imcomp)==iproc)
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
                      end select
                   end select
                   select type(atcfloc)
                   class is (atom_config_e)
                      select type (atcfcomp)
                      class is (atom_config_e)
                         if((atcfcomp%lxpp).and.(atcfloc%lxpp))then
                            atcfloc%xpp(:,iloc)=atcfcomp%xpp(:,i)
                         end if
                         if((atcfcomp%lsigat).and.(atcfloc%lsigat))then
                            atcfloc%sigat(:,:,iloc)=atcfcomp%sigat(:,:,i)
                         endif
                         if((atcfcomp%lprteat).and.(atcfloc%lprteat))then
                            atcfloc%eat(iloc)=atcfcomp%eat(i)
                         endif
                         if((atcfcomp%lax).and.(atcfloc%lax))then
                            atcfloc%ax(:,iloc)=atcfcomp%ax(:,i)
                         endif
                         if((atcfcomp%llangevin).and.(atcfloc%llangevin))then
                            atcfloc%glangv(:,iloc)=atcfcomp%glangv(:,i)
                         endif
                      end select
                   end select
                   select type(atcfloc)
                   class is (atom_config_arps)
                      select type (atcfcomp)
                      class is (atom_config_arps)
                         atcfloc%mov(iloc)=atcfcomp%mov(i)
                      end select
                   end select
                      
                end if
             end do
             atcfloc%im=ns
          else
             if (allocated(ibuffer)) then
                deallocate(ibuffer);deallocate(rbuffer);deallocate(lbuffer)
             end if
             call buffersizes (atcfcomp,sizeI,sizeR,sizel,IposF,Rposf,Lposf,carac,ns) ! ns car on transfère seulement ns atomes vers atcfloc
             allocate(ibuffer(sizeI));     allocate(Lbuffer(sizeL));     allocate(Rbuffer(sizeR));
             call buildbuff(atcfcomp,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf,mask,ns)
             csT(1)=csi
             csT(2)=csl
             csT(3)=csR
             call mpic%SEND(ns,iproc,20001)
             call mpic%send (cst,iproc,212)
             if (csi.ne.0)call mpic%send(ibuffer,iproc,214)
             if (csl.ne.0)call mpic%send(lbuffer,iproc,215)
             if (csR.ne.0)call mpic%send(Rbuffer,iproc,216)

          end if
       end do

    else
     
       call mpic%recv(imrecv, idmaster, 20001)
       call mpic%recv (cst,idmaster,212)
       call buffersizes (atcfloc,sizeI,sizeR,sizel,IposF,Rposf,Lposf,carac,atcfloc%im) ! pas de ns car on reçoit tous les atomes de atloc
       if (atcfloc%im.ne.imrecv) then
          write(6,*)'ERREUR M2L', atcfloc%im,imrecv
          call arret_ndm
       end if
       if (cst(1).ne.sizeI) then
          write(6,*)'erreur CST1B ',sizeI,cst(1)
          call arret_ndm
       end if
       if (cst(2).ne.sizel) then
          write(6,*)'erreur CST2B ',sizel,cst(2)
          call arret_ndm
       end if
       if (cst(3).ne.sizeR) then
          write(6,*)'erreur CST3B ',sizer,cst(3)
          call arret_ndm
       end if
       
       allocate(Rbuffer(sizeR));allocate(ibuffer(sizeI));allocate(lbuffer(sizeL))
       if (cst(1).ne.0)call mpic%recv(ibuffer,idmaster,214)
       if (cst(2).ne.0)call mpic%recv(lbuffer,idmaster,215)
       if (cst(3).ne.0)call mpic%recv(Rbuffer,idmaster,216)
       call copybuff (atcfloc,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf,atcfloc%im)

    end if

    return
#endif

  end subroutine master2loc_atom
 

  subroutine buffersizes (atcf,sizeI,sizeR,sizel,IposF,Rposf,Lposf,carac,nmask)
    integer, intent(out)::sizeI,sizel,sizer
    integer, intent(out),dimension (0:26):: Iposf,Rposf,Lposf
    class(atom_config),intent(in)::atcf
     character(len=26),intent(in)::carac
    integer,intent(in),optional::nmask
    integer::nmaskV
    integer:: size1,size3,sizeV,size9
    integer::nvi,nvr,nvl
    !    integer, dimension (0:26):: Iposf,Rposf,Lposf
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   

    if (present(nmask))then
       nmaskV=nmask
    else
       nmaskV=atcf%imm
    end if

    size1=nmaskV;size3=3*size1; size9=3*size3
!    call atcf%print_type('in buffersize')
    nvi=0; sizeI=0
    Iposf(:)=0
    nvR=0; sizeR=0
    Rposf(:)=0
    nvl=0; sizel=0
    lposf(:)=0
    if(scan('n',carac).ne.0)  then
       nvi=nvi+1
       Iposf(nvi)=Iposf(nvi-1)+size1
!       call MPI_SEND(atcf%num_at_glob, size1, MPI_INTEGER, rgcib,104,comm,ierr)
    end if
    if(scan('i',carac).ne.0) then
       !call MPI_SEND(atcf%ityp, size1, MPI_INTEGER, rgcib,105,comm,ierr)
       nvi=nvi+1
       Iposf(nvi)=Iposf(nvi-1)+size1
 !      write(6,*)'nvi iposf ityp', nvi, iposf(nvi)
    end if
    if(scan('e',carac).ne.0) then
       !       call MPI_SEND(atcf%ielat, size1, MPI_INTEGER, rgcib,106,comm,ierr)
       nvi=nvi+1
       Iposf(nvi)=Iposf(nvi-1)+size1
!       write(6,*)'nvi iposf ielat', nvi, iposf(nvi)
    end if
    if(scan('p',carac).ne.0) then
       nvi=nvi+1
       Iposf(nvi)=Iposf(nvi-1)+size1
!       write(6,*)'nvi iposf proc_at', nvi, iposf(nvi)
       !       call MPI_SEND(atcf%proc_at, size1, MPI_INTEGER, rgcib,102,comm,ierr)
    end if
    if(scan('l',carac).ne.0)  then
          !call MPI_SEND(atcf%lgul, size1, MPI_LOGICAL, rgcib,103,comm,ierr)
       nvl=nvl+1
       lposf(nvl)=lposf(nvl-1)+size1
    end if
    if(scan('x',carac).ne.0)   then
       nvR=nvR+1
       Rposf(nvR)=Rposf(nvR-1)+size3
!       call MPI_SEND(atcf%xp, size3, NDM_MPI_REAL_DOUBLE, rgcib,100,comm,ierr)
    end if
    if(scan('f',carac).ne.0) then
       nvR=nvR+1
       Rposf(nvR)=Rposf(nvR-1)+size3
       !call MPI_SEND(atcf%fp, size3, NDM_MPI_REAL_DOUBLE, rgcib,101,comm,ierr)
    end if
    
    if (atcf%ltabvois) then
       sizeV=size(atcf%indi)
       if(scan('w',carac).ne.0) then
          nvi=nvi+1
        Iposf(nvi)=Iposf(nvi-1)+size1
!          call MPI_SEND(atcf%iwmax, size1, MPI_INTEGER, rgcib,107,comm,ierr)
       end if
       if(scan('d',carac).ne.0) then
          nvi=nvi+1
          Iposf(nvi)=Iposf(nvi-1)+sizeV
!          call MPI_SEND(atcf%indi, sizeV, MPI_INTEGER, rgcib,108,comm,ierr)
       end if
    end if

    select type (atcf)
    class is  (atom_config_d)
       if(scan('v',carac).ne.0) then
          nvR=nvR+1
       Rposf(nvR)=Rposf(nvR-1)+size3
          !call MPI_SEND(atcf%vp, size3, NDM_MPI_REAL_DOUBLE, rgcib,109,comm,ierr)
       end if
    end select
    select type (atcf)
    class is  (atom_config_e)
       if (atcf%lxpp)then
          if(scan('r',carac).ne.0) then
             nvR=nvR+1
             Rposf(nvR)=Rposf(nvR-1)+size3
             !call MPI_SEND(atcf%xpp, size3, NDM_MPI_REAL_DOUBLE, rgcib,110,comm,ierr)
          end if
       end if
       if (atcf%lprteat)then
          if(scan('u',carac).ne.0) then
             nvR=nvR+1
             Rposf(nvR)=Rposf(nvR-1)+size1
             !call MPI_SEND(atcf%eat, size1, NDM_MPI_REAL_DOUBLE, rgcib,111,comm,ierr)
          end if
       end if
       if (atcf%llangevin)then
          if(scan('g',carac).ne.0) then
             nvR=nvR+1
             Rposf(nvR)=Rposf(nvR-1)+size3
             !call MPI_SEND(atcf%glangv, size3, NDM_MPI_REAL_DOUBLE, rgcib,112,comm,ierr)
          end if
       end if
       if (atcf%lax)then
          if(scan('a',carac).ne.0) then
             nvR=nvR+1
             Rposf(nvR)=Rposf(nvR-1)+size3
             !call MPI_SEND(atcf%ax, size3, NDM_MPI_REAL_DOUBLE, rgcib,113,comm,ierr)
          end if
       end if
       if (atcf%lsigat)then
          if(scan('s',carac).ne.0) then
             nvR=nvR+1
             Rposf(nvR)=Rposf(nvR-1)+size9
             !call MPI_SEND(atcf%sigat, 3*size3, NDM_MPI_REAL_DOUBLE, rgcib,114,comm,ierr)
          end if
       end if
    end select
    select type (atcf)
    class is  (atom_config_arps)
       if(scan('m',carac).ne.0) then
          nvi=nvi+1
          Iposf(nvi)=iposf(nvi-1)+size1
          !call MPI_SEND(atcf%vp, size3, NDM_MPI_REAL_DOUBLE, rgcib,109,comm,ierr)
       end if
    end select

!****************************************************
    sizeI=Iposf(nvI);    sizeR=Rposf(nvR);    sizel=Lposf(nvl)


    return
  end subroutine buffersizes

  subroutine buildbuff(atcf,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf,maskR,nmask)
     integer, intent(in)::sizeI,sizel,sizer
    integer, intent(in),dimension (0:26):: Iposf,Rposf,Lposf
    class(atom_config),intent(in)::atcf
     character(len=26),intent(in)::carac
    integer,allocatable,intent(inout):: ibuffer(:)
    logical,allocatable,intent(inout)::lbuffer(:)
    integer,intent(out)::csi,csr,csl
    real(double),allocatable,intent(inout)::rbuffer(:)
    logical, optional,intent(in)::maskR(:)

   logical,allocatable :: mask(:)
    integer:: size1,size3,sizeV,size9,iat
    integer::ivi,ivr,ivl,ibi,ibr,ibl,ic,ic2,ib,ip
    integer,intent(in),optional::nmask
    integer::nmaskV

    allocate(mask(atcf%imm))
    
    if (present(nmask))then
       nmaskV=nmask
    else
       nmaskV=atcf%imm
    end if
    if (present(maskR))then
       mask(:)=maskR(:)
    else
       mask=.true.
    end if
    

    size1=nmaskV;size3=3*size1; size9=3*size3
   


    ibi=0;ibl=0;ibr=0
    ivi=0;ivl=0;ivR=0
    csi=0;csl=0;csr=0
    !****************************************************
!    write(6,*)'CARAC',carac
    if(scan('n',carac).ne.0)  then
       ivi=ivi+1
       ip=0
       do iat=1,atcf%imm
          if (mask(iat).eqv..true.) then
             ip=ip+1
             ib=Iposf(ivi-1)+ip
             ibuffer(ib)=atcf%num_at_glob(iat)
             csi=csi+1
          end if
       end do
!       call MPI_SEND(atcf%num_at_glob, size1, MPI_INTEGER, rgcib,104,comm,ierr)
    end if
    if(scan('i',carac).ne.0) then
       !call MPI_SEND(atcf%ityp, size1, MPI_INTEGER, rgcib,105,comm,ierr)
       ivi=ivi+1
       ip=0
       do iat=1,atcf%imm
          if (mask(iat).eqv..true.) then
             ip=ip+1
             ib=Iposf(ivi-1)+ip
             ibuffer(ib)=atcf%ityp(iat)
             csi=csi+1
          end if
       end do
    end if
    if(scan('e',carac).ne.0) then
       !       call MPI_SEND(atcf%ielat, size1, MPI_INTEGER, rgcib,106,comm,ierr)
       ivi=ivi+1
       ip=0
       do iat=1,atcf%imm
          if (mask(iat).eqv..true.) then
             ip=ip+1
             ib=Iposf(ivi-1)+ip
             ibuffer(ib)=atcf%ielat(iat)
             csi=csi+1
          end if
       end do
    end if
#ifdef PARA
    if(scan('p',carac).ne.0) then
       ivi=ivi+1
       ip=0
       do iat=1,atcf%imm
          if (mask(iat).eqv..true.) then
             ip=ip+1
             ib=Iposf(ivi-1)+ip
             ibuffer(ib)=atcf%proc_at(iat)
             csi=csi+1
          end if
       end do
       !       call MPI_SEND(atcf%proc_at, size1, MPI_INTEGER, rgcib,102,comm,ierr)
    end if
#endif
    if(scan('l',carac).ne.0)  then
       !call MPI_SEND(atcf%lgul, size1, MPI_LOGICAL, rgcib,103,comm,ierr)
       ivl=ivl+1
       ip=0
       do iat=1,atcf%imm
          if (mask(iat).eqv..true.) then
             ip=ip+1
             ib=Lposf(ivl-1)+ip
             Lbuffer(ib)=atcf%lgul(iat)
             csl=csl+1
          end if
       end do
    end if
    if(scan('x',carac).ne.0)   then
       ivR=ivR+1
       ip=0
       do iat=1,atcf%imm
          if (mask(iat).eqv..true.) then
             do ic=1,3
                ip=ip+1
                ib=Rposf(ivR-1)+ip
                rbuffer(ib)=atcf%xp(ic,iat)
                csR=csR+1
             end do
          end if
       end do
!       call MPI_SEND(atcf%xp, size3, NDM_MPI_REAL_DOUBLE, rgcib,100,comm,ierr)
    end if
    if(scan('f',carac).ne.0) then
       ivR=ivR+1
       ip=0
       do iat=1,atcf%imm
          if (mask(iat).eqv..true.) then
             do ic=1,3
                ip=ip+1
                ib=Rposf(ivR-1)+ip
                rbuffer(ib)=atcf%fp(ic,iat)
                csR=csR+1
             end do
          end if
       end do
       !call MPI_SEND(atcf%fp, size3, NDM_MPI_REAL_DOUBLE, rgcib,101,comm,ierr)
    end if
    
    if (atcf%ltabvois) then

       if(scan('w',carac).ne.0) then
          if (any(mask(1:atcf%im).eqv..false.))then
             write(6,*)'trf latbavois et conf incomplète stop'
             call arret_ndm
          end if
          ivi=ivi+1
          do ip=1,size1
             ib=Iposf(ivi-1)+ip
             ibuffer(ib)=atcf%iwmax(ip)
             csi=csi+1
          end do
      
!          call MPI_SEND(atcf%iwmax, size1, MPI_INTEGER, rgcib,107,comm,ierr)
       end if
       if(scan('d',carac).ne.0) then
          sizeV=size(atcf%indi)
          ivi=ivi+1
          do ip=1,sizeV
             ib=Iposf(ivi-1)+ip
             ibuffer(ib)=atcf%indi(ip)
             csi=csi+1
          end do
!          call MPI_SEND(atcf%indi, sizeV, MPI_INTEGER, rgcib,108,comm,ierr)
       end if
    end if

    select type (atcf)
    class is  (atom_config_d)
       if(scan('v',carac).ne.0) then
          ivR=ivR+1
          ip=0
          do iat=1,atcf%imm
             if (mask(iat).eqv..true.) then
                do ic=1,3
                   ip=ip+1
                   ib=Rposf(ivR-1)+ip
                   rbuffer(ib)=atcf%vp(ic,iat)
                   csR=csR+1
                end do
             end if
          end do
          !call MPI_SEND(atcf%vp, size3, NDM_MPI_REAL_DOUBLE, rgcib,109,comm,ierr)
       end if

    end select
    select type (atcf)
    class is  (atom_config_e)
       if (atcf%lxpp)then
          if(scan('r',carac).ne.0) then
             ivR=ivR+1
             ip=0
             do iat=1,atcf%imm
                if (mask(iat).eqv..true.) then
                   do ic=1,3
                      ip=ip+1
                      ib=Rposf(ivR-1)+ip
                      rbuffer(ib)=atcf%xpp(ic,iat)
                      csR=csR+1
                   end do
                end if
             end do
             !call MPI_SEND(atcf%xpp, size3, NDM_MPI_REAL_DOUBLE, rgcib,110,comm,ierr)
          end if
       end if

       if (atcf%lprteat)then
          if(scan('u',carac).ne.0) then
          ivR=ivR+1
          ip=0
          do iat=1,atcf%imm
             if (mask(iat).eqv..true.) then
                ip=ip+1
                ib=Rposf(ivR-1)+ip
                rbuffer(ib)=atcf%eat(iat)
                csR=csR+1
             end if
          end do
             !call MPI_SEND(atcf%eat, size1, NDM_MPI_REAL_DOUBLE, rgcib,111,comm,ierr)
          end if
       end if
       if (atcf%llangevin)then
          if(scan('g',carac).ne.0) then
             ivR=ivR+1
             ip=0
             do iat=1,atcf%imm
                if (mask(iat).eqv..true.) then
                   do ic=1,3
                      ip=ip+1
                      ib=Rposf(ivR-1)+ip
                      rbuffer(ib)=atcf%glangv(ic,iat)
                      csR=csR+1
                   end do
                end if
             end do
             !call MPI_SEND(atcf%glangv, size3, NDM_MPI_REAL_DOUBLE, rgcib,112,comm,ierr)
          end if
       end if
       if (atcf%lax)then
          if(scan('a',carac).ne.0) then
             !call MPI_SEND(atcf%ax, size3, NDM_MPI_REAL_DOUBLE, rgcib,113,comm,ierr)
             ivR=ivR+1
             ip=0
             do iat=1,atcf%imm
                if (mask(iat).eqv..true.) then
                   do ic=1,3
                      ip=ip+1
                      ib=Rposf(ivR-1)+ip
                      rbuffer(ib)=atcf%ax(ic,iat)
                      csR=csR+1
                   end do
                end if
             end do
          end if
       end if
       if (atcf%lsigat)then
          if(scan('s',carac).ne.0) then
             ivR=ivR+1
             ip=0
             do iat=1,atcf%imm
                if (mask(iat).eqv..true.) then
                   do ic=1,3
                      do ic2=1,3
                         ip=ip+1
                         ib=Rposf(ivR-1)+ip
                         rbuffer(ib)=atcf%sigat(ic,ic2,iat)
                         csR=csR+1
                      end do
                   end do
                end if
             end do
             !call MPI_SEND(atcf%sigat, 3*size3, NDM_MPI_REAL_DOUBLE, rgcib,114,comm,ierr)
          end if
       end if
       
    end select
    select type (atcf)
    class is  (atom_config_arps)
       if(scan('m',carac).ne.0) then
          ivi=ivi+1
          ip=0
          do iat=1,atcf%imm
             if (mask(iat).eqv..true.) then
                ip=ip+1
                ib=Iposf(ivi-1)+ip
                ibuffer(ib)=atcf%mov(iat)
                csi=csi+1
             end if
          end do
       end if
    end select

    if (csi.ne.sizeI) then
       write(6,*)'erreur CSI 2',sizeI,csi
!       call endmpi
       call arret_ndm
    end if
    if (csr.ne.sizer) then
       write(6,*)'erreur CSR ',sizeR,csr
!       call endmpi
       call arret_ndm
    end if
    if (csl.ne.sizel) then
       write(6,*)'erreur CSL ',sizel,csl
!       call endmpi
       call arret_ndm
    end if

    return
  end subroutine buildbuff

  
  subroutine addatim(atcf,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf,immax)
     integer, intent(in)::sizeI,sizel,sizer,immax
    integer, intent(in),dimension (0:26):: Iposf,Rposf,Lposf
    class(atom_config),intent(inout)::atcf
     character(len=26),intent(in)::carac
    integer,allocatable,intent(in):: ibuffer(:)
    logical,allocatable,intent(in)::lbuffer(:)
    integer,intent(out)::csi,csr,csl
    real(double),allocatable,intent(in)::rbuffer(:)

    integer:: size1,size3,sizeV,size9
    integer::ivi,ivr,ivl,ibi,ibr,ibl,ic,ic2,ib,imp


    size1=1;size3=3*size1; size9=3*size3
    ibi=0;ibl=0;ibr=0
    ivi=0;ivl=0;ivR=0
    csi=0;csl=0;csr=0

    atcf%im=atcf%im+1
    imp=atcf%im
    if(scan('n',carac).ne.0)  then
       ivi=ivi+1
       ib=Iposf(ivi-1)+1
       atcf%num_at_glob(imp)=ibuffer(ib)
       csi=csi+1
    end if
    if(scan('i',carac).ne.0) then
       ivi=ivi+1
       ib=Iposf(ivi-1)+1
       atcf%ityp(imp)=ibuffer(ib)
       csi=csi+1
    end if
    if(scan('e',carac).ne.0) then
       ivi=ivi+1
       ib=Iposf(ivi-1)+1
       atcf%ielat(imp)=ibuffer(ib)
       csi=csi+1
    end if

#ifdef PARA
    if(scan('p',carac).ne.0) then
       ivi=ivi+1
          ib=Iposf(ivi-1)+1
          atcf%proc_at(imp)=ibuffer(ib)
          csi=csi+1
    end if
#endif
    
    if(scan('l',carac).ne.0)  then
       ivl=ivl+1
          ib=Lposf(ivi-1)+1
          atcf%lgul(imp)=Lbuffer(ib)
          csl=csl+1
    end if
    if(scan('x',carac).ne.0)   then
       ivR=ivR+1
          do ic=1,3
             ib=Rposf(ivR-1)+ic
             atcf%xp(ic,imp)=rbuffer(ib)
             csR=csR+1
       end do
    end if
    if(scan('f',carac).ne.0) then
       ivR=ivR+1
       do ic=1,3
          ib=Rposf(ivR-1)+ic
          atcf%fp(ic,imp)=rbuffer(ib)
          csR=csR+1
       end do
    end if

    if (atcf%ltabvois) then
       sizeV=size(atcf%indi)
       if(scan('w',carac).ne.0) then
          ivi=ivi+1
          ib=Iposf(ivi-1)+1
          atcf%iwmax(imp)=ibuffer(ib)
          csi=csi+1
       end if
       if(scan('d',carac).ne.0) then
          ivi=ivi+1
             ib=Iposf(ivi-1)+1
             atcf%indi(imp)=ibuffer(ib)
             csi=csi+1
       end if
    end if

    select type (atcf)
       class is  (atom_config_d)
       if(scan('v',carac).ne.0) then
          ivR=ivR+1
          do ic=1,3
             ib=Rposf(ivR-1)+ic
             atcf%vp(ic,imp)=rbuffer(ib)
             csR=csR+1
          end do
       end if

    end select
    select type (atcf)
    class is  (atom_config_e)
       if (atcf%lxpp)then
          if(scan('r',carac).ne.0) then
             ivR=ivR+1
             do ic=1,3
                ib=Rposf(ivR-1)+ic
                atcf%xpp(ic,imp)= rbuffer(ib)
                csR=csR+1
                end do
          end if
       end if

       if (atcf%lprteat)then
          if(scan('u',carac).ne.0) then
             ivR=ivR+1
             ib=Rposf(ivR-1)+1
             atcf%eat(imp)=Rbuffer(ib)
             csR=csR+1
          end if
       end if
       if (atcf%llangevin)then
          if(scan('g',carac).ne.0) then
             ivR=ivR+1
             do ic=1,3
                ib=Rposf(ivR-1)+ic
                atcf%glangv(ic,imp)=rbuffer(ib)
                csR=csR+1
             end do
          end if
       end if
       if (atcf%lax)then
          if(scan('a',carac).ne.0) then
             ivR=ivR+1
             do ic=1,3
                ib=Rposf(ivR-1)+ic
                   atcf%ax(ic,imp)= rbuffer(ib)
                   csR=csR+1
                end do
          end if
       end if
       if (atcf%lsigat)then
          if(scan('s',carac).ne.0) then
             ivR=ivR+1
             do ic=1,3
                do ic2=1,3
                   ib=Rposf(ivR-1)+(ic-1)*3+ic2
                   atcf%sigat(ic,ic2,imp)=rbuffer(ib)
                   csR=csR+1
                end do
             end do
          end if
       end if
    end select
    select type (atcf)
       class is  (atom_config_arps)

        if(scan('m',carac).ne.0) then
           !call MPI_SEND(atcf%ityp, size1, MPI_INTEGER, rgcib,105,comm,ierr)
           ivi=ivi+1
           ib=Iposf(ivi-1)+1
           atcf%mov(imp)=ibuffer(ib)
           csi=csi+1
        end if
     end select

    if (csi.ne.sizeI) then
       write(6,*)'erreur CSI 1 ',sizeI,csi
!       call endmpi
       call arret_ndm
    end if
    if (csr.ne.sizer) then
       write(6,*)'erreur CSR ',sizeR,csr
!       call endmpi
       call arret_ndm
    end if
    if (csl.ne.sizel) then
       write(6,*)'erreur CSL ',sizel,csl
!       call endmpi
       call arret_ndm
    end if
  end subroutine addatim
  subroutine copybuff(atcf,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf,immax)
     integer, intent(in)::sizeI,sizel,sizer,immax
    integer, intent(in),dimension (0:26):: Iposf,Rposf,Lposf
    class(atom_config),intent(inout)::atcf
     character(len=26),intent(in)::carac
    integer,allocatable,intent(in):: ibuffer(:)
    logical,allocatable,intent(in)::lbuffer(:)
    integer,intent(out)::csi,csr,csl
    real(double),allocatable,intent(in)::rbuffer(:)

    integer:: size1,size3,sizeV,size9
    integer::ivi,ivr,ivl,ibi,ibr,ibl,ic,ic2,ib,ip


    size1=immax;size3=3*size1; size9=3*size3
    ibi=0;ibl=0;ibr=0
    ivi=0;ivl=0;ivR=0
    csi=0;csl=0;csr=0


    if(scan('n',carac).ne.0)  then
       ivi=ivi+1
       do ip=1,size1
          ib=Iposf(ivi-1)+ip
          atcf%num_at_glob(ip)=ibuffer(ib)
          csi=csi+1
       end do
       !       call MPI_SEND(atcf%num_at_glob, size1, MPI_INTEGER, rgcib,104,comm,ierr)
    end if
    if(scan('i',carac).ne.0) then
       !call MPI_SEND(atcf%ityp, size1, MPI_INTEGER, rgcib,105,comm,ierr)
       ivi=ivi+1
       do ip=1,size1
          ib=Iposf(ivi-1)+ip
          atcf%ityp(ip)=ibuffer(ib)
          csi=csi+1
       end do
    end if
    if(scan('e',carac).ne.0) then
       !       call MPI_SEND(atcf%ielat, size1, MPI_INTEGER, rgcib,106,comm,ierr)
       ivi=ivi+1
       do ip=1,size1
          ib=Iposf(ivi-1)+ip
          atcf%ielat(ip)=ibuffer(ib)
          csi=csi+1
       end do
    end if
#ifdef PARA
    if(scan('p',carac).ne.0) then
       ivi=ivi+1
       do ip=1,size1
          ib=Iposf(ivi-1)+ip
          atcf%proc_at(ip)=ibuffer(ib)
          csi=csi+1
       end do
       !       call MPI_SEND(atcf%proc_at, size1, MPI_INTEGER, rgcib,102,comm,ierr)
    end if
#endif
    
    if(scan('l',carac).ne.0)  then
       !call MPI_SEND(atcf%lgul, size1, MPI_LOGICAL, rgcib,103,comm,ierr)
       ivl=ivl+1
       do ip=1,size1
          ib=Lposf(ivi-1)+ip
          atcf%lgul(ip)=Lbuffer(ib)
          csl=csl+1
       end do
    end if
    if(scan('x',carac).ne.0)   then
       ivR=ivR+1
       do ip=1,size1
          do ic=1,3
             ib=Rposf(ivR-1)+3*(ip-1)+ic
             atcf%xp(ic,ip)=rbuffer(ib)
             csR=csR+1
          end do
       end do
       !       call MPI_SEND(atcf%xp, size3, NDM_MPI_REAL_DOUBLE, rgcib,100,comm,ierr)
    end if
    if(scan('f',carac).ne.0) then
       ivR=ivR+1
       do ip=1,size1
          do ic=1,3
             ib=Rposf(ivR-1)+3*(ip-1)+ic
             atcf%fp(ic,ip)=rbuffer(ib)
             csR=csR+1
          end do
       end do
       !call MPI_SEND(atcf%fp, size3, NDM_MPI_REAL_DOUBLE, rgcib,101,comm,ierr)
    end if

    if (atcf%ltabvois) then
       sizeV=size(atcf%indi)
       if(scan('w',carac).ne.0) then
          ivi=ivi+1
          do ip=1,size1
             ib=Iposf(ivi-1)+ip
             atcf%iwmax(ip)=ibuffer(ib)
             csi=csi+1
          end do

          !          call MPI_SEND(atcf%iwmax, size1, MPI_INTEGER, rgcib,107,comm,ierr)
       end if
       if(scan('d',carac).ne.0) then
          ivi=ivi+1
          do ip=1,sizeV
             ib=Iposf(ivi-1)+ip
             atcf%indi(ip)=ibuffer(ib)
             csi=csi+1
          end do
          !          call MPI_SEND(atcf%indi, sizeV, MPI_INTEGER, rgcib,108,comm,ierr)
       end if
    end if

    select type (atcf)
       class is  (atom_config_d)
       if(scan('v',carac).ne.0) then
          ivR=ivR+1
          do ip=1,size1
             do ic=1,3
                ib=Rposf(ivR-1)+3*(ip-1)+ic
                atcf%vp(ic,ip)=rbuffer(ib)
                csR=csR+1
             end do
          end do
          !call MPI_SEND(atcf%vp, size3, NDM_MPI_REAL_DOUBLE, rgcib,109,comm,ierr)
       end if

    end select
    select type (atcf)
    class is  (atom_config_e)
       if (atcf%lxpp)then
          if(scan('r',carac).ne.0) then
             ivR=ivR+1
             do ip=1,size1
                do ic=1,3
                   ib=Rposf(ivR-1)+3*(ip-1)+ic
                   atcf%xpp(ic,ip)= rbuffer(ib)
                   csR=csR+1
                end do
             end do
             !call MPI_SEND(atcf%xpp, size3, NDM_MPI_REAL_DOUBLE, rgcib,110,comm,ierr)
          end if
       end if

       if (atcf%lprteat)then
          if(scan('u',carac).ne.0) then
             ivR=ivR+1
             do ip=1,size1
                ib=Rposf(ivR-1)+ip
                atcf%eat(ip)=Rbuffer(ib)
                csR=csR+1
             end do
             !call MPI_SEND(atcf%eat, size1, NDM_MPI_REAL_DOUBLE, rgcib,111,comm,ierr)
          end if
       end if
       if (atcf%llangevin)then
          if(scan('g',carac).ne.0) then
             ivR=ivR+1
             do ip=1,size1
                do ic=1,3
                   ib=Rposf(ivR-1)+3*(ip-1)+ic
                   atcf%glangv(ic,ip)=rbuffer(ib)
                   csR=csR+1
                end do
             end do
             !call MPI_SEND(atcf%glangv, size3, NDM_MPI_REAL_DOUBLE, rgcib,112,comm,ierr)
          end if
       end if
       if (atcf%lax)then
          if(scan('a',carac).ne.0) then
             ivR=ivR+1
             do ip=1,size1
                do ic=1,3
                   ib=Rposf(ivR-1)+3*(ip-1)+ic
                   atcf%ax(ic,ip)= rbuffer(ib)
                   csR=csR+1
                end do
             end do
             !call MPI_SEND(atcf%ax, size3, NDM_MPI_REAL_DOUBLE, rgcib,113,comm,ierr)
          end if
       end if
       if (atcf%lsigat)then
          if(scan('s',carac).ne.0) then
             ivR=ivR+1
             do ip=1,size1
                do ic=1,3
                   do ic2=1,3
                      ib=Rposf(ivR-1)+(ip-1)*9+(ic-1)*3+ic2
                      atcf%sigat(ic,ic2,ip)=rbuffer(ib)
                      csR=csR+1
                   end do
                end do
             end do
             !call MPI_SEND(atcf%sigat, 3*size3, NDM_MPI_REAL_DOUBLE, rgcib,114,comm,ierr)
          end if
       end if
    end select
    select type (atcf)
       class is  (atom_config_arps)

        if(scan('m',carac).ne.0) then
           !call MPI_SEND(atcf%ityp, size1, MPI_INTEGER, rgcib,105,comm,ierr)
           ivi=ivi+1
           do ip=1,size1
              ib=Iposf(ivi-1)+ip
              atcf%mov(ip)=ibuffer(ib)
              csi=csi+1
           end do
        end if
     end select

    if (csi.ne.sizeI) then
       write(6,*)'erreur CSI 1 ',sizeI,csi
!       call endmpi
       call arret_ndm
    end if
    if (csr.ne.sizer) then
       write(6,*)'erreur CSR ',sizeR,csr
!       call endmpi
       call arret_ndm
    end if
    if (csl.ne.sizel) then
       write(6,*)'erreur CSL ',sizel,csl
!       call endmpi
       call arret_ndm
    end if

  end subroutine copybuff


  subroutine distribnag(nag,inag,immax,atcf,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf)
    integer,intent(in),allocatable::nag(:),inag(:)  ! nag (1:imloc) liste des num_at_glob ; inag(1:max(num_at_glob))pointeur inverse de num_at_glob vers icomp
    integer, intent(in)::sizeI,sizel,sizer,immax
    integer, intent(in),dimension (0:26):: Iposf,Rposf,Lposf
    class(atom_config),intent(inout)::atcf
     character(len=26),intent(in)::carac
    integer,allocatable,intent(in):: ibuffer(:)
    logical,allocatable,intent(in)::lbuffer(:)
    integer,intent(out)::csi,csr,csl
    real(double),allocatable,intent(in)::rbuffer(:)

    integer:: size1,size3,sizeV,size9
    integer::ivi,ivr,ivl,ibi,ibr,ibl,ic,ic2,ib,ip,icomp

    size1=immax;size3=3*size1; size9=3*size3
    ibi=0;ibl=0;ibr=0
    ivi=0;ivl=0;ivR=0
    csi=0;csl=0;csr=0

    if(scan('n',carac).ne.0)  then
       ivi=ivi+1
       do ip=1,size1
          ib=Iposf(ivi-1)+ip
          icomp=inag(nag(ip))
          atcf%num_at_glob(icomp)=ibuffer(ib)
          csi=csi+1
       end do
       !       call MPI_SEND(atcf%num_at_glob, size1, MPI_INTEGER, rgcib,104,comm,ierr)
    end if
    if(scan('i',carac).ne.0) then
       !call MPI_SEND(atcf%ityp, size1, MPI_INTEGER, rgcib,105,comm,ierr)
       ivi=ivi+1
       do ip=1,size1
          ib=Iposf(ivi-1)+ip
          atcf%ityp(inag(nag(ip)))=ibuffer(ib)
          csi=csi+1
       end do
    end if
    if(scan('e',carac).ne.0) then
       !       call MPI_SEND(atcf%ielat, size1, MPI_INTEGER, rgcib,106,comm,ierr)
       ivi=ivi+1
       do ip=1,size1
          ib=Iposf(ivi-1)+ip
          atcf%ielat(inag(nag(ip)))=ibuffer(ib)
          csi=csi+1
       end do
    end if
#ifdef PARA
    if(scan('p',carac).ne.0) then
       ivi=ivi+1
       do ip=1,size1
          ib=Iposf(ivi-1)+ip
          atcf%proc_at(inag(nag(ip)))=ibuffer(ib)
          csi=csi+1
       end do
       !       call MPI_SEND(atcf%proc_at, size1, MPI_INTEGER, rgcib,102,comm,ierr)
    end if
#endif
    
    if(scan('l',carac).ne.0)  then
       !call MPI_SEND(atcf%lgul, size1, MPI_LOGICAL, rgcib,103,comm,ierr)
       ivl=ivl+1
       do ip=1,size1
          ib=Lposf(ivl-1)+ip
          atcf%lgul(inag(nag(ip)))=Lbuffer(ib)
          csl=csl+1
       end do
    end if
    if(scan('x',carac).ne.0)   then
       ivR=ivR+1
       do ip=1,size1
          do ic=1,3
             ib=Rposf(ivR-1)+3*(ip-1)+ic
             atcf%xp(ic,inag(nag(ip)))=rbuffer(ib)
             csR=csR+1
          end do
       end do
       !       call MPI_SEND(atcf%xp, size3, NDM_MPI_REAL_DOUBLE, rgcib,100,comm,ierr)
    end if
    if(scan('f',carac).ne.0) then
       ivR=ivR+1
       do ip=1,size1
          do ic=1,3
             ib=Rposf(ivR-1)+3*(ip-1)+ic
             atcf%fp(ic,inag(nag(ip)))=rbuffer(ib)
             csR=csR+1
          end do
       end do
       !call MPI_SEND(atcf%fp, size3, NDM_MPI_REAL_DOUBLE, rgcib,101,comm,ierr)
    end if

    if (atcf%ltabvois) then
       sizeV=size(atcf%indi)
       if(scan('w',carac).ne.0) then
          ivi=ivi+1
          do ip=1,size1
             ib=Iposf(ivi-1)+ip
             atcf%iwmax(inag(nag(ip)))=ibuffer(ib)
             csi=csi+1
          end do

          !          call MPI_SEND(atcf%iwmax, size1, MPI_INTEGER, rgcib,107,comm,ierr)
       end if
       if(scan('d',carac).ne.0) then
          ivi=ivi+1
          do ip=1,sizeV
             ib=Iposf(ivi-1)+ip
             atcf%indi(inag(nag(ip)))=ibuffer(ib)
             csi=csi+1
          end do
          !          call MPI_SEND(atcf%indi, sizeV, MPI_INTEGER, rgcib,108,comm,ierr)
       end if
    end if

    select type (atcf)
       class is  (atom_config_d)
       if(scan('v',carac).ne.0) then
          ivR=ivR+1
          do ip=1,size1
             do ic=1,3
                ib=Rposf(ivR-1)+3*(ip-1)+ic
                atcf%vp(ic,inag(nag(ip)))=rbuffer(ib)
                csR=csR+1
             end do
          end do
          !call MPI_SEND(atcf%vp, size3, NDM_MPI_REAL_DOUBLE, rgcib,109,comm,ierr)
       end if

    end select
    select type (atcf)
    class is  (atom_config_e)
       if (atcf%lxpp)then
          if(scan('r',carac).ne.0) then
             ivR=ivR+1
             do ip=1,size1
                do ic=1,3
                   ib=Rposf(ivR-1)+3*(ip-1)+ic
                   atcf%xpp(ic,inag(nag(ip)))= rbuffer(ib)
                   csR=csR+1
                end do
             end do
             !call MPI_SEND(atcf%xpp, size3, NDM_MPI_REAL_DOUBLE, rgcib,110,comm,ierr)
          end if
       end if
       if (atcf%lprteat)then
          if(scan('u',carac).ne.0) then
             ivR=ivR+1
             do ip=1,size1
                ib=Rposf(ivR-1)+ip
                atcf%eat(inag(nag(ip)))=Rbuffer(ib)
                csR=csR+1
             end do
             !call MPI_SEND(atcf%eat, size1, NDM_MPI_REAL_DOUBLE, rgcib,111,comm,ierr)
          end if
       end if
       if (atcf%llangevin)then
          if(scan('g',carac).ne.0) then
             ivR=ivR+1
             do ip=1,size1
                do ic=1,3
                   ib=Rposf(ivR-1)+3*(ip-1)+ic
                   atcf%glangv(ic,inag(nag(ip)))=rbuffer(ib)
                   csR=csR+1
                end do
             end do
             !call MPI_SEND(atcf%glangv, size3, NDM_MPI_REAL_DOUBLE, rgcib,112,comm,ierr)
          end if
       end if
       if (atcf%lax)then
          if(scan('a',carac).ne.0) then
             ivR=ivR+1
             do ip=1,size1
                do ic=1,3
                   ib=Rposf(ivR-1)+3*(ip-1)+ic
                   atcf%ax(ic,inag(nag(ip)))= rbuffer(ib)
                   csR=csR+1
                end do
             end do
             !call MPI_SEND(atcf%ax, size3, NDM_MPI_REAL_DOUBLE, rgcib,113,comm,ierr)
          end if
       end if
       if (atcf%lsigat)then
          if(scan('s',carac).ne.0) then
             ivR=ivR+1
             do ip=1,size1
                do ic=1,3
                   do ic2=1,3
                      ib=Rposf(ivR-1)+(ip-1)*9+(ic-1)*3+ic2
                      atcf%sigat(ic,ic2,inag(nag(ip)))=rbuffer(ib)
                      csR=csR+1
                   end do
                end do
             end do
             !call MPI_SEND(atcf%sigat, 3*size3, NDM_MPI_REAL_DOUBLE, rgcib,114,comm,ierr)
          end if
       end if
    end select
    select type (atcf)
    class is  (atom_config_arps)

       if(scan('m',carac).ne.0) then
          !call MPI_SEND(atcf%ityp, size1, MPI_INTEGER, rgcib,105,comm,ierr)
          ivi=ivi+1
          do ip=1,size1
             ib=Iposf(ivi-1)+ip
             atcf%mov(inag(nag(ip)))=ibuffer(ib)
             csi=csi+1
          end do
       end if
    end select

    if (csi.ne.sizeI) then
       write(6,*)'erreur CSIC 1 ',sizeI,csi
!       call endmpi
       call arret_ndm
    end if
    if (csr.ne.sizer) then
       write(6,*)'erreur CSRC ',sizeR,csr
!       call endmpi
       call arret_ndm
    end if
    if (csl.ne.sizel) then
       write(6,*)'erreur CSLC ',sizel,csl
!       call endmpi
       call arret_ndm
    end if

  end subroutine distribnag
  
  subroutine lgcheck(atconf,message)
    class(atom_config)::atconf
    character(*)::message
    if (atconf%lglock) then
       write(6,*)' lgul is already used correct code',rang
       write(6,*) message
       call arret_ndm
    end if
  end subroutine lgcheck
end module atomconfig


