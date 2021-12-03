module atomconfig
  USE T_kind_param_m,only:double,long
  USE Mat_utils_mod,only: fillbuffer3D,fillbuffer1D,fillbuffer9D
#ifdef PARA
  use mpi
  USE Tpara,only:NDM_MPI_REAL_DOUBLE,ierr,mpi_communicator,endmpi
  use gen_com_m,only:rang


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
     logical, allocatable:: lgul(:)
     integer,allocatable::num_at_glob(:)
#ifdef PARA
     integer,allocatable::proc_at(:) ! tableau de taille im_glob total indiquant le numéro du proc qui gère l'atome
#endif     
   contains
     procedure, pass::init=>init_atom_config
     procedure, pass::Eegal
     procedure, pass::copy_atom=>copy_atom_b
     procedure, pass::dealloc=>dealloc_atom_config
     procedure, pass::vers_master=>vers_master_atom !(atcfloc,atcfcomp,div)
     procedure, pass::master2loc=>master2loc_atom !(atcfcomp,atcfloc,div)
     procedure, pass::copy_config
     procedure, pass::print
     procedure, pass::pack
     procedure, pass::fab
     procedure, pass::backto
     procedure, pass::add2conf
     procedure, pass::extend
     procedure, pass::deftype
     procedure, pass::send2proc=>s2p_atom
     procedure, pass::send2all=>s2a_atom
     procedure, pass::recv=>rcv_atom
     procedure, pass::zero=>zero_atom
     procedure, pass::switch_atom

     !
  end type atom_config

  type, extends (atom_config):: atom_config_d ! type dynamique des configurations atomiques(+vp/+xpp). vp et xpp seront toujours allouées
     real(double),allocatable ::xpp(:,:)
     real(double),allocatable::vp(:,:)
   contains
     procedure, pass::copy_atom=>copy_atom_d
     procedure, pass::dealloc=>dealloc_atom_config_d
     !#ifdef PARA
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
     procedure, pass::zero=>zero_atom_e
   
     !#endif     
  end type atom_config_e

  private ::buffersizes
contains
  !initialisations
  
  subroutine init_atom_config(atconf,imin,immin,ltabvois,nvois,rvois,lreallocate,im_glob,imm_glob)
    class(atom_config),intent(inout)::atconf
    integer,intent(in):: imin
    logical,optional, intent(in)::ltabvois,lreallocate
    integer, optional::nvois,immin,im_glob,imm_glob
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
    type is (atom_config_d)
       if ((lrealloc).and.(allocated(atconf%vp)))then
          deallocate(atconf%vp); deallocate(atconf%xpp)
       end if
       if (.not.allocated(atconf%vp))then
          allocate(atconf%vp(3,atconf%imm));allocate(atconf%xpp(3,atconf%imm))
       end if
       atconf%vp=0;atconf%xpp=0
    type is (atom_config_e)
       if ((lrealloc).and.(allocated(atconf%vp)))then
          deallocate(atconf%vp); deallocate(atconf%xpp)
       end if
       if (.not.allocated(atconf%vp))then
          allocate(atconf%vp(3,atconf%imm));allocate(atconf%xpp(3,atconf%imm))
       end if
       atconf%vp=0;atconf%xpp=0
       
       !       if ((lrealloc).and.(allocated(atconf%vp)))then
!          deallocate(atconf%vp); deallocate(atconf%xpp)
!       end if
!       if (.not.allocated(atconf%vp))then
!          allocate(atconf%vp(3,atconf%imm));allocate(atconf%xpp(3,atconf%imm))
!       end if
!       atconf%vp=0;atconf%xpp=0
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
!!$    select type(atcf)
!!$    type is (atom_config)
!!$       call attemp_b%init(imcn,immcn,atcf%ltabvois,nvois,rvois)
!!$       attemp=>attemp_b
!!$    type is (atom_config_d)
!!$       call attemp_d%init(imcn,immcn,atcf%ltabvois,nvois,rvois)
!!$       attemp=>attemp_d
!!$    type is (atom_config_e)
!!$       call atcf%Eegal(attemp_e)
!!$       call attemp_e%init(imcn,immcn,atcf%ltabvois,nvois,rvois)
!!$       attemp=>attemp_e
!!$    end select
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

  subroutine s2p_atom (atcf, rgcib,mpic,caracT)
    class(atom_config):: atcf
    type(mpi_communicator),intent(in)::mpic
    integer,intent(in)::rgcib
    integer:: size1,size3,sizeV,size9
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    integer::nvi,nvr,sizeI,sizeR,nvl,sizel,ivi,ivr,ivl,ibi,ibr,ibl,ic,ic2,csi,csr,csl,ib,ip
    integer, dimension (0:26):: Iposf,Rposf,Lposf
    integer,allocatable:: ibuffer(:)
    logical,allocatable::lbuffer(:)
    real(double),allocatable::rbuffer(:)
    integer:: cst(3)
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   
    size1=atcf%imm;size3=3*size1; size9=3*size3
#ifdef PARA
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugas'
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
    integer:: size1,size3,sizeV,size9
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    integer::nvi,nvr,sizeI,sizeR,nvl,sizel,ivi,ivr,ivl,ibi,ibr,ibl,ic,ic2,csi,csr,csl,ib,ip
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
       carac='xfniewdlpvrugas'
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
       stop
    end if
    if (cst(2).ne.sizel) then
       write(6,*)'erreur CST2 ',sizel,cst(2)
       stop
    end if
    if (cst(3).ne.sizeR) then
       write(6,*)'erreur CST3 ',sizer,cst(3)
       stop
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
    integer:: size1,size3,sizeV,size9
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at
    ! v=vp,r=xpp
    ! u=eat,g=glangv;a=ax;s=sigat    
    integer::nvi,nvr,sizeI,sizeR,nvl,sizel,ivi,ivr,ivl,ibi,ibr,ibl,ic,ic2,csi,csr,csl,ib,ip
    integer, dimension (0:26):: Iposf,Rposf,Lposf
    integer,allocatable:: ibuffer(:)
    logical,allocatable::lbuffer(:)
    real(double),allocatable::rbuffer(:)
    integer:: cst(3)


    size1=atcf%imm;size3=3*size1; size9=3*size3
#ifdef PARA
    if (.not.present(caracT)) then
       carac='xfniewdlpvrugas'
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
  
  subroutine Eegal(atsource,atcible)
    class (atom_config),intent(in)::atsource
    class(atom_config)::atcible
    select type (atcible)
    type is (atom_config_e)
       select type (atsource)
       type is (atom_config_e)
          atcible%lax=atsource%lax
          atcible%lprteat=atsource%lprteat
          atcible%llangevin=atsource%llangevin
          atcible%lsigat=atsource%lsigat
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
    atcible%im_glob=atsource%im_glob
    atcible%imm_glob=atsource%imm_glob
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
!          write(6,*)'cible',atcible%imm,atcible%im,size(atcible%vp)
!          write(6,*)'source',atsource%imm,atsource%im,size(atsource%vp)
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
       end select
       call atsource%Eegal(atcible)
     end subroutine deftype
     
  !on inverse deux atomes dans la configuration
  subroutine switch_atom(atsource,ind_switch_1, ind_switch_2)
    class(atom_config)::atsource
    class(atom_config),allocatable:: intermediaire
    integer :: ind_switch_1, ind_switch_2,nag1,nag2
    logical :: lex = .true.
    call atsource%deftype(intermediaire)
    call intermediaire%init(imin=2)
    nag1=atsource%num_at_glob(ind_switch_1)
    nag2=atsource%num_at_glob(ind_switch_2)
    call atsource%copy_atom(ind_switch_1,intermediaire,1)
    call atsource%copy_atom(ind_switch_2,intermediaire,2)
    call intermediaire%copy_atom(2,atsource,ind_switch_1)
    call intermediaire%copy_atom(1,atsource,ind_switch_2)
    atsource%num_at_glob(ind_switch_1)=nag1
    atsource%num_at_glob(ind_switch_2)=nag2
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

    call at2pack%deftype(at)
    call at%init(imn,immn,at2pack%ltabvois,nvois,rvois,im_glob=at2pack%im_glob,imm_glob=at2pack%imm_glob)
    do i=1,immn
       call at2pack%copy_atom(i,at,i)
    end do
    call at2pack%dealloc
    call at%copy_config(at2pack,lrescl=.true.)
    

  end subroutine pack

  subroutine fab (atsource,atcible,lback,lrescl) ! construit atsource à partir de lgul de atcible , ecrase atcible
    class(atom_config),intent(in)::atsource
    class(atom_config),intent(out)::atcible
    logical::lback
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

       call atsource%Eegal(atcible)
       call atcible%init(imtrf,imtrf,atsource%ltabvois,nvois,rvois)
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
       stop
    end if

  end subroutine fab

  subroutine backto(atfab,atback)
    class(atom_config),intent(in)::atfab
    class(atom_config)::atback
    integer::i2,imtrf,i
    do i=1,atfab%im
       i2=atfab%num_at_glob(i)
       imtrf=atback%num_at_glob(i2)
       call atfab%copy_atom(i,atback,i2,lextend=.false.)
       atback%num_at_glob(i2)=imtrf
    end do
  end subroutine backto

    
  subroutine add2conf (atsource,atcible,lextend,ldealloc)
    class(atom_config),intent(inout)::atsource
    class(atom_config),intent(inout)::atcible
    logical,optional,intent(in)::ldealloc,lextend



    logical :: ldal ! par defaut pas de destruction de atsource
    logical :: lext ! par defaut on etend atcible si besoin

    integer::i2,i,imcib,imnew,immcib,imsrc,immsrc,immnew
    class(atom_config),allocatable:: atcor
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
             stop
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



  subroutine print(atin,i1,i2,unit,natg1,natg2,caracT)
    class(atom_config), intent(in)::atin
    integer,optional::i1,i2,unit,natg1,natg2
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac
    integer::i,im,ifin,ideb,ist,ifn,natpr,ig,iprt,unitw
    class (atom_config),allocatable::atprt
    unitw=6
    if (.not.present(caracT)) then
       carac='xfniewdlpvfrugas'
    else
       carac=caracT
    end if
    if (present(unit))unitw=unit

    call atin%deftype(atprt)
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
          
    atprt%im_glob=atin%im_glob
    atprt%imm_glob=atin%imm_glob
    
    write(unitw,*)'im = ',atprt%im
    write(unitw,*)'imm = ',atprt%imm
    write(unitw,*)'im_glob = ',atprt%im_glob
    write(unitw,*)'imm_glob = ',atprt%imm_glob
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
       atndm%llangevin=llangevin
       atndm%lax=lax
       atndm%lprteat=lprteat
       atndm%lsigat=lsigat
       call atndm%init(im,imm,ltabvois,nvois,rvois=0.d0)
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


  
  subroutine vers_master_atom(atcfloc,atcfcomp,div,caracT)
    class(atom_config),intent(in)::atcfloc
    class(atom_config)::atcfcomp
    type(para_config)::div
    character(len=*),optional,intent(in)::caracT
    character(len=26)::carac

#ifdef PARA
    integer::idmaster,idloc,icomm ! proc master,proclocal ,communicateur

    integer::nvi,nvr,sizeI,sizeR,nvl,sizel,ivi,ivr,ivl,ibi,ibr,ibl,ic,ic2,csi,csr,csl,ib,ip,iag,itot
    integer, dimension (0:26):: Iposf,Rposf,Lposf
    integer,allocatable:: ibuffer(:)
    logical,allocatable::lbuffer(:)
    real(double),allocatable::rbuffer(:)
    integer:: cst(3)
    type(mpi_communicator)::mpic

    logical,allocatable::lnag(:)
    integer::natgM
    integer,allocatable::inag(:)

    integer,allocatable::nag(:)
    integer::imloc,imloc3,iproc,imrecv,icomp,imrecv3,imrecv9,imloc9
    integer::imcomp,imtot,proc_source,npim,iloc
    logical,allocatable::mask(:)

    if (.not.present(caracT)) then
       carac='xfniewdlpvrugas'
    else
       carac=caracT//'np'
    end if
    idmaster=0
    idloc=div%mpi_image%rank
    npim=div%mpi_image%nproc
    mpic=div%mpi_image

    natgM=maxval(atcfloc%num_at_glob(1:atcfloc%im))
    call mpic%max(natgM)
    allocate(lnag(natgM))
    allocate(inag(natgM))
    inag=0
    lnag=.false.

    if (idloc==idmaster) then
       imtot=0
       do iproc=0,npim-1
          if (iproc==idmaster) then
             imrecv=atcfloc%im
             imtot=imtot+imrecv
             do iloc=1,imrecv
                lnag(atcfloc%num_at_glob(iloc))=.true.
             end do
          else
             if (allocated(ibuffer)) then
                deallocate(rbuffer);deallocate(ibuffer);deallocate(lbuffer)
             end if
             call  mpic%probe(11011,sourceout=proc_source)
             call mpic%RECV(imrecv,proc_source, 11011)
             imtot=imtot+imrecv
             allocate(nag(imrecv))
             call mpic%RECV(nag(1:imrecv), proc_source,10014)
             do iloc=1,imrecv
                lnag(nag(iloc))=.true.
             end do
             deallocate(nag)
          end if
       end do
       if (imtot.ne.atcfcomp%im) then
          write(6,*)'atomes perdus ?',idloc, imtot,atcfcomp%im,div%mpi_orig%rank
          call MPI_finalize(ierr)
          stop
       end if
       itot=0
       do iag=1,natgM
          if (lnag(iag)) then
             itot=itot+1
             inag(iag)=itot
          end if
       end do
       if (itot.ne.imtot) then  
          write(6,*)'atomes perdus 2 ?',itot,imtot
          call MPI_finalize(ierr)
          stop
       end if


    else
       imloc=atcfloc%im
       allocate(nag(imloc))
       call mpic%SEND(imloc,idmaster,11011)
       nag(1:imloc)=atcfloc%num_at_glob(1:imloc)
       call mpic%SEND(nag, idmaster, 10014)
       deallocate(nag)
    end if
       

    
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
                atcfcomp%num_at_glob(inag(nag(iloc)))=nag(iloc)
                atcfcomp%xp(1:3,inag(nag(iloc)))=atcfloc%xp(1:3,iloc)
                atcfcomp%fp(1:3,inag(nag(iloc)))=atcfloc%fp(1:3,iloc)
                atcfcomp%ityp(inag(nag(iloc)))=atcfloc%ityp(iloc)
                if( allocated(atcfloc%lgul))  atcfcomp%lgul(inag(nag(iloc)))=atcfloc%lgul(iloc)
                atcfcomp%proc_at(inag(nag(iloc)))=idmaster 
                select type(atcfloc)
                class is (atom_config_d)
                   select type (atcfcomp)
                   class is (atom_config_d)
                      atcfcomp%vp(1:3,inag(nag(iloc)))=atcfloc%vp(1:3,iloc)
                      atcfcomp%xpp(1:3,inag(nag(iloc)))=atcfloc%xpp(1:3,iloc)
                   end select
                end select
                select type(atcfloc)
                class is (atom_config_e)
                   select type (atcfcomp)
                   class is (atom_config_e)
                      if((atcfcomp%lsigat).and.(atcfloc%lsigat))then
                         atcfcomp%sigat(1:3,1:3,inag(nag(iloc)))=atcfloc%sigat(1:3,1:3,iloc)
                      endif
                      if((atcfcomp%lprteat).and.(atcfloc%lprteat))then
                         atcfcomp%eat(inag(nag(iloc)))=atcfloc%eat(iloc)
                      endif
                      if((atcfcomp%llangevin).and.(atcfloc%llangevin))then
                         atcfcomp%glangv(:,inag(nag(iloc)))=atcfloc%glangv(:,iloc)
                      endif
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
             do iloc=1,imrecv
                atcfcomp%num_at_glob(inag(nag(iloc)))=nag(iloc)
             end do

             call buffersizes (atcfloc,sizeI,sizeR,sizel,IposF,Rposf,Lposf,carac,imrecv)
             if (cst(1).ne.sizeI) then
                write(6,*)'erreur CST1B ',sizeI,cst(1)
                stop
             end if
             if (cst(2).ne.sizel) then
                write(6,*)'erreur CST2B ',sizel,cst(2)
                stop
             end if
             if (cst(3).ne.sizeR) then
                write(6,*)'erreur CST3B ',sizer,cst(3)
                stop
             end if
             allocate(Rbuffer(sizeR));allocate(ibuffer(sizeI));allocate(lbuffer(sizeL))
             if (cst(1).ne.0)call mpic%recv(ibuffer,proc_source,314)
             if (cst(2).ne.0)call mpic%recv(lbuffer,proc_source,315)
             if (cst(3).ne.0)call mpic%recv(Rbuffer,proc_source,316)
             call distribnag(nag,inag,imrecv,atcfcomp,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf)
             do iloc=1,imrecv
                atcfcomp%proc_at(inag(nag(iloc)))=proc_source !  ibuffer(1:imrecv)
             end do

             deallocate(nag)
             

          end if
       end do


       if (imtot.ne.atcfcomp%im) then
          write(6,*)'atomes perdus ?',idloc, imtot,atcfcomp%im,div%mpi_orig%rank
          call MPI_finalize(ierr)
          stop
       end if
       do icomp=2,imtot
          if (atcfcomp%num_at_glob(icomp).lt.atcfcomp%num_at_glob(icomp-1)) then
             write(6,*)'atomes mal rangés L2M ?', icomp,atcfcomp%num_at_glob(icomp),atcfcomp%num_at_glob(icomp-1)
          end if
       end do
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
    integer::idmaster,idloc,icomm,npim ! proc master,proclocal ,communicateur
    integer::imloc,imloc3,iproc,imrecv,ideb,ifin,imrecv3,imrecv9,icomp
    integer::imcomp,imtot,proc_source,ns,iloc,i,iu
    logical,allocatable::mask(:)

    integer::nvi,nvr,sizeI,sizeR,nvl,sizel,ivi,ivr,ivl,ibi,ibr,ibl,ic,ic2,csi,csr,csl,ib,ip
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

!!$    call atcfcomp%print(unit=500+mpic%rank)
!!$    call atcfloc%print(unit=600+mpic%rank)

    if (.not.present(caracT)) then
       carac='xfniewdlpvrugas'
    else
       carac=caracT//'np'
       
    end if



    if (idloc==idmaster) then
       allocate(mask(atcfcomp%imm))
       do icomp=2,imtot
          if (atcfcomp%num_at_glob(icomp).lt.atcfcomp%num_at_glob(icomp-1)) then
             write(6,*)'atomes mal rangés M2L ?', icomp,atcfcomp%num_at_glob(icomp),atcfcomp%num_at_glob(icomp-1)
          end if
       end do

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
                         atcfloc%xpp(:,iloc)=atcfcomp%xpp(:,i)
                      end select
                   end select
                   select type(atcfloc)
                   class is (atom_config_e)
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
             if (allocated(ibuffer)) then
                deallocate(ibuffer);deallocate(rbuffer);deallocate(lbuffer)
             end if
             call buffersizes (atcfcomp,sizeI,sizeR,sizel,IposF,Rposf,Lposf,carac,ns) ! ns car on transfère seulement ns atomes vers atcfloc
             allocate(ibuffer(sizeI));     allocate(Lbuffer(sizeL));     allocate(Rbuffer(sizeR));
             call buildbuff(atcfcomp,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf,mask,ns)
             csT(1)=csi
             csT(2)=csl
             csT(3)=csR
!             write(6,*)'SEND',div%image,iproc,ns,cst
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
!       write(6,*)'RECV',div%image,idloc,imrecv,cst
!       write(6,*)'atcfloc',atcfloc%im,atcfloc%imm,sizeI,sizel,sizer
       if (atcfloc%im.ne.imrecv) then
          write(6,*)'ERREUR M2L', atcfloc%im,imrecv
          stop
       end if
       if (cst(1).ne.sizeI) then
          write(6,*)'erreur CST1B ',sizeI,cst(1)
          stop
       end if
       if (cst(2).ne.sizel) then
          write(6,*)'erreur CST2B ',sizel,cst(2)
          stop
       end if
       if (cst(3).ne.sizeR) then
          write(6,*)'erreur CST3B ',sizer,cst(3)
          stop
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
    integer:: size1,size3,sizeV,size9,iat
    integer::nvi,nvr,nvl,ivi,ivr,ivl,ibi,ibr,ibl,ic,ic2,csi,csr,csl,ib,ip
    !    integer, dimension (0:26):: Iposf,Rposf,Lposf
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at   

    if (present(nmask))then
       nmaskV=nmask
    else
       nmaskV=atcf%imm
    end if

    size1=nmaskV;size3=3*size1; size9=3*size3

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
       if(scan('r',carac).ne.0) then
       nvR=nvR+1
        Rposf(nvR)=Rposf(nvR-1)+size3
       !call MPI_SEND(atcf%xpp, size3, NDM_MPI_REAL_DOUBLE, rgcib,110,comm,ierr)
       end if
    end select
    select type (atcf)
    class is  (atom_config_e)
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
    integer::nvi,nvr,nvl,ivi,ivr,ivl,ibi,ibr,ibl,ic,ic2,ib,ip
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
             stop
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
    end select
    select type (atcf)
    class is  (atom_config_e)
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
                      rbuffer(ib)=atcf%glangv(ic,ip)
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
                      rbuffer(ib)=atcf%ax(ic,ip)
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
    if (csi.ne.sizeI) then
       write(6,*)'erreur CSI 2',sizeI,csi
!       call endmpi
       stop
    end if
    if (csr.ne.sizer) then
       write(6,*)'erreur CSR ',sizeR,csr
!       call endmpi
       stop
    end if
    if (csl.ne.sizel) then
       write(6,*)'erreur CSL ',sizel,csl
!       call endmpi
       stop
    end if

    return
  end subroutine buildbuff

  
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
    integer::nvi,nvr,nvl,ivi,ivr,ivl,ibi,ibr,ibl,ic,ic2,ib,ip


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
    end select
    select type (atcf)
       class is  (atom_config_e)
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

    if (csi.ne.sizeI) then
       write(6,*)'erreur CSI 1 ',sizeI,csi
!       call endmpi
       stop
    end if
    if (csr.ne.sizer) then
       write(6,*)'erreur CSR ',sizeR,csr
!       call endmpi
       stop
    end if
    if (csl.ne.sizel) then
       write(6,*)'erreur CSL ',sizel,csl
!       call endmpi
       stop
    end if

  end subroutine copybuff


  subroutine distribnag(nag,inag,immax,atcf,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf)
    integer,intent(in),allocatable::nag(:),inag(:)
    integer, intent(in)::sizeI,sizel,sizer,immax
    integer, intent(in),dimension (0:26):: Iposf,Rposf,Lposf
    class(atom_config),intent(inout)::atcf
     character(len=26),intent(in)::carac
    integer,allocatable,intent(in):: ibuffer(:)
    logical,allocatable,intent(in)::lbuffer(:)
    integer,intent(out)::csi,csr,csl
    real(double),allocatable,intent(in)::rbuffer(:)

    integer:: size1,size3,sizeV,size9
    integer::nvi,nvr,nvl,ivi,ivr,ivl,ibi,ibr,ibl,ic,ic2,ib,ip


    size1=immax;size3=3*size1; size9=3*size3
    ibi=0;ibl=0;ibr=0
    ivi=0;ivl=0;ivR=0
    csi=0;csl=0;csr=0


    if(scan('n',carac).ne.0)  then
       ivi=ivi+1
       do ip=1,size1
          ib=Iposf(ivi-1)+ip
          atcf%num_at_glob(inag(nag(ip)))=ibuffer(ib)
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
    end select
    select type (atcf)
       class is  (atom_config_e)
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

    if (csi.ne.sizeI) then
       write(6,*)'erreur CSIC 1 ',sizeI,csi
!       call endmpi
       stop
    end if
    if (csr.ne.sizer) then
       write(6,*)'erreur CSRC ',sizeR,csr
!       call endmpi
       stop
    end if
    if (csl.ne.sizel) then
       write(6,*)'erreur CSLC ',sizel,csl
!       call endmpi
       stop
    end if

  end subroutine distribnag
  
end module atomconfig


