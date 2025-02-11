# NDM atoms, cells and boxes

## Principles
NDM stores configurations in derived types for atoms, cells (see below) and boxes. This document describes the derived types of these objects.

## ATOMIC CONFIGURATIONS
Atomic configurations are stored in atom_config(_*) derived types. There can be many different atomic configurations within a single run.

There are 4 levels of atom_config derived types :
atom_config for positions, types and forces ;
atom_config_d with velocities added ;
atom_config_e with optional quantities added (previous position, atomic energy, atomic stress, etc...)  ;
atom_config_* for local extensions of atom_config_e; e;g. atom_confih_neb, atom_config_mcgc, etc..


###atom_config type
This is the minimal type for atomic configurations. All components are always allocated, except for the neighbour list.
type(atom_config):: atcf
atcf can be a global configurations, a local(space) configurations or a specific configurations (i.e. only type 1 atoms).  

####  type atom_config 
     integer::im=0,imm=0,im_glob,imm_glob    
     !im is the number of atoms in the config. For paralell calculations and local confs this is the number of atoms in the current space process nb d'atomes sur le proc,  
     imm is the size of the arrays (imm>=im)   
     im_glob is the global number of atoms (over all procs)   
     imm_glob is the size of the global arrays (i;e. the imm of the global conf if it exists)   


     real(double),allocatable:: xp(:,:) ! %xp(1:3,1:imm) positions of the atoms   
     real(double),allocatable::fp(:,:)   ! %fp(1:3,1:imm) forces on  the atoms   
     integer,allocatable::ityp(:)  !%ityp(1:imm) types of the atoms   
     integer,allocatable::ielat(:)  ! index of the cell the atom belongs to  
     logical, allocatable:: lgul(:)  ! A general utility logical   
     integer,allocatable::num_at_glob(:)  ! a global index that is preserved when the atom changes processor   
     logical :: ltabvois  ! if T a neighbour table can be built   
     integer,allocatable:: iwmax(:) !  %iwmax (i) Index of the last neighbour of i   
     integer:: nvois ! size of the neighbour table (>=total number of neighbours)   
     real(double)::rvois ! radius for the neighbour table   
     integer,allocatable:: indi(:) ! %indi(1:nvois) index of neighbouring atoms  
     integer(long)::icaltabt   ! This is the index of the cell repartition must be identical for the corresponding cell configuration   
     integer,allocatable::proc_at(:) ! Array of size im_glob that indicates the processor that deals with the corresponding atom   
     integer::imf=0 ! index of the last ghost atom between om+1 and imf atoms are ghosts from other processors   
**    contains  **  
     procedure, pass::init=>init_atom_config  ! initilization 
     procedure, pass::copy_atom=>copy_atom_b  ! copy one atom from a config to another
     procedure, pass::dealloc=>dealloc_atom_config  ! deallcoates the config
     procedure, pass::vers_master=>vers_master_atom !(atcfloc,atcfcomp,div)   ! transfers a local config to a master process glocabl config  
     procedure, pass::master2loc=>master2loc_atom !(atcfcomp,atcfloc,div)    ! transfers  glocal config to a local config  
     procedure, pass::copy_config copy a config into another config   
     procedure, pass::print  !prints the config (for debug)  
     procedure, pass::send2proc=>s2p_atom  ! sends a config to a proc  
     procedure, pass::send2all=>s2a_atom  !sends a cofig to all procs  
     procedure, pass::recv=>rcv_atom  rceived a config from a proc  
     procedure, pass::zero=>zero_atom  ! sets all arraus to O but does not deallocate  
     procedure, pass::switch_atom  switch two atoms inside a config  
     procedure, pass::print_type   print the type of the config  
     procedure, pass::pack  ! rescales a config with many empty elements (when imm is way too large)  
     procedure, pass::fab  ! builds an extract config from another one based on a logical  
     procedure, pass::sort  ! builds a config from another one with true lgul followed by false lgul   
     procedure, pass::backto  ! transfers back a selection to an origin config (post%fab)  
     procedure, pass::add2conf  adds a conf to a nother one   
     procedure, pass::extend  extends the size of a config (with add2conf)   
     procedure, pass::deftype  !  ! initialize atcible to the type of atsource, inluding the values of lax, lpreeat, etc.   
procedure, pass::Eegal  stes the type of config based on a model  
  end type atom_config  


### atom_config extensions   
####   type, extends (atom_config):: atom_config_d ! with velocities. The most used config type in the code  
     real(double),allocatable::vp(:,:)  
   contains  
     procedure, pass::copy_atom=>copy_atom_d  
     procedure, pass::dealloc=>dealloc_atom_config_d  
     procedure, pass::zero=>zero_atom_d  
  end type atom_config_d  
  
####  type, extends (atom_config_d):: atom_config_e ! extended types with optional quantities  
     logical::lxpp   
     real(double),allocatable ::xpp(:,:)   ! previous positions (t-delta-t)    
     logical::lprteat 
     real(double),allocatable ::eat(:)  ! atomic energy
     logical::lsigat  
     real(double),allocatable ::sigat(:,:,:)   atomic stress (1:3,1:3;1:imm)   
     logical::lLangevin  
     real(double),allocatable ::Glangv(:,:)  ! Langevin factor  (1:3,1:imm)  
     logical::lax  
     real(double),allocatable ::ax(:,:)  ! initial positions   
   contains  
     procedure, pass::copy_atom=>copy_atom_e  
     procedure, pass::dealloc=>dealloc_atom_config_e  
     procedure, pass::zero=>zero_atom_e  
  end type atom_config_e
  
####  type, extends (atom_config_e)::atom_config_arps**  An example of extension
     real(double),allocatable::rho(:),fpr(:,:),fpg(:,:)  
     integer, allocatable::mov(:)  
     contains  
       procedure, pass::copy_atom=>copy_atom_arps   
       procedure, pass::dealloc=>dealloc_atom_config_arps   
       procedure, pass::zero=>zero_atom_arps  
  end type atom_config_arps


### atom_config routines
####  subroutine print_type(atcf,mess)  
Prints the actual type of the config  

class(atom_config),intent(in)::atcf  
    character(len=*),optional::mess  


####subroutine init_atom_config(atconf,imin,immin,ltabvois,nvois,rvois,lreallocate,im_glob,imm_glob)  
initilization of atconf

    class(atom_config),intent(inout)::atconf  
    integer,intent(in):: imin  ! the number of atoms   
    logical,optional, intent(in)::ltabvois   with a neighbour table (default=false)   
    logical,optional, intent(in)::lreallocate	forces a reallocation (default= false, i. e. stops il already allocated)  
    integer, optional::nvois,immin,im_glob,imm_glob  ! default values for nvois imm, im_glob and imm_glob (default set to zero except imm=imin   )
    real(double),optional::rvois  !
  
####  subroutine copy_atom_b(atsource,i,atcible,j,lextend,caracT)   

copies atom i from atsource to atom j of atcible   
class(atom_config), intent(in)::atsource  
    integer,intent(in):: i  
    class(atom_config), intent(inout)::atcible  
    integer,intent(in):: j  
    character(len=\*),optional,intent(in)::caracT  
    logical , optional, intent(in) :: lextend (default = false)  ! extends the target conf if needed  
*** caracT***:
 string driving what is to be copied  
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at, v=vp,r=xpp , u=eat,g=glangv;a=ax;s=sigat; m=mov(arps)   
  end subroutine copy_atom_b  

####  subroutine extend(atcf,iadd)


  extends the size of atconf arrays  of iadd elements
  class(atom_config),intent(inout)::atcf   
    integer,intent(in)::iadd    
  end subroutine extend  


####  subroutine send2proc (atcf, rgcib,mpic,caracT)   
sends mpi copy of  atcf to process rgcib of communicator mic(see para help)

    class(atom_config):: atcf  
    type(mpi_communicator),intent(in)::mpic  
    integer,intent(in)::rgcib  
    character(len=*),optional,intent(in)::caracT  
    if (.not.present(caracT)) then  
       carac='xfniewdlpvrugasm'  
  end subroutine send2proc



####  subroutine recv (atcf, rgem,mpic,caracT)
receives mpi copy of  atcf from process rgem of communicator mic(see para help)  
    class(atom_config):: atcf  
    type(mpi_communicator),intent(in)::mpic  
    integer,intent(in)::rgem  
    character(len=*),optional,intent(in)::caracT   
  end subroutine recv_atom  



####  subroutine zero_atom (atcf)  
sets arrays at 0, no deallocation  
class(atom_config):: atcf  

###  subroutine deftype(atsource,atcible)  
 initialize atcible to the type of atsource, inluding the values of lax, lpreeat, etc.
 
    class(atom_config),intent(in)::atsource  
    class(atom_config),allocatable::atcible  
     end subroutine deftype  

####subroutine Eegal(atsource,atcible)  
initialzes optional flags of  atcible to source (continuation of the preceeing routine)  
    class (atom_config),intent(in)::atsource  
    class(atom_config)::atcible  
  end subroutine Eegal  
  
####  subroutine copy_config (atsource,atcible,lrescl)
copies the whole atsource to atcible

    class(atom_config),intent(in)::atsource  
    class(atom_config)::atcible  
    logical,intent(in)::lrescl  ! if True escales atcible to match atsource  
  end subroutine copy_config  

####  subroutine dealloc_atom_config(atconf)
deallocates atconf

    class(atom_config), intent(inout)::atconf  
  end subroutine dealloc_atom_config  

     
####  subroutine switch_atom(atsource,ind_switch_1, ind_switch_2)  
switches ind_switch_1 and inswitch_2

class(atom_config)::atsource  
    integer :: ind_switch_1, ind_switch_2   
  end subroutine

####subroutine pack(at2pack,imm_in)  
reallocates the arrays to their minimal size  

    class(atom_config),intent(inout):: at2pack  
    integer,optional, intent(in):: imm_in  if present segs the size to imin
  end subroutine pack  

####  subroutine fab(atsource,atcible,lback,lrescl,commsp)
  builds atcible from elements of atsource with %lgul(i)=True  
  lback true indicates that the modified atcible will be returned to atsource
  
    class(atom_config),intent(in)::atsource  
    class(atom_config),intent(out)::atcible  
    logical::lback  
    logical, optional::lrescl  if true (default) rescales atcible  
    type(mpi_communicator),optional::commsp  a communicator to put the total number of elements pf atcible in im_glob  
  end subroutine fab  



####  subroutine backto(atfab,atback)
copies back atfab builty with fab to atback

    class(atom_config),intent(in)::atfab
    class(atom_config)::atback
  end subroutine backto

  
####  subroutine sort (atsource,atcible) ! construit atsource à partir de lgul de atcible , ecrase atcible

builds atcible from atsource 
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
    write(6,*)'AAAAA',rang,atin%im,atin%imf
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
    write(6,*)'borders',i1l,i2l,natg1l,natg2l
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
             if(scan('r',carac).ne.0)then
                   write(unitw,'(A,2i9,3E15.7)')'%xpp= ', i,atin%num_at_glob(i),atin%xpp(:,i)
             end if
          class is (atom_config_e)
             write(unitw,*)'prt_e'
             if(scan('v',carac).ne.0)then
                   write(unitw,'(A,2i9,3E15.7)')'%vp= ', i,atin%num_at_glob(i),atin%vp(:,i)
             end if
             if(scan('r',carac).ne.0)then
                   write(unitw,'(A,2i9,3E15.7)')'%xpp= ', i,atin%num_at_glob(i),atin%xpp(:,i)
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
                 class is (atom_config_e)
                    select type (atcfcomp)
                    class is (atom_config_e)
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
  

  
end module atomconfig







