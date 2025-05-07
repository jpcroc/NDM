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
#####    contains        
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
      
####  type, extends (atom_config_e)::atom_config_arps   An example of extension    
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
    
#### subroutine Eegal(atsource,atcible)      
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
  end subroutine sort      
      
     
####  subroutine add2conf (atsource,atcible,lextend,ldealloc)    
  adds atsource to atcible      
    class(atom_config),intent(inout)::atsource      
    class(atom_config),intent(inout)::atcible      
    logical,optional,intent(in)::ldealloc  ! false by default = no deallocation of atsource       
    logical,optional,intent(in)::lextend  ! true by default= atcible allocations are extended if need to accomodate atsource      
  end subroutine add2conf    
    
    
    
####  subroutine print(atin,i1,i2,unit,natg1,natg2,caracT,mess,iter)  prints elements of atin configurations. For debugging only.      
    class(atom_config), intent(in)::atin  the configuration to be (partially) printed.      
    integer,optional::i1,i2 first and last elements of the arrays to be printed      
    integer,optional::natg1,natg2 first and last num_at_glob indexes  to be printed        
    integer,optional::unit unit to print to      
    integer,optional::iter  an integer to denote the "ietratio" the code is in          
    integer::i1l,i2l,natg1l,natg2l      
    character(len=*),optional,intent(in)::caracT    
    *** caracT***:    
 string driving what is to be copied      
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at, v=vp,r=xpp , u=eat,g=glangv;a=ax;s=sigat; m=mov(arps)       
    character(len=*),optional::mess a message to print      
  end subroutine print    
      
####  subroutine vers_master_atom(atcfloc,atcfcomp,div,caracT)    
sends the local configuration atcfloc to the global configuration atcfcomp on the ran 0 of "div" (div%mpi_image%rank==0)      
    class(atom_config),intent(in)::atcfloc      
    class(atom_config)::atcfcomp      
    type(para_config)::div      
    character(len=*),optional,intent(in)::caracT      
  end subroutine vers_master_atom    
    
      
####  subroutine master2loc_atom(atcfcomp,atcfloc,div,caracT)      
sends (from div%mpi_image%rank==0) the elements of the global atcfcomp configurations to the local atcfloc      
The receiving proc of each atom is atcfcomp%proc_at(1:imcomp)==iproc      
    class(atom_config)::atcfloc      
    class(atom_config)::atcfcomp      
    type(para_config)::div      
    character(len=*),optional,intent(in)::caracT      
    !x=xp;f=fp,n=num_at_glob,,i=ityp,e=ielat,w=iwmax,d=indi,l=lgul p=proc_at      
    ! v=vp,r=xpp      
    ! u=eat,g=glangv;a=ax;s=sigat      
    
  end subroutine master2loc_atom      
     
    
####  subroutine buffersizes (atcf,sizeI,sizeR,sizel,IposF,Rposf,Lposf,carac,nmask) is private      
  end subroutine buffersizes      
    
####  subroutine buildbuff(atcf,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf,maskR,nmask) should be private      
  end subroutine buildbuff      
    
      
####  subroutine copybuff(atcf,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf,immax) should be private      
    
  end subroutine copybuff      
    
    
####  subroutine distribnag(nag,inag,immax,atcf,sizeI,sizel,sizer,ibuffer,lbuffer,rbuffer,csi,csl,csr,carac,iposf,rposf,lposf) should be private      
  end subroutine distribnag    
      
    
      
end module atomconfig    
    
    
## BOX CONFIGURATIONS    
    
  ###type box_config    
     real(double):: at(3,3) ! at(:,1) is the first basis vector    
     real(double):: bg(3,3) ! recip from at at(:,i).bg(:,j)=delta(i,j)    
     real(double):: as(3,3) ! = b/norm2(b)     
     real(double):: volume of the box    
     real(double):: normat(3) norms of the box vectors    
     real(double):: normbg(3)  norms of the reciprocal vectors of the box    
     integer(long)::icaltabt index of the last cell re  partition should be indentical for the cell config    
     logical::lperiod positions are all "inside" the box  or not      
     logical::islarge=.true. ! useless for now     
     integer::ipbc(3) ! if 1 (default)the periodic boundary conditions are applied on planes  b-c,a-c,a-b      
   contains    
     procedure, pass::print=>boxprint   !prints the box   
     procedure, pass::init=>initbox  !inits the box componenents from at(3,3)  
     procedure, pass::showtype=>boxshowtype  ! idicates wether type is box_config or box_config_lpr  
     procedure, pass::master2slave=>boxmaster2slave    
     procedure, pass::send2proc=>boxsend2proc    
     procedure, pass::recv=>boxrecv    
  end type box_config    
    
  ###type, extends (box_config):: box_config_lpr  
  derived type which contains variables in the Parinello-Rhaman format, the only interesting variable is hdot (3,3) time derivative (velocity) of the box   
     real(double), dimension(3,3)  :: h, hDot    
     real(double), dimension(3,3)  :: trh, invh, invtrh, Gmat, invGmat, Gdot    
     real(double) :: invVolu,wbox    
         
  
  end type box_config_lpr    
  
  
###contains  
  
####  subroutine boxrecv(box,rgem,mpic)  
  receives in present rank the box rgem communicator in mpic  
    type(mpi_communicator),intent(in)::mpic  
    class(box_config)::box  
    integer,intent(in)::rgem  
  end subroutine boxrecv  
    
####  subroutine boxsend2proc(box,rgcib,mpic)  
  sends  the box from present rank to rgcib communicator in mpic  
    type(mpi_communicator),intent(in)::mpic  
    class(box_config),intent(in)::box  
    integer,intent(in)::rgcib  
  end subroutine boxsend2proc  
  
  
####  subroutine boxmaster2slave(box,rgem,mpic)  
  broadcasts the box from rank rgem to all communicators in mpic  
  
    type(mpi_communicator),intent(in)::mpic  
    class(box_config)::box  
    integer,intent(in)::rgem  
  end subroutine boxmaster2slave  
     
####  subroutine boxshowtype(box) ! shows the box type of the object  
    class(box_config)::box  
  end subroutine boxshowtype  
            
####   subroutine initbox(boxnew,at,ipbc,zl)  
initilizes the box according to at(3,3) (or zl(3) tetragonal leght, defective) and ipbc  
    class(box_config),intent(inout)::boxnew  
    real(double),intent(in),optional::at(3,3)  
    real(double),optional,intent(in)::zl(3)  
  end subroutine initbox  
  
####  subroutine updatebox(boxnew,at,zl,check)  
  updates the components of the box when at has changed  
    USE Mat_utils_mod,only:  MatInv  
    class(box_config),intent(inout)::boxnew  
    real(double),intent(in),optional::at(3,3)  
    real(double),optional,intent(in)::zl(3)  
    integer,optional::check  
  end subroutine updatebox  
  
  
####  subroutine boxprint(boxprt,unit,mess)  
    class(box_config),intent(in)::boxprt  
   end subroutine boxprint  
  
####   subroutine periodbox(box,atcf)  
   applies periodid boundary conditions (according to ipbc(3) on atcf positions, incluidng if present xpp and ax  
     class(box_config),intent(in)::box  
     class(atom_config)::atcf  
  
   end subroutine periodbox  
  
  
 end module boxconfig  
  
  
  
    
    
    
