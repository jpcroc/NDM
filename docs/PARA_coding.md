
# NDM PARALELIZATION

## Principles
NDM is paralelizeed in MPI. It uses an internal module for MPI calls. The basic paralelization (communicator comm_space) is a domain/force calculation decomposition. Additional levels of paralelization can be "easily" added above this level.


##MPI routines  
MPI uses an internal module for MPI calls. There should be no MPI calls outside his module.
###Name of the module : TPara
Object types within the module:
type mpi_communicator
type para_space_config

###mpi_communiccator
this module is a ganeral utility MPI encapsulator.

####Example of use
 call comm_space%send(atdml%xp(1:3,1:im),0,11003)   ! sends xp to rank 0 of comm_space
 call mpi_world%bcast(0,iseed)  brodcast iseed from 0  
 call comm_space%sum(imtot)  sums imtot over all ranks  
 call comm_space%max(natgm)  returns max value of natgm  
 call mpi_world%init(MPI_COMM_WORLD) !initializes the NDM  mpi_communicator object from the actual MPI comm  

Note that there is no need to specify the size, nor the shape, nor the type of the quantity, which can be an array or a value ; BUT of intrinsic fortran types

####  type mpi_communicator
     integer    :: comm       ! MPI communicator index
     integer    :: nproc      ! number of procs in the communicator 
     integer    :: rank       ! index           in the communicator 
     integer    :: group       ! group          of the communicator 
####   contains ! procedures are all "pass" 
     procedure :: init => mpic_init
     procedure :: init0 => mpic_init0
     procedure :: probe => mpic_probe
     procedure :: barrier => mpic_barrier
     procedure :: print ! gives info
     generic :: send  => mpic_send_dp,mpic_send_cdp,mpic_send_i,mpic_send_l
     generic :: recv  => mpic_recv_dp,mpic_recv_cdp,mpic_recv_i,mpic_recv_l
     ! sum
     generic :: sum  => mpic_sum_dp
     generic :: sum  => mpic_sum_cdp
     generic :: sum  => mpic_sum_i
     procedure :: mpic_sum_dp
     procedure :: mpic_sum_cdp
     procedure :: mpic_sum_i
     procedure:: mpic_send_dp,mpic_send_cdp,mpic_send_i,mpic_send_l
     procedure:: mpic_recv_dp,mpic_recv_cdp,mpic_recv_i,mpic_recv_l
     ! min ALLREDUCE !!
     generic :: min  => mpic_min_dp
     generic :: min  => mpic_min_i
     procedure :: mpic_min_dp
     procedure :: mpic_min_i
     ! max ALLREDUCE !!
     generic :: max  => mpic_max_dp
     generic :: max  => mpic_max_i
     generic::maxloc=>mpic_maxloc_dp
     procedure::mpic_maxloc_dp
     procedure :: mpic_max_dp
     procedure :: mpic_max_i
     ! and
     procedure :: and => mpic_and_l
     ! broadcast
     generic :: bcast  => mpic_bcast_dp
     generic :: bcast  => mpic_bcast_l
     generic :: bcast  => mpic_bcast_i
     generic :: bcast  => mpic_bcast_cdp
     procedure :: mpic_bcast_dp,mpic_bcast_l
     procedure:: mpic_bcast_i
     procedure :: mpic_bcast_cdp
     generic :: build=>build_i,build_dp,build_cdp,build_l
     procedure :: build_i,build_dp,build_cdp,build_l

  end type mpi_communicator


### PARA SPACE
This type defines the values of the space communicator. Regular programmers should not have to deal with this type; See Tpara if you have to.


##multi-level Division of comm_world

### Principles
comm_world can be divided into multiple communicators which can be further divided down to comm_space

### Paraconfig type and module 
objects of type paraconfig deal with the division of a communicator into smaller communicators.  
The outer (larger) communicator is mpi_orig  
The inner communicators are mpi_image. they all have a "master" process.  
The master processes are joined by the mpi_master communicator  

  type :: para_config  
     type(mpi_communicator)::mpi_image,mpi_master,mpi_orig   
     integer:: image ! rank of the image (from 0 to nimage)  
     integer::nimage ! = %mpi_master%nproc  
     logical :: lmaster ! .true. if the process is the  master process of the mpi_image comm  
   contains  
     procedure, pass::print  
  end type para_config  

####Note :
para_config%nimage=mpi_master%nproc  
para_config%nimage*mpi_image%nproc=mpi_orig%nproc  

### Example of a division :
with parapath of type para_config  
####Building of the mpi_orig communicator  

parapath%mpi_orig%nproc=nprocs 
 parapath%mpi_orig%rank=rang  
       call MPI_COMM_DUP(MPI_COMM_WORLD,parapath%mpi_orig%comm,ierr)  
       call MPI_COMM_GROUP(parapath%mpi_orig%comm,parapath%mpi_orig%group,ierr)  
#### definition of the number of divisions  
 parapath_mpi_orig is to be divided in nparapath subcommunicators:  
   parapath%nimage=nparapath  
   call commconstr(parapath)  


That last line builds the division and initializes the parapath mpi_comminicators objects.  

####One can go further down :
with paramcgc a subdivison of parapath  
paramcgc%mpi_orig%nproc= parapath%mpi_image%nproc ! =parapath%mpi_orig%nproc/nparapath  
    paramcgc%mpi_orig%comm= parapath%mpi_image%comm ! =parapath%mpi_orig%nproc/nparapath   
    paramcgc%mpi_orig%rank=parapath%mpi_image%rank  
    paramcgc%nimage=2  
    call MPI_COMM_DUP(parapath%mpi_image%comm,paramcgc%mpi_orig%comm,ierr)  
    call MPI_COMM_GROUP(paramcgc%mpi_orig%comm,paramcgc%mpi_orig%group,ierr)  
    call commconstr(paramcgc)


####Finally one must redefine comm_space ad the paramcgc_mpi_image:  
    call MPI_COMM_free(mpi_comm_space,ierr)  
    MPI_COMM_space=paramcgc%mpi_image%comm  
    call comm_space%init(MPI_COMM_SPACE)  
    nprocspace=paramcgc%mpi_image%nproc  


One thus has 3 levels of para  
    lmegamaster=.false. ;     if (parapath%mpi_orig%rank==0) lmegamaster=.true.  
lbigmaster=parapath%lmaster  
    lmaster=paramcgc%lmaster  






##Foces calculations with multiple levels of para : "master and partner"  
In reguler MD runs there is no real master. However, some algorithms (e.g. conjugate gradient, NEB, montecarlo) require global interventions on atomic positions. This is done by a master which has all the positions and plays with them. The other processes (partners)  are needed just to calulate the forces in the space/force  paralelization level. This is done with:
These calls should be done for all processes (master and partners).  


call initloc(atconf,cells,atmcgcloc,cellmcgcloc,boxmcgc,paramcgc,rumax,lperiod&  
            &,psc=pscgc,ldistrib=ldistrib,lcalcvois=lcalcvois,lboxchange=lprahman) !initloc contient caltabtc sur atloc


call pointer_caltabt_calfo(sig,potist,atconf,cells,boxmcgc,atmcgcloc,cellmcgcloc,paramcgc,&  
            &lperiod,lupdate=lchange,psc=pscgc,lcalcvois=lcalcvois,lboxchange=lprahman)  

Above atconf and cells are global atomic configurations (see atoms and cells) atmcgcloc and cellmcgcloc are pointers of the same types for local (with only local atoms) configurations. They do not have to be initialized or dimensionned, initloc does that. 
  class(atom_config),pointer::atmcgcloc  
  type(cell_config),pointer::cellmcgcloc
  
On return atconf contains the forces. pointer_caltact_calfo deals with cell repartition, and force calculations.

### transfer of atomic configurations from one process to another (see atoms and cells):  
 call atconf%send2proc(rgcib,commexample,'x')  
 This send the atconf configuration  (actually only the positions ('x') from local process to rgcib process (e.g. rgcib=0 <==> master) in communicator "commexample"  


## Master/slave scheme
A master/slave scheme is an algo where slaves wait in a storage routine and are activated by a call from a master, perform a task (force calculation), return data to the master and return to the storage routine. The master can call the slaves from many different places in the master code. that allows to program the master as if its process was sequential. Only from time to time does it call the slaves to work.  An example is in NGC.F90:


    gcpara being a division (not a communicator)		
type(para_config)::gcpara     
    USE parautils,only:WORKER_TAG,tolstoi,STOP_TAG  
    USE parautils,only:depeche_mode	    

    call set_pointers_gc ! initilisations of  pointers for  tolstoi and calfo
    if (gcpara%lmaster.neqv..true.) then
       call tolstoi (WORKER_TAG,gcpara)
*** the slaves are sent to storage routine***   
    else   
*** this is the master which will perform analyses       ***   
    .  
    .  
    .  


At some point the master enters "depeche_mode":

    call depeche_mode (gcpara)

That triggers the force calculations by master and slaves. The depeche_mode call can be done in any routine. 
    
*** masters calls back slaves to stop the code***  
       call tolstoi (STOP_TAG,gcpara) ! make servants return   
    end if

Some pointers must be defined to get the data see     call set_pointers_gc


