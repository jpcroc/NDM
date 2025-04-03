#PATH Monte Carlo calculations

## Principles

### algorithm
NDM implements Path Monte carlo Calculations as explained in
Work-biased path-sampling calculations of chemical potentials: Principles and applications to uranium oxide,  
O. Barbour, J. P. Crocombette, T. Beigbedder, J. Tranchida, E. Bourasseau and M. Athènes,  
Journal of Chemical Physics 162, (2025) 

The idea is to step by step insert/delete particles and make a Monte-carlo sampling on the insertion/deletion paths based on the work needed to perform the insertion.

The standard insertion proceeds at random positions in the box (possibly with a cutoff inter-atomic distanc, see below).

Advanced insertion considers a spring interaction to either :
-a point;  
-a sphere;
-a line,;
-a plane.
This allows to perform calculations close to vacancy/bubble/dislocation/ grain boundary respectively.  

Alternatively, a bias can be input for insertions or delations; In this case there is no spring but a bias  in the positions of insertions or the choice of deletions.  


A semi grand canonical version has been coded but has not been fully tested.

### paralelization
The code can run in serial. However multiple paralelization levels make it rather parallel.  
The lowest level is the space division of the simulation box  
The intermediate  level is the division between the system with N and N+1 atoms ;  
the upper level is the simulatneous insertion of multiple (nparapath)  atoms.  

By default, the intermediate level is activated first, then the upper one provided lparapath is set to .true., and finally the space parallelization is activated.
The number of tasks must be a multiple of 2*nparapath.

##din variables

    pas_lambda_mc : number of intermediate steps in the insertions/deletions , MUST be specified  
    n_path : number of insertion/desinsertion paths MUST be specified  
    distminat ! minimum interatomic distance of insertion distance , default 1Ang  
    itypcalc  type of the inserted:deleted atoms MUST be specified  
    nparapath,  performs multiple simultaneous inerstions/desinsertions  and choose one (or stays on the previous path)  
    lparapath, runs the code in parallel along naprapath  
    idirectionmcgc,   default 0 , starts with insertions  

    ins_typ  ! 0 everywhere, with bias) : 1 in a sphere, 3 in a slice ; with a  spring : 11 in a site, 33 in a plane , 44 in a line , 55 in a bubble  
    lspring=.false. ! MCC calculation with a slowly vanishing vabishing string
    k_spring  spring stifness in Ev/ANG**2 (default 1)
    
    R0mcgc  radius of the sphere of biased insertion
    fdfactmcgc by d"efault 18, pilots the sharpness of the fermi-dirac function of insertgion in the sphere
    bublcenter position of the center of the bubble OR the spring attachment point for a 3D spring. Expressed in box reduced coordinates  
    izlins specifies the direction of  the plane/line to attach the spring.  
    zlcenter(3) specifies the position of the plane/line to attach the spring.  
    --For ins_type= 3/33 only one non zero component (along izlins) indicating the poition of the plane perpendicaular to a or b or c.  
    --For ins_type= 44 two  non zero component (NOT izlins) indicating the intersection of the line with the plane perpendicaular to a or b or c (indicated by the zero component)  

    protocol_mcc  = "MCP" or "cos" MCP performs the insertion linearly. "cos" performs the insertion as in (1-cos((pi/2)*step/pas_lambda_mc)  

    lbiais_inser  default .false. indicates that a bias will be applied on the choice of deletions  
    lbiais_retrait  default .false. indicates that a bias will be applied on the positions of  insertions   
    fdmc_1 numbers piloting the delation bias; As of the time of writing, this works only for fluorite UO2 with a bias on the deletion of O interstitials; needs addition coding to work in other cases  
    fdmc_2  
    nbatplus  insertions of nbatplus atoms at the same time default=1