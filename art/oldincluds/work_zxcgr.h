 module work_cg     
       type (int_vect3d)  , dimension(:)  , allocatable :: geometric_constraint             
       double precision   energ_tot 
       double precision   unit_vect(3,3)        
       double precision , dimension(:)    , allocatable ::  rho_site,V_site     
       double precision , dimension(:,:)  , allocatable ::  tau,force,force_relax
       double precision , dimension(:,:,:)  ,	allocatable :: Rn
       integer  nat_up,iflag,NN_MAX(3)
       integer         , dimension(:)    , allocatable ::  ndir
       integer         , dimension(:,:)  , allocatable ::  neigh_type       
       COMMON /param_work_geom_cg/unit_vect,tau,Rn,Rmax,NN_MAX,ndir,neigh_type         
       COMMON /param_work_dim_cg/nat_up,ndir_max          
       COMMON /param_work_pot_cg/V_site,rho_site
       COMMON /param_energ_cg/energ_tot
       COMMON /param_forces_cg/force,force_relax
       COMMON /param_constraint_cg/geometric_constraint 
end module work_cg       	
