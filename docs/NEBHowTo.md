# NEB calculations guide

## Principles
NDM can perform NEB/climbing NEB calculations in parallel. The number of procs must be a multiple of the number of intermediate images. 

## PREPARATION  
###Initial and final configurations  
Initial and final configurations must be preliminarily relaxed/quenched with LPERIOD=.FALSE. !!!  

Considering a "neb.din" run 
Initial and Final  configurations are in deb_neb.gin and fin_neb.gin respectively. They will NOT be evolved by NDM.  
One thus needs the following files :
neb.din, deb_neb.gin, fin_neb.gin, "potential".potin files (lammps potential files), name.in.  

Be careful to respect a soncistent order of atoms between deb and fin structures !  

### din file
The din file must contain  
dmtype=9 !this triggers NEB calculation  
fpstop=0.01 !stopping energy criterion (eV/ang)  
lperiod=.false. no periodic conditions on atomic positions  VERY IMPORTANT !

It is good to specify of number of images:  
npath = 15 ! 15 images of the neb is the default. That means 13 INTERMEDIATE configurations.  
  
The other options are optional.  
nebtype=2 NEB search is the default, nebtype=1 is for DRAG search  
maxneb = 700  the MAXimum of NEB steps  
lclimb=.false. ! set to true for climbin NEB  
nwclimb=3 ! starts the climbing at the second evaluation of forces (in VASP =1, in Henkelmann's paper is said to be "a few iterations").  
    
lPathFromGin = .FALSE.  if T : read initial path in gin files *.1.gin, *.2.gin, ... (NEB calculaion). useful for restart, not sure it still works.     
nebrelaxation=2   We relax all the atoms; if nebrelaxation==1 only the most "deplaced" atoms (deplaced more than deltaRmax Ang) are relaxed, the other are kept fixed. No real interest.  
deltaRmax=1.d-2  
neb_noise=0  0 without noise, 1 with noise. Useful to break symmetries.  
neb_noise_scale=0.001   this will affect the 4th digit  ! inputs some noise on the positions. 
i_neb_drag =5 makes 5 in hyperplan (drag like) relaxation between succesive calculation of the spring forces. This is Cosmin's version. "Pure" NEB is i_neb_drag=1. test to check if this is a better choice..  

## RUN
One can run either in sequential or in multiple of the number of intermediate images, provided the space decomposition is possible.  
Ewample with npath=15 :
mpirun -np 13 rundm90 or mpirun -np 26 rundm90 or mpirun -np 52 rundm90 etc.




