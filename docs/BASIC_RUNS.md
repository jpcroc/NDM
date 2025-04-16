# NDM basic Runs  

## Atomic Configurations 
Configurations can be read from formatted (.gin) or unformatted (.cin)  files.
In the .din file :
igen =0 implies reading of formatted .gin file
igen =1 implies reading of formatted .cin file	
See INPUT_FILES for details about gin/cin files.
imm must be larger or equal to the number of atoms in the simulation


## Quenching/Relaxation

In the din file specify one of the below values of dmtype. I suggest **dmtype = 23 or 24 (Fire algorithm)**.  
   2 -> fast quenching algorithm (==21) 
   21 -> fast quenching constant volume  
   22 -> fast quench constant pressure   
   23 ->  fire quench algorithm constant volume  
   24 -> fire quench constant pressure  

It is useful to set in the din file the format of the output for atomic positions:  
lwgin =.true. will output a new gin file and/or
ivisu=4 output format for atomic configurations ivisu=4= CFG ; ivisu=5 new gin file


Consider setting :
lperiod=.false. That allows the atomic positions to "get out" of the box. This is only apparent as periodic Boundary Conditions on forces are always enforced, but this can be useful if one wants to compare two configurations (e.g. before and after relaxation; initial or final NEB images). NEB preparation MUST be performed with lperiod=.false.  


## Thermalization/equilibration/NVE/NPH/NPT runs
### DMTYPE
In the din file specify one of the below values of dmtype. I suggest **dmtype = 4**  
   1 -> Standard Verlet  algorithm
   4 -> Velocity Verlet  algorithm
  
Consider setting :
Tinit=300.0 ! will initialize the velocities at 300K feither from random distribution (gin file or cin file with lvpread=.false.) or rescaled to 300K (cin file with lvpread=.true. which is the default).  

### NVE
NVE runs : no additional variable.

### NVT
NVT runs :
In the din file  
**Text=300.0**  the external temperature, mandatory for cst temperature algorithms  
Choice of the constant T  algorithm. Select one of :  
lTcon=.false. if true isokinetic algo   
lTberendsen=.false.  T: Berendsen cst temp    
lTNose=.false. T: Nose  cst temp   
lTHoover=.false.  T: NoseHoover  cst temp     
lTandersen=.false. temperature constante a la Andersen   
lLangevin=.false. Langevin MD for cst temperature
  
See DINfile.md for additional optional options pertaining to each algo.  

### NPH
In dyn file specify:  
**lpr=.true.**  
By default the target stress/pressure is zero. Otherwise :  
sigext = 0.0 Symetric tensor of  the external stress  in kbar  

Additionnal variables can be :  
lpcon2 = .FALSE. if true a damping term is added on the dynamics of the box. Proves useful to stabilize the oscillation of the  box.  
lpconxyz = .FALSE. T: the relaxation is allowed only along the X, Y and Z axis  (an orthorhombic box remain orthorhombic).  
lpconx = .FALSE.  T: the relaxation is allowed only along the X axis  
lpcony = .FALSE.  T: the relaxation is allowed only along the  Y  axis  
lpconz = .FALSE.  T:  the relaxation is allowed only along the  Z axis  
lpcube=.false. !  T: volume relaxation only (a cube remains a cube).  
sigstop= 0.05 kbar is tbe consergence criterion for stress.  

### NPT
Combine lpr and a cst temperature algo. Any combination should work (hopefully).
Example:  
dmtype=4  
ltlangevin=.true.  
text=300.0  
lpr=.true.  
  
Again, it is useful to set specify in the din file :  
lwgin =.true. will output a new gin file and/or  
ivisu=4 putput format for atomic configurationsivisu=4 CFG, ivisu=5 = newgin file.


