
## NDM input files

### name.in  
This file contains just one line with a string which is the head of the calculation specific input and output files. Example ;  
\> more name.in  
calc  
\>  
calc.gin and calc.gin/name.cin will then be read.  

### .gin file (e.g. name.gin)  
formatted/readable configuration file,   
L1 : repetition of the cell structure in the a, b, c directions  
Na Nb Nc  
L2-4 : cell parameters in Angstroms  
ax ay az ! first cell vector in an orthonoral direct base (X,Y,Z)  
bx by bz ! second cell vector in an orthonoral direct base (X,Y,Z)  
cx cy cz ! third cell vector in an orthonoral direct base (X,Y,Z)  
L5: number of atoms in the above cell  
Nat  
L6-L6+Nat : positions of the atoms in a,b,c reduced coodinates and type of a	atom (integer) to be connected to an actual type in the .potin file possible additional columns are not read  
u1 v1 w1 t1  
u2 v2 w2 t2  
.  
uNat vNat wNat tNat  
  
```  
  
### .cin file (e.g. name.cin)  
unformatted configuration produced by ndm as old_name.cout to be renamed as name.cin . Contains the box structure, positions, types and velocities of atoms, as well as additional parameters allowing a direct restart/continuation of a calculation. Often negative (default ) values indicate that the feature is not activated

### .din file (e.g. name.din)
General driving file for the MD calculation
This file is a fortran namelist with a few mandatory and many optional variables.
See specific help file 

### .potin file (e.g. CRG.potin, eamtab.potin, etc.)
Potential file. Eah potential form has its own input format.  
See specific help file : potential.potin.md.  

