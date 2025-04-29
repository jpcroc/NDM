# NDM potential files   
Potential files are stored in NDM_ROOT/potentials with many examples.  
The type of the potentials are indicated in the .din file :  
  
## in the din file:  
For a single type of potential :  
ipotential=  
>**7 tabulated pair (with or without charge)   **    
 **10 EAM;**  
 **-10=LAMMPS atom style atomic;**  
 **-11 LAMMPS atom style charge (changes only simple.potin)**  
 **16= CRG**  
 **20=MILADY**  
 0=Born-Mayer-Huggins,  
1=Buckingham,  
2=watanabe,  
3=buck8,  
4=Morelon UO2,  
5  Morse term,  
9 Basak  
6=Stillinger Weber a  la Vashista  
  12 ZrC JuLi(+Tersoff Doan)  ;  
 13 Tersoff coupure COS;  
 14 Tersoff coupure FD ;  
 15 tersoff coupure SIN (original) ;  
 11 Ercollesi ;;  
 8 bandura 2017= Bukingham +Morse+Fermi-Dirac+Inverse gaussian  
   
  
###For multiple potential types, set :  
npotentiel = n : number of potentials if >1 Lpotential is used instead of ipotential   
lpotentiel(:)=.false. if npotential >1 lpotential indicates the used potential    
<div style="page-break-after: always;"></div>  
ntyp=-1 the number of types in the simulation must be specified ONLY if npotentiel >1. Other wise it is dealt with in the .potin file.    
  
  
  
## eamtab.potin  
this is a tabulated eam input file.  
### example  
../potentials/eamtab.potin_NiAl  
for ntyp types of atoms, pairs are enumerated as follows :  
1-1  
1-2  
.  
.  
1-n  
2-2  
2-3  
.  
.  
2-n  
3-3  
.  
.  
.  
n-n  
  
  
### line by line explanation  
 2 ! number of atom types. It is the actual number of types in the simulation if npotential=1 OR the number of types in this file if npotential>1. The explanation below is for npotential=1  
   6.30000000000       ! cutoff radius  
56.693 28 'Ni' ! mass ; atomic number of symbol of type 1   
26.982 13 'Al'  
  -1.0 -1.5 ! Ziegler joining points (in Ang) for pair 1-1 , negative values mean no Ziegler potential connection.  
  -1.0 -1.5  !  Ziegler joining points (in Ang) for pair 1-2 , negative values mean no Ziegler potential connection.  
  -1.0 -1.5  !  Ziegler joining points (in Ang) for pair 2-2 , negative values mean no Ziegler potential connection.  
  6300  ! number of points in the potential grids  
 1 ! type 1 GLUE/EAM embedding function  
 6300 7.957506165079365E-004 ! one restates explicitly the number of points and step size (ssz).  
  7.957506165079365E-004 -5.521848303000000E-002 !One starts at 1*ssz.  
  1.591501233015873E-003 -7.217284523679759E-002 !First columns is in arb.u., second in eV.  
  2.387251849523810E-003 -9.279867925370675E-002  
.  
.  
5.01243313338349        29.1296427056877     
   5.01322888400000        29.1408414200000     ! There are 6300 lines for type 1  
1 ! type 1 density  
 6300 1.000000000000000E-003   ! one restates explicitly the number of points and step size.  
  1.000000000000000E-003  0.157126463100000  !units are Ang and arb.u.   
  2.000000000000000E-003  0.157126463100000     
  3.000000000000000E-003  0.157126463100000     
.  
.  
6.29800000000000       0.000000000000000E+000  
   6.29900000000000       0.000000000000000E+000  
   6.30000000000000       0.000000000000000E+000  
2 ! type 2 GLUE/EAM embedding function  
6300 7.957506165079365E-004  
  7.957506165079365E-004 -1.132095561000000E-002  
  1.591501233015873E-003 -1.794826230533940E-002  
.  
.  
5.01322888400000        6.28007349000000     
2 ! type 2 density  
 6300 1.000000000000000E-003   
  1.000000000000000E-003  0.123929132300000     
  2.000000000000000E-003  0.123929132300000     
.  
.  
6.29900000000000       1.226568485036625E-010  
   6.30000000000000       1.567585409656496E-010  
1  ! pair repulsion 1-1  
 6300 1.000000000000000E-003   
  1.000000000000000E-003   4.82001809700000  !units are Ang and eV     
  2.000000000000000E-003   11164133.6668597     
.  
.  
6.29900000000000       0.000000000000000E+000  
   6.30000000000000       0.000000000000000E+000  
2  ! pair repulsion 1-2  
 6300 1.000000000000000E-003   
  1.000000000000000E-003   478.747452500000     
  2.000000000000000E-003   5187873.21607976     
.  
.  
6.30000000000000       0.000000000000000E+000  
3  ! pair repulsion 2-2   
 6300 1.000000000000000E-003  
  1.000000000000000E-003   6.89225935500000     
  2.000000000000000E-003   2410747.26079528     
.  
.  
<div style="page-break-after: always;"></div>6.30000000000000       4.497386693493672E-008  
  
  
##pair_tab.potin  
This is for tabulated pair potential, with or without additional Coulombic interaction.  
Ewald summation driving with ieawald :  
iewald=0 no coulombic term  
iewald=1 Coulombic term with regular Ewald Summation   
iewald=3 Coulombic term with regular Wolf Summation  
iewald=2 Coulombic term with regular Particle Mesh Ewald Summation (pretty sure it does not work)  
  
### ewald namelist  
It is assumed you know what a fortran namelist is.  
    namelist /ewald/ rue, alpha, precis, ncouc3, ncoucx, ncoucy, ncoucz, ecrue,ipotrep,RcWolf  
rue is the Cutoff radius in ANg.  MANDATORY  
alpha is the alpha parameter of the EWald summation.  
precis is the precision of the Ewald Summation  
ncouc3, ncoucx, ncoucy, ncoucz are thesize of the summations in reciprocal space with ncoucx=ncoucy=ncouz=ncouc3 if only this last parameter is specified.  
Only two of the last three parameters must be specifid, the last one is deduced from the first two.  
Rcwolf = rue by default is the cutoff for Wolf Summation.  
ipotrep=0 by default : if 0 no ZBL input . if 2 ZBL connection radii will be read  
  
Example for a non coulombic potential:  
 &ewald  
   rue=11.00 Cutoff radius in ANg.   
 &end  
Example for a Coulombic potential  
 &ewald  
   rue=11.00 Cutoff radius in ANg  
   iewald=3  
   alpha=0.4  
 &end  
  
  
### line by line  
LINE 1:  
0 .false. !integer input for Ewald summation and three-body terms. three body potentials have not been used in the 21st century...  
LINE 2:  
&ewald  
   rue=11.00  
   iewald=3,   alpha=0.4  
 &end  
LINE "3":  
2 ! number of atomic types  
LINE "4-5" WITHOUT  Coulombic term (no charge):  
15.9994 8.0 'O ' 1      !  mass; atomic number, symbol, type number  
238.03  92.0 'U ' 2    !   mass; atomic number, symbol, type number  
LINE "4-5" WITH Coulombic term/charge:  
-1.1104 15.9994 8.0 'O ' ! charge  mass; atomic number, symbol, type number  
2.2208 238.03  92.0 'U ' ! charge  mass; atomic number, symbol, type number  
LINE 6:  
3 50000   ! number of pairs to be read grid size (ngr). The number of pairs is specified as it may differ from ntyp*(ntyp+1)/2 in case of npotential>1  
LINE 7  if ipoteprep=0  :  
2 2 ! types of atoms forming the pair to be read  
LINE 7  if ipoteprep=1  :  
2 2 0.8 1.0 ! types of atoms forming the pair to be read  connexion distances for ZBL connection.  
LINE 8 to 8+ngr (see line 6)  
  0.234000000000E-03   519155925.182 distance (ANG) ; pair interaction (eV)  
  0.468000000000E-03   258716171.282  
  .  
  .    
1 2 ! types of atoms forming the pair to be read  
  0.234000000000E-03   519155925.182 distance (ANG) ; pair interaction (eV)  
  0.468000000000E-03   258716171.282  
  .  
  .    
1 1 ! types of atoms forming the pair to be read  
  0.234000000000E-03   519155925.182 distance (ANG) ; pair interaction (eV)  
  0.468000000000E-03   258716171.282  
  .  
  .    
<div style="page-break-after: always;"></div> EOF  
  
##Other NDM potentials  
Other potentials are extensions of the pair_tab format with analytic tabulation perfiomrd inside NDM with parameters read in the files.  
###CRG.potin  
ipotential=16  
 2 ! nb de types  
  &ewald  
   iewald=1, rue=11.00  
   alpha=0.420560 ,   ncouc3=4  
   ipotrep=0  
&end  
-1.1104 15.9994 8.0 'O ' 1  ! charge  mass; atomic number, symbol, type number    
2.2208 238.03  92.0 'U ' 2  ! charge  mass; atomic number, symbol, type number  
11001  ! Number of pojnts in the grid to build  
1 !First type  
0.69 106.856 1000.0 !G n rhomax= estimation of the maximum density on atoms of that type from other atoms   
2 !type U  
1.806 3450.995 100.0 !G n rhomax  
1 !pair1= O O  
0 0 820.283 0.3529 3.8843 0 ! D gam A rho C r0  
2 !pair 2=  O U  
0.6608 2.058 448.779 0.3878 0 2.381 ! D gam A rho C r0  
3 !pair 3 = U U  
0 0 18600 0.2747 0 0 ! D gam A rho C r0  
###Basak.potin  
3  .false.  
 &ewald  
   rue=11.0,ipotrep=2, alpha=0.420560  
  
  &end  
 2 !2 types d'atomes  
 -1.1104 15.9994 8.0 'O '      ! Q CM (masse);  numero atomique ; symbole  
2.2208 238.03  92.0 'U '    ! Q CM (masse);  numero atomique ; symbole  
3  
1 1 883.12 0.3422 3.996 0. 0. 0. 0.6 0.8 : read(lupotin,*) tt1,tt2, abasakr,rhor, cbasakr, Dr,betar,rstar  
1 2 432.18 0.3422 0. 0.5055 1.864 2.378 0.6 0.8  
2 2 187.03 0.3422 0. 0. 0. 0. 0.6 0.8  
<div style="page-break-after: always;"></div> EOF  
  
##LAMMPS and Milady potentials  
See LAMMPS_MILADY.md for details  
These external potentials need a simple;potin file to input basic data to NDM. No mixing between LAMMPS and/or Milady and/or NDM potentials is possible.  
###simple.potin  
 3 ! nb of types  
11.0 rue  ! cutoff distance;  this is not necessary but I input it anyway to have cells consistent with LAMMPS/MLD cutoffs.  
15.9994 8.0 'O ' -1.1104     !  mass; atomic number, symbol, type number  
238.03  92.0 'U'  2.2208     !  mass; atomic number, symbol, type number  
131.4  54.0 'Xe'  0     !  mass; atomic number, symbol, type number  
These last lines must be consistent with the LAMMPS/MLD data. no internal check is performed.  
  
  
