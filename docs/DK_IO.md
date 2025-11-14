## DynamicsKit Input/Output library



### Summary:

NDM can be compiled with the DynamicsKit IO library allowing new input/output formats to be read and exported. Formats can be selected by specifying the associated **igen** and **ivisu** numbers in the `.din` file. The file that NDM will read must then be named : `name.extension` where `name` is read from `name.in` and `extension` corresponds to the correct igen.



### Compile NDM with DK_IO:

```bash
#Works only with GNU compilers for now
mkdir build
cd build
cmake .. #Chech that NDM_OPT_DK-IO=ON or add -D NDM_OPT_DK-IO=ON
make
```



### Available formats:

| Format  | Input (igen)  | Input with speeds* ** (igen)  | Output (ivisu) | Output with speeds (ivisu) | Extension |
| :-----: | :-----------: | :---------------------------: | :------------: | :------------------------: | --------- |
| AtomEye |      11       |              21               |       11       |             21             | .cfg      |
| AtomEye |      12       |              22               |       12       |             22             | .xfg      |
| CASTEP  |      13       |              23               |       13       |             23             | .cell     |
|   CIF   |      14       |            &cross;            |       14       |          &cross;           | .cif      |
| DL_POLY |      15       |              25               |       15       |             25             | .CONFIG   |
|  GULP   |      16       |            &cross;            |       16       |          &cross;           | .gulp     |
| LAMMPS  |      17       |            &cross;            |       17       |          &cross;           | .lmp      |
|  VASP   |      18       |            &cross;            |       18       |          &cross;           | .POSCAR   |
|   XYZ   |      19       |            &cross;            |       19       |          &cross;           | .xyz      |



(*) dmtype 30, 32, 34, 33, 19, 35, 12 are not compatible with reading atoms' speeds from the input file. (See src/prog.F90)

(**) For dmtype 3, 30, 5, 11, 31, 32, 33, 21, 22, 2 the input speeds will be reset. (See src/initspeed.F90)



If `tinit` is specified in the `.din` file, the atom velocities will be read then rescaled to match the temperature `tinit` (K). If `tinit` <= 0, atoms' speeds will remain unchanged.

