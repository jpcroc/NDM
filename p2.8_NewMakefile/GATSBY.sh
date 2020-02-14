# to compile in serial on gatsby 
# with gcc compiler 
# no module are need and no module will conflict 
module purge
module load openmpi/gcc/4.0.2
make ndm_serial_gfortran

# with lammps and gfortran 
# one need to compile its own lammps library first

make ndm_lammps_serial_gfortran LAMMPS=1

# with openmpi and gcc
# module openmpi/gcc/4.0.2, the module openmpi/intel/4.0.2 will conflict if not removed first.

make ndm_mpi_gfortran PARA=1





# with intel compiler 
# intel/219.5.281 module is needed and no module will conflict
module purge
module load openmpi/intel/4.0.2
make ndm_serial 

# with lammps and intel
# one need to compile its own lammps library first

make ndm_lammps_serial LAMMPS=1

# with openmpi and intel 
# module openmpi/intel/4.0.2, the module openmpi/intel/4.0.2 will conflict if not removed first.

make ndm_mpi PARA=1


# LAST but not least NDM + ML + intel
make ndm_ml ML=1





