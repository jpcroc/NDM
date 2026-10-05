## Compilation babylon



### gfortran Dynamics kit
module load mkl/latest
 module load openmpi/5.0.9-gcc11.5.0
 cmake ..  -D NDM_OPT_DK-IO=ON
 make -j 12


### Compile NDM LAMMPS


git clone -b patch_4Jul2026 --depth=1 https://github.com/lammps/lammps.git

cd lammps/

mkdir build

cd build/

module purge

module load mpi/latest mkl/latest compiler/latest

cmake  -C ../cmake/presets/basic.cmake -DPKG_REAXFF=ON -DCMAKE_CXX_COMPILER=icpx -DCMAKE_C_COMPILER=icx -DBUILD_SHARED_LIBS=ON ../cmake

OR 

with  SMTBQ :

cmake -C ../cmake/presets/basic.cmake   -DPKG_REAXFF=ON   -DCMAKE_CXX_COMPILER=icpx -DCMAKE_C_COMPILER=icx   -DBUILD_SHARED_LIBS=ON -DPKG_SMTBQ=yes   -DCMAKE_C_FLAGS="-qmkl=sequential"   -DCMAKE_CXX_FLAGS="-qmkl=sequential"   -DCMAKE_EXE_LINKER_FLAGS="-qmkl=sequential"   -DCMAKE_SHARED_LINKER_FLAGS="-qmkl=sequential"   ../cmake

make -j 12


module load mpi/latest mkl/latest compiler/latest

export LAMMPS_HOME=$PWD

cd ../..

git clone -b ndm_lammps_2026 --depth=1 https://github.com/jpcroc/NDM.git

cd NDM

mkdir build

cd build

module load mkl

cmake ..  -D NDM_PACKAGE_LIST=LAMMPS  -DCMAKE_Fortran_COMPILER=mpiifx

make -j 12



