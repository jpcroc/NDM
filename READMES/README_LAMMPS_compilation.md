# How to compile LAMMPS as a library to use with NDM

On Gatsby one version is installed already and accessible with:
```
module purge
module load lammps/29Sep2021-u2
mkdir build_lammps ou rm -rf build_lammps/*
cd build_lammps
cmake ../ D NDM_PACKAGE_LIST="LAMMPS"
make -j 


```
It is also possible to compile LAMMPS from sources to use different versions or packages as explained below.

## Local LAMMPS compilation

To compile LAMMPS start with downloading the sources:
```
wget https://download.lammps.org/tars/lammps-stable.tar.gz
tar -xvf lammps-stable.tar.gz
cd lammps-2Aug2023
```
Older stable versions of LAMMPS are also available at this adress: https://download.lammps.org/tars/

LAMMPS compilation can then be configured with command line options but it is often more convenient to use cmake preset files.
See below examples of preset files to compile LAMMPS as a library. LAMMPS normally provides more examples in the cmake directory.

To compile LAMMPS with cmake, replace <preset.cmake> with the name of a preset file and run the following command lines in the lammps directory :

```
mkdir build; cd build
cmake ../cmake/. -C <preset.cmake>
make -j

```

For more information see LAMMPS documentation to build with CMake: https://docs.lammps.org/Build_cmake.html

## Exemple 1 of preset.cmake to compile LAMMPS 3Mar20 serial

```
set(BUILD_LIB ON CACHE BOOL "compile lammps as a library" FORCE)
set(BUILD_MPI OFF CACHE BOOL "no mpi" FORCE)
set(BUILD_OMP OFF CACHE BOOL "no omp" FORCE)

# Packages
set(PKG_MANYBODY ON CACHE BOOL "" FORCE)
set(PKG_USER-ATC ON CACHE BOOL "" FORCE)
set(PKG_CORESHELL ON CACHE BOOL "" FORCE)
set(PKG_USER-DIFFRACTION ON CACHE BOOL "" FORCE)
set(PKG_USER-FEP ON CACHE BOOL "" FORCE)
set(PKG_USER-INTEL ON CACHE BOOL "" FORCE)
set(PKG_KSPACE ON CACHE BOOL "" FORCE)
set(PKG_MC ON CACHE BOOL "" FORCE)
set(PKG_USER-MEAM ON CACHE BOOL "" FORCE)
set(PKG_USER-MOLFILE ON CACHE BOOL "" FORCE)
set(PKG_OPT ON CACHE BOOL "" FORCE)
set(PKG_USER-PHONON ON CACHE BOOL "" FORCE)
set(PKG_QEQ ON CACHE BOOL "" FORCE)
set(PKG_USER-QMMM ON CACHE BOOL "" FORCE)
set(PKG_USER-QTB ON CACHE BOOL "" FORCE)
set(PKG_USER-REAXFF ON CACHE BOOL "" FORCE)
set(PKG_REPLICA ON CACHE BOOL "" FORCE)
set(PKG_USER-SMTBQ ON CACHE BOOL "" FORCE)
set(PKG_SPIN ON CACHE BOOL "" FORCE)
set(PKG_VORONOI ON CACHE BOOL "" FORCE)
```

## Exemple 2 of preset.cmake to compile LAMMPS 2Aug23 parallel

Note that for this later version of LAMMPS some variables have changed names like the `BUILD_SHARED_LIBS` or `PKG_DIFFRACTION` instead of `PKG_USER-DIFFRACTION`.

```
set(BUILD_SHARED_LIBS ON CACHE BOOL "compile lammps as a library" FORCE)
set(BUILD_MPI ON CACHE BOOL "no mpi" FORCE)
set(BUILD_OMP OFF CACHE BOOL "no omp" FORCE)

# Packages
set(PKG_MANYBODY ON CACHE BOOL "" FORCE)
set(PKG_ATC ON CACHE BOOL "" FORCE)
set(PKG_CORESHELL ON CACHE BOOL "" FORCE)
set(PKG_DIFFRACTION ON CACHE BOOL "" FORCE)
set(PKG_FEP ON CACHE BOOL "" FORCE)
set(PKG_INTEL ON CACHE BOOL "" FORCE)
set(PKG_KSPACE ON CACHE BOOL "" FORCE)
set(PKG_MC ON CACHE BOOL "" FORCE)
set(PKG_MEAM ON CACHE BOOL "" FORCE)
set(PKG_MOLFILE ON CACHE BOOL "" FORCE)
set(PKG_OPT ON CACHE BOOL "" FORCE)
set(PKG_PHONON ON CACHE BOOL "" FORCE)
set(PKG_QEQ ON CACHE BOOL "" FORCE)
set(PKG_QMMM ON CACHE BOOL "" FORCE)
set(PKG_QTB ON CACHE BOOL "" FORCE)
set(PKG_REAXFF ON CACHE BOOL "" FORCE)
set(PKG_REPLICA ON CACHE BOOL "" FORCE)
set(PKG_SMTBQ ON CACHE BOOL "" FORCE)
set(PKG_SPIN ON CACHE BOOL "" FORCE)
set(PKG_VORONOI ON CACHE BOOL "" FORCE)
```

# How to compile NDM using LAMMPS as a library
To compile NDM with LAMMPS add `LAMMPS` to the list of NDM packages:
```
cmake .. -D NDM_PACKAGE_LIST=LAMMPS
```
For more detail, please refer to the main README.md file describing how to compile NDM with cmake.
