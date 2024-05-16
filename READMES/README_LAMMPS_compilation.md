# How to compile LAMMPS as a library to use with NDM

On Gatsby on version is installed already.

## Custom lammps


LAMMPS compilation can be configured with command line options but it is often more convenient to use cmake preset files.
See below examples of preset files to compile LAMMPS as a library. LAMMPS normally provides more examples in the cmake directory.

```
cmake ../cmake/. -C <preset.cmake>
make -j
make install
```






## Exemple 1 of preset.cmake to compile LAMMPS 3Mar20 serial
```
set(LAMMPS_MACHINE "serial" CACHE STRING "serial compilation" FORCE)
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

## Exemple 2 of preset.cmake to compile LAMMPS 3Mar20 serial


# How to compile NDM using LAMMPS as a library
Turn the option on and pray your preferred divinity.
