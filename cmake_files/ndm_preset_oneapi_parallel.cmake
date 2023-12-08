###################################################
#       CMake cache file for NDM
#       preset to compute NDM with oneapi parallel
###################################################

# preset that turns on just oneapi parallel
# this will be compiled quickly and handle a lot of common inputs.
message("PRESET:   use ndm_preset_oneapi_parallel.cmake")

# Define compilers
set(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "$FC user choice fortran compiler f95 gfortran ifort ... mpifort mpiifort ..." FORCE)
set(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "TODO" FORCE) # required for LAMMPS

# NDM compilation mode
set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "NDM for $CMAKE_BUILD_TYPE RELEASE or D" FORCE)
set(NDM_LIBRARY_TYPE "STATIC" CACHE STRING "create libndm as a static/shared library" FORCE)

# MPI
set(NDM_WITH_MPI ON CACHE BOOL "Compilation serial as 'off', parallel as 'on'" FORCE)

# Lammps
set(NDM_WITH_LAMMPS ON CACHE BOOL "use lammps" FORCE)

# MKL configuration, see MKLConfig.cmake file for more info on available parameters
set(NDM_WITH_MKL ON CACHE BOOL "use MKL" FORCE)

#set(ENABLE_SCALAPACK ON CACHE STRING "Scalapack" FORCE)
set(MKL_THREADING "sequential" CACHE STRING "Threading" FORCE)
set(MKL_INTERFACE "lp64" CACHE STRING "Interface" FORCE)
#set(MKL_MPI "openmpi" CACHE STRING "OpenMPI" FORCE)



# -------------- printing infos ----------------------

message("PRESET:   CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
message("PRESET:   CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
message("PRESET:   CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
message("PRESET:   MKL_THREADING=${MKL_THREADING}")
message("PRESET:   MKL_INTERFACE=${MKL_INTERFACE}")
get_cmake_property(_variableNames VARIABLES)
foreach( _variableName ${_variableNames} )
  if ( _variableName MATCHES "NDM_." )
    message("PRESET:   ${_variableName}=${${_variableName}}")
  endif()
endforeach()


