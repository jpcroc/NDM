# launch cmake
# https://cmake.org/cmake/help/latest/manual/cmake.1.html#generate-a-project-buildsystem

# preset that turns on just gnu serial
# this will be compiled quickly and handle a lot of common inputs.

message("PRESET:   use ndm_preset_gnu_serial.cmake")

# CMAKE usual preset variables
# https://cmake.org/cmake/help/latest/envvar/FC.html
set(CMAKE_Fortran_COMPILER "f95" CACHE STRING "$FC user choice fortran compiler f95 gfortran ifort ... mpifort ..." FORCE)

# NDM usual preset variables
set(NDM_TYPE "GNU" CACHE STRING "GNU or MIX or INTEL" FORCE)
set(NDM_BUILD_TYPE "Release" CACHE STRING "NDM for $CMAKE_BUILD_TYPE Release or Debug" FORCE)
set(NDM_WITH_MPI "off" CACHE BOOL "Compilation serial as 'off', parallel as 'on'" FORCE)
set(NDM_MPI_INSDIR "useless" CACHE STRING "NDM with MPI /usr/lib64/openmpi" FORCE)


message("PRESET:   CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
get_cmake_property(_variableNames VARIABLES)
foreach(_variableName ${_variableNames})
  if ( _variableName MATCHES "NDM_." )
    message("PRESET:   ${_variableName}=${${_variableName}}")
  endif()
endforeach()
