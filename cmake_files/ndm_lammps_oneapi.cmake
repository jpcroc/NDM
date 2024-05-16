###################################################
#       CMake cache file for NDM
#       preset to compile NDM with LAMMPS and oneapi parallel
###################################################
set(NDM_PRESET_TITLE "Compile NDM parallel with LAMMPS and intel" CACHE STRING "" FORCE)

# ------------------- Define compilers ----------------------
set(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "TODO" FORCE) # required for LAMMPS

# ------------------- LAMMPS ----------------------
#set(LAMMPS_HOME "/volatile/catB/jd270899/local_softwares/lammps/lammps-3Mar20/build_library" CACHE STRING "Custom LAMMPS location" FORCE)

# ---------------- Choose packages to use -------------------
list(APPEND TMP_PACKAGE_LIST "LAMMPS")
set(NDM_PACKAGE_LIST ${TMP_PACKAGE_LIST} CACHE INTERNAL "List of packages" FORCE)


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
