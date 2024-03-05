###################################################
#       CMake cache file for NDM
#       preset to compute NDM with GNU serial
###################################################
set(NDM_PRESET_TITLE "Compile NDM GNU serial : gfortran" CACHE STRING "" FORCE)

# ------------------- Define compilers ----------------------
set(CMAKE_Fortran_COMPILER "gfortran" CACHE STRING "$FC user choice fortran compiler f95 gfortran ifort ... mpifort mpiifort ..." FORCE)
#set(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "TODO" FORCE) # required for LAMMPS

# ---------------- preprocessor definitions ------------------
list(APPEND TMP_COMPILE_DEFINITION "MKL")
set(NDM_COMPILE_DEFINITION ${TMP_COMPILE_DEFINITION} CACHE INTERNAL "List of definitions for preprocessor" FORCE)

# ------ Choose compilation mode and corresponding flags ------
set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "RELEASE or DEBUG" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE " -O3 -ffree-line-length-512 -funroll-all-loops" CACHE STRING "" FORCE) 
set(CMAKE_Fortran_FLAGS_DEBUG " -O0 -g -C -fpe-all=0 -ffree-line-length-512" CACHE STRING "" FORCE)

# --------------------- cmake options ------------------------
set( NDM_OPT_TRACE OFF CACHE BOOL "trace all variables for cmake debug" FORCE )
set( NDM_OPT_COMPILE_DOC OFF CACHE BOOL "Compile documentation" FORCE )

# ------------------- MKL configuration ----------------------
set(ENABLE_SCALAPACK OFF CACHE STRING "Scalapack" FORCE)
set(MKL_THREADING "sequential" CACHE STRING "Threading" FORCE)
set(MKL_INTERFACE "lp64" CACHE STRING "Interface" FORCE)

# ---------------- Choose packages to use -------------------
list(APPEND TMP_PACKAGE_LIST "MKL")
#list(APPEND TMP_PACKAGE_LIST "LAMMPS")
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
