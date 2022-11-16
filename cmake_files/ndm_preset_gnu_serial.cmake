# preset that turns on just gnu serial
# this will be compiled quickly and handle a lot of common inputs.

message(STATUS "!!!! preset on ndm_gnu_serial.cmake")

set(NDM_TYPE "GNU" CACHE STRING "GNU or MIX or INTEL" FORCE)
set(NDM_BUILD_TYPE "Release" CACHE STRING "NDM for $CMAKE_BUILD_TYPE Release or Debug" FORCE)
set(NDM_WITH_MPI "off" CACHE BOOL "Compilation serial as 'off', parallel as 'on'" FORCE)
set(NDM_FC_USER "f95" CACHE STRING "NDM for $FC user choice fortran compiler f95 gfortran ..." FORCE)
set(NDM_CMAKE_VERBOSE "off" CACHE STRING "NDM for cmake more verbosity  'off|on'" FORCE)

message(STATUS "!!!! NDM_TYPE : ${NDM_TYPE}")
message(STATUS "!!!! NDM_FC_USER : ${NDM_FC_USER}")
message(STATUS "!!!! NDM_BUILD_TYPE : ${NDM_BUILD_TYPE}")
