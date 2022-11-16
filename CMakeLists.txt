################################################################################################
# thanks to
# https://stackoverflow.com/questions/52712416/cmake-fortran-module-directory-to-be-used-with-add-library
# https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/languages/fortran/ForFortranExample
# https://gitlab.com/CLIUtils/modern-cmake/-/blob/master/examples/extended-project/CMakeLists.txt
# https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html
#
# export NDM_TYPE=MIX # or INTEL, GNU
# cmake3 -S<path-to-source> -B<path-to-build>
#
# all NDM_xxx vars are for EZ cmake resume debug with something as
# cmake3 -S${SRCDIR} -B${BUIDIR} -DCMAKE_BUILD_TYPE=Release -DNDM_OPT_TRACE=ON | grep "TRACE " | grep "NDM_"
################################################################################################

CMAKE_minimum_required( VERSION 3.17 )

project( NDM
         VERSION 2.0.0
         DESCRIPTION "NDM project with CMake"
         HOMEPAGE_URL "https://codev-tuleap.intra.cea.fr/plugins/git/ndm"
         LANGUAGES Fortran C )

set( PROJECT_AUTHOR "Jean-Paul Crocombette" )

# get NDM cmake macro etc
list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_files" )


################################################
# utilities
################################################

include( TraceNewMacros )
include( CTest )

setcolors()
log_info( "CMAKE_MODULE_PATH : ${CMAKE_MODULE_PATH}" )

################################################
# user options
# begin with "NDM_"
################################################

option( NDM_OPT_TRACE "Debug trace mode for cmake" OFF )
option( NDM_OPT_COMPILE_ONLY_CTESTS "Compile ONLY ctests NDM" OFF )
option( NDM_OPT_COMPILE_DOC "Compile doc NDM" OFF )

# from environment
set( NDM_ROODIR $ENV{NDM_ROODIR} ) # .../NDM_ROODIR
set( NDM_TESDIR $ENV{NDM_TESDIR} ) # .../NDM_ROODIR/ndm_tests
set( NDM_SRCDIR $ENV{NDM_SRCDIR} ) # .../NDM_ROODIR/NDM
set( NDM_BUIDIR $ENV{NDM_BUIDIR} ) # .../NDM_ROODIR/ndm_build...
set( NDM_INSDIR $ENV{NDM_INSDIR} ) # .../NDM_ROODIR/ndm_install...

set( NDM_OAP_INSDIR $ENV{NDM_OAP_INSDIR} ) # .../oneAPI
set( NDM_MPI_INSDIR $ENV{NDM_MPI_INSDIR} ) # ...openmpi


#set( VT_ROOT $ENV{VT_ROOT} ) # ...itac
#set( CMPLR_ROOT $ENV{CMPLR_ROOT} ) # ...oneapi/compiler/


# CMAKE usual variables
set( CMAKE_BUILD_TYPE ${NDM_BUILD_TYPE} ) # Release or Debug
set( CMAKE_COMPILE_MODE ${NDM_BUILD_TYPE}$ )
set( FC ${NDM_FC_USER} )
set( Fortran_COMPILER_NAME ${FC} )

set( MPI_HOME ${NDM_MPI_INSDIR} )
set( CMAKE_INSTALL_PREFIX ${NDM_INSDIR} )

# check
if ( NOT ${CMAKE_BINARY_DIR} STREQUAL ${NDM_BUIDIR} )
  log_fatal( "$NDM_BUIDIR and $CMAKE_BINARY_DIR (from cmake -B option) not equals, fix it.\n ${NDM_BUIDIR}  !=  ${CMAKE_BINARY_DIR}" )
endif()

if (NOT (NDM_TYPE STREQUAL "MIX" OR  NDM_TYPE STREQUAL "INTEL" OR NDM_TYPE STREQUAL "GNU"))
  log_fatal( "NDM_TYPE environ variable unknown, fix it as MIX or INTEL or GNU" )
endif()


if ( NDM_WITH_MPI )
  log_info( "NDM parallel as NDM_WITH_MPI ${NDM_WITH_MPI}" )
  find_package(MPI REQUIRED)
else()
  log_info( "NDM serial as NDM_WITH_MPI ${NDM_WITH_MPI}" )
endif()

# find_package(MKL REQUIRED) # because not found on irene...
##find_package(MKL)
##include_directories( SYSTEM ${MKL_INCLUDE} )


# If you want yo use Intel MPI with Intel Fortran compiler, the wrapper is mpiifort.
# If you want to use an other MPI (such as Open MPI) with use mpi
# then you likely need to rebuild the MPI library because Fortran modules
# are not interoperable between vendors,
# and sometimes not even different compiler versions of the same vendor.
# I strongly discourage you to use ifort and manually set the include and library path.


if( NDM_OPT_COMPILE_DOC )
  find_package( Doxygen )  # TODO Sphinx instead
  if (Doxygen_FOUND)
    # TODO to make doc ?
    log_warning( "add_subdirectory docs" )
    add_subdirectory( docs )
  else()
    log_info( "Doxygen not found, not building docs" )
  endif()
endif()

# before all cmake default installation
# get_filename_component( default_prefix ".." ABSOLUTE )
set( default_prefix ${NDM_ROODIR} )
log_info( "cmake default_prefix : ${default_prefix}" )

# where make install do
set( NDM_DEFAULT_PREFIX ${default_prefix} )

################################################
# trace user options
# begin with "NDM_"
################################################

log_info( "cmake with trace : ${NDM_OPT_TRACE}" )
log_info( "compile ONLY ctests NDM : ${NDM_OPT_COMPILE_ONLY_CTESTS}" )

# may be debug to fix it
# set( CMAKE_BUILD_TYPE Debug )
# set( CMAKE_BUILD_TYPE Release )

# TODO verify if not defined if(DEFINED ENV{<name>})
# set( MKLROOT $ENV{MKLROOT} ) # MKLROOT and MKL_ROOT are equals in fine
set( MKL_ROOT $ENV{MKL_ROOT} )
log_info( "CMAKE MKL_ROOT : ${MKL_ROOT}" )

if ( NDM_WITH_MPI )
  set( NDM_MPI_ROODIR $ENV{NDM_MPI_ROODIR} )
  log_info( "NDM_MPI_ROODIR : ${NDM_MPI_ROODIR}" )
  log_info( "MPI_HOME :${MPI_HOME}" )
  log_info( "MPI_DETERMINE_Fortran_CAPABILITIES : ${MPI_DETERMINE_Fortran_CAPABILITIES}" )
  log_info( "MPI_Fortran_MODULE_DIR : ${MPI_Fortran_MODULE_DIR}" )
  log_info( "MPI_Fortran_HAVE_F08_MODULE : ${MPI_Fortran_HAVE_F08_MODULE}" )
  log_info( "MPI_Fortran_HAVE_F90_MODULE: ${MPI_Fortran_HAVE_F90_MODULE}" )
  log_info( "MPI_Fortran_VERSION : ${MPI_Fortran_VERSION}" )
endif()

################################################
# Other options
# as CMAKE usage needs
################################################

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  log_warning( "CMAKE_BUILD_TYPE set Release BY DEFAULT" )
  set( CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE )
endif()

set( CMAKE_INSTALL_PREFIX ${NDM_INSDIR} CACHE STRING
     "Choose the installation directory; by default it installs in the NDM_INSTALL_DIR."
     FORCE )

# FFLAGS depend on the compiler
get_filename_component( Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME )


log_info( "CMAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}" )
log_info( "CMAKE_INSTALL_PREFIX : ${CMAKE_INSTALL_PREFIX}" )
log_info( "CMAKE_Fortran_COMPILER full path : ${CMAKE_Fortran_COMPILER}" )
log_info( "Fortran_COMPILER_NAME : ${Fortran_COMPILER_NAME}" )
log_info( "MKL_ROOT : ${MKL_ROOT}" )


if (DEFINED ENV{NDM_MKL_LIB})
  list( APPEND NDM_MKL_LIB $ENV{NDM_MKL_LIB} )
  log_info( "NDM_MKL_LIB as reached environ variable :\n         ${NDM_MKL_LIB}" )
else()
  log_warning( "NDM_MKL_LIB environ variable unknown, we just try BY DEFAULT" )
  # set it for 'gf' as 'mpifort-gnu' and 'intel' as 'mpiifort-intel'
  if ( NDM_TYPE STREQUAL "MIX")
    list( APPEND NDM_MKL_LIB mkl_scalapack_lp64;mkl_intel_lp64;mkl_sequential;mkl_core;mkl_blacs_openmpi_lp64;pthread;iomp5;m;dl;mpi;mpi_mpifh)
  endif()
  if ( NDM_TYPE STREQUAL "INTEL")
    list( APPEND NDM_MKL_LIB  mkl_scalapack_lp64;mkl_intel_lp64;mkl_sequential;mkl_core;mkl_blacs_intelmpi_lp64;pthread;m;dl)
  endif()
  if ( NDM_TYPE STREQUAL "GNU")
    list( APPEND NDM_MKL_LIB mkl_scalapack_lp64;mkl_gf_lp64;mkl_sequential;mkl_core;mkl_blacs_openmpi_lp64;pthread;m;dl)
  endif()
  log_warning( "NDM_MKL_LIB as DEFAULT :\n         ${NDM_MKL_LIB}")
endif()

set( NDM_MKL_LIB "" )
log_warning( "TODO NDM_MKL_LIB reset as debug")

set(NDM_RESUME "(${NDM_TYPE}) NDM_WITH_MPI ${NDM_WITH_MPI} compilation ${CMAKE_BUILD_TYPE}")
log_title( ${NDM_RESUME} )

if (NDM_TYPE STREQUAL "INTEL")
  if (NOT (Fortran_COMPILER_NAME STREQUAL "mpiifort" OR Fortran_COMPILER_NAME STREQUAL "ifort"))
    log_fatal( "NDM_TYPE INTEL needs Fortran_COMPILER_NAME FC=mpiifort or ifort, fix it" )
  endif()

  # mpiifort as all oneAPI compiled products
  add_compile_definitions(IFORT)

  # MAGICAL DETECT OVERFLOW ... TODO choose previous for production, -C is check bounds
  set( TTMP "${MKL_ROOT}/lib/intel64" )  # /volatile2/intel/oneapi/mkl/2022.0.1/lib/intel64
  set( CMAKE_Fortran_FLAGS_RELEASE " -O3  -traceback -L${TMP_MKL} -L${TMP_MPI} -L${TTMP}" )


  set( CMAKE_Fortran_FLAGS_DEBUG " -O0 -debug all -g -fpe0 -traceback -C " )
  # could be from ${MKL_SCALAPACK} ${MKL_LIBRARIES} ${MKL_BLACS} :
  # gives mkl_scalapack_ilp64 mkl_intel_ilp64 mkl_sequential mkl_core mkl_intel_thread mkl_blacs_intelmpi_ilp64 m dl
  # The Intel MKL ILP64 libraries use the 64-bit integer type
  # (necessary for indexing large arrays, with more than 231-1 elements),
  # whereas the LP64 libraries index arrays with the 32-bit integer type.
  # at usage ILP64 have problems as ctest... 200224
  # http://cali2.unilim.fr/intel-xe/mkl/mkl_userguide/GUID-87821148-338B-4022-8C90-F24C866F2878.htm
  # may be Compiling for ILP64 'ifort -i8 -I<mkl directory>/include ...' 220224



elseif (NDM_TYPE STREQUAL "MIX")
  if (NOT (Fortran_COMPILER_NAME STREQUAL "mpiifort" OR Fortran_COMPILER_NAME STREQUAL ifort))
    log_fatal( "NDM_TYPE MIX needs Fortran_COMPILER_NAME FC=mpiiifort or ifort, fix it" )
  endif()

  # ifort oneAPI as cosmin preferred config 220203 with openmpi ifort compiled
  add_compile_definitions(IFORT)

  set( TMP_MKL "${MKL_ROOT}/lib/intel64" )  # ...intel/oneapi/mkl/2022.0.1/lib/intel64
  set( TMP_MPI ${MPI_Fortran_MODULE_DIR})  # .../mpi/install_412/lib

  # MAGICAL DETECT OVERFLOW ... # -C is check bounds
  #set( CMAKE_Fortran_FLAGS_RELEASE "  -C -fpe-all=0 -no-ftz -traceback -L${TMP_MKL} -L${TMP_MPI}")
  # set( CMAKE_Fortran_FLAGS_RELEASE " -C -fpe0 -traceback -L${TMP_MKL} -L${TMP_MPI}" )

  set( CMAKE_Fortran_FLAGS_RELEASE " -O3 -traceback -L${TMP_MKL} -L${TMP_MPI} " )
  set( CMAKE_Fortran_FLAGS_DEBUG " -O0 -debug all -g -fpe0 -traceback -C -L${TMP_MKL} -L${TMP_MPI} " )

  # .../mpi/...412/lib .../mpi/...316/lib )
  include_directories( SYSTEM  ${MPI_Fortran_MODULE_DIR} )  # get *.mod ifort
  # LD_LIBRARY_PATH contains .../mpi/install_316/lib

  add_compile_definitions(MPI_STATUSES_IGNORE=1)



elseif (NDM_TYPE STREQUAL "GNU")
  if (NOT (Fortran_COMPILER_NAME STREQUAL "mpifort" OR Fortran_COMPILER_NAME STREQUAL "f95" OR Fortran_COMPILER_NAME STREQUAL "gfortran"))
    log_fatal( "NDM_TYPE GNU needs Fortran_COMPILER_NAME FC=mpifort or f95 or gfortran, fix it: ${Fortran_COMPILER_NAME}" )
  endif()

  # gfortran gnu 10 mini 220203
  add_compile_definitions(GNU)

  set( TMP_MKL "${MKL_ROOT}/lib/intel64" )  # /volatile2/intel/oneapi/mkl/2022.0.1/lib/intel64
  set( TMP_MPI ${MPI_Fortran_MODULE_DIR})

  # set( TMP_MPI ${MPI_Fortran_MODULE_DIR})   # .../mpi/install_316/lib
  # from module load gcc/9.3.0 & module load openmpi/gcc_9.3.0/4.0.1
  # set( TMP_MPI /product/centos7-x86_64/apps/OpenMPI/4.0.1-GCC-9.3.0/lib ) # from module load gcc/9.3.0

  # set (CMAKE_Fortran_FLAGS_RELEASE " -funroll-all-loops -fno-f2c -O3 -ffree-form -ffree-line-length-512 -L${TMP_MKL} -L${TMP_MPI}")

  set (CMAKE_Fortran_FLAGS_RELEASE " -funroll-all-loops -O3 -ffree-line-length-512 ")
  set (CMAKE_Fortran_FLAGS_DEBUG   " - fno-f2c -O0 -g -ffree-form  -ffree-line-length-512 ")

  include_directories(SYSTEM  ${MPI_Fortran_MODULE_DIR} )  # get *.mod ifort


else()
  # others general else are fatal houston problem
  log_fatal( "No optimized Fortran compiler flags are known, we just try NOTHING!" )
endif()


log_info( "NDM_MKL_LIB : ${NDM_MKL_LIB}")


################################################
# preprocessor definitions used in NDM 2022
# cf in code '#if(PHONDY || PARAPH || MAB || ML || PARAML)'
# au cas ou, pour une cible specifique,
# vous pouvez effectuer les operations suivantes :
# target_compile_definitions(my_target PRIVATE FOO=1 BAR=1)
################################################

# useless from cosmin 220210
# add_compile_definitions(PHONDY="")
# add_compile_definitions(PARAPH="")
# add_compile_definitions(MAB="")
# add_compile_definitions(PARA=1)
# add_compile_definitions(PARAM=1)

# TODO if NDM_WITH_MPI
#add_compile_definitions(ML=1)
#add_compile_definitions(PARAML=0)
#add_compile_definitions(PARA=)

log_info( "CMAKE_Fortran_FLAGS_RELEASE :\n       ${CMAKE_Fortran_FLAGS_RELEASE}" )
log_info( "CMAKE_Fortran_FLAGS_DEBUG :\n       ${CMAKE_Fortran_FLAGS_DEBUG}" )

# TODO if windows
configure_file(ndm_setenv.bash.in ndm_setenv.bash @ONLY NEWLINE_STYLE UNIX )

################################################
# set paths mod executables
################################################
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod )
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin )


if( NOT NDM_OPT_COMPILE_ONLY_CTESTS )
  log_warning( "add_subdirectory src" )
  add_subdirectory( src )
endif()

log_warning( "add_subdirectory tests" )
add_subdirectory( tests )

add_custom_target(copy_all ALL)

list( APPEND NDM_COPY_DIRS "scripts;READMES;images" )
# TODO later list( APPEND NDM_COPY_DIRS "ndmpy" )

foreach( adir ${NDM_COPY_DIRS} )
  # to build dir
  add_custom_command(TARGET copy_all POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_CURRENT_SOURCE_DIR}/${adir}
                ${CMAKE_CURRENT_BINARY_DIR}/${adir} )
  # to install dir
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${adir} DESTINATION ${CMAKE_INSTALL_PREFIX}/${adir})
endforeach()



# resume TRACE log
# "NDM_." is regExp
printCmakeTraceVariable( "NDM_." )
# printCmakeTraceVariable( "CMAKE_INSTALL." )


# optionally get trace of all cmake main variables in log
if( NDM_OPT_TRACE )
  printCmakeTrace()
endif()

# test_colormessage()
