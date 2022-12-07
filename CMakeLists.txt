
CMAKE_minimum_required( VERSION 3.17 )

project( NDM
         VERSION 3.0.0
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


log_info( "cmake in ${PROJECT_SOURCE_DIR}" )
log_info( "CMAKE_MODULE_PATH : ${CMAKE_MODULE_PATH}" )


################################################
# user options
# begin with "NDM_"
################################################

option( NDM_OPT_TRACE "Debug trace mode for cmake" OFF )
option( NDM_OPT_COMPILE_ONLY_CTESTS "Compile ONLY ctests NDM" OFF )
option( NDM_OPT_COMPILE_DOC "Compile doc NDM" OFF )

# from user environment
if( NOT DEFINED ENV{NDM_ROODIR} )
  log_fatal( "Env var $NDM_ROODIR not set, fix it (usually parent folder of NDM)." )
endif()

set( NDM_ROODIR $ENV{NDM_ROODIR} ) # .../NDM_ROODIR
set( NDM_TESDIR $ENV{NDM_TESDIR} ) # .../NDM_ROODIR/ndm_tests
set( NDM_SRCDIR $ENV{NDM_SRCDIR} ) # .../NDM_ROODIR/NDM
set( NDM_BUIDIR $ENV{NDM_BUIDIR} ) # .../NDM_ROODIR/ndm_build...
set( NDM_INSDIR $ENV{NDM_INSDIR} ) # .../NDM_ROODIR/ndm_install...

set( NDM_OAP_INSDIR $ENV{NDM_OAP_INSDIR} ) # .../oneAPI
set( NDM_MPI_INSDIR $ENV{NDM_MPI_INSDIR} ) # ...openmpi


# may be future useful
# set( VT_ROOT $ENV{VT_ROOT} )       # ...itac
# set( CMPLR_ROOT $ENV{CMPLR_ROOT} ) # ...oneapi/compiler/2022.0.1


################################################
# CMAKE usual variables
################################################

# cmake version 3.17+
set( CMAKE_MESSAGE_LOG_LEVEL STATUS CACHE STRING "cmake verbosity  'ERROR|WARNING|NOTICE|STATUS|VERBOSE|DEBUG|TRACE'")
set( CMAKE_VERBOSE_MAKEFILE off CACHE STRING "build verbosity  'off|on'")
log_info ( "CMAKE_VERBOSE_MAKEFILE : ${CMAKE_VERBOSE_MAKEFILE}" )


set( CMAKE_BUILD_TYPE ${NDM_BUILD_TYPE} ) # Release or Debug
# may be debug to fix something
# set( CMAKE_BUILD_TYPE Debug )
# set( CMAKE_BUILD_TYPE Release )


set( CMAKE_INSTALL_PREFIX ${NDM_INSDIR} CACHE STRING
     "Choose the installation directory; by default it installs in the NDM_INSTALL_DIR."
     FORCE )

# check
if ( NOT ${CMAKE_BINARY_DIR} STREQUAL ${NDM_BUIDIR} )
  log_fatal( "$NDM_BUIDIR and $CMAKE_BINARY_DIR (from cmake -B option) not equals, fix it.\n ${NDM_BUIDIR}  !=  ${CMAKE_BINARY_DIR}" )
endif()

#if (NOT (NDM_TYPE STREQUAL "MIX" OR  NDM_TYPE STREQUAL "INTEL" OR NDM_TYPE STREQUAL "GNU"))
if (NOT (NDM_TYPE STREQUAL "INTEL" OR NDM_TYPE STREQUAL "GNU"))
  log_fatal( "NDM_TYPE environ variable unknown, fix it as INTEL or GNU" )
endif()


if ( NDM_WITH_MPI )
  log_info( "NDM parallel as NDM_WITH_MPI ${NDM_WITH_MPI}" )
  set( MPI_HOME ${NDM_MPI_INSDIR} )   # CMAKE usual variable
  find_package(MPI REQUIRED)
  set( NDM_MPI_RESUME "PARALLEL")
else()
  log_info( "NDM serial as NDM_WITH_MPI ${NDM_WITH_MPI}" )
  set( NDM_MPI_RESUME "SERIAL")
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


if ( NDM_WITH_MPI )
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

# make sure that the default is a Release
if (NOT CMAKE_BUILD_TYPE)
  log_warning( "CMAKE_BUILD_TYPE set Release BY DEFAULT" )
  set( CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are None Debug Release."
       FORCE )
endif()

get_filename_component( Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME )
log_info( "CMAKE_Fortran_FLAGS_RELEASE default : ${CMAKE_Fortran_FLAGS_RELEASE}" )
log_info( "CMAKE_Fortran_FLAGS_DEBUG   default : ${CMAKE_Fortran_FLAGS_DEBUG}" )

log_info( "CMAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}" )
log_info( "CMAKE_INSTALL_PREFIX : ${CMAKE_INSTALL_PREFIX}" )
log_info( "CMAKE_Fortran_COMPILER full path : ${CMAKE_Fortran_COMPILER}" )
log_info( "Fortran_COMPILER_NAME : ${Fortran_COMPILER_NAME}" )


set(NDM_RESUME "(${NDM_TYPE}) NDM ${NDM_MPI_RESUME}, compilation ${CMAKE_BUILD_TYPE}")
log_title( ${NDM_RESUME} )


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

# may be #ifdef or #ifndef or #if  TODO have to be homogenous with LAMMPS or MILADY
if ( NDM_WITH_MPI STREQUAL "on")
  #add_compile_definitions(ML=0)
  add_compile_definitions(PARAML=1)
  add_compile_definitions(PARA=1)
  log_warning( "add_compile_definitions : PARALLEL" )
else()
  #add_compile_definitions(ML=1)
  log_warning( "add_compile_definitions : SERIAL" )
endif()
log_info( "CMAKE_Fortran_FLAGS_RELEASE current : ${CMAKE_Fortran_FLAGS_RELEASE}" )
log_info( "CMAKE_Fortran_FLAGS_DEBUG   current : ${CMAKE_Fortran_FLAGS_DEBUG}" )

# setenv.bash to set environ to use code ouside cmake/make
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

if( NDM_OPT_COMPILE_DOC )
  find_program( NDM_JUPYTER_BOOK NAMES jupyter-book DOC "jupyter-book build executable")
  if ( NDM_JUPYTER_BOOK )
    log_warning( "add_subdirectory docs" )
    add_subdirectory( docs )
  else()
    log_warning( "prerequisite jupyter-book not found, not building docs" )
  endif()
else()
  log_info( "cmake option NDM_OPT_COMPILE_DOC ${NDM_OPT_COMPILE_DOC}, not building docs" )
endif()

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
