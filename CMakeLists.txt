##############################################
#   Project configuration
##############################################

CMAKE_minimum_required( VERSION 3.28 )

project( NDM
         VERSION 3.0.0
         DESCRIPTION "NDM project with CMake"
         HOMEPAGE_URL "https://codev-tuleap.intra.cea.fr/plugins/git/ndm"
         LANGUAGES Fortran CXX)

set( PROJECT_AUTHOR "Jean-Paul Crocombette" )

################################################
# utilities
################################################
# get NDM cmake macro etc
get_property(${PROJECT_NAME}_LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)   # Get the languages list for this project


list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_files" )
include( TraceNewMacros )   # useful for colorful logs.
include( import_lammps )    # macro for importing lammps
setcolors()                 # defines log colors

include( CTest )

log_info( "cmake in ${PROJECT_SOURCE_DIR}" )                # log
log_info( "CMAKE_MODULE_PATH : ${CMAKE_MODULE_PATH}" )      # log

################################################
# set user options and cache variables
################################################

option(${PROJECT_NAME}_WITH_MPI "Compile ${PROJECT_NAME} in Parallel using MPI" ON)
option(${PROJECT_NAME}_WITH_LAMMPS "Compile ${PROJECT_NAME} using lammps as a library" ON)

# MKL configuration
set(ENABLE_SCALAPACK OFF CACHE BOOL "using scalapack")
set(MKL_THREADING "" CACHE STRING "description threading")
set(MKL_INTERFACE "" CACHE STRING "description interface")
set(MKL_MPI "" CACHE STRING "description MKL MPI")


################################################
# CMAKE usual variables
################################################

# cmake verbosity
set( CMAKE_MESSAGE_LOG_LEVEL STATUS CACHE STRING "cmake verbosity  'ERROR|WARNING|NOTICE|STATUS|VERBOSE|DEBUG|TRACE'")
set( CMAKE_VERBOSE_MAKEFILE off CACHE STRING "build verbosity  'off|on'")
log_info ( "CMAKE_VERBOSE_MAKEFILE : ${CMAKE_VERBOSE_MAKEFILE}" )

###################################################
#   Loading packages
###################################################

# MPI
if ( ${PROJECT_NAME}_WITH_MPI )
    log_info( "${PROJECT_NAME} parallel as ${PROJECT_NAME}_WITH_MPI ${${PROJECT_NAME}_WITH_MPI}" )
    add_compile_definitions(PARA=1)         # Preprocessor definition
    add_compile_definitions(PARAML=1)       # Preprocessor definition
    find_package(MPI REQUIRED)
    foreach (LANGUAGE IN ITEMS ${${PROJECT_NAME}_LANGUAGES})
        link_libraries(${MPI_${LANGUAGE}_LIBRARIES}) # add MPI libraries to all later-defined targets
    endforeach()
    set( ${PROJECT_NAME}_MPI_RESUME "PARALLEL")
else()
    log_info( "${PROJECT_NAME} serial as ${PROJECT_NAME}_WITH_MPI ${${PROJECT_NAME}_WITH_MPI}" )
    set( ${PROJECT_NAME}_MPI_RESUME "SERIAL")
endif()

if ( ${PROJECT_NAME}_WITH_MKL )
    find_package(MKL CONFIG REQUIRED)
    link_libraries($<LINK_ONLY:MKL::MKL>)                       # add MKL libraries to all later-defined targets
endif()
if (${PROJECT_NAME}_WITH_LAMMPS)
    log_info( "Compile ${PROJECT_NAME} with Lammps as a library" )
    import_lammps()
    add_compile_definitions(LAMMPS_VERSION=1)   # Preprocessor definition
endif()


################################################
# Other options
# as CMAKE usage needs
################################################

# make sure that the default is a Release
if (NOT CMAKE_BUILD_TYPE)
  log_warning( "CMAKE_BUILD_TYPE set Release BY DEFAULT" )
  set( CMAKE_BUILD_TYPE "RELEASE" CACHE STRING
       "Choose the type of build, options are None Debug Release."
       FORCE )
endif()

log_info("Build type : ${CMAKE_BUILD_TYPE}")
foreach(LANGUAGE IN ITEMS ${${PROJECT_NAME}_LANGUAGES})
    log_info("${LANGUAGE} compiler : ${CMAKE_${LANGUAGE}_COMPILER}")
    log_info("${LANGUAGE} compiler flags : ${CMAKE_${LANGUAGE}_FLAGS} ${CMAKE_${LANGUAGE}_FLAGS_${CMAKE_BUILD_TYPE}}")
endforeach()

set(${PROJECT_NAME}_RESUME "(${${PROJECT_NAME}_TYPE}) ${PROJECT_NAME} ${${PROJECT_NAME}_MPI_RESUME}, compilation ${CMAKE_BUILD_TYPE}")
log_title( ${${PROJECT_NAME}_RESUME} )


################################################
# set paths mod executables
################################################
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod )
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin )

if( NOT ${PROJECT_NAME}_OPT_COMPILE_ONLY_CTESTS )
  log_warning( "add_subdirectory src" )
  add_subdirectory( src )
endif()

log_warning( "add_subdirectory tests" )
add_subdirectory( tests )

#if( NDM_OPT_COMPILE_DOC )
#  find_program( NDM_JUPYTER_BOOK NAMES jupyter-book DOC "jupyter-book build executable")
#  if ( NDM_JUPYTER_BOOK )
#    log_warning( "add_subdirectory docs" )
#    add_subdirectory( docs )
#  else()
#    log_warning( "prerequisite jupyter-book not found, not building docs" )
#  endif()
#else()
#  log_info( "cmake option NDM_OPT_COMPILE_DOC ${NDM_OPT_COMPILE_DOC}, not building docs" )
#endif()

#add_custom_target(copy_all ALL)

#list( APPEND NDM_COPY_DIRS "scripts;READMES;images" )
# TODO later list( APPEND NDM_COPY_DIRS "ndmpy" )

#foreach( adir ${NDM_COPY_DIRS} )
#  # to build dir
#  add_custom_command(TARGET copy_all POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E copy_directory
#                ${CMAKE_CURRENT_SOURCE_DIR}/${adir}
#                ${CMAKE_CURRENT_BINARY_DIR}/${adir} )
#  # to install dir
#  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${adir} DESTINATION ${CMAKE_INSTALL_PREFIX}/${adir})
#endforeach()

# resume TRACE log
# "NDM_." is regExp
#printCmakeTraceVariable( "NDM_." )
# printCmakeTraceVariable( "CMAKE_INSTALL." )
# optionally get trace of all cmake main variables in log
#if( NDM_OPT_TRACE )
#printCmakeTrace()
#endif()

# test_colormessage()
