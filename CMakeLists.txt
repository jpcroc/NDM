CMAKE_minimum_required( VERSION 3.20 )
################################################
# Auto configuration
# in case no configuration file is used
################################################
# Cache variables will not be reset if already defined with cmake .. -C <cache file>

find_program(VAR_MPIIFORT mpiifort) # check if mpiifort is available for intel parallel
find_program(VAR_IFORT ifort)       # else check if ifort is available for intel serial
find_program(VAR_MPIFORT mpifort)   # else check if mpifort is available for GNU parallel
find_program(VAR_GFORTRAN gfortran) # else check if gfortran is available for GNU serial
set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "RELEASE or DEBUG")

if (NOT ${VAR_MPIIFORT} STREQUAL "VAR_MPIIFORT-NOTFOUND") # if mpiifort not not found, then INTEL parallel
    set(NDM_PRESET_TITLE "Compile NDM intel parallel : mpiifort" CACHE STRING "")
    set(CMAKE_Fortran_COMPILER "${VAR_MPIIFORT}" CACHE STRING "")    
    set(NDM_PARA ON)
elseif (NOT ${VAR_IFORT} STREQUAL "VAR_IFORT-NOTFOUND") # if mpiifort not not found, then INTEL serial
    set(NDM_PRESET_TITLE "Compile NDM intel serial : ifort" CACHE STRING "")
    set(CMAKE_Fortran_COMPILER "${VAR_IFORT}" CACHE STRING "")
elseif (NOT ${VAR_MPIFORT} STREQUAL "VAR_MPIFORT-NOTFOUND") # if mpifort not not found, then GNU parallel
    set(NDM_PRESET_TITLE "Compile NDM GNU parallel : gfortran" CACHE STRING "")
    set(CMAKE_Fortran_COMPILER "${VAR_GFORTRAN}" CACHE STRING "")
    set(NDM_GNU ON)
    set(NDM_PARA ON)
elseif (NOT ${VAR_GFORTRAN} STREQUAL "VAR_GFORTRAN-NOTFOUND") # if mpifort not not found, then GNU serial
    set(NDM_PRESET_TITLE "Compile NDM GNU serial : gfortran" CACHE STRING "")
    set(CMAKE_Fortran_COMPILER "${VAR_GFORTRAN}" CACHE STRING "")
    set(NDM_GNU ON)
endif()

if (NDM_GNU)
    set(CMAKE_Fortran_FLAGS_RELEASE " -O3 -ffree-line-length-512 -funroll-all-loops" CACHE STRING "")
    set(CMAKE_Fortran_FLAGS_DEBUG " -O0 -g -C -fpe-all=0 -ffree-line-length-512" CACHE STRING "")
    set(MKL_MPI "openmpi" CACHE STRING "MPI")
endif()

if (NDM_PARA)
    list(APPEND TMP_PACKAGE_LIST "MPI")
    list(APPEND TMP_COMPILE_DEFINITION "PARA")
    list(APPEND TMP_COMPILE_DEFINITION "PARAML")
endif()

set(CMAKE_Fortran_FLAGS_RELEASE " -O3 " CACHE STRING "")
set(CMAKE_Fortran_FLAGS_DEBUG " -O0 -g -C -fpe-all=0" CACHE STRING "")

list(APPEND TMP_PACKAGE_LIST "MKL")
set(NDM_PACKAGE_LIST ${TMP_PACKAGE_LIST} CACHE STRING "List of packages")

set(MKL_MPI "intelmpi" CACHE STRING "MPI")
set(ENABLE_SCALAPACK OFF CACHE STRING "Scalapack" )
set(MKL_THREADING "sequential" CACHE STRING "Threading")
set(MKL_INTERFACE "lp64" CACHE STRING "Interface")
set(MKL_ARCH "intel64" CACHE STRING "Architecture")

list(APPEND TMP_COMPILE_DEFINITION "MKL")
set(NDM_COMPILE_DEFINITION ${TMP_COMPILE_DEFINITION} CACHE STRING "List of definitions for preprocessor")

################################################
# Project definition
################################################
project( NDM
         VERSION 3.0.0
         DESCRIPTION "NDM project with CMake"
         HOMEPAGE_URL "https://codev-tuleap.intra.cea.fr/plugins/git/ndm"
         LANGUAGES Fortran ) #CXX)

set( PROJECT_AUTHOR "Jean-Paul Crocombette" )

################################################
# default user options
################################################

option( ${PROJECT_NAME}_OPT_TRACE "trace all variables for cmake debug" OFF )
option( ${PROJECT_NAME}_OPT_COMPILE_TESTS "Allows test compilation" OFF )
option( ${PROJECT_NAME}_OPT_COMPILE_DOC "Compile documentation" OFF )


################################################
# utilities
################################################

# Get the languages list for this project
get_property(${PROJECT_NAME}_LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)   

# cmake scripts location
list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_files" )        
message(STATUS "CMAKE_MODULE_PATH : ${CMAKE_MODULE_PATH}" )

# Import cmake scripts
include(doxygen_doc)        # to compile documentation with Doxygen
include( CTest )            # testing with CTest
include( TraceNewMacros )   # colorful logs
setcolors()                 # setting log colors


################################################
# Packages, ex: MPI, MKL, etc...
################################################
log_title(${${PROJECT_NAME}_PRESET_TITLE})

foreach( PACKAGE IN ITEMS ${${PROJECT_NAME}_PACKAGE_LIST})
    log_info("Compile ${PROJECT_NAME} with ${PACKAGE} ")
    include(compile_with_${PACKAGE})
    cmake_language(CALL compile_with_${PACKAGE})        #cmake 3.18
endforeach()


################################################
# preprocessor definitions 
# cf in code '#if(PHONDY || PARAPH || MAB || ML || PARAML)'
# au cas ou, pour une cible specifique,
# vous pouvez effectuer les operations suivantes :
# target_compile_definitions(my_target PRIVATE FOO=1 BAR=1)
################################################

foreach( VARIABLE IN ITEMS ${${PROJECT_NAME}_COMPILE_DEFINITION})
    add_compile_definitions(${VARIABLE}=1)
endforeach()


################################################
# set paths to modules and executables
################################################

# setting locations in the build directory
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )

################################################
#   set targets and sources
################################################

# Define main library to add source files to in subdirectories
set(MAIN_LIBRARY_NAME "ndm")
set(MAIN_EXECUTABLE_NAME "rundm90_ndm")
if ("MPI" IN_LIST ${PROJECT_NAME}_PACKAGE_LIST)
    string(APPEND MAIN_EXECUTABLE_NAME "_mpi")
else()
    string(APPEND MAIN_EXECUTABLE_NAME "_serial")
endif()
add_library( ${MAIN_LIBRARY_NAME} STATIC )

# adding required libraries to main library
foreach( PACKAGE IN ITEMS ${${PROJECT_NAME}_PACKAGE_LIST})
    target_link_libraries(${MAIN_LIBRARY_NAME} ${${PROJECT_NAME}_${PACKAGE}_LIBRARIES})
    target_link_directories(${MAIN_LIBRARY_NAME} PUBLIC ${${PROJECT_NAME}_${PACKAGE}_LIBRARY_DIRS})
    target_include_directories( ${MAIN_LIBRARY_NAME} PUBLIC ${${PROJECT_NAME}_${PACKAGE}_INCLUDE_DIR} )
endforeach()

# reading src/CMakeLists.txt
add_subdirectory( src )

if (${PROJECT_NAME}_OPT_COMPILE_TESTS)
    log_warning( "add_subdirectory tests" )
    add_subdirectory( tests )
endif()

# Defines target to install
install( TARGETS ${MAIN_LIBRARY_NAME} ${MAIN_EXECUTABLE_NAME})

################################################
# PKG-Config
# generates .pc file at installation
################################################

install(CODE "
    set(PROJECT_NAME ${PROJECT_NAME})
    set(MAIN_LIBRARY_NAME ${MAIN_LIBRARY_NAME})
    set(PROJECT_HOMEPAGE_URL ${PROJECT_HOMEPAGE_URL})
    set(PROJECT_VERSION ${PROJECT_VERSION})
    set(PROJECT_DESCRIPTION \"${PROJECT_DESCRIPTION}\")
    set(PC_FILE_PATH \"\${CMAKE_INSTALL_PREFIX}/pkgconfig\")
    configure_file(${CMAKE_SOURCE_DIR}/cmake_files/liblibrary.pc.in \${PC_FILE_PATH}/lib${MAIN_LIBRARY_NAME}.pc @ONLY)
    message(STATUS \"pkg-config file: \${PC_FILE_PATH}/lib${MAIN_LIBRARY_NAME}.pc\" )
")

################################################
# DOXYGEN Documentation
################################################
# Must be after source files definition

if( ${PROJECT_NAME}_OPT_COMPILE_DOC )
  doxygen_doc()
endif()

################################################
# Trace some logs
################################################

# Display options
get_cmake_property(${PROJECT_NAME}_VARIABLE_NAMES VARIABLES)                
foreach( VARIABLE_NAME in ${${PROJECT_NAME}_VARIABLE_NAMES})
    if(VARIABLE_NAME MATCHES "._OPT_.")
        log_info("${VARIABLE_NAME} : ${${VARIABLE_NAME}}")
    endif()
endforeach()

# Display compilers and flags
foreach (LANGUAGE IN ITEMS ${${PROJECT_NAME}_LANGUAGES})
    log_info("CMAKE_${LANGUAGE}_COMPILER : ${CMAKE_${LANGUAGE}_COMPILER}")
    log_info("CMAKE_${LANGUAGE}_FLAGS_${CMAKE_BUILD_TYPE} : ${CMAKE_${LANGUAGE}_FLAGS_${CMAKE_BUILD_TYPE}}")
endforeach()

# Display preprocessor definitions
get_directory_property(COMPILE_DEFINITIONS_LIST COMPILE_DEFINITIONS)
foreach (COMP_DEF IN ITEMS ${COMPILE_DEFINITIONS_LIST})
    log_warning("Using preprocessor definition : ${COMP_DEF}")
endforeach()

# Display all variables in cmake, optionnal
if( ${PROJECT_NAME}_OPT_TRACE )
  printCmakeTrace()
endif()

