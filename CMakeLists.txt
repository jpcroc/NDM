CMAKE_minimum_required( VERSION 3.20 )
################################################
# Auto configuration
# in case no configuration file is used
################################################
# Cache variables will not be reset if already defined with cmake .. -C <cache file>


set(languages "Fortran")
if ("LAMMPS" IN_LIST NDM_PACKAGE_LIST)
    list(APPEND languages "CXX")
endif()

set(compile_configuration_list INTEL GNU)

# List fortran mpi wrappers to look for
set(INTEL_Fortran_WRAPPERS mpiifx mpiifort mpifort)
set(INTEL_Fortran_COMPILERS ifx ifort)
set(INTEL_CXX_WRAPPERS mpiicpc)
set(INTEL_CXX_COMPILERS icpc)

set(GNU_Fortran_WRAPPERS mpifort mpifc mpif90 mpif77)
set(GNU_Fortran_COMPILERS gfortran f95 f77)
set(GNU_CXX_WRAPPERS mpicxx mpiCC mpic++ mpigxx mpigcc mpicc )
set(GNU_CXX_COMPILERS g++ gcc)

# define configurations
macro(CONFIG_INTEL)
    set(CMAKE_${lang}_COMPILER "${VAR_${config}_${lang}}" CACHE STRING "")
    cmake_path(GET CMAKE_${lang}_COMPILER FILENAME compiler_name)     
    set(NDM_PRESET_TITLE "Compile NDM INTEL: ${compiler_name}" CACHE STRING "")
endmacro(CONFIG_INTEL)

macro(CONFIG_GNU)
    set(CMAKE_${lang}_COMPILER "${VAR_${config}_${lang}}" CACHE STRING "")
    cmake_path(GET CMAKE_${lang}_COMPILER FILENAME compiler_name)
    set(NDM_PRESET_TITLE "Compile NDM GNU : ${compiler_name}" CACHE STRING "")
    set(CMAKE_Fortran_FLAGS_RELEASE " -O3 -ffree-line-length-512 -funroll-all-loops" CACHE STRING "")
    set(CMAKE_Fortran_FLAGS_DEBUG " -O0 -g -C -fpe-all=0 -ffree-line-length-512" CACHE STRING "")
    set(MKL_MPI "openmpi" CACHE STRING "MPI")
endmacro(CONFIG_GNU)

foreach(lang IN ITEMS ${languages})

# Find mpirun directory to look for wrapper
set(list_MPIEXEC mpiexec mpirun)
set(VAR_MPIEXEC "VAR_MPIEXEC-NOTFOUND" CACHE STRING "")
foreach(mpi_prog IN ITEMS ${list_MPIEXEC})
	if (VAR_MPIEXEC STREQUAL "VAR_MPIEXEC-NOTFOUND")
	    if (DEFINED ENV{OMP_ROOT})
	        find_program(VAR_MPIEXEC ${mpi_prog} PATHS $ENV{OMP_ROOT} NO_DEFAULT_PATH)
	    elseif(DEFINED ENV{MPI_HOME})
	        find_program(VAR_MPIEXEC ${mpi_prog} PATHS $ENV{MPI_HOME} NO_DEFAULT_PATH)
	    else()
	        find_program(VAR_MPIEXEC ${mpi_prog})	# search anywhere
            endif()
	endif()
endforeach()

if (VAR_MPIEXEC STREQUAL "VAR_MPIEXEC-NOTFOUND")
	set(Project_PARA OFF)
else()
	cmake_path(GET VAR_MPIEXEC PARENT_PATH mpi_path)
	set(Project_PARA ON)
endif()

# For each possible configuration, look for a wrapper in the mpirun directory
set(WRAPPER_FOUND_${lang} False)
foreach (config IN ITEMS ${compile_configuration_list})
    set(VAR_${config}_${lang} "VAR_${config}_${lang}-NOTFOUND" CACHE STRING "")
    foreach (wrapper IN ITEMS ${${config}_${lang}_WRAPPERS})
        if (VAR_${config}_${lang} STREQUAL "VAR_${config}_${lang}-NOTFOUND")
            find_program(VAR_${config}_${lang} ${wrapper} PATHS ${mpi_path} NO_DEFAULT_PATH)
        endif()
    endforeach()
    if (NOT VAR_${config}_${lang} STREQUAL "VAR_${config}_${lang}-NOTFOUND") # if a wrapper is found, get the underlying compiler        
	set(WRAPPER_FOUND_${lang} True)
	execute_process(
            COMMAND ${VAR_${config}_${lang}} -show
            OUTPUT_VARIABLE SHOW_UNDERLYING)
        if (SHOW_UNDERLYING)
            string(REPLACE " " ";" SHOW_UNDERLYING2 ${SHOW_UNDERLYING})
            list(GET SHOW_UNDERLYING2 0 ${config}_${lang}_underlying_compiler)
        endif()
    endif()
endforeach()

# If no wrapper was found, look wider
if (NOT ${WRAPPER_FOUND_${lang}})
	foreach (config IN ITEMS ${compile_configuration_list})
	    foreach (compiler IN ITEMS ${${config}_${lang}_COMPILERS})
		if (VAR_${config}_${lang} STREQUAL "VAR_${config}_${lang}-NOTFOUND")
		    find_program(VAR_${config}_${lang} ${compiler})
		endif()
	    endforeach()
	    if (NOT VAR_${config}_${lang} STREQUAL "VAR_${config}_${lang}-NOTFOUND") # if a compiler is found, maybe it's a wrapper and we can get the underlying compiler
		execute_process(
		    COMMAND ${VAR_${config}_${lang}} -show
		    OUTPUT_VARIABLE SHOW_UNDERLYING)
		if (SHOW_UNDERLYING)
		    string(REPLACE " " ";" SHOW_UNDERLYING2 ${SHOW_UNDERLYING})
		    list(GET SHOW_UNDERLYING2 0 ${config}_${lang}_underlying_compiler)
		endif()
	    endif()
	endforeach()
endif()

# Choose previously defined configuration
set(config_${lang}_set OFF)
foreach (config IN ITEMS ${compile_configuration_list})
    if (NOT ${config_${lang}_set})
	# if the underlying compiler is defined it means we have found a working wrapper		
	if (DEFINED ${config}_${lang}_underlying_compiler)
	    if (${${config}_${lang}_underlying_compiler} IN_LIST ${config}_${lang}_COMPILERS)
                cmake_language(CALL CONFIG_${config})        #cmake 3.18
                set(config_${lang}_set ON)
            endif()
	elseif(DEFINED CMAKE_${lang}_COMPILER)
	    cmake_path(GET CMAKE_${lang}_COMPILER FILENAME compiler_name)
            if (${compiler_name} IN_LIST ${config}_${lang}_COMPILERS)
		cmake_language(CALL CONFIG_${config})        #cmake 3.18
		set(config_${lang}_set ON)
            endif()
	elseif (NOT ${Project_PARA} AND NOT VAR_${config}_${lang} STREQUAL "VAR_${config}_${lang}-NOTFOUND")
		cmake_language(CALL CONFIG_${config})        #cmake 3.18
                set(config_${lang}_set ON)
        endif()
    endif()
endforeach()

endforeach()

if (${Project_PARA})
    list(APPEND NDM_PACKAGE_LIST "MPI")
    list(APPEND NDM_COMPILE_DEFINITION "PARA")
    list(APPEND NDM_COMPILE_DEFINITION "PARAML")
endif()

set(CMAKE_Fortran_FLAGS_RELEASE " -O3 " CACHE STRING "")
set(CMAKE_Fortran_FLAGS_DEBUG " -O0 -g -C -fpe-all=0" CACHE STRING "")

if(NOT NDM_WITHOUT_MKL)
    list(APPEND NDM_PACKAGE_LIST "MKL")
    list(APPEND NDM_COMPILE_DEFINITION "MKL")
endif()

set(MKL_MPI "intelmpi" CACHE STRING "MPI")
set(ENABLE_SCALAPACK OFF CACHE STRING "Scalapack" )
set(ENABLE_LAPACK95 ON CACHE STRING "lapack fortran 95 api" )
set(MKL_THREADING "sequential" CACHE STRING "Threading")
set(MKL_INTERFACE "lp64" CACHE STRING "Interface")
set(MKL_ARCH "intel64" CACHE STRING "Architecture")

if ("LAMMPS" IN_LIST NDM_PACKAGE_LIST)
    set(CMAKE_CXX_FLAGS_RELEASE " -no-multibyte-chars " CACHE STRING "")
endif()

################################################
# Project definition
################################################
project( NDM
         VERSION 3.0.0
         DESCRIPTION "NDM project with CMake"
         HOMEPAGE_URL "https://codev-tuleap.intra.cea.fr/plugins/git/ndm"
         LANGUAGES ${languages})

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

