CMAKE_minimum_required( VERSION 3.21 )

# cmake scripts location
list( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_files" )
include( utils )

# met CMAKE_BUILD_TYPE en majuscules
if(CMAKE_BUILD_TYPE)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _tmp)
    set(CMAKE_BUILD_TYPE "${_tmp}" CACHE STRING "RELEASE or DEBUG" FORCE)
endif()

# release mode by default
set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "RELEASE or DEBUG")

# List fortran compilers to look for
# set(INTEL_Fortran_COMPILERS ifx ifort)
set(GNU_Fortran_COMPILERS gfortran)
# set(INTEL_C_COMPILERS icc icx)
set(GNU_C_COMPILERS gcc)

# list configurations, priority to first ones
set(compile_configuration_list GNU)

# Define possible configurations
# Compilation with INTEL not possible for now

# macro(CONFIG_INTEL) # INTEL
#     set(CMAKE_Fortran_COMPILER "${VAR_INTEL_Fortran}" CACHE STRING "")    
#     cmake_path(GET CMAKE_Fortran_COMPILER FILENAME compiler_name)     
#     set(DK-Core_PRESET_TITLE "Compile DK-Core INTEL : ${compiler_name} + intelmpi" CACHE STRING "")
#     set(CMAKE_Fortran_FLAGS_RELEASE " -O3 -fpp -diag-disable=10448 " CACHE STRING "")
#     set(CMAKE_Fortran_FLAGS_DEBUG " -O0 -g -C -fpe-all=0" CACHE STRING "")
# endmacro(CONFIG_INTEL)

macro(CONFIG_GNU) # GNU
    set(CMAKE_Fortran_COMPILER "${VAR_GNU_Fortran}" CACHE STRING "")    
    cmake_path(GET CMAKE_Fortran_COMPILER FILENAME compiler_name)     
    set(DK-Core_PRESET_TITLE "Compile DK-Core GNU : ${compiler_name} + openmpi" CACHE STRING "")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -g -march=native -Wimplicit-procedure" CACHE STRING "")
    set(CMAKE_Fortran_FLAGS_DEBUG " -O0 -g -C -Wimplicit-procedure" CACHE STRING "")
    set(CMAKE_C_FLAGS_RELEASE "-g -Ofast -march=native" CACHE STRING "")
    set(CMAKE_C_FLAGS_DEBUG "-g -O0" CACHE STRING "")
endmacro(CONFIG_GNU)

# Look for available compilers
# 1. Check if environment variable OMPI_FC, MPI_FC or FC is defined
# 2. Check for each configuration (INTEL, GNU) if a compiler exist
set_compiler("DK-Core" "FC" "Fortran")
set_compiler("DK-Core" "CC" "C")

# Determine the configuration
# if a suitable compiler has been found, run the macro CONFIG_<config>
set_config("DK-Core" "Fortran")
set_config("DK-Core" "C")

################################################
# Project definition
################################################

project( DK-Core
         VERSION 0.1
         DESCRIPTION "DK-Core project compiled with CMake"
         HOMEPAGE_URL ""
         LANGUAGES "Fortran" "C")

set( PROJECT_AUTHOR "Paul Fossati" )

################################################
# default user options
################################################
if (PROJECT_IS_TOP_LEVEL) # cmake 3.21
    #option( ${PROJECT_NAME}_OPT_TEST "Compile tests utilities" OFF )
    #option( ${PROJECT_NAME}_OPT_COVERAGE "compile for code coverage analysis" OFF )
    option( ${PROJECT_NAME}_OPT_TRACE "trace all variables for cmake debug" OFF )
endif()

################################################
#   set targets and sources
################################################

# listing sources
add_subdirectory(src)

# Defines target to install
if (PROJECT_IS_TOP_LEVEL)
    install(TARGETS dk_core 
        EXPORT DKCore
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
    install(DIRECTORY ${CMAKE_BINARY_DIR}/mod/dk_core/
        DESTINATION mod/dk_core
        FILES_MATCHING PATTERN "*.mod"
    )
    install(EXPORT DKCore
        FILE DKCoreConfig.cmake
        NAMESPACE DKCore::
        DESTINATION lib/cmake/DKCore/
    )
endif()
################################################
# Trace some logs
################################################

# Get the languages list for this project
get_property(${PROJECT_NAME}_LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)   
include( TraceNewMacros )   # colorful logs
setcolors()                 # setting log colors

if (PROJECT_IS_TOP_LEVEL) # cmake 3.21
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
        log_info("CMAKE_${LANGUAGE}_COMPILER_ID : ${CMAKE_${LANGUAGE}_COMPILER_ID}")
        log_info("CMAKE_${LANGUAGE}_FLAGS_${CMAKE_BUILD_TYPE} : ${CMAKE_${LANGUAGE}_FLAGS_${CMAKE_BUILD_TYPE}}")
    endforeach()

    # Display preprocessor definitions
    get_target_property(COMPILE_DEFINITIONS_LIST dk_core COMPILE_DEFINITIONS)
    foreach (COMP_DEF IN ITEMS ${COMPILE_DEFINITIONS_LIST})
        log_info("Using preprocessor definition : ${COMP_DEF}")
    endforeach()

    # Display all variables in cmake, optionnal
    if( ${PROJECT_NAME}_OPT_TRACE )
    printCmakeTrace()
    endif()
endif()