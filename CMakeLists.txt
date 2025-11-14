CMAKE_minimum_required( VERSION 3.21 )

# cmake scripts location
list( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_files" )        
include( utils )

################################################
# Auto configuration
# in case no configuration file is used
################################################
# Cache variables will not be reset if already defined with cmake .. -C <cache file>

# met CMAKE_BUILD_TYPE en majuscules
if(CMAKE_BUILD_TYPE)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _tmp)
    set(CMAKE_BUILD_TYPE "${_tmp}" CACHE STRING "RELEASE or DEBUG" FORCE)
endif()

# release mode by default
set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "RELEASE or DEBUG")

# List fortran compilers to look for
#set(INTEL_Fortran_COMPILERS ifx ifort)
set(GNU_Fortran_COMPILERS gfortran)
#set(INTEL_C_COMPILERS icc)
set(GNU_C_COMPILERS gcc)
#set(INTEL_CXX_COMPILERS i++)
set(GNU_CXX_COMPILERS g++)

# list configurations, priority to first ones
#set(compile_configuration_list GNU INTEL)
set(compile_configuration_list GNU)

# Define possible configurations
# Compilation with INTEL not possible for now

# macro(CONFIG_INTEL) # INTEL
#     set(CMAKE_Fortran_COMPILER "${VAR_INTEL_Fortran}" CACHE STRING "")    
#     cmake_path(GET CMAKE_Fortran_COMPILER FILENAME compiler_name)     
#     set(DK-IO_PRESET_TITLE "Compile DK-IO INTEL : ${compiler_name} + intelmpi" CACHE STRING "")
#     set(CMAKE_Fortran_FLAGS_RELEASE " -O3 -fpp -diag-disable=10448 " CACHE STRING "")
#     set(CMAKE_Fortran_FLAGS_DEBUG " -O0 -g -C -fpe-all=0" CACHE STRING "")
# endmacro(CONFIG_INTEL)

macro(CONFIG_GNU) # GNU
    set(CMAKE_Fortran_COMPILER "${VAR_GNU_Fortran}" CACHE STRING "")    
    cmake_path(GET CMAKE_Fortran_COMPILER FILENAME compiler_name)     
    set(DK-IO_PRESET_TITLE "Compile DK-IO GNU : ${compiler_name} + openmpi" CACHE STRING "")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -g -march=native -Wimplicit-procedure" CACHE STRING "")
    set(CMAKE_Fortran_FLAGS_DEBUG " -O0 -g -C -Wimplicit-procedure" CACHE STRING "")
    set(CMAKE_C_FLAGS_RELEASE "-g -Ofast -march=native" CACHE STRING "")
    set(CMAKE_C_FLAGS_DEBUG "-g -O0" CACHE STRING "")
    set(CMAKE_CXX_FLAGS_RELEASE "-g -Ofast -march=native" CACHE STRING "")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "")
endmacro(CONFIG_GNU)

# Look for available compilers
# 1. Check if environment variable OMPI_FC, MPI_FC or FC is defined
# 2. Check for each configuration (INTEL, GNU) if a compiler exist
set_compiler("DK-IO" "FC" "Fortran")
set_compiler("DK-IO" "CC" "C")
set_compiler("DK-IO" "CC" "CXX")

# Determine the configuration
# if a suitable compiler has been found, run the macro CONFIG_<config>
set_config("DK-IO" "Fortran")
set_config("DK-IO" "C")
set_config("DK-IO" "CXX")

################################################
# Project definition
################################################

project( DK-IO
         VERSION 0.1
         DESCRIPTION "DK-IO project compiled with CMake"
         HOMEPAGE_URL ""
         LANGUAGES "Fortran" "C" "CXX" )

set( PROJECT_AUTHOR "Paul Fossati" )

################################################
# default user options
################################################
if (PROJECT_IS_TOP_LEVEL) # cmake 3.21
    #option( ${PROJECT_NAME}_OPT_TEST "Compile tests utilities" OFF )
    #option( ${PROJECT_NAME}_OPT_COVERAGE "compile for code coverage analysis" OFF )
    option( ${PROJECT_NAME}_OPT_TRACE "trace all variables for cmake debug" OFF )
    option( DK-IO_OPT_COMPILE_BIN "compile cfginfo, cfgconvert, and cfgdiff" ON )
endif()

################################################
#   set targets and sources
################################################

# Add the dk_core project here:
find_package(DKCore QUIET)
if(DKCore_FOUND)
    set(DKCore_HINT_DIR "${DKCore_DIR}") #Pour garder en mémoire l'emplacement de DKCore dans DKIOConfig.cmake
else()
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../dk_core/CMakeLists.txt")
        add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../dk_core dk_core-build)
        add_library(DKCore::dk_core ALIAS dk_core)
    else()
        message(FATAL_ERROR "La bibliothèque dk_core n'a pas été trouvée")
    endif()
endif()

add_subdirectory(src)

# Defines target to install
if (PROJECT_IS_TOP_LEVEL)
    install(TARGETS dk_io cfginfo cfgconvert cfgdiff 
        EXPORT DKIO
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)
    install(DIRECTORY ${CMAKE_BINARY_DIR}/mod/dk_io/
        DESTINATION mod/dk_io
        FILES_MATCHING PATTERN "*.mod"
    )
    install(EXPORT DKIO
        FILE DKIOTargets.cmake
        NAMESPACE DKIO::
        DESTINATION lib/cmake/DKIO/
    )
    include(CMakePackageConfigHelpers)
    configure_package_config_file(
        cmake_files/DKIOConfig.cmake.in
        ${CMAKE_BINARY_DIR}/cmake/DKIOConfig.cmake
        INSTALL_DESTINATION lib/cmake/DKIO/
    )
    install(FILES
        ${CMAKE_BINARY_DIR}/cmake/DKIOConfig.cmake
        DESTINATION lib/cmake/DKIO/
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
    get_target_property(COMPILE_DEFINITIONS_LIST dk_io COMPILE_DEFINITIONS)
    foreach (COMP_DEF IN ITEMS ${COMPILE_DEFINITIONS_LIST})
        log_info("Using preprocessor definition : ${COMP_DEF}")
    endforeach()

    # Display all variables in cmake, optionnal
    if( ${PROJECT_NAME}_OPT_TRACE )
    printCmakeTrace()
    endif()
endif()