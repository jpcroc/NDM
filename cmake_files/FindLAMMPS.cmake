#[[
FindLAMMPS

Tries to find lammps using pkg-config and some more hazardous ways.
Called with Find_package(LAMMPS).

PKG-CONFIG
----------
To use pkg-config successfully make sure the .pc file (usually liblammps.pc) is included in the PKG_CONFIG_PATH.
The .pc file name depends on lammps library name, see below.

LAMMPS_LIBRARY_NAME
-------------------
Variable to set in CMake or as an environment variable. By default "lammps".
Depends on <LAMMPS_MACHINE> custom parameter which may be set at lammps compilation and is related to the executable name:
executable name : lmp_<LAMMPS_MACHINE>
library name: lammps_<LAMMPS_MACHINE>
For example "lammps_mpi" is popular for parallel lammps but there is no convention.

pkg-config will be looking for lib<LAMMPS_LIBRARY_NAME>.pc file.

LAMMPS_HOME
-----------
Variable to set in CMake or as an environment variable, undefined by default.
Path to lammps library directory or to lammps .pc file directory.
Typically this will be the path to LAMMPS build directory if LAMMPS has been compiled with CMake or the src/ directory if compiled with make in the sources.

#]]
cmake_minimum_required(VERSION 3.20)
################################################
# Functions and Macros definitions
################################################

# Function to extract definition variables from a list of cflags
# Definitions variables are stored in RESULT_GET_DEFINITION variable
function (get_definition TMP_FLAGS)
    set(RESULT_GET_DEFINITION "")
    foreach(TMP_FLAG IN ITEMS ${TMP_FLAGS})
        if ( ${TMP_FLAG} MATCHES "^-D" ) # if flag starts with -D
            string(SUBSTRING ${TMP_FLAG} 2 -1 TMP_DEF) # ignore the two first characters
	    string(FIND ${TMP_DEF} "$" INCORRECT_DEFINITION) # Check for illegal character
	    if (INCORRECT_DEFINITION EQUAL -1)
		list(APPEND RESULT_GET_DEFINITION ${TMP_DEF})
	    endif()
        endif()
    endforeach()
    # return result
    set(RESULT_GET_DEFINITION ${RESULT_GET_DEFINITION} PARENT_SCOPE)
endfunction(get_definition)

# Macro to try looking for lammps using pkg-config
macro(pkg_search)
    find_package(PkgConfig QUIET)
    if (${PKG_CONFIG_FOUND})
        pkg_check_modules(LAMMPS "lib${LAMMPS_LIBRARY_NAME}")
        if(LAMMPS_FOUND)
            get_definition("${LAMMPS_CFLAGS}")
            list(APPEND LAMMPS_DEFINITIONS ${RESULT_GET_DEFINITION})
        endif()
        
    endif()
endmacro(pkg_search)

# Macro to try finding lammps library in LD_LIBRARY_PATH environment variable
macro (ld_library_search)
    string(REPLACE ":" ";" TMP_LD_PATH $ENV{LD_LIBRARY_PATH})
    find_library(TMP_PATH_TO_LIBRARY ${LAMMPS_LIBRARY_NAME} HINTS ${TMP_LD_PATH})
    if (NOT TMP_PATH_TO_LIBRARY STREQUAL "TMP_PATH_TO_LIBRARY-NOTFOUND")
        set(LAMMPS_FOUND ON)
        cmake_path(REMOVE_FILENAME TMP_PATH_TO_LIBRARY)
        append_if_new_library_dirs(${TMP_PATH_TO_LIBRARY})
        foreach (TMP_DIR IN ITEMS ${POSSIBLE_INCLUDE_DIRECTORIES})
            append_if_new_include_dirs("${TMP_PATH_TO_LIBRARY}${TMP_DIR}")
        endforeach()
    endif()
endmacro(ld_library_search)

macro (append_if_exists DIRECTORY DIRECTORY_LIST)
    if (EXISTS ${DIRECTORY})
        set(NEW_DIRECTORY_LIST "")
        set(TMP_DIRECTORY_LIST "")
        list(APPEND TMP_DIRECTORY_LIST ${DIRECTORY_LIST} ${DIRECTORY})
        foreach(TMP_DIRECTORY IN ITEMS ${TMP_DIRECTORY_LIST})
            cmake_path(SET TMP_DIRECTORY NORMALIZE "${TMP_DIRECTORY}/")
            if (NOT ${TMP_DIRECTORY} IN_LIST NEW_DIRECTORY_LIST)
                list(APPEND NEW_DIRECTORY_LIST ${TMP_DIRECTORY})
            endif()
        endforeach()
    else()
        set(NEW_DIRECTORY_LIST ${DIRECTORY_LIST})
    endif()
endmacro()

function (append_if_new_include_dirs DIRECTORY)
    append_if_exists(${DIRECTORY} "${LAMMPS_INCLUDE_DIRS}")
    set(LAMMPS_INCLUDE_DIRS ${NEW_DIRECTORY_LIST} PARENT_SCOPE)
endfunction(append_if_new_include_dirs)

function (append_if_new_library_dirs DIRECTORY)
    append_if_exists(${DIRECTORY} "${LAMMPS_LIBRARY_DIRS}")
    set(LAMMPS_LIBRARY_DIRS ${NEW_DIRECTORY_LIST} PARENT_SCOPE)
endfunction(append_if_new_library_dirs)

################################################
# Sequential script
################################################

# Standard arguments for Find_package function
include(FindPackageHandleStandardArgs)

# Initialize LAMMPS_LIBRARY_NAME
set(LAMMPS_LIBRARY_NAME "lammps" CACHE STRING "Name of the lammps library")
if (DEFINED ENV{LAMMPS_LIBRARY_NAME})
    set(LAMMPS_LIBRARY_NAME $ENV{LAMMPS_LIBRARY_NAME})
endif()
set(LAMMPS_DEFINITIONS "")
set(POSSIBLE_INCLUDE_DIRECTORIES 
    "/../include" 
    "/includes" 
    "/../src")

# Try to find Lammps with pkg-config
if (DEFINED LAMMPS_HOME)
    set(ENV{PKG_CONFIG_PATH} "${LAMMPS_HOME}:$ENV{PKG_CONFIG_PATH}")
    list(PREPEND CMAKE_PREFIX_PATH "${LAMMPS_HOME}")
    set(pc_file "${LAMMPS_HOME}/lib${LAMMPS_LIBRARY_NAME}.pc")
    if(EXISTS ${pc_file})
        pkg_search()
        if (LAMMPS_FOUND)
            list(APPEND LAMMPS_LIBRARY_DIRS "${LAMMPS_HOME}")
            append_if_new_library_dirs(${LAMMPS_HOME})
            foreach (TMP_DIR IN ITEMS ${POSSIBLE_INCLUDE_DIRECTORIES})
                append_if_new_include_dirs("${LAMMPS_HOME}${TMP_DIR}")
            endforeach()
        endif()
    endif()
else()
    pkg_search()
endif()

# if pkg-config didn't work, try with LD_LIBRARY_PATH
if (NOT LAMMPS_FOUND)
    # look for lammps libraries in LD_LIBRARY_PATH
    ld_library_search()  
    list(APPEND LAMMPS_LIBRARIES "${LAMMPS_LIBRARY_NAME}")
    list(APPEND LAMMPS_DEFINITIONS "LAMMPS_SMALLBIG")
endif()

# correct some include directories
foreach (TMP_DIR IN ITEMS ${LAMMPS_INCLUDE_DIRS})
    append_if_new_include_dirs("${TMP_DIR}/lammps")
endforeach()

# search for extra libraries
set(LAMMPS_MODULES_LIST     # List of possible extra libraries for lammps
    atc
    awpmd
    colvars
    gomp
    gsl
    jpeg
    lammps_atc_mpi
    lammps_qmmm_mpi
    meam
    milady
    mpi_stubs
    plumed
    png
    poems
    qmmm
    reax
    voro++
    z
)

list(APPEND POSSIBLE_LIB_DIRS ${TMP_LD_PATH})
list(PREPEND POSSIBLE_LIB_DIRS
#	"/usr/lib64/"
)
list(PREPEND POSSIBLE_SUFFIX_DIRS 
	"/STUBS"
	"voro_build-prefix/src/voro_build/src/"
)

if (LAMMPS_FOUND)
    foreach (EXTRA_LIB IN ITEMS ${LAMMPS_MODULES_LIST})
        set(TMP_PATH_TO_EXTRA_LIB "TMP_PATH_TO_EXTRA_LIB-NOTFOUND")
        find_library(TMP_PATH_TO_EXTRA_LIB ${EXTRA_LIB} HINTS ${POSSIBLE_LIB_DIRS} PATH_SUFFIXES ${POSSIBLE_SUFFIX_DIRS})
        if (NOT TMP_PATH_TO_EXTRA_LIB STREQUAL "TMP_PATH_TO_EXTRA_LIB-NOTFOUND")
            list(APPEND LAMMPS_LIBRARIES ${EXTRA_LIB})

            cmake_path(REMOVE_FILENAME TMP_PATH_TO_EXTRA_LIB)
            append_if_new_library_dirs(${TMP_PATH_TO_EXTRA_LIB})
        endif()
    endforeach()
endif()

# required variables for find_package()
find_package_handle_standard_args( LAMMPS 
    REQUIRED_VARS LAMMPS_LIBRARY_DIRS LAMMPS_LIBRARIES 
    VERSION_VAR LAMMPS_VERSION
)
# optionnal output variables :
set(LAMMPS_LINK_LIBRARIES ${LAMMPS_LINK_LIBRARIES} PARENT_SCOPE)
set(LAMMPS_CFLAGS ${LAMMPS_CFLAGS} PARENT_SCOPE)
set(LAMMPS_LDFLAGS ${LAMMPS_LDFLAGS} PARENT_SCOPE)

set(LAMMPS_INCLUDE_DIRS ${LAMMPS_INCLUDE_DIRS} PARENT_SCOPE)
set(LAMMPS_DEFINITIONS ${LAMMPS_DEFINITIONS} PARENT_SCOPE)
