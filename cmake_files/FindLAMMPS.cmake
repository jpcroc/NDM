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
Variable to set in Cmake or as an environment variable. By default "lammps".
Depends on <LAMMPS_MACHINE> custom parameter which may be set at lammps compilation and is related to the executable name:
executable name : lmp_<LAMMPS_MACHINE>
library name: lammps_<LAMMPS_MACHINE>
For example "lammps_mpi" is popular for parallel lammps but there is no convention.

#]]
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
            list(APPEND RESULT_GET_DEFINITION ${TMP_DEF})
        endif()
    endforeach()
endfunction(get_definition)

# Macro to try looking for lammps using pkg-config
macro(pkg_search)
    find_package(PkgConfig QUIET)
    if (${PKG_CONFIG_FOUND})
        pkg_check_modules(LAMMPS "lib${LAMMPS_LIBRARY_NAME}")
        get_definition("${LAMMPS_CFLAGS}")
        list(APPEND LAMMPS_DEFINITIONS ${RESULT_GET_DEFINITION})
    endif()
endmacro(pkg_search)

# Macro to try finding lammps library in LD_LIBRARY_PATH environment variable
macro (ld_library_search)
    string(REPLACE ":" ";" TMP_LD_PATH $ENV{LD_LIBRARY_PATH})
    find_library(TMP_PATH_TO_LIBRARY ${LAMMPS_LIBRARY_NAME} HINTS ${TMP_LD_PATH})
    if (NOT TMP_PATH_TO_LIBRARY STREQUAL "TMP_PATH_TO_LIBRARY-NOTFOUND")
        set(LAMMPS_FOUND ON)
        cmake_path(REMOVE_FILENAME TMP_PATH_TO_LIBRARY)
        list(APPEND LAMMPS_LIBRARY_DIRS ${TMP_PATH_TO_LIBRARY})
        list(APPEND LAMMPS_INCLUDE "${TMP_PATH_TO_LIBRARY}/../include" )
    endif()
endmacro(ld_library_search)

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

# Try to find Lammps with pkg-config
pkg_search()

if (NOT LAMMPS_FOUND)
    # look for lammps libraries in LD_LIBRARY_PATH
    ld_library_search()  
    list(APPEND LAMMPS_LIBRARIES "${LAMMPS_LIBRARY_NAME}")
    list(APPEND LAMMPS_DEFINITIONS "LAMMPS_SMALLBIG")
endif()

# correct some include directories
foreach (TMP_DIR IN ITEMS ${LAMMPS_INCLUDE_DIRS})
    if (EXISTS "${TMP_DIR}/lammps")
        list(APPEND LAMMPS_INCLUDE_DIRS "${TMP_DIR}/lammps")
    endif()
endforeach()

# search for extra libraries
set(LAMMPS_MODULES_LIST     # List of possible extra libraries for lammps
    meam
    voro++
    gsl
    plumed
    milady
    lammps_atc_mpi
    lammps_qmmm_mpi
    colvars
    reax
    poems
    awpmd
)

if (NOT TMP_PATH_TO_LIBRARY STREQUAL "TMP_PATH_TO_LIBRARY-NOTFOUND")
    foreach (EXTRA_LIB IN ITEMS ${LAMMPS_MODULES_LIST})
        find_library(TMP_PATH_TO_EXTRA_LIB ${EXTRA_LIB} HINTS ${TMP_PATH_TO_LIBRARY})
        if (NOT TMP_PATH_TO_EXTRA_LIB STREQUAL "TMP_PATH_TO_EXTRA_LIB-NOTFOUND")
            list(APPEND LAMMPS_LIBRARIES ${EXTRA_LIB})
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
