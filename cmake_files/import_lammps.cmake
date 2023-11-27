###########
#
#   MODULE TO FIND LAMMPS LIBRARY
#
#   Make sure lammps.pc is included in the PKG_CONFIG_PATH or
#   indicate lammps path manually using LAMMPS_PATH and lammps
#   ldflags and cflags using LAMMPS_FLAGS.
#
#
############


function (import_lammps)
    
    add_compile_definitions(LAMMPS_LIB_MPI=1)
    message(STATUS "test lammps --------------------------------------------------")
    
    if(NOT LAMMPS_FLAGS OR NOT LAMMPS_PATH)                 # if information isn't provided, call pkgconfig
        find_package(PkgConfig QUIET)
        pkg_check_modules(LAMMPS liblammps)
    else()
        message(STATUS "go on.")
    endif()

    set(LAMMPS_FLAGS "" CACHE STRING "necessary flags to include lammps")
    set(LAMMPS_PATH "" CACHE STRING "lammps library path, liblammps.so or liblammps.a")
    
    
    if (LAMMPS_PATH AND EXISTS ${LAMMPS_PATH})      # a valid path is given
        message(STATUS "Using the given LAMMPS_PATH=${LAMMPS_PATH}")
    elseif (LAMMPS_PATH AND NOT EXISTS ${LAMMPS_PATH})
        message(STATUS "The provided path does not exist: LAMMPS_PATH=${LAMMPS_PATH}")
        message(STATUS "Please provide the full path to liblammps.so or liblammps.a with '-DLAMMPS_PATH=path/to/library/liblammps.so' including the file name or do not define LAMMPS_PATH to let pkg-config search for lammps.")
    elseif(${LAMMPS_FOUND}) # Find path automatically with pkg-config
        set(LAMMPS_PATH_SHARED "${LAMMPS_LIBDIR}/lib${LAMMPS_LIBRARIES}.so")
        set(LAMMPS_PATH_STATIC "${LAMMPS_LIBDIR}/lib${LAMMPS_LIBRARIES}.a")
    
        if (EXISTS ${LAMMPS_PATH_SHARED})
            message(STATUS "Found Lammps: ${LAMMPS_PATH_SHARED}")
            set(LAMMPS_PATH ${LAMMPS_PATH_SHARED})
        elseif (EXISTS ${LAMMPS_PATH_STATIC})
            message(STATUS "Found Lammps: ${LAMMPS_PATH_STATIC}")
            set(LAMMPS_PATH ${LAMMPS_PATH_STATIC})
        else()
            message(STATUS "Couldn't find liblammps.so or liblammps.a. Please provide the correct path with '-DLAMMPS_PATH=path/to/library/liblammps.so' including the file name.")
        endif()
    else()
        message(STATUS "Couldn't find lammps using pkg-config.")
        message(STATUS "Please verify that lammps.pc can be found in PKG_CONFIG_PATH or provide the correct path to lammps with '-DLAMMPS_PATH=path/to/library/liblammps.so' including the file name.")
    endif()
    
    add_library(lammps_library SHARED IMPORTED)
    set_target_properties(lammps_library PROPERTIES IMPORTED_LOCATION ${LAMMPS_PATH})
    
    if(LAMMPS_FLAGS)
        set_target_properties(lammps_library PROPERTIES SHARED_LIBRARY_FLAGS ${LAMMPS_FLAGS})
    elseif(${LAMMPS_FOUND})
        set_target_properties(lammps_library PROPERTIES SHARED_LIBRARY_FLAGS "${LAMMPS_LDFLAGS}${LAMMPS_CFLAGS}")
    else()
        message(STATUS "Couldn't find lammps using pkg-config.")
        message(STATUS "Please verify that lammps.pc can be found in PKG_CONFIG_PATH or provide lammps flags with -DLAMMPS_FLAGS=\"-l<libraries>-L<link_directories>-D<VARIABLE>-I<include_directories>\" ")
    endif()

    #set_target_properties(lammps_library PROPERTIES SHARED_LIBRARY_FLAGS "${LAMMPS_LDFLAGS}${LAMMPS_CFLAGS}")
    #set_target_properties(lammps_library PROPERTIES SHARED_LIBRARY_FLAGS ${LAMMPS_FLAGS})
    #target_link_libraries(lammps_library INTERFACE lammps)
    #target_link_directories(lammps_library INTERFACE /home/catB/jd270899/.local/lib64)
    link_libraries( lammps_library )

endfunction(import_lammps)
