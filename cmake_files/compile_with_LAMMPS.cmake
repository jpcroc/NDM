###########
#
#   MODULE TO FIND LAMMPS LIBRARY
#
#   LAMMPS_ROOT
#
############

macro(determine_library_name)
    # Determine library name
    find_program(VAR_LMP_MPI lmp_mpi)   # check if lmp_mpi executable is available
    find_program(VAR_LMP lmp)
    find_program(VAR_LMP_SERIAL lmp_serial) # check if lmp_serial executable is available
    if (${VAR_LMP_MPI}) # lammps mpi is available
        set(lammps_mpi ON)
    else()
        set(lammps_mpi OFF)
    endif()
    if (${VAR_LMP} OR ${VAR_LMP_SERIAL}) # lammps serial is available
        set(lammps_serial ON)
    else()
        set(lammps_serial OFF)
    endif()
    if ("MPI" IN_LIST ${PROJECT_NAME}_PACKAGE_LIST) # project uses mpi
        set(project_mpi ON)
    else()
        set(project_mpi OFF)
    endif()
    if (${project_mpi} AND ${lammps_mpi}) # if project uses mpi and lammps_mpi is available
        set(LAMMPS_LIBRARY_NAME "lammps_mpi")
    elseif(${lammps_mpi} AND NOT ${lammps_serial}) # project doesn't use mpi but only lammps_mpi is available
        set(LAMMPS_LIBRARY_NAME "lammps_mpi")
    else()  # lammps_mpi is not available
        if (${VAR_LMP})
            set(LAMMPS_LIBRARY_NAME "lammps")
        elseif (${VAR_LMP_SERIAL})
            set(LAMMPS_LIBRARY_NAME "lammps_serial")
        endif()
    endif()
endmacro()

function (compile_with_LAMMPS)
    
    if (NOT DEFINED LAMMPS_LIBRARY_NAME)
        # Determine lammps library name from available executables and MPI
        determine_library_name()
    endif()
    
    # Check environment variables
    if (DEFINED ENV{LAMMPS_HOME})
        set(LAMMPS_HOME $ENV{LAMMPS_HOME})
    endif()
    if (DEFINED ENV{LAMMPS_EXTRA_LIBRARIES})
        set(LAMMPS_EXTRA_LIBRARIES $ENV{LAMMPS_EXTRA_LIBRARIES})
    endif()

    # Try to find lammps using FindLAMMPS.cmake file. Custom because not included in cmake 3.29
    find_package(LAMMPS)
    
    # Add custom libraries if necessary
    if (DEFINED LAMMPS_EXTRA_LIBRARIES)
        foreach(EXTRA_LIB IN ITEMS ${LAMMPS_EXTRA_LIBRARIES})
            list(APPEND LAMMPS_LIBRARIES ${EXTRA_LIB})
        endforeach()
    endif()
    
    # Send relevant variables to the parent scope
    set(${PROJECT_NAME}_LAMMPS_LIBRARIES ${LAMMPS_LIBRARIES} PARENT_SCOPE)
    set(${PROJECT_NAME}_LAMMPS_LIBRARY_DIRS ${LAMMPS_LIBRARY_DIRS} PARENT_SCOPE)
    set(${PROJECT_NAME}_LAMMPS_INCLUDE_DIR ${LAMMPS_INCLUDE_DIRS} PARENT_SCOPE)

    list(APPEND LAMMPS_DEFINITIONS "LAMMPS_VERSION")
    list(APPEND LAMMPS_DEFINITIONS "LAMMPS_LIB_MPI")
    Foreach( DEFINITION IN ITEMS ${LAMMPS_DEFINITIONS})
        list(APPEND ${PROJECT_NAME}_COMPILE_DEFINITION ${DEFINITION})
    endforeach()
    set(${PROJECT_NAME}_COMPILE_DEFINITION ${${PROJECT_NAME}_COMPILE_DEFINITION} PARENT_SCOPE)
endfunction(compile_with_LAMMPS)
