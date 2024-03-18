# no subdirectory

log_warning( "NDM gets fortran sources in ${CMAKE_CURRENT_SOURCE_DIR}" )
file( GLOB NDM_LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/[^.]*.F90" )
file( GLOB NDM_MAIN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/main.F90" )
log_warning( "NDM gets fortran sources in ${CMAKE_CURRENT_SOURCE_DIR}/ARTN" )
file( GLOB NDM_LIB_SOURCES_ART "${CMAKE_CURRENT_SOURCE_DIR}/ARTN/[^.]*.F90" )
list( REMOVE_ITEM NDM_LIB_SOURCES ${NDM_MAIN_SOURCES})  # remove main.F90 from the library set
list( APPEND NDM_LIB_SOURCES ${NDM_LIB_SOURCES_ART} )

if (${PROJECT_NAME}_WITH_LAMMPS)
  log_warning( "NDM gets fortran and c++ sources in ${CMAKE_CURRENT_SOURCE_DIR}/lammps" )
  file( GLOB NDM_LIB_SOURCES_LAMMPS "${CMAKE_CURRENT_SOURCE_DIR}/lammps/[^.]*.F90" )
  list( APPEND NDM_LIB_SOURCES ${NDM_LIB_SOURCES_LAMMPS} )
  list( APPEND NDM_LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/lammps/LAMMPS-wrapper.cpp")
  include_directories("/home/catB/jd270899/.local/include/lammps")
endif()

if (NOT "MPI" IN_LIST NDM_PACKAGE_LIST)
    file( GLOB NDM_REMOVE "${CMAKE_CURRENT_SOURCE_DIR}/ARTN/gen_init_MPI.F90" ) # parallel source to remove
    list( REMOVE_ITEM NDM_LIB_SOURCES ${NDM_REMOVE})  
endif()

set_target_properties(${MAIN_LIBRARY_NAME} PROPERTIES LINKER_LANGUAGE Fortran)  # Important when c++ is enabled for lammps, without this line c++ linker takes priority over Fortran
target_sources( ${MAIN_LIBRARY_NAME} PRIVATE ${NDM_LIB_SOURCES} )

#foreach( TMP ${NDM_MAIN_SOURCES} )
#  log_warning( "add_executable ${TMP} ndm_main.exe"}

add_executable( ndm_main.exe ${NDM_MAIN_SOURCES} )
set_target_properties(ndm_main.exe PROPERTIES LINKER_LANGUAGE Fortran)  # Important when c++ is enabled for lammps, without this line c++ linker takes priority over Fortran

#endforeach()

#add_dependencies(ndm_main.exe ${MAIN_LIBRARY_NAME})
target_link_libraries( ndm_main.exe PRIVATE ${MAIN_LIBRARY_NAME} )
