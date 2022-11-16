
# order matter
# no sub dirs yet
# add_subdirectory(???)

enable_language( Fortran )

log_warning( "NDM get fortran sources in ${CMAKE_CURRENT_SOURCE_DIR}" )

file( GLOB NDM_LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/[^.]*.F90 )
file( GLOB NDM_MAIN_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/main.F90 )
list(REMOVE_ITEM NDM_LIB_SOURCES ${MAIN_SOURCES})

list( SORT NDM_LIB_SOURCES )
list( SORT NDM_MAIN_SOURCES )

if ( NDM_CMAKE_VERBOSE STREQUAL "on")
  log_info( "NDM_LIB_SOURCES" )
  foreach( TMP ${NDM_LIB_SOURCES} )
    get_filename_component(TMP1 ${TMP} NAME )
    message(STATUS "${TMP1}" )
  endforeach()

  log_info( "NDM_MAIN_SOURCES" )
  message(STATUS "${NDM_MAIN_SOURCES}" )
endif()

add_library( libndm STATIC ${NDM_LIB_SOURCES} )
target_sources( libndm PRIVATE ${NDM_LIB_SOURCES} )

#foreach( TMP ${NDM_MAIN_SOURCES} )
#  log_warning( "add_executable ${TMP} ndm_main.exe"}

add_executable( ndm_main.exe ${NDM_MAIN_SOURCES} )
#endforeach()

add_dependencies(ndm_main.exe libndm)
target_link_libraries( ndm_main.exe PRIVATE libndm ${MLD_MKL_LIB} )


install( TARGETS libndm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
install( TARGETS ndm_main.exe DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )
