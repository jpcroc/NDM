# https://gitlab.com/CLIUtils/modern-cmake/-/blob/master/examples/extended-project/tests/CMakeLists.txt

# Tests need to be added as executables first
# list and foreach ...

# set( CMAKE_LIBRARY_PATH ${MKLROOT}/lib/intel64)
# include_directories( ${NDM_MPI}/include )

list( APPEND simpletests # as *.f90
  test_010_hello_world
  test_010_hello_world_openmpi
)

list( APPEND preprosimpletests # as *.F90
  test_315_prepro
)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
# Command can be a target (simple) BUT with mpirun it is a bash command
# current directory during tests seems ...build/tests !
# TODO append integration tests cosmin (python3 ....)

list( SORT simpletests )

if (NOT "MPI" IN_LIST NDM_PACKAGE_LIST)
  list( FILTER simpletests EXCLUDE REGEX "_openmpi$" )
  list( FILTER simpletests EXCLUDE REGEX "_mpi_" )
endif()

add_custom_target(build_tests)

foreach( atest IN LISTS simpletests )
  log_info( "test ${atest}")
  add_executable( ${atest}.exe EXCLUDE_FROM_ALL ./${atest}.f90 )
  add_dependencies(build_tests ${atest}.exe)
  set_target_properties(${atest}.exe PROPERTIES LINKER_LANGUAGE Fortran) # Important when c++ is enabled for lammps, without this line c++ linker takes priority over Fortran

  # cmake regex see https://cmake.org/cmake/help/latest/command/string.html
  if ( ${atest} MATCHES "_openmpi$" )   # ending "_openmpi" -> mpirun 4
    add_test( NAME ${atest} COMMAND mpirun -np 4 ../bin/${atest}.exe )
  elseif ( ${atest} MATCHES "_mpi_" )   # with "_mpi_" -> mpirun 2
    add_test( NAME ${atest} COMMAND mpirun -np 2 ../bin/${atest}.exe )
  else()
    add_test( NAME ${atest} COMMAND ${atest}.exe )
  endif()

  set_tests_properties(${atest} PROPERTIES TIMEOUT 60)  # 60s timeout

  #target_link_libraries( ${atest}.exe PRIVATE ${MAIN_LIBRARY_NAME} )
endforeach()

foreach( atest IN LISTS preprosimpletests )
  log_info( "test ${atest}")
  add_executable( ${atest}.exe EXCLUDE_FROM_ALL ./${atest}.F90 )
  add_dependencies(build_tests ${atest}.exe)
  set_target_properties(${atest}.exe PROPERTIES LINKER_LANGUAGE Fortran) # Important when c++ is enabled for lammps, without this line c++ linker takes priority over Fortran
  add_test( NAME ${atest} COMMAND ${atest}.exe )
  set_tests_properties(${atest} PROPERTIES TIMEOUT 30)

  #target_link_libraries( ${atest}.exe PRIVATE ${MAIN_LIBRARY_NAME} )
endforeach()




# TODO append integration tests crocombette (python3 ....) on other list pythontests
