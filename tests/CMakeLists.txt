# https://gitlab.com/CLIUtils/modern-cmake/-/blob/master/examples/extended-project/tests/CMakeLists.txt

# Tests need to be added as executables first
# list and foreach ...

# set( CMAKE_LIBRARY_PATH ${MKLROOT}/lib/intel64)
# include_directories( ${MLD_MPI}/include )

list( APPEND simpletests # as *.f90
  test_010_hello_world
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

foreach( atest IN LISTS simpletests )
  log_info( "test ${atest}")
  add_executable( ${atest}.exe ./${atest}.f90 )

  # cmake regex see https://cmake.org/cmake/help/latest/command/string.html
  if (${atest} MATCHES "_openmpi$")   # ending "_openmpi" -> mpirun 4
    add_test( NAME ${atest} COMMAND mpirun -np 4 ../bin/${atest}.exe )
  elseif (${atest} MATCHES "_mpi_")   # with "_mpi_" -> mpirun 2
    add_test( NAME ${atest} COMMAND mpirun -np 2 ../bin/${atest}.exe )
  else()
    add_test( NAME ${atest} COMMAND ${atest}.exe )
  endif()
  set_tests_properties(${atest} PROPERTIES TIMEOUT 30)

  target_link_libraries( ${atest}.exe PRIVATE libndm ${NDM_MKL_LIB} )
endforeach()

foreach( atest IN LISTS preprosimpletests )
  log_info( "test ${atest}")
  add_executable( ${atest}.exe ./${atest}.F90 )
  add_test( NAME ${atest} COMMAND ${atest}.exe )
  set_tests_properties(${atest} PROPERTIES TIMEOUT 30)

  target_link_libraries( ${atest}.exe PRIVATE libndm ${MLD_MKL_LIB} )
endforeach()



# TODO append integration tests crocombette (python3 ....) on other list pythontests
