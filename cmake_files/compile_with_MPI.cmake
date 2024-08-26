# Compile_with_MPI
#   Provides all necessary functions to compile with MPI
#   using the following variables:
#   
#   PROJECT_NAME : set by cmake
#   MPI_HOME : can be set to specify MPI location 

function(compile_with_MPI)
    # set ENV{MPI_HOME} to get it
    find_package(MPI)
    set(${PROJECT_NAME}_MPI_INCLUDE_DIR ${MPI_Fortran_MODULE_DIR} PARENT_SCOPE)
    foreach (LANGUAGE IN ITEMS ${${PROJECT_NAME}_LANGUAGES})
        #log_info("MPI for ${LANGUAGE}")
        #log_info("${MPI_${LANGUAGE}_LIBRARIES}")
        list(APPEND TMP_MPI_LIBRARIES ${MPI_${LANGUAGE}_LIBRARIES})
        list(APPEND TMP_MPI_INCLUDE_DIR ${MPI_${LANGUAGE}_INCLUDE_DIRS})
    endforeach()
    set(${PROJECT_NAME}_MPI_LIBRARIES ${TMP_MPI_LIBRARIES} PARENT_SCOPE)
 #   set(${PROJECT_NAME}_MPI_LIBRARIES ${MPI_C_LIBRARIES} PARENT_SCOPE)
    set(${PROJECT_NAME}_FOUND_MPI ${FOUND_MPI} PARENT_SCOPE)
    set(${PROJECT_NAME}_MPI_Fortran_MODULE_DIR ${MPI_Fortran_MODULE_DIR} PARENT_SCOPE)
    set(${PROJECT_NAME}_MPI_INCLUDE_DIR ${TMP_MPI_INCLUDE_DIR} PARENT_SCOPE)
endfunction()
