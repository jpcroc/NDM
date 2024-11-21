# Compile_with_MILADY

function(compile_with_MILADY)
    #search for HDF5
    #find_package(MILADY)
    
    # If MILADY is not found, compile a local version
    if (NOT MILADY_FOUND)
        include(ExternalProject)
        set(milady_build_dir "${CMAKE_BINARY_DIR}/ml")
        ExternalProject_Add(
            milady
            PREFIX ${milady_build_dir}
            SOURCE_DIR "/home/catB/jd270899/git_rep/ml"
            CMAKE_ARGS ""
            #"-DMLD_NDM=ON"
            INSTALL_COMMAND ""
            TEST_COMMAND ""       
            )
    endif()
    
    set(ENV{PATH} "${CMAKE_BINARY_DIR}/bin:$ENV{PATH}")
    set(ENV{LD_LIBRARY_PATH} "${CMAKE_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH}")
    
    #set(${PROJECT_NAME}_MILADY_LIBRARIES ${MILADY_LIBRARIES} ${MILADY_HL_LIBRARIES} PARENT_SCOPE)
    #set(${PROJECT_NAME}_MILADY_INCLUDE_DIR ${MILADY_INCLUDE_DIRS} PARENT_SCOPE)
endfunction()
