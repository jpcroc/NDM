# Compile_with_MILADY

function(compile_with_MILADY)

    #find_package(MILADY)
    if (DEFINED ENV{MILADY_ROOT})
        set(MLD_SOURCES "$ENV{MILADY_ROOT}")
    endif()
    
    # If MILADY is not found, compile a local version
    if (NOT MILADY_FOUND)
        include(ExternalProject)
        set(milady_build_dir "${CMAKE_BINARY_DIR}/ml")
        ExternalProject_Add(
            milady
            PREFIX ${milady_build_dir}
            SOURCE_DIR "${MLD_SOURCES}"
            CMAKE_ARGS ""
            #"-DMLD_NDM=ON"
            INSTALL_COMMAND ""
            TEST_COMMAND ""       
            )
    endif()
    
    set(${PROJECT_NAME}_MILADY_LIBRARIES "MILADY" PARENT_SCOPE)
    set(${PROJECT_NAME}_MILADY_LIBRARY_DIRS "${milady_build_dir}/src/milady-build/lib" PARENT_SCOPE)
    set(${PROJECT_NAME}_MILADY_INCLUDE_DIR "${milady_build_dir}/src/milady-build/mod" PARENT_SCOPE)

endfunction()
