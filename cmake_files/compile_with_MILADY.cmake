# Compile_with_MILADY

function(compile_with_MILADY)

    #find_package(MILADY)
    if (NOT DEFINED ENV{MILADY_ROOT})
        message(FATAL_ERROR "Missing MILADY_ROOT environment variable")
    endif()
    
    if (EXISTS "$ENV{MILADY_ROOT}/MILADYConfig.cmake")
        find_package(MILADY CONFIG HINTS $ENV{MILADY_ROOT})
        set(MILADY_FOUND ON)
    else()
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
        set(MILADY_LIBRARIES "MILADY")
        set(MILADY_LIBRARY_DIRS "${milady_build_dir}/src/milady-build/lib")
        set(MILADY_INCLUDE_DIR "${milady_build_dir}/src/milady-build/mod" )
    endif()
    
    set(${PROJECT_NAME}_MILADY_LIBRARIES "${MILADY_LIBRARIES}" PARENT_SCOPE)
    set(${PROJECT_NAME}_MILADY_LIBRARY_DIRS "${MILADY_LIBRARY_DIRS}" PARENT_SCOPE)
    set(${PROJECT_NAME}_MILADY_INCLUDE_DIR "${MILADY_INCLUDE_DIR}" PARENT_SCOPE)

endfunction()
