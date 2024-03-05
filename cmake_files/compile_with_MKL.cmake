# Compile_with_MKL
#   Provides all necessary functions to compile with MKL
#   using the following variables:
#   
#   PROJECT_NAME : set by cmake
#   MKL_DIR : directory of MKLConfig.cmake

function(compile_with_MKL)
    # Using custom FindMKL.cmake
    if (NOT DEFINED ENV{MKL_DIR})
    find_package(MKL)
    else()
        set(MKL_CONFIG ON)
        find_package(MKL CONFIG)
    endif()
    if(NOT MKL_FOUND)
        set(MKL_CONFIG ON)
        find_package(MKL CONFIG)
    endif()
    
    # excluding libraries which are already in CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES
    set(TMP_MKL_LIBRARIES "")
    foreach(MKL_LIBRARY IN ITEMS ${MKL_LIBRARIES})
        if (NOT ("${MKL_LIBRARY}" IN_LIST CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES))
            list(APPEND TMP_MKL_LIBRARIES ${MKL_LIBRARY})
        endif()
    endforeach()
    
    set(${PROJECT_NAME}_MKL_LIBRARIES ${TMP_MKL_LIBRARIES} PARENT_SCOPE)
    set(${PROJECT_NAME}_MKL_INCLUDE_DIR ${MKL_INCLUDE} PARENT_SCOPE)
    set(${PROJECT_NAME}_MKL_FOUND ${MKL_FOUND} PARENT_SCOPE)
    set(${PROJECT_NAME}_MKL_CFLAGS ${MKL_CFLAGS} PARENT_SCOPE)

    # already defined by MKLConfig.cmake
    if (NOT MKL_CONFIG)
        #message(STATUS "MKL::MKL is set")
        #message(STATUS "MKL:${MKL_LDFLAGS}")
        add_library(MKL::MKL INTERFACE IMPORTED GLOBAL)
        target_compile_options(MKL::MKL INTERFACE ${MKL_CFLAGS} )
        target_link_libraries(MKL::MKL INTERFACE ${MKL_LDFLAGS})
        target_link_libraries(MKL::MKL INTERFACE ${MKL_LIBRARIES})
    endif()
   
endfunction()

