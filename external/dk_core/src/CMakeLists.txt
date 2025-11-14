file(GLOB SRC_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.F90"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.c"
)

add_library(dk_core ${SRC_FILES})

target_compile_definitions(dk_core PRIVATE PACKAGE_NAME="dk_core")
target_compile_definitions(dk_core PRIVATE PACKAGE_STRING="dk_core 0.1")
target_compile_definitions(dk_core PRIVATE PACKAGE_VERSION="0.1")

#Build properties only
set_target_properties(dk_core PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/dk_core
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

#Pour transmettre l'emplacement des modules fortran:
target_include_directories(dk_core INTERFACE
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/mod/dk_core>
    $<INSTALL_INTERFACE:mod/dk_core>
)

#Plus portable qu'un if CMAKE_BUILD_TYPE=RELEASE ou DEBUG
target_compile_options(dk_core PRIVATE
    $<$<CONFIG:Debug>:
        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-512>
        $<$<COMPILE_LANGUAGE:C>:-fopenmp>
    >
    $<$<CONFIG:Release>:
        $<$<COMPILE_LANGUAGE:Fortran>:-funroll-loops -finline-limit=600 -ffree-line-length-512>
        $<$<COMPILE_LANGUAGE:C>:-funroll-loops -fopenmp>
    >
)