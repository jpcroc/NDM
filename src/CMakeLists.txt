set(LIB_FILES
    config.F90
    io_wrappers.c
    posix_io.f90
    array_io.c
    ryu/s2d.c
    fast_double_parser.cc
    external/dtoa.c
    external/grisu3.c
    FileReader.F90
    zlib.f90
    GzipFileReader.f90
    TextFileReader.f90
    FortranFileReader.f90
    bzlib.f90
    Bzip2FileReader.f90
    zstd.f90
    ZstdFileReader.f90
    MmapFileReader.F90
    FileObject.F90
    Token.F90
    Structure.F90
    structure_io.F90
    StructureReader.F90
    structure_io_dlpoly.F90
    structure_io_lammps.F90
    SpaceGroup.F90
    structure_io_cif.F90
    structure_io_castep.F90
    structure_io_abinit.F90
    structure_io_xfg.F90
    structure_io_xyz.F90
    structure_io_xsf.F90
    structure_io_gulp.F90
    structure_io_vasp.F90
)

add_executable(generate_mman generate_mman.c)
set_target_properties(generate_mman PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tmp
)

add_custom_command(
    OUTPUT mman.f90
    COMMAND generate_mman > mman.f90
    DEPENDS generate_mman
    COMMENT "Exécution du générateur pour produire mman.f90"
)

add_library(dk_io ${LIB_FILES} mman.f90)
target_compile_definitions(dk_io PRIVATE PACKAGE_STRING="dk_io 0.1")
target_compile_definitions(dk_io PRIVATE PACKAGE_VERSION="0.1")
target_compile_definitions(dk_io PRIVATE IEEE_8087=1)
target_compile_definitions(dk_io PRIVATE INQUIRE_DIR=file)
target_include_directories(dk_io PRIVATE external)
target_include_directories(dk_io PRIVATE .)

#-DHAVE_LIBZ=1 -DHAVE_LIBBZ2=1 -DHAVE_LIBZSTD=1 -lzstd -lbz2 -lz
find_library(LIBZ z)
find_library(LIBBZ2 bz2)
find_library(LIBZSTD zstd)
set(EXTRA_LIBS "")

if(LIBZ)
    target_compile_definitions(dk_io PRIVATE HAVE_LIBZ=1)
    list(APPEND EXTRA_LIBS ${LIBZ})
endif()

if(LIBBZ2)
    target_compile_definitions(dk_io PRIVATE HAVE_LIBBZ2=1)
    list(APPEND EXTRA_LIBS ${LIBBZ2})
endif()

if(LIBZSTD)
    target_compile_definitions(dk_io PRIVATE HAVE_LIBZSTD=1)
    list(APPEND EXTRA_LIBS ${LIBZSTD})
endif()

target_link_libraries(dk_io PUBLIC DKCore::dk_core ${EXTRA_LIBS})

set_target_properties(dk_io PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/dk_io
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

#Pour transmettre l'emplacement des modules fortran:
target_include_directories(dk_io INTERFACE
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/mod/dk_io>
    $<INSTALL_INTERFACE:mod/dk_io>
)

#Flags de compilations additionnels:
target_compile_options(dk_io PRIVATE
    $<$<CONFIG:Debug>:
        $<$<COMPILE_LANGUAGE:C>:-fopenmp>
        $<$<COMPILE_LANGUAGE:CXX>:-fopenmp>
        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-512>
    >
    $<$<CONFIG:Release>:
        $<$<COMPILE_LANGUAGE:C>:-funroll-loops -fopenmp>
        $<$<COMPILE_LANGUAGE:CXX>:-funroll-loops -fopenmp>
        $<$<COMPILE_LANGUAGE:Fortran>:-funroll-loops -finline-limit=600 -ffree-line-length-512>
    >
)

if(DK-IO_OPT_COMPILE_BIN)
    add_executable(cfginfo cfginfo.F90)
    target_link_libraries(cfginfo PRIVATE dk_io)
    target_compile_definitions(cfginfo PRIVATE PACKAGE_VERSION="0.1")
    set_target_properties(cfginfo PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dk_io
    )

    add_executable(cfgconvert cfgconvert.F90)
    target_link_libraries(cfgconvert PRIVATE dk_io)
    target_compile_definitions(cfgconvert PRIVATE PACKAGE_VERSION="0.1")
    set_target_properties(cfgconvert PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dk_io
    )

    add_executable(cfgdiff cfgdiff.F90)
    target_link_libraries(cfgdiff PRIVATE dk_io)
    target_compile_definitions(cfgdiff PRIVATE PACKAGE_VERSION="0.1")
    set_target_properties(cfgdiff PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dk_io
    )
    foreach(t cfginfo cfgconvert cfgdiff)
        target_compile_options(${t} PRIVATE
            $<$<CONFIG:Debug>:
                $<$<COMPILE_LANGUAGE:C>:-fopenmp>
                $<$<COMPILE_LANGUAGE:CXX>:-fopenmp>
                $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-512>
            >
            $<$<CONFIG:Release>:
                $<$<COMPILE_LANGUAGE:C>:-funroll-loops -fopenmp>
                $<$<COMPILE_LANGUAGE:CXX>:-funroll-loops -fopenmp>
                $<$<COMPILE_LANGUAGE:Fortran>:-funroll-loops -finline-limit=600 -ffree-line-length-512>
            >
        )
    endforeach()
endif()