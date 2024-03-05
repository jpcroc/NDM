#[[
FindMKL

Find a Math Kernel Library (MKL) implementation.

Using pkg-config.

#]]
include(FindPackageHandleStandardArgs)

find_package(PkgConfig QUIET)

if (${PKG_CONFIG_FOUND})
    # library name
    if (MKL_LINK STREQUAL "static")
        set(pc_link_name "static")
    elseif(MKL_LINK STREQUAL "dynamic")
        set(pc_link_name "dynamic")
    elseif(MKL_LINK STREQUAL "sdl")
        set(pc_link_name "sdl")
    else()
        set(pc_link_name "dynamic")
    endif()

    if (MKL_INTERFACE STREQUAL "lp64")
        set(pc_interface_name "lp64")
    elseif(MKL_INTERFACE STREQUAL "ilp64")
        set(pc_interface_name "ilp64")
    else()
        set(pc_interface_name "ilp64")
    endif()

    if (MKL_THREADING STREQUAL "sequential")
        set(pc_threading_name "seq")
    elseif(MKL_THREADING STREQUAL "intel_thread")
        set(pc_threading_name "iomp")
    elseif(MKL_THREADING STREQUAL "gnu_thread")
        set(pc_threading_name "gomp")
    elseif(MKL_THREADING STREQUAL "tbb_thread")
        set(pc_threading_name "tbb")
    else()
        set(pc_threading_name "iomp")
    endif()
    set(pc_name "mkl-${pc_link_name}-${pc_interface_name}-${pc_threading_name}")

    # searching for the package
    pkg_check_modules(PKGC_MKL ${pc_name})

    set(MKL_LIBRARIES ${PKGC_MKL_LIBRARIES})
    set(MKL_LINK_LIBRARIES ${PKGC_MKL_LINK_LIBRARIES})
    set(MKL_LIBRARY_DIRS ${PKGC_MKL_LIBRARY_DIRS})
    set(MKL_LDFLAGS ${PKGC_MKL_LDFLAGS})
    set(MKL_CFLAGS ${PKGC_MKL_CFLAGS})
    set(MKL_INCLUDE /home/catB/jd270899/softwares/oneAPI/oneapi/mkl/2023.2.0/include)
endif()

if (ENABLE_SCALAPACK)
    list(APPEND MKL_LIBRARIES mkl_scalapack_${MKL_INTERFACE} mkl_blacs_${MKL_MPI}_${MKL_INTERFACE})
    list(APPEND MKL_LDFLAGS -lmkl_scalapack_${MKL_INTERFACE} -lmkl_blacs_${MKL_MPI}_${MKL_INTERFACE})
endif()


# required variables for find_package()
find_package_handle_standard_args(MKL REQUIRED_VARS MKL_LIBRARIES MKL_LINK_LIBRARIES MKL_CFLAGS MKL_LDFLAGS MKL_INCLUDE MKL_LIBRARY_DIRS)
