#******************************************************************************#
#                              configure.ac
#------------------------------------------------------------------------------#
#> Build system for the dk_io library.
#>
#> Version 0.1
#
#  Written by Paul Fossati, <paul.fossati@cea.fr>
#  Copyright (c) 2024-2025 CEA
#
#------------------------------------------------------------------------------#
# Redistribution and use in source and binary forms, with or without           #
# modification, are permitted provided that the following conditions are met:  #
#                                                                              #
#     * Redistributions of source code must retain the above copyright notice, #
#       this list of conditions and the following disclaimer.                  #
#                                                                              #
#     * Redistributions in binary form must reproduce the above copyright      #
#       notice, this list of conditions and the following disclaimer in the    #
#       documentation and/or other materials provided with the distribution.   #
#                                                                              #
#     * The name of the author may not be used to endorse or promote products  #
#      derived from this software without specific prior written permission    #
#      from the author.                                                        #
#                                                                              #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"  #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE   #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE     #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR          #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF         #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS     #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)      #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE   #
# POSSIBILITY OF SUCH DAMAGE.                                                  #
#******************************************************************************#

#******************************************************************************#
# Autoconf startup
#******************************************************************************#
AC_PREREQ([2.71])
AC_INIT([dk_io],[0.1],[paul.fossati@cea.fr])
AC_CONFIG_MACRO_DIRS([m4])
AC_PREFIX_DEFAULT([$HOME])

CD_MSG_SECTION([Common programs])
AM_INIT_AUTOMAKE([color-tests subdir-objects foreign])
AC_PROG_RANLIB

#AC_CONFIG_SRCDIR([src/io_wrappers.c])
#AC_CONFIG_HEADERS([config.h])

#******************************************************************************#
# Check the Fortran compiler
#******************************************************************************#
CD_MSG_SECTION([Fortran compiler])

AC_PROG_FC

# Set current language to Fortran
AC_LANG([Fortran])

# Check for useful Fortran compiler flags
AC_OPENMP
AC_FC_FREEFORM
AC_FC_LINE_LENGTH([unlimited])

# Detect the compiler type
AC_MSG_CHECKING([the compiler vendor])
AM_CONDITIONAL([USE_IFORT], [$FC --version | grep 'ifort' > /dev/null])
AM_CONDITIONAL([USE_GFORTRAN], [$FC --version | grep 'GNU Fortran' > /dev/null])
AM_CONDITIONAL([USE_FLANG], [$FC --version | grep 'flang'  > /dev/null])
AM_CONDITIONAL([USE_IFX], [$FC --version | grep 'ifx'  > /dev/null])

# Set the compiler flags
AM_COND_IF([USE_GFORTRAN], [
  AC_MSG_RESULT([gfortran])
  AC_SUBST([FORTRAN_COMPILEFLAGS], ['-g -Ofast -mcpu=native -funroll-loops -Wimplicit-procedure -finline-limit=600'])
  AC_SUBST([C_COMPILEFLAGS], ['-g -Ofast -mcpu=native -funroll-loops'])
])
AM_COND_IF([USE_FLANG], [
  AC_MSG_RESULT([flang])
  AC_SUBST([FORTRAN_COMPILEFLAGS], ['-g -O3 -funroll-loops'])
  AC_SUBST([C_COMPILEFLAGS], ['-g -O3 -funroll-loops'])
])
AM_COND_IF([USE_IFX], [
  AC_MSG_RESULT([ifx])
dnl  AC_SUBST([FORTRAN_COMPILEFLAGS], ['-g -O3 -march=native -mtune=native -funroll-loops -traceback -check all -debug all'])
  AC_SUBST([FORTRAN_COMPILEFLAGS], ['-g -O3 -march=native -mtune=native -funroll-loops'])
  AC_SUBST([C_COMPILEFLAGS], ['-g -O3 -march=native -mtune=native -funroll-loops'])
])
AM_COND_IF([USE_IFORT], [
  AC_MSG_RESULT([ifort])
  AC_SUBST([FORTRAN_COMPILEFLAGS], ['-g -diag-disable 5142,7416 -O3 -fast -xhost'])
  AC_SUBST([C_COMPILEFLAGS], ['-g -diag-disable 5142,7416 -O3 -fast -xhost'])
])

#
AC_MSG_CHECKING([the INQUIRE argument for directories])
AC_RUN_IFELSE(
  [
    program test
        logical :: exists
        inquire(file=".", exist=exists)
        if (.not.exists) error stop
    end program
  ],
  [
    AC_DEFINE([INQUIRE_DIR], [file])
    AC_MSG_RESULT([file])
  ],
  [
    AC_RUN_IFELSE(
      [
      program test
        logical :: exists
        inquire(dir=".", exist=exists)
        if (.not.exists) error stop
      end program
      ],
      [
        AC_DEFINE([INQUIRE_DIR], [dir])
        AC_MSG_RESULT([dir])
      ],
      [
        AC_RUN_IFELSE(
          [
          program test
            logical :: exists
            inquire(directory=".", exist=exists)
            if (.not.exists) error stop
          end program
          ],
          [
            AC_DEFINE([INQUIRE_DIR], [directory])
            AC_MSG_RESULT([directory])
          ],
          [
            AC_MSG_RESULT([none])
          ]
        )
      ]
    )
  ]
)

#******************************************************************************#
# Check the C++ compiler
#******************************************************************************#
CD_MSG_SECTION([C++ compiler])

AC_PROG_CXX

#******************************************************************************#
# Check the C compiler
#******************************************************************************#
CD_MSG_SECTION([C compiler])

AC_PROG_CC


#******************************************************************************#
# Check libraries
#******************************************************************************#
CD_MSG_SECTION([Libraries])

# Set current language to C for the following tests
AC_LANG([C])

# Add the prefix to the linker flags
if test \( "x$prefix" != "xNONE" \)
then
  if test \( -d "$prefix/lib" \)
  then
    LDFLAGS="$LDFLAGS -L$prefix/lib"
  fi
  if test \( -d "$prefix/lib64" \)
  then
    LDFLAGS="$LDFLAGS -L$prefix/lib64"
  fi
fi


AC_SEARCH_LIBS([inflateEnd], [z], [
  AC_DEFINE([HAVE_LIBZ])
  AM_CONDITIONAL([HAVE_LIBZ], [true])
], [AM_CONDITIONAL([HAVE_LIBZ], [false])])
AC_SEARCH_LIBS([BZ2_bzlibVersion], [bz2], [
  AC_DEFINE([HAVE_LIBBZ2])
  AM_CONDITIONAL([HAVE_LIBBZ2], [true])
], [AM_CONDITIONAL([HAVE_LIBBZ2], [false])])
AC_SEARCH_LIBS([ZSTD_versionNumber], [zstd], [
  AC_DEFINE([HAVE_LIBZSTD])
  AM_CONDITIONAL([HAVE_LIBZSTD], [true])
], [AM_CONDITIONAL([HAVE_LIBZSTD], [false])])
AC_SEARCH_LIBS([spg_get_major_version], [symspg], [
  AC_CHECK_FILE(
    [$prefix/include/spglib_f08.F90],
    [
      AC_SUBST([SPGLIB_F08], [$prefix/include/spglib_f08.F90])
      AC_DEFINE([HAVE_SPGLIB])
      AM_CONDITIONAL([HAVE_SPGLIB], [true])
    ],
    [
      AC_MSG_WARN(["⚠️ The spglib library is installed but the spglib_f08.f90 could not be found."])
      AM_CONDITIONAL([HAVE_SPGLIB], [false])
    ]
  )
], [AM_CONDITIONAL([HAVE_SPGLIB], [false])])
AC_SEARCH_LIBS([dk_core], [dk_core])

# Set current language to Fortran for the following tests
AC_LANG([Fortran])

# Add the path to the dk_core module files to the compiler flags
comp=$(echo $FC | awk '{print $1}')
AC_MSG_CHECKING([the location of the dk_core modules])
if test \( "x$prefix" != "xNONE" \)
then
  if test \( -d "$prefix/lib/modules/${comp}/dk_core" \)
  then
    dk_core_moduledir=$prefix/lib/modules/${comp}/dk_core
    dk_io_moduledir=$prefix/lib/modules/${comp}/dk_io
  fi
  if test \( -d "$prefix/lib64/modules/${comp}/dk_core" \)
  then
    dk_core_moduledir=$prefix/lib64/modules/${comp}/dk_core
    dk_io_moduledir=$prefix/lib64/modules/${comp}/dk_io
  fi
fi
if test \( "x$dk_core_moduledir" != "x" \)
then
  FCFLAGS="$FCFLAGS -I$dk_core_moduledir"
  AC_SUBST([dk_core_moduledir], [$dk_core_moduledir])
  AC_SUBST([dk_io_moduledir], [$dk_io_moduledir])
  AC_MSG_RESULT([$dk_core_moduledir])
else
  AC_MSG_RESULT([not found])
fi
AC_MSG_CHECKING([whether we can use the dk_core modules])
AC_COMPILE_IFELSE(
  [
    program test
    use ext_character
    end program
  ],
  [
    have_dk_core=yes
    AC_MSG_RESULT([yes])
  ], [
    AC_MSG_RESULT([no])
  ]
)


#******************************************************************************#
# Tests
#******************************************************************************#
CD_MSG_SECTION([Tests])

CODE_COVERAGE_LCOV_OPTIONS="--config-file ${srcdir}/doc/lcovrc --exclude **/spglib_f08.F90"
AC_SUBST([CODE_COVERAGE_LCOV_OPTIONS])
CODE_COVERAGE_GENHTML_OPTIONS="--config-file ${srcdir}/doc/lcovrc "
AC_SUBST([CODE_COVERAGE_GENHTML_OPTIONS])
CODE_COVERAGE_OUTPUT_DIRECTORY="${srcdir}/doc/coverage"
AC_SUBST([CODE_COVERAGE_OUTPUT_DIRECTORY])

# Check whether we are in coverage mode
AX_CODE_COVERAGE

# Set up the test suite
AC_CONFIG_TESTDIR([tests])


#******************************************************************************#
# Output
#******************************************************************************#
CD_MSG_SECTION([Output])

# Write Makefiles
AC_CONFIG_FILES([
  Makefile
  src/Makefile
  doc/Makefile
  tests/Makefile
])
AC_OUTPUT
