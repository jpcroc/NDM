\#******************************************************************************#
# Autoconf startup
#******************************************************************************#
AC_PREREQ([2.71])
AC_INIT([dk_core],[0.1],[paul.fossati@gmail.com])
AC_CONFIG_MACRO_DIRS([m4])
AC_PREFIX_DEFAULT([$HOME])

CD_MSG_SECTION([Common programs])
AM_INIT_AUTOMAKE([color-tests subdir-objects foreign])
AC_PROG_RANLIB

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
  AC_SUBST([FORTRAN_COMPILEFLAGS], ['-g -O3 -march=native -mtune=native -funroll-loops -traceback -check all -debug all'])
  AC_SUBST([FORTRAN_COMPILEFLAGS], ['-g -O3 -march=native -mtune=native -funroll-loops -traceback'])
  AC_SUBST([C_COMPILEFLAGS], ['-g -O3 -march=native -mtune=native -funroll-loops'])
])
AM_COND_IF([USE_IFORT], [
  AC_MSG_RESULT([ifort])
  AC_SUBST([FORTRAN_COMPILEFLAGS], ['-g -diag-disable 5142,7416 -O3 -fast -xhost'])
  AC_SUBST([C_COMPILEFLAGS], ['-g -diag-disable 5142,7416 -O3 -fast -xhost'])
])

AC_MSG_CHECKING([whether $FC supports real(real128) variables])
AC_COMPILE_IFELSE(
  [
    program test
        use iso_fortran_env, only: real128
        real(real128) :: x
        x = asin(0.0_real128)
        write(*,*) x
    end program
  ],
  [
    AC_DEFINE([HAS_REAL128])
    AC_MSG_RESULT([yes])
  ], [
    AC_MSG_RESULT([no])
  ]
)

AC_MSG_CHECKING([whether $FC provides the backtrace subroutine])
AC_LINK_IFELSE(
  [
    program test
        implicit none(type)
        call backtrace
    end program
  ],
  [
    AC_DEFINE([HAS_BACKTRACE], [1])
    AC_MSG_RESULT([yes])
  ], [
    AC_MSG_RESULT([no])
  ]
)


#******************************************************************************#
# Check the C compiler
#******************************************************************************#
CD_MSG_SECTION([C compiler])

AC_PROG_CC

# Set current language to Fortran
AC_LANG([C])

#******************************************************************************#
# Tests
#******************************************************************************#
CD_MSG_SECTION([Tests])

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
  tests/Makefile
])
AC_OUTPUT
