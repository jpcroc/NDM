include $(top_srcdir)/aminclude_static.am

EXTRA_DIST= \
  bzip2/incomplete.bz2 bzip2/corrupted.bz2 bzip2/empty.bz2 bzip2/magic.bz2 \
  gzip/incomplete.gz   gzip/corrupted.gz   gzip/empty.gz   gzip/magic.gz   \
  zstd/incomplete.zst  zstd/corrupted.zst  zstd/empty.zst  zstd/magic.zst  \
  structure_basic_checks.inc

check_PROGRAMS= fileobject textfilereader gzipfilereader bzip2filereader   \
  fortranfilereader mmapfilereader zstdfilereader structure_castep         \
  structure_cif structure_dlpoly structure_gulp structure_lammps           \
  structure_vasp structure_xfg structure_xsf structure_xyz structure_io

structure_io_SOURCES= structure_io.f90
structure_io.o: $(abs_top_builddir)/src/libdk_io.a

fileobject_SOURCES= fileobject.f90
fileobject.o: $(abs_top_builddir)/src/libdk_io.a

textfilereader_SOURCES= textfilereader.f90
textfilereader.o: $(abs_top_builddir)/src/libdk_io.a

gzipfilereader_SOURCES= gzipfilereader.F90
gzipfilereader.o: $(abs_top_builddir)/src/libdk_io.a

bzip2filereader_SOURCES= bzip2filereader.F90
bzip2filereader.o: $(abs_top_builddir)/src/libdk_io.a

fortranfilereader_SOURCES= fortranfilereader.f90
fortranfilereader.o: $(abs_top_builddir)/src/libdk_io.a

mmapfilereader_SOURCES= mmapfilereader.f90
mmapfilereader.o: $(abs_top_builddir)/src/libdk_io.a

structure_castep_SOURCES= structure_castep.f90 structure_basic_checks.inc
structure_castep.o: $(abs_top_builddir)/src/libdk_io.a

structure_cif_SOURCES= structure_cif.f90 structure_basic_checks.inc
structure_cif.o: $(abs_top_builddir)/src/libdk_io.a

structure_dlpoly_SOURCES= structure_dlpoly.f90 structure_basic_checks.inc
structure_dlpoly.o: $(abs_top_builddir)/src/libdk_io.a

structure_lammps_SOURCES= structure_lammps.f90 structure_basic_checks.inc
structure_lammps.o: $(abs_top_builddir)/src/libdk_io.a

structure_gulp_SOURCES= structure_gulp.f90 structure_basic_checks.inc
structure_gulp.o: $(abs_top_builddir)/src/libdk_io.a

structure_vasp_SOURCES= structure_vasp.f90 structure_basic_checks.inc
structure_vasp.o: $(abs_top_builddir)/src/libdk_io.a

structure_xfg_SOURCES= structure_xfg.f90 structure_basic_checks.inc
structure_xfg.o: $(abs_top_builddir)/src/libdk_io.a

structure_xsf_SOURCES= structure_xsf.f90 structure_basic_checks.inc
structure_xsf.o: $(abs_top_builddir)/src/libdk_io.a

structure_xyz_SOURCES= structure_xyz.f90 structure_basic_checks.inc
structure_xyz.o: $(abs_top_builddir)/src/libdk_io.a

zstdfilereader_SOURCES= zstdfilereader.F90 structure_basic_checks.inc
zstdfilereader.o: $(abs_top_builddir)/src/libdk_io.a

#******************************************************************************#
# Flags                                                                        #
#******************************************************************************#

#     my_program_LIBS = ... $(CODE_COVERAGE_LIBS) ...
#     my_program_CPPFLAGS = ... $(CODE_COVERAGE_CPPFLAGS) ...
#     my_program_CFLAGS = ... $(CODE_COVERAGE_CFLAGS) ...
#     my_program_CXXFLAGS = ... $(CODE_COVERAGE_CXXFLAGS) ...
#
#     clean-local: code-coverage-clean
#     distclean-local: code-coverage-dist-clean
#
#

# Set compiler flags
AM_FCFLAGS= @FORTRAN_COMPILEFLAGS@ $(OPENMP_FCFLAGS) $(CODE_COVERAGE_CFLAGS) \
 -I$(abs_top_builddir)/src

# Set linker flags
AM_LDFLAGS= $(CODE_COVERAGE_LDFLAGS) \
 -L$(abs_top_builddir)/src

LDADD= -ldk_io -ldk_core

#clean-local: code-coverage-clean
#distclean-local: code-coverage-dist-clean

$(srcdir)/package.m4: $(top_srcdir)/configure.ac
	:;{ \
	  echo '# Signature of the current package.' && \
	  echo 'm4_define([AT_PACKAGE_NAME],      [$(PACKAGE_NAME)])' && \
	  echo 'm4_define([AT_PACKAGE_TARNAME],   [$(PACKAGE_TARNAME)])' && \
	  echo 'm4_define([AT_PACKAGE_VERSION],   [$(PACKAGE_VERSION)])' && \
	  echo 'm4_define([AT_PACKAGE_STRING],    [$(PACKAGE_STRING)])' && \
	  echo 'm4_define([AT_PACKAGE_BUGREPORT], [$(PACKAGE_BUGREPORT)])'; \
	  echo 'm4_define([AT_PACKAGE_URL],       [$(PACKAGE_URL)])'; \
	} > $@-t
	mv $@-t $@

EXTRA_DIST += testsuite.at package.m4 $(TESTSUITE)
TESTSUITE = $(srcdir)/testsuite

config.inc:
	echo "character(*), parameter :: TEST_PATH = \"$(abs_top_srcdir)/tests\"" > config.inc

check-local: atconfig $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' $(TESTSUITEFLAGS)

installcheck-local: atconfig $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' $(TESTSUITEFLAGS)


## FIXME: how to autodetect the path to autom4te?
AUTOM4TE = autom4te
AUTOTEST = $(AUTOM4TE) --language=autotest
# Always generate testsuite into the source directory, not into the
# build directory, since it must be distributed, along with
# package.m4, configure, etc.
$(TESTSUITE): $(srcdir)/testsuite.at \
              $(srcdir)/package.m4
	$(AUTOTEST) -I '$(srcdir)' -o $@.tmp $@.at
	mv $@.tmp $@

DISTCLEANFILES = atconfig

# Remove various files when cleaning
clean-local:
	test ! -f '$(TESTSUITE)' || \
	  $(SHELL) '$(TESTSUITE)' --clean
	rm -f *.mod
	rm -f *.oo
	rm -f *.fh
	rm -f *.gcno
	rm -f *.gcda
	rm -f *.gcov
	rm -rf *.dSYM
	rm -f *.g90
	rm -f bzip2filereader
	rm -f fileobject
	rm -f fortranfilereader
	rm -f mmapfilereader
	rm -f gzipfilereader
	rm -f structure_xfg
	rm -f textfilereader
