include $(top_srcdir)/aminclude_static.am

EXTRA_DIST=

check_PROGRAMS=DataTable

DataTable_SOURCES= DataTable.f90
DataTable.o: $(abs_top_builddir)/src/libdk_core.a


# Include Makefiles for subdirectories

#******************************************************************************#
# Flags                                                                        #
#******************************************************************************#

# Set compiler flags
AM_FCFLAGS= @FORTRAN_COMPILEFLAGS@ $(OPENMP_FCFLAGS) $(CODE_COVERAGE_CFLAGS) \
	-I$(abs_top_builddir)/src

# Set linker flags
AM_LDFLAGS= $(CODE_COVERAGE_LDFLAGS) \
	-L$(abs_top_builddir)/src

LDADD= -ldk_core

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
