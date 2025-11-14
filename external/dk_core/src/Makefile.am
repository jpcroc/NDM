include $(top_srcdir)/aminclude_static.am

lib_LIBRARIES= libdk_core.a

libdk_core_a_SOURCES= alloc.f90 config.F90 parameters.f90 Event.f90 ext_character.f90 \
  Exception.f90 core.F90 String.f90 Console.F90 DataColumn.f90 DataTable.f90 grisu3.c \
  test.f90 pointer.f90 pointer_c.c math.f90

EXTRA_DIST= almostEqual.inc characterWithReal.inc almostEqual_wrapper.inc

alloc.o:
dk_alloc.mod: alloc.o

Console.o: Console.F90 dk_string.mod core.mod ext_character.mod
dk_console.mod: Console.o

DataColumn.o: dk_alloc.mod dk_string.mod dk_exception.mod
dk_datacolumn.mod: DataColumn.o

DataTable.o: dk_datacolumn.mod dk_string.mod
dk_datatable.mod: DataTable.o

Exception.o: Exception.f90 dk_console.mod event_levels.mod
dk_exception.mod: Exception.o

Event.o: Event.f90 core.mod
dk_event.mod: Event.o
event_levels.mod: Event.o

String.o: String.f90
dk_string.mod: String.o

parameters.o: parameters.f90
dk_parameters.mod: parameters.o

core.o: dk_parameters.mod config.mod
core.mod: core.o

config.o: config.F90
config.mod: config.o

ext_character.o: core.mod
ext_character.mod: ext_character.o

math.o: dk_math.mod dk_parameters.mod
dk_math.mod: math.f90

dk_test.mod: test.o
test.o: ext_character.mod dk_exception.mod

# Set compiler flags
AM_FCFLAGS= @FORTRAN_COMPILEFLAGS@ $(OPENMP_FCFLAGS) $(CODE_COVERAGE_CFLAGS) \
	-I$(abs_top_builddir)/src

AM_CFLAGS= @C_COMPILEFLAGS@ $(OPENMP_FCFLAGS) $(CODE_COVERAGE_CFLAGS)

# Remove the .mod files when cleaning
clean-local:
	rm -f *.mod
	rm -f compiler.inc
	rm -f *.oo
	rm -rf *.dSYM

install-exec-local: install-modules

comp := $(shell echo $(FC) | awk '{print $1}')

install-modules: $(libdir)/modules/dk_core
	mkdir -p $(libdir)/modules/$(comp)/dk_core
	$(INSTALL_SCRIPT) $(abs_top_builddir)/src/*.mod $(libdir)/modules/$(comp)/dk_core

$(libdir)/modules/dk_core:
	mkdir -p $(libdir)/modules/$(FC)/dk_core
