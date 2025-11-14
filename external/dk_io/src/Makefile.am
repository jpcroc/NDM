#******************************************************************************#
#                            src/Makefile.am
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
include $(top_srcdir)/aminclude_static.am

EXTRA_DIST= external/grisu3.h generate_mman.c external/fast_double_parser.h

bin_PROGRAMS= cfginfo cfgconvert cfgdiff
lib_LIBRARIES= libdk_io.a

cfginfo_SOURCES= cfginfo.F90
cfginfo.o: libdk_io.a

cfgconvert_SOURCES= cfgconvert.F90
cfgconvert.o: libdk_io.a

cfgdiff_SOURCES= cfgdiff.F90
cfgdiff.o: libdk_io.a

libdk_io_a_SOURCES= config.F90 io_wrappers.c posix_io.f90 array_io.c ryu/s2d.c \
 fast_double_parser.cc external/dtoa.c external/grisu3.c FileObject.F90        \
 FileReader.F90 FortranFileReader.f90 TextFileReader.f90 MmapFileReader.F90    \
 Token.F90 Structure.F90 structure_io.F90 StructureReader.F90                  \
 structure_io_dlpoly.F90 structure_io_lammps.F90 structure_io_cif.F90          \
 structure_io_castep.F90 structure_io_abinit.F90 structure_io_xfg.F90          \
 structure_io_xyz.F90 structure_io_xsf.F90 structure_io_gulp.F90               \
 structure_io_vasp.F90 SpaceGroup.F90


dk_fileobject.mod: FileObject.o
FileObject.o: FileObject.F90 dk_filereader.mod dk_gzipfilereader.mod           \
 dk_textfilereader.mod dk_fortranfilereader.mod dk_bzip2filereader.mod         \
 dk_zstdfilereader.mod dk_mmapfilereader.mod

dk_filereader.mod: FileReader.o

if HAVE_LIBZ
libdk_io_a_SOURCES+= zlib.f90 GzipFileReader.f90
dk_gzipfilereader.mod: GzipFileReader.o
GzipFileReader.o: GzipFileReader.f90 zlib.mod dk_filereader.mod dk_posix_io.mod
else
dk_gzipfilereader.mod:
	touch dk_gzipfilereader.mod
endif

if HAVE_LIBBZ2
libdk_io_a_SOURCES+= bzlib.f90 Bzip2FileReader.f90
dk_bzip2filereader.mod: Bzip2FileReader.o
Bzip2FileReader.o: Bzip2FileReader.f90 bzlib.mod dk_filereader.mod dk_posix_io.mod
else
dk_bzip2filereader.mod:
	touch dk_bzip2filereader.mod
endif

dk_fortranfilereader.mod: FortranFileReader.o
FortranFileReader.o: FortranFileReader.f90 dk_filereader.mod

dk_mmapfilereader.mod: mman.mod MmapFileReader.o
MmapFileReader.o: mman.mod MmapFileReader.F90 dk_filereader.mod

dk_textfilereader.mod: TextFileReader.o
TextFileReader.o: TextFileReader.f90 dk_filereader.mod dk_posix_io.mod

if HAVE_LIBZSTD
libdk_io_a_SOURCES+= zstd.f90 ZstdFileReader.f90
dk_zstdfilereader.mod: ZstdFileReader.o
ZstdFileReader.o: ZstdFileReader.f90 zstd.mod dk_filereader.mod
else
dk_zstdfilereader.mod:
	touch dk_zstdfilereader.mod
endif

bzlib.mod: bzlib.o
zlib.mod: zlib.o
zstd.mod: zstd.o

dk_posix_io.mod: posix_io.o

dk_structure.mod: Structure.o
Structure.o: Structure.F90

dk_structurereader.mod: StructureReader.o
StructureReader.o: StructureReader.F90 dk_structure_io.mod dk_structure.mod

dk_structure_io.mod: structure_io.o
structure_io.o: structure_io.F90 dk_fileobject.mod dk_structure.mod

dk_structure_io_abinit.mod: structure_io_abinit.o
structure_io_abinit.o: structure_io_abinit.F90 dk_structure_io.mod

dk_structure_io_castep.mod: structure_io_castep.o
structure_io_castep.o: structure_io_castep.F90 dk_structure_io.mod

dk_structure_io_dlpoly.mod: structure_io_dlpoly.o
structure_io_dlpoly.o: structure_io_dlpoly.F90 dk_structure_io.mod

dk_structure_io_lammps.mod: structure_io_lammps.o
structure_io_lammps.o: structure_io_lammps.F90 dk_structure_io.mod

dk_structure_io_vasp.mod: structure_io_vasp.o
structure_io_vasp.o: structure_io_vasp.F90 dk_token.mod dk_structure_io.mod

dk_structure_io_xfg.mod: structure_io_xfg.o
structure_io_xfg.o: structure_io_xfg.F90 dk_structure_io.mod

dk_structure_io_xsf.mod: structure_io_xsf.o
structure_io_xsf.o: structure_io_xsf.o dk_structure_io.mod

dk_structure_io_xyz.mod: structure_io_xyz.o
structure_io_xyz.o: structure_io_xyz.F90 dk_structure_io.mod

dk_structure_io_cif.mod: structure_io_cif.o

dk_structure_io_gulp.mod: structure_io_gulp.o

dk_token.mod: Token.o
Token.o: Token.F90

dk_spacegroup.mod: SpaceGroup.o

if HAVE_SPGLIB

libdk_io_a_SOURCES+= spglib_f08.F90
SpaceGroup.o: SpaceGroup.F90 spglib_f08.mod

structure_io_cif.o: structure_io_cif.F90 dk_structure_io.mod dk_spacegroup.mod \
 dk_token.mod

structure_io_gulp.o: structure_io_gulp.F90 dk_token.mod dk_structure_io.mod    \
 dk_spacegroup.mod

spglib_f08.mod: spglib_f08.o
spglib_f08.o: spglib_f08.F90
spglib_f08.F90: $(prefix)/include/spglib_f08.F90
	cp $(prefix)/include/spglib_f08.F90 spglib_f08.F90

else

SpaceGroup.o: SpaceGroup.F90

structure_io_cif.o: structure_io_cif.F90 dk_fileobject.mod dk_structure_io.mod dk_token.mod dk_spacegroup.mod

structure_io_gulp.o: structure_io_cif.F90 dk_fileobject.mod dk_structure_io.mod dk_token.mod dk_spacegroup.mod

endif

mman.mod: mman.f90
	$(FC) -c -o mman.o mman.f90

mman.f90: generate_mman
	./generate_mman > mman.f90

generate_mman: generate_mman.c
	$(CC) -o generate_mman generate_mman.c

# Set compiler flags
AM_FCFLAGS= @FORTRAN_COMPILEFLAGS@ $(OPENMP_FCFLAGS) $(CODE_COVERAGE_CFLAGS) \
 -I$(abs_top_builddir)/src

AM_CFLAGS = @C_COMPILEFLAGS@ $(OPENMP_FCFLAGS) $(CODE_COVERAGE_CFLAGS) \
  -I$(abs_top_builddir)/src/external

# Necessary for dtoa.cming IEEE floating point numbers on a x86 platform
AM_CFLAGS += -DIEEE_8087

AM_CXXFLAGS = @C_COMPILEFLAGS@ $(OPENMP_FCFLAGS) $(CODE_COVERAGE_CFLAGS) \
  -I$(abs_top_builddir)/src/external

# Set linker flags
AM_LDFLAGS= $(CODE_COVERAGE_LDFLAGS) \
 -L$(abs_top_builddir)/src

LDADD= -ldk_io

# Remove the .mod files when cleaning
clean-local:
	rm -f *.mod
	rm -f *.smod
	rm -f *.o
	rm -f *.oo
	rm -f compiler.inc
	rm -rf *.dSYM
	rm -f *.gcda
	rm -f external/*.gcda
	rm -f *.gcno
	rm -f external/*.gcno
	rm -f generate_mman

install-exec-local: all install-modules

install-modules: $(dk_io_moduledir)
	$(INSTALL_SCRIPT) $(abs_top_builddir)/src/*.mod $(dk_io_moduledir)

$(dk_io_moduledir):
	mkdir -p $(dk_io_moduledir)

