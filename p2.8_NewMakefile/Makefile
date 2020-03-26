# MCE Makefile
#########################################
#
# Makefile - makefile for the MCE code
#
# DEN/DMN/SRMP
# CEA Saclay
#
# Copyright Mickaël Trochet
#
#########################################

SHELL = /bin/bash

#sinclude Objectlist

# Definitions

f77files = $(wildcard *.f)
f90files = $(wildcard *.F90)
cfiles   = $(wildcard *.c)
#cppfiles = $(wildcard *.cpp)
INC_SRC  = $(wildcard *.h)
OBJ_SRC   = $(f77files:%.f=%.o) $(f90files:%.F90=%.o)
#MOD_SRC = $(f90files:%.F90=%.mod)

#f77mod   = $(wildcard mod/*.f)
f90mod 	  = $(wildcard mod/*.F90)
#cmod     = $(wildcard mod/*.c)
#cppmod   = $(wildcard mod/*.cpp)
INC_MOD  = $(wildcard mod/*.h)
OBJ_MOD   = $(f90mod:mod/%.F90=%.o) 
#MOD_MOD =   = $(f90mod:mod/%.F90=%.mod) 

#f77para  = $(wildcard para/*.f)
f90para   = $(wildcard para/*.F90)
#cpara    = $(wildcard para/*.c)
#cpppara   = $(wildcard para/*.cpp)
INC_PARA  = $(wildcard para/*.h)
OBJ_PARA   = $(f90para:para/%.F90=%.o)
#MOD_PARA= $(f90para:para/%.F90=%.mod)

#f77ml  = $(wildcard ml/*.f)
f90ml = $(wildcard ml/*.F90)   
#cml    = $(wildcard ml/*.c)
#cppml   = $(wildcard ml/*.cpp)
INC_ML  = $(wildcard ml/*.h)
OBJ_ML   = $(f90ml:ml/%.F90=%.o)
#MOD_ML   = $(f90ml:ml/%.F90=%.mod) 

#f77lammps  = $(wildcard lammps/*.f)
f90lammps   = $(wildcard lammps/*.F90)
#clammps     = $(wildcard lammps/*.c)
cpplammps   = $(wildcard lammps/*.cpp)
INC_LAMMPS  = $(wildcard lammps/*.h)
OBJ_LAMMPS  = $(f90lammps:lammps/%.F90=%.o) $(cpplammps:lammps/%.cpp=%.o)
#MOD_LAMMPS  = $(f90lammps:lammps/%.F90=%.mod) $(cpplammps:lammps/%.cpp=%.mod) 



checkfiles = $(f77files) $(f90files) $(f90mod)  
OBJ = $(OBJ_SRC) $(OBJ_MOD)
INC = $(INC_SRC) $(INC_MOD)
#MOD = 
ROOT = rundm90
EXT = $@
FLAG = # BECAREFUL THIS LINE IS IMPORTANT BECAUSE WE APPEND TO IT AFTER IF MACRO ARE DEFINED 
EXE = $(ROOT)_$@    #_$(shell ./version.sh)


# For dependencies here are some variables :
MAKEDEPEND=fortdepend
# dependence file name with arch extension
DEP_FILE= .dependencies_$@ 


#EXE = $(ROOT)_$@_v$(shell ./version.sh) # if git and version.sh present

# MACRO DECLARATION SECTION

ifdef PARA
OBJ +=  $(OBJ_PARA)
checkfiles +=  $(f90para) 
INC += $(INC_PARA)
FLAG += -DPARA
endif

ifdef LAMMPS
OBJ +=  $(OBJ_LAMMPS)
checkfiles += $(f90lammps) 
INC += $(INC_LAMMPS)
FLAG += -DLAMMPS_VERSION 
endif


ifdef ML # this is ML + PARAML, serial version of ML is not considered  
OBJ +=$(OBJ_ML) 
checkfiles +=  $(f90ml) 
INC += $(INC_ML) 
FLAG += -DML -DPARAML #-DPARA 
endif

# List of all targets


help:
	@echo ''
	@echo 'make arch ARGS=1          build NDM where arch is one of:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
          for file in $$files; do head -1 $$file; done
	@echo ''
	@echo 'make clean-arch          delete object folder and executable of arch'
	@echo 'make softclean-arch      only delete executable of arch'
	@echo 'make clean-all           delete object folder and executable of all arch'
	@echo 'make softclean-all       only delete executable of all arch'
	@echo 'make tar                 create NDM.tar.gz of src dir and packages'
	@echo ''



# Create a tarball of src dir and packages

.DEFAULT:
	@test -f MAKE/Makefile.$@
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@./version.sh
	@cat version.h
	@echo \#define COMPILE_TARGET \"$@\" >> version.h
	@#fortdepend -w $(MPIFLAG) -i mpi -f  *.F90 mod/* para/*  *.f  -o .dependencies
	$(MAKEDEPEND) -w $(FLAG)  -i mpi iso_c_binding -f  $(checkfiles)  -o $(DEP_FILE)
	@chmod +x MAKE/Makefile.$@ 
	@cp MAKE/Makefile.$@ Obj_$@/Makefile
	@cp name.in Obj_$@/.
	@echo "sinclude ../$(DEP_FILE)" >> Obj_$@/Makefile 
	@#echo "sinclude ../ml/milady.depends" >> Obj_$@/Makefile 
	cd Obj_$@; \
        $(MAKE) $(MFLAGS) "OBJ = $(OBJ)" "INC = $(INC)"  "EXE = ../$(EXE)" ../$(EXE)
#	@ln -sf ./$(EXE) ./$(ROOT)_exec
	#cp $(ROOT)_* ~/bin/.



tar:
	@$(MAKE) clean-all
	@tar cvzf ../$(ROOT)_src.tar.gz ./*

clean-all:
	rm -rf Obj_* $(ROOT)_*

clean-%:
	rm -rf Obj_$(@:clean-%=%) $(ROOT)_$(@:clean-%=%)

softclean-all:
	rm -rf $(ROOT)_*

softclean-%:
	rm -rf $(ROOT)_$(@:softclean-%=%)

