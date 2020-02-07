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
OBJ_SRC   = $(f77files:%.f=%.o) $(f90files:%.F90=%.o) #$(cfiles:%.c=%.o)  $(cppfiles:%.cpp=%.o)

#f77mod   = $(wildcard mod/*.f)
f90mod 	  = $(wildcard mod/*.F90)
#cmod     = $(wildcard mod/*.c)
#cppmod   = $(wildcard mod/*.cpp)
INC_MOD  = $(wildcard mod/*.h)
OBJ_MOD   = $(f90mod:mod/%.F90=%.o) #$(f77files:mod/%.f=%.o) $(cfiles:mod/%.c=%.o)  $(cppfiles:mod/%.cpp=%.o)


#f77para  = $(wildcard para/*.f)
f90para   = $(wildcard para/*.F90)
#cpara    = $(wildcard para/*.c)
#cpppara   = $(wildcard para/*.cpp)
#INC_PARA  = $(wildcard para/*.h)
OBJ_PARA   = $(f90para:para/%.F90=%.o) #$(f77files:para/%.f=%.o) $(cfiles:para/%.c=%.o)  $(cppfiles:para/%.cpp=%.o)

#f77lammps  = $(wildcard lammps/*.f)
f90lammps   = $(wildcard lammps/*.F90)
#clammps     = $(wildcard lammps/*.c)
cpplammps   = $(wildcard lammps/*.cpp)
INC_LAMMPS  = $(wildcard lammps/*.h)
OBJ_LAMMPS  = $(f90lammps:lammps/%.F90=%.o) $(cpplammps:lammps/%.cpp=%.o) #$(cfiles:para/%.c=%.o)   $(f77files:lammps/%.f=%.o)



checkfiles = $(f77files) $(f90files) $(f90mod)  
OBJ = $(OBJ_SRC) $(OBJ_MOD)
INC = $(INC_SRC) $(INC_MOD) 
ROOT = rundm90
EXT = $@
FLAG = # BECAREFUL THIS LINE IS IMPORTANT BECAUSE WE APPEND TO IT AFTER IF MACRO ARE DEFINED 
EXE = $(ROOT)_$@    #_$(shell ./version.sh)
#EXE = $(ROOT)_$@_v$(shell ./version.sh) # if git and version.sh present

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



# List of all targets

help:
	@echo ''
	@echo 'make machine             build MCE where machine is one of:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
          for file in $$files; do head -3 $$file; done
	@echo ''
	@echo 'make clean               delete all object files'
	@echo 'make tar                 create MCE.tar.gz of src dir and packages'
	@echo ''



# Create a tarball of src dir and packages

.DEFAULT:
	@test -f MAKE/Makefile.$@
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@./version.sh
	@cat version.h
	@echo \#define COMPILE_TARGET \"$@\" >> version.h
	@#fortdepend -w ${MPIFLAG} -i mpi -o .depend -f  *.f90 *F90 
	@#fortdepend -w ${MPIFLAG} -i mpi -f  *.F90 mod/* para/*  *.f  -o .dependencies
	fortdepend -w ${FLAG} -i "mpi" -f  ${checkfiles}  -o .dependencies
	@chmod +x MAKE/Makefile.$@ 
	@cp MAKE/Makefile.$@ Obj_$@/Makefile
	@cp name.in Obj_$@/.
	@echo "sinclude ../.dependencies" >> Obj_$@/Makefile 
	cd Obj_$@; \
        $(MAKE) $(MFLAGS) "OBJ = $(OBJ)" "INC = ${INC}"  "EXE = ../$(EXE)" ../$(EXE)
	@ln -sf ./$(EXE) ./$(ROOT)_exec
	#cp $(ROOT)_* ~/bin/.


tar:
	@$(MAKE) clean-all
	@cd ..; tar cvzf p2.8/$(ROOT)_src.tar.gz p2.8/* \

clean-all:
	rm -rf Obj_*

clean-%:
	rm -rf Obj_$(@:clean-%=%)

