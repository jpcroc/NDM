
OBJ_MODULE = mab_in_ndm_module.o

OBJ_MODULE_LAMMPS=LAMMPS-wrapper.o \
	LAMMPS.o 



OBJ   = read_mab_file.o\
        tools.o\
	reaction_neb_module.o\
        einstein.o\
	calfo_atomic_forces.o\
        test_minimum.o\
        test_minimum_abf.o\
        langevin.o\
        langevin_overdamped.o\
        langevin_overdamped_zeta.o\
	reaction.o\
	reaction_neb.o\
        histo.o\
        calfo_BIN_OMEGA_GAUSSIEN.o\
        calfo_ABFee.o\
        calfo_ABFee_const_iter.o\
        calfo_ABF.o\
	data_processing.o\
	error_stat.o\
	mab.o 



include ../linux.arch

OBJ_LAMMPS= ndm_lammps.o

EXEC = libmab.a

all : $(EXEC)

ifdef LAMMPS_VERSION
$(EXEC) : $(OBJ_MODULE) $(OBJ_MODULE_LAMMPS)  $(OBJ_LAMMPS) $(OBJ)
	ar rcsv  $(EXEC)  $(OBJ_MODULE) $(OBJ_MODULE_LAMMPS)  $(OBJ_LAMMPS) $(OBJ)
else
$(EXEC) :$(OBJ_MODULE) $(OBJ) 
	ar rcsv  $(EXEC) $(OBJ_MODULE) $(OBJ)
endif



.PRECIOUS: libmab.a

.SUFFIXES: .f .f90 .F90 .o .a .cpp


# Compilation rules
%.o: %.c
	$(CC)  -c  ${EXTRA_INC} $(CFLAGS)            $< -o $@  
%.o: %.cpp
	$(C++)  -c ${EXTRA_INC} $(C++FLAGS)          $< -o $@ 

.f.o:
	$(F77) $(F77FLAGS_MAB) -c $<

.f.a:
	$(F77) $(F77FLAGS_MAB) -c $<
	ar rv $@ $%
	rm -f $%

.f90.o:
	$(F90) $(F90FLAGS_MAB) -c $<

.f90.a:
	$(F90) $(F90FLAGS_MAB) -c $<
	ar rv $@ $%
	rm -f $%
.F90.o :
	$(CPP)  $(CPPFLAGS_MAB) $< > $*.tmp.f90
#ifort style	$(F90) -c  $(F90FLAGS_PHONDY) -I${MKLINCLUDE}  -mkl  $*.tmp.f90
	$(F90) -c  ${EXTRA_INC} $(F90FLAGS_MAB)   $*.tmp.f90
	mv $*.tmp.o $*.o
.F90.a :
	$(CPP)  $(CPPFLAGS_MAB) $< > $*.tmp.f90
	$(F90) -c  ${EXTRA_INC} $(CPPFLAGS_MAB) $(F90FLAGS_MAB) $*.tmp.f90
	mv $*.tmp.o $*.o
	ar rv $@ $%
	rm -rf $%





clean:
	rm -f *.o *.mod *.tmp.f90  $(EXEC) 
