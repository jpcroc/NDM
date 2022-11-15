
OBJ_MODULE= phondy_in_ndm_module.o \
	thermo_module.o derived_types_ph.o

OBJ_MODULE_LAMMPS=LAMMPS.o \
	LAMMPS-wrapper.o


OBJ   =	read_phondy_file.o\
	calfoeamtabvois_ph.o\
        neighbours_phondy.o\
        calfo_phondy.o\
	test_minimum.o \
        partial_index.o \
        partial_hessian.o \
	force_constant.o \
        dos.o \
        diago.o \
	small_boxes_q.o \
	dump_lammps.o \
	phondy01.o

#OBJ_ML = phondy_and_ml_module.o\
#					dyn_matrix_q_points.o


OBJ_ML = phondy_and_ml_module.o



#OBJ_ML = phondy_and_ml_module.o\


include ${NDMPLACE}/linux.arch

OBJ_LAMMPS= ndm_lammps.o



EXEC = libphondy.a

all : $(EXEC)

ifdef LAMMPS_VERSION
$(EXEC) : $(OBJ_MODULE) $(OBJ_MODULE_LAMMPS)  $(OBJ_LAMMPS) $(OBJ)
	ar rcsv  $(EXEC)  $(OBJ_MODULE) $(OBJ_MODULE_LAMMPS)  $(OBJ_LAMMPS) $(OBJ)
else
ifdef ML
$(EXEC) :$(OBJ_MODULE) $(OBJ_ML) $(OBJ)
	ar rcsv  $(EXEC) $(OBJ_MODULE) $(OBJ_ML) $(OBJ)
else
$(EXEC) :$(OBJ_MODULE) $(OBJ)
	ar rcsv  $(EXEC) $(OBJ_MODULE) $(OBJ)
endif
endif


.PRECIOUS: libphondy.a

.SUFFIXES: .f .f90 .F90 .o .a .c .cpp


# Compilation rules
#%.o: %.c
#	$(CC)  -c  ${EXTRA_INC} $(CFLAGS)            $< -o $@
#%.o: %.cpp
#	$(C++)  -c ${EXTRA_INC} $(C++FLAGS)          $< -o $@



.cpp.o:
	$(C++)  -c -I${LAMMPSPLACE}/src    ${EXTRA_INC} $(C++FLAGS)          $< -o $@
.cpp.a:
	$(C++)  -c ${EXTRA_INC} $(C++FLAGS)          $< -o $@
	ar rv $@ $%
	rm -f $%

.c.o:
	$(CC)  -c  ${EXTRA_INC} $(CFLAGS)            $< -o $@

.c.a:
	$(CC)  -c  ${EXTRA_INC} $(CFLAGS)            $< -o $@
	ar rv $@ $%
	rm -f $%




.f.o:
	$(F77) ${EXTRA_INC} $(F77FLAGS_PHONDY) -c $<

.f.a:
	$(F77) ${EXTRA_INC} $(F77FLAGS_PHONDY) -c $<
	ar rv $@ $%
	rm -f $%

.f90.o:
	$(F90) ${EXTRA_INC} $(F90FLAGS_PHONDY) -c $<

.f90.a:
	$(F90) ${EXTRA_INC} $(F90FLAGS_PHONDY) -c $<
	ar rv $@ $%
	rm -f $%
.F90.o :
	$(CPP)  $(CPPFLAGS_PHONDY) $< > $*.tmp.f90
#ifort style	$(F90) -c  $(F90FLAGS_PHONDY) -I${MKLINCLUDE}  -mkl  $*.tmp.f90
	$(F90) -c  ${EXTRA_INC} $(F90FLAGS_PHONDY)   $*.tmp.f90
	mv $*.tmp.o $*.o
.F90.a :
	$(CPP)  $(CPPFLAGS_PHONDY) $< > $*.tmp.f90
	$(F90) -c  ${EXTRA_INC} $(CPPFLAGS_PHONDY) $(F90FLAGS_PHONDY) $*.tmp.f90
	mv $*.tmp.o $*.o
	ar rv $@ $%
	rm -rf $%


clean:
	rm -f *.o *.mod *.tmp.f90  $(EXEC)
