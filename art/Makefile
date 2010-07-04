OBJ   = defs.o \
	random_art.o \
	art_in_ndm_module.o\
	lanczos.o \
	art01.o \
	initialize.o \
	write_refconfig.o \
	utils.o \
	find_saddle.o  \
	saddle_converge.o \
	read_parameters.o \
	restart.o \
	locate_energy.o \
	art_defects.o 

SW_pot =  min_converge.o calcfo_ndm.o 


include ../linux.arch

EXEC = libart.a

all : $(EXEC)

$(EXEC) :$(OBJ) $(SW_pot)
	ar rcsv  $(EXEC) $(OBJ) $(SW_pot)



.PRECIOUS: libart.a

.SUFFIXES: .f .f90 .o .a

.f.o:
	$(F77) $(F77FLAGS_ART) -c $<

.f.a:
	$(F77) $(F77FLAGS_ART) -c $<
	ar rv $@ $%
	rm -f $%

.f90.o:
	$(F90) $(F90FLAGS_ART) -c $<

.f90.a:
	$(F90) $(F90FLAGS_ART) -c $<
	ar rv $@ $%
	rm -f $%


clean:
	rm -f *.o *.mod $(EXEC) 
