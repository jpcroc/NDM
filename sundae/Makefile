OBJ   = lanczosdef.o\
	tools.o \
	sundae_module.o\
	sundae.o
#	lanczos.o


PREPROC=-fpp -DV2

include ../linux.arch

EXEC = libsundae.a

all : $(EXEC)

$(EXEC) :$(OBJ) $(SW_pot)
	ar rcsv  $(EXEC) $(OBJ) $(SW_pot)



.PRECIOUS: libsundae.a

.SUFFIXES: .f .f90 .o .a

.f.o:
	$(F77) $(F77FLAGS_SUNDAE) -c $<

.f.a:
	$(F77) $(F77FLAGS_SUNDAE) -c $<
	ar rv $@ $%
	rm -f $%

.f90.o:
	$(F90) $(F90FLAGS_SUNDAE) $(PREPROC) -c $<

.f90.a:
	$(F90) $(F90FLAGS_SUNDAE) -c $<
	ar rv $@ $%
	rm -f $%


clean:
	rm -f *.o *.mod $(EXEC) 
