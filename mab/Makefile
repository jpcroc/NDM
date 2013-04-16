OBJ   = mab_in_ndm_module.o\
        tools.o\
        langevin.o\
	mab.o 



include ../linux.arch

EXEC = libmab.a

all : $(EXEC)

$(EXEC) :$(OBJ) 
	ar rcsv  $(EXEC) $(OBJ)



.PRECIOUS: libmab.a

.SUFFIXES: .f .f90 .o .a

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


clean:
	rm -f *.o *.mod $(EXEC) 
