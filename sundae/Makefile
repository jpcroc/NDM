OBJ =	lanczosdef.o\
		tools.o \
		sundae_module.o\
        calfo_sundae.o\
        init_sundae.o\
        read_sundae.o\
        langevin.o\
        lanczos.o\
		index_voisins.o\
		sundae.o


PREPROC = -fpp -DPARASUN

include ../linux.arch

EXEC = libsundae.a

all : $(EXEC)

$(EXEC) :$(OBJ) $(SW_pot)
	ar rcsv  $(EXEC) $(OBJ) $(SW_pot)



.PRECIOUS: libsundae.a

.SUFFIXES: .f .f90 .F90 .o .a

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
.F90.o :
	$(CPP) $(CPPFLAGS_SUNDAE) $< > $*.tmp.f90
	$(F90) -c  $(F90FLAGS_SUNDAE) -I${MKLINCLUDE}  -mkl  $*.tmp.f90
	mv $*.tmp.o $*.o
.F90.a :
	$(CPP) $(CPPFLAGS_SUNDAE) $< > $*.tmp.f90
	$(F90) -c  $(F90FLAGS_SUNDAE) $*.tmp.f90
	mv $*.tmp.o $*.o
	ar rv $@ $%
	rm -rf $%


clean:
	rm -f *.o *.mod $(EXEC) 
