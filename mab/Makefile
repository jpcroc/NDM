OBJ   = mab_in_ndm_module.o\
        read_mab_file.o\
        tools.o\
        einstein.o\
        test_minimum_abf.o\
        langevin.o\
        langevin_overdamped.o\
        langevin_overdamped_csi.o\
	reaction.o\
        histo.o\
        calfo_ABF.o\
	data_processing.o\
	error_stat.o\
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
