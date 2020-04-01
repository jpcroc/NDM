OBJ   = ml_in_ndm_module.o\
	derived_types.o\
	angular_functions.o\
	kcross.o\
	def_kernels.o\
	partial_index_ml.o\
	toy_models.o\
	stats_tools.o\
        read_poscar.o\
        build_subdata.o\
	read_ml_file.o\
	neighbours.o\
        compute_g2.o\
        compute_g3.o\
        compute_afs.o\
        compute_pow_so3.o\
        compute_bispectrum_so3.o\
        compute_pow_so4.o\
        compute_bispectrum_so4.o\
        compute_mtp.o\
        compute_soap.o\
        compute_acd.o\
        compute_acd_local.o\
				compute_magnetic.o\
        descriptors.o\
        alloc_typ_ml.o\
        read_desc_data.o\
	build_kernel_matrix.o\
	print_matrix.o\
	cholesky.o\
	extrapolation.o\
	ml_lbfgs.o\
	marginal_likelihood.o\
        svd.o\
	sparsification.o\
	wesley_style.o\
        nr.o\
	snap.o\
	main_compute_descriptors.o\
  regularization.o\
  snap_weights_train.o\
  calfo_ml.o\
	ml.o



include ${HOME}/Flow/NDM/linux.arch

#include milady.depends

EXEC = libml.a

all : $(EXEC)

$(EXEC) :$(OBJ)
	ar rcsv  $(EXEC) $(OBJ)



.PRECIOUS: libml.a

.SUFFIXES: .f .f90 .F90 .o .a

.f.o:
	$(F77) ${EXTRA_INC} $(F77FLAGS_ML) -c $<

.f.a:
	$(F77) ${EXTRA_INC} $(F77FLAGS_ML) -c $<
	ar rv $@ $%
	rm -f $%

.f90.o:
	$(F90) ${EXTRA_INC} $(F90FLAGS_ML) -c $<

.f90.a:
	$(F90) ${EXTRA_INC} $(F90FLAGS_ML) -c $<
	ar rv $@ $%
	rm -f $%
.F90.o :
	$(CPP) $(CPPFLAGS_ML) $< > $*.tmp.f90
#	#echo ${EXTRA_INC}
#	#pwd
	$(F90) -c  ${EXTRA_INC} $(F90FLAGS_ML)      $*.tmp.f90
	mv $*.tmp.o $*.o
.F90.a :
	$(CPP) $(CPPFLAGS_ML) $< > $*.tmp.f90
	$(F90) -c ${EXTRA_INC}  $(F90FLAGS_ML) $*.tmp.f90
	mv $*.tmp.o $*.o
	ar rv $@ $%
	rm -rf $%


clean:
	rm -f *.o *.mod *.tmp.f90 $(EXEC)
