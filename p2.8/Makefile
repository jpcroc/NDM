# LES OPTIONS DE COMPILATIONS SONT DANS LINUX.ARCH
# VOIR CE FICHIER


 BIN=.
#BIN=../../bin

SHELL=/bin/sh

# Suffixes rule
.SUFFIXES :
.SUFFIXES : .o .f90 .F90 .f


include ./linux.arch

MODULS =


                
OBJS_MOD=   T_kind_param_m.o\
	gen_com_m.o \
        pos_com_m.o \
        fft_com_m.o \
	SMjuli.o \
	eamtab.o \
        contrainte.o \
        neb_module.o \
        jqmod.o \
        tab_imm_m.o\
        caltabtcr.o \
        anasic.o 



                
OBJS_PARA=   T_kind_param_mpi_m.o\
	mod_mpi.o\
	init_voisinage.o\
	init_MPI.o \


OBJS_SEQ=	anapos.o \
	work_cgII.o \
        defcdp.o \
        fcc.o \
        cfg.o \
        suivinonpbc_module.o \
        neb_controle.o \
        pr2.o \
        prNose2.o \
	DynamicPME.o\
        epme.o\
        transf.o\
        moduli.o\
        param_det.o\
	aluerco.o \
	correl.o \
        calccoordo.o \
        calcangle.o \
	calfoeamcel.o \
	calfoeamtabvois.o \
	julitabvois.o \
	calpoeam.o \
	adf.o\
        endrun.o \
        controle.o \
        period.o \
	notperiod.o \
        trempe.o \
        input_pair.o \
        sauvegarde.o \
        neigcel.o \
        analyse.o \
	gcII.o \
	gcmodII.o \
        init.o \
        calcdigr.o \
        tempinstcyl.o \
	Hcyl.o \
        initspeed.o \
        layer.o \
        recips.o \
        calcdepla2.o \
        zieg2.o \
        erf.o \
        rdf.o \
        calctemp.o \
        tempinst.o \
        dmloop.o \
        scalebox.o \
        calfo.o\
        initcasca.o \
        prog.o \
        caltabi.o \
        cryst_to_car.o \
        calpo.o \
        rasmol.o \
        spline.o \
        caltabt.o \
        deftimestep.o \
        calfo2ctabvois.o \
        calfo2ccel.o \
        calfoew.o \
        sauveposition.o \
        calfo3c.o \
        dyn.o \
        main.o \
        calcdepla.o \
        config.o \
        readdm.o \
        DynamicAlloc.o\
        calfow.o\
        dmloop_vverlet.o\
        dmloop_lpr.o\
        dyn_vverlet.o\
        Mat_util.o\
        divid.o\
        vanbeest.o\
        potrep.o\
        loopforcetest.o\
        Tersoff.o\
        calfoberend.o\
        configcr.o\
        fft_inter.o\
        force_tersoff_facteurs.o\
        inputtersoff.o\
        force_tersoff.o\
        force_tersoff_cel.o\
        initcdp.o\
        creadp.o\
        jqbh.o\
        bondval.o\
        decoup3D.o \
        coord_to_cell.o\
        prtplz.o\
        tr_fire.o\
        neb.o\
        extension.o\
        tersoff_zbl.o\
        test_position_He.o\
        desinteg_insert.o\
        dislo.o\
        arret_ndm.o \

ifdef PARA
 # Version Parallele
 OBJS =  ${OBJS_MOD} ${OBJS_PARA} ${OBJS_SEQ} 
else
 # Version Sequentielle
 OBJS = ${OBJS_MOD} ${OBJS_SEQ}
endif


 DEP=${OBJS}
 CLEANART=clean_art

ifdef ART
 # Version ART
 OBJS = ${OBJS_MOD} ${OBJS_SEQ}
 DEP=${OBJS}  ${LIBART}
 CLEANART=clean_art

endif
 
.f.o :
	${F90} -c  ${F90FLAGS} $<

.f90.o :
	${F90} -c  ${F90FLAGS} $<

.F90.o :
	$(CPP) $(CPPFLAGS) $< > $*.tmp.f90
	$(F90) -c  $(F90FLAGS) $*.tmp.f90
	mv $*.tmp.o $*.o


LIBART = libart.a




${BIN}/rundm90 : ${DEP} 
	${LD} -o ${BIN}/rundm90${EXT}  ${F90FLAGS} ${LDFLAGS} ${OBJS} ${LIBRARY}
	rm -rf art


$(LIBART): 
	(cp -rp ../art ./ ; cd art	; make ; cp *.o *.a ../../art)

clean : clean_seq clean_para

clean_seq :
	@rm -f  *.o *.tmp.f90 *.mod rundm90_seq*  *.l *_genmod.f90


clean_para:
	@rm -f  *.o *.tmp.f90 *.mod rundm90_para* *.l


obj_clean:
	@rm -f  *.o *.tmp.f90 *.mod *.l



fft_inter.o : fft_inter.F90
	$(CPP) $(CPPFLAGS) fft_inter.F90 > fft_inter.tmp.f90
	$(F90) -c  -O0   fft_inter.tmp.f90
	mv fft_inter.tmp.o fft_inter.o


clean_art:
	rm -rf art ; cd ../art ; make clean   
#
