# LES OPTIONS DE COMPILATIONS SONT DANS LINUX.ARCH
# VOIR CE FICHIER


BIN=./


SHELL=/bin/sh

# Suffixes rule
.SUFFIXES :
.SUFFIXES : .o .f90 .F90 .f


include ./linux.arch

MODULS =


                
OBJS_MOD=   T_kind_param_m.o\
	gen_com_m.o \
	var_pot.o \
        pos_com_m.o \
        fft_com_m.o \
	SMjuli.o \
	eamerco.o \
	eamtab.o \
        contrainte.o \
        cfg.o \
	gin.o \
        neb_module.o \
        jqmod.o \
        tab_imm_m.o\
        caltabtcr.o \
        calpo_ew.o \
        anasic.o \
        modules.o 



                
OBJS_PARA=   T_kind_param_mpi_m.o\
	mod_mpi.o\
	init_voisinage.o\
	init_MPI.o 


OBJS_PARAPH=   init_MPI_phondy.o
OBJS_PARAML=   init_MPI_ml.o



OBJS_SEQ=	anapos.o \
	work_cgII.o \
        DynamicAlloc.o\
        initeloss.o \
	alloc_typ.o \
        defcdp.o \
        fcc.o \
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
        tr_fire.o\
        input_pair.o \
        sauvegarde.o \
        neigcel.o \
        analyse.o \
	 spebc_fin.o\
	gcII.o \
	gcmodII.o \
        init.o \
	init_spebc.o \
	writing.o \
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
	strain_bc.o \
	stress_bc.o \
        scalebox.o \
        calfo.o\
	calfo_decalage.o\
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
        sauveforce.o \
        calfo3c.o \
        dyn.o \
        main.o \
        calcdepla.o \
        config.o \
        readdm.o \
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
        neb.o\
        extension.o\
        tersoff_zbl.o\
        test_position_He.o\
        desinteg_insert.o\
        dislo.o\
        heat.o\
        julicel.o\
        arret_ndm.o \

ifdef PARA
 # Version Parallele
 OBJS =  ${OBJS_MOD} ${OBJS_PARA} ${OBJS_SEQ} 
 EXT=_para
else
 # Version Sequentielle
 OBJS = ${OBJS_MOD} ${OBJS_SEQ}
endif


 DEP=${OBJS}
 CLEANART=clean_art
 CLEANPHONDY=clean_phondy

ifdef ART
 # Version ART
 OBJS = ${OBJS_MOD} ${OBJS_SEQ}
 DEP=${OBJS}  ${LIBART}
 CLEANART=clean_art
 DIREXT=art
 EXT=_art
endif

ifdef PHONDY
 # Version PHONDY
ifdef PARAPH
 OBJS = ${OBJS_MOD} ${OBJS_PARAPH} ${OBJS_SEQ}
 EXT=_phondy_para 
else
 OBJS = ${OBJS_MOD} ${OBJS_SEQ} 
 EXT=_phondy
endif
 DEP=${OBJS}  ${LIBPHONDY}
 CPPFLAGS=${CPPFLAGS_PHONDY}
 CLEANPHONDY=clean_phondy
 DIREXT=phondy
endif

ifdef ML
 # Version Machine Learning
ifdef PARAML
 OBJS = ${OBJS_MOD} ${OBJS_PARAML} ${OBJS_SEQ}
 EXT=_ml_para 
else
 OBJS = ${OBJS_MOD} ${OBJS_SEQ} 
 EXT=_ml
endif
 DEP=${OBJS}  ${LIBML}
 CPPFLAGS=${CPPFLAGS_ML}
 CLEANML=clean_ml
 DIREXT=ml
endif


ifdef SUNDAE
 # Version SUNDAE
ifdef PARASUN
OBJS = ${OBJS_MOD} ${OBJS_SEQ}
 EXT=_sundae_para
else
OBJS = ${OBJS_MOD} ${OBJS_SEQ} 
 EXT=_sundae
endif
 OBJS = ${OBJS_MOD} ${OBJS_SEQ}
 DEP=${OBJS}  ${LIBSUNDAE}
 CPPFLAGS=${CPPFLAGS_SUNDAE}
 CLEANSUNDAE=clean_sundae
 DIREXT=sundae
 EXT=_sundae
endif

ifdef  MAB
 # Version MAB
 OBJS = ${OBJS_MOD} ${OBJS_SEQ}
 DEP=${OBJS}  ${LIBMAB}
 CLEANMAB=clean_mab
 DIREXT=mab
 EXT=_mab
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
LIBPHONDY=libphondy.a
LIBML=libml.a
LIBSUNDAE=libsundae.a
LIBMAB=libmab.a


${BIN}/rundm90 : ${DEP} 
	${LD} -o ${BIN}/rundm90${EXT}${EXT2}  ${F90FLAGS} ${LDFLAGS} ${OBJS} ${LIBRARY}
#	rm -rf ${DIREXT}


$(LIBART): 
	(cp -rp ../art ./ ; cd art	; make ; cp *.o *.a ../../art)
$(LIBPHONDY): 
	(cp -rp ../phondy ./ ; cd phondy	; make ; cp *.o *.a ../../phondy)
$(LIBML): 
	(cp -rp ../ml ./ ; cd ml	; make ; cp *.o *.a ../../ml)
$(LIBSUNDAE): 
	(cp -rp ../sundae ./ ; cd sundae	; make ; cp *.o *.a *.tmp.f90 ../../sundae)
$(LIBMAB): 
	(cp -rp ../mab ./ ; cd mab	; make ; cp *.o *.a ../../mab)
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

configcr.o : configcr.F90
	$(CPP) $(CPPFLAGS) configcr.F90 > configcr.tmp.f90
	$(F90) -c  -O0   configcr.tmp.f90
	mv configcr.tmp.o configcr.o


clean_art:
	rm -rf art ; cd ../art ; make clean   
clean_phondy:
	rm -rf phondy ; cd ../phondy ; make clean   
clean_ml:
	rm -rf ml ; cd ../ml ; make clean   
clean_sundae:
	rm -rf sundae ; cd ../sundae ; make clean   
clean_mab:
	rm -rf mab ; cd ../mab ; make clean   
#
