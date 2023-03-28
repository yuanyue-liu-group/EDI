include ../../make.inc

MODFLAGS= $(BASEMOD_FLAGS) \
    	  $(MOD_FLAG)../../PW/src \

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a\
    	 ../../upflib/libupf.a\
    	 ../../KS_Solvers/libks_solvers.a \
         ../../FFTXlib/libqefft.a\
    	 ../../UtilXlib/libutil.a\
    	 ../../XClib/xc_lib.a\
    	 ../../LAXlib/libqela.a\

MODULES = $(PWOBJS) $(QEMODS)

TLDEPS= pwlibs
EDICMODS = edic_mod.f90

EDICSRC = edic.f90\
          get_betavkb.f90\
          get_qeh_eps_data.f90\
          get_chi_data.f90\
          ml-rs.f90\
          mnl-ks.f90\
          mcharge-lfa.f90\
          mcharge.f90\
          ml-nc.f90\
          mnl-nc.f90\
          mnl-soc.f90\
          get_vloc_dat.f90\
          get_gind_rhoandpsi_gw.f90\
          get_wt_data.f90\
          gw_eps_read.f90\
          gw_eps_bcast.f90\
          gw_eps_init.f90\
          postprocess.f90\
  
OBJ = $(EDICMODS:.f90=.o) $(EDICSRC:.f90=.o)

EXE = edic.x

default: $(OBJ)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJ) $(MODULES) $(LIBOBJS) $(QELIBS)

test: $(OBJ)
	$(LD) $(LDFLAGS) -o $(test) $(OBJ) $(MODULES) $(LIBOBJS) $(QELIBS)

clean:
	rm -f *.o *~ *.d *.mod *.i *.x

