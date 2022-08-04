
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
EDICMODS = \
edic_mod.f90


#  v_type.f90\
#  input_parameters_calcmdefect.f90\
#  wavefunctions_calcmdefect.f90\
#  becmod_perturb.f90\
#  gvect_perturb.f90

#EDICSRC = \
#  edic.f90\
#  initialization.f90\
#  read_perturb_file.f90\
#  calc_m.f90\
#  init_us_2_sc.f90\

EDICSRC = \
  edic.f90\
  init_us_2_sc.f90\
getepsdata.f90\
ml-rs.f90\
mnl-ks.f90\
mcharge-lfa.f90\
mcharge.f90\
ml-nc.f90\
mnl-nc.f90\
mnl-soc.f90\
read_perturb_file.f90\
get_gind_rhoandpsi_gw.f90\
gw_eps_init.f90\

#h5gw_read.f90\

#mnl-nc.f90\
get_vloc_colin.f90\
getvrsc.f90\
  #get_vloc_colin.f90\
  #calcmdefect_noncolin.f90\
#  init_us_2_perturb.f90

OBJ = $(EDICMODS:.f90=.o) $(EDICSRC:.f90=.o)

EXE = edic.x

default: $(OBJ)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJ) $(MODULES) $(LIBOBJS) $(QELIBS)

clean:
	rm -f *.o *~ *.d *.mod *.i *.x

