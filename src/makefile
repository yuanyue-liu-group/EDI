
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
TESTMODS = \
  v_type.f90\
  input_parameters_calcmdefect.f90\
  wavefunctions_calcmdefect.f90\
  becmod_perturb.f90\
  gvect_perturb.f90

TESTSRC = \
  main.f90\
  initialization.f90\
  calcmdefect_noncolin.f90\
  read_perturb_file.f90\
  get_vloc_colin.f90\
  calc_m.f90\
  init_us_2_perturb.f90

OBJ = $(TESTMODS:.f90=.o) $(TESTSRC:.f90=.o)

EXE = edsc.x

default: $(OBJ)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJ) $(MODULES) $(LIBOBJS) $(QELIBS)

clean:
	rm -f *.o *~ *.d *.mod *.i

