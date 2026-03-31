include ../../make.inc

MODFLAGS= $(BASEMOD_FLAGS) \
    	  $(MOD_FLAG)../../PW/src \
    	  $(MOD_FLAG)../../external/wannier90/src/obj \

PWOBJS = ../../PW/src/libpw.a
W90LIB = ../../external/wannier90/libwannier.a
QEMODS = ../../Modules/libqemod.a\
    	 ../../upflib/libupf.a\
    	 ../../KS_Solvers/libks_solvers.a \
         ../../FFTXlib/src/libqefft.a\
    	 ../../UtilXlib/libutil.a\
    	 ../../XClib/xc_lib.a\
    	 ../../LAXlib/libqela.a\
    	 ../../MBD/libmbd.a\
    	 ../../dft-d3/libdftd3qe.a\

MODULES = $(PWOBJS) $(QEMODS) $(W90LIB)
PWLIBS  = $(PWOBJS) $(QEMODS) $(W90LIB)

F90FLAGS_W90 = $(F90FLAGS) -D__WANLIB

TLDEPS= pwlibs

# EDI core modules
EDIMODS = edic_mod.f90 delta_weights.f90 bz_symmetry.f90 \
          wigner_seitz.f90 wan2bloch.f90 wannierize.f90 edi_read_hr.f90

# Wannier/transport modules
W90MODS = ep_constants.f90 io_var_edi.f90 input_edi.f90 wann_common.f90 \
          global_var.f90 low_lvl_edi.f90 parallelism_edi.f90 \
          openfilepw_edi.f90 io_edi.f90 kfold_edi.f90 edi_input.f90 \
          edbloch2wan.f90 ed_coarse.f90 transport.f90

# Subroutines (matrix elements, potentials)
SUBS = get_betavkb.f90 get_vloc_onthefly.f90

OBJ = $(EDIMODS:.f90=.o) $(W90MODS:.f90=.o) $(SUBS:.f90=.o)

# Default: build both edi.x and extract_pot.x
.DEFAULT_GOAL := all

all: edi.x extract_pot.x

# Main EDI executable
pw2wan_edi.o: pw2wan_edi.f90
	$(MPIF90) $(F90FLAGS_W90) -c pw2wan_edi.f90 -o pw2wan_edi.o

edi_pw2wan.o: edi_pw2wan.f90
	$(MPIF90) $(F90FLAGS_W90) -c edi_pw2wan.f90 -o edi_pw2wan.o

edi.o: edi.f90
	$(MPIF90) $(F90FLAGS_W90) -c edi.f90 -o edi.o

edi.x: $(OBJ) pw2wan_edi.o edi_pw2wan.o edi.o
	$(LD) $(LDFLAGS) -o edi.x $(OBJ) pw2wan_edi.o edi_pw2wan.o edi.o $(MODULES) $(LIBOBJS) $(QELIBS)

# Potential extraction tool (serial)
extract_pot.o: extract_pot.f90
	$(MPIF90) $(F90FLAGS) -c extract_pot.f90 -o extract_pot.o

extract_pot.x: extract_pot.o
	$(LD) $(LDFLAGS) -o extract_pot.x extract_pot.o $(PWLIBS) $(LIBOBJS) $(QELIBS)

# Test utilities
test_sum_rule.o: test_sum_rule.f90
	$(MPIF90) $(F90FLAGS) -c test_sum_rule.f90 -o test_sum_rule.o

test_potential.o: test_potential.f90
	$(MPIF90) $(F90FLAGS) -c test_potential.f90 -o test_potential.o

edi_test_sr.x: test_sum_rule.o
	$(LD) $(LDFLAGS) -o edi_test_sr.x test_sum_rule.o $(PWLIBS) $(LIBOBJS) $(QELIBS)

edi_test_pot.x: test_potential.o
	$(LD) $(LDFLAGS) -o edi_test_pot.x test_potential.o $(PWLIBS) $(LIBOBJS) $(QELIBS)

clean:
	rm -f *.o *~ *.d *.mod *.i edi.x extract_pot.x edi_test_sr.x edi_test_pot.x
