# Makefile for nulike.  
#
# Author: Pat Scott
# p.scott@imperial.ac.uk
#
# This is a pretty simple program, so just change the
# makefile by hand to suit your system, or call
# it from another makefile, overriding the variables
# FF, FOPT and MODULE.
#
# Note that the test program requires DarkSUSY,
# and the preparatory program requires nusigma.
# Just compiling the library requires neither.

# Define fortran compiler and options: intel
#FF=ifort
#FOPT=-O2 -extend_source # -warn all -check all #(contributed numerical routines cause warnings)
#MODULE=module
# Define fortran compiler and options: gnu
FF=gfortran
FOPT=-O2 -ffixed-line-length-none # -Wall -fcheck=all #(contributed numerical routines cause warnings)
MODULE=J

# DarkSUSY location, library name and include path
DSLIBDIR = ../gambit/extras/DarkSUSY/DarkSUSY/lib
DSLIBINC = ../gambit/extras/DarkSUSY/DarkSUSY/include
DSLIBNAME = darksusy

# nusigma location and library name
NUSIGDIR = ../nusigma-1.17-pyr/lib
NUSIGINC = ../nusigma-1.17-pyr/inc
NUSIGNAME = nusigma

# Define library-making options
RANLIB = ranlib
AR = ar
ARFLAGS = rvs

# Directories
INC=include
BUILD=build
SRC=src
LIB=lib
PROGS=programs
TSPACK=contrib/TSPACK
CUBPACK=contrib/CUBPACK

# Set internal compile commands
FC=$(FF)
FFLAGS=$(FOPT) -I$(INC) -$(MODULE) $(BUILD) -fPIC

# Headers
INC_DEP = nulike.h nulike_internal.h nuprep.h nuversion.h nucommon.h nuconst.h constants.h
vpath %.h $(INC)

# Sources/objects
F90SOURCES = cons.f90 incgam.f90 marcumq.f90

SOURCES = besi0.f dgamic.f lngamma.f lambertw.f lambertw2.f \
lambertwln.f flagblocks.f lnpin.f lnpiln.f lnpinsum.f lnpilnsum.f \
lnpoisint.f anglike.f angres.f bgangpdf.f bginit.f bglikeprecomp.f \
bgpredinit.f bgspec.f bounds.f dP_SdE.f dPhi_SdE.f sensinit.f \
edisp.f edispinit.f sens.f eventinit.f init.f nlike.f psf.f pval.f \
analysis_map.f sigintegrand.f signal.f specintegrand.f speclike.f \
utils.f partials.f specanglike.f specangintegrand.f specanginit.f \
tabulated_weight.f preparse_files.f partintegrand.f d1mach.f \
offctrpsf.f

OBJ = $(BUILD)/tspack.o $(patsubst %.f90,$(BUILD)/%.o,$(F90SOURCES)) $(patsubst %.f,$(BUILD)/%.o,$(SOURCES)) 

TSPACK_SOURCES = ENDSLP.f SIGS.f SNHCSH.f STORE.f \
YPCOEF.f YPC1.f YPC1P.f YPC2.f YPC2P.f TSPSI.f \
INTRVL.f HVAL.f HPVAL.f TSINTL.f HPPVAL.f TSVAL1.f
TSPACK_FULL_SOURCES = $(patsubst %.f,$(TSPACK)/%.f,$(TSPACK_SOURCES))

LOCAL_ROOT = $(PWD)

CUBPACK_OBJ_BARE = buckley.o internal_types.o ds_routines.o divide.o \
	rule_tn.o rule_t3.o rule_t2.o rule_c2.o rule_c3.o rule_cn.o  \
	rule_1.o rule_general.o region_processor.o volume.o \
	check.o global_all.o error_handling.o cui.o
CUBPACK_OBJ = $(CUBPACK_OBJ_BARE:%.o=$(BUILD)/%.o)

default: libnulike.a

all : libnulike.a libnulike.so nulike_prep nulike_test

$(BUILD)/%.o : $(SRC)/%.f $(INC_DEP)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o : $(SRC)/%.f90 $(INC_DEP)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o : $(CUBPACK)/src/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/tspack.o : $(TSPACK_FULL_SOURCES)
	cat $(TSPACK_FULL_SOURCES) > tspack.f
	$(FC) $(FFLAGS) -c tspack.f -o $(BUILD)/tspack.o
	rm tspack.f

libnulike.a : $(CUBPACK_OBJ) $(OBJ)
	$(AR) $(ARFLAGS) $(LIB)/$@ $(OBJ) $(CUBPACK_OBJ)

libnulike.so : $(CUBPACK_OBJ) $(OBJ)
	$(FC) -shared -o $(LIB)/$@ $(OBJ) $(CUBPACK_OBJ)

nulike_prep : libnulike.a $(PROGS)/nulike_prep.f
	$(FC) $(FFLAGS) -I$(NUSIGINC) -o $@ $(PROGS)/nulike_prep.f -L$(NUSIGDIR) -l$(NUSIGNAME) -L$(LIB) -lnulike -lgfortran

nulike_test : libnulike.a $(PROGS)/nulike_test.f
	$(FF) $(FFLAGS) -I$(DSLIBINC) -o $@ $(PROGS)/nulike_test.f -L$(DSLIBDIR) -l$(DSLIBNAME) -lHB -lFH -L$(LIB) -lnulike

clean : 
	rm -f $(BUILD)/* tspack.f Key.dat nulike_test nulike_prep

distclean : clean
	rm -f $(LIB)/libnulike.a $(LIB)/libnulike.so
