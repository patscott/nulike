# Makefile for nulike.
#
# Author: Pat Scott
# p.scott@imperial.ac.uk
#
# This is a pretty simple program, so just change
# this makefile by hand to suit your system, or
# override the variables FF, FOPT and MODULE
# from the command line.
#
# Note that the test programs require DarkSUSY,
# and the preparatory program requires nusigma.
# Just compiling the library requires neither.
# As with the compiler options, you can modify
# the DSLIB* and/or NUSIG* variables below
# to suit your installation, either by hand here
# or via the command line when calling this makefile.

# Define fortran compiler and options: intel
#FF=ifort
#FOPT=-O2 -extend_source # -g -warn all -check all #(contributed numerical routines cause warnings)
#MODULE=module
# Define fortran compiler and options: gnu
FF=gfortran
FOPT=-O2 -ffixed-line-length-none # -g -Wall -fcheck=all #(contributed numerical routines cause warnings)
MODULE=J

# DarkSUSY location, library name and include path
#DSLIBDIR = ../darksusy-5.1.3/lib
#DSLIBINC = ../darksusy-5.1.3/include
DSLIBDIR = ../../gambit/Backends/installed/darksusy/5.1.3/lib
DSLIBINC = ../../gambit/Backends/installed/darksusy/5.1.3/include
DSLIBNAME = darksusy

# nusigma location, library name and include path
NUSIGDIR = ../nusigma-1.18-pyr/lib
NUSIGINC = ../nusigma-1.18-pyr/inc
NUSIGNAME = nusigma

# Define library-making options
RANLIB = ranlib
AR = ar
ARFLAGS = rvs
SHARFLAGS = -shared -fopenmp

##### Users should not need to modify anything below this line. #####

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
FFLAGS=$(FOPT) -I$(INC) -$(MODULE) $(BUILD) -fPIC -fopenmp

# Headers
INC_DEP_BARE = nulike.h nulike_internal.h nuprep.h nuversion.h nucommon.h nuconst.h constants.h
INC_DEP = $(patsubst %.h,$(INC)/%.h,$(INC_DEP_BARE))

# Sources/objects
F90SOURCES = cons.f90 incgam.f90 marcumq.f90

SOURCES = expbesi0.f dgamic.f lngamma.f lambertw.f lambertw2.f \
lambertwln.f flagblocks.f lnpin.f lnpiln.f lnpinsum.f lnpilnsum.f \
lnpoisint.f anglike.f angres.f bgangpdf.f bginit.f bglikeprecomp.f \
bgpredinit.f bgspec.f bounds.f dP_SdE.f dPhi_SdE.f sensinit.f \
edisp.f edispinit.f sens.f eventinit.f init.f nlike.f psf.f pval.f \
analysis_map.f sigintegrand.f signal.f specintegrand.f speclike.f \
utils.f partials.f specanglike.f specangintegrand.f specanginit.f \
tabulated_weight.f preparse_files.f partintegrand.f d1mach.f \
offctrpsf.f bias.f biasinit.f

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

all : libnulike.a libnulike.so nulike_prep nulike_test nulike_test_mssm25 nulike_test_wimp


$(BUILD)/tspack.o : $(TSPACK_FULL_SOURCES)
	cat $(TSPACK_FULL_SOURCES) > tspack.f
	$(FC) $(FFLAGS) -c tspack.f -o $(BUILD)/tspack.o
	rm tspack.f


$(BUILD)/buckley.o : $(CUBPACK)/src/buckley.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/error_handline.o : $(CUBPACK)/src/error_handline.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o : $(CUBPACK)/src/%.f90 $(BUILD)/buckley.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/divide.o : $(CUBPACK)/src/divide.f90 $(BUILD)/buckley.o $(BUILD)/internal_types.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/check.o : $(CUBPACK)/src/check.f90 $(BUILD)/buckley.o $(BUILD)/internal_types.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/rule_general.o : $(CUBPACK)/src/rule_general.f90 $(BUILD)/buckley.o $(BUILD)/internal_types.o $(BUILD)/rule_1.o $(BUILD)/rule_t2.o $(BUILD)/rule_t3.o $(BUILD)/rule_tn.o $(BUILD)/rule_c2.o $(BUILD)/rule_c3.o $(BUILD)/rule_cn.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/region_processor.o : $(CUBPACK)/src/region_processor.f90 $(BUILD)/buckley.o $(BUILD)/internal_types.o $(BUILD)/rule_general.o $(BUILD)/divide.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/global_all.o : $(CUBPACK)/src/global_all.f90 $(BUILD)/buckley.o $(BUILD)/internal_types.o $(BUILD)/region_processor.o $(BUILD)/volume.o $(BUILD)/ds_routines.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/cui.o : $(CUBPACK)/src/cui.f90 $(BUILD)/buckley.o $(BUILD)/ds_routines.o $(BUILD)/rule_general.o $(BUILD)/global_all.o $(BUILD)/internal_types.o $(BUILD)/check.o $(BUILD)/error_handling.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/cons.o : $(SRC)/cons.f90 $(INC_DEP)
	$(FC) $(FFLAGS) -c $< -o $@


$(BUILD)/partials.o : $(SRC)/partials.f $(INC_DEP) $(BUILD)/buckley.o $(BUILD)/cui.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/signal.o : $(SRC)/signal.f $(INC_DEP) $(BUILD)/buckley.o $(BUILD)/cui.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/specanglike.o : $(SRC)/specanglike.f $(INC_DEP) $(BUILD)/buckley.o $(BUILD)/cui.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/speclike.o : $(SRC)/speclike.f $(INC_DEP) $(BUILD)/buckley.o $(BUILD)/cui.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/incgam.o : $(SRC)/incgam.f90 $(INC_DEP) $(BUILD)/cons.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/marcumq.o : $(SRC)/marcumq.f90 $(INC_DEP) $(BUILD)/incgam.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/offctrpsf.o : $(SRC)/offctrpsf.f $(INC_DEP) $(BUILD)/marcumq.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/partintegrand.o : $(SRC)/partintegrand.f $(INC_DEP) $(BUILD)/marcumq.o
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o : $(SRC)/%.f $(INC_DEP)
	$(FC) $(FFLAGS) -c $< -o $@


libnulike.a : $(CUBPACK_OBJ) $(OBJ)
	$(AR) $(ARFLAGS) $(LIB)/$@ $(OBJ) $(CUBPACK_OBJ)

libnulike.so : $(CUBPACK_OBJ) $(OBJ)
	$(FC) $(SHARFLAGS) -o $(LIB)/$@ $(OBJ) $(CUBPACK_OBJ)

#Note the link order of the libraries in the following executables! -lnulike must always come before -ldarksusy,
#as the nulike versions of the contributed TSPACK routines are threadsafe, whereas DarkSUSY's are not.

nulike_prep : libnulike.a $(PROGS)/nulike_prep.f
	$(FC) $(FFLAGS) -I$(NUSIGINC) -o $@ $(PROGS)/nulike_prep.f -L$(LIB) -lnulike -lgfortran -L$(NUSIGDIR) -l$(NUSIGNAME)

nulike_test : libnulike.a $(PROGS)/nulike_test.f
	$(FF) $(FFLAGS) -I$(DSLIBINC) -o $@ $(PROGS)/nulike_test.f -L$(LIB) -lnulike -L$(DSLIBDIR) -l$(DSLIBNAME) -lHB -lFH

nulike_test_mssm25 : libnulike.a $(PROGS)/nulike_test_mssm25.f
	$(FF) $(FFLAGS) -I$(DSLIBINC) -o $@ $(PROGS)/nulike_test_mssm25.f -L$(LIB) -lnulike -L$(DSLIBDIR) -l$(DSLIBNAME) -lHB -lFH

nulike_test_wimp : libnulike.a $(PROGS)/nulike_test_wimp.f
	$(FF) $(FFLAGS) -I$(DSLIBINC) -o $@ $(PROGS)/nulike_test_wimp.f -L$(LIB) -lnulike -L$(DSLIBDIR) -l$(DSLIBNAME) -lHB -lFH

clean :
	rm -f $(BUILD)/* tspack.f Key.dat nulike_prep nulike_test nulike_test_wimp nulike_test_mssm25

distclean : clean
	rm -f $(LIB)/libnulike.a $(LIB)/libnulike.so
