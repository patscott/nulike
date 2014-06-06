# Makefile for nulike.  
# Author: Pat Scott, patscott@physics.mcgill.ca
# This is a pretty simple program, so just change the
# makefile by hand to suit your system.
# Note that the test program requires DarkSUSY,
# and the preparatory program requires nusigma.
# Just compiling the library requires neither.

# Directories
INC=include
BUILD=build
SRC=src
LIB=lib
TEST=test
TSPACK=contrib/TSPACK

# Define fortran compiler and options
FF=ifort
FC=$(FF)
FOPT=-O -extend_source -I$(INC) -module $(BUILD) -fPIC -warn all -check all
FFLAGS=$(FOPT) -c

# DarkSUSY location, library name and include path
DSLIBDIR = ../iclike2/lib
DSLIBINC = ../iclike2/include
DSLIBNAME = darksusy

# nusigma location and library name
NUSIGDIR = ../nusigma/lib
NUSIGNAME = nusigma

# Define library-making options
RANLIB = ranlib
AR = ar
ARFLAGS = rvs

# Headers
INC_DEP = nulike.h
vpath %.h $(INC)

# Sources/objects
SOURCES = lngamma.f lambertw.f lambertw2.f lambertwln.f \
dgamic.f lnpin.f lnpiln.f lnpinsum.f lnpilnsum.f \
simpson.f lnpoisint.f anglike.f angres.f bgangpdf.f \
bginit.f bglikeprecomp.f bgpredinit.f bgspec.f bounds.f \
dP_SdE.f dPhi_SdE.f eainit.f edisp.f \
edispcheckout.f edispinit.f effarea.f \
eventinit.f init.f nlike.f psf.f pval.f analysis_map.f \
sigintegrand.f signal.f specintegrand.f speclike.f
OBJ = $(patsubst %.f,$(BUILD)/%.o,$(SOURCES)) $(BUILD)/tspack.o
TSPACK_SOURCES = ENDSLP.f SIGS.f SNHCSH.f STORE.f \
YPCOEF.f YPC1.f YPC1P.f YPC2.f YPC2P.f TSPSI.f \
INTRVL.f HVAL.f HPVAL.f TSINTL.f HPPVAL.f TSVAL1.f
TSPACK_FULL_SOURCES = $(patsubst %.f,$(TSPACK)/%.f,$(TSPACK_SOURCES))

default: libnulike.a

all : $(OBJ) libnulike.a libnulike.so nulike_prep nulike_test

$(BUILD)/%.o : $(SRC)/%.F $(INC_DEP)
	$(FC) $(FFLAGS) src/$< -o $@

$(BUILD)/%.o : $(SRC)/%.f $(INC_DEP)
	$(FC) $(FFLAGS) $< -o $@

$(BUILD)/tspack.o : $(TSPACK_FULL_SOURCES)
	cat $(TSPACK_FULL_SOURCES) > tspack.f
	$(FC) $(FFLAGS) tspack.f -o $(BUILD)/tspack.o
	rm tspack.f

libnulike.a : $(OBJ)
	$(AR) $(ARFLAGS) $(LIB)/$@ $(OBJ)

libnulike.so : $(OBJ)
	$(FC) -shared -o $(LIB)/$@ $(OBJ)

nulike_prep : libnulike.a

nulike_test : libnulike.a $(TEST)/nulike_test.f
	$(FF) $(FOPT) -I$(DSLIBINC) -o $@ $(TEST)/nulike_test.f -L$(DSLIBDIR) -static -l$(DSLIBNAME) -lHB -lFH -L$(LIB) -static -lnulike -lisospin

clean : 
	rm -f $(BUILD)/* tspack.f Key.dat nulike_test nulike_prep

distclean : clean
	rm -f $(LIB)/libnulike.a $(LIB)/libnulike.so
