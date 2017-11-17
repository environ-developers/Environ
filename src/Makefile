# Makefile for Solvent (dielectric solvent effects)

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG).

ENVIRONOBJS = \
boundary.o                \
cavity.o                  \
charges.o                 \
dielectric.o              \
electrostatic.o           \
electrostatic_init.o      \
electrostatic_mod.o       \
electrostatic_types.o     \
environ_init.o            \
environ_input.o           \
environ_main.o            \
environ_makov_payne.o     \
environ_mod.o             \
environ_output.o          \
environ_types.o           \
extcharge.o               \
functions.o               \
generalized.o             \
generate_boundary.o       \
ioncc.o                   \
ions.o                    \
linearized_pb.o           \
periodic.o                \
poisson.o                 \
pressure.o                \
solvent_tddfpt.o          \
util.o

all : libenviron.a

libenviron.a: $(ENVIRONOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.o *~ *.a *.d *.i *.L *.mod *.F90

include make.depend

# DO NOT DELETE
