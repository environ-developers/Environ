# Makefile for Solvent (dielectric solvent effects)

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG).

ENVIRONOBJS = \
environ_types.o                       \
solvent.o                             \
cavity.o                              \
pressure.o                            \
epsregion.o                           \
environ_mod.o                         \
environ_init.o                        \
environ_info.o                        \
environ_main.o                        \
environ_input.o                       \
environ_makov_payne.o                 \
generate_dielectric.o                 \
environ_debug.o                       \
electrostatic.o                       \
electrostatic_mod.o                   \
boundary.o                            \
functions.o                           \
poisson.o 

all : libenviron.a

libenviron.a: $(ENVIRONOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.o *~ *.a *.d *.i *.L *.mod *.F90

include make.depend

# DO NOT DELETE
