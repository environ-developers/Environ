# Makefile for Solvent (dielectric solvent effects)

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG).

ENVIRONOBJS = \
solvent.o                             \
solvent_tddfpt.o                      \
calc_rhopol.o                         \
calc_vpol.o                           \
cavity.o                              \
pressure.o                            \
periodic.o                            \
ioncc.o                               \
extcharge.o                           \
epsregion.o                           \
environ_mod.o                         \
environ_init.o                        \
environ_info.o                        \
environ_main.o                        \
environ_input.o                       \
environ_makov_payne.o                 \
generate_dielectric.o                 \
environ_debug.o                       \
util.o                                \
lbfgs.o

all : libenviron.a

libenviron.a: $(ENVIRONOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.o *~ *.a *.d *.i *.L *.mod *.F90

include make.depend

# DO NOT DELETE
