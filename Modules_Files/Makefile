# Makefile for Modules

include ../../make.inc

# location of needed modules
MODFLAGS=$(BASEMOD_FLAGS) $(MOD_FLAG)../UtilXlib \
	$(MOD_FLAG)../FFTXlib $(MOD_FLAG).

# list of modules

MODULES = \
cell_base.o \
constants.o \
control_flags.o \
fft_base.o \
invmat.o \
io_files.o \
io_global.o \
ions_base.o \
kind.o \
latgen.o \
mp_bands.o \
mp_images.o \
parameters.o \
parser.o \
paw_variables.o \
pseudo_types.o \
radial_grids.o \
random_numbers.o \
recips.o \
recvec.o \
uspp.o \
volume.o \
wrappers.o \
ws_base.o

TLDEPS=libenvutil libenvfft

all : libenvmod.a

## The following is needed only for lapack compiled from sources

dlamch.o : dlamch.f
	$(F77) $(FFLAGS_NOOPT) -c $<

libenvmod.a: 	$(MODULES) $(OBJS)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi


clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L

include make.depend
