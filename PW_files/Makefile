# Makefile for PW

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) $(MOD_FLAG)../UtilXlib \
	$(MOD_FLAG)../FFTXlib $(MOD_FLAG).

PWLIBS = \
esm.o \
extfield.o \
lsda_mod.o \
martyna_tuckerman.o \
scf_mod.o \
vlocal.o \
v_of_rho.o

QEMODS=../Modules_Files/libenvmod.a \
       ../FFTXlib/libenvfft.a ../UtilXlib/libenvutil.a

TLDEPS= libenvutil libenvfft

#all : tldeps pw.x generate_vdW_kernel_table.x generate_rVV10_kernel_table.x
all : libenvpw.a

libenvpw.a : $(PWLIBS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L


include make.depend

# DO NOT DELETE
