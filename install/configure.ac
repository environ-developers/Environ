# Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
# Copyright (C) 2001-2021 Quantum ESPRESSO group
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

AC_INIT(ENVIRON, 3.1.1, , environ)

AC_PREREQ(2.64)
AC_CONFIG_MACRO_DIR([m4/])

# Initialize variables, filling with one comes from the environment...
X_AC_ENV_DEFAULT_ENV()

# configure for current directory by default
if test "$topdir" = "" ; then topdir="`pwd`" ; fi

# check system type (no cross-compilation for now)
AC_CANONICAL_BUILD

# Checking Architecture...
X_AC_ENV_ARCH()

# Add all needed -D options to try_dflags
try_dflags=""

# Checking archiver...
X_AC_ENV_AR()

# Checking OpenMP...
X_AC_ENV_OPENMP()

# Checking MPIF90...
X_AC_ENV_MPIF90()

# Checking CC...
X_AC_ENV_CC()

# Checking F90...
X_AC_ENV_F90()

# Checking preprocessor...
X_AC_ENV_CPP()

# Checking linker...
X_AC_ENV_LD()

# Checking CUDA...
X_AC_ENV_CUDA()

# Checking F90 rule...
X_AC_ENV_F90RULE()

AC_LANG_PUSH(Fortran 77)
F77=$f90 # use Fortran 90 actually
FFLAGS="$test_fflags"
LDFLAGS="$test_ldflags"

# Checking BLAS...
X_AC_ENV_BLAS()

# Checking LAPACK...
# X_AC_ENV_LAPACK()

# Checking for FFT...
X_AC_ENV_FFT()

# Checking for IBM MASS library...
X_AC_ENV_MASS()

# check for MPI library...
X_AC_ENV_MPI()

# Checking Quantum ESPRESSO coupling...
X_AC_ENV_QE()

# Checking for ranlib...
X_AC_ENV_RANLIB()

# Checking wget or curl...
X_AC_ENV_WGET()

# Show dflags before adding $(MANUAL_DFLAGS) and adapt to XLF (if needed)
if test "$dflags" = "" ; then dflags="$try_dflags" ; fi
echo setting DFLAGS... $try_dflags

# xlf compilers (AIX and powerpc) want comma-separated -D directives
X_AC_ENV_AIX_DFLAGS()

if test "$iflags" = "" ; then iflags="$try_iflags" ; fi
echo setting IFLAGS... $iflags

# export additional settings to generated files
AC_SUBST(dflags)
AC_SUBST(fdflags)
AC_SUBST(iflags)
AC_SUBST(ld_libs)
AC_SUBST(topdir)

AC_CONFIG_FILES(make.inc)
AC_CONFIG_FILES(configure.msg)

AC_OUTPUT

# final messages
sed '/@delete@/d' configure.msg
echo configure: success
