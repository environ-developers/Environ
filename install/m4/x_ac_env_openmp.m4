# Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_ENV_OPENMP], [

AC_ARG_ENABLE(openmp,
   [AS_HELP_STRING([--enable-openmp],
       [compile for openmp execution if possible (default: no)])],
   [if   test "$enableval" = "yes" ; then
      use_openmp=1
   else
      use_openmp=0
   fi],
   [use_openmp=0])
   
])
