# Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_ENV_WGET], [

  AC_CHECK_PROG(wget, wget, wget -O)
  if test "$wget" = ""; then
    AC_CHECK_PROG(wget, curl, curl -o)
  fi
  echo setting WGET... $wget
  
  AC_SUBST(wget)
  
  ]
)
