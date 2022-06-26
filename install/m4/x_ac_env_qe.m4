# Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
#
AC_DEFUN([X_AC_ENV_QE], [

	AC_MSG_CHECKING([for QE])

	AC_ARG_WITH(qe,
		[AS_HELP_STRING([--with-qe], [absolute path to QE root directory])],
		[
			if test "$withval" != yes ;
		 	then
		 	 	qedir="$withval"
		 	else
				qedir="$topdir/../"
		 	fi
		 	
			if test -d $qedir; then
				AC_MSG_RESULT(found at $qedir)
			else
			    AC_MSG_ERROR([$qedir is not a valid path])
			fi

		], 
		[AC_MSG_RESULT(not coupled)]
	)

	sed "s|^PREFIX=|PREFIX=$qedir|" examples/qe/environment.in > examples/qe/environment
	sed "s|^export ESPRESSO_ROOT=|export ESPRESSO_ROOT=$qedir|" tests/environment.in > tests/environment

	AC_SUBST(qedir)
	]
)
