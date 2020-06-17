#!/bin/bash
if grep -q "version_number = 6.5" ../Modules/version.f90 ; then internal=1; else internal=0 ; fi
if [ $internal -eq 0 ] ; then
	fils=$( ls ./src/*.f90 | sed 's/^\.\/src\///g' | sed 's/.f90/.o/' )
	for fil in $fils
	do
		sed "s/$fil : ..\/UtilXlib/$fil : ..\/Environ\/libs\/UtilXlib/" ../Modules/make.depend > ../Modules/make.depend.tmp
		cp ../Modules/make.depend.tmp ../Modules/make.depend
		sed "s/$fil : ..\/FFTXlib/$fil : ..\/Environ\/libs\/FFTXlib/" ../Modules/make.depend > ../Modules/make.depend.tmp
		cp ../Modules/make.depend.tmp ../Modules/make.depend
	done
	sed "s/libla libfft/libla libfft envlibfft envlibutil/" ../Makefile > ../Makefile.tmp
	cp ../Makefile.tmp ../Makefile
fi
