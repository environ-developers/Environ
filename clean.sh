#!/bin/bash
#First install should be 0, 0, 1, 1 for the four variables below.
revertpatch=0
removesteptwo=0
initialinstall=1
steptwo=1
####################################################
#Removing created files
if [ $removesteptwo -eq 1 ] ; then
	rm Modules/make.depend
	rm PW/src/Makefile
	rm Modules/Makefile
	rm Makefile
#Renaming old files
	mv Modules/make.depend.preEnviron2 Modules/make.depend
	mv Modules/Makefile.preEnviron2 Modules/Makefile
	mv PW/src/Makefile.preEnviron PW/src/Makefile
	mv Makefile.preEnviron Makefile
fi
if [ $revertpatch -eq 1 ] ; then
	./install/addsonpatch.sh Environ Environ/src Modules -revert
	./Environ/patches/environpatch.sh -revert
	./install/makedeps.sh
	make clean
fi
####################################################
#Reinstalling
if [ $initialinstall -eq 1 ] ; then
	./configure
	make pw
	./install/addsonpatch.sh Environ Environ/src Modules -patch
	./Environ/patches/environpatch.sh -patch
	./install/makedeps.sh
fi
if [ $steptwo -eq 1 ] ; then
#Copying Makefiles
	mv Makefile Makefile.preEnviron
	mv Modules/Makefile Modules/Makefile.preEnviron2
	mv PW/src/Makefile PW/src/Makefile.preEnviron
#Copying make.depend files
	mv Modules/make.depend Modules/make.depend.preEnviron2
#Adding corrections to Makefiles
	python Environ/main_make_corrections.py Makefile.preEnviron > Makefile
	python Environ/Modules_makefile_corrections.py Modules/Makefile.preEnviron2 > Modules/Makefile
	python Environ/PW_makefile_corrections.py PW/src/Makefile.preEnviron > PW/src/Makefile
#Adding corrections to make.depend files
	python Environ/make.depend_corrections.py Modules/make.depend.preEnviron2 Modules/Environ.inc > Modules/make.depend
	make pw
fi
