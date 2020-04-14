#Cleaning first
#Removing created files
rm Modules/make.depend
rm PW/src/Makefile
rm Modules/Makefile
rm Makefile
#Renaming old files
mv Modules/make.depend.preEnviron2 Modules/make.depend
mv Modules/Makefile.preEnviron2 Modules/Makefile
mv PW/src/Makefile.preEnviron PW/src/Makefile
mv Makefile.preEnviron Makefile
./install/addsonpatch.sh Environ Environ/src Modules -revert
./Environ/patches/environpatch.sh -revert
./install/makedeps.sh
make clean
#Reinstalling
./configure
make pw
./install/addsonpatch.sh Environ Environ/src Modules -patch
./Environ/patches/environpatch.sh -patch
./install/makedeps.sh
#Copying Makefiles
mv Makefile Makefile.preEnviron
mv Modules/Makefile Modules/Makefile.preEnviron2
mv PW/src/Makefile PW/src/Makefile.preEnviron
#Copying make.depend files
mv Modules/make.depend Modules/make.depend.preEnviron2
#mv PW/src/make.depend PW/src/make.depend.preEnviron
#Adding corrections to Makefiles
python Environ/main_make_corrections.py Makefile.preEnviron > Makefile
python Environ/Modules_makefile_corrections.py Modules/Makefile.preEnviron2 > Modules/Makefile
python Environ/PW_makefile_corrections.py PW/src/Makefile.preEnviron > PW/src/Makefile
#Adding corrections to make.depend files
python Environ/make.depend_corrections.py Modules/make.depend.preEnviron2 Modules/Environ.inc > Modules/make.depend
make pw
