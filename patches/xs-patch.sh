#!/bin/bash
#
# Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
#
#    This file is part of Environ version 1.0
#
#    Environ 1.0 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    Environ 1.0 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more detail, either the file
#    `License' in the root directory of the present distribution, or
#    online at <http://www.gnu.org/licenses/>.
#
# PATCH script for files in XSpectra/src
#
# Authors: Francesco Nattino  (THEOS and NCCR-MARVEL, Ecole Polytechnique Federale de Lausanne)
#          Oliviero Andreussi (Department of Physics, University of North Thexas)
#

cd $XS_SRC

patch_makefile

check_src_patched
if test "$PATCHED" == 1; then 
   return
else
   message "Patching"
fi

echo "#Please do not remove or modify this file"                          >  Environ_PATCH
echo "#It keeps track of patched versions of the Environ addson package" >> Environ_PATCH
echo "$ENVIRON_VERSION"                                                  >> Environ_PATCH

# The following patches were generated using the command
# git diff XSpectra/src/read_input_and_bcast.f90 

cat > tmp.XS.1 <<EOF
diff --git a/XSpectra/src/read_input_and_bcast.f90 b/XSpectra/src/read_input_and_bcast.f90
index dd5c371..3997245 100644
--- a/XSpectra/src/read_input_and_bcast.f90
+++ b/XSpectra/src/read_input_and_bcast.f90
@@ -153,6 +153,8 @@ subroutine read_input_and_bcast(filerecon, r_paw)
 
   ENDIF
 
+  CALL plugin_arguments_bcast( ionode_id, world_comm )
+
   ! \$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$
   ! \$   Variables broadcasting
   ! \$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$
EOF

cat > tmp.XS.2 <<EOF
diff --git a/XSpectra/src/reset_k_points_and_reinit.f90 b/XSpectra/src/reset_k_points_and_reinit.f90
index 3f61ed0..8e6bbda 100644
--- a/XSpectra/src/reset_k_points_and_reinit.f90
+++ b/XSpectra/src/reset_k_points_and_reinit.f90
@@ -35,6 +35,9 @@ SUBROUTINE reset_k_points_and_reinit_nscf()
   CALL read_k_points()
   
   nkstot=nks
+
+  CALL plugin_read_input("XS")
+  CALL plugin_summary()
   
   CALL divide_et_impera( nkstot, xk, wk, isk, nks )
 

EOF

cat > tmp.XS.3 <<EOF
diff --git a/XSpectra/src/xspectra.f90 b/XSpectra/src/xspectra.f90
index c8f9685..37f2708 100644
--- a/XSpectra/src/xspectra.f90
+++ b/XSpectra/src/xspectra.f90
@@ -107,6 +107,8 @@ PROGRAM X_Spectra
 #endif
   CALL environment_start ( 'XSpectra' )
 
+  CALL plugin_arguments()
+
   CALL banner_xspectra()
 
   CALL read_input_and_bcast(filerecon, r_paw)
@@ -408,6 +410,7 @@ PROGRAM X_Spectra
 
   CALL stop_clock( calculation  )
   CALL print_clock( calculation )
+  CALL plugin_clock()
 
   WRITE (stdout, 1000)
   WRITE (stdout,'(5x,a)') '                           END JOB XSpectra'
EOF


patch -b -z PreENVIRON --ignore-whitespace -i tmp.XS.1 >/dev/null
patch -b -z PreENVIRON --ignore-whitespace -i tmp.XS.2 >/dev/null
patch -b -z PreENVIRON --ignore-whitespace -i tmp.XS.3 >/dev/null

rm tmp.XS.? 

printf " done!\n"

cd $QE_DIR
