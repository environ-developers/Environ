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

if test -e "Environ_PATCH" ; then
    echo "-- File Environ_PATCH exists in XSpectra/src directory"
    echo "-- I guess you have already patched XSpectra/src with Environ $(tail -1 Environ_PATCH)"
    echo "-- Please unpatch it first, or start from a clean source tree"
    echo "-- See you later..."
    echo "* ABORT"
    exit
fi

echo "* I will try to patch XSpectra/src with Environ version $ENVIRON_VERSION ..."
echo "#Please do not remove or modify this file"                          >  Environ_PATCH
echo "#It keeps track of patched versions of the Environ addson package" >> Environ_PATCH
echo "$ENVIRON_VERSION"                                                  >> Environ_PATCH

# The following patches were generated using the command
# git diff XSpectra/src/read_input_and_bcast.f90 

cat > tmp.XS.1 <<EOF
diff --git a/XSpectra/src/read_input_and_bcast.f90 b/XSpectra/src/read_input_and_bcast.f90
index 4b8efc5..3228cb5 100644
--- a/XSpectra/src/read_input_and_bcast.f90
+++ b/XSpectra/src/read_input_and_bcast.f90
@@ -155,6 +155,8 @@ subroutine read_input_and_bcast(filerecon, r_paw)
 
   ENDIF
 
+  CALL plugin_arguments_bcast( ionode_id, world_comm )
+
   ! \$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$
   ! \$   Variables broadcasting
   ! \$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$\$
EOF

cat > tmp.XS.2 <<EOF
diff --git a/XSpectra/src/reset_k_points_and_reinit.f90 b/XSpectra/src/reset_k_points_and_reinit.f90
index 08940b3..ae6035e 100644
--- a/XSpectra/src/reset_k_points_and_reinit.f90
+++ b/XSpectra/src/reset_k_points_and_reinit.f90
@@ -34,6 +34,9 @@ SUBROUTINE reset_k_points_and_reinit_nscf()
   
   CALL read_k_points()
   
+  CALL plugin_read_input("XS")
+  CALL plugin_summary()
+
   nkstot=nks
   
   CALL divide_et_impera( nkstot, xk, wk, isk, nks )
EOF

cat > tmp.XS.3 <<EOF
diff --git a/XSpectra/src/xspectra.f90 b/XSpectra/src/xspectra.f90
index 223e8bb..3a2f449 100644
--- a/XSpectra/src/xspectra.f90
+++ b/XSpectra/src/xspectra.f90
@@ -110,6 +110,8 @@ PROGRAM X_Spectra
 
   CALL banner_xspectra()
 
+  CALL plugin_arguments()
+
   CALL read_input_and_bcast(filerecon, r_paw)
 
   call write_sym_param_to_stdout()
@@ -409,6 +411,7 @@ PROGRAM X_Spectra
 
   CALL stop_clock( calculation  )
   CALL print_clock( calculation )
+  CALL plugin_clock()
 
   WRITE (stdout, 1000)
   WRITE (stdout,'(5x,a)') '                           END JOB XSpectra'
EOF


patch -b -z PreENVIRON --ignore-whitespace -i tmp.XS.1
patch -b -z PreENVIRON --ignore-whitespace -i tmp.XS.2
patch -b -z PreENVIRON --ignore-whitespace -i tmp.XS.3

rm tmp.XS.? 

echo "- DONE!"

cd $QE_DIR
