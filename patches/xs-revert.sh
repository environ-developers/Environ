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
# PATCH REVERT script for files in XSpectra/src
#
# Authors: Francesco Nattino  (THEOS and NCCR-MARVEL, Ecole Polytechnique Federale de Lausanne)
#          Oliviero Andreussi (Department of Physics, University of North Thexas)
#

cd $XS_SRC

if test ! -e Environ_PATCH ; then
    echo "-- File Environ_PATCH is not there"
    echo "-- I guess you never patched, so there is nothing to revert"
    echo "* ABORT"
    exit
fi

echo "* I will try to revert XSpectra/src with Environ version $ENVIRON_VERSION ..."
rm "Environ_PATCH"

if [ -e read_input_and_bcast.f90PreENVIRON ]; then
    mv -f read_input_and_bcast.f90PreENVIRON read_input_and_bcast.f90 
else
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
    patch -R --ignore-whitespace -i tmp.XS.1
    rm tmp.XS.1
fi

if [ -e reset_k_points_and_reinit.f90PreENVIRON ]; then
    mv -f reset_k_points_and_reinit.f90PreENVIRON reset_k_points_and_reinit.f90
else
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
    patch -R --ignore-whitespace -i tmp.XS.2
    rm tmp.XS.2
fi

if [ -e xspectra.f90PreENVIRON ]; then
    mv -f xspectra.f90PreENVIRON xspectra.f90
else
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
    patch -R --ignore-whitespace -i tmp.XS.3
    rm tmp.XS.3
fi

echo "* DONE!"

cd $QE_DIR
