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
# PATCH script for files in TDDFPT/src
#
# Authors: Oliviero Andreussi (Department of Physics, University of North Thexas)
#          Iurii Timrov       (THEOS and NCCR-MARVEL, Ecole Polytechnique Federale de Lausanne)
#

#!/bin/bash

cd $TD_SRC

if test -e "Environ_PATCH" ; then
    echo "-- File Environ_PATCH exists in TDDFPT/src directory"
    echo "-- I guess you have already patched TDDFPT/src with Environ $(tail -1 Environ_PATCH)"
    echo "-- Please unpatch it first, or start from a clean source tree"
    echo "-- See you later..."
    echo "* ABORT"
    exit
fi

echo "* I will try to patch TDDFPT/src with Environ version $ENVIRON_VERSION ..."
echo "#Please do not remove or modify this file"                          >  Environ_PATCH
echo "#It keeps track of patched versions of the Environ addson package" >> Environ_PATCH
echo "$ENVIRON_VERSION"                                                  >> Environ_PATCH

# The following patches were generated using the command
# diff -u lr_readin.f90.original lr_readin.f90.new > tmp

cat > tmp.6.1 <<EOF
--- lr_readin.f90	2018-01-29 15:53:52.000000000 -0600
+++ lr_readin.f90	2018-01-29 15:53:58.000000000 -0600
@@ -59,18 +59,16 @@
   USE constants,           ONLY : eps4
   USE control_lr,          ONLY : lrpa
 #if defined(__ENVIRON)
-  USE environ_base,        ONLY : environ_base_init, ir_end
+  USE plugin_flags,        ONLY : use_environ
+  USE environ_output,      ONLY : set_environ_output
   USE environ_input,       ONLY : read_environ
-  USE environ_base,        ONLY : ifdtype, nfdpoint
-  USE ions_base,           ONLY : nsp, ityp, zv, tau, nat
+  USE ions_base,           ONLY : nsp, atm, ityp, zv, tau, nat
   USE cell_base,           ONLY : at, alat, omega, ibrav
-  USE solvent_tddfpt,      ONLY : solvent_initbase_tddfpt
-  USE environ_init,        ONLY : environ_initions, environ_initcell,      &
-                                  environ_clean, environ_initbase,         &
-                                  environ_initions_allocate
+  USE environ_init,        ONLY : environ_initions, environ_initcell, &
+                                  environ_clean, environ_initbase,    &
+                                  environ_initelectrons, environ_initpotential
   USE environ_main,        ONLY : calc_venviron
-  USE mp_bands,            ONLY : me_bgrp
-  USE plugin_flags,        ONLY : use_environ
+  USE mp_bands,            ONLY : intra_bgrp_comm, me_bgrp, root_bgrp_id
 #endif
   
 
@@ -82,7 +80,7 @@
   ! Fine control of beta_gamma_z file
   CHARACTER(LEN=80) :: disk_io
   ! Specify the amount of I/O activities
-  INTEGER :: ios, iunout, ierr, ipol
+  INTEGER :: ios, iunout, ierr, ipol, ir_end
   LOGICAL :: auto_rs
   CHARACTER(LEN=6) :: int_to_char
   !
@@ -439,33 +437,35 @@
      !
      ! Copied from PW/src/input.f90
      !
-     CALL read_environ( nat, nsp, assume_isolated, ibrav )
+     CALL set_environ_output("TD", ionode, ionode_id, intra_image_comm, stdout)
+     CALL read_environ( "TD", NINT(nelec), nspin, nat, nsp, atm, assume_isolated )
      !
      ! Taken from PW/src/init_run.f90
      !
      ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
-     CALL environ_initbase( dfftp%nnr )
+     CALL environ_initbase( dfftp%nr1, dfftp%nr2, dfftp%nr3, ibrav, alat, omega, at, &
+          & dfftp%nnr, ir_end, intra_bgrp_comm, me_bgrp, root_bgrp_id )
      !
      ! Taken from PW/src/electrons.f90
      !
-     CALL environ_initions( dfftp%nnr, nat, nsp, ityp, zv, tau, alat )
-     CALL environ_initcell( dfftp%nnr, dfftp%nr1, dfftp%nr2, dfftp%nr3, ibrav, omega, alat, at )
+     CALL environ_initions( dfftp%nnr, nat, nsp, ityp, zv, tau )
+     CALL environ_initcell( omega, at )
      !
      ! Compute additional unperturbed potentials due to the presence of the
      ! environment and add them to the SCF (HXC) potential.
      !
      WRITE( stdout, '(/5x,"Computing and adding the polarization potentials to HXC potential")' )
      !
-     CALL calc_venviron( .TRUE., dfftp%nnr, nspin, 0.0d0, rho%of_r, v%of_r(:,1))
+     CALL environ_initelectrons( nspin, dfftp%nnr, rho%of_r, nelec )
      !
-     ! Now, once the Environ potentials were computed we can deallocate numerous 
-     ! Environ arrays because they are not needed any longer.
+     CALL environ_initpotential( dfftp%nnr, v%of_r(:,1) )
      !
-     CALL environ_clean( .TRUE. )
+     CALL calc_venviron( .TRUE., dfftp%nnr, v%of_r(:,1) )
      !
-     ! Allocations for TDDFPT
+     ! Now, once the Environ potentials were computed we can deallocate numerous
+     ! Environ arrays because they are not needed any longer.
      !
-     CALL solvent_initbase_tddfpt(ifdtype, nfdpoint, dfftp%nnr)
+     CALL environ_clean( .TRUE. )
      !
   ENDIF
   !
EOF

cat > tmp.6.2 <<EOF
--- lr_readin.f90	2018-02-08 15:27:58.000000000 -0600
+++ lr_readin.f90	2018-02-08 15:37:21.000000000 -0600
@@ -59,17 +59,16 @@
   USE constants,           ONLY : eps4
   USE control_lr,          ONLY : lrpa
 #if defined(__ENVIRON)
-  USE environ_base,        ONLY : environ_base_init, ir_end
+  USE plugin_flags,        ONLY : use_environ
+  USE environ_output,      ONLY : set_environ_output
   USE environ_input,       ONLY : read_environ
-  USE environ_base,        ONLY : ifdtype, nfdpoint
-  USE ions_base,           ONLY : nsp, ityp, zv, tau, nat
+  USE ions_base,           ONLY : nsp, atm, ityp, zv, tau, nat
   USE cell_base,           ONLY : at, alat, omega, ibrav
-  USE solvent_tddfpt,      ONLY : solvent_initbase_tddfpt
-  USE environ_init,        ONLY : environ_initions, environ_initcell,      &
-                                  environ_clean, environ_initbase,         &
-                                  environ_initions_allocate
+  USE environ_init,        ONLY : environ_initions, environ_initcell, &
+                                  environ_clean, environ_initbase,    &
+                                  environ_initelectrons, environ_initpotential
   USE environ_main,        ONLY : calc_venviron
-  USE plugin_flags,        ONLY : use_environ
+  USE mp_bands,            ONLY : intra_bgrp_comm, me_bgrp, root_bgrp_id
 #endif
   
 
@@ -438,33 +437,35 @@
      !
      ! Copied from PW/src/input.f90
      !
-     CALL read_environ( nat, nsp, assume_isolated, ibrav )
+     CALL set_environ_output("TD", ionode, ionode_id, intra_image_comm, stdout)
+     CALL read_environ( "TD", NINT(nelec), nspin, nat, nsp, atm, assume_isolated )
      !
      ! Taken from PW/src/init_run.f90
      !
      ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p)
-     CALL environ_initbase( dfftp%nnr )
+     CALL environ_initbase( dfftp%nr1, dfftp%nr2, dfftp%nr3, ibrav, alat, omega, at, &
+          & dfftp%nnr, ir_end, intra_bgrp_comm, me_bgrp, root_bgrp_id )
      !
      ! Taken from PW/src/electrons.f90
      !
-     CALL environ_initions( dfftp%nnr, nat, nsp, ityp, zv, tau, alat )
-     CALL environ_initcell( dfftp%nnr, dfftp%nr1, dfftp%nr2, dfftp%nr3, ibrav, omega, alat, at )
+     CALL environ_initions( dfftp%nnr, nat, nsp, ityp, zv, tau )
+     CALL environ_initcell( omega, at )
      !
      ! Compute additional unperturbed potentials due to the presence of the
      ! environment and add them to the SCF (HXC) potential.
      !
      WRITE( stdout, '(/5x,"Computing and adding the polarization potentials to HXC potential")' )
      !
-     CALL calc_venviron( .TRUE., dfftp%nnr, nspin, 0.0d0, rho%of_r, v%of_r(:,1))
+     CALL environ_initelectrons( nspin, dfftp%nnr, rho%of_r, nelec )
      !
-     ! Now, once the Environ potentials were computed we can deallocate numerous 
-     ! Environ arrays because they are not needed any longer.
+     CALL environ_initpotential( dfftp%nnr, v%of_r(:,1) )
      !
-     CALL environ_clean( .TRUE. )
+     CALL calc_venviron( .TRUE., dfftp%nnr, v%of_r(:,1) )
      !
-     ! Allocations for TDDFPT
+     ! Now, once the Environ potentials were computed we can deallocate numerous
+     ! Environ arrays because they are not needed any longer.
      !
-     CALL solvent_initbase_tddfpt(ifdtype, nfdpoint, dfftp%nnr)
+     CALL environ_clean( .TRUE. )
      !
   ENDIF
   !
EOF

# BACKWARD COMPATIBILITY
# Compatible with QE-5.X QE-6.1.X
#patch -b -z PreENVIRON --ignore-whitespace -i tmp.6.1
# Compatible with QE-6.2, QE-6.2.1 and QE-GIT
patch -b -z PreENVIRON --ignore-whitespace -i tmp.6.2
# END BACKWARD COMPATIBILITY

rm tmp.6.1 tmp.6.2

echo "- DONE!"

cd $QE_DIR
