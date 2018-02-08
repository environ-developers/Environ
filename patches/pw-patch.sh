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
# PATCH script for plugin files in PW/src
#
# Authors: Oliviero Andreussi (Department of Physics, University of North Thexas)
#          Francesco Nattino  (THEOS and NCCR-MARVEL, Ecole Polytechnique Federale de Lausanne)
#          Ismaila Dabo       (Department of Materials Science and Engineering, Penn State)
#

#!/bin/bash

cd $PW_SRC

if test -e "Environ_PATCH" ; then
    echo "-- File Environ_PATCH exists in PW/src directory"
    echo "-- I guess you have already patched PW/src with Environ $(tail -1 Environ_PATCH)"
    echo "-- Please unpatch it first, or start from a clean source tree"
    echo "-- See you later..."
    echo "* ABORT"
    exit
fi

echo "* I will try to patch PW/src with Environ version $ENVIRON_VERSION ..."
echo "#Please do not remove or modify this file"                          >  Environ_PATCH
echo "#It keeps track of patched versions of the Environ addson package" >> Environ_PATCH
echo "$ENVIRON_VERSION"                                                  >> Environ_PATCH

# plugin_int_forces

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
  USE environ_main,  ONLY : calc_fenviron\
!Environ patch
' plugin_int_forces.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch \
  REAL(DP), ALLOCATABLE :: force_environ(:,:)\
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
  IF (use_environ) THEN\
    !\
    ALLOCATE(force_environ(3,nat))\
    !\
    force_environ=0.0_dp\
    !\
    ! ... Add environment contributions\
    !\
    CALL calc_fenviron( nat, force_environ )\
    !\
    IF ( iverbosity > 0 ) THEN\
      WRITE( stdout, 9001 )\
      DO na = 1, nat\
         WRITE( stdout, 9002 ) na, ityp(na), ( force_environ(ipol,na), ipol = 1, 3 )\
      END DO\
      WRITE( stdout, * )\
    ENDIF\
    !\
    force = force_environ\
    !\
    DEALLOCATE(force_environ)\
    !\
  END IF\
  !\
9001 FORMAT(5x,"The global environment contribution to forces")\
9002 FORMAT(5X,"atom ",I4," type ",I2,"   force = ",3F14.8)\
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_int_forces.f90

# plugin_read_input

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE io_global,  ONLY : ionode, ionode_id, stdout\
USE mp_images,  ONLY : intra_image_comm\
USE input_parameters, ONLY : nspin\
USE input_parameters, ONLY : nat, ntyp, atom_label\
USE input_parameters, ONLY : assume_isolated, ibrav\
USE environ_input,    ONLY : read_environ\
USE environ_output,   ONLY : set_environ_output\
!Environ patch
' plugin_read_input.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   IF (use_environ) THEN\
      CALL set_environ_output("PW", ionode, ionode_id, intra_image_comm, stdout)\
      CALL read_environ("PW",1, nspin, nat, ntyp, atom_label, assume_isolated)\
   ENDIF\
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_read_input.f90

# plugin_clean

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_init, ONLY : environ_clean \
USE    control_flags, ONLY : tddfpt \
!Environ patch
' plugin_clean.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ.AND..NOT.tddfpt) CALL environ_clean(lflag) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_clean.f90

# plugin_summary

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    io_global,      ONLY : stdout \
USE    environ_output, ONLY : environ_summary, & \
                              update_output_program_unit \
!Environ patch
' plugin_summary.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL update_output_program_unit( stdout ) \
   if(use_environ) CALL environ_summary() \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_summary.f90

# plugin_initbase

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    mp_bands,     ONLY : intra_bgrp_comm, me_bgrp, root_bgrp_id \
USE    cell_base,    ONLY : at, alat, omega, ibrav \
USE    environ_init, ONLY : environ_initbase \
!Environ patch
' plugin_initbase.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch \
INTEGER :: ir_end \
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
! BACKWARD COMPATIBILITY \
! Compatible with QE-5.X QE-6.1.X \
!  ir_end = dfftp%nr1x*dfftp%nr2x*dfftp%npl \
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT \
#if defined (__MPI) \
    ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p) \
#else \
    ir_end = dfftp%nnr \
#endif \
! END BACKWARD COMPATIBILITY \
  IF ( use_environ ) CALL environ_initbase( dfftp%nr1, dfftp%nr2, dfftp%nr3, ibrav, alat, omega, at, & \
                                          & dfftp%nnr, ir_end, intra_bgrp_comm, me_bgrp, root_bgrp_id ) \
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_initbase.f90

# plugin_clock

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_output, ONLY : environ_clock \
!Environ patch
' plugin_clock.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL environ_clock() \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_clock.f90

# plugin_print_energies

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    control_flags,  ONLY : conv_elec \
USE    klist,          ONLY : lgauss, ltetra \
USE    environ_output, ONLY : environ_print_energies, & \
                              environ_print_fermi_shift \
!Environ patch
' plugin_print_energies.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if (use_environ) then \
     CALL environ_print_energies() \
     if (conv_elec.and.(lgauss.or.ltetra)) then \
       CALL environ_print_fermi_shift() \
     end if \
   end if \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_print_energies.f90

# plugin_init_ions

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE ions_base,            ONLY : zv, nat, nsp, ityp, tau\
USE environ_init,         ONLY : environ_initions\
!Environ patch
' plugin_init_ions.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
  IF ( use_environ ) call environ_initions( dfftp%nnr, nat, nsp, ityp, zv, tau )\
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_ions.f90

# plugin_init_cell

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE cell_base,            ONLY : at, omega\
USE environ_init,         ONLY : environ_initcell\
!Environ patch
' plugin_init_cell.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
  IF ( use_environ ) call environ_initcell( omega, at )\
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_cell.f90

# plugin_scf_energy

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_base,          ONLY : deenviron, eelectrostatic, & \
                                  esurface, evolume, eelectrolyte \
USE environ_main,          ONLY : calc_eenviron \
!Environ patch
' plugin_scf_energy.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) THEN \
        ! \
        ! compute environ contributions to total energy \
        ! \
        CALL calc_eenviron( deenviron, eelectrostatic, esurface, evolume, eelectrolyte ) \
        ! \
        plugin_etot = plugin_etot + deenviron + eelectrostatic + esurface + evolume + eelectrolyte \
        ! \
  END IF \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_scf_energy.f90

# plugin_init_potential

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_init,         ONLY : environ_initpotential \
!Environ patch
' plugin_init_potential.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) CALL environ_initpotential( dfftp%nnr, vltot ) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_potential.f90

# plugin_scf_potential

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE klist,                 ONLY : nelec \
USE environ_base,          ONLY : update_venviron, environ_thr, & \
                                  environ_restart \
USE environ_init,          ONLY : environ_initelectrons \
USE environ_main,          ONLY : calc_venviron \
!Environ patch
' plugin_scf_potential.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
     IF(use_environ) THEN \
        ! \
        ! update electrons-related quantities in environ \
        ! \
        CALL environ_initelectrons( nspin, dfftp%nnr, rhoin%of_r, nelec ) \
        ! \
        ! environ contribution to the local potential \
        ! \
        IF ( dr2 .GT. 0.0_dp ) THEN \
           update_venviron = .NOT. conv_elec .AND. dr2 .LT. environ_thr \
        ! \
        ELSE \
           update_venviron = environ_restart \
           ! for subsequent steps of optimization or dynamics, compute \
           ! environ contribution during initialization \
           IF ( .NOT. environ_restart ) environ_restart = .TRUE. \
        ENDIF \
        ! \
        IF ( update_venviron ) WRITE( stdout, 9200 ) \
        CALL calc_venviron( update_venviron, dfftp%nnr, vltot ) \
        ! \
9200 FORMAT(/"     add environment contribution to local potential") \
     ENDIF \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_scf_potential.f90

# makov_payne

# sed '/Environ MODULES BEGIN/ a\
# !Environ patch \
# USE environ_mp,  ONLY : environ_makov_payne \
# !Environ patch
# ' makov_payne.f90 > tmp.1
# 
# sed '/Environ CALLS BEGIN/ a\
# !Environ patch \
#      IF(use_environ) THEN \
#        CALL environ_makov_payne( dfftp%nnr, nspin, rho%of_r, x0, etot ) \
#      ENDIF \
# !Environ patch
# ' tmp.1 > tmp.2
# 
# mv tmp.2 makov_payne.f90

rm tmp.1

# force_lc

cat >> force_lc.f90 <<EOF
!Environ patch
subroutine external_force_lc( rhor, force )

  use kinds,            only : DP
  use cell_base,        only : at, bg, alat, omega
  use ions_base,        only : nat, ntyp => nsp, ityp, tau, zv, amass
  use fft_base,         only : dfftp
  use fft_interfaces,   only : fwfft
  use gvect,            only : ngm, gstart, ngl, nl, igtongl, g, gg, gcutm
  use lsda_mod,         only : nspin
  use vlocal,           only : strf, vloc
  use control_flags,    only : gamma_only
  use martyna_tuckerman, only: do_comp_mt, wg_corr_force
  implicit none

  real( dp ), intent(in) ::  rhor (dfftp%nnr, nspin)
  real( dp ), intent(out) :: force (3, nat)

  real( dp ), allocatable :: force_tmp(:,:)
  complex( dp ), allocatable :: auxg(:), auxr(:)

  force = 0.0_dp

  allocate(force_tmp(3,nat))

  if ( do_comp_mt) then
     force_tmp = 0.0_dp
     allocate(auxr(dfftp%nnr))
     allocate(auxg(ngm))
     auxg = cmplx(0.0_dp,0.0_dp)
     auxr = cmplx(rhor(:,1),0.0_dp)
     if ( nspin .eq. 2 ) auxr = auxr + cmplx(rhor(:,2),0.0_dp)
     call fwfft ("Dense", auxr, dfftp)
     auxg(:)=auxr(nl(:))
     call wg_corr_force(.false.,omega, nat, ntyp, ityp, ngm, g, tau, zv, strf, &
                        1, auxg, force_tmp)
     deallocate(auxr,auxg)
     force = force + force_tmp
  endif

  force_tmp = 0.0_dp
  call force_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, &
       g, rhor, nl, nspin, gstart, gamma_only, vloc, force_tmp )
  force = force + force_tmp

  return
end subroutine external_force_lc
!Environ patch
EOF

echo "- DONE!"

cd $QE_DIR
