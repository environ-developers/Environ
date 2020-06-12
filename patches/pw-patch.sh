#!/bin/bash
#
# Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
#
#    This file is part of Environ version 1.1
#
#    Environ 1.1 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    Environ 1.1 is distributed in the hope that it will be useful,
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
! BACKWARD COMPATIBILITY\
! Compatible with QE-5.X QE-6.1.X, QE-6.2.X\
!  USE input_parameters, ONLY : nspin\
!  USE input_parameters, ONLY : nat, ntyp, atm => atom_label\
!  USE input_parameters, ONLY : ibrav\
! Compatible with QE-6.3.X\
!  USE lsda_mod,   ONLY : nspin\
!  USE ions_base,  ONLY : nat, ntyp => nsp, atm\
!  USE cell_base,  ONLY : ibrav\
! Compatible with QE-6.4.X QE-GIT\
  USE ions_base,  ONLY : nat, ntyp => nsp, atm\
  USE cell_base,  ONLY : ibrav\
! END BACKWARD COMPATIBILITY\
  USE martyna_tuckerman, ONLY : do_comp_mt\
  USE environ_input,     ONLY : read_environ\
  USE environ_output,    ONLY : set_environ_output\
!Environ patch
' plugin_read_input.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
! BACKWARD COMPATIBILITY \
! Compatible with QE-5.X QE-6.1.X, QE-6.2.X\
!  CHARACTER( LEN = 2 ) :: prog\
! Compatible with QE-6.3.X and QE-GIT\
! END BACKWARD COMPATIBILITY\
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   IF (use_environ) THEN\
! BACKWARD COMPATIBILITY\
! Compatible with QE-5.X QE-6.1.X QE-6.2.X\
!      prog = "PW"\
! Compatible with QE-6.3.X QE-6.4.X QE-GIT\
! END BACKWARD COMPATIBILITY\
      CALL set_environ_output(prog, ionode, ionode_id, intra_image_comm, stdout)\
! BACKWARD COMPATIBILITY\
! Compatible with QE-5.X QE-6.1.X QE-6.2.X QE-6.3.X\
!      CALL read_environ(prog, 1, nspin, nat, ntyp, atm, do_comp_mt)\
! Compatible with QE-6.4.X QE-GIT\
      CALL read_environ(prog, 1, nat, ntyp, atm, do_comp_mt)\
! END BACKWARD COMPATIBILITY\
   ENDIF\
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_read_input.f90

# plugin_clean

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE environ_base, ONLY : ltddfpt\
USE environ_init, ONLY : environ_clean, environ_clean_pw, &\
                         environ_clean_tddfpt\
!Environ patch
' plugin_clean.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
! BACKWARD COMPATIBILITY \
! Compatible with QE-5.X QE-6.1.X, QE-6.2.X\
!  CHARACTER( LEN = 2 ) :: prog\
! Compatible with QE-6.3.X and QE-GIT\
! END BACKWARD COMPATIBILITY\
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   IF (use_environ) THEN\
      !\
! BACKWARD COMPATIBILITY\
! Compatible with QE-5.X QE-6.1.X, QE-6.2.X\
!       prog = "PW"\
! Compatible with QE-6.3.X and QE-GIT\
! END BACKWARD COMPATIBILITY\
      IF ( prog(1:2) == "PW" ) THEN\
         !\
         ! When called by PW, but inside a TD calculation\
         ! do not clean environ variables, they have been\
         ! already cleaned by TD. The lflag input is used\
         ! to fully clean the variable or to only clean\
         ! variables initialized during the PW run and not the\
         ! ones initialized while processing the input:\
         ! this allows NEB simulations\
         !\
         IF ( .NOT. ltddfpt ) CALL environ_clean(lflag)\
         !\
      ELSE IF ( prog(1:2) == "TD" ) THEN\
         !\
         ! When called by TD, use the flag input variable to\
         ! specify whether to clean the PW variables or\
         ! the TD variables. In both cases, the variables are\
         ! fully cleaned (no NEB with TD).\
         !\
	 IF ( .NOT. lflag ) THEN\
            CALL environ_clean_pw(.TRUE.)\
         ELSE\
            CALL environ_clean_tddfpt(.TRUE.)\
	 END IF\
         !\
      END IF\
      !\
   END IF\
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_clean.f90

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
USE    mp_bands,     ONLY : intra_bgrp_comm, me_bgrp, root_bgrp\
USE    cell_base,    ONLY : at, alat, ibrav\
USE    environ_init, ONLY : environ_initbase\
USE    gvect,        ONLY : gstart, gcutm\
!Environ patch
' plugin_initbase.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
INTEGER :: ir_end, idx0, j0, k0\
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
! BACKWARD COMPATIBILITY\
! Compatible with QE-5.X QE-6.0.X QE-6.1.X\
!  idx0 = dfftp%nr1x*dfftp%nr2x*dfftp%ipp(me_bgrp+1)\
!  ir_end = dfftp%nr1x*dfftp%nr2x*dfftp%npl\
! Compatible with QE-6.2.X QE-6.3.X QE-6.4.X QE-GIT\
#if defined (__MPI)\
    j0 = dfftp%my_i0r2p ; k0 = dfftp%my_i0r3p\
    ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p)\
#else\
    j0 = 0; k0 = 0;\
    ir_end = dfftp%nnr\
#endif\
! END BACKWARD COMPATIBILITY\
  IF ( use_environ ) CALL environ_initbase( ibrav, alat, at, &\
                             & intra_bgrp_comm, me_bgrp, root_bgrp, &\
                             & gcutm, gstart )\
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
USE    environ_output, ONLY : environ_print_energies, & \
                              environ_print_potential_warning \
!Environ patch
' plugin_print_energies.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if (use_environ) then \
     CALL environ_print_energies() \
     if (conv_elec) then \
       CALL environ_print_potential_warning() \
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
USE fft_interfaces,       ONLY : invfft\
USE gvect,                ONLY : igtongl\
USE control_flags,        ONLY : gamma_only\
USE vlocal,               ONLY : vloc\
!Environ patch
' plugin_init_ions.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
INTEGER :: i\
COMPLEX( DP ), DIMENSION(:), ALLOCATABLE :: aux\
REAL( DP ), DIMENSION(:,:), ALLOCATABLE :: vloc_of_r\
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
IF ( use_environ ) THEN\
   ALLOCATE(aux(dfftp%nnr))\
   ALLOCATE(vloc_of_r(dfftp%nnr,nsp))\
   DO i = 1, nsp\
      aux(dfftp%nl(:)) = CMPLX( vloc(igtongl(:),i), 0.0_DP, KIND=DP )\
      IF ( gamma_only ) aux(dfftp%nlm(:)) = CONJG( aux(dfftp%nl(:)) )\
      CALL invfft( "Rho", aux, dfftp )\
      vloc_of_r(:,i) = DBLE( aux(:) )\
   ENDDO\
   CALL environ_initions( dfftp%nnr, nat, nsp, ityp, zv, tau, vloc_of_r )\
   DEALLOCATE(aux)\
   DEALLOCATE(vloc_of_r)\
ENDIF\
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_init_ions.f90

# plugin_init_cell

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE cell_base,            ONLY : at\
USE environ_init,         ONLY : environ_initcell\
!Environ patch
' plugin_init_cell.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
  IF ( use_environ ) call environ_initcell( at )\
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_cell.f90

# plugin_scf_energy

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_base,          ONLY : deenviron, eelectrostatic, & \
                                  esurface, evolume, econfine, eelectrolyte \
USE environ_main,          ONLY : calc_eenviron \
!Environ patch
' plugin_scf_energy.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) THEN \
        ! \
        ! compute environ contributions to total energy \
        ! \
        CALL calc_eenviron( deenviron, eelectrostatic, esurface, evolume, econfine, eelectrolyte ) \
        ! \
        plugin_etot = plugin_etot + deenviron + eelectrostatic + esurface + evolume + econfine + eelectrolyte \
        ! \
  END IF \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_scf_energy.f90

# plugin_init_potential

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_init,         ONLY : environ_initpotential\
!Environ patch
' plugin_init_potential.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) CALL environ_initpotential( dfftp%nnr, vltot )\
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_potential.f90

# plugin_scf_potential

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE klist,                 ONLY : nelec\
USE control_flags,         ONLY : lscf\
USE environ_base,          ONLY : update_venviron, environ_thr, &\
                                  environ_restart, ltddfpt\
USE environ_init,          ONLY : environ_initelectrons\
USE environ_main,          ONLY : calc_venviron\
USE environ_output,        ONLY : environ_print_potential_shift\
!Environ patch
' plugin_scf_potential.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
INTEGER :: local_verbose\
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
     IF(use_environ) THEN\
        !\
        ! reduce output at each scf iteration\
        !\
        local_verbose = 0\
        IF ( .NOT. lscf .OR. conv_elec ) local_verbose = 1\
        !\
        ! update electrons-related quantities in environ\
        !\
! BACKWARD COMPATIBILITY\
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X\
!        CALL environ_initelectrons( nspin, dfftp%nnr, rhoin%of_r, nelec )\
! Compatible with QE-6.4.X QE-GIT\
        CALL environ_initelectrons( dfftp%nnr, rhoin%of_r(:,1), nelec )\
! END BACKWARD COMPATIBILITY\
        !\
        ! environ contribution to the local potential\
        !\
        IF ( dr2 .GT. 0.0_dp ) THEN\
           update_venviron = .NOT. conv_elec .AND. dr2 .LT. environ_thr\
        !\
        ELSE\
           update_venviron = environ_restart .OR. ltddfpt\
           ! for subsequent steps of optimization or dynamics, compute\
           ! environ contribution during initialization\
           IF ( .NOT. environ_restart ) environ_restart = .TRUE.\
        ENDIF\
        !\
        IF ( update_venviron ) WRITE( stdout, 9200 )\
        CALL calc_venviron( update_venviron, dfftp%nnr, vltot, local_verbose )\
        !\
        IF ( .NOT. lscf .OR. conv_elec ) THEN\
          CALL environ_print_potential_shift()\
        END IF\
        !\
9200 FORMAT(/"     add environment contribution to local potential")\
     ENDIF\
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_scf_potential.f90

# plugin_check

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
IF (use_environ) CALL errore( calling_subroutine, &\
   & "Calculation not compatible with Environ embedding", 1)\
!Environ patch
' plugin_check.f90 > tmp.1

mv tmp.1 plugin_check.f90

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

# force_lc

cat > tmp.1 <<EOF
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
# BACKWARD COMPATIBILITY
# Compatible with QE-5.X QE-6.1.X, QE-6.2.X
# cat tmp.1 >> force_lc.f90
# Compatible with QE-6.3.X and QE-GIT
# END BACKWARD COMPATIBILITY
rm tmp.1

echo "- DONE!"

cd $QE_DIR
