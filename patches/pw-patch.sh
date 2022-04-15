#!/bin/bash
#----------------------------------------------------------------------------------------
#
# Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
#
#----------------------------------------------------------------------------------------
#
#     This file is part of Environ version 2.0
#
#     Environ 2.0 is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 2 of the License, or
#     (at your option) any later version.
#
#     Environ 2.0 is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more detail, either the file
#     `License' in the root directory of the present distribution, or
#     online at <http://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------------------------
#
# Authors: Oliviero Andreussi (Department of Physics, UNT)
#          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
#          Ismaila Dabo       (DMSE, Penn State)
#          Edan Bainglass     (Department of Physics, UNT)
#
#----------------------------------------------------------------------------------------
#
# PATCH script for plugin files and Makefile in PW/src
#
#----------------------------------------------------------------------------------------

cd $PW_SRC

patch_makefile

check_src_patched
if test "$PATCHED" == 1; then
   return
else
   message "Patching"
fi

echo "#Please do not remove or modify this file" >Environ_PATCH
echo "#It keeps track of patched versions of the Environ addson package" >>Environ_PATCH
echo "$ENVIRON_VERSION" >>Environ_PATCH

# plugin_int_forces

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
  USE environ_api, ONLY : environ\
!Environ patch
' plugin_int_forces.f90 >tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
  REAL(DP), ALLOCATABLE :: force_environ(:,:)\
!Environ patch
' tmp.1 >tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
  IF (use_environ) THEN\
    !\
    ALLOCATE(force_environ(3, nat))\
    !\
    ! ... Add environment contributions\
    !\
    CALL environ%calc%force( nat, force_environ )\
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
  END IF\
  !\
9001 FORMAT(5x,"The global environment contribution to forces")\
9002 FORMAT(5X,"atom ",I4," type ",I2,"   force = ",3F14.8)\
!Environ patch
' tmp.2 >tmp.1

mv tmp.1 plugin_int_forces.f90

# plugin_read_input

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
  USE io_global,         ONLY : ionode, ionode_id, stdout\
  USE mp_images,         ONLY : intra_image_comm\
  USE environ_api,       ONLY : environ\
!Environ patch
' plugin_read_input.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   IF (use_environ) THEN\
      !\
      CALL environ%init_interface()\
      !\
      CALL environ%init_io(ionode, ionode_id, intra_image_comm, stdout, ionode)\
      !\
      CALL environ%read_input()\
      !\
      CALL environ%setup%init()\
      !\
      IF (prog == "TD") CALL environ%setup%set_tddfpt(.TRUE.)\
      !\
   ENDIF\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_read_input.f90

# plugin_clean

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_clean.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   IF (use_environ) THEN\
      !\
      IF (prog(1:2) == "PW") THEN\
         !\
         ! When called by PW, but inside a TD calculation\
         ! do not clean environ variables, they have been\
         ! already cleaned by TD. The lflag input is used\
         ! to fully clean the variable or to only clean\
         ! variables initialized during the PW run and not the\
         ! ones initialized while processing the input:\
         ! this allows NEB simulations\
         !\
         IF (lflag) THEN\
            CALL environ%destroy(1) ! NEB image reading phase\
         ELSE IF (.NOT. environ%setup%is_tddfpt()) THEN\
            CALL environ%destroy(2)\
         END IF\
         !\
      ELSE IF ( prog(1:2) == "TD" ) THEN\
         !\
         ! When called by TD, use the flag input variable to\
         ! specify whether to clean the PW variables or\
         ! the TD variables. In both cases, the variables are\
         ! fully cleaned (no NEB with TD).\
         !\
         IF (.NOT. lflag) THEN\
            CALL environ%destroy(3)\
         ELSE\
            CALL environ%destroy(4)\
         END IF\
         !\
      END IF\
      !\
   END IF\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_clean.f90

# plugin_summary

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE io_global,   ONLY : stdout\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_summary.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   IF (use_environ) CALL environ%setup%print_summary(stdout)\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_summary.f90

# plugin_initbase

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE kinds,             ONLY : DP\
USE mp_bands,          ONLY : intra_bgrp_comm, me_bgrp, root_bgrp\
USE cell_base,         ONLY : at, alat\
USE ions_base,         ONLY : nat, nsp, ityp, atm, zv\
USE gvect,             ONLY : gcutm\
USE martyna_tuckerman, ONLY : do_comp_mt\
USE environ_api,       ONLY : environ\
!Environ patch
' plugin_initbase.f90 >tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
REAL(DP), ALLOCATABLE :: at_scaled(:, :)\
REAL(DP) :: gcutm_scaled\
INTEGER :: nr(3)\
CHARACTER(LEN=80) :: sub_name = "plugin_initbase"\
!Environ patch
' tmp.1 >tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
  IF (use_environ) THEN\
      !\
      IF (alat < 1.D-8) CALL errore(sub_name, "Wrong alat", 1)\
      !\
      ALLOCATE (at_scaled(3, 3))\
      at_scaled = at * alat\
      !\
      gcutm_scaled = gcutm / alat**2\
      !\
      nr(1) = dfftp%nr1\
      nr(2) = dfftp%nr2\
      nr(3) = dfftp%nr3\
      !\
      CALL environ%setup%init_cell(intra_bgrp_comm, at_scaled, gcutm=gcutm_scaled, nr=nr)\
      !\
      DEALLOCATE (at_scaled)\
      !\
      CALL environ%setup%init_numerical(do_comp_mt)\
      !\
      CALL environ%main%init(nat, nsp, atm, ityp, zv)\
      !\
  END IF\
!Environ patch
' tmp.2 >tmp.1

mv tmp.1 plugin_initbase.f90

# plugin_clock

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_clock.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   if(use_environ) CALL environ%setup%print_clocks()\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_clock.f90

# plugin_print_energies

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE control_flags, ONLY : conv_elec\
USE environ_api,   ONLY : environ\
!Environ patch
' plugin_print_energies.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   if (use_environ) then\
     CALL environ%main%print_energies("PW")\
     if (conv_elec) then\
       CALL environ%setup%print_potential_warning()\
     end if\
   end if\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_print_energies.f90

# plugin_init_ions

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE cell_base,   ONLY : alat\
USE ions_base,   ONLY : nat, tau\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_init_ions.f90 >tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
REAL(DP), ALLOCATABLE :: tau_scaled(:, :)\
!Environ patch
' tmp.1 >tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
IF (use_environ) THEN\
   ALLOCATE (tau_scaled(3, nat))\
   tau_scaled = tau * alat\
   !\
   CALL environ%main%update_ions(nat, tau_scaled)\
   !\
   DEALLOCATE (tau_scaled)\
END IF\
!Environ patch
' tmp.2 >tmp.1

mv tmp.1 plugin_init_ions.f90

# plugin_init_cell

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE cell_base,   ONLY : at, alat\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_init_cell.f90 >tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
REAL(DP), ALLOCATABLE :: at_scaled(:, :)\
!Environ patch
' tmp.1 >tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
IF ( use_environ ) THEN\
   ALLOCATE (at_scaled(3, 3))\
   at_scaled = at * alat\
   !\
   CALL environ%update_cell(at_scaled)\
   !\
   DEALLOCATE (at_scaled)\
END IF\
!\
!Environ patch
' tmp.2 >tmp.1

mv tmp.1 plugin_init_cell.f90

# plugin_scf_energy

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_scf_energy.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
IF(use_environ) THEN\
   !\
   ! compute environ contributions to total energy\
   !\
   ! Note: plugin_etot is set to 0.0_dp right before\
   !       this routine is called\
   !\
   CALL environ%calc%denergy(plugin_etot)\
   !\
   CALL environ%calc%energy(plugin_etot)\
   !\
END IF\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_scf_energy.f90

# plugin_init_potential

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_init_potential.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
  IF(use_environ) CALL environ%main%update_potential( dfftp%nnr, vltot )\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_init_potential.f90

# plugin_scf_potential

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE kinds,          ONLY : DP\
USE global_version, ONLY : version_number\
USE klist,          ONLY : nelec\
USE control_flags,  ONLY : lscf\
USE lsda_mod,       ONLY : nspin\
USE environ_api,    ONLY : environ\
!Environ patch
' plugin_scf_potential.f90 >tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
LOGICAL :: update_venviron\
INTEGER :: local_verbose\
REAL(DP), ALLOCATABLE :: rhoaux(:)\
!Environ patch
' tmp.1 >tmp.2

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
        ALLOCATE ( rhoaux(dfftp%nnr) )\
        rhoaux(:) = rhoin%of_r(:, 1)\
        !\
        IF ( version_number == "6.3" ) THEN\
            IF ( nspin == 2 ) rhoaux(:) = rhoaux(:) + rhoin%of_r(:, 2)\
        END IF\
        !\
        CALL environ%update_electrons( rhoaux, nelec )\
        !\
        ! environ contribution to the local potential\
        !\
        IF ( dr2 .GT. 0.0_dp ) THEN\
           update_venviron = .NOT. conv_elec .AND. dr2 .LT. environ%setup%get_threshold()\
        !\
        ELSE\
           update_venviron = environ%setup%is_restart() .OR. environ%setup%is_tddfpt()\
           ! for subsequent steps of optimization or dynamics, compute\
           ! environ contribution during initialization\
           CALL environ%setup%set_restart(.TRUE.)\
        ENDIF\
        !\
        IF ( update_venviron ) WRITE( stdout, 9200 )\
        !\
        CALL environ%calc%potential(update_venviron, local_verbose)\
        !\
        vltot = environ%main%get_vzero(dfftp%nnr) + environ%main%get_dvtot(dfftp%nnr)\
        !\
        IF ( .NOT. lscf .OR. conv_elec ) CALL environ%main%print_potential_shift()\
        !\
9200 FORMAT(/"     add environment contribution to local potential")\
     ENDIF\
!Environ patch
' tmp.2 >tmp.1

mv tmp.1 plugin_scf_potential.f90

# plugin_check

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
IF (use_environ) CALL errore( calling_subroutine, &\
   & "Calculation not compatible with Environ embedding", 1)\
!Environ patch
' plugin_check.f90 >tmp.1

mv tmp.1 plugin_check.f90

rm tmp.2

# plugin initialization
# Note, when I tried this from a fresh compilation, it didn't actually patch in
# may need a different spot to place this and plugin_ext_forces

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE klist,            ONLY : tot_charge\
USE control_flags,    ONLY : lbfgs, lforce => tprnfor\
USE control_flags,    ONLY : nstep\
USE environ_api,      ONLY : environ\
!Environ patch
' plugin_initialization.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
!\
\
! *****************************************************************************\
!\
! This checks on whether semiconductor optimization is used and either starts \
! the initial calculation of flatband potential or reads flatband potential from \
! file according to user input \
! \
! ***************************************************************************** \
 \
IF (use_environ) THEN \
 \
IF (environ%setup%lmsgcs) THEN \
CALL start_clock( "semiconductor" ) \
lforce = .TRUE. \
lbfgs = .FALSE. \
nstep = 100 \
tot_charge = 0.0 \
CALL stop_clock( "semiconductor" ) \
 \
END IF \
 \
END IF \
 \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_initialization.f90

#plugin_ext_forces (where I'm hiding all the semiconductor shit)


sed '/Environ MODULES BEGIN/ a\
!Environ patch \
!------------------------------------------------ \
! \
!Note: I am using the forces plugin as a backdoor \
!for the semiconductor loop. Its kinda off, but it works \
!If youre actually interested in plugin forces, check \
!the plugin_int_forces module \
! \
!------------------------------------------------ \
 \
\
USE class_io,           ONLY : io\
 \
USE mp,             ONLY: mp_bcast, mp_barrier, mp_sum \
USE mp_world,       ONLY: world_comm \
USE mp_images,      ONLY: intra_image_comm \
USE mp_bands,       ONLY: intra_bgrp_comm \
USE klist,            ONLY : tot_charge, nelec \
USE cell_base,        ONLY : omega \
USE lsda_mod,         ONLY : nspin \
USE control_flags,    ONLY : conv_ions, nstep, istep \
USE ener,             ONLY : ef \
USE constants,        ONLY : rytoev \
USE fft_base,         ONLY : dfftp \
USE ions_base,        ONLY : nat, ityp, zv \
USE extrapolation,    ONLY : update_pot \
USE qexsd_module,     ONLY:   qexsd_set_status \
USE environ_api,      ONLY : environ\
!Environ patch
' plugin_ext_forces.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch \
\
SAVE \
REAL(DP)                  ::   cur_chg \
REAL(DP)                  ::   prev_chg, prev_chg2 \
REAL(DP)                  ::   cur_dchg \
REAL(DP)                  ::   cur_fermi \
REAL(DP)                  ::   prev_dchg \
REAL(DP)                  ::   gamma_mult \
REAL(DP)                  ::   prev_step_size \
REAL(DP)                  ::   ss_chg, charge \
INTEGER                   ::   chg_step, na \
REAL(DP)                  ::   surf_area \
REAL(DP)                  :: chg_per_area \
REAL(DP)                  :: ss_chg_per_area \
REAL(DP)                  :: ss_potential, total_potential \
REAL(DP)                  :: dft_chg_max, dft_chg_min \
REAL(DP)                  :: change_vec \
REAL(DP)                  :: v_cut, bulk_potential \
REAL(DP)                  :: ionic_charge \
LOGICAL                   :: converge \
! !Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
 \
!************************************************* \
! \
! This section designed to run after a call to electrons. Basically, it checks \
! whether the semiconductor charge has converged and then updates the relevant \
! quantities (tot_charge) accordingly \
! \
!************************************************* \
 \
gamma_mult = 0.15 \
\
\
converge = .TRUE. \
\
! calculating ionic charge \
ionic_charge = 0._DP \
DO na = 1, nat \
ionic_charge = ionic_charge + zv( ityp(na) ) \
END DO \
\
 \
 \
IF (use_environ .AND. environ%setup%lmsgcs) THEN \
CALL start_clock( "semiconductor" ) \
\
chg_step = istep \
! Initializing the constraints of possible DFT charges \
! Should probably be initialized at chg_step =1 but that seemed to be \
! creating some trouble possibly \
\
! Set dft charge boundaries \
IF (chg_step == 1) THEN \
! this is an option that feels like it should be useful to edit in the future \
\
IF (environ%main%semiconductor%base%electrode_charge > 0.0) THEN \
dft_chg_max = environ%main%semiconductor%base%electrode_charge \
dft_chg_min = 0.0 \
ELSE \
dft_chg_min = environ%main%semiconductor%base%electrode_charge \
dft_chg_max = 0.0 \
END IF \
\
END IF \
 \
 \
! \
IF (chg_step == 0) THEN \
! After first scf step, extract fermi, set charge \
environ%main%semiconductor%base%flatband_fermi = ef!*rytoev \
tot_charge = 0.7*environ%main%semiconductor%base%electrode_charge \
environ%main%semiconductor%base%slab_charge = tot_charge \
environ%main%environment_charges%externals%functions%array(1)%volume = -(environ%main%semiconductor%base%electrode_charge - tot_charge) \
environ%main%environment_charges%externals%functions%array(2)%volume = environ%main%semiconductor%base%electrode_charge \
conv_ions = .FALSE. \
istep =  istep + 1 \
WRITE( stdout, 1001) environ%main%semiconductor%base%flatband_fermi*rytoev,tot_charge \
!  \
! ... re-initialize atomic position-dependent quantities  \
!  \
nelec = ionic_charge - tot_charge \
CALL update_pot() \
CALL hinit1() \
ELSE \
 \
! Calculate steepest descent for changing charge at each scf step \
 \
cur_fermi = ef!*rytoev  \
! for now, will try to keep everything in Ry, should basically work the same  \
 \
 \
! not calling it dfermi because difference is only used as a guide to changing \
! charge. Calling it dchg to conform with steepest descent model  \
cur_dchg = environ%main%semiconductor%base%bulk_sc_fermi - cur_fermi \
bulk_potential = (environ%main%semiconductor%base%bulk_sc_fermi - environ%main%semiconductor%base%flatband_fermi)*rytoev \
ss_chg = environ%main%semiconductor%base%ss_chg \
 \
! making sure constraints are updated  \
! 1) if electrode chg positive and dft is negative stop it \
! 2) if electrode chg negative and dft is positive stop it   \
IF (environ%main%semiconductor%base%electrode_charge > 0) THEN \
IF (ss_chg < 0.0) THEN \
dft_chg_min = tot_charge \
converge = .FALSE. \
ELSE \
prev_chg2 = tot_charge \
END IF \
ELSE \
IF (ss_chg > 0.0) THEN \
dft_chg_max = tot_charge \
converge = .FALSE. \
ELSE \
prev_chg2 = tot_charge \
END IF \
END IF \
CALL mp_bcast(dft_chg_min, ionode_id,intra_image_comm) \
CALL mp_bcast(dft_chg_max, ionode_id,intra_image_comm) \
 \
! Updating the steepest descent parameter, gamma_mult \
IF (chg_step > 1 )THEN \
gamma_mult = (cur_chg - prev_chg)/(cur_dchg - prev_dchg) \
END IF \
 \
change_vec = -gamma_mult*cur_dchg \
prev_chg = tot_charge \
 \
! Making sure that charge doesnt go above dft_chg_max or below dft_chg_min    \
! Updating tot_charge \
IF ((tot_charge + change_vec) > dft_chg_max ) THEN \
IF (tot_charge >= dft_chg_max) THEN \
tot_charge = prev_chg2 + 0.7*(dft_chg_max-prev_chg2) \
ELSE \
tot_charge = tot_charge + 0.7*(dft_chg_max-tot_charge) \
END IF \
ELSE IF ((tot_charge + change_vec) < dft_chg_min) THEN \
IF (tot_charge <= dft_chg_min) THEN \
tot_charge = prev_chg2 - 0.7*(prev_chg2-dft_chg_min) \
ELSE \
tot_charge = tot_charge - 0.7*(tot_charge-dft_chg_min) \
END IF \
 \
ELSE \
tot_charge = tot_charge + change_vec \
END IF \
 \
 \
CALL mp_bcast(tot_charge, ionode_id,intra_image_comm) \
!updating variables based on new_tot_charge \
cur_chg = tot_charge \
prev_step_size = ABS(cur_chg - prev_chg) \
prev_dchg = cur_dchg \
CALL mp_bcast(converge,ionode_id, intra_image_comm) \
CALL mp_bcast(prev_step_size,ionode_id,intra_image_comm) \
 \
! decide if loop has converged based on change in charge \
IF (((prev_step_size > environ%main%semiconductor%base%charge_threshold) .OR. (.NOT. converge)) & \
& .AND. (chg_step < nstep-1))  THEN \
! not converged \
conv_ions = .FALSE. \
WRITE( STDOUT, 1002)& \
&chg_step,cur_fermi*rytoev,ss_chg,prev_step_size,cur_dchg,tot_charge \
istep =  istep + 1 \
nelec = ionic_charge - tot_charge \
environ%main%semiconductor%base%slab_charge = tot_charge \
environ%main%environment_charges%externals%functions%array(1)%volume=-(environ%main%semiconductor%base%electrode_charge - tot_charge) \
CALL mp_bcast(nelec, ionode_id,intra_image_comm) \
CALL update_pot() \
CALL hinit1() \
ELSE \
! converged  \
! case where about to exceed max number of steps \
IF (chg_step == nstep -1) THEN \
WRITE(STDOUT,*)NEW_LINE("a")//"   Exceeded Max number steps!"//& \
&NEW_LINE("a")//"   Results probably out of accurate range"//& \
&NEW_LINE("a")//"   Larger chg_thr recommended."//& \
&NEW_LINE("a")//"   Writing current step to q-v.dat." \
END IF \
 \
!writing output for successful converge \
WRITE(STDOUT, 1003)chg_step,prev_step_size,ss_chg,cur_dchg,& \
&bulk_potential \
OPEN(21,file = "q-v.dat", status = "unknown") \
WRITE(21, *)"Potential (V-V_fb)  Surface State Potential (V-V_cut)",& \
&"  Electrode Charge (e)",& \
&"  Surface States Charge (e)    ",& \
&"Electrode Charge per surface area (e/cm^2)     ",& \
&"Surface State Charge per surface area (e/cm^2)" \
surf_area = environ%main%semiconductor%base%surf_area_per_sq_cm \
chg_per_area = environ%main%semiconductor%base%electrode_charge/surf_area \
ss_chg_per_area = ss_chg/surf_area \
ss_potential = environ%main%semiconductor%base%ss_v_cut \
CALL mp_bcast(ss_potential, ionode_id, intra_image_comm) \
WRITE(21, 1004) -bulk_potential, ss_potential,& \
&environ%main%semiconductor%base%electrode_charge, ss_chg,& \
&chg_per_area,ss_chg_per_area \
CLOSE(21) \
END IF \
END IF \
 \
CALL stop_clock( "semiconductor" ) \
END IF  \
 \
1001 FORMAT(/, 5X, 44("*"), //, & \
               5X, "flatband potential        = ", F16.8, //, & \
               5X, "initial DFT charge        = ", F16.8, //, & \
               5X, 44("*"), /) \
! \
1002 FORMAT(/, 5X, 44("*"), //, & \
               5X, "charge convergence step   = ", I16, //, & \
               5X, "DFT fermi level           = ", F16.8, /, & \
               5X, "charge in surface states  = ", F16.8, /, & \
               5X, "charge accuracy           < ", F16.8, /, & \
               5X, "bulk/DFT fermi difference = ", F16.8, /, & \
               5X, "DFT charge                = ", F16.8, //, & \
               5X, 44("*"), /) \
! \
1003 FORMAT(/, 5X, 44("*"), //, & \
               5X, "charge convergence step   = ", I16, //, & \
               5X, "converged charge accuracy < ", F16.8, /, & \
               5X, "charge in surface states  = ", F16.8, /, & \
               5X, "bulk/DFT fermi difference = ", F16.8, /, & \
               5X, "final potential (V)       = ", F16.8, /, & \
               5X, "output written to q-v.dat", //, & \
               5X, 44("*"), /) \
! \
1004 FORMAT(1X, 4F14.8, 2ES12.5) \
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_ext_forces.f90


rm tmp.2

printf " done!\n"

cd $QE_DIR
