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
! Compatible with QE-6.3.X and QE-GIT\
  USE lsda_mod,   ONLY : nspin\
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
! Compatible with QE-5.X QE-6.1.X, QE-6.2.X\
!      prog = "PW"\
! Compatible with QE-6.3.X and QE-GIT\
! END BACKWARD COMPATIBILITY\
      CALL set_environ_output(prog, ionode, ionode_id, intra_image_comm, stdout)\
      CALL read_environ(prog,1, nspin, nat, ntyp, atm, do_comp_mt)\
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
! BACKWARD COMPATIBILITY \
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
! Compatible with QE-6.2.X, QE-6.3.X and QE-GIT \
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
USE environ_base,          ONLY : update_venviron, environ_thr, &\
                                  environ_restart, ltddfpt\
USE environ_init,          ONLY : environ_initelectrons\
USE environ_main,          ONLY : calc_venviron\
!Environ patch
' plugin_scf_potential.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
     IF(use_environ) THEN\
        !\
        ! update electrons-related quantities in environ\
        !\
        CALL environ_initelectrons( nspin, dfftp%nnr, rhoin%of_r, nelec )\
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
        CALL calc_venviron( update_venviron, dfftp%nnr, vltot )\
        !\
9200 FORMAT(/"     add environment contribution to local potential")\
     ENDIF\
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_scf_potential.f90

# plugin initialization

sed '! /***Environ MODULES BEGIN***/ a\
!Environ patch\
USE klist,            ONLY : tot_charge\
USE force_mod,        ONLY : lforce\
USE control_flags,    ONLY : lbfgs\
USE environ_base,  ONLY : lsemiconductor\
USE control_flags,    ONLY : nstep\
!Environ patch
' plugin_initialization.f90 > tmp.1

sed '! ***Environ CALLS BEGIN***\
!Environ patch\
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
IF (lsemiconductor) THEN \
CALL start_clock( "semiconductor" ) \
lforce = .TRUE. \
lbfgs = .FALSE. \
nstep = 100 \
tot_charge = 0.0 \
!WRITE( stdout, 1000) \
WRITE( stdout, 1002) tot_charge \
CALL stop_clock( "semiconductor" ) \
 \
END IF \
 \
END IF \
 \
1000 FORMAT(5X,//"*******************************************"//,& \
&"  Please cite                              "//,& \
&"  Q. Campbell, D. Fisher and I. Dabo, Phys. Rev. Mat. 3, 015404 (2019)."//,& \
&"  doi: 10.1103/PhysRevMaterials.3.015404   "//,& \
&"  In any publications resulting from this work.") \
 \
1002 FORMAT(5x,//"*******************************************"//, \
&"     Running initial calculation for flatband."//& \
&   "     Using charge of: ",F14.8,//& \
&"*******************************************") \
 \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_initialization.f90

#plugin ext_forces (where I'm hiding all the semiconductor shit)


sed '! ***Environ MODULES BEGIN*** \
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
USE environ_base,  ONLY : lsemiconductor, semiconductor, cell \
USE environ_output,  ONLY : environ_unit \
 \
USE mp,             ONLY: mp_bcast, mp_barrier, mp_sum \
USE mp_world,       ONLY: world_comm \
USE mp_images,      ONLY: intra_image_comm \
USE mp_bands,       ONLY: intra_bgrp_comm \
USE klist,            ONLY : tot_charge, nelec \
USE cell_base,        ONLY : omega \
USE lsda_mod,         ONLY : nspin \
USE scf,              ONLY : rho \
USE control_flags,    ONLY : conv_ions, nstep, istep \
USE ener,             ONLY : ef \
USE constants,        ONLY : rytoev \
USE fft_base,         ONLY : dfftp \
USE ions_base,        ONLY : nat, ityp, zv \
USE extrapolation,    ONLY : update_pot \
USE qexsd_module,     ONLY:   qexsd_set_status \
!Environ patch
' plugin_ext_forces.f90 > tmp.1

sed '! ***Environ VARIABLES BEGIN*** \
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

sed '  ! ! ***Environ CALLS BEGIN*** \
! !Environ patch \
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
ionic_charge = 0._DP \
DO na = 1, nat \
ionic_charge = ionic_charge + zv( ityp(na) ) \
END DO \
 \
 \
 \
IF (use_environ .AND. lsemiconductor) THEN \
CALL start_clock( "semiconductor" ) \
 \
chg_step = istep-1 \
!! Initializing the constraints of possible DFT charges \
! Should probably be initialized at chg_step =1 but that seemed to be \
! creating some trouble possibly \
IF (chg_step == 1) THEN \
! this is an option that feels like it should be useful to edit in the future \
IF (semiconductor%electrode_charge > 0.0) THEN \
dft_chg_max = 2.0*semiconductor%electrode_charge \
dft_chg_min = 0.0 \
ELSE \
dft_chg_min = 2.0*semiconductor%electrode_charge \
dft_chg_max = 0.0 \
END IF \
 \
END IF \
 \
 \
IF (chg_step == 0) THEN \
tot_charge = 0.7*semiconductor%electrode_charge \
semiconductor%flatband_fermi = ef!*rytoev \
conv_ions = .FALSE. \
! CALL qexsd_set_status(255) \
! CALL punch( "config" ) \
! CALL add_qexsd_step(istep) \
istep =  istep + 1 \
!CALL save_flatband_pot(dfftp%nnr) \
WRITE( stdout, 1001) semiconductor%flatband_fermi*rytoev,tot_charge \
! \
! ... re-initialize atomic position-dependent quantities \
! \
nelec = ionic_charge - tot_charge \
CALL update_pot() \
CALL hinit1() \
ELSE \
cur_fermi = ef!*rytoev \
! for now, will try to keep everything in Ry, should basically work the same \
 \
!CALL save_current_pot(dfftp%nnr,cur_fermi,cur_dchg,ss_chg,v_cut,chg_step) \
cur_dchg = semiconductor%bulk_sc_fermi - cur_fermi \
bulk_potential = (semiconductor%bulk_sc_fermi - semiconductor%flatband_fermi)*rytoev \
ss_chg = tot_chg \
!IF (ionode) THEN \
! making sure constraints are updated \
IF (semiconductor%electrode_charge > 0) THEN \
IF (v_cut >= 0.0 ) THEN \
dft_chg_max = tot_charge \
converge = .FALSE. \
ELSE IF (ss_chg < 0.0) THEN \
dft_chg_min = tot_charge \
converge = .FALSE. \
ELSE \
prev_chg2 = tot_charge \
END IF \
ELSE \
IF (v_cut <= 0.0 ) THEN \
dft_chg_min = tot_charge \
converge = .FALSE. \
ELSE IF (ss_chg > 0.0) THEN \
dft_chg_max = tot_charge \
converge = .FALSE. \
ELSE \
prev_chg2 = tot_charge \
END IF \
END IF \
CALL mp_bcast(dft_chg_min, ionode_id,intra_image_comm) \
CALL mp_bcast(dft_chg_max, ionode_id,intra_image_comm) \
IF (chg_step > 1 )THEN \
gamma_mult = (cur_chg - prev_chg)/(cur_dchg - prev_dchg) \
END IF \
WRITE(environ_unit,*)"cur_chg: ",cur_chg \
WRITE(environ_unit,*)"prev_chg: ",prev_chg \
WRITE(environ_unit,*)"cur_dchg: ",cur_dchg \
WRITE(environ_unit,*)"prev_dchg: ",prev_dchg \
WRITE(environ_unit,*)"Using gamma of ",gamma_mult \
change_vec = -gamma_mult*cur_dchg \
prev_chg = tot_charge \
! This is my way of trying to impose limited constraints with an \
! unknown constraining function. Theres almost certainly a more \
! efficient way to do this but I havent thought of it yet \
 \
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
 \
END IF \
WRITE(environ_unit,*)"DFT_min ",dft_chg_min \
WRITE(environ_unit,*)"DFT_max ",dft_chg_max \
CALL mp_bcast(tot_charge, ionode_id,intra_image_comm) \
!print *,"DFT_max",dft_chg_max \
cur_chg = tot_charge \
prev_step_size = ABS(cur_chg - prev_chg) \
prev_dchg = cur_dchg \
WRITE(environ_unit,*)"Convergeable? ",converge \
CALL mp_bcast(converge,ionode_id, intra_image_comm) \
CALL mp_bcast(prev_step_size,ionode_id,intra_image_comm) \
IF (((prev_step_size > semiconductor%charge_threshold) .OR. (.NOT. converge)) & \
& .AND. (chg_step < nstep-1))  THEN \
conv_ions = .FALSE. \
WRITE( STDOUT, 1002)& \
&chg_step,cur_fermi*rytoev,ss_chg,prev_step_size,cur_dchg,tot_charge \
!CALL qexsd_set_status(255) \
!CALL punch( "config" ) \
!CALL add_qexsd_step(istep) \
istep =  istep + 1 \
nelec = ionic_charge - tot_charge \
CALL mp_bcast(nelec, ionode_id,intra_image_comm) \
CALL update_pot() \
CALL hinit1() \
ELSE \
IF (chg_step == nstep -1) THEN \
WRITE(STDOUT,*)NEW_LINE("a")//"   Exceeded Max number steps!"//& \
&NEW_LINE("a")//"   Results probably out of accurate range"//& \
&NEW_LINE("a")//"   Smaller chg_thr recommended."//& \
&NEW_LINE("a")//"   Writing current step to q-v.dat." \
END IF \
WRITE(STDOUT, 1003)chg_step,prev_step_size,ss_chg,cur_dchg,& \
&bulk_potential \
OPEN(21,file = "q-v.dat", status = "unknown") \
WRITE(37, *)"Potential (V-V_fb)  Surface State Potential (V-V_cut)",& \
&"  Electrode Charge (e)",& \
&"  Surface States Charge (e)    ",& \
&"Electrode Charge per surface area (e/cm^2)     ",& \
&"Surface State Charge per surface area (e/cm^2)" \
surf_area = semiconductor%surf_area_per_sq_cm \
chg_per_area = semiconductor%electrode_charge/surf_area \
ss_chg_per_area = ss_chg/surf_area \
ss_potential = -bulk_potential \
CALL mp_bcast(ss_potential, ionode_id, intra_image_comm) \
!print *, bulk_potential,ss_potential \
WRITE(37, 1004)total_potential, ss_potential,& \
&semiconductor%electrode_charge, ss_chg,& \
&chg_per_area,ss_chg_per_area \
CLOSE(21) \
END IF \
END IF \
 \
CALL stop_clock( "semiconductor" ) \
END IF \
 \
 \
 \
1001 FORMAT(5x,//"***************************************************",//& \
&"     Flatband potential calculated as ",F14.8,// & \
&"     Now using initial charge of:  ",F14.8,// & \
"***************************************************") \
! \
1002 FORMAT(5x,//"***************************************************",//& \
&"     Finished Charge convergence step : ",I3,//& \
&"     DFT Fermi level calculated as ",F14.8,// & \
&"     Charge trapped in surface states: ",F14.8," e",//& \
&"     Charge Accuracy < ",F14.8,// & \
&"     Difference between bulk and DFT fermi: ",F14.8,//& \
&"     Now using DFT charge of:  ",F14.8,// & \
"***************************************************") \
1003 FORMAT(5x,//"***************************************************",//& \
&"     Finished charge convergence step : ",I3,//& \
&"     Convergence of charge with accuracy < ",F14.8," e",// & \
&"     Charge trapped in surface states: ",F14.8,//& \
&"     Difference between bulk and DFT fermi: ",F14.8,//& \
&"     Final Potential: ",F14.8," V", //& \
&"     Output written to q-v.dat       ",//& \
"***************************************************") \
1004 FORMAT(1x,4F14.8,2ES12.5) \
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_ext_forces.f90


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
