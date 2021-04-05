
! Copyright (C) 2017 Environ group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by Q. Campbell
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE correction_ms
    !------------------------------------------------------------------------------------
    !
    USE environ_types
    USE electrostatic_types
    USE environ_output, ONLY: environ_unit
    USE modules_constants, ONLY: e2, fpi, k_boltzmann_ry, pi, tpi
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: calc_vms, calc_gradvms
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Given the total explicit charge, the value of the field at the boundary is obtained by Gauss's law
    !!
    !! (1) ez = - tpi * e2 * charge * axis_length / omega / env_static_permittivity
    !!
    !! By integrating the Mott Schottky model, we can find that the schottky barrier is equal to
    !!
    !! (2) vms = fact*(ez)**2 + kbt
    !!
    !! where
    !!
    !! (3) fact = sc_permittivity/tpi / e2 /2.D0 /carrier_density
    !!
    !!  You can then find the potential as a parabolic function of distance
    !!
    !! (4) v_analytic = (distance)**2/fact/4 - ez*(distance)
    !!
    !!  This parabolic function is only accurate until the maximum is reached
    !!  We can quantify when this will happen as the depletion length
    !!
    !! (5) depletion_length = 2.D0 *fact*ez
    !!
    !!  After the depletion length, the potential will remain flat
    !!  For more details and derivation see the appendix in Ch11 of Schmickler
    !!  Interfacial Electrochemistry book
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_vms(oned_analytic, semiconductor, charges, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(oned_analytic_core), TARGET, INTENT(IN) :: oned_analytic
        TYPE(environ_semiconductor), TARGET, INTENT(IN) :: semiconductor
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        INTEGER, POINTER :: nnr
        TYPE(environ_cell), POINTER :: cell
        !
        INTEGER, POINTER :: env_periodicity, slab_axis
        REAL(DP), POINTER :: omega, axis_length
        REAL(DP), POINTER :: origin(:), axis(:, :)
        !
        REAL(DP), POINTER :: permittivity, xstern, carrier_density
        REAL(DP) :: kbt, invkbt
        !
        TYPE(environ_density), TARGET :: local
        REAL(DP), POINTER :: v(:)
        !
        INTEGER :: i, icount
        !
        REAL(DP) :: ez, fact, vms
        REAL(DP) :: arg, const, depletion_length
        REAL(DP) :: area, vtmp, vbound, distance
        REAL(DP) :: dipole(0:3), quadrupole(3)
        REAL(DP) :: tot_charge, tot_dipole(3), tot_quadrupole(3)
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_vms'
        !
        !--------------------------------------------------------------------------------
        !
        CALL start_clock('calc_vms')
        !
        IF (.NOT. ASSOCIATED(potential%cell, charges%cell)) &
            CALL errore(sub_name, 'Missmatch in domains of potential and charges', 1)
        !
        IF (potential%cell%nnr .NE. oned_analytic%n) &
            CALL errore(sub_name, 'Missmatch in domains of potential and solver', 1)
        !
        cell => potential%cell
        nnr => cell%nnr
        !
        omega => cell%omega
        env_periodicity => oned_analytic%d
        slab_axis => oned_analytic%axis
        axis_length => oned_analytic%size
        origin => oned_analytic%origin
        axis => oned_analytic%x
        !
        !--------------------------------------------------------------------------------
        ! Get parameters of semiconductor to compute analytic correction
        !
        permittivity => semiconductor%permittivity
        carrier_density => semiconductor%carrier_density
        xstern => semiconductor%simple%width
        !
        WRITE (environ_unit, *) "carrier density: ", carrier_density
        !
        WRITE (environ_unit, *) "xstern: ", xstern
        WRITE (environ_unit, *) "carrier density: ", carrier_density
        FLUSH (environ_unit)
        !
        !--------------------------------------------------------------------------------
        ! Set Boltzmann factors
        !
        kbt = semiconductor%temperature * k_boltzmann_ry
        invkbt = 1.D0 / kbt
        !
        IF (env_periodicity .NE. 2) &
            CALL errore(sub_name, &
                        'Option not yet implemented: 1D Poisson-Boltzmann &
                        &solver only for 2D systems', 1)
        !
        CALL init_environ_density(cell, local)
        !
        v => local%of_r
        !
        !--------------------------------------------------------------------------------
        ! Compute multipoles of the system wrt the chosen origin
        !
        CALL compute_dipole(nnr, charges%of_r, origin, dipole, quadrupole)
        !
        tot_charge = dipole(0)
        !
        tot_dipole = dipole(1:3)
        tot_quadrupole = quadrupole
        area = omega / axis_length
        !
        !--------------------------------------------------------------------------------
        ! First apply parabolic correction
        !
        fact = e2 * tpi / omega
        !
        const = -pi / 3.D0 * tot_charge / axis_length * e2 - &
                fact * tot_quadrupole(slab_axis)
        !
        v(:) = -tot_charge * axis(1, :)**2 + 2.D0 * tot_dipole(slab_axis) * axis(1, :)
        v(:) = fact * v(:) + const
        !
        !--------------------------------------------------------------------------------
        ! Compute the physical properties of the interface
        !
        WRITE (environ_unit, *) "charge: ", tot_charge
        ez = -tpi * e2 * tot_charge / area
        WRITE (environ_unit, *) "electric field: ", ez
        fact = 1.D0 / tpi / e2 / 2.D0 / carrier_density
        WRITE (environ_unit, *) "Prefactor: ", fact
        arg = fact * (ez**2.D0)
        vms = arg ! +kbt ! #TODO figure this out
        !
        !--------------------------------------------------------------------------------
        ! Finds the total length of the depletion region
        !
        depletion_length = 2.D0 * fact * ez
        WRITE (environ_unit, *) "depletion length: ", depletion_length
        WRITE (environ_unit, *) "vms: ", vms
        !
        vbound = 0.D0
        icount = 0
        !
        DO i = 1, nnr
            !
            IF (ABS(axis(1, i)) .GE. xstern) THEN
                icount = icount + 1
                !
                vbound = vbound + potential%of_r(i) + v(i) - &
                         ez * (ABS(axis(1, i)) - xstern)
                !
            END IF
            !
        END DO
        !
        CALL mp_sum(icount, cell%dfft%comm)
        !
        CALL mp_sum(vbound, cell%dfft%comm)
        !
        vbound = vbound / DBLE(icount)
        WRITE (environ_unit, *) "vbound: ", vbound
        !
        !--------------------------------------------------------------------------------
        ! Compute the analytic potential and charge
        !
        v = v - vbound - vms
        !
        DO i = 1, nnr
            !
            IF (ABS(axis(1, i)) .GE. xstern) THEN
                distance = ABS(axis(1, i)) - xstern
                !
                !------------------------------------------------------------------------
                ! Only applies parabolic equation if still within the depletion width
                !
                IF (distance <= depletion_length) THEN
                    !
                    vtmp = (distance)**2.D0 / fact / 4.D0 + ez * (distance)
                    ! Mott-Schottky analytic solution on the outside
                    !
                ELSE
                    vtmp = 0.D0
                END IF
                !
                v(i) = v(i) + vtmp - ez * distance
                ! remove source potential (linear) and add analytic one
                !
            END IF
            !
        END DO
        !
        potential%of_r = potential%of_r + v
        !
        CALL destroy_environ_density(local)
        !
        CALL stop_clock('calc_vms')
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_vms
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_gradvms(oned_analytic, semiconductor, charges, gradv)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(oned_analytic_core), TARGET, INTENT(IN) :: oned_analytic
        TYPE(environ_semiconductor), TARGET, INTENT(IN) :: semiconductor
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradv
        !
        INTEGER, POINTER :: nnr
        TYPE(environ_cell), POINTER :: cell
        !
        INTEGER, POINTER :: env_periodicity, slab_axis
        REAL(DP), POINTER :: omega, axis_length
        REAL(DP), POINTER :: origin(:), axis(:, :)
        !
        REAL(DP), POINTER :: permittivity, xstern, carrier_density
        REAL(DP) :: kbt, invkbt
        !
        TYPE(environ_gradient), TARGET :: glocal
        REAL(DP), POINTER :: gvms(:, :)
        !
        INTEGER :: i
        !
        REAL(DP) :: ez, fact, vms
        REAL(DP) :: arg
        REAL(DP) :: area, dvtmp_dx
        REAL(DP) :: dipole(0:3), quadrupole(3)
        REAL(DP) :: tot_charge, tot_dipole(3), tot_quadrupole(3)
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_gradvms'
        !
        !--------------------------------------------------------------------------------
        !
        CALL start_clock('calc_gvms')
        !
        IF (.NOT. ASSOCIATED(gradv%cell, charges%cell)) &
            CALL errore(sub_name, 'Missmatch in domains of potential and charges', 1)
        !
        IF (gradv%cell%nnr .NE. oned_analytic%n) &
            CALL errore(sub_name, 'Missmatch in domains of potential and solver', 1)
        !
        cell => gradv%cell
        nnr => cell%nnr
        !
        omega => cell%omega
        env_periodicity => oned_analytic%d
        slab_axis => oned_analytic%axis
        axis_length => oned_analytic%size
        origin => oned_analytic%origin
        axis => oned_analytic%x
        !
        !--------------------------------------------------------------------------------
        ! Get parameters of semiconductor to compute analytic correction
        !
        permittivity => semiconductor%permittivity
        carrier_density => semiconductor%carrier_density
        xstern => semiconductor%simple%width
        !
        !--------------------------------------------------------------------------------
        ! Set Boltzmann factors
        !
        kbt = semiconductor%temperature * k_boltzmann_ry
        invkbt = 1.D0 / kbt
        !
        IF (env_periodicity .NE. 2) &
            CALL errore(sub_name, &
                        'Option not yet implemented: 1D Poisson-Boltzmann &
                        &solver only for 2D systems', 1)
        !
        CALL init_environ_gradient(cell, glocal)
        !
        gvms => glocal%of_r
        !
        !--------------------------------------------------------------------------------
        ! Compute multipoles of the system wrt the chosen origin
        !
        CALL compute_dipole(nnr, charges%of_r, origin, dipole, quadrupole)
        !
        tot_charge = dipole(0)
        tot_dipole = dipole(1:3)
        !
        !--------------------------------------------------------------------------------
        ! First compute the gradient of parabolic correction
        !
        fact = e2 * fpi / omega
        gvms(slab_axis, :) = tot_dipole(slab_axis) - tot_charge * axis(1, :)
        gvms = gvms * fact
        !
        gradv%of_r = gradv%of_r + gvms
        !
        CALL destroy_environ_gradient(glocal)
        !
        CALL stop_clock('calc_gvms')
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_gradvms
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE correction_ms
!----------------------------------------------------------------------------------------
