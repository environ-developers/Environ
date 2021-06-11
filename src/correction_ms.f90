
! Copyright (C) 2017 Environ group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------------------
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
    USE env_mp, ONLY: env_mp_sum
    !
    USE environ_param, ONLY: DP, e2, fpi, k_boltzmann_ry, pi, tpi
    !
    USE types_core, ONLY: oned_analytic_core
    USE types_cell, ONLY: environ_cell
    USE types_physical, ONLY: environ_semiconductor
    USE types_representation, ONLY: environ_density, environ_gradient
    !
    USE utils_density, ONLY: init_environ_density, destroy_environ_density
    USE utils_gradient, ONLY: init_environ_gradient, destroy_environ_gradient
    !
    USE tools_math, ONLY: multipoles_environ_density
    !
    USE environ_output, ONLY: ionode, verbose, environ_unit
    !
    !------------------------------------------------------------------------------------
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
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_vms'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        IF (.NOT. ASSOCIATED(potential%cell, charges%cell)) &
            CALL env_errore(sub_name, 'Mismatch in domains of potential and charges', 1)
        !
        IF (potential%cell%nnr /= oned_analytic%n) &
            CALL env_errore(sub_name, 'Mismatch in domains of potential and solver', 1)
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
        !--------------------------------------------------------------------------------
        ! Set Boltzmann factors
        !
        kbt = semiconductor%temperature * k_boltzmann_ry
        invkbt = 1.D0 / kbt
        !
        IF (env_periodicity /= 2) &
            CALL env_errore(sub_name, &
                            'Option not yet implemented: 1D Poisson-Boltzmann solver &
                            &only for 2D systems', 1)
        !
        CALL init_environ_density(cell, local)
        !
        v => local%of_r
        area = omega / axis_length
        !
        CALL multipoles_environ_density(charges, origin, charge, dipole, quadrupole)
        ! compute multipoles of the system w.r.t the chosen origin
        !
        !--------------------------------------------------------------------------------
        ! First apply parabolic correction
        !
        fact = e2 * tpi / omega
        const = -pi / 3.D0 * charge / axis_length * e2 - fact * quadrupole(slab_axis)
        !
        v(:) = -charge * axis(1, :)**2 + 2.D0 * dipole(slab_axis) * axis(1, :)
        v(:) = fact * v(:) + const
        !
        !--------------------------------------------------------------------------------
        ! Compute the physical properties of the interface
        !
        ez = -tpi * e2 * charge / area
        fact = 1.D0 / tpi / e2 / 2.D0 / carrier_density
        arg = fact * (ez**2.D0)
        vms = arg ! +kbt ! #TODO figure this out
        !
        !--------------------------------------------------------------------------------
        ! Finds the total length of the depletion region
        !
        depletion_length = 2.D0 * fact * ez
        !
        IF (ionode .AND. verbose > 0) THEN
            WRITE (environ_unit, *) "depletion length: ", depletion_length
            WRITE (environ_unit, *) "vms: ", vms
        END IF
        !
        vbound = 0.D0
        icount = 0
        !
        DO i = 1, nnr
            !
            IF (ABS(axis(1, i)) >= xstern) THEN
                icount = icount + 1
                !
                vbound = vbound + potential%of_r(i) + v(i) - &
                         ez * (ABS(axis(1, i)) - xstern)
                !
            END IF
            !
        END DO
        !
        CALL env_mp_sum(icount, cell%dfft%comm)
        !
        CALL env_mp_sum(vbound, cell%dfft%comm)
        !
        vbound = vbound / DBLE(icount)
        !
        !--------------------------------------------------------------------------------
        ! Compute the analytic potential and charge
        !
        v = v - vbound - vms
        !
        DO i = 1, nnr
            !
            IF (ABS(axis(1, i)) >= xstern) THEN
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
        CALL env_stop_clock(sub_name)
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
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_gradvms'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        IF (.NOT. ASSOCIATED(gradv%cell, charges%cell)) &
            CALL env_errore(sub_name, 'Mismatch in domains of potential and charges', 1)
        !
        IF (gradv%cell%nnr /= oned_analytic%n) &
            CALL env_errore(sub_name, 'Mismatch in domains of potential and solver', 1)
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
        IF (env_periodicity /= 2) &
            CALL env_errore(sub_name, &
                            'Option not yet implemented: 1D Poisson-Boltzmann solver &
                            &only for 2D systems', 1)
        !
        CALL init_environ_gradient(cell, glocal)
        !
        gvms => glocal%of_r
        !
        CALL multipoles_environ_density(charges, origin, charge, dipole, quadrupole)
        ! compute multipoles of the system w.r.t the chosen origin
        !
        !--------------------------------------------------------------------------------
        ! First compute the gradient of parabolic correction
        !
        fact = e2 * fpi / omega
        gvms(slab_axis, :) = dipole(slab_axis) - charge * axis(1, :)
        gvms = gvms * fact
        !
        gradv%of_r = gradv%of_r + gvms
        !
        CALL destroy_environ_gradient(glocal)
        !
        CALL env_stop_clock(sub_name)
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
