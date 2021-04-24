! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.1
!
!    Environ 1.1 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.1 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------------------
!>
!! Module to include a real-space correction of the potential in order to
!! remove periodic boundary conditions on a partially periodic system.
!! The correction is computed for the time being using the planar average
!! approximation and is of quadratic nature: the first order proportional
!! to the dipole of the system along the direction perpendicular to the slab,
!! the second order proportional to the total charge of the system.
!!
!! The variables needed to correct periodic boundary conditions for a partially periodic
!! system. Real space correction with planar average approximation.
!!
!----------------------------------------------------------------------------------------
MODULE correction_periodic
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: e2, pi, tpi, fpi
    USE environ_types
    USE core_types
    USE environ_output
    !
    IMPLICIT NONE
    !
    REAL(DP), PARAMETER :: madelung(3) = (/2.837297479D0, 2.8883D0, 2.885D0/)
    !
    PRIVATE
    !
    PUBLIC :: calc_vperiodic, calc_fperiodic, calc_gradvperiodic
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_vperiodic(oned_analytic, charges, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(oned_analytic_core), TARGET, INTENT(IN) :: oned_analytic
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        REAL(DP), POINTER :: vperiodic(:)
        TYPE(environ_cell), POINTER :: cell
        !
        INTEGER, POINTER :: env_periodicity, slab_axis
        REAL(DP), POINTER :: alat, omega, axis_length
        REAL(DP), POINTER :: origin(:), axis(:, :)
        !
        INTEGER :: icor
        !
        REAL(DP) :: fact
        REAL(DP) :: const, charge, dipole(3), quadrupole(3)
        !
        TYPE(environ_density), TARGET :: local
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_vperiodic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL start_clock('calc_vpbc')
        !
        IF (.NOT. ASSOCIATED(potential%cell, charges%cell)) &
            CALL errore(sub_name, 'Missmatch in domains of potential and charges', 1)
        !
        IF (.NOT. ASSOCIATED(potential%cell, oned_analytic%cell)) &
            CALL errore(sub_name, 'Missmatch in domains of potential and solver', 1)
        !
        cell => potential%cell
        alat => cell%alat
        omega => cell%omega
        !
        env_periodicity => oned_analytic%d
        slab_axis => oned_analytic%axis
        axis_length => oned_analytic%size
        origin => oned_analytic%origin
        axis => oned_analytic%x
        !
        IF (env_periodicity == 0 .AND. .NOT. cell%cubic) &
            CALL errore(sub_name, &
                        'Parabolic correction in 0D is only for cubic cells', 1)
        !
        CALL init_environ_density(cell, local)
        !
        vperiodic => local%of_r
        !
        CALL multipoles_environ_density(charges, origin, charge, dipole, quadrupole)
        ! compute multipoles of the system with respect to the chosen origin
        !
        !--------------------------------------------------------------------------------
        ! Compute quadratic PBC correction
        !
        fact = e2 * tpi / omega
        !
        SELECT CASE (env_periodicity)
        CASE (0)
            !
            const = madelung(1) * charge / alat * e2 &
                    - fact * SUM(quadrupole(:)) / 3.D0
            !
            vperiodic = 0.D0
            !
            DO icor = 1, 3
                !
                vperiodic(:) = vperiodic(:) - charge * axis(icor, :)**2 + &
                               2.D0 * dipole(icor) * axis(icor, :)
                !
            END DO
            !
            vperiodic = fact / 3.D0 * vperiodic + const
        CASE (1)
            CALL errore(sub_name, 'Option not yet implemented', 1)
        CASE (2)
            !
            const = -pi / 3.D0 * charge / axis_length * e2 - &
                    fact * quadrupole(slab_axis)
            !
            vperiodic(:) = -charge * axis(1, :)**2 + &
                           2.D0 * dipole(slab_axis) * axis(1, :)
            !
            vperiodic = fact * vperiodic + const
        CASE (3)
            const = 0.D0
            vperiodic = 0.D0
        CASE DEFAULT
            !
            CALL errore(sub_name, &
                        'Unexpected option in dimensionality of PBC correction', 1)
            !
        END SELECT
        !
        potential%of_r = potential%of_r + vperiodic
        !
        CALL destroy_environ_density(local)
        !
        CALL stop_clock('calc_vpbc')
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_vperiodic
    !------------------------------------------------------------------------------------
    !>
    !! Note that in this subroutine the ionic density is implicit (and thus
    !! spread gaussians). Since the gradient of the corrective potential does
    !! not depend on the quadrupole moment of rhotot, it should be independent
    !! on the shape of the ionic density
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_gradvperiodic(oned_analytic, charges, gvtot)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(oned_analytic_core), TARGET, INTENT(IN) :: oned_analytic
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gvtot
        !
        REAL(DP), POINTER :: gvperiodic(:, :)
        TYPE(environ_cell), POINTER :: cell
        !
        INTEGER, POINTER :: env_periodicity, slab_axis
        REAL(DP), POINTER :: omega
        REAL(DP), POINTER :: origin(:), axis(:, :)
        !
        INTEGER :: icor
        !
        REAL(DP) :: fact
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        !
        TYPE(environ_gradient), TARGET :: glocal
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_gradvperiodic'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(gvtot%cell, charges%cell)) &
            CALL errore(sub_name, 'Missmatch in domains of gradient and charges', 1)
        !
        IF (.NOT. ASSOCIATED(gvtot%cell, oned_analytic%cell)) &
            CALL errore(sub_name, 'Missmatch in domains of gradient and solver', 1)
        !
        cell => gvtot%cell
        omega => cell%omega
        !
        env_periodicity => oned_analytic%d
        slab_axis => oned_analytic%axis
        origin => oned_analytic%origin
        axis => oned_analytic%x
        !
        CALL init_environ_gradient(cell, glocal)
        !
        gvperiodic => glocal%of_r
        !
        CALL multipoles_environ_density(charges, origin, charge, dipole, quadrupole)
        ! compute multipoles of the system with respect to the chosen origin
        !
        !--------------------------------------------------------------------------------
        ! Compute gradient of periodic images correction
        !
        fact = e2 * fpi / omega
        !
        SELECT CASE (env_periodicity)
        CASE (0)
            !
            DO icor = 1, 3
                gvperiodic(icor, :) = (dipole(icor) - charge * axis(icor, :)) / 3.D0
            END DO
            !
        CASE (1)
            CALL errore(sub_name, 'Option not yet implemented', 1)
        CASE (2)
            gvperiodic(slab_axis, :) = dipole(slab_axis) - charge * axis(1, :)
        CASE (3)
            gvperiodic = 0.D0
        CASE DEFAULT
            CALL errore(sub_name, 'Unexpected option', 1)
        END SELECT
        !
        gvperiodic = gvperiodic * fact
        !
        gvtot%of_r = gvtot%of_r + gvperiodic
        ! sum the periodic contribution to the total gradient of the potential
        !
        CALL destroy_environ_gradient(glocal)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_gradvperiodic
    !------------------------------------------------------------------------------------
    !>
    !! Computes the contribution to the atomic forces
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_fperiodic(oned_analytic, natoms, charges, auxiliary, f)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(oned_analytic_core), TARGET, INTENT(IN) :: oned_analytic
        INTEGER, INTENT(IN) :: natoms
        TYPE(environ_charges), TARGET, INTENT(IN) :: charges
        TYPE(environ_density), INTENT(IN) :: auxiliary
        !
        REAL(DP), INTENT(OUT) :: f(3, natoms)
        !
        INTEGER, POINTER :: ityp(:)
        REAL(DP), POINTER :: tau(:, :)
        !
        INTEGER, POINTER :: env_periodicity, slab_axis
        REAL(DP), POINTER :: alat, omega
        REAL(DP), POINTER :: origin(:)
        !
        INTEGER :: i
        !
        REAL(DP) :: fact, pos(3)
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        REAL(DP) :: ftmp(3, natoms)
        !
        TYPE(environ_density) :: local
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_fperiodic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL start_clock('calc_fpbc')
        !
        IF (.NOT. ASSOCIATED(charges%density%cell, oned_analytic%cell)) &
            CALL errore(sub_name, 'Missmatch in domains of charges and solver', 1)
        !
        IF (natoms /= charges%ions%number) &
            CALL errore(sub_name, &
                        'Missmatch in numbers of atoms passed in input and stored', 1)
        !
        alat => charges%density%cell%alat
        omega => charges%density%cell%omega
        tau => charges%ions%tau
        ityp => charges%ions%ityp
        !
        env_periodicity => oned_analytic%d
        slab_axis => oned_analytic%axis
        origin => oned_analytic%origin
        !
        CALL init_environ_density(charges%density%cell, local)
        !
        local%of_r = charges%density%of_r + auxiliary%of_r
        !
        CALL multipoles_environ_density(local, origin, charge, dipole, quadrupole)
        ! compute multipoles of the system with respect to the chosen origin
        !
        !--------------------------------------------------------------------------------
        ! Interatomic forces, quadrupole is not needed, thus the same
        ! expression holds for point-like and gaussian nuclei
        !
        fact = e2 * fpi / omega
        ftmp = 0.D0
        !
        DO i = 1, natoms
            pos(:) = (tau(:, i) - origin(:)) * alat
            !
            SELECT CASE (env_periodicity)
            CASE (0)
                ftmp(:, i) = (charge * pos(:) - dipole(:)) / 3.D0
            CASE (1)
                CALL errore(sub_name, 'Option not yet implemented', 1)
            CASE (2)
                ftmp(slab_axis, i) = charge * pos(slab_axis) - dipole(slab_axis)
            CASE (3)
                ftmp = 0.D0
            CASE DEFAULT
                CALL errore(sub_name, 'Unexpected', 1)
            END SELECT
            !
            ftmp(:, i) = ftmp(:, i) * fact * charges%ions%iontype(ityp(i))%zv
            !
        END DO
        !
        f = f + ftmp
        !
        CALL stop_clock('calc_fpbc')
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_fperiodic
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE correction_periodic
!----------------------------------------------------------------------------------------
