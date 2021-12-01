!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 2.0
!
!     Environ 2.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 2.0 is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more detail, either the file
!     `License' in the root directory of the present distribution, or
!     online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Quinn Campbell     (Sandia National Laboratories)
!          Ismaila Dabo       (DMSE, Penn State)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_core_1da_electrostatics
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_mp, ONLY: env_mp_sum
    !
    USE environ_param, ONLY: DP, e2, K_BOLTZMANN_RY, pi, tpi, fpi, madelung
    !
    USE class_density
    USE class_gradient
    USE class_function
    USE class_function_gaussian
    !
    USE class_core_1da
    !
    USE class_electrolyte_base
    USE class_semiconductor_base
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE, EXTENDS(core_1da), PUBLIC :: core_1da_electrostatics
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: calc_1da_vperiodic, calc_1da_vgcs, calc_1da_vms
        PROCEDURE :: calc_1da_gradvperiodic, calc_1da_gradvgcs, calc_1da_gradvms
        PROCEDURE :: calc_1da_fperiodic
        !
        !--------------------------------------------------------------------------------
    END TYPE core_1da_electrostatics
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                 CORRECTION METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_1da_vperiodic(this, charges, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        CLASS(core_1da_electrostatics), TARGET, INTENT(INOUT) :: this
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        REAL(DP), POINTER :: vperiodic(:)
        !
        INTEGER :: icor
        !
        REAL(DP) :: fact
        REAL(DP) :: const, charge, dipole(3), quadrupole(3)
        !
        TYPE(environ_density), TARGET :: local
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_1da_vperiodic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(potential%cell, charges%cell)) &
            CALL io%error(sub_name, 'Mismatch in domains of potential and charges', 1)
        !
        IF (.NOT. ASSOCIATED(potential%cell, this%cell)) &
            CALL io%error(sub_name, 'Mismatch in domains of potential and solver', 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => potential%cell, &
                   omega => potential%cell%omega, &
                   env_periodicity => this%dim, &
                   slab_axis => this%axis, &
                   axis_length => this%size, &
                   origin => this%origin, &
                   axis => this%x)
            !
            !----------------------------------------------------------------------------
            ! Check correction availability
            !
            IF (env_periodicity == 0 .AND. .NOT. cell%cubic) &
                CALL io%error(sub_name, &
                              'Parabolic correction in 0D is only for cubic cells', 1)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL local%init(cell)
            !
            vperiodic => local%of_r
            !
            !----------------------------------------------------------------------------
            !
            CALL charges%multipoles(origin, charge, dipole, quadrupole)
            ! compute multipoles of the system with respect to the chosen origin
            !
            !----------------------------------------------------------------------------
            ! Compute quadratic PBC correction
            !
            fact = e2 * tpi / omega
            !
            SELECT CASE (env_periodicity)
                !
            CASE (0)
                !
                const = madelung(1) * charge * e2 / omega**(1 / 3.D0) - &
                        fact * SUM(quadrupole(:)) / 3.D0
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
                !
            CASE (1)
                CALL io%error(sub_name, 'Option not yet implemented', 1)
                !
            CASE (2)
                !
                const = -pi / 3.D0 * charge / axis_length * e2 - &
                        fact * quadrupole(slab_axis)
                !
                vperiodic(:) = -charge * axis(1, :)**2 + &
                               2.D0 * dipole(slab_axis) * axis(1, :)
                !
                vperiodic = fact * vperiodic + const
                !
            CASE (3)
                const = 0.D0
                vperiodic = 0.D0
                !
            CASE DEFAULT
                !
                CALL io%error(sub_name, &
                              'Unexpected option in dimensionality of PBC correction', 1)
                !
            END SELECT
            !
            potential%of_r = potential%of_r + vperiodic
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL local%destroy()
            !
        END ASSOCIATE
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_1da_vperiodic
    !------------------------------------------------------------------------------------
    !>
    !! Note that in this subroutine the ionic density is implicit (and thus
    !! spread gaussians). Since the gradient of the corrective potential does
    !! not depend on the quadrupole moment of rhotot, it should be independent
    !! on the shape of the ionic density
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_1da_gradvperiodic(this, charges, gvtot)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        CLASS(core_1da_electrostatics), TARGET, INTENT(INOUT) :: this
        TYPE(environ_gradient), INTENT(INOUT) :: gvtot
        !
        REAL(DP), POINTER :: gvperiodic(:, :)
        !
        INTEGER :: icor
        !
        REAL(DP) :: fact
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        !
        TYPE(environ_gradient), TARGET :: glocal
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_1da_gradvperiodic'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(gvtot%cell, charges%cell)) &
            CALL io%error(sub_name, 'Mismatch in domains of gradient and charges', 1)
        !
        IF (.NOT. ASSOCIATED(gvtot%cell, this%cell)) &
            CALL io%error(sub_name, 'Mismatch in domains of gradient and solver', 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => gvtot%cell, &
                   omega => gvtot%cell%omega, &
                   env_periodicity => this%dim, &
                   slab_axis => this%axis, &
                   origin => this%origin, &
                   axis => this%x)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL glocal%init(cell)
            !
            gvperiodic => glocal%of_r
            !
            !----------------------------------------------------------------------------
            !
            CALL charges%multipoles(origin, charge, dipole, quadrupole)
            ! compute multipoles of the system with respect to the chosen origin
            !
            !----------------------------------------------------------------------------
            ! Compute gradient of periodic images correction
            !
            fact = e2 * fpi / omega
            !
            SELECT CASE (env_periodicity)
                !
            CASE (0)
                !
                DO icor = 1, 3
                    gvperiodic(icor, :) = (dipole(icor) - charge * axis(icor, :)) / 3.D0
                END DO
                !
            CASE (1)
                CALL io%error(sub_name, 'Option not yet implemented', 1)
                !
            CASE (2)
                gvperiodic(slab_axis, :) = dipole(slab_axis) - charge * axis(1, :)
                !
            CASE (3)
                gvperiodic = 0.D0
                !
            CASE DEFAULT
                CALL io%error(sub_name, 'Unexpected option', 1)
                !
            END SELECT
            !
            gvperiodic = gvperiodic * fact
            !
            gvtot%of_r = gvtot%of_r + gvperiodic
            ! sum the periodic contribution to the total gradient of the potential
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL glocal%destroy()
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_1da_gradvperiodic
    !------------------------------------------------------------------------------------
    !>
    !! Computes the contribution to the atomic forces
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_1da_fperiodic(this, natoms, ions, auxiliary, f)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: natoms
        CLASS(environ_function), TARGET, INTENT(IN) :: ions(:)
        TYPE(environ_density), INTENT(IN) :: auxiliary
        !
        CLASS(core_1da_electrostatics), TARGET, INTENT(INOUT) :: this
        REAL(DP), INTENT(INOUT) :: f(3, natoms)
        !
        REAL(DP), POINTER :: Z
        !
        INTEGER :: i
        !
        REAL(DP) :: fact, pos(3)
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        REAL(DP) :: ftmp(3, natoms)
        !
        TYPE(environ_density) :: local
        !
        TYPE(environ_function_gaussian), POINTER :: local_ions(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_1da_fperiodic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(auxiliary%cell, this%cell)) &
            CALL io%error(sub_name, 'Mismatch in domains of charges and solver', 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => auxiliary%cell, &
                   omega => auxiliary%cell%omega, &
                   env_periodicity => this%dim, &
                   slab_axis => this%axis, &
                   origin => this%origin)
            !
            !----------------------------------------------------------------------------
            !
            SELECT TYPE (ions)
                !
            TYPE IS (environ_function_gaussian)
                local_ions => ions
                !
            CLASS DEFAULT
                CALL io%error(sub_name, "Unexpected function type", 1)
                !
            END SELECT
            !
            IF (natoms /= SIZE(local_ions)) &
                CALL io%error(sub_name, &
                              'Mismatch between input and stored number of ions', 1)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL local%init(cell)
            !
            local%of_r = auxiliary%of_r
            !
            !----------------------------------------------------------------------------
            !
            CALL local%multipoles(origin, charge, dipole, quadrupole)
            ! compute multipoles of the system with respect to the chosen origin
            !
            !----------------------------------------------------------------------------
            ! Interatomic forces, quadrupole is not needed, thus the same
            ! expression holds for point-like and gaussian nuclei
            !
            fact = e2 * fpi / omega
            ftmp = 0.D0
            !
            DO i = 1, natoms
                pos = local_ions(i)%pos - origin
                Z => local_ions(i)%volume
                !
                SELECT CASE (env_periodicity)
                    !
                CASE (0)
                    ftmp(:, i) = (charge * pos - dipole) / 3.D0
                    !
                CASE (1)
                    CALL io%error(sub_name, 'Option not yet implemented', 1)
                    !
                CASE (2)
                    ftmp(slab_axis, i) = charge * pos(slab_axis) - dipole(slab_axis)
                    !
                CASE (3)
                    ftmp = 0.D0
                    !
                CASE DEFAULT
                    CALL io%error(sub_name, 'Unexpected', 1)
                    !
                END SELECT
                !
                ftmp(:, i) = ftmp(:, i) * fact * Z
            END DO
            !
            f = f + ftmp
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL local%destroy()
            !
        END ASSOCIATE
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_1da_fperiodic
    !------------------------------------------------------------------------------------
    !>
    !! Function that calculates the Gouy-Chapman correction and adds it to the potential
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_1da_vgcs(this, electrolyte, charges, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_electrolyte_base), TARGET, INTENT(IN) :: electrolyte
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        CLASS(core_1da_electrostatics), TARGET, INTENT(INOUT) :: this
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        REAL(DP), POINTER :: v(:)
        !
        INTEGER :: i, icount
        !
        REAL(DP) :: kbt, invkbt
        REAL(DP) :: ez, fact, vstern, const, constl, constr
        REAL(DP) :: dv, vbound, zion
        REAL(DP) :: arg, asinh, coth, acoth
        REAL(DP) :: f1, f2
        REAL(DP) :: area, vtmp
        REAL(DP) :: lin_k, lin_e, lin_c
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        !
        TYPE(environ_density), TARGET :: local
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_1da_vgcs'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(potential%cell, charges%cell)) &
            CALL io%error(sub_name, 'Mismatch in domains of potential and charges', 1)
        !
        IF (.NOT. ASSOCIATED(potential%cell, this%cell)) &
            CALL io%error(sub_name, 'Mismatch in domains of potential and solver', 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => potential%cell, &
                   nnr => potential%cell%nnr, &
                   omega => potential%cell%omega, &
                   env_periodicity => this%dim, &
                   slab_axis => this%axis, &
                   axis_length => this%size, &
                   origin => this%origin, &
                   axis => this%x, &
                   cion => electrolyte%ioncctype(1)%cbulk, &
                   permittivity => electrolyte%permittivity, &
                   xstern => electrolyte%distance)
            !
            !----------------------------------------------------------------------------
            ! Get parameters of electrolyte to compute analytic correction
            !
            IF (electrolyte%ntyp /= 2) &
                CALL io%error(sub_name, &
                              'Unexpected number of counterionic species, &
                              &different from two', 1)
            !
            zion = ABS(electrolyte%ioncctype(1)%z)
            !
            !----------------------------------------------------------------------------
            ! Set Boltzmann factors
            !
            kbt = electrolyte%temperature * K_BOLTZMANN_RY
            invkbt = 1.D0 / kbt
            !
            IF (env_periodicity /= 2) &
                CALL io%error(sub_name, &
                              'Option not yet implemented: 1D Poisson-Boltzmann solver &
                              &only for 2D systems', 1)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL local%init(cell)
            !
            v => local%of_r
            area = omega / axis_length
            !
            !----------------------------------------------------------------------------

            CALL charges%multipoles(origin, charge, dipole, quadrupole)
            ! compute multipoles of the system with respect to the chosen origin
            !
            !----------------------------------------------------------------------------
            ! First apply parabolic correction
            !
            fact = e2 * tpi / omega
            const = -pi / 3.D0 * charge / axis_length * e2 - fact * quadrupole(slab_axis)
            v(:) = -charge * axis(1, :)**2 + 2.D0 * dipole(slab_axis) * axis(1, :)
            v = fact * v + const
            !
            !----------------------------------------------------------------------------
            ! Compute the physical properties of the interface
            !
            ez = -tpi * e2 * charge / area
            ! the input charge density includes explicit and polarization charges, so
            ! tot_charge already accounts for the dielectric screening. permittivity
            ! needs not to be included
            !
            fact = -e2 * SQRT(8.D0 * fpi * cion * kbt / e2 / permittivity)
            arg = ez / fact
            asinh = LOG(arg + SQRT(arg**2 + 1))
            vstern = 2.D0 * kbt / zion * asinh
            dv = 2.D0 * fpi * dipole(slab_axis) / area
            !
            !----------------------------------------------------------------------------
            ! Compute left/right conditions for GCS potential
            !
            arg = (vstern - dv * 0.5D0) * 0.25D0 * invkbt * zion
            coth = (EXP(2.D0 * arg) + 1.D0) / (EXP(2.D0 * arg) - 1.D0)
            constl = coth * EXP(zion * fact * invkbt * 0.5D0 * xstern)
            arg = (vstern + dv * 0.5D0) * 0.25D0 * invkbt * zion
            coth = (EXP(2.D0 * arg) + 1.D0) / (EXP(2.D0 * arg) - 1.D0)
            constr = coth * EXP(zion * fact * invkbt * 0.5D0 * xstern)
            !
            !----------------------------------------------------------------------------
            ! Compute linearized quantities
            ! note that this is simplified due to same c_ion
            !
            lin_k = SQRT(electrolyte%k2)
            lin_e = SQRT(permittivity)
            lin_c = -1.D0 * ez * lin_e / lin_k * EXP(lin_k * xstern / lin_e)
            !
            IF (electrolyte%linearized) &
                vstern = lin_c * EXP(-1.D0 * lin_k * xstern / lin_e)
            !
            !----------------------------------------------------------------------------
            ! Compute value of the reference potential at the boundary with electrolyte
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
            v = v - vbound + vstern
            !
            !----------------------------------------------------------------------------
            ! Compute some constants needed for the calculation
            !
            f1 = -fact * zion * invkbt * 0.5D0
            f2 = 4.D0 * kbt / zion
            !
            !----------------------------------------------------------------------------
            ! Compute the analytic potential
            !
            IF (electrolyte%linearized) THEN
                !
                DO i = 1, nnr
                    !
                    IF (ABS(axis(1, i)) >= xstern) THEN
                        !
                        vtmp = lin_c * EXP(-1.D0 * lin_k * ABS(axis(1, i)) / lin_e)
                        ! linearized Gouy-Chapmann-Stern analytic solution on the outside
                        !
                        v(i) = vtmp - potential%of_r(i)
                        ! remove source potential and add analytic one
                        !
                    END IF
                    !
                END DO
                !
            ELSE
                !
                DO i = 1, nnr
                    !
                    !--------------------------------------------------------------------
                    ! Gouy-Chapmann-Stern analytic solution on the outside
                    !
                    IF (axis(1, i) <= -xstern) THEN
                        !
                        !----------------------------------------------------------------
                        ! Left solution
                        !
                        arg = constl * EXP(ABS(axis(1, i)) * f1)
                        !
                        IF (ABS(arg) > 1.D0) THEN
                            acoth = 0.5D0 * LOG((arg + 1.D0) / (arg - 1.D0))
                        ELSE
                            acoth = 0.D0
                        END IF
                        !
                        vtmp = f2 * acoth
                        !
                        v(i) = vtmp - potential%of_r(i)
                        ! remove source potential and add analytic one
                        !
                    ELSE IF (axis(1, i) >= xstern) THEN
                        !
                        !----------------------------------------------------------------
                        ! Right solution
                        !
                        arg = constr * EXP(ABS(axis(1, i)) * f1)
                        !
                        IF (ABS(arg) > 1.D0) THEN
                            acoth = 0.5D0 * LOG((arg + 1.D0) / (arg - 1.D0))
                        ELSE
                            acoth = 0.D0
                        END IF
                        !
                        vtmp = f2 * acoth
                        !
                        v(i) = vtmp - potential%of_r(i)
                        ! remove source potential and add analytic one
                        !
                    END IF
                    !
                END DO
                !
            END IF
            !
            potential%of_r = potential%of_r + v
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL local%destroy()
            !
        END ASSOCIATE
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_1da_vgcs
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_1da_gradvgcs(this, electrolyte, charges, gradv)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_electrolyte_base), TARGET, INTENT(IN) :: electrolyte
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        CLASS(core_1da_electrostatics), TARGET, INTENT(INOUT) :: this
        TYPE(environ_gradient), INTENT(INOUT) :: gradv
        !
        REAL(DP), POINTER :: gvstern(:, :)
        !
        INTEGER :: i
        !
        REAL(DP) :: zion, kbt, invkbt
        REAL(DP) :: ez, fact, vstern, const, constl, constr
        REAL(DP) :: arg, asinh, coth, acoth
        REAL(DP) :: lin_k, lin_e, lin_c
        REAL(DP) :: f1, f2
        REAL(DP) :: area, dvtmp_dx
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        !
        TYPE(environ_gradient), TARGET :: glocal
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_1da_gradvgcs'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(gradv%cell, charges%cell)) &
            CALL io%error(sub_name, 'Mismatch in domains of potential and charges', 1)
        !
        IF (.NOT. ASSOCIATED(gradv%cell, this%cell)) &
            CALL io%error(sub_name, 'Mismatch in domains of potential and solver', 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => gradv%cell, &
                   nnr => gradv%cell%nnr, &
                   omega => gradv%cell%omega, &
                   env_periodicity => this%dim, &
                   slab_axis => this%axis, &
                   axis_length => this%size, &
                   origin => this%origin, &
                   axis => this%x, &
                   cion => electrolyte%ioncctype(1)%cbulk, &
                   permittivity => electrolyte%permittivity, &
                   xstern => electrolyte%distance)
            !
            !----------------------------------------------------------------------------
            ! Get parameters of electrolyte to compute analytic correction
            !
            IF (electrolyte%ntyp /= 2) &
                CALL io%error(sub_name, &
                              'Unexpected number of counterionic species, &
                              &different from two', 1)
            !
            zion = ABS(electrolyte%ioncctype(1)%z)
            !
            !----------------------------------------------------------------------------
            ! Set Boltzmann factors
            !
            kbt = electrolyte%temperature * K_BOLTZMANN_RY
            invkbt = 1.D0 / kbt
            !
            IF (env_periodicity /= 2) &
                CALL io%error(sub_name, &
                              'Option not yet implemented: 1D Poisson-Boltzmann solver &
                              &only for 2D systems', 1)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL glocal%init(cell)
            !
            gvstern => glocal%of_r
            area = omega / axis_length
            !
            !----------------------------------------------------------------------------
            !
            CALL charges%multipoles(origin, charge, dipole, quadrupole)
            ! compute multipoles of the system with respect to the chosen origin
            !
            !----------------------------------------------------------------------------
            ! First compute the gradient of parabolic correction
            !
            fact = e2 * fpi / omega
            gvstern(slab_axis, :) = dipole(slab_axis) - charge * axis(1, :)
            gvstern = gvstern * fact
            !
            !----------------------------------------------------------------------------
            ! Compute the physical properties of the interface
            !
            ez = -tpi * e2 * charge / area !/ permittivity
            ! the input charge density includes explicit and polarization charges, so
            ! tot_charge already accounts for the dielectric screening. permittivity
            ! needs not to be included
            !
            fact = -e2 * SQRT(8.D0 * fpi * cion * kbt / e2 / permittivity)
            arg = ez / fact
            asinh = LOG(arg + SQRT(arg**2 + 1))
            vstern = 2.D0 * kbt / zion * asinh
            arg = vstern * 0.25D0 * invkbt * zion
            coth = (EXP(2.D0 * arg) + 1.D0) / (EXP(2.D0 * arg) - 1.D0)
            const = coth * EXP(zion * fact * invkbt * 0.5D0 * xstern)
            !
            !----------------------------------------------------------------------------
            ! Compute linearized quantities
            ! Note that this is simplified due to same c_ion
            !
            lin_k = SQRT(electrolyte%k2)
            lin_e = SQRT(permittivity)
            lin_c = -1.D0 * ez * lin_e / lin_k * EXP(lin_k * xstern / lin_e)
            !
            !----------------------------------------------------------------------------
            ! Compute the analytic gradient of potential
            ! Note that the only contribution different from the parabolic
            ! correction is in the region of the diffuse layer
            !
            IF (electrolyte%linearized) THEN
                !
                !------------------------------------------------------------------------
                ! Compute some constants needed for the calculation
                !
                f1 = -1.D0 * lin_k / lin_e
                !
                DO i = 1, nnr
                    !
                    IF (ABS(axis(1, i)) >= xstern) THEN
                        !
                        !----------------------------------------------------------------
                        ! Linearized Gouy-Chapmann-Stern analytic solution on the outside
                        !
                        arg = f1 * ABS(axis(1, i))
                        dvtmp_dx = lin_c * f1 * EXP(arg)
                        !
                        ! remove source potential and add analytic one
                        gvstern(slab_axis, i) = &
                            -gradv%of_r(slab_axis, i) + &
                            (dvtmp_dx - ez) * ABS(axis(1, i)) / axis(1, i)
                        !
                    END IF
                    !
                END DO
                !
            ELSE
                !
                !------------------------------------------------------------------------
                ! Compute some constants needed for the calculation
                !
                f1 = -fact * zion * invkbt * 0.5D0
                f2 = 4.D0 * kbt / zion
                !
                DO i = 1, nnr
                    !
                    IF (axis(1, i) <= -xstern) THEN
                        !
                        !----------------------------------------------------------------
                        ! Gouy-Chapmann-Stern analytic solution on the outside
                        !
                        arg = constl * EXP(ABS(axis(1, i)) * f1)
                        dvtmp_dx = f1 * f2 * arg / (1.D0 - arg**2)
                        !
                        ! remove source potential (linear) and add analytic one
                        gvstern(slab_axis, i) = &
                            -gradv%of_r(slab_axis, i) + &
                            (dvtmp_dx - ez) * ABS(axis(1, i)) / axis(1, i)
                        !
                    ELSE IF (axis(1, i) >= xstern) THEN
                        !
                        !----------------------------------------------------------------
                        ! Gouy-Chapmann-Stern analytic solution on the outside
                        !
                        arg = constr * EXP(ABS(axis(1, i)) * f1)
                        dvtmp_dx = f1 * f2 * arg / (1.D0 - arg**2)
                        !
                        ! remove source potential (linear) and add analytic one
                        gvstern(slab_axis, i) = &
                            -gradv%of_r(slab_axis, i) + &
                            (dvtmp_dx - ez) * ABS(axis(1, i)) / axis(1, i)
                        !
                    END IF
                    !
                END DO
                !
            END IF
            !
            gradv%of_r = gradv%of_r + gvstern
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL glocal%destroy()
            !
        END ASSOCIATE
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_1da_gradvgcs
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_1da_vms(this, semiconductor, charges, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_semiconductor_base), TARGET, INTENT(IN) :: semiconductor
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        CLASS(core_1da_electrostatics), TARGET, INTENT(INOUT) :: this
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        REAL(DP), POINTER :: v(:)
        !
        INTEGER :: i, icount
        !
        REAL(DP) :: kbt, invkbt
        REAL(DP) :: ez, fact, vms
        REAL(DP) :: arg, const, depletion_length
        REAL(DP) :: area, vtmp, vbound, distance
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        !
        TYPE(environ_density), TARGET :: local
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_1da_vms'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(potential%cell, charges%cell)) &
            CALL io%error(sub_name, 'Mismatch in domains of potential and charges', 1)
        !
        IF (potential%cell%nnr /= this%nnr) &
            CALL io%error(sub_name, 'Mismatch in domains of potential and solver', 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => potential%cell, &
                   nnr => potential%cell%nnr, &
                   omega => potential%cell%omega, &
                   env_periodicity => this%dim, &
                   slab_axis => this%axis, &
                   axis_length => this%size, &
                   origin => this%origin, &
                   axis => this%x, &
                   permittivity => semiconductor%permittivity, &
                   carrier_density => semiconductor%carrier_density, &
                   xstern => semiconductor%sc_distance)
            !
            !----------------------------------------------------------------------------
            ! Set Boltzmann factors
            !
            kbt = semiconductor%temperature * K_BOLTZMANN_RY
            invkbt = 1.D0 / kbt
            !
            IF (env_periodicity /= 2) &
                CALL io%error(sub_name, &
                              'Option not yet implemented: 1D Poisson-Boltzmann solver &
                              &only for 2D systems', 1)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL local%init(cell)
            !
            v => local%of_r
            area = omega / axis_length
            !
            !----------------------------------------------------------------------------
            !
            CALL charges%multipoles(origin, charge, dipole, quadrupole)
            ! compute multipoles of the system w.r.t the chosen origin
            !
            !----------------------------------------------------------------------------
            ! First apply parabolic correction
            !
            fact = e2 * tpi / omega
            const = -pi / 3.D0 * charge / axis_length * e2 - fact * quadrupole(slab_axis)
            !
            v = -charge * axis(1, :)**2 + 2.D0 * dipole(slab_axis) * axis(1, :)
            v = fact * v + const
            !
            !----------------------------------------------------------------------------
            ! Compute the physical properties of the interface
            !
            ez = -tpi * e2 * charge / area
            fact = 1.D0 / tpi / e2 / 2.D0 / carrier_density
            arg = fact * (ez**2.D0)
            vms = arg ! +kbt ! #TODO figure this out
            !
            !----------------------------------------------------------------------------
            ! Finds the total length of the depletion region
            !
            depletion_length = 2.D0 * fact * ez
            !
            IF (io%lnode .AND. io%verbosity > 0) THEN
                WRITE (io%debug_unit, *) "depletion length: ", depletion_length
                WRITE (io%debug_unit, *) "vms: ", vms
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
            !----------------------------------------------------------------------------
            ! Compute the analytic potential and charge
            !
            v = v - vbound - vms
            !
            DO i = 1, nnr
                !
                IF (ABS(axis(1, i)) >= xstern) THEN
                    distance = ABS(axis(1, i)) - xstern
                    !
                    !--------------------------------------------------------------------
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
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL local%destroy()
            !
        END ASSOCIATE
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_1da_vms
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_1da_gradvms(this, semiconductor, charges, gradv)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_semiconductor_base), TARGET, INTENT(IN) :: semiconductor
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        CLASS(core_1da_electrostatics), TARGET, INTENT(INOUT) :: this
        TYPE(environ_gradient), INTENT(INOUT) :: gradv
        !
        REAL(DP), POINTER :: gvms(:, :)
        !
        INTEGER :: i
        !
        REAL(DP) :: kbt, invkbt
        REAL(DP) :: ez, fact, vms
        REAL(DP) :: arg
        REAL(DP) :: area, dvtmp_dx
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        !
        TYPE(environ_gradient), TARGET :: glocal
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_1da_gradvms'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(gradv%cell, charges%cell)) &
            CALL io%error(sub_name, 'Mismatch in domains of potential and charges', 1)
        !
        IF (gradv%cell%nnr /= this%nnr) &
            CALL io%error(sub_name, 'Mismatch in domains of potential and solver', 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => gradv%cell, &
                   omega => gradv%cell%omega, &
                   env_periodicity => this%dim, &
                   slab_axis => this%axis, &
                   axis_length => this%size, &
                   origin => this%origin, &
                   axis => this%x, &
                   permittivity => semiconductor%permittivity, &
                   carrier_density => semiconductor%carrier_density, &
                   xstern => semiconductor%sc_distance)
            !
            !----------------------------------------------------------------------------
            ! Set Boltzmann factors
            !
            kbt = semiconductor%temperature * K_BOLTZMANN_RY
            invkbt = 1.D0 / kbt
            !
            IF (env_periodicity /= 2) &
                CALL io%error(sub_name, &
                              'Option not yet implemented: 1D Poisson-Boltzmann solver &
                              &only for 2D systems', 1)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL glocal%init(cell)
            !
            gvms => glocal%of_r
            !
            !----------------------------------------------------------------------------
            !
            CALL charges%multipoles(origin, charge, dipole, quadrupole)
            ! compute multipoles of the system w.r.t the chosen origin
            !
            !----------------------------------------------------------------------------
            ! First compute the gradient of parabolic correction
            !
            fact = e2 * fpi / omega
            gvms(slab_axis, :) = dipole(slab_axis) - charge * axis(1, :)
            gvms = gvms * fact
            !
            gradv%of_r = gradv%of_r + gvms
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL glocal%destroy()
            !
        END ASSOCIATE
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_1da_gradvms
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_core_1da_electrostatics
!----------------------------------------------------------------------------------------
