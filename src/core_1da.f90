!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 3.0
!
!     Environ 3.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 3.0 is distributed in the hope that it will be useful,
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
MODULE class_core_1da
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_mp, ONLY: env_mp_sum
    !
    USE environ_param, ONLY: DP, e2, K_BOLTZMANN_RY, pi, tpi, fpi, madelung
    !
    USE class_cell
    USE class_density
    USE class_gradient
    USE class_functions
    !
    USE class_core
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
    TYPE, EXTENDS(environ_core), PUBLIC :: core_1da
        !--------------------------------------------------------------------------------
        !
        INTEGER :: nnr = 0
        INTEGER :: dim = 0
        INTEGER :: pdim = 0
        INTEGER :: axis = 0
        !
        REAL(DP) :: size = 0.D0
        REAL(DP) :: origin(3) = 0.D0 ! set later to system center of mass
        !
        REAL(DP), ALLOCATABLE :: x(:, :)
        !
        !--------------------------------------------------------------------------------
        ! Reduced arrays for optimization
        !
        INTEGER, ALLOCATABLE :: ir(:) ! indices of points of interest
        REAL(DP), ALLOCATABLE :: r(:, :) ! displacements w.r.t. origin
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_core_1da
        PROCEDURE :: init => init_core_1da
        PROCEDURE :: update_cell => update_core_1da_cell
        PROCEDURE :: update_origin => update_core_1da_origin
        PROCEDURE :: destroy => destroy_core_1da
        !
        PROCEDURE :: calc_vperiodic => calc_1da_vperiodic
        PROCEDURE :: calc_grad_vperiodic => calc_1da_grad_vperiodic
        PROCEDURE :: calc_fperiodic => calc_1da_fperiodic
        !
        PROCEDURE :: calc_vgcs => calc_1da_vgcs
        PROCEDURE :: calc_grad_vgcs => calc_1da_grad_vgcs
        !
        PROCEDURE :: calc_vms => calc_1da_vms
        PROCEDURE :: calc_grad_vms => calc_1da_grad_vms
        !
        PROCEDURE :: calc_vms_gcs => calc_1da_vms_gcs
        PROCEDURE :: calc_grad_vms_gcs => calc_1da_grad_vms_gcs

        !--------------------------------------------------------------------------------
    END TYPE core_1da
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   ADMIN METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_core_1da(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_1da), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_core_1da'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%cell)) CALL io%create_error(routine)
        !
        IF (ALLOCATED(this%x)) CALL io%create_error(routine)
        !
        IF (ALLOCATED(this%ir)) CALL io%create_error(routine)
        !
        IF (ALLOCATED(this%r)) CALL io%create_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        this%core_type = '1d-analytic'
        this%nnr = 0
        this%dim = 0
        this%pdim = 0
        this%axis = 0
        this%size = 0.D0
        this%origin = 0.D0
        !
        NULLIFY (this%cell)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_core_1da
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_core_1da(this, dim, axis, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, axis
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        CLASS(core_1da), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'init_core_1da'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        IF (dim == 3 .OR. dim < 0) &
            CALL io%error(routine, &
                          "Wrong dimensions for analytic one dimensional core", 1)
        !
        this%dim = dim
        this%pdim = 3 - dim
        !
        IF ((dim == 1 .OR. dim == 2) .AND. (axis > 3 .OR. axis < 1)) &
            CALL io%error(routine, &
                          "Wrong choice of axis for analytic one dimensional core", 1)
        !
        this%axis = axis
        this%nnr = cell%nnr
        this%cell => cell
        ALLOCATE (this%x(this%pdim, this%nnr))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_core_1da
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_core_1da_cell(this, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        CLASS(core_1da), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%cell => cell
        !
        IF (this%dim == 0) THEN
            this%size = cell%omega
        ELSE IF (this%dim == 1) THEN
            this%size = cell%omega / cell%at(this%axis, this%axis)
        ELSE IF (this%dim == 2) THEN
            this%size = cell%at(this%axis, this%axis)
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_core_1da_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_core_1da_origin(this, origin)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: origin(3)
        !
        CLASS(core_1da), INTENT(INOUT) :: this
        !
        LOGICAL :: physical
        INTEGER :: i, count
        REAL(DP) :: r(3), r2
        !
        INTEGER, ALLOCATABLE :: ir(:)
        REAL(DP), ALLOCATABLE :: disps(:, :)
        !
        CHARACTER(LEN=80) :: routine = 'update_core_1da_origin'
        !
        !--------------------------------------------------------------------------------
        !
        this%origin = origin
        !
        ALLOCATE (ir(this%cell%ir_end))
        ALLOCATE (disps(3, this%cell%ir_end))
        ir = 0
        disps = 0.D0
        count = 0
        !
        ASSOCIATE (cell => this%cell, &
                   dim => this%dim)
            !
            IF (dim == 0) THEN
                !
                DO i = 1, cell%ir_end
                    !
                    CALL cell%get_min_distance(i, 0, 0, origin, r, r2, physical)
                    !
                    IF (.NOT. physical) CYCLE
                    !
                    count = count + 1
                    ir(count) = i
                    disps(:, count) = r
                    !
                    this%x(:, i) = r
                END DO
                !
            ELSE IF (dim == 1) THEN
                CALL io%error(routine, "Option not yet implemented", 1)
            ELSE IF (dim == 2) THEN
                !
                DO i = 1, cell%ir_end
                    !
                    CALL cell%get_min_distance(i, 0, 0, origin, r, r2, physical)
                    !
                    IF (.NOT. physical) CYCLE
                    !
                    count = count + 1
                    ir(count) = i
                    disps(:, count) = r
                    !
                    this%x(1, i) = r(this%axis)
                END DO
                !
            END IF
            !
        END ASSOCIATE
        !
        IF (.NOT. ALLOCATED(this%ir)) ALLOCATE (this%ir(count - 2))
        !
        IF (.NOT. ALLOCATED(this%r)) ALLOCATE (this%r(3, count - 2))
        !
        this%ir = ir(:count)
        this%r = disps(:, :count)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_core_1da_origin
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_core_1da(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_1da), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'destroy_core_1da'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%cell)) RETURN
        !
        IF (.NOT. ALLOCATED(this%x)) CALL io%destroy_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        DEALLOCATE (this%x)
        !
        IF (ALLOCATED(this%ir)) DEALLOCATE (this%ir)
        !
        IF (ALLOCATED(this%r)) DEALLOCATE (this%r)
        !
        NULLIFY (this%cell)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_core_1da
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
    SUBROUTINE calc_1da_vperiodic(this, charges, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_1da), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        !
        REAL(DP), POINTER :: vperiodic(:)
        !
        INTEGER :: i
        !
        REAL(DP) :: fact
        REAL(DP) :: const, charge, dipole(3), quadrupole(3)
        !
        TYPE(environ_density), TARGET :: local
        !
        CHARACTER(LEN=80) :: routine = 'calc_1da_vperiodic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(routine)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(v%cell, charges%cell)) &
            CALL io%error(routine, "Mismatch in domains of potential and charges", 1)
        !
        IF (.NOT. ASSOCIATED(v%cell, this%cell)) &
            CALL io%error(routine, "Mismatch in domains of potential and solver", 1)
        !
        IF (this%dim == 0 .AND. .NOT. v%cell%cubic) &
            CALL io%error(routine, &
                          "Parabolic correction in 0D is only for cubic cells", 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (omega => v%cell%omega, &
                   slab_axis => this%axis, &
                   axis => this%x)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL local%init(v%cell)
            !
            vperiodic => local%of_r
            !
            !----------------------------------------------------------------------------
            ! Compute multipoles of the system with respect to the chosen origin
            !
            CALL charges%multipoles(this%origin, charge, dipole, quadrupole, this%ir, &
                                    this%r)
            !
            !----------------------------------------------------------------------------
            ! Compute quadratic PBC correction
            !
            fact = e2 * tpi / omega
            !
            SELECT CASE (this%dim)
                !
            CASE (0)
                !
                const = madelung(1) * charge * e2 / omega**(1 / 3.D0) - &
                        fact * SUM(quadrupole) / 3.D0
                !
                vperiodic = 0.D0
                !
                DO i = 1, 3
                    !
                    vperiodic = vperiodic - charge * axis(i, :)**2 + &
                                2.D0 * dipole(i) * axis(i, :)
                    !
                END DO
                !
                vperiodic = fact / 3.D0 * vperiodic + const
                !
            CASE (1)
                CALL io%error(routine, "Option not yet implemented", 1)
                !
            CASE (2)
                !
                const = -pi / 3.D0 * charge / this%size * e2 - &
                        fact * quadrupole(slab_axis)
                !
                vperiodic = -charge * axis(1, :)**2 + &
                            2.D0 * dipole(slab_axis) * axis(1, :)
                !
                vperiodic = fact * vperiodic + const
                !
            CASE (3)
                const = 0.D0
                vperiodic = 0.D0
                !
            CASE DEFAULT
                CALL io%error(routine, "Unexpected system dimensions", 1)
                !
            END SELECT
            !
            v%of_r = v%of_r + vperiodic
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL local%destroy()
            !
        END ASSOCIATE
        !
        CALL env_stop_clock(routine)
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
    SUBROUTINE calc_1da_grad_vperiodic(this, charges, grad_v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_1da), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad_v
        !
        REAL(DP), POINTER :: gvperiodic(:, :)
        !
        INTEGER :: i
        !
        REAL(DP) :: fact
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        !
        TYPE(environ_gradient), TARGET :: glocal
        !
        CHARACTER(LEN=80) :: routine = 'calc_1da_grad_vperiodic'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(grad_v%cell, charges%cell)) &
            CALL io%error(routine, "Mismatch in domains of gradient and charges", 1)
        !
        IF (.NOT. ASSOCIATED(grad_v%cell, this%cell)) &
            CALL io%error(routine, "Mismatch in domains of gradient and solver", 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (slab_axis => this%axis, &
                   axis => this%x)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL glocal%init(grad_v%cell)
            !
            gvperiodic => glocal%of_r
            !
            !----------------------------------------------------------------------------
            ! Compute multipoles of the system with respect to the chosen origin
            !
            CALL charges%multipoles(this%origin, charge, dipole, quadrupole, this%ir, &
                                    this%r)
            !
            !----------------------------------------------------------------------------
            ! Compute gradient of periodic images correction
            !
            fact = e2 * fpi / grad_v%cell%omega
            !
            SELECT CASE (this%dim)
                !
            CASE (0)
                !
                DO i = 1, 3
                    gvperiodic(i, :) = (dipole(i) - charge * axis(i, :)) / 3.D0
                END DO
                !
            CASE (1)
                CALL io%error(routine, "Option not yet implemented", 1)
                !
            CASE (2)
                gvperiodic(slab_axis, :) = dipole(slab_axis) - charge * axis(1, :)
                !
            CASE (3)
                gvperiodic = 0.D0
                !
            CASE DEFAULT
                CALL io%error(routine, "Unexpected system dimensions", 1)
                !
            END SELECT
            !
            gvperiodic = gvperiodic * fact
            !
            grad_v%of_r = grad_v%of_r + gvperiodic
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
    END SUBROUTINE calc_1da_grad_vperiodic
    !------------------------------------------------------------------------------------
    !>
    !! Computes the contribution to the atomic forces
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_1da_fperiodic(this, nat, ions, auxiliary, force)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_1da), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: nat
        TYPE(environ_functions), INTENT(IN) :: ions
        TYPE(environ_density), INTENT(IN) :: auxiliary
        !
        REAL(DP), INTENT(INOUT) :: force(3, nat)
        !
        INTEGER :: i
        !
        REAL(DP) :: fact
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        REAL(DP) :: ftmp(3, nat)
        !
        TYPE(environ_density) :: local
        !
        CHARACTER(LEN=80) :: routine = 'calc_1da_fperiodic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(routine)
        !
        !--------------------------------------------------------------------------------
        !
        IF (nat /= ions%number) &
            CALL io%error(routine, &
                          "Mismatch between input and stored number of ions", 1)
        !
        IF (.NOT. ASSOCIATED(auxiliary%cell, this%cell)) &
            CALL io%error(routine, "Mismatch in domains of charges and solver", 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (slab_axis => this%axis, &
                   origin => this%origin)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL local%init(auxiliary%cell)
            !
            local%of_r = auxiliary%of_r
            !
            !----------------------------------------------------------------------------
            ! Compute multipoles of the system with respect to the chosen origin
            !
            CALL local%multipoles(origin, charge, dipole, quadrupole, this%ir, this%r)
            !
            !----------------------------------------------------------------------------
            ! Interatomic forces, quadrupole is not needed, thus the same
            ! expression holds for point-like and gaussian nuclei
            !
            fact = e2 * fpi / auxiliary%cell%omega
            ftmp = 0.D0
            !
            DO i = 1, nat
                !
                ASSOCIATE (pos => ions%array(i)%pos - origin, &
                           Z => ions%array(i)%volume)
                    !
                    SELECT CASE (this%dim)
                        !
                    CASE (0)
                        ftmp(:, i) = (charge * pos - dipole) / 3.D0
                        !
                    CASE (1)
                        CALL io%error(routine, "Option not yet implemented", 1)
                        !
                    CASE (2)
                        ftmp(slab_axis, i) = charge * pos(slab_axis) - dipole(slab_axis)
                        !
                    CASE (3)
                        ftmp = 0.D0
                        !
                    CASE DEFAULT
                        CALL io%error(routine, "Unexpected system dimensions", 1)
                        !
                    END SELECT
                    !
                    ftmp(:, i) = ftmp(:, i) * fact * Z
                    !
                END ASSOCIATE
                !
            END DO
            !
            force = force + ftmp
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL local%destroy()
            !
        END ASSOCIATE
        !
        CALL env_stop_clock(routine)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_1da_fperiodic
    !------------------------------------------------------------------------------------
    !>
    !! Function that calculates the Gouy-Chapman correction and adds it to the potential
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_1da_vgcs(this, electrolyte, charges, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_1da), INTENT(IN) :: this
        TYPE(environ_electrolyte_base), INTENT(IN) :: electrolyte
        TYPE(environ_density), INTENT(IN) :: charges
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        !
        REAL(DP), POINTER :: vgcs(:)
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
        CHARACTER(LEN=80) :: routine = 'calc_1da_vgcs'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(routine)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(v%cell, charges%cell)) &
            CALL io%error(routine, "Mismatch in domains of potential and charges", 1)
        !
        IF (.NOT. ASSOCIATED(v%cell, this%cell)) &
            CALL io%error(routine, "Mismatch in domains of potential and solver", 1)
        !
        IF (electrolyte%ntyp /= 2) &
            CALL io%error(routine, &
                          "Unexpected number of counterionic species, different from two", 1)
        !
        IF (this%dim /= 2) &
            CALL io%error(routine, &
                          "Option not yet implemented: 1D Poisson-Boltzmann solver only for 2D systems", 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (comm => v%cell%dfft%comm, &
                   nnr => v%cell%nnr, &
                   omega => v%cell%omega, &
                   slab_axis => this%axis, &
                   axis_length => this%size, &
                   axis => this%x, &
                   cion => electrolyte%ioncctype(1)%cbulk, &
                   zion => ABS(electrolyte%ioncctype(1)%z), &
                   permittivity => electrolyte%permittivity, &
                   xstern => electrolyte%distance)
            !
            !----------------------------------------------------------------------------
            ! Set Boltzmann factors
            !
            kbt = electrolyte%temperature * K_BOLTZMANN_RY
            invkbt = 1.D0 / kbt
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL local%init(v%cell)
            !
            vgcs => local%of_r
            area = omega / axis_length
            !
            !----------------------------------------------------------------------------
            ! Compute multipoles of the system with respect to the chosen origin
            !
            CALL charges%multipoles(this%origin, charge, dipole, quadrupole, this%ir, &
                                    this%r)
            !
            !----------------------------------------------------------------------------
            ! First apply parabolic correction
            !
            fact = e2 * tpi / omega
            const = -pi / 3.D0 * charge / axis_length * e2 - fact * quadrupole(slab_axis)
            vgcs = -charge * axis(1, :)**2 + 2.D0 * dipole(slab_axis) * axis(1, :)
            vgcs = fact * vgcs + const
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
                    vbound = vbound + v%of_r(i) + vgcs(i) - &
                             ez * (ABS(axis(1, i)) - xstern)
                    !
                END IF
                !
            END DO
            !
            CALL env_mp_sum(icount, comm)
            !
            CALL env_mp_sum(vbound, comm)
            !
            vbound = vbound / DBLE(icount)
            vgcs = vgcs - vbound + vstern
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
                        vgcs(i) = vtmp - v%of_r(i)
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
                        vgcs(i) = vtmp - v%of_r(i)
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
                        vgcs(i) = vtmp - v%of_r(i)
                        ! remove source potential and add analytic one
                        !
                    END IF
                    !
                END DO
                !
            END IF
            !
            v%of_r = v%of_r + vgcs
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL local%destroy()
            !
        END ASSOCIATE
        !
        CALL env_stop_clock(routine)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_1da_vgcs
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_1da_grad_vgcs(this, electrolyte, charges, grad_v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_1da), INTENT(IN) :: this
        TYPE(environ_electrolyte_base), INTENT(IN) :: electrolyte
        TYPE(environ_density), INTENT(IN) :: charges
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad_v
        !
        REAL(DP), POINTER :: grad_vgcs(:, :)
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
        CHARACTER(LEN=80) :: routine = 'calc_1da_grad_vgcs'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(routine)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(grad_v%cell, charges%cell)) &
            CALL io%error(routine, "Mismatch in domains of potential and charges", 1)
        !
        IF (.NOT. ASSOCIATED(grad_v%cell, this%cell)) &
            CALL io%error(routine, "Mismatch in domains of potential and solver", 1)
        !
        IF (electrolyte%ntyp /= 2) &
            CALL io%error(routine, &
                          "Unexpected number of counterionic species, different from two", 1)
        !
        IF (this%dim /= 2) &
            CALL io%error(routine, &
                          "Option not yet implemented: 1D Poisson-Boltzmann solver only for 2D systems", 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (nnr => grad_v%cell%nnr, &
                   omega => grad_v%cell%omega, &
                   slab_axis => this%axis, &
                   axis => this%x, &
                   cion => electrolyte%ioncctype(1)%cbulk, &
                   zion => ABS(electrolyte%ioncctype(1)%z), &
                   permittivity => electrolyte%permittivity, &
                   xstern => electrolyte%distance)
            !
            !----------------------------------------------------------------------------
            ! Set Boltzmann factors
            !
            kbt = electrolyte%temperature * K_BOLTZMANN_RY
            invkbt = 1.D0 / kbt
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL glocal%init(grad_v%cell)
            !
            grad_vgcs => glocal%of_r
            area = omega / this%size
            !
            !----------------------------------------------------------------------------
            ! Compute multipoles of the system with respect to the chosen origin
            !
            CALL charges%multipoles(this%origin, charge, dipole, quadrupole, this%ir, &
                                    this%r)
            !
            !----------------------------------------------------------------------------
            ! First compute the gradient of parabolic correction
            !
            fact = e2 * fpi / omega
            grad_vgcs(slab_axis, :) = dipole(slab_axis) - charge * axis(1, :)
            grad_vgcs = grad_vgcs * fact
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
                        grad_vgcs(slab_axis, i) = &
                            -grad_v%of_r(slab_axis, i) + &
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
                        grad_vgcs(slab_axis, i) = &
                            -grad_v%of_r(slab_axis, i) + &
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
                        grad_vgcs(slab_axis, i) = &
                            -grad_v%of_r(slab_axis, i) + &
                            (dvtmp_dx - ez) * ABS(axis(1, i)) / axis(1, i)
                        !
                    END IF
                    !
                END DO
                !
            END IF
            !
            grad_v%of_r = grad_v%of_r + grad_vgcs
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL glocal%destroy()
            !
        END ASSOCIATE
        !
        CALL env_stop_clock(routine)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_1da_grad_vgcs
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_1da_vms(this, semiconductor, charges, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_1da), INTENT(IN) :: this
        TYPE(environ_semiconductor_base), INTENT(IN) :: semiconductor
        TYPE(environ_density), INTENT(IN) :: charges
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        !
        REAL(DP), POINTER :: vms(:)
        !
        INTEGER :: i, icount
        !
        REAL(DP) :: kbt, invkbt
        REAL(DP) :: ez, fact, delta_vms
        REAL(DP) :: arg, const, depletion_length
        REAL(DP) :: area, vtmp, vbound, distance
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        !
        TYPE(environ_density), TARGET :: local
        !
        CHARACTER(LEN=80) :: routine = 'calc_1da_vms'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(routine)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(v%cell, charges%cell)) &
            CALL io%error(routine, "Mismatch in domains of potential and charges", 1)
        !
        IF (v%cell%nnr /= this%nnr) &
            CALL io%error(routine, "Mismatch in domains of potential and solver", 1)
        !
        IF (this%dim /= 2) &
            CALL io%error(routine, &
                          "Option not yet implemented: 1D Poisson-Boltzmann solver only for 2D systems", 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (comm => v%cell%dfft%comm, &
                   nnr => v%cell%nnr, &
                   omega => v%cell%omega, &
                   slab_axis => this%axis, &
                   axis_length => this%size, &
                   axis => this%x, &
                   xstern => semiconductor%sc_distance)
            !
            !----------------------------------------------------------------------------
            ! Set Boltzmann factors
            !
            kbt = semiconductor%temperature * K_BOLTZMANN_RY
            invkbt = 1.D0 / kbt
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL local%init(v%cell)
            !
            vms => local%of_r
            area = omega / axis_length
            !
            !----------------------------------------------------------------------------
            ! Compute multipoles of the system w.r.t the chosen origin
            !
            CALL charges%multipoles(this%origin, charge, dipole, quadrupole, this%ir, &
                                    this%r)
            !
            !----------------------------------------------------------------------------
            ! First apply parabolic correction
            !
            fact = e2 * tpi / omega
            const = -pi / 3.D0 * charge / axis_length * e2 - fact * quadrupole(slab_axis)
            !
            vms = -charge * axis(1, :)**2 + 2.D0 * dipole(slab_axis) * axis(1, :)
            vms = fact * vms + const
            !
            !----------------------------------------------------------------------------
            ! Compute the physical properties of the interface
            !
            ez = -tpi * e2 * charge / area
            fact = 1.D0 / tpi / e2 / 2.D0 / semiconductor%carrier_density
            arg = fact * (ez**2.D0)
            delta_vms = arg
            !
            !----------------------------------------------------------------------------
            ! Finds the total length of the depletion region
            !
            depletion_length = 2.D0 * fact * ez
            !
            IF (io%lnode .AND. io%verbosity > 0) THEN
                WRITE (io%debug_unit, *) "depletion length: ", depletion_length
                WRITE (io%debug_unit, *) "delta_vms: ", delta_vms
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
                    vbound = vbound + v%of_r(i) + vms(i) - &
                             ez * (ABS(axis(1, i)) - xstern)
                    !
                END IF
                !
            END DO
            !
            CALL env_mp_sum(icount, comm)
            !
            CALL env_mp_sum(vbound, comm)
            !
            vbound = vbound / DBLE(icount)
            !
            !----------------------------------------------------------------------------
            ! Compute the analytic potential and charge
            !
            vms = vms - vbound - delta_vms
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
                    vms(i) = vms(i) + vtmp - ez * distance
                    ! remove source potential (linear) and add analytic one
                    !
                END IF
                !
            END DO
            !
            v%of_r = v%of_r + vms
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL local%destroy()
            !
        END ASSOCIATE
        !
        CALL env_stop_clock(routine)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_1da_vms
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_1da_grad_vms(this, semiconductor, charges, grad_v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_1da), INTENT(IN) :: this
        TYPE(environ_semiconductor_base), INTENT(IN) :: semiconductor
        TYPE(environ_density), INTENT(IN) :: charges
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad_v
        !
        REAL(DP), POINTER :: grad_vms(:, :)
        !
        INTEGER :: i
        !
        REAL(DP) :: kbt, invkbt
        REAL(DP) :: ez, fact, delta_vms
        REAL(DP) :: arg
        REAL(DP) :: area, dvtmp_dx
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        !
        TYPE(environ_gradient), TARGET :: glocal
        !
        CHARACTER(LEN=80) :: routine = 'calc_1da_grad_vms'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(routine)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(grad_v%cell, charges%cell)) &
            CALL io%error(routine, "Mismatch in domains of potential and charges", 1)
        !
        IF (grad_v%cell%nnr /= this%nnr) &
            CALL io%error(routine, "Mismatch in domains of potential and solver", 1)
        !
        IF (this%dim /= 2) &
            CALL io%error(routine, &
                          "Option not yet implemented: 1D Poisson-Boltzmann solver only for 2D systems", 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (slab_axis => this%axis, &
                   permittivity => semiconductor%permittivity)
            !
            !----------------------------------------------------------------------------
            ! Set Boltzmann factors
            !
            kbt = semiconductor%temperature * K_BOLTZMANN_RY
            invkbt = 1.D0 / kbt
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL glocal%init(grad_v%cell)
            !
            grad_vms => glocal%of_r
            !
            !----------------------------------------------------------------------------
            ! Compute multipoles of the system w.r.t the chosen origin
            !
            CALL charges%multipoles(this%origin, charge, dipole, quadrupole, this%ir, &
                                    this%r)
            !
            !----------------------------------------------------------------------------
            ! First compute the gradient of parabolic correction
            !
            fact = e2 * fpi / grad_v%cell%omega
            grad_vms(slab_axis, :) = dipole(slab_axis) - charge * this%x(1, :)
            grad_vms = grad_vms * fact
            !
            grad_v%of_r = grad_v%of_r + grad_vms
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL glocal%destroy()
            !
        END ASSOCIATE
        !
        CALL env_stop_clock(routine)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_1da_grad_vms
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_1da_vms_gcs(this, electrolyte, semiconductor, charges, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_1da), INTENT(IN) :: this
        TYPE(environ_electrolyte_base), INTENT(IN) :: electrolyte
        TYPE(environ_density), INTENT(IN) :: charges
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_semiconductor_base), INTENT(INOUT) :: semiconductor
        !
        REAL(DP), POINTER :: vms_gcs(:)
        !
        INTEGER :: i, icount, ir, j, k, naxis, z_cut_idx, avg_window
        !
        REAL(DP) :: kbt, invkbt
        REAL(DP) :: ez_ms, fact, vms
        REAL(DP) :: arg, const, depletion_length
        REAL(DP) :: v_cut, v_edge, z_cut, z_val
        REAL(DP) :: asinh, coth, acoth
        REAL(DP) :: max_axis
        REAL(DP) :: area
        REAL(DP) :: charge, dipole(3), quadrupole(3)
        !
        REAL(DP), DIMENSION(:), ALLOCATABLE :: current_pot, subtracted_pot, avgd_pot
        REAL(DP), DIMENSION(:), ALLOCATABLE :: flatband_pot
        !
        LOGICAL :: physical
        !
        TYPE(environ_density), TARGET :: local
        !
        CHARACTER(LEN=80) :: routine = 'calc_1da_vms_gcs'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(routine)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(v%cell, charges%cell)) &
            CALL io%error(routine, 'Missmatch in domains of potential and charges', 1)
        !
        IF (.NOT. ASSOCIATED(v%cell, this%cell)) &
            CALL io%error(routine, "Mismatch in domains of potential and solver", 1)
        !
        IF (electrolyte%ntyp /= 2) &
            CALL io%error(routine, &
                          "Unexpected number of counterionic species, different from two", 1)
        !
        IF (this%dim /= 2) &
            CALL io%error(routine, &
                          "Option not yet implemented: 1D Poisson-Boltzmann solver only for 2D systems", 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (comm => v%cell%dfft%comm, &
                   nnr => v%cell%nnr, &
                   omega => v%cell%omega, &
                   slab_axis => this%axis, &
                   axis_length => this%size, &
                   axis => this%x, &
                   permittivity_ms => semiconductor%permittivity, &
                   carrier_density => semiconductor%carrier_density, &
                   electrode_charge => semiconductor%electrode_charge, &
                   xstern_ms => semiconductor%sc_distance, &
                   slab_charge => semiconductor%slab_charge)
            !
            !----------------------------------------------------------------------------
            ! Set Boltzmann factors
            !
            kbt = semiconductor%temperature * k_boltzmann_ry
            invkbt = 1.D0 / kbt
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL local%init(v%cell)
            !
            vms_gcs => local%of_r
            area = omega / axis_length
            semiconductor%surf_area_per_sq_cm = area * 2.8002D-17
            ! value of 1 square bohr in cm^2
            !
            !----------------------------------------------------------------------------
            !
            CALL charges%multipoles(this%origin, charge, dipole, quadrupole)
            ! compute multipoles of the system w.r.t the chosen origin
            !
            !----------------------------------------------------------------------------
            ! First apply parabolic correction
            !
            fact = e2 * tpi / omega
            const = -pi / 3.D0 * charge / axis_length * e2 - fact * quadrupole(slab_axis)
            vms_gcs = -charge * axis(1, :)**2 + 2.D0 * dipole(slab_axis) * axis(1, :)
            vms_gcs = fact * vms_gcs(:) + const
            !
            !
            !----------------------------------------------------------------------------
            ! Find the edge of the slab
            !
            max_axis = 0.0
            !
            DO i = 1, nnr
                IF (axis(1, i) > max_axis) max_axis = axis(1, i)
            END DO
            !
            v_edge = 0.D0
            icount = 0
            !
            DO i = 1, nnr
                !
                IF (max_axis - axis(1, i) <= 0.1) THEN
                    icount = icount + 1
                    v_edge = v_edge + v%of_r(i) + vms_gcs(i)
                END IF
                !
            END DO
            !
            CALL env_mp_sum(icount, comm)
            !
            CALL env_mp_sum(v_edge, comm)
            !
            v_edge = v_edge / DBLE(icount)
            vms_gcs = vms_gcs - v_edge
            !
            !----------------------------------------------------------------------------
            ! Save in periodic changes to potential
            !
            v%of_r = v%of_r + vms_gcs
            !
            !----------------------------------------------------------------------------
            ! Compute physical properties - MS system
            !
            naxis = v%cell%nr(3)
            !
            avg_window = INT(semiconductor%sc_spread / 2.0 / v%cell%at(3, 3) * naxis)
            !
            IF (io%lnode .AND. io%verbosity >= 1) WRITE (io%debug_unit, 1000)
            !
            IF (slab_charge == 0.D0) THEN
                !
                !------------------------------------------------------------------------
                ! Save planar average of flatband potential
                !
                ez_ms = 0.D0
                !
                IF (.NOT. ALLOCATED(avgd_pot)) ALLOCATE (avgd_pot(naxis))
                !
                IF (.NOT. ALLOCATED(flatband_pot)) ALLOCATE (flatband_pot(naxis))
                !
                CALL v%cell%planar_average(nnr, naxis, 3, 0, .FALSE., v%of_r, &
                                           flatband_pot)
                !
                CALL v%cell%running_average(naxis, avg_window, flatband_pot, avgd_pot)
                !
                semiconductor%flatband_pot_planar_avg = avgd_pot
            ELSE
                !
                !------------------------------------------------------------------------
                ! Generate subtracted potential and corresponding ez_ms
                !
                IF (.NOT. ALLOCATED(current_pot)) ALLOCATE (current_pot(naxis))
                !
                IF (.NOT. ALLOCATED(subtracted_pot)) ALLOCATE (subtracted_pot(naxis))
                !
                IF (.NOT. ALLOCATED(avgd_pot)) ALLOCATE (avgd_pot(naxis))
                !
                CALL v%cell%planar_average(nnr, naxis, 3, 0, .FALSE., v%of_r, current_pot)
                !
                CALL v%cell%running_average(naxis, avg_window, current_pot, avgd_pot)
                !
                current_pot = avgd_pot
                subtracted_pot = current_pot - semiconductor%flatband_pot_planar_avg
                z_cut = this%origin(3) - xstern_ms
                !
                IF (io%lnode .AND. io%verbosity > 1) WRITE (io%debug_unit, 1001) z_cut
                !
                z_cut_idx = INT(z_cut / v%cell%at(3, 3) * naxis)
                !
                v_cut = 0.0
                ez_ms = 0.0
                icount = 0
                !
                !------------------------------------------------------------------------
                ! Averaging v_cut and ez_ms for stability
                !
                DO i = 1, naxis
                    !
                    IF (ABS(z_cut_idx - i) < INT(avg_window / 2)) THEN
                        v_cut = v_cut + subtracted_pot(z_cut_idx)
                        !
                        ez_ms = ez_ms + &
                                (subtracted_pot(z_cut_idx + 1) - &
                                 subtracted_pot(z_cut_idx)) / (v%cell%at(3, 3) / naxis)
                        !
                        icount = icount + 1
                    END IF
                    !
                END DO
                !
                CALL env_mp_sum(icount, comm)
                !
                CALL env_mp_sum(v_cut, comm)
                !
                CALL env_mp_sum(ez_ms, comm)
                !
                v_cut = v_cut / DBLE(icount)
                ez_ms = ez_ms / DBLE(icount)
                !
                !------------------------------------------------------------------------
                ! Calculate what charge ez_ms corresponds to with Gauss law
                !
                semiconductor%ss_v_cut = v_cut
                semiconductor%ss_chg = ez_ms * area / tpi / e2
                !
                IF (io%lnode .AND. io%verbosity > 1) THEN
                    WRITE (io%debug_unit, 1002) v_cut
                    WRITE (io%debug_unit, 1003) ez_ms
                END IF
                !
            END IF
            !
            fact = 1.D0 / carrier_density * permittivity_ms
            arg = fact * (ez_ms**2.D0)
            !
            IF (ez_ms < 0) THEN
                vms = arg
            ELSE
                vms = -arg
            END IF
            !
            IF (io%lnode .AND. io%verbosity >= 1) WRITE (io%debug_unit, 1004) vms
            !
            depletion_length = ABS(2.D0 * fact * ez_ms)
            !
            !----------------------------------------------------------------------------
            ! Set bulk fermi
            !
            semiconductor%bulk_sc_fermi = vms + semiconductor%flatband_fermi + v_cut
            !
            !----------------------------------------------------------------------------
            ! Write potentials
            !
            IF (io%lnode .AND. io%verbosity > 1) THEN
                !
                OPEN (93, file='subtracted_pot.dat', status='replace')
                OPEN (94, file='current_pot.dat', status='replace')
                OPEN (95, file='flatband_pot.dat', status='replace')
                !
                DO i = 1, naxis
                    z_val = i * v%cell%at(3, 3) / naxis
                    WRITE (95, *) z_val, semiconductor%flatband_pot_planar_avg(i)
                    !
                    IF (slab_charge /= 0.D0) THEN
                        WRITE (93, *) z_val, subtracted_pot(i)
                        WRITE (94, *) z_val, current_pot(i)
                    END IF
                    !
                END DO
                !
                CLOSE (93)
                CLOSE (94)
                CLOSE (95)
            END IF
            !
            !----------------------------------------------------------------------------
            !
            CALL local%destroy()
            !
        END ASSOCIATE
        !
        CALL env_stop_clock(routine)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), " MS-GCS ", 68('%'),/)
        !
1001    FORMAT(1X, "MS/DFT cutoff (bohr)       = ", F14.6,/)
        !
1002    FORMAT(1X, "Potential at cutoff        = ", G18.10,/)
        !
1003    FORMAT(1X, "Electric field at cutoff   = ", G18.10,/)
        !
1004    FORMAT(1X, "MS potential               = ", F14.6)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_1da_vms_gcs
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_1da_grad_vms_gcs(this, electrolyte, semiconductor, charges, grad_v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_1da), INTENT(IN) :: this
        TYPE(environ_electrolyte_base), INTENT(IN) :: electrolyte
        TYPE(environ_semiconductor_base), INTENT(IN) :: semiconductor
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_gradient), INTENT(INOUT) :: grad_v
        !
        REAL(DP), POINTER :: grad_vms_gcs(:, :)
        !
        !
        REAL(DP) :: kbt, invkbt
        REAL(DP) :: fact
        REAL(DP) :: charge, dipole(0:3), quadrupole(3)
        !
        TYPE(environ_gradient), TARGET :: glocal
        !
        CHARACTER(LEN=80) :: routine = 'calc_1da_grad_vms_gcs'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(routine)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(grad_v%cell, charges%cell)) &
            CALL io%error(routine, "Mismatch in domains of potential and charges", 1)
        !
        IF (.NOT. ASSOCIATED(grad_v%cell, this%cell)) &
            CALL io%error(routine, "Mismatch in domains of potential and solver", 1)
        !
        IF (electrolyte%ntyp /= 2) &
            CALL io%error(routine, &
                          "Unexpected number of counterionic species, different from two", 1)
        !
        IF (this%dim /= 2) &
            CALL io%error(routine, &
                          "Option not yet implemented: 1D Poisson-Boltzmann solver only for 2D systems", 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (slab_axis => this%axis, &
                   omega => grad_v%cell%omega, &
                   axis => this%x)
            !
            !----------------------------------------------------------------------------
            !
            kbt = semiconductor%temperature * k_boltzmann_ry
            invkbt = 1.D0 / kbt
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL glocal%init(grad_v%cell)
            !
            grad_vms_gcs => glocal%of_r
            !
            !----------------------------------------------------------------------------
            !
            CALL charges%multipoles(this%origin, charge, dipole, quadrupole)
            ! compute multipoles of the system with respect to the chosen origin
            !
            !----------------------------------------------------------------------------
            ! compute the gradient of parabolic correction
            !
            fact = e2 * fpi / omega
            grad_vms_gcs(slab_axis, :) = dipole(slab_axis) - charge * axis(1, :)
            grad_vms_gcs = grad_vms_gcs * fact
            grad_v%of_r = grad_v%of_r + grad_vms_gcs
            !
            CALL glocal%destroy()
            !
        END ASSOCIATE
        !
        CALL env_stop_clock(routine)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_1da_grad_vms_gcs
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_core_1da
!----------------------------------------------------------------------------------------
