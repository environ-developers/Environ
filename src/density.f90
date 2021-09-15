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
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_density
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_mp, ONLY: env_mp_sum
    !
    USE env_base_scatter, ONLY: env_gather_grid
    USE env_types_fft, ONLY: env_fft_type_descriptor
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
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
    TYPE, PUBLIC :: environ_density
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lupdate = .FALSE. ! optionally have an associated logical status
        !
        CHARACTER(LEN=80) :: label = 'density'
        ! optionally have an associated label, used for printout and debugs
        !
        TYPE(environ_cell), POINTER :: cell => NULL()
        ! each quantity in real-space is associated with its definition domain
        !
        REAL(DP), ALLOCATABLE :: of_r(:)
        ! the quantity in real-space, local to each processor
        !
        !--------------------------------------------------------------------------------
        ! Multipole moments of the quantity
        !
        REAL(DP) :: charge = 0.D0
        REAL(DP) :: dipole(3) = 0.D0
        REAL(DP) :: quadrupole(3) = 0.D0
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_density
        PROCEDURE :: init => init_environ_density
        PROCEDURE :: copy => copy_environ_density
        PROCEDURE :: destroy => destroy_environ_density
        !
        PROCEDURE :: multipoles => multipoles_environ_density
        PROCEDURE :: integrate => integrate_environ_density
        PROCEDURE :: euclidean_norm => euclidean_norm_environ_density
        PROCEDURE :: quadratic_mean => quadratic_mean_environ_density
        PROCEDURE :: scalar_product => scalar_product_environ_density
        !
        PROCEDURE :: dipole_of_origin, quadrupole_of_origin
        !
        PROCEDURE :: printout => print_environ_density
        PROCEDURE :: write_cube => write_cube_density
        !
        PROCEDURE :: write_cube_no_ions
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_density
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
    SUBROUTINE create_environ_density(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_density), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%cell)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%of_r)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_density(this, cell, label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        CHARACTER(LEN=80), INTENT(IN), OPTIONAL :: label
        !
        CLASS(environ_density), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        IF (PRESENT(label)) this%label = label
        !
        this%cell => cell
        !
        ALLOCATE (this%of_r(this%cell%nnr))
        this%of_r = 0.D0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE copy_environ_density(this, copy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_density), INTENT(IN) :: this
        !
        TYPE(environ_density), INTENT(OUT) :: copy
        !
        CHARACTER(LEN=80) :: sub_name = 'copy_environ_density'
        !
        INTEGER :: n
        !
        !--------------------------------------------------------------------------------
        !
        copy%cell => this%cell
        !
        copy%lupdate = this%lupdate
        copy%label = this%label
        copy%charge = this%charge
        copy%dipole = this%dipole
        copy%quadrupole = this%quadrupole
        !
        IF (ALLOCATED(this%of_r)) THEN
            n = SIZE(this%of_r)
            !
            IF (ALLOCATED(copy%of_r)) DEALLOCATE (copy%of_r)
            !
            ALLOCATE (copy%of_r(n))
            copy%of_r = this%of_r
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE copy_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_density(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_density), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%cell)) CALL io%destroy_error(sub_name)
        !
        IF (.NOT. ALLOCATED(this%of_r)) CALL io%destroy_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (this%cell)
        !
        DEALLOCATE (this%of_r)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_density
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  GENERAL METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE multipoles_environ_density(this, origin, monopole, dipole, quadrupole)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_density), TARGET, INTENT(IN) :: this
        REAL(DP), INTENT(IN) :: origin(3)
        !
        REAL(DP), INTENT(OUT) :: monopole
        REAL(DP), DIMENSION(3), INTENT(OUT) :: dipole, quadrupole
        !
        TYPE(environ_cell), POINTER :: cell
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: r(3), rhoir, r2
        INTEGER :: dim, axis
        !
        !--------------------------------------------------------------------------------
        !
        cell => this%cell
        !
        monopole = 0.D0
        dipole = 0.D0
        quadrupole = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL cell%get_min_distance(ir, 0, 3, origin, r, r2, physical)
            ! compute minimum distance using minimum image convention
            !
            IF (.NOT. physical) CYCLE
            !
            rhoir = this%of_r(ir)
            !
            !----------------------------------------------------------------------------
            ! Multipoles
            !
            monopole = monopole + rhoir
            dipole = dipole - rhoir * r
            quadrupole = quadrupole + rhoir * r**2
            !
        END DO
        !
        CALL env_mp_sum(monopole, cell%dfft%comm)
        !
        CALL env_mp_sum(dipole, cell%dfft%comm)
        !
        CALL env_mp_sum(quadrupole, cell%dfft%comm)
        !
        monopole = monopole * cell%domega
        dipole = dipole * cell%domega
        quadrupole = quadrupole * cell%domega
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE multipoles_environ_density
    !------------------------------------------------------------------------------------
    !>
    !! #TODO unused
    !!
    !------------------------------------------------------------------------------------
    FUNCTION dipole_of_origin(this, origin) RESULT(dipole)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_density), INTENT(IN) :: this
        !
        REAL(DP), DIMENSION(3) :: origin, dipole
        !
        !--------------------------------------------------------------------------------
        !
        dipole = this%dipole + this%charge * (this%cell%origin - origin)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION dipole_of_origin
    !------------------------------------------------------------------------------------
    !>
    !! #TODO unused
    !!
    !------------------------------------------------------------------------------------
    FUNCTION quadrupole_of_origin(this, origin) RESULT(quadrupole)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_density), INTENT(IN) :: this
        !
        REAL(DP), DIMENSION(3) :: origin, quadrupole
        !
        !--------------------------------------------------------------------------------
        !
        quadrupole = this%quadrupole + &
                     this%charge * (this%cell%origin - origin)**2 + &
                     2.D0 * this%dipole * (this%cell%origin - origin)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION quadrupole_of_origin
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION integrate_environ_density(this) RESULT(integral)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_density), INTENT(IN) :: this
        !
        REAL(DP) :: integral
        !
        !--------------------------------------------------------------------------------
        !
        integral = SUM(this%of_r(1:this%cell%ir_end))
        !
        CALL env_mp_sum(integral, this%cell%dfft%comm)
        !
        integral = integral * this%cell%domega
        !
        !--------------------------------------------------------------------------------
    END FUNCTION integrate_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION euclidean_norm_environ_density(this) RESULT(euclidean_norm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_density), INTENT(IN) :: this
        !
        INTEGER, POINTER :: ir_end
        !
        REAL(DP) :: euclidean_norm
        !
        !--------------------------------------------------------------------------------
        !
        ir_end => this%cell%ir_end
        euclidean_norm = DOT_PRODUCT(this%of_r(1:ir_end), this%of_r(1:ir_end))
        !
        CALL env_mp_sum(euclidean_norm, this%cell%dfft%comm)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION euclidean_norm_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION quadratic_mean_environ_density(this) RESULT(quadratic_mean)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_density), INTENT(IN) :: this
        !
        INTEGER, POINTER :: ir_end
        !
        REAL(DP) :: quadratic_mean
        !
        !--------------------------------------------------------------------------------
        !
        ir_end => this%cell%ir_end
        quadratic_mean = DOT_PRODUCT(this%of_r(1:ir_end), this%of_r(1:ir_end))
        !
        CALL env_mp_sum(quadratic_mean, this%cell%dfft%comm)
        !
        quadratic_mean = SQRT(quadratic_mean / this%cell%ntot)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION quadratic_mean_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION scalar_product_environ_density(this, density2) RESULT(scalar_product)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_density), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: density2
        !
        !
        INTEGER, POINTER :: ir_end
        REAL(DP) :: scalar_product
        !
        CHARACTER(LEN=80) :: fun_name = 'scalar_product_environ_density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%cell, density2%cell)) &
            CALL io%error(fun_name, 'Operation on fields with inconsistent domains', 1)
        !
        ir_end => this%cell%ir_end
        scalar_product = DOT_PRODUCT(this%of_r(1:ir_end), density2%of_r(1:ir_end))
        !
        CALL env_mp_sum(scalar_product, this%cell%dfft%comm)
        !
        scalar_product = scalar_product * this%cell%domega
        !
        !--------------------------------------------------------------------------------
    END FUNCTION scalar_product_environ_density
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Prints the details of the density
    !!
    !! If called by a parent object, prints details in block format
    !!
    !! @param verbose       : (INTEGER) adds verbosity to global verbose
    !! @param debug_verbose : (INTEGER) replaces global verbose for debugging
    !! @param unit          : (INTEGER) output target (default = io%debug_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_density(this, verbose, debug_verbose, unit, lcube)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_density), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: verbose, debug_verbose, unit
        LOGICAL, INTENT(IN), OPTIONAL :: lcube
        !
        INTEGER :: base_verbose, local_verbose, local_unit
        !
        LOGICAL :: print_cube = .TRUE.
        REAL(DP) :: integral
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(debug_verbose)) THEN
            base_verbose = debug_verbose
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = verbose
            ELSE
                local_verbose = debug_verbose
            END IF
            !
        ELSE IF (io%verbosity > 0) THEN
            base_verbose = io%verbosity
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = base_verbose + verbose
            ELSE
                local_verbose = base_verbose
            END IF
            !
        ELSE
            RETURN
        END IF
        !
        IF (PRESENT(unit)) THEN
            local_unit = unit
        ELSE
            local_unit = io%debug_unit
        END IF
        !
        IF (PRESENT(lcube)) print_cube = lcube
        !
        IF (local_verbose >= 1) THEN
            !
            IF (io%lnode) THEN
                !
                IF (local_verbose >= base_verbose) THEN ! header
                    WRITE (local_unit, 1000)
                ELSE
                    !
                    CALL io%block_divider(local_verbose, base_verbose, local_unit)
                    !
                    WRITE (local_unit, 1001)
                END IF
                !
                WRITE (local_unit, 1002) ADJUSTL(this%label)
            END IF
            !
            integral = this%integrate()
            !
            IF (io%lnode) WRITE (local_unit, 1003) integral
            !
            ! #TODO ADD MAXVAL AND MINVAL
            !
            IF (local_verbose >= 3 .AND. print_cube) CALL this%write_cube_no_ions()
            !
            IF (local_verbose < base_verbose) &
                CALL io%block_divider(local_verbose, base_verbose, local_unit)
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), ' DENSITY ', 67('%'))
1001    FORMAT(/, ' DENSITY', /, ' =======')
        !
1002    FORMAT(/, ' density label              = ', A50)
        !
1003    FORMAT(/, ' integral of density        = ', G18.10)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE write_cube_no_ions(this, idx, label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_density), TARGET, INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: idx
        CHARACTER(LEN=100), INTENT(IN), OPTIONAL :: label
        !
        CHARACTER(LEN=100) :: filename, filemod, local_label
        !
        !--------------------------------------------------------------------------------
        ! Set filename
        !
        IF (PRESENT(idx)) THEN
            WRITE (filemod, '(i4.4)') idx
        ELSE
            filemod = ""
        END IF
        !
        IF (PRESENT(label)) THEN
            local_label = label
        ELSE
            local_label = this%label
        END IF
        !
        filename = TRIM(ADJUSTL(local_label))//TRIM(filemod)//'.cube'
        !
        !--------------------------------------------------------------------------------
        ! Write cube
        !
        OPEN (300, file=TRIM(filename), status='unknown')
        !
        CALL this%cell%write_cube(1) ! write cube cell data
        !
        WRITE (300, '(i5,4f12.6)') 1, 0.D0, 0.D0, 0.D0, 0.D0 ! hydrogen ion
        !
        CALL this%write_cube()
        !
        CLOSE (300)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE write_cube_no_ions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE write_cube_density(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_density), TARGET, INTENT(IN) :: this
        !
        INTEGER :: ir, ir1, ir2, ir3
        INTEGER :: count
        INTEGER :: nr1x, nr2x, nr3x
        INTEGER :: nr1, nr2, nr3
        !
        REAL(DP) :: tmp
        REAL(DP), ALLOCATABLE :: flocal(:)
        !
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        CHARACTER(LEN=80) :: sub_name = 'write_cube_density'
        !
        !--------------------------------------------------------------------------------
        !
        dfft => this%cell%dfft
        !
        nr1x = dfft%nr1x
        nr2x = dfft%nr2x
        nr3x = dfft%nr3x
        !
        nr1 = dfft%nr1
        nr2 = dfft%nr2
        nr3 = dfft%nr3
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE (flocal(nr1x * nr2x * nr3x))
#if defined(__MPI)
        flocal = 0.D0
        !
        CALL env_gather_grid(dfft, this%of_r, flocal)
        !
        CALL env_mp_sum(flocal, dfft%comm)
        !
#else
        flocal = this%of_r
#endif
        !
        !--------------------------------------------------------------------------------
        ! Write cube density data
        !
        count = 0
        !
        DO ir1 = 1, nr1
            !
            DO ir2 = 1, nr2
                !
                DO ir3 = 1, nr3
                    count = count + 1
                    ir = ir1 + (ir2 - 1) * nr1 + (ir3 - 1) * nr1 * nr2
                    tmp = DBLE(flocal(ir))
                    !
                    IF (ABS(tmp) < 1.D-99) tmp = 0.D0
                    !
                    IF (MOD(count, 6) == 0) THEN
                        WRITE (300, '(e13.6,1x)') tmp
                    ELSE
                        WRITE (300, '(e13.6,1x)', advance='no') tmp
                    END IF
                    !
                END DO
                !
            END DO
            !
        END DO
        !
        DEALLOCATE (flocal)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE write_cube_density
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_density
!----------------------------------------------------------------------------------------
