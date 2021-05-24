!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE utils_cell
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP, tpi
    !
    USE env_stick_base, ONLY: env_sticks_map, env_sticks_map_deallocate
    !
    USE env_fft_types, ONLY: env_fft_type_descriptor, env_fft_type_init, &
                             env_fft_type_deallocate
    !
    USE types_cell, ONLY: environ_cell
    !
    USE tools_cell, ONLY: volume, recips
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: init_environ_cell, update_environ_cell, destroy_environ_cell
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_cell(gcutm, comm, alat, at, cell, nr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm
        REAL(DP), INTENT(IN) :: gcutm, alat, at(3, 3)
        INTEGER, OPTIONAL, INTENT(IN) :: nr(3)
        !
        TYPE(environ_cell), INTENT(INOUT) :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_cell'
        !
        !--------------------------------------------------------------------------------
        !
        IF (alat < 1.D-8) CALL env_errore(sub_name, 'Wrong alat', 1)
        !
        cell%alat = alat
        !
        !--------------------------------------------------------------------------------
        ! Calculate units of reciprocal space
        !
        cell%tpiba = tpi / alat
        cell%tpiba2 = cell%tpiba**2
        !
        !--------------------------------------------------------------------------------
        ! Create fft descriptor for system cell
        !
        IF (PRESENT(nr)) THEN
            cell%dfft%nr1 = nr(1)
            cell%dfft%nr2 = nr(2)
            cell%dfft%nr3 = nr(3)
        END IF
        !
        CALL init_dfft(gcutm, comm, at, cell%dfft)
        !
        cell%in1 = 1.D0 / DBLE(cell%dfft%nr1)
        cell%in2 = 1.D0 / DBLE(cell%dfft%nr2)
        cell%in3 = 1.D0 / DBLE(cell%dfft%nr3)
        !
        !--------------------------------------------------------------------------------
        ! Real space grid, local dimensions (processor-specific)
        !
        cell%nnr = cell%dfft%nnr
#if defined(__MPI)
        cell%j0 = cell%dfft%my_i0r2p
        cell%k0 = cell%dfft%my_i0r3p
        !
        cell%ir_end = MIN(cell%nnr, &
                          cell%dfft%nr1x * cell%dfft%my_nr2p * cell%dfft%my_nr3p)
#else
        cell%j0 = 0
        cell%k0 = 0
        cell%ir_end = cell%nnr
#endif
        !
        cell%ntot = cell%dfft%nr1 * cell%dfft%nr2 * cell%dfft%nr3
        ! total number of physical points
        !
        !--------------------------------------------------------------------------------
        ! Set basic cell properties
        !
        CALL update_environ_cell(at, cell)
        !
        cell%origin = 0.D0
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_cell(at, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: at(3, 3)
        !
        TYPE(environ_cell), INTENT(INOUT) :: cell
        !
        INTEGER :: ic, ix, iy, iz
        REAL(DP) :: dx, dy, dz
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_cell'
        !
        !--------------------------------------------------------------------------------
        !
        cell%at = at
        !
        CALL volume(cell%alat, cell%at(1, 1), cell%at(1, 2), cell%at(1, 3), cell%omega)
        ! calculate cell volume
        !
        cell%cubic = iscubic(cell%at) ! Check if the cell is cubic
        !
        ! calculate reciprocal cell
        CALL recips(cell%at(1, 1), cell%at(1, 2), cell%at(1, 3), &
                    cell%bg(1, 1), cell%bg(1, 2), cell%bg(1, 3))
        !
        cell%domega = cell%omega / cell%ntot ! Set volume element
        !
        !--------------------------------------------------------------------------------
        ! Calculate corners for minimum image convention
        !
        ic = 0
        !
        DO ix = 0, 1
            dx = DBLE(-ix)
            !
            DO iy = 0, 1
                dy = DBLE(-iy)
                !
                DO iz = 0, 1
                    dz = DBLE(-iz)
                    ic = ic + 1
                    cell%corners(1, ic) = dx * at(1, 1) + dy * at(1, 2) + dz * at(1, 3)
                    cell%corners(2, ic) = dx * at(2, 1) + dy * at(2, 2) + dz * at(2, 3)
                    cell%corners(3, ic) = dx * at(3, 1) + dy * at(3, 2) + dz * at(3, 3)
                END DO
                !
            END DO
            !
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_cell(lflag, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(environ_cell), INTENT(INOUT) :: cell
        !
        !--------------------------------------------------------------------------------
        !
        CALL destroy_dfft(lflag, cell%dfft)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_dfft(gcutm, comm, at, dfft)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm
        REAL(DP), INTENT(IN) :: gcutm, at(3, 3)
        !
        TYPE(env_fft_type_descriptor), INTENT(INOUT) :: dfft
        !
        TYPE(env_sticks_map) :: smap
        REAL(DP) :: bg(3, 3)
        !
        !--------------------------------------------------------------------------------
        ! Calculate the reciprocal lattice vectors
        !
        CALL recips(at(1, 1), at(1, 2), at(1, 3), bg(1, 1), bg(1, 2), bg(1, 3))
        !
        CALL env_fft_type_init(dfft, smap, .TRUE., .TRUE., comm, at, bg, gcutm, &
                               nyfft=1, nmany=1)
        !
        dfft%rho_clock_label = 'fft'
        !
        CALL env_sticks_map_deallocate(smap)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_dfft
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_dfft(lflag, dfft)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(env_fft_type_descriptor), INTENT(INOUT) :: dfft
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_fft_type_deallocate(dfft)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_dfft
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION iscubic(at)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), PARAMETER :: tol = 1.D-8
        REAL(DP) :: at(3, 3)
        LOGICAL :: iscubic
        !
        INTEGER :: ipol, jpol
        REAL(DP) :: tmp
        !
        !--------------------------------------------------------------------------------
        ! If at(3, 3) is a cubic cell, at(1, 1) = at(2, 2) = at(3, 3) and
        ! the other elements are equal to 0.D0
        !
        iscubic = .FALSE.
        !
        tmp = 0.D0
        !
        DO ipol = 1, 3
            !
            DO jpol = 1, 3
                !
                IF (ipol == jpol) THEN
                    tmp = tmp + ABS(at(ipol, ipol) - at(1, 1))
                ELSE
                    tmp = tmp + ABS(at(ipol, jpol))
                END IF
                !
            END DO
            !
        END DO
        !
        iscubic = tmp < tol
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION iscubic
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_cell
!----------------------------------------------------------------------------------------
