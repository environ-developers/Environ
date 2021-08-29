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
MODULE class_core_fd
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    !
    USE class_core_numerical
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
    TYPE, EXTENDS(numerical_core), PUBLIC :: core_fd
        !--------------------------------------------------------------------------------
        !
        INTEGER :: ifdtype
        INTEGER :: nfdpoint
        INTEGER, ALLOCATABLE :: icfd(:)
        INTEGER :: ncfd
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_core_fd
        PROCEDURE :: init => init_core_fd
        PROCEDURE :: destroy => destroy_core_fd
        !
        PROCEDURE :: set_coefficients
        !
        !--------------------------------------------------------------------------------
    END TYPE core_fd
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
    SUBROUTINE create_core_fd(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fd), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_core_fd'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%cell)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        IF (ALLOCATED(this%icfd)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (this%cell)
        !
        this%core_type = 'fd'
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_core_fd
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_core_fd(this, ifdtype, nfdpoint, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ifdtype, nfdpoint
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        CLASS(core_fd), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        ! Set finite differences tools
        !
        CALL this%create()
        !
        this%ifdtype = ifdtype
        this%nfdpoint = nfdpoint
        ALLOCATE (this%icfd(-nfdpoint:nfdpoint))
        !
        CALL this%set_coefficients()
        !
        this%cell => cell
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_core_fd
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_core_fd(this, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(core_fd), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_core_fd'
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (this%cell)
        !
        IF (lflag) THEN
            !
            IF (.NOT. ALLOCATED(this%icfd)) &
                CALL env_errore(sub_name, 'Trying to deallocate an empty object', 1)
            !
            DEALLOCATE (this%icfd)
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_core_fd
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               PRIVATE HELPER METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_coefficients(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fd), TARGET, INTENT(INOUT) :: this
        !
        INTEGER, POINTER :: ifdtype, nfdpoint, ncfd
        INTEGER, POINTER :: icfd(:)
        !
        INTEGER :: in
        !
        !--------------------------------------------------------------------------------
        !
        ifdtype => this%ifdtype
        nfdpoint => this%nfdpoint
        ncfd => this%ncfd
        icfd => this%icfd
        !
        ncfd = 0
        icfd = 0
        !
        SELECT CASE (ifdtype)
            !
        CASE (1) ! (2N+1)-point Central Differences
            !
            IF (nfdpoint == 1) THEN
                ncfd = 2
                icfd(1) = 1
            ELSE IF (nfdpoint == 2) THEN
                ncfd = 12
                icfd(2) = -1
                icfd(1) = 8
            ELSE IF (nfdpoint == 3) THEN
                ncfd = 60
                icfd(3) = 1
                icfd(2) = -9
                icfd(1) = 45
            ELSE IF (nfdpoint == 4) THEN
                ncfd = 840
                icfd(4) = -3
                icfd(3) = 32
                icfd(2) = -168
                icfd(1) = 672
            ELSE
                !
                !------------------------------------------------------------------------
                ! DEBUGGING
                !
                WRITE (*, *) 'ERROR: wrong number of points', nfdpoint, &
                    ' for finite difference type ', ifdtype
                !
                STOP
                !
            END IF
            !
        CASE (2) ! low-Noise Lanczos Differentiators ( M = 2 )
            !
            IF (nfdpoint >= 2) THEN
                ncfd = (nfdpoint) * (nfdpoint + 1) * (2 * nfdpoint + 1) / 3
                !
                DO in = 1, nfdpoint
                    icfd(in) = in
                END DO
                !
            ELSE
                !
                !------------------------------------------------------------------------
                ! DEBUGGING
                !
                WRITE (*, *) 'ERROR: wrong number of points', nfdpoint, &
                    ' for finite difference type ', ifdtype
                !
                STOP
                !
            END IF
            !
        CASE (3) ! super Lanczos Low-Noise Differentiators ( M = 4 )
            !
            IF (nfdpoint == 3) THEN
                ncfd = 252
                icfd(3) = -22
                icfd(2) = 67
                icfd(1) = 58
            ELSE IF (nfdpoint == 4) THEN
                ncfd = 1188
                icfd(4) = -86
                icfd(3) = 142
                icfd(2) = 193
                icfd(1) = 126
            ELSE IF (nfdpoint == 5) THEN
                ncfd = 5148
                icfd(5) = -300
                icfd(4) = 294
                icfd(3) = 532
                icfd(2) = 503
                icfd(1) = 296
            ELSE
                !
                !------------------------------------------------------------------------
                ! DEBUGGING
                !
                WRITE (*, *) 'ERROR: wrong number of points', nfdpoint, &
                    ' for finite difference type ', ifdtype
                !
                STOP
                !
            END IF
            !
        CASE (4) ! smooth Noise-Robust Differentiators  ( n = 2 )
            !
            IF (nfdpoint == 2) THEN
                ncfd = 8
                icfd(2) = 1
                icfd(1) = 2
            ELSE IF (nfdpoint == 3) THEN
                ncfd = 32
                icfd(3) = 1
                icfd(2) = 4
                icfd(1) = 5
            ELSE IF (nfdpoint == 4) THEN
                ncfd = 128
                icfd(4) = 1
                icfd(3) = 6
                icfd(2) = 14
                icfd(1) = 14
            ELSE IF (nfdpoint == 5) THEN
                ncfd = 512
                icfd(5) = 1
                icfd(4) = 8
                icfd(3) = 27
                icfd(2) = 48
                icfd(1) = 42
            ELSE
                !
                !------------------------------------------------------------------------
                ! DEBUGGING
                !
                WRITE (*, *) 'ERROR: wrong number of points', nfdpoint, &
                    ' for finite difference type ', ifdtype
                !
                STOP
                !
            END IF
            !
        CASE (5) ! smooth Noise-Robust Differentiators  ( n = 4 )
            !
            IF (nfdpoint == 3) THEN
                ncfd = 96
                icfd(3) = -5
                icfd(2) = 12
                icfd(1) = 39
            ELSE IF (nfdpoint == 4) THEN
                ncfd = 96
                icfd(4) = -2
                icfd(3) = -1
                icfd(2) = 16
                icfd(1) = 27
            ELSE IF (nfdpoint == 5) THEN
                ncfd = 1536
                icfd(5) = -11
                icfd(4) = -32
                icfd(3) = 39
                icfd(2) = 256
                icfd(1) = 322
            ELSE
                !
                !------------------------------------------------------------------------
                ! DEBUGGING
                !
                WRITE (*, *) 'ERROR: wrong number of points', nfdpoint, &
                    ' for finite difference type ', ifdtype
                !
                STOP
                !
            END IF
            !
        CASE DEFAULT
            !
            !----------------------------------------------------------------------------
            ! DEBUGGING
            !
            WRITE (*, *) 'ERROR: finite difference type unknown, ifdtype=', ifdtype
            !
            STOP
            !
        END SELECT
        !
        DO in = 1, nfdpoint
            icfd(-in) = -icfd(in)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_coefficients
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_core_fd
!----------------------------------------------------------------------------------------
