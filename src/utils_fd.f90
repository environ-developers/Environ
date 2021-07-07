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
MODULE utils_fd
    !------------------------------------------------------------------------------------
    !
    USE types_core, ONLY: fd_core
    USE types_cell, ONLY: environ_cell
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: init_fd_core_first, init_fd_core_second, destroy_fd_core
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! #TODO unused
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_fd_core(fd)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(fd_core), INTENT(INOUT) :: fd
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (fd%cell) ! create empty finite difference core
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_fd_core
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_fd_core_first(ifdtype, nfdpoint, fd)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ifdtype, nfdpoint
        !
        TYPE(fd_core), INTENT(INOUT) :: fd
        !
        !--------------------------------------------------------------------------------
        ! Set finite differences tools
        !
        fd%ifdtype = ifdtype
        fd%nfdpoint = nfdpoint
        ALLOCATE (fd%icfd(-nfdpoint:nfdpoint))
        !
        CALL set_fd_coefficients(fd)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_fd_core_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_fd_core_second(cell, fd)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        TYPE(fd_core), INTENT(INOUT) :: fd
        !
        !--------------------------------------------------------------------------------
        !
        fd%cell => cell
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_fd_core_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_fd_core(lflag, fd)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(fd_core), INTENT(INOUT) :: fd
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_fd_core'
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (fd%cell)
        !
        IF (lflag) THEN
            !
            IF (.NOT. ALLOCATED(fd%icfd)) &
                CALL env_errore(sub_name, &
                                'Trying to deallocate a non-allocated object', 1)
            !
            DEALLOCATE (fd%icfd)
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_fd_core
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_fd_coefficients(fd)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(fd_core), TARGET, INTENT(INOUT) :: fd
        !
        INTEGER, POINTER :: ifdtype, nfdpoint, ncfd
        INTEGER, POINTER :: icfd(:)
        !
        INTEGER :: in
        !
        !--------------------------------------------------------------------------------
        !
        ifdtype => fd%ifdtype
        nfdpoint => fd%nfdpoint
        ncfd => fd%ncfd
        icfd => fd%icfd
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
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_fd_coefficients
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_fd
!----------------------------------------------------------------------------------------
