!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
! Copyright (C) 2011 Quantum ESPRESSO group
!
!----------------------------------------------------------------------------------------
!
! This file is part of Environ version 2.0
!
! Environ 2.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
!
! Environ 2.0 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more detail, either the file
! `License' in the root directory of the present distribution, or
! online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors:
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE env_fftw_interfaces
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
    INTERFACE
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE env_create_plan_1d(plan, nz, i) &
            BIND(C, name="env_create_plan_1d")
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_PTR, C_INT
            TYPE(C_PTR) :: plan
            INTEGER(KIND=C_INT) :: nz, i
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE env_create_plan_1d
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE env_destroy_plan_1d(plan) &
            BIND(C, name="env_destroy_plan_1d")
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_PTR
            TYPE(C_PTR) :: plan
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE env_destroy_plan_1d
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE env_create_plan_2d(plan, nx, ny, i) &
            BIND(C, name="env_create_plan_2d")
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_PTR, C_INT
            TYPE(C_PTR) :: plan
            INTEGER(KIND=C_INT) :: nx, ny, i
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE env_create_plan_2d
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE env_destroy_plan_2d(plan) &
            BIND(C, name="env_destroy_plan_2d")
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_PTR
            TYPE(C_PTR) :: plan
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE env_destroy_plan_2d
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE env_create_plan_3d(plan, nx, ny, nz, i) &
            BIND(C, name="env_create_plan_3d")
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_PTR, C_INT
            TYPE(C_PTR) :: plan
            INTEGER(KIND=C_INT) :: nx, ny, nz, i
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE env_create_plan_3d
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE env_destroy_plan_3d(plan) &
            BIND(C, name="env_destroy_plan_3d")
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_PTR
            TYPE(C_PTR) :: plan
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE env_destroy_plan_3d
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE env_fft_x_stick(plan, xy, nx, ny, nz, ldx, ldy) &
            BIND(C, name="env_fft_x_stick")
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE
            TYPE(C_PTR) :: plan
            INTEGER(KIND=C_INT) :: nx, ny, nz, ldx, ldy
            COMPLEX(KIND=C_DOUBLE) :: xy
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE env_fft_x_stick
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE env_fft_x_stick_single(plan, xy, nx, ny, nz, ldx, ldy) &
            BIND(C, name="env_fft_x_stick_single")
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE
            TYPE(C_PTR) :: plan
            INTEGER(KIND=C_INT) :: nx, ny, nz, ldx, ldy
            COMPLEX(KIND=C_DOUBLE) :: xy
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE env_fft_x_stick_single
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE env_fft_y_stick(plan, xy, ny, ldx) &
            BIND(C, name="env_fft_y_stick")
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE
            TYPE(C_PTR) :: plan
            INTEGER(KIND=C_INT) :: ny, ldx
            COMPLEX(KIND=C_DOUBLE) :: xy
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE env_fft_y_stick
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE env_fft_z_stick(plan, z, ldz, nzl) &
            BIND(C, name="env_fft_z_stick")
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE
            TYPE(C_PTR) :: plan
            INTEGER(KIND=C_INT) :: ldz, nzl
            COMPLEX(KIND=C_DOUBLE) :: z
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE env_fft_z_stick
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE env_fft_z_stick_single(plan, z, ldz) &
            BIND(C, name="env_fft_z_stick_single")
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE
            TYPE(C_PTR) :: plan
            INTEGER(KIND=C_INT) :: ldz
            COMPLEX(KIND=C_DOUBLE) :: z
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE env_fft_z_stick_single
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE env_fftw_inplace_drv_1d(plan, m, z, inc1, inc2) &
            BIND(C, name="env_fftw_inplace_drv_1d")
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE
            TYPE(C_PTR) :: plan
            INTEGER(KIND=C_INT) :: m, inc1, inc2
            COMPLEX(KIND=C_DOUBLE) :: z
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE env_fftw_inplace_drv_1d
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE env_fftw_inplace_drv_3d(plan, m, z, inc1, inc2) &
            BIND(C, name="env_fftw_inplace_drv_3d")
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE
            TYPE(C_PTR) :: plan
            INTEGER(KIND=C_INT) :: m, inc1, inc2
            COMPLEX(KIND=C_DOUBLE) :: z
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE env_fftw_inplace_drv_3d
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    END INTERFACE
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_fftw_interfaces
!----------------------------------------------------------------------------------------
