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
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE init_core
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    USE types_cell, ONLY: environ_cell
    !
    USE base_core
    !
    USE utils_fd
    USE utils_fft
    USE utils_oned_analytic
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: set_core_base, core_initbase, core_initcell, core_initions, core_clean
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_core_base(ifdtype, nfdpoint, use_internal_pbc_corr, dim, axis)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: use_internal_pbc_corr
        INTEGER, INTENT(IN) :: ifdtype, nfdpoint, dim, axis
        !
        !--------------------------------------------------------------------------------
        ! Set up active numerical cores
        !
        IF (lfd) CALL init_fd_core_first(ifdtype, nfdpoint, fd)
        !
        IF (loned_analytic) CALL init_oned_analytic_core_first(dim, axis, oned_analytic)
        !
        IF (lfft_system) CALL init_fft_core_first(system_fft, use_internal_pbc_corr)
        !
        IF (lfft_environment) &
            CALL init_fft_core_first(environment_fft, use_internal_pbc_corr)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_core_base
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE core_initbase(gcutm, environment_cell, system_cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: gcutm
        TYPE(environ_cell), INTENT(IN) :: system_cell, environment_cell
        !
        !--------------------------------------------------------------------------------
        !
        IF (loned_analytic) &
            CALL init_oned_analytic_core_second(environment_cell, oned_analytic)
        !
        IF (lfd) CALL init_fd_core_second(environment_cell, fd)
        !
        IF (lfft_environment) &
            CALL init_fft_core_second(gcutm, environment_cell, environment_fft)
        !
        IF (lfft_system) CALL init_fft_core_second(gcutm, system_cell, system_fft)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE core_initbase
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE core_initcell(system_cell, environment_cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: system_cell, environment_cell
        !
        !--------------------------------------------------------------------------------
        !
        IF (loned_analytic) &
            CALL update_oned_analytic_core_cell(environment_cell, oned_analytic)
        !
        IF (lfft_environment) &
            CALL update_fft_core_cell(environment_cell, environment_fft)
        !
        IF (lfft_system) CALL update_fft_core_cell(system_cell, system_fft)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE core_initcell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE core_initions(pos)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: pos(3)
        !
        !--------------------------------------------------------------------------------
        !
        IF (loned_analytic) CALL update_oned_analytic_core_origin(pos, oned_analytic)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE core_initions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE core_clean(lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        !--------------------------------------------------------------------------------
        !
        IF (lfft_environment) CALL destroy_fft_core(lflag, environment_fft)
        !
        IF (lfft_system) CALL destroy_fft_core(lflag, system_fft)
        !
        IF (lfd) CALL destroy_fd_core(lflag, fd)
        !
        IF (loned_analytic) CALL destroy_oned_analytic_core(lflag, oned_analytic)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE core_clean
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE init_core
!----------------------------------------------------------------------------------------
