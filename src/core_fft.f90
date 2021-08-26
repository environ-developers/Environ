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
MODULE class_core_fft
    !------------------------------------------------------------------------------------
    !
    USE env_mp, ONLY: env_mp_sum
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    !
    USE class_core_numerical
    !
    USE generate_gvectors, ONLY: env_ggen
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
    TYPE, EXTENDS(numerical_core), PUBLIC :: core_fft
        !--------------------------------------------------------------------------------
        !
        INTEGER :: ngm ! local  number of G vectors (on this processor)
        ! with gamma tricks, only vectors in G>
        !
        REAL(DP) :: gcutm ! ecutrho/(2 pi/a)^2, cut-off for |G|^2
        !
        INTEGER :: gstart ! index of the first G vector whose module is > 0
        ! needed in parallel execution:
        ! gstart=2 for the proc that holds G=0
        ! gstart=1 for all others
        !
        REAL(DP), ALLOCATABLE :: gg(:)
        ! G^2 in increasing order (in units of (2pi/a)^2)
        !
        REAL(DP), ALLOCATABLE :: g(:, :)
        ! G-vectors cartesian components (in units 2pi/a)
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_core_fft
        PROCEDURE :: init_first => init_core_fft_first
        PROCEDURE :: init_second => init_core_fft_second
        PROCEDURE :: update_cell => update_core_fft_cell
        PROCEDURE :: destroy => destroy_core_fft
        !
        PROCEDURE, PRIVATE :: init_gvect => env_gvect_init
        PROCEDURE, PRIVATE :: deallocate_gvect => env_deallocate_gvect
        !
        !--------------------------------------------------------------------------------
    END TYPE core_fft
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
    SUBROUTINE create_core_fft(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_core_fft'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%cell)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        IF (ALLOCATED(this%g)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        IF (ALLOCATED(this%gg)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (this%cell)
        !
        this%core_type = 'fft'
        !
        this%ngm = 0
        this%gcutm = 0.0_DP
        this%gstart = 2
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_core_fft
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_core_fft_first(this, use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN), OPTIONAL :: use_internal_pbc_corr
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_core_fft_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_core_fft_second(this, gcutm, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: gcutm
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        INTEGER :: ngm_g ! global number of G vectors (summed on all procs)
        ! in serial execution, ngm_g = ngm
        !
        !--------------------------------------------------------------------------------
        !
        this%gcutm = gcutm
        this%cell => cell
        !
        this%ngm = cell%dfft%ngm
        !
        !--------------------------------------------------------------------------------
        ! #TODO The following routines are in generate_gvectors
        ! and may need to be simplified
        !
        CALL this%init_gvect(ngm_g, cell%dfft%comm)
        !
        CALL env_ggen(this%cell%dfft, cell%dfft%comm, cell%at, cell%bg, this%gcutm, &
                      ngm_g, this%ngm, this%g, this%gg, this%gstart, .TRUE.)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_core_fft_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_core_fft_cell(this, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%cell => cell
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_core_fft_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_core_fft(this, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_core_fft'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%cell)) &
                CALL env_errore(sub_name, 'Trying to destroy an empty object', 1)
        !
        NULLIFY (this%cell)
        !
        CALL this%deallocate_gvect()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_core_fft
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               PRIVATE HELPER METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Set local and global dimensions, allocate arrays
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_gvect_init(this, ngm_g, comm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        INTEGER, INTENT(INOUT) :: ngm_g
        !
        INTEGER :: ngm
        !
        INTEGER, INTENT(IN) :: comm
        ! communicator of the group on which g-vecs are distributed
        !
        !--------------------------------------------------------------------------------
        ! Calculate sum over all processors
        !
        ngm = this%ngm
        ngm_g = ngm
        !
        CALL env_mp_sum(ngm_g, comm)
        !
        !--------------------------------------------------------------------------------
        ! Allocate arrays - only those that are always kept until the end
        !
        ALLOCATE (this%gg(ngm))
        ALLOCATE (this%g(3, ngm))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_gvect_init
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_deallocate_gvect(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(this%gg)) DEALLOCATE (this%gg)
        !
        IF (ALLOCATED(this%g)) DEALLOCATE (this%g)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_deallocate_gvect
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_core_fft
!----------------------------------------------------------------------------------------
