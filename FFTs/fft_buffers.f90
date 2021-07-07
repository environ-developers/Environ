!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
! Copyright (C) Quantum ESPRESSO (www.quantum-espresso.org)
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
! Authors: Modified by Edan Bainglass
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE env_fft_buffers
    !------------------------------------------------------------------------------------
    !
    USE env_fft_param, ONLY: DP
    !
    USE env_types_fft, ONLY: env_fft_type_descriptor
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    SAVE
    !
    INTEGER :: current_size = 0
    !
    COMPLEX(DP), ALLOCATABLE :: dev_space_fftparallel(:)
    COMPLEX(DP), ALLOCATABLE :: dev_space_scatter_dblbuffer(:)
    COMPLEX(DP), ALLOCATABLE :: pin_space_scatter_dblbuffer(:)
    COMPLEX(DP), ALLOCATABLE :: pin_space_scatter_in(:)
    COMPLEX(DP), ALLOCATABLE :: pin_space_scatter_out(:)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: dev_space_fftparallel
    attributes(DEVICE) :: dev_space_scatter_dblbuffer
    attributes(PINNED) :: pin_space_scatter_in
    attributes(PINNED) :: pin_space_scatter_out
    attributes(PINNED) :: pin_space_scatter_dblbuffer
    !
#endif
    !------------------------------------------------------------------------------------
    ! Keeping imports private
    !
    PRIVATE :: DP, env_fft_type_descriptor
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_check_fft_buffers_size(desc, howmany)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), INTENT(IN) :: desc
        INTEGER, OPTIONAL, INTENT(IN) :: howmany
        !
        INTEGER :: howmany_
        INTEGER :: info
        !
        CHARACTER(LEN=80) :: sub_name = 'env_check_fft_buffers_size'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(howmany)) THEN
            howmany_ = howmany
        ELSE
            howmany_ = 1
        END IF
        !
        IF (current_size < desc%nnr * howmany_) THEN
            !
            current_size = desc%nnr * howmany_
            !
            IF (ALLOCATED(dev_space_fftparallel)) DEALLOCATE (dev_space_fftparallel)
            !
            IF (ALLOCATED(pin_space_scatter_in)) DEALLOCATE (pin_space_scatter_in)
            !
            IF (ALLOCATED(pin_space_scatter_out)) DEALLOCATE (pin_space_scatter_out)
            !
            ALLOCATE (dev_space_fftparallel(current_size), STAT=info)
            !
            IF (info /= 0) CALL env_errore(sub_name, 'Allocation failed', 1)
            !
            ALLOCATE (pin_space_scatter_in(current_size), STAT=info)
            !
            IF (info /= 0) CALL env_errore(sub_name, 'Allocation failed', 2)
            !
            ALLOCATE (pin_space_scatter_out(current_size), STAT=info)
            !
            IF (info /= 0) CALL env_errore(sub_name, 'Allocation failed', 3)
            !
            !----------------------------------------------------------------------------
            ! Slab decomposition implements double buffering
            !
            IF ((.NOT. desc%use_pencil_decomposition)) THEN
                !
                IF (ALLOCATED(dev_space_scatter_dblbuffer)) &
                    DEALLOCATE (dev_space_scatter_dblbuffer)
                !
                IF (ALLOCATED(pin_space_scatter_dblbuffer)) &
                    DEALLOCATE (pin_space_scatter_dblbuffer)
                !
                ALLOCATE (dev_space_scatter_dblbuffer(current_size), STAT=info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'Allocation failed', 4)
                !
                ALLOCATE (pin_space_scatter_dblbuffer(current_size), STAT=info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'Allocation failed', 5)
                !
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_check_fft_buffers_size
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_deallocate_fft_buffers()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        current_size = 0
        !
        IF (ALLOCATED(dev_space_fftparallel)) DEALLOCATE (dev_space_fftparallel)
        !
        IF (ALLOCATED(pin_space_scatter_in)) DEALLOCATE (pin_space_scatter_in)
        !
        IF (ALLOCATED(pin_space_scatter_out)) DEALLOCATE (pin_space_scatter_out)
        !
        IF (ALLOCATED(dev_space_scatter_dblbuffer)) &
            DEALLOCATE (dev_space_scatter_dblbuffer)
        !
        IF (ALLOCATED(pin_space_scatter_dblbuffer)) &
            DEALLOCATE (pin_space_scatter_dblbuffer)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_deallocate_fft_buffers
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_fft_buffers
!----------------------------------------------------------------------------------------
