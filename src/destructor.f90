!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
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
! Authors: Edan Bainglass (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_destructor
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE env_global_objects, ONLY: setup, env
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
    TYPE, PUBLIC :: environ_destructor
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, NOPASS :: all => environ_clean
        PROCEDURE, NOPASS :: first => environ_clean_first
        PROCEDURE, NOPASS :: second => environ_clean_second
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_destructor
    !------------------------------------------------------------------------------------
    !
    TYPE(environ_destructor), PUBLIC, SAVE :: clean
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Clean up all the Environ related allocated variables, and call
    !! clean up subroutines of specific Environ modules.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_clean()
        !--------------------------------------------------------------------------------
        !
        CALL environ_clean_first()
        !
        CALL environ_clean_second()
        !
        CALL env_deallocate_mp_buffers()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_clean
    !------------------------------------------------------------------------------------
    !>
    !! Clean up all the Environ related allocated variables, and call clean up
    !! subroutines of specific Environ modules.
    !!
    !! The structure of this subroutine mirrors the one of init_environ subroutines
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_clean_first()
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
        ! Deallocate environment variables
        !
        IF (ASSOCIATED(env%vzero%cell)) CALL env%vzero%destroy()
        !
        IF (ASSOCIATED(env%dvtot%cell)) CALL env%dvtot%destroy()
        !
        !--------------------------------------------------------------------------------
        ! base_environ variables
        !
        IF (setup%lelectrostatic .AND. ASSOCIATED(env%vreference%cell)) &
            CALL env%vreference%destroy()
        !
        IF (setup%lsoftcavity .AND. ASSOCIATED(env%vsoftcavity%cell)) &
            CALL env%vsoftcavity%destroy()
        !
        IF (setup%lconfine .AND. ASSOCIATED(env%vconfine%cell)) &
            CALL env%vconfine%destroy()
        !
        !--------------------------------------------------------------------------------
        ! Destroy derived types which were allocated in input
        !
        IF (setup%lelectrostatic .OR. setup%lconfine) THEN
            !
            CALL env%system_charges%destroy()
            !
            CALL env%environment_charges%destroy()
            !
        END IF
        !
        IF (setup%lexternals) CALL env%externals%destroy()
        !
        IF (setup%lstatic) CALL env%static%destroy()
        !
        IF (setup%lelectrolyte) CALL env%electrolyte%destroy()
        !
        IF (setup%lsemiconductor) CALL env%semiconductor%destroy()
        !
        CALL env%system_electrons%destroy()
        !
        CALL env%system_ions%destroy()
        !
        CALL env%system_system%destroy()
        !
        CALL env%environment_electrons%destroy()
        !
        CALL env%environment_ions%destroy()
        !
        CALL env%environment_system%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_clean_first
    !------------------------------------------------------------------------------------
    !>
    !! Clean up all the Environ-related allocated variables and call clean up
    !! subroutines of specific Environ modules. These are quantities that may
    !! be needed by TDDFPT, thus may need to be cleaned later
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_clean_second()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL :: opnd
        !
        !--------------------------------------------------------------------------------
        !
        INQUIRE (unit=io%debug_unit, opened=opnd)
        !
        IF (opnd) CLOSE (unit=io%debug_unit)
        !
        !--------------------------------------------------------------------------------
        ! base_environ variables
        !
        IF (setup%lelectrostatic) THEN
            !
            IF (ASSOCIATED(env%velectrostatic%cell)) CALL env%velectrostatic%destroy()
            !
            CALL electrostatic_clean()
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Destroy derived types which were allocated in input
        !
        IF (setup%loptical) THEN
            !
            CALL env%environment_response_charges%destroy()
            !
            CALL env%environment_response_electrons%destroy()
            !
            CALL env%system_response_charges%destroy()
            !
            CALL env%system_response_electrons%destroy()
            !
            CALL env%optical%destroy()
            !
        END IF
        !
        IF (setup%lsolvent) CALL env%solvent%destroy()
        !
        IF (setup%lboundary) CALL setup%derivatives%destroy()
        !
        IF (setup%ldoublecell) THEN
            !
            CALL setup%mapping%destroy()
            !
            CALL setup%environment_cell%destroy()
            !
        END IF
        !
        CALL setup%system_cell%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_clean_second
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
    SUBROUTINE electrostatic_clean()
        !--------------------------------------------------------------------------------
        !
        CALL setup%outer%destroy()
        !
        CALL setup%reference%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE electrostatic_clean
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_destructor
!----------------------------------------------------------------------------------------