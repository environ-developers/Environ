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
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Edan Bainglass     (Department of Physics, UNT)
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
    USE class_environ
    USE class_setup
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
        TYPE(environ_main), POINTER :: main => NULL()
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: everything => environ_clean
        PROCEDURE :: first => environ_clean_first
        PROCEDURE :: second => environ_clean_second
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_destructor
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Clean up all the Environ related allocated variables, and call
    !! clean up subroutines of specific Environ modules.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_clean(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_destructor), INTENT(IN) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'environ_clean'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%first()
        !
        CALL this%second()
        !
        CALL env_deallocate_buffers()
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
    SUBROUTINE environ_clean_first(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_destructor), TARGET, INTENT(IN) :: this
        !
        TYPE(environ_main), POINTER :: main
        TYPE(environ_setup), POINTER :: setup
        !
        CHARACTER(LEN=80) :: sub_name = 'environ_clean_first'
        !
        !--------------------------------------------------------------------------------
        !
        main => this%main
        setup => main%setup
        !
        !--------------------------------------------------------------------------------
        ! Deallocate environment variables
        !
        IF (ASSOCIATED(main%vzero%cell)) CALL main%vzero%destroy()
        !
        IF (ASSOCIATED(main%dvtot%cell)) CALL main%dvtot%destroy()
        !
        !--------------------------------------------------------------------------------
        ! base_environ variables
        !
        IF (setup%lelectrostatic .AND. ASSOCIATED(main%vreference%cell)) &
            CALL main%vreference%destroy()
        !
        IF (setup%lsoftcavity .AND. ASSOCIATED(main%vsoftcavity%cell)) &
            CALL main%vsoftcavity%destroy()
        !
        IF (setup%lconfine .AND. ASSOCIATED(main%vconfine%cell)) &
            CALL main%vconfine%destroy()
        !
        !--------------------------------------------------------------------------------
        ! Destroy derived types which were allocated in input
        !
        IF (setup%lelectrostatic .OR. setup%lconfine) THEN
            !
            CALL main%system_charges%destroy()
            !
            CALL main%environment_charges%destroy()
            !
        END IF
        !
        IF (setup%lexternals) CALL main%externals%destroy()
        !
        IF (setup%lstatic) CALL main%static%destroy()
        !
        IF (setup%lelectrolyte) CALL main%electrolyte%destroy()
        !
        IF (setup%lsemiconductor) CALL main%semiconductor%destroy()
        !
        IF (setup%laddcharges) CALL main%additional_charges%destroy()
        !
        CALL main%system_electrons%destroy()
        !
        CALL main%system_ions%destroy()
        !
        CALL main%system_system%destroy()
        !
        CALL main%environment_electrons%destroy()
        !
        CALL main%environment_ions%destroy()
        !
        CALL main%environment_system%destroy()
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
    SUBROUTINE environ_clean_second(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_destructor), INTENT(IN) :: this
        !
        LOGICAL :: opnd
        !
        TYPE(environ_main), POINTER :: main
        TYPE(environ_setup), POINTER :: setup
        !
        CHARACTER(LEN=80) :: sub_name = 'environ_clean_second'
        !
        !--------------------------------------------------------------------------------
        !
        main => this%main
        setup => main%setup
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
            IF (ASSOCIATED(main%velectrostatic%cell)) &
                CALL main%velectrostatic%destroy()
            !
            CALL setup%outer%destroy()
            !
            CALL setup%reference%destroy()
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Destroy derived types which were allocated in input
        !
        IF (setup%loptical) THEN
            !
            CALL main%environment_response_charges%destroy()
            !
            CALL main%environment_response_electrons%destroy()
            !
            CALL main%system_response_charges%destroy()
            !
            CALL main%system_response_electrons%destroy()
            !
            CALL main%optical%destroy()
            !
        END IF
        !
        IF (setup%lsolvent) CALL main%solvent%destroy()
        !
        !--------------------------------------------------------------------------------
        ! Destroy cells
        !
        CALL setup%system_cell%destroy()
        !
        CALL setup%environment_cell%destroy()
        !
        CALL setup%mapping%destroy()
        !
        !--------------------------------------------------------------------------------
        !
        main%initialized = .FALSE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_clean_second
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_destructor
!----------------------------------------------------------------------------------------
