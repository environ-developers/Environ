!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 3.0
!
!     Environ 3.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 3.0 is distributed in the hope that it will be useful,
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
    !
    USE env_base_input
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
        PROCEDURE :: everything => clean_all
        !
        PROCEDURE :: first => clean_physical_first
        PROCEDURE :: second => clean_physical_second
        PROCEDURE :: numerical => clean_numerical
        !
        PROCEDURE, NOPASS :: input => deallocate_input_registers
        PROCEDURE, NOPASS :: debug => close_debug_file
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
    SUBROUTINE clean_all(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_destructor), INTENT(IN) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%first()
        !
        CALL this%second()
        !
        CALL this%numerical()
        !
        CALL this%input()
        !
        CALL this%debug()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE clean_all
    !------------------------------------------------------------------------------------
    !>
    !! Clean up all the Environ related allocated variables, and call clean up
    !! subroutines of specific Environ modules.
    !!
    !! The structure of this subroutine mirrors the one of init_environ subroutines
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE clean_physical_first(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_destructor), TARGET, INTENT(IN) :: this
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (main => this%main, &
                   setup => this%main%setup)
            !
            !----------------------------------------------------------------------------
            !
            IF (ASSOCIATED(main%vzero%cell)) CALL main%vzero%destroy()
            !
            IF (ASSOCIATED(main%dvtot%cell)) CALL main%dvtot%destroy()
            !
            !----------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------
            !
            CALL main%system_charges%destroy()
            !
            CALL main%environment_charges%destroy()
            !
            IF (setup%lexternals) CALL main%externals%destroy()
            !
            IF (setup%lstatic) CALL main%static%destroy()
            !
            IF (setup%lelectrolyte) THEN
                !
                CALL main%electrolyte%destroy()
                !
                DEALLOCATE (main%electrolyte%boundary)
            END IF
            !
            IF (setup%lsemiconductor) CALL main%semiconductor%destroy()
            !
            IF (setup%laddcharges) CALL main%additional_charges%destroy()
            !
            !----------------------------------------------------------------------------
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
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE clean_physical_first
    !------------------------------------------------------------------------------------
    !>
    !! Clean up all the Environ-related allocated variables and call clean up
    !! subroutines of specific Environ modules. These are quantities that may
    !! be needed by TDDFPT, thus may need to be cleaned later
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE clean_physical_second(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_destructor), INTENT(IN) :: this
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (main => this%main, &
                   setup => this%main%setup)
            !
            !----------------------------------------------------------------------------
            !
            IF (ASSOCIATED(main%velectrostatic%cell)) CALL main%velectrostatic%destroy()
            !
            !----------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------
            !
            IF (setup%lsolvent) THEN
                !
                CALL main%solvent%destroy()
                !
                DEALLOCATE (main%solvent)
            END IF
            !
            !----------------------------------------------------------------------------
            !
            main%initialized = .FALSE.
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE clean_physical_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE clean_numerical(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_destructor), INTENT(IN) :: this
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (setup => this%main%setup)
            !
            IF (setup%lelectrostatic) THEN
                !
                CALL setup%outer%destroy()
                !
                CALL setup%reference%destroy()
                !
            END IF
            !
            IF (setup%has_numerical_setup) THEN
                !
                IF (setup%lelectrostatic) CALL setup%reference_container%destroy()
                !
                CALL setup%outer_container%destroy()
                !
                IF (setup%need_inner) CALL setup%inner_container%destroy()
                !
            END IF
            !
            CALL setup%system_cell%destroy()
            !
            CALL setup%environment_cell%destroy()
            !
            CALL setup%mapping%destroy()
            !
        END ASSOCIATE
        !
        CALL env_deallocate_buffers()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE clean_numerical
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE close_debug_file()
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
    END SUBROUTINE close_debug_file
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE deallocate_input_registers()
        !--------------------------------------------------------------------------------
        !
        DEALLOCATE (atomicspread, solvationrad, corespread)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE deallocate_input_registers
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_destructor
!----------------------------------------------------------------------------------------
