! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.1
!
!    Environ 1.1 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.1 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------------------
!>
!! This module provides output subroutines for Environ, including summary
!! of input parameters, timings, print of energy contributions, and print
!! subroutines for Environ derived data types
!!
!----------------------------------------------------------------------------------------
MODULE environ_output
    !------------------------------------------------------------------------------------
    !
    USE env_io, ONLY: env_find_free_unit
    USE env_mp, ONLY: env_mp_sum
    !
    USE env_fft_types, ONLY: env_fft_type_descriptor
    USE env_scatter_mod, ONLY: env_gather_grid
    !
    USE modules_constants, ONLY: DP, amu_si, bohr_radius_si, rydberg_si, RYTOEV
    !
    USE cell_types, ONLY: environ_cell
    USE representation_types
    USE physical_types
    !
    USE environ_base, ONLY: lelectrostatic, eelectrostatic, lsurface, esurface, &
                            lvolume, evolume, lelectrolyte, eelectrolyte, lconfine, &
                            econfine, deenviron, lsmearedions, potential_shift, &
                            environ_thr, lsolvent, solvent, env_static_permittivity, &
                            env_optical_permittivity, env_surface_tension, &
                            env_pressure, lelectrostatic, ltddfpt, derivatives, &
                            lsemiconductor, system_ions
    !
    USE electrostatic_base, ONLY: need_pbc_correction, outer
    USE core_base, ONLY: lfd, fd
    !
    USE utils_density, ONLY: init_environ_density, destroy_environ_density
    !
    USE tools_math, ONLY: integrate_environ_density
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    SAVE
    !
    LOGICAL :: ionode = .TRUE.
    INTEGER :: ionode_id
    !
    LOGICAL :: lstdout ! whether environ can print on standard output
    !
    INTEGER :: comm ! WE MAY NEED A SECOND COMMUNICATOR FOR IMAGE PARALLELIZATION
    !
    INTEGER :: program_unit
    INTEGER :: environ_unit
    INTEGER :: verbose
    INTEGER :: depth = 1
    !
    CHARACTER(LEN=2) :: prog
    !
    INTEGER, PARAMETER :: nbibliography = 6
    CHARACTER(LEN=80) :: bibliography(nbibliography)
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: ionode, ionode_id, lstdout, comm, program_unit, environ_unit, verbose, &
              prog, set_environ_output, update_output_program_unit, &
              print_environ_density, print_environ_gradient, print_environ_hessian, &
              print_environ_functions, print_environ_iontype, print_environ_ions, &
              print_environ_electrons, print_environ_externals, print_environ_charges, &
              print_environ_system, print_environ_boundary, print_environ_dielectric, &
              print_environ_electrolyte, environ_print_energies, &
              environ_print_potential_shift, environ_print_potential_warning, &
              environ_summary, environ_clock, write_cube
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Populate bibliography entries
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_bibliography()
        !--------------------------------------------------------------------------------
        !
        bibliography(1) = "O. Andreussi, I. Dabo and N. Marzari, J. &
                          &Chem. Phys. 136, 064102 (2012)"
        bibliography(2) = "I. Timrov, O. Andreussi, A. Biancardi, N. &
                          &Marzari, and S. Baroni, J. Chem. Phys. &
                          &142, 034111 (2015)"
        bibliography(3) = "O. Andreussi, I. Dabo and N. Marzari, J. &
                          &Chem. Phys. 136, 064102 (2012)"
        bibliography(4) = "O. Andreussi, I. Dabo and N. Marzari, J. &
                          &Chem. Phys. 136, 064102 (2012)"
        bibliography(5) = "O. Andreussi, N.G. Hoermann, F. Nattino, &
                          &G. Fisicaro, S. Goedecker, and N. Marzari, &
                          &J. Chem. Theory Comput. 15, 1996 (2019)"
        bibliography(6) = "F. Nattino, M. Truscott, N. Marzari, and &
                          &O. Andreussi, J. Chem. Phys. 150, 041722 &
                          &(2019)"
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_bibliography
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_environ_output(prog_, ionode_, ionode_id_, comm_, program_unit_)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: prog_
        LOGICAL, INTENT(IN) :: ionode_
        INTEGER, INTENT(IN) :: ionode_id_
        INTEGER, INTENT(IN) :: comm_
        INTEGER, INTENT(IN) :: program_unit_
        !
        !--------------------------------------------------------------------------------
        !
        CALL set_bibliography()
        !
        ionode = ionode_
        ionode_id = ionode_id_
        comm = comm_
        !
        program_unit = program_unit_
        environ_unit = env_find_free_unit()
        !
        prog = prog_(1:2)
        !
        SELECT CASE (prog)
        CASE ('PW', 'pw')
            lstdout = .TRUE.
        CASE ('CP', 'cp')
            lstdout = .TRUE.
        CASE DEFAULT
            lstdout = .FALSE.
        END SELECT
        !
        lstdout = lstdout .AND. ionode
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_environ_output
    !------------------------------------------------------------------------------------
    !>
    !! Sets the output file target
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_output_program_unit(program_unit_)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: program_unit_
        !
        !--------------------------------------------------------------------------------
        !
        program_unit = program_unit_
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_output_program_unit
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_cell(cell, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN ! environ output file has not been opened
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN ! nothing to output
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose .AND. ionode) WRITE (UNIT=environ_unit, FMT=1000)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1001) cell%alat, cell%omega
            !
            IF (verbosity >= 3) THEN
                !
                IF (ionode) WRITE (UNIT=environ_unit, FMT=1002) cell%at
                !
                IF (ionode) &
                    WRITE (UNIT=environ_unit, FMT=1003) &
                    cell%dfft%nr1, cell%dfft%nr2, cell%dfft%nr3
                !
                IF (ionode) &
                    WRITE (UNIT=environ_unit, FMT=1004) &
                    cell%ntot, cell%nnr, cell%domega
                !
            END IF
            !
        END IF
        !
        RETURN
        !
1000    FORMAT(/, 4('%'), ' CELL ', 70('%'))
        !
1001    FORMAT(1X, 'lattice spacing            = ', F12.6, ' ' &
               /, 1X, 'cell volume                = ', F12.6, ' ')
        !
1002    FORMAT(1X, 'simulation cell axes       = ', 3(F12.6), ' ' &
               /, 1X, '                             ', 3(F12.6), ' ' &
               /, 1X, '                             ', 3(F12.6), ' ')
        !
1003    FORMAT(1X, 'real space grid dim.s      = ', 3I4, ' ')
        !
1004    FORMAT(1X, 'total size of grid         = ', I10, ' ' &
               /, 1X, 'size of r-space per proc.  = ', I10, ' ' &
               /, 1X, 'finite element volume      = ', F12.6, ' ')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_density(density, local_verbose, local_depth, local_ions)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: density
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        TYPE(environ_ions), INTENT(IN), OPTIONAL :: local_ions
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        REAL(DP) :: integral
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN ! environ output file has not been opened
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN ! nothing to output
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose .AND. ionode) WRITE (UNIT=environ_unit, FMT=1100)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1101) ADJUSTL(density%label)
            !
            integral = integrate_environ_density(density)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1102) integral
            !
            ! #TODO ADD MAXVAL AND MINVAL
            !
            IF (verbosity >= 3) THEN
                !
                CALL print_environ_cell(density%cell, passed_verbosity, passed_depth)
                !
                IF (PRESENT(local_ions)) THEN
                    !
                    IF (ionode) WRITE (environ_unit, *) 'present local ions'
                    !
                    CALL write_cube(density, local_ions)
                    !
                ELSE IF (system_ions%initialized) THEN
                    !
                    IF (ionode) WRITE (environ_unit, *) 'using stored ions'
                    !
                    CALL write_cube(density, system_ions)
                    !
                ELSE
                    !
                    IF (ionode) WRITE (environ_unit, *) 'no ions'
                    !
                    CALL write_cube(density)
                    !
                END IF
                !
            END IF
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
1100    FORMAT(/, 4('%'), ' DENSITY ', 67('%'))
1101    FORMAT(1X, 'density label              = ', A80)
1102    FORMAT(1X, 'integral of density        = ', G20.10)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_gradient(gradient, local_verbose, local_depth, local_ions)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_gradient), INTENT(IN) :: gradient
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        TYPE(environ_ions), INTENT(IN), OPTIONAL :: local_ions
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density) :: dens
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        REAL(DP) :: integral
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN ! environ output file has not been opened
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN ! nothing to output
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose .AND. ionode) WRITE (UNIT=environ_unit, FMT=1200)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1201) ADJUSTL(gradient%label)
            !
            integral = integrate_environ_density(gradient%modulus)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1202) integral
            !
            ! #TODO ADD MAXVAL AND MINVAL
            !
            IF (verbosity >= 3) THEN
                !
                CALL print_environ_cell(gradient%cell, passed_verbosity, passed_depth)
                !
                IF (PRESENT(local_ions)) THEN
                    CALL write_cube(gradient%modulus, local_ions)
                ELSE IF (ASSOCIATED(system_ions%tau)) THEN
                    CALL write_cube(gradient%modulus, system_ions)
                ELSE
                    CALL write_cube(gradient%modulus)
                END IF
                !
            END IF
            !
            IF (verbosity >= 4) THEN
                cell => gradient%cell
                !
                CALL init_environ_density(cell, dens)
                !
                dens%label = TRIM(ADJUSTL(gradient%label))//'_x'
                dens%of_r(:) = gradient%of_r(1, :)
                !
                CALL print_environ_density(dens, passed_verbosity, passed_depth)
                !
                dens%label = TRIM(ADJUSTL(gradient%label))//'_y'
                dens%of_r(:) = gradient%of_r(2, :)
                !
                CALL print_environ_density(dens, passed_verbosity, passed_depth)
                !
                dens%label = TRIM(ADJUSTL(gradient%label))//'_z'
                dens%of_r(:) = gradient%of_r(3, :)
                !
                CALL print_environ_density(dens, passed_verbosity, passed_depth)
                !
                CALL destroy_environ_density(dens)
                !
            END IF
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
1200    FORMAT(/, 4('%'), ' GRADIENT ', 66('%'))
1201    FORMAT(1X, 'gradient label             = ', A80)
1202    FORMAT(1X, 'integral of square modul.  = ', G20.10)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_hessian(hessian, local_verbose, local_depth, local_ions)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_hessian), INTENT(IN) :: hessian
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        TYPE(environ_ions), INTENT(IN), OPTIONAL :: local_ions
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density) :: dens
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        REAL(DP) :: integral
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN ! environ output file has not been opened
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN ! nothing to output
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose .AND. ionode) WRITE (UNIT=environ_unit, FMT=1250)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1251) ADJUSTL(hessian%label)
            !
            integral = integrate_environ_density(hessian%laplacian)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1252) integral
            !
            ! #TODO ADD MAXVAL AND MINVAL
            !
            IF (verbosity >= 3) THEN
                !
                CALL print_environ_cell(hessian%cell, passed_verbosity, passed_depth)
                !
                IF (PRESENT(local_ions)) THEN
                    CALL write_cube(hessian%laplacian, local_ions)
                ELSE IF (ASSOCIATED(system_ions%tau)) THEN
                    CALL write_cube(hessian%laplacian, system_ions)
                ELSE
                    CALL write_cube(hessian%laplacian)
                END IF
                !
            END IF
            !
            IF (verbosity >= 4) THEN
                cell => hessian%cell
                !
                CALL init_environ_density(cell, dens)
                !
                dens%label = TRIM(ADJUSTL(hessian%label))//'_xx'
                dens%of_r(:) = hessian%of_r(1, 1, :)
                !
                CALL print_environ_density(dens, passed_verbosity, passed_depth)
                !
                dens%label = TRIM(ADJUSTL(hessian%label))//'_xy'
                dens%of_r(:) = hessian%of_r(1, 2, :)
                !
                CALL print_environ_density(dens, passed_verbosity, passed_depth)
                !
                dens%label = TRIM(ADJUSTL(hessian%label))//'_xz'
                dens%of_r(:) = hessian%of_r(1, 3, :)
                !
                CALL print_environ_density(dens, passed_verbosity, passed_depth)
                !
                dens%label = TRIM(ADJUSTL(hessian%label))//'_yy'
                dens%of_r(:) = hessian%of_r(2, 2, :)
                !
                CALL print_environ_density(dens, passed_verbosity, passed_depth)
                !
                dens%label = TRIM(ADJUSTL(hessian%label))//'_yz'
                dens%of_r(:) = hessian%of_r(2, 3, :)
                !
                CALL print_environ_density(dens, passed_verbosity, passed_depth)
                !
                dens%label = TRIM(ADJUSTL(hessian%label))//'_zz'
                dens%of_r(:) = hessian%of_r(3, 3, :)
                !
                CALL print_environ_density(dens, passed_verbosity, passed_depth)
                !
                CALL destroy_environ_density(dens)
                !
            END IF
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
1250    FORMAT(/, 4('%'), ' HESSIAN ', 67('%'))
1251    FORMAT(1X, 'hessian label              = ', A80)
1252    FORMAT(1X, 'integral of laplacian      = ', G20.10)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_hessian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_functions(nfunctions, functions, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nfunctions
        TYPE(environ_functions), INTENT(IN) :: functions(nfunctions)
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        INTEGER :: ifunctions
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_functions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN ! environ output file has not been opened
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN ! nothing to output
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose .AND. ionode) WRITE (UNIT=environ_unit, FMT=1300)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1301) nfunctions
            !
            DO ifunctions = 1, nfunctions
                !
                SELECT CASE (functions(ifunctions)%type_)
                CASE (1)
                    !
                    IF (ionode) &
                        WRITE (UNIT=environ_unit, FMT=1302) ifunctions, &
                        functions(ifunctions)%dim, functions(ifunctions)%axis, &
                        functions(ifunctions)%spread, functions(ifunctions)%volume
                    !
                CASE (2)
                    !
                    IF (ionode) &
                        WRITE (UNIT=environ_unit, FMT=1303) ifunctions, &
                        functions(ifunctions)%dim, functions(ifunctions)%axis, &
                        functions(ifunctions)%width, functions(ifunctions)%spread, &
                        functions(ifunctions)%volume
                    !
                CASE (3)
                    !
                    IF (ionode) &
                        WRITE (UNIT=environ_unit, FMT=1304) ifunctions, &
                        functions(ifunctions)%spread, functions(ifunctions)%volume
                    !
                CASE (4)
                    !
                    IF (ionode) &
                        WRITE (UNIT=environ_unit, FMT=1305) ifunctions, &
                        functions(ifunctions)%dim, functions(ifunctions)%axis, &
                        functions(ifunctions)%width, functions(ifunctions)%spread, &
                        functions(ifunctions)%volume
                    !
                CASE (5)
                    !
                    IF (ionode) &
                        WRITE (UNIT=environ_unit, FMT=1306) ifunctions, &
                        functions(ifunctions)%dim, functions(ifunctions)%axis, &
                        functions(ifunctions)%width, functions(ifunctions)%spread, &
                        functions(ifunctions)%volume
                    !
                CASE DEFAULT
                    CALL env_errore(sub_name, 'Unexpected function type', 1)
                END SELECT
                !
                IF (verbosity >= 3 .AND. ionode) &
                    WRITE (UNIT=environ_unit, FMT=1307) functions(ifunctions)%pos
                !
            END DO
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
1300    FORMAT(/, 4('%'), ' FUNCTIONS ', 65('%'))
1301    FORMAT(1X, 'number of functions        = ', I10)
        !
1302    FORMAT(1X, 'Gaussian function, type    = 1 ' &
               /, 1X, 'number                     = ', I3, ' ' &
               /, 1X, 'dimensionality             = ', I1, ' ' &
               /, 1X, 'axis                       = ', I1, ' ' &
               /, 1X, 'spread                     = ', F12.6, ' ' &
               /, 1X, 'integral                   = ', F12.6, ' ')
        !
1303    FORMAT(1X, 'ERFC function, type        = 2 ' &
               /, 1X, 'number                     = ', I3, ' ' &
               /, 1X, 'dimensionality             = ', I1, ' ' &
               /, 1X, 'axis                       = ', I1, ' ' &
               /, 1X, 'width                      = ', F12.6, ' ' &
               /, 1X, 'spread                     = ', F12.6, ' ' &
               /, 1X, 'integral                   = ', F12.6, ' ')
        !
1304    FORMAT(1X, 'exponential function, type = 3 ' &
               /, 1X, 'number                     = ', I3, ' ' &
               /, 1X, 'spread                     = ', F12.6, ' ' &
               /, 1X, 'integral                   = ', F12.6, ' ')
        !
1305    FORMAT(1X, 'Scaled ERFC function, type = 4 ' &
               /, 1X, 'number                     = ', I3, ' ' &
               /, 1X, 'dimensionality             = ', I1, ' ' &
               /, 1X, 'axis                       = ', I1, ' ' &
               /, 1X, 'width                      = ', F12.6, ' ' &
               /, 1X, 'spread                     = ', F12.6, ' ' &
               /, 1X, 'max value                  = ', F12.6, ' ')
        !
1306    FORMAT(1X, 'Scaled ERF function, type  = 5 ' &
               /, 1X, 'number                     = ', I3, ' ' &
               /, 1X, 'dimensionality             = ', I1, ' ' &
               /, 1X, 'axis                       = ', I1, ' ' &
               /, 1X, 'width                      = ', F12.6, ' ' &
               /, 1X, 'spread                     = ', F12.6, ' ' &
               /, 1X, 'max value                  = ', F12.6, ' ')
        !
1307    FORMAT(1X, 'position                   = ', 3(F14.7), ' ')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_functions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_iontype(ntyp, iontype, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ntyp
        TYPE(environ_iontype), INTENT(IN) :: iontype(ntyp)
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        INTEGER :: ityp
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_iontype'
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN ! environ output file has not been opened
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN ! nothing to output
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose .AND. ionode) WRITE (UNIT=environ_unit, FMT=1400)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1401) ntyp
            !
            DO ityp = 1, ntyp
                !
                IF (ionode) &
                    WRITE (UNIT=environ_unit, FMT=1402) ityp, iontype(ityp)%index, &
                    iontype(ityp)%label, iontype(ityp)%atmnum, iontype(ityp)%zv
                !
                IF (verbosity >= 3 .AND. ionode) &
                    WRITE (UNIT=environ_unit, FMT=1403) &
                    iontype(ityp)%atomicspread, iontype(ityp)%corespread, &
                    iontype(ityp)%solvationrad
                !
            END DO
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
1400    FORMAT(/, 4('%'), ' IONTYPES ', 66('%'))
1401    FORMAT(1X, 'number of ionic types      = ', I10)
        !
1402    FORMAT(1X, 'iontype = ', I2, ' index = ', I2, ' label = ', A3, &
               ' atomic num = ', I2, ' charge = ', F7.2)
        !
1403    FORMAT(1X, 'atomic spread = ', F7.2, ' core spread = ', F7.2, &
               ' solvation radius = ', F7.2)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_iontype
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_ions(ions, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_ions), INTENT(IN) :: ions
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_ions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN ! environ output file has not been opened
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN ! nothing to output
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose .AND. ionode) WRITE (UNIT=environ_unit, FMT=1500)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1501) ions%number
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1502) ions%center
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1503) ions%charge
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1505) ions%quadrupole_pc
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1504) ions%dipole
            !
            IF (ions%use_smeared_ions .AND. ionode) &
                WRITE (UNIT=environ_unit, FMT=1506) ions%quadrupole_gauss
            !
            DO i = 1, ions%number
                !
                IF (ionode) &
                    !
                    WRITE (UNIT=environ_unit, FMT=1507) &
                    i, ions%ityp(i), ions%tau(:, i) * ions%alat
                !
            END DO
            !
            IF (verbosity >= 3) THEN
                !
                CALL print_environ_iontype(ions%ntyp, ions%iontype, passed_verbosity, &
                                           passed_depth)
                !
                IF (ions%use_smeared_ions) THEN
                    !
                    CALL print_environ_density(ions%density, passed_verbosity, &
                                               passed_depth)
                    !
                    IF (verbosity >= 4) &
                        CALL print_environ_functions(ions%number, ions%smeared_ions, &
                                                     passed_verbosity, passed_depth)
                    !
                END IF
                !
                IF (ions%use_core_electrons) THEN
                    !
                    IF (verbosity >= 4) &
                        CALL print_environ_density(ions%core, passed_verbosity, &
                                                   passed_depth)
                    !
                    IF (verbosity >= 5) &
                        CALL print_environ_functions(ions%number, ions%core_electrons, &
                                                     passed_verbosity, passed_depth)
                    !
                END IF
                !
            END IF
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
1500    FORMAT(/, 4('%'), ' IONS ', 70('%'))
1501    FORMAT(1X, 'number of ions             = ', I10)
1502    FORMAT(1X, 'ionic center of charge     = ', 3(F14.7))
1503    FORMAT(1X, 'total ionic charge         = ', F14.7)
1504    FORMAT(1X, 'ionic dipole               = ', 3(F14.7))
1505    FORMAT(1X, 'ionic quadrupole (pc)      = ', 3(F14.7))
1506    FORMAT(1X, 'ionic quadrupole (gauss)   = ', 3(F14.7))
1507    FORMAT(1X, 'ion ', I3, ' type = ', I2, ' coordinates = ', 3(F14.7))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_ions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_electrons(electrons, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_electrons), INTENT(IN) :: electrons
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_electrons'
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN ! environ output file has not been opened
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN ! nothing to output
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose .AND. ionode) WRITE (UNIT=environ_unit, FMT=1600)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1601) electrons%number
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1602) electrons%charge
            !
            IF (verbosity >= 3) &
                CALL print_environ_density(electrons%density, passed_verbosity, &
                                           passed_depth)
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
1600    FORMAT(/, 4('%'), ' ELECTRONS ', 65('%'))
1601    FORMAT(1X, 'number of electrons        = ', I10)
1602    FORMAT(1X, 'total electronic charge    = ', F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_electrons
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_externals(externals, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_externals), INTENT(IN) :: externals
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_externals'
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN ! environ output file has not been opened
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN ! nothing to output
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose .AND. ionode) WRITE (UNIT=environ_unit, FMT=1700)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1701) externals%number
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1702) externals%charge
            !
            IF (verbosity >= 3) THEN
                !
                CALL print_environ_functions(externals%number, externals%functions, &
                                             passed_verbosity, passed_depth)
                !
                CALL print_environ_density(externals%density, passed_verbosity, &
                                           passed_depth)
                !
            END IF
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
1700    FORMAT(/, 4('%'), ' EXTERNALS ', 65('%'))
1701    FORMAT(1X, 'number of external charges = ', I10)
1702    FORMAT(1X, 'total external charge      = ', F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_externals
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_charges(charges, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_charges), INTENT(IN) :: charges
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_charges'
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN ! environ output file has not been opened
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN ! nothing to output
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose .AND. ionode) WRITE (UNIT=environ_unit, FMT=1800)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1801) charges%number
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1802) charges%charge
            !
            IF (verbosity >= 3) THEN
                !
                CALL print_environ_density(charges%density, passed_verbosity, &
                                           passed_depth)
                !
                IF (charges%include_ions) &
                    CALL print_environ_ions(charges%ions, passed_verbosity, &
                                            passed_depth)
                !
                IF (charges%include_electrons) &
                    CALL print_environ_electrons(charges%electrons, passed_verbosity, &
                                                 passed_depth)
                !
                IF (charges%include_externals) &
                    CALL print_environ_externals(charges%externals, passed_verbosity, &
                                                 passed_depth)
                !
                IF (charges%include_dielectric) &
                    CALL print_environ_dielectric(charges%dielectric, &
                                                  passed_verbosity, passed_depth)
                !
                IF (charges%include_electrolyte) &
                    CALL print_environ_electrolyte(charges%electrolyte, &
                                                   passed_verbosity, passed_depth)
                !
            END IF
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
1800    FORMAT(/, 4('%'), ' CHARGES ', 67('%'))
1801    FORMAT(1X, 'total number of charges    = ', I10)
1802    FORMAT(1X, 'total charge               = ', F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_system(system, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_system), INTENT(IN) :: system
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_system'
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN ! environ output file has not been opened
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN ! nothing to output
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose .AND. ionode) WRITE (UNIT=environ_unit, FMT=1900)
            !
            IF (system%ntyp == 0) THEN
                IF (ionode) WRITE (UNIT=environ_unit, FMT=1901)
            ELSE
                IF (ionode) WRITE (UNIT=environ_unit, FMT=1902) system%ntyp
            END IF
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1903) system%dim, system%axis
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=1904) system%pos, system%width
            !
            IF (verbosity >= 3) &
                CALL print_environ_ions(system%ions, passed_verbosity, passed_depth)
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
1900    FORMAT(/, 4('%'), ' SYSTEM ', 68('%'))
1901    FORMAT(1X, 'system is built from all present ionic types')
1902    FORMAT(1X, 'system is built from the first ', I3, ' ionic types')
        !
1903    FORMAT(1X, 'system defined dimension   = ', I2, ' ' &
               /, 1X, 'system defined axis        = ', I2, ' ')
        !
1904    FORMAT(1X, 'system center              = ', 3F14.7, ' ' &
               /, 1X, 'system width               = ', F14.7, ' ')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_system
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_boundary(boundary, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_boundary), INTENT(IN) :: boundary
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN ! environ output file has not been opened
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN ! nothing to output
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose .AND. ionode) WRITE (UNIT=environ_unit, FMT=2000)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=2001) boundary%label
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=2002) boundary%mode
            !
            IF (boundary%need_electrons) THEN
                !
                IF (ionode) WRITE (UNIT=environ_unit, FMT=2003) boundary%type_
                !
                SELECT CASE (boundary%type_)
                CASE (0)
                    !
                    IF (ionode) &
                        WRITE (UNIT=environ_unit, FMT=2004) &
                        boundary%rhozero, boundary%tbeta
                    !
                CASE (1)
                    !
                    IF (ionode) &
                        WRITE (UNIT=environ_unit, FMT=2005) &
                        boundary%rhomax, boundary%rhomin
                    !
                    IF (verbosity >= 3 .AND. ionode) &
                        WRITE (UNIT=environ_unit, FMT=2006) boundary%fact
                    !
                CASE (2)
                    !
                    IF (ionode) &
                        WRITE (UNIT=environ_unit, FMT=2007) &
                        boundary%rhomax, boundary%rhomin
                    !
                END SELECT
                !
                IF (verbosity >= 4) THEN
                    !
                    CALL print_environ_density(boundary%density, passed_verbosity, &
                                               passed_depth)
                    !
                    IF (boundary%need_ions) THEN
                        !
                        IF (ionode) WRITE (UNIT=environ_unit, FMT=2008)
                        !
                        IF (verbosity >= 4) &
                            CALL print_environ_density(boundary%ions%core)
                        !
                    END IF
                    !
                    IF (verbosity >= 4) &
                        CALL print_environ_electrons(boundary%electrons, &
                                                     passed_verbosity, passed_depth)
                    !
                END IF
                !
                IF (verbosity >= 5) &
                    CALL print_environ_density(boundary%dscaled, passed_verbosity, &
                                               passed_depth)
                !
                IF (verbosity >= 5) &
                    CALL print_environ_density(boundary%d2scaled, passed_verbosity, &
                                               passed_depth)
                !
            ELSE IF (boundary%need_ions) THEN
                !
                IF (ionode) &
                    WRITE (UNIT=environ_unit, FMT=2009) &
                    boundary%alpha, boundary%softness
                !
                IF (verbosity >= 3) &
                    CALL print_environ_functions(boundary%ions%number, &
                                                 boundary%soft_spheres, &
                                                 passed_verbosity, passed_depth)
                !
            ELSE IF (boundary%need_system) THEN
                !
                IF (ionode) &
                    WRITE (UNIT=environ_unit, FMT=2010) &
                    boundary%simple%pos, boundary%simple%width, &
                    boundary%simple%spread, boundary%simple%dim, &
                    boundary%simple%axis
                !
            END IF
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=2011) boundary%volume
            !
            IF (boundary%deriv >= 1 .AND. ionode) &
                WRITE (UNIT=environ_unit, FMT=2012) boundary%surface
            !
            IF (verbosity >= 4) &
                CALL print_environ_density(boundary%scaled, passed_verbosity, &
                                           passed_depth)
            !
            IF (boundary%solvent_aware) THEN
                !
                IF (ionode) &
                    WRITE (UNIT=environ_unit, FMT=2013) &
                    boundary%filling_threshold, boundary%filling_spread, &
                    boundary%solvent_probe%width, boundary%solvent_probe%spread
                !
                IF (verbosity >= 4) &
                    CALL print_environ_density(boundary%local, passed_verbosity, &
                                               passed_depth)
                !
                IF (verbosity >= 4) &
                    CALL print_environ_density(boundary%filling, passed_verbosity, &
                                               passed_depth)
                !
                IF (verbosity >= 5) &
                    CALL print_environ_density(boundary%dfilling, passed_verbosity, &
                                               passed_depth)
                !
                IF (verbosity >= 5) &
                    CALL print_environ_density(boundary%probe, passed_verbosity, &
                                               passed_depth)
                !
            END IF
            !
            IF (verbosity >= 5 .AND. boundary%deriv >= 1) &
                CALL print_environ_gradient(boundary%gradient, passed_verbosity, &
                                            passed_depth)
            IF (verbosity >= 5 .AND. boundary%deriv >= 2) &
                CALL print_environ_density(boundary%laplacian, passed_verbosity, &
                                           passed_depth)
            IF (verbosity >= 5 .AND. boundary%deriv == 3) &
                CALL print_environ_density(boundary%dsurface, passed_verbosity, &
                                           passed_depth)
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
2000    FORMAT(/, 4('%'), ' BOUNDARY ', 66('%'))
2001    FORMAT(1X, 'boundary label             = ', A20, ' ')
2002    FORMAT(1X, 'boundary mode              = ', A20, ' ')
        !
2003    FORMAT(1X, 'boundary is built as a function of a smooth density' &
               /, 1X, 'function type              = ', I2, ' ')
        !
2004    FORMAT(1X, 'using the Fattebert-Gygi function with parameters ' &
               /, 1X, 'rhozero                    = ', F14.7, ' ' &
               /, 1X, '2*beta                     = ', F14.7, ' ')
        !
2005    FORMAT(1X, 'using the optimal SCCS function with parameters ' &
               /, 1X, 'rhomax                     = ', F14.7, ' ' &
               /, 1X, 'rhomin                     = ', F14.7, ' ')
        !
2006    FORMAT(1X, 'log(rhomax/rhomin)         = ', F14.7, ' ')
        !
2007    FORMAT(1X, 'using the modified SCCS function with parameters ' &
               /, 1X, 'rhomax                     = ', F14.7, ' ' &
               /, 1X, 'rhomin                     = ', F14.7, ' ')
        !
2008    FORMAT(1X, 'adding fictitious core-electrons')
        !
2009    FORMAT(1X, 'boundary is built from soft-spheres centered on ionic positions' &
               /, 1X, 'solvent-dependent scaling  = ', F14.7, ' ' &
               /, 1X, 'softness parameter         = ', F14.7, ' ')
        !
2010    FORMAT(1X, 'boundary is built as an analytic function centered on system position' &
               /, 1X, 'center of the boundary     = ', 3F14.7, ' ' &
               /, 1X, 'distance from the center   = ', F14.7, ' ' &
               /, 1X, 'spread of the interface    = ', F14.7, ' ' &
               /, 1X, 'dimensionality             = ', I2, ' ' &
               /, 1X, 'axis                       = ', I2, ' ')
        !
2011    FORMAT(1X, 'volume of the QM region    = ', F14.7, ' ')
2012    FORMAT(1X, 'surface of the QM region   = ', F14.7, ' ')
        !
2013    FORMAT(1X, 'using solvent-aware boundary           ' &
               /, 1X, 'filling threshold          = ', F14.7, ' ' &
               /, 1X, 'filling spread             = ', F14.7, ' ' &
               /, 1X, 'solvent radius x rad scale = ', F14.7, ' ' &
               /, 1X, 'spread of solvent probe    = ', F14.7, ' ')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_dielectric(dielectric, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_dielectric), INTENT(IN) :: dielectric
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_dielectric'
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN ! environ output file has not been opened
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN ! nothing to output
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose .AND. ionode) WRITE (UNIT=environ_unit, FMT=2100)
            !
            IF (dielectric%nregions == 0) THEN
                IF (ionode) WRITE (UNIT=environ_unit, FMT=2101) dielectric%constant
            ELSE
                !
                IF (ionode) &
                    WRITE (UNIT=environ_unit, FMT=2102) &
                    dielectric%constant, dielectric%nregions
                !
                CALL print_environ_functions(dielectric%nregions, dielectric%regions, &
                                             passed_verbosity, passed_depth)
                !
                IF (verbosity >= 4) &
                    CALL print_environ_density(dielectric%background, &
                                               passed_verbosity, passed_depth)
                !
            END IF
            !
            IF (verbosity >= 3) &
                CALL print_environ_density(dielectric%density, passed_verbosity, &
                                           passed_depth)
            !
            IF (verbosity >= 3) &
                CALL print_environ_density(dielectric%epsilon, passed_verbosity, &
                                           passed_depth)
            !
            IF (verbosity >= 5) &
                CALL print_environ_density(dielectric%depsilon, passed_verbosity, &
                                           passed_depth)
            !
            IF (ionode) &
                WRITE (UNIT=environ_unit, FMT=2103) &
                dielectric%need_gradient, dielectric%need_factsqrt
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=2104) dielectric%charge
            !
            IF (verbosity >= 5) THEN
                !
                CALL print_environ_gradient(dielectric%gradlog, passed_verbosity, &
                                            passed_depth)
                !
                IF (dielectric%need_gradient) &
                    CALL print_environ_gradient(dielectric%gradient, passed_verbosity, &
                                                passed_depth)
                !
                IF (dielectric%need_factsqrt) &
                    CALL print_environ_density(dielectric%factsqrt, passed_verbosity, &
                                               passed_depth)
                !
            END IF
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
2100    FORMAT(/, 4('%'), ' DIELECTRIC ', 65('%'))
        !
2101    FORMAT(1X, 'dielectric build on homogeneous background' &
               /, 1X, 'environment bulk permitt.  = ', F14.7, ' ')
        !
2102    FORMAT(1X, 'dielectric build in the presence of dielectric regions' &
               /, 1X, 'environment bulk permitt.  = ', F14.7, ' ' &
               /, 1X, 'number of dielec. regions  = ', I4, ' ')
        !
2103    FORMAT(1X, 'dielectric flags' &
               /, 1X, 'need gradient              = ', L2, ' ' &
               /, 1X, 'need factor depend. sqrt   = ', L2, ' ')
        !
2104    FORMAT(1X, 'total dielectric charge    = ', F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_dielectric
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_electrolyte(electrolyte, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_electrolyte), INTENT(IN) :: electrolyte
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth, ityp
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_electrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN ! environ output file has not been opened
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN ! nothing to output
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose .AND. ionode) WRITE (UNIT=environ_unit, FMT=3100)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=3101) electrolyte%ntyp
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=3102) electrolyte%temperature
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=3103) 1.D0 / SQRT(electrolyte%k2)
            !
            IF (electrolyte%cionmax > 0.D0) THEN
                IF (ionode) WRITE (UNIT=environ_unit, FMT=3104) electrolyte%cionmax
            END IF
            !
            DO ityp = 1, electrolyte%ntyp
                !
                IF (ionode) &
                    WRITE (UNIT=environ_unit, FMT=3105) &
                    electrolyte%ioncctype(ityp)%index
                !
                IF (ionode) &
                    WRITE (UNIT=environ_unit, FMT=3106) &
                    electrolyte%ioncctype(ityp)%cbulk
                !
                IF (ionode) &
                    WRITE (UNIT=environ_unit, FMT=3107) &
                    electrolyte%ioncctype(ityp)%cbulk * amu_si / bohr_radius_si**3
                !
                IF (ionode) &
                    WRITE (UNIT=environ_unit, FMT=3108) electrolyte%ioncctype(ityp)%z
                !
                IF (verbosity >= 5) &
                    CALL print_environ_density(electrolyte%ioncctype(ityp)%c, &
                                               passed_verbosity, passed_depth)
                !
                IF (verbosity >= 5) &
                    CALL print_environ_density(electrolyte%ioncctype(ityp)%cfactor, &
                                               passed_verbosity, passed_depth)
                !
            END DO
            !
            IF (verbosity >= 3) &
                CALL print_environ_density(electrolyte%density, passed_verbosity, &
                                           passed_depth)
            !
            IF (verbosity >= 3) &
                CALL print_environ_density(electrolyte%gamma, passed_verbosity, &
                                           passed_depth)
            !
            IF (verbosity >= 5) &
                CALL print_environ_density(electrolyte%dgamma, passed_verbosity, &
                                           passed_depth)
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=3109) electrolyte%linearized
            !
            IF (ionode) WRITE (UNIT=environ_unit, FMT=3110) electrolyte%charge
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
3100    FORMAT(/, 4('%'), ' ELECTROLYTE ', 64('%'))
3101    FORMAT(1X, 'number electrol. species   = ', I4)
3102    FORMAT(1X, 'solvent temperature        = ', F7.1)
3103    FORMAT(1X, 'Debye length / sqrt(eps)   = ', F14.7)
        !
3104    FORMAT(1X, 'modified Poisson-Boltzmann ' &
               /, 1X, 'maximum concentration      = ', F14.7)
        !
3105    FORMAT(1X, 'electrolyte species ', I4)
3106    FORMAT(1X, 'bulk concentration  (a.u.) = ', E15.4)
3107    FORMAT(1X, 'bulk concentration (mol/L) = ', F14.7)
3108    FORMAT(1X, 'ionic charge               = ', F7.2)
        !
3109    FORMAT(1X, 'electrolyte flags' &
               /, 1X, 'linearized                 = ', L2)
        !
3110    FORMAT(1X, 'total electrolyte charge   = ', F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_electrolyte
    !------------------------------------------------------------------------------------
    !>
    !! Write out the different Environ contributions to the energy.
    !! Called by electrons.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_print_energies()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80) :: sub_name = 'environ_print_energies'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ionode) THEN
            !
            IF (prog == 'PW') THEN
                !
                IF (lelectrostatic) WRITE (program_unit, 9201) eelectrostatic
                !
                IF (lsurface) WRITE (program_unit, 9202) esurface
                !
                IF (lvolume) WRITE (program_unit, 9203) evolume
                !
                IF (lelectrolyte) WRITE (program_unit, 9205) eelectrolyte
                !
                IF (lconfine) WRITE (program_unit, 9206) econfine
                !
                WRITE (program_unit, 9204) deenviron
            ELSE IF (prog == 'CP') THEN
                !
                IF (lelectrostatic) WRITE (program_unit, 9301) eelectrostatic
                !
                IF (lsurface) WRITE (program_unit, 9302) esurface
                !
                IF (lvolume) WRITE (program_unit, 9303) evolume
                !
                IF (lelectrolyte) WRITE (program_unit, 9305) eelectrolyte
                !
                IF (lconfine) WRITE (program_unit, 9306) econfine
                !
                WRITE (program_unit, 9304) deenviron
            ELSE
                CALL env_errore(sub_name, 'Wrong program calling Environ', 1)
            END IF
            !
        END IF
        !
        RETURN
        !
9201    FORMAT('     electrostatic embedding   =', F17.8, ' Ry')
9202    FORMAT('     cavitation energy         =', F17.8, ' Ry')
9203    FORMAT('     PV energy                 =', F17.8, ' Ry')
9206    FORMAT('     confinement energy        =', F17.8, ' Ry')
9205    FORMAT('     electrolyte free energy   =', F17.8, ' Ry')
9204    FORMAT('     correction to one-el term =', F17.8, ' Ry')
9301    FORMAT('     electrostatic embedding = ', F14.5, ' Hartree a.u.')
9302    FORMAT('           cavitation energy = ', F14.5, ' Hartree a.u.')
9303    FORMAT('                   PV energy = ', F14.5, ' Hartree a.u.')
9305    FORMAT('     electrolyte free energy = ', F14.5, ' Hartree a.u.')
9306    FORMAT('          confinement energy = ', F14.5, ' Hartree a.u.')
9304    FORMAT('   correction to one-el term = ', F14.5, ' Hartree a.u.')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_print_energies
    !------------------------------------------------------------------------------------
    !>
    !! If Gaussian nuclei are used, write out the corresponding potential shift
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_print_potential_shift()
        !--------------------------------------------------------------------------------
        !
        IF (lsmearedions) WRITE (program_unit, 9400) potential_shift * RYTOEV
        !
9400    FORMAT(/, 5(' '), &
                'the potential shift due to the Gaussian-smeared nuclei is ', &
                F10.4, ' ev')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_print_potential_shift
    !------------------------------------------------------------------------------------
    !
    !>
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_print_potential_warning()
        !--------------------------------------------------------------------------------
        !
        IF (need_pbc_correction) WRITE (program_unit, 9401)
        !
9401    FORMAT(/, &
                5(' '), 'WARNING: you are using the parabolic pbc correction;', /, &
                5(' '), '         the potential shift above must be added to ', /, &
                5(' '), '         band and Fermi energies.')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_print_potential_warning
    !------------------------------------------------------------------------------------
    !
    !> Write out the main parameters of Environ calculations, summarizing
    !! the input keywords (some info also on internal vs input units).
    !! Called by summary.f90
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_summary()
        !--------------------------------------------------------------------------------
        !
        IF (ionode) THEN
            !
            IF (prog == 'PW') THEN
                WRITE (program_unit, *)
                WRITE (UNIT=program_unit, FMT=9000)
                WRITE (UNIT=program_unit, FMT=8000) bibliography(1)
                !
                !------------------------------------------------------------------------
                ! Environ Summary
                !
                WRITE (UNIT=program_unit, FMT=9001) environ_thr
                !
                IF (lsolvent) THEN
                    !
                    IF (solvent%type_ == 0) THEN
                        WRITE (UNIT=program_unit, FMT=9002) 'Fatteber-Gygi'
                        !
                        WRITE (UNIT=program_unit, FMT=9003) &
                            solvent%rhozero, solvent%tbeta
                        !
                    ELSE IF (solvent%type_ == 1) THEN
                        WRITE (UNIT=program_unit, FMT=9002) 'SCCS'
                        !
                        WRITE (UNIT=program_unit, FMT=9004) &
                            solvent%rhomax, solvent%rhomin
                        !
                    END IF
                    !
                    IF (solvent%solvent_aware) WRITE (UNIT=program_unit, FMT=9013)
                    !
                    IF (solvent%field_aware) THEN
                        WRITE (UNIT=program_unit, FMT=9014)
                        !
                        WRITE (UNIT=program_unit, FMT=9015) &
                            solvent%field_factor, solvent%charge_asymmetry
                        !
                        WRITE (UNIT=program_unit, FMT=9016) &
                            solvent%field_min, solvent%field_max
                        !
                    END IF
                    !
                END IF
                !
                IF (env_static_permittivity > 1.D0) THEN
                    WRITE (UNIT=program_unit, FMT=9005) env_static_permittivity
                    !
                    IF (ltddfpt) &
                        WRITE (UNIT=program_unit, FMT=9006) env_optical_permittivity
                    !
                    WRITE (UNIT=program_unit, FMT=9007) TRIM(solvent%mode)
                END IF
                !
                IF (env_surface_tension > 0.D0) &
                    WRITE (UNIT=program_unit, FMT=9010) &
                    env_surface_tension / 1.D-3 / bohr_radius_si**2 * rydberg_si, &
                    env_surface_tension
                !
                IF (env_pressure /= 0.D0) &
                    WRITE (UNIT=program_unit, FMT=9011) &
                    env_pressure * rydberg_si / bohr_radius_si**3 * 1.D-9, &
                    env_pressure
                !
                !------------------------------------------------------------------------
                ! Electrostatic Summary
                !
                IF (lelectrostatic) THEN
                    !
                    WRITE (UNIT=program_unit, FMT=9100)
                    !
                    WRITE (UNIT=program_unit, FMT=9101) &
                        outer%problem, outer%solver%type_, outer%solver%auxiliary
                    !
                    WRITE (UNIT=program_unit, FMT=9102) outer%core%type_
                    !
                    IF (lfd) THEN
                        !
                        IF (fd%ifdtype == 1) THEN
                            !
                            WRITE (UNIT=program_unit, FMT=9103) &
                                'central diff.', fd%nfdpoint
                            !
                        ELSE IF (fd%ifdtype == 2 .OR. fd%ifdtype == 3) THEN
                            !
                            WRITE (UNIT=program_unit, FMT=9103) &
                                'lanczos diff.', fd%nfdpoint
                            !
                        ELSE IF (fd%ifdtype == 4 .OR. fd%ifdtype == 5) THEN
                            !
                            WRITE (UNIT=program_unit, FMT=9103) &
                                'noise-robust diff.', fd%nfdpoint
                            !
                        END IF
                        !
                    END IF
                    !
                END IF
                !
                WRITE (UNIT=program_unit, FMT=8001)
                !
            END IF
            !
        END IF
        !
8000    FORMAT(/, 5X, 'Plese cite', /, 9X, A80, &
                /, 5X, 'in publications or presentations arising from this work.',/)
        !
8001    FORMAT(/)
9000    FORMAT(/, 5X, 'Environ Module', /, 5X, '==============')
9001    FORMAT('     compensation onset threshold      = ', E24.4, ' ')
9002    FORMAT('     switching function adopted        = ', A24, ' ')
        !
9003    FORMAT('     solvation density threshold       = ', E24.4, ' ' &
               /'     smoothness exponent (2 x beta)    = ', F24.2, ' ')
        !
9004    FORMAT('     density limit for vacuum region   = ', E24.4, ' ' &
               /'     density limit for bulk solvent    = ', E24.4, ' ')
        !
9005    FORMAT('     static permittivity               = ', F24.2, ' ')
9006    FORMAT('     optical permittivity              = ', F24.4, ' ')
9007    FORMAT('     epsilon calculation mode          = ', A24, ' ')
        !
9010    FORMAT('     surface tension in input (dyn/cm) = ', F24.2, ' ' &
               /'     surface tension in internal units = ', E24.4, ' ')
        !
9011    FORMAT('     external pressure in input (GPa)  = ', F24.2, ' ' &
               /'     external pressure in inter. units = ', E24.4, ' ')
        !
9012    FORMAT('     correction slab geom. along axis  = ', I24, ' ')
9013    FORMAT('     interface is solvent aware            ')
9014    FORMAT('     interface is field aware            ')
        !
9015    FORMAT('     field aware factor                = ', F24.2, ' ' &
               /'     asymmetry of field-awareness      = ', F24.2, ' ')
        !
9016    FORMAT('     field limit for no correction     = ', F24.2, ' ' &
               /'     field limit for full correction   = ', F24.2, ' ')
        !
9100    FORMAT(/, 5X, 'Electrostatic Setup', /, 5X, '-------------------')
        !
9101    FORMAT('     electrostatic problem to solve    = ', A24, ' ' &
               /'     numerical solver adopted          = ', A24, ' ' &
               /'     type of auxiliary density adopted = ', A24, ' ')
        !
9102    FORMAT('     type of core tool for poisson     = ', A24, ' ')
        !
9103    FORMAT('     type of numerical differentiator  = ', A24, ' ' &
               /'     number of points in num. diff.    = ', I24, ' ')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_summary
    !------------------------------------------------------------------------------------
    !>
    !! Write out the time informations of the Environ dependent calculations.
    !! Called by print_clock_pw.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_clock(passed_unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN), OPTIONAL :: passed_unit
        !
        INTEGER :: actual_unit
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(passed_unit)) THEN
            actual_unit = passed_unit
        ELSE
            actual_unit = program_unit
        END IF
        !
        WRITE (actual_unit, *)
        WRITE (actual_unit, '(5X,"Environ routines")')
        !
        !--------------------------------------------------------------------------------
        ! Dielectric subroutines
        !
        IF (lelectrostatic) THEN
            !
            CALL env_print_clock('calc_eelect')
            !
            CALL env_print_clock('calc_velect')
            !
            CALL env_print_clock('calc_vgcs')
            !
            CALL env_print_clock('dielectric')
            !
            CALL env_print_clock('electrolyte')
            !
            CALL env_print_clock('calc_felect')
            !
        END IF
        !
        IF (lsemiconductor) CALL env_print_clock('calc_vms')
        !
        !--------------------------------------------------------------------------------
        ! TDDFT
        !
        IF (ltddfpt) CALL env_print_clock('calc_vsolvent_tddfpt')
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_clock
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE write_cube(f, ions, idx, label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), TARGET, INTENT(IN) :: f
        TYPE(environ_ions), TARGET, OPTIONAL, INTENT(IN) :: ions
        INTEGER, INTENT(IN), OPTIONAL :: idx
        CHARACTER(LEN=100), INTENT(IN), OPTIONAL :: label
        !
        INTEGER :: ir, ir1, ir2, ir3, num
        INTEGER :: ipol, iat, typ, count
        INTEGER :: nr1x, nr2x, nr3x
        INTEGER :: nr1, nr2, nr3
        !
        REAL(DP) :: tmp, scale
        REAL(DP), ALLOCATABLE :: flocal(:)
        !
        CHARACTER(LEN=100) :: filename
        REAL(DP), POINTER :: alat
        REAL(DP), POINTER :: at(:, :), origin(:)
        !
        INTEGER :: nat
        INTEGER, POINTER :: ityp(:)
        REAL(DP), POINTER :: tau(:, :)
        !
        CHARACTER(LEN=100) :: filemod
        !
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        CHARACTER(LEN=80) :: sub_name = 'write_cube'
        !
        !--------------------------------------------------------------------------------
        !
        dfft => f%cell%dfft
        !
        nr1x = dfft%nr1x
        nr2x = dfft%nr2x
        nr3x = dfft%nr3x
        !
        nr1 = dfft%nr1
        nr2 = dfft%nr2
        nr3 = dfft%nr3
        !
        IF (PRESENT(idx)) THEN
            WRITE (filemod, '(i4.4)') idx
        ELSE
            filemod = ""
        END IF
        !
        IF (PRESENT(label)) THEN
            filename = TRIM(ADJUSTL(label))//TRIM(filemod)//".cube"
        ELSE
            filename = TRIM(ADJUSTL(f%label))//TRIM(filemod)//".cube"
        END IF
        !
        alat => f%cell%alat
        at => f%cell%at
        origin => f%cell%origin
        !
        IF (PRESENT(ions)) THEN
            nat = ions%number
            ityp => ions%ityp
            tau => ions%tau
        ELSE
            nat = 1
        END IF
        !
        ALLOCATE (flocal(nr1x * nr2x * nr3x))
#if defined(__MPI)
        flocal = 0.D0
        !
        CALL env_gather_grid(dfft, f%of_r, flocal)
        !
        CALL env_mp_sum(flocal, dfft%comm)
        !
#else
        flocal = f%of_r
#endif
        IF (ionode) THEN
            !
            OPEN (300, file=TRIM(filename), status='unknown')
            !
            scale = alat !*0.52917720859d0
            WRITE (300, *) 'CUBE FILE GENERATED BY PW.X'
            WRITE (300, *) 'OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z'
            WRITE (300, '(i5,3f12.6)') nat, origin * scale
            WRITE (300, '(i5,3f12.6)') nr1, (at(ipol, 1) / DBLE(nr1) * scale, ipol=1, 3)
            WRITE (300, '(i5,3f12.6)') nr2, (at(ipol, 2) / DBLE(nr2) * scale, ipol=1, 3)
            WRITE (300, '(i5,3f12.6)') nr3, (at(ipol, 3) / DBLE(nr3) * scale, ipol=1, 3)
            !
            IF (PRESENT(ions)) THEN
                !
                DO iat = 1, nat
                    typ = ityp(iat)
                    num = ions%iontype(typ)%atmnum
                    !
                    WRITE (300, '(i5,4f12.6)') &
                        num, 0.D0, tau(1, iat) * scale, &
                        tau(2, iat) * scale, tau(3, iat) * scale
                    !
                END DO
                !
            ELSE
                WRITE (300, '(i5,4f12.6)') 1, 0.D0, 0.D0, 0.D0, 0.D0
            END IF
            !
            count = 0
            !
            DO ir1 = 1, nr1
                !
                DO ir2 = 1, nr2
                    !
                    DO ir3 = 1, nr3
                        count = count + 1
                        ir = ir1 + (ir2 - 1) * nr1 + (ir3 - 1) * nr1 * nr2
                        tmp = DBLE(flocal(ir))
                        !
                        IF (ABS(tmp) < 1.D-99) tmp = 0.D0
                        !
                        IF (MOD(count, 6) == 0) THEN
                            WRITE (300, '(e13.6,1x)') tmp
                        ELSE
                            WRITE (300, '(e13.6,1x)', advance='no') tmp
                        END IF
                        !
                    END DO
                    !
                END DO
                !
            END DO
            !
            CLOSE (300)
            !
        END IF
        !
        DEALLOCATE (flocal)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE write_cube
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE environ_output
!----------------------------------------------------------------------------------------
