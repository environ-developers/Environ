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
!          Matthew Truscott   (Department of Physics, UNT)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE environ_debugging
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
    USE class_function_gaussian
    USE class_functions
    USE class_gradient
    !
    USE class_boundary
    USE class_electrons
    USE class_ions
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: test_de_dboundary
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! (DEBUGGING) Test functional derivative of energy w.r.t local boundary
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE test_de_dboundary(boundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_boundary), INTENT(IN), TARGET :: boundary
        !
        TYPE(environ_boundary) :: localbound
        TYPE(environ_density) :: de_dboundary
        TYPE(environ_function_gaussian) :: test_function
        TYPE(environ_density) :: delta
        TYPE(environ_gradient) :: graddelta
        !
        TYPE(environ_cell), POINTER :: cell
        !
        INTEGER :: i
        REAL(DP) :: localpressure, localsurface_tension
        REAL(DP) :: eplus, eminus, de_fd, de_analytic, epsilon
        !
        !--------------------------------------------------------------------------------
        !
        cell => boundary%scaled%cell
        !
        CALL boundary%copy(localbound)
        !
        CALL de_dboundary%init(cell)
        !
        CALL localbound%solvent_aware_boundary()
        !
        localpressure = 100.D0
        !
        CALL localbound%devolume_dboundary(localpressure, de_dboundary)
        !
        localsurface_tension = 100.D0
        !
        CALL localbound%desurface_dboundary(localsurface_tension, de_dboundary)
        !
        CALL localbound%sa_de_dboundary(de_dboundary)
        !
        test_function%f_type = 1
        test_function%dim = 0
        test_function%axis = 3
        test_function%spread = 0.3D0
        test_function%width = 0.D0
        test_function%volume = 1.D0
        !
        epsilon = 0.000001
        !
        CALL delta%init(cell)
        !
        CALL graddelta%init(cell)
        !
        ALLOCATE (test_function%pos(3))
        test_function%pos(1) = 11.79D0
        test_function%pos(2) = 12.05D0
        !
        DO i = 1, cell%dfft%nr3
            !
            test_function%pos(3) = DBLE(i - 1) * cell%at(3, 3) / DBLE(cell%dfft%nr3)
            !
            CALL test_function%density(delta, .TRUE.)
            !
            CALL test_function%gradient(graddelta, .TRUE.)
            !
            de_fd = 0.D0
            !
            CALL boundary%copy(localbound)
            !
            localbound%scaled%of_r = localbound%scaled%of_r + epsilon * delta%of_r
            localbound%volume = localbound%scaled%integrate()
            !
            IF (localbound%deriv .GE. 1) THEN
                !
                localbound%gradient%of_r = localbound%gradient%of_r + &
                                           epsilon * graddelta%of_r
                !
                CALL localbound%gradient%update_modulus()
                !
                localbound%surface = localbound%gradient%modulus%integrate()
            END IF
            !
            CALL localbound%solvent_aware_boundary()
            !
            CALL localbound%evolume(localpressure, eplus)
            !
            de_fd = de_fd + eplus
            !
            CALL localbound%esurface(localsurface_tension, eplus)
            !
            de_fd = de_fd + eplus
            !
            CALL boundary%copy(localbound)
            !
            localbound%scaled%of_r = localbound%scaled%of_r - epsilon * delta%of_r
            localbound%volume = localbound%scaled%integrate()
            !
            IF (localbound%deriv .GE. 1) THEN
                !
                localbound%gradient%of_r = localbound%gradient%of_r - &
                                           epsilon * graddelta%of_r
                !
                CALL localbound%gradient%update_modulus()
                !
                localbound%surface = localbound%gradient%modulus%integrate()
            END IF
            !
            CALL localbound%solvent_aware_boundary()
            !
            CALL localbound%evolume(localpressure, eminus)
            !
            de_fd = de_fd - eminus
            !
            CALL localbound%esurface(localsurface_tension, eminus)
            !
            de_fd = de_fd - eminus
            !
            de_fd = 0.5D0 * de_fd / epsilon
            !
            de_analytic = de_dboundary%scalar_product(delta)
            !
            IF (io%lnode) &
                WRITE (io%debug_unit, '(1X,a,f20.10,3f20.10)') &
                ' z = ', test_function%pos(3), de_analytic, de_fd, de_analytic - de_fd
            !
            FLUSH (io%debug_unit)
            !
        END DO
        !
        CALL delta%destroy()
        !
        CALL graddelta%destroy()
        !
        CALL de_dboundary%destroy()
        !
        CALL localbound%destroy()
        !
        STOP
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE test_de_dboundary
    !------------------------------------------------------------------------------------
    !>
    !! (DEBUGGING)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE test_ion_field_derivatives(ideriv, bound)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ideriv ! .EQ. 1/2 for electronic/ionic
        !
        TYPE(environ_boundary), INTENT(INOUT) :: bound
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density) :: rho
        TYPE(environ_gradient) :: field
        !
        INTEGER :: i, ipol
        REAL(DP) :: dx, x0, epsilon
        REAL(DP), ALLOCATABLE :: fd_partial_of_ion_field(:)
        REAL(DP), ALLOCATABLE :: analytic_partial_of_ion_field(:, :, :)
        !
        INTEGER :: j
        REAL(DP) :: tmp
        REAL(DP), ALLOCATABLE :: fd_dion_field_drho(:)
        TYPE(environ_density), ALLOCATABLE :: analytic_dion_field_drho(:)
        TYPE(environ_function_gaussian) :: test_function
        TYPE(environ_density) :: delta
        TYPE(environ_electrons) :: localelectrons
        !
        !--------------------------------------------------------------------------------
        !
        cell => bound%scaled%cell
        !
        !--------------------------------------------------------------------------------
        ! Recompute total charge density and field
        !
        CALL rho%init(cell)
        !
        rho%of_r = bound%electrons%density%of_r + bound%ions%density%of_r
        !
        CALL field%init(cell)
        !
        ! CALL gradv_h_of_rho_r(rho%of_r, field%of_r)
        !
        !--------------------------------------------------------------------------------
        ! Print out individual and global fluxes
        !
        DO i = 1, bound%ions%number
            !
            WRITE (io%debug_unit, '(a,i3,a,f14.7)') &
                'flux through soft-sphere number ', i, &
                ' equal to ', bound%ion_field(i)
            !
        END DO
        !
        rho%of_r = rho%of_r * bound%scaled%of_r
        !
        WRITE (io%debug_unit, '(a,f14.7,a,f14.7)') &
            'total_charge = ', rho%integrate(), &
            ' total flux throught !soft-spheres = ', SUM(bound%ion_field(:))
        !
        CALL bound%gradient%scalar_product(field, rho)

        WRITE (io%debug_unit, '(a,f14.7)') 'actual total flux = ', rho%integrate()
        !
        IF (ideriv .EQ. 1) THEN
            !
            !----------------------------------------------------------------------------
            ! Test functional derivative wrt electronic density of each flux
            !
            ALLOCATE (analytic_dion_field_drho(bound%ions%number))
            !
            DO i = 1, bound%ions%number
                CALL bound%dion_field_drho(i)%copy(analytic_dion_field_drho(i))
            END DO
            !
            ALLOCATE (fd_dion_field_drho(bound%ions%number))
            !
            test_function%f_type = 1
            test_function%dim = 0
            test_function%axis = 3
            test_function%spread = 0.3D0
            test_function%width = 0.D0
            test_function%volume = 1.D0
            !
            epsilon = 0.000001
            !
            CALL delta%init(cell)
            !
            CALL localelectrons%density%init(cell)
            !
            !----------------------------------------------------------------------------
            ! We are only going to check delta functions along the z axis
            ! passing throught the O atom
            !
            ALLOCATE (test_function%pos(3))
            test_function%pos(1) = 11.79D0
            test_function%pos(2) = 12.05D0
            !
            DO i = 1, cell%dfft%nr3
                !
                test_function%pos(3) = DBLE(i - 1) * cell%at(3, 3) / DBLE(cell%dfft%nr3)
                !
                CALL test_function%density(delta, .TRUE.)
                !
                fd_dion_field_drho = 0.D0
                !
                localelectrons%density%of_r = bound%electrons%density%of_r - &
                                              epsilon * delta%of_r
                !
                ! CALL compute_ion_field(bound%ions%number, bound%local_spheres, &
                !                        bound%ions, localelectrons, bound%ion_field)
                !
                fd_dion_field_drho = bound%ion_field
                !
                localelectrons%density%of_r = bound%electrons%density%of_r + &
                                              epsilon * delta%of_r
                !
                ! CALL compute_ion_field(bound%ions%number, bound%local_spheres, &
                !                        bound%ions, localelectrons, bound%ion_field)
                !
                fd_dion_field_drho = bound%ion_field - fd_dion_field_drho
                fd_dion_field_drho = fd_dion_field_drho * 0.5D0 / epsilon
                !
                IF (io%lnode) &
                    WRITE (io%debug_unit, *) ' z = ', test_function%pos(3)
                !
                DO j = 1, bound%ions%number
                    tmp = analytic_dion_field_drho(j)%scalar_product(delta)
                    !
                    IF (io%lnode) &
                        WRITE (io%debug_unit, '(a,i3,3f20.10)') 'ion = ', j, &
                        tmp, fd_dion_field_drho(j), tmp - fd_dion_field_drho(j)
                    !
                END DO
                !
            END DO
            !
            DEALLOCATE (test_function%pos)
            !
            CALL delta%destroy()
            !
            DO i = 1, bound%ions%number
                CALL analytic_dion_field_drho(i)%destroy()
            END DO
            !
            DEALLOCATE (analytic_dion_field_drho)
            DEALLOCATE (fd_dion_field_drho)
            !
        ELSE IF (ideriv .EQ. 2) THEN
            !
            !----------------------------------------------------------------------------
            ! Test derivative wrt atomic positions with finite differences
            !
            ! CALL compute_ion_field_partial(bound%ions%number, bound%local_spheres, &
            !                                bound%ions, bound%electrons, &
            !                                bound%ion_field, bound%partial_of_ion_field, &
            !                                bound%core%fft)
            !
            ALLOCATE (analytic_partial_of_ion_field(3, bound%ions%number, &
                                                    bound%ions%number))
            !
            analytic_partial_of_ion_field = bound%partial_of_ion_field
            !
            dx = 2.0D-3 ! dx expected units: BOHR
            !
            ALLOCATE (fd_partial_of_ion_field(bound%ions%number))
            !
            DO i = 1, bound%ions%number
                !
                DO ipol = 1, 3
                    x0 = bound%ions%tau(ipol, i)
                    bound%ions%tau(ipol, i) = x0 - dx
                    !
                    CALL bound%ions%update(bound%ions%number, bound%ions%tau)
                    !
                    rho%of_r = bound%electrons%density%of_r + bound%ions%density%of_r
                    !
                    ! CALL gradv_h_of_rho_r(rho%of_r, field%of_r)
                    !
                    ! CALL compute_ion_field(bound%ions%number, bound%local_spheres, &
                    !                        bound%ions, bound%electrons, &
                    !                        bound%ion_field)
                    !
                    fd_partial_of_ion_field = bound%ion_field
                    !
                    bound%ions%tau(ipol, i) = x0 + dx
                    !
                    CALL bound%ions%update(bound%ions%number, bound%ions%tau)
                    !
                    rho%of_r = bound%electrons%density%of_r + bound%ions%density%of_r
                    !
                    ! CALL gradv_h_of_rho_r(rho%of_r, field%of_r)
                    !
                    ! CALL compute_ion_field(bound%ions%number, bound%local_spheres, &
                    !                        bound%ions, bound%electrons, &
                    !                        bound%ion_field)
                    !
                    fd_partial_of_ion_field = bound%ion_field - fd_partial_of_ion_field
                    !
                    fd_partial_of_ion_field = fd_partial_of_ion_field / 2.D0 / dx
                    !
                    WRITE (io%debug_unit, *) ' i  = ', i, ' ipol = ', ipol
                    !
                    WRITE (io%debug_unit, '(a,10f20.10)') 'analytic     = ', &
                        analytic_partial_of_ion_field(ipol, :, i)
                    !
                    WRITE (io%debug_unit, '(a,10f20.10)') 'finite-diff  = ', &
                        fd_partial_of_ion_field(:)
                    !
                    WRITE (io%debug_unit, *) ' '
                    !
                    bound%ions%tau(ipol, i) = x0
                END DO
                !
            END DO
            !
            DEALLOCATE (analytic_partial_of_ion_field)
            DEALLOCATE (fd_partial_of_ion_field)
        END IF
        !
        STOP
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE test_ion_field_derivatives
    !------------------------------------------------------------------------------------
    !>
    !! (DEBUGGING)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_test_boundary(bound, electrons)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_electrons), INTENT(IN) :: electrons
        !
        TYPE(environ_boundary), INTENT(INOUT) :: bound
        !
        INTEGER :: debugcubes = 0
        !
        CHARACTER(len=80) :: sub_name = 'update_test_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT CASE (bound%mode)
            !
        CASE ('fa-electronic')
            !
            ! CALL compute_normal_field(bound%ions, electrons, bound%normal_field)
            !
            ! CALL field_aware_density(electrons, bound)
            !
            CALL bound%boundary_of_density()
            !
            SELECT CASE (debugcubes)
                !
            CASE (2)
                CALL bound%scaled%write_cube_no_ions(label='standard')
                !
            CASE (1)
                CALL bound%density%write_cube_no_ions(label='standard')
                !
            CASE DEFAULT
                !
            END SELECT
            !
        CASE ('fa-ionic')
            !
            ! CALL compute_ion_field(bound%ions%number, bound%local_spheres, &
            !                        bound%ions, electrons, bound%ion_field)
            !
            ! CALL bound%set_soft_spheres()
            !
            CALL bound%boundary_of_functions()
            !
        CASE ('electronic')
            CALL bound%boundary_of_density(electrons%density)
            !
        CASE DEFAULT
            CALL io%error(sub_name, 'Unrecognized boundary mode', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_test_boundary
    !------------------------------------------------------------------------------------
    !>
    !! (DEBUGGING)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_ionic(ions)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_ions), INTENT(IN) :: ions
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        DO i = 1, ions%number
            PRINT *, "ATOM", i, ions%tau(1, i), ions%tau(2, i), ions%tau(3, i)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_ionic
    !------------------------------------------------------------------------------------
    !>
    !! (DEBUGGING)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE test_energy_derivatives(ideriv, bound)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ideriv
        !
        TYPE(environ_boundary), INTENT(INOUT) :: bound
        !
        TYPE(environ_cell), POINTER :: cell
        !
        INTEGER :: i, ipol
        REAL(DP) :: epsilon
        REAL(DP) :: localpressure, localsurface_tension
        REAL(DP) :: de_fd, de_analytic, etmp
        REAL(DP) :: dx, x0, force(3), ssforce(3)
        REAL(DP) :: flux
        REAL(DP), ALLOCATABLE :: tau0(:, :)
        !
        TYPE(environ_density) :: de_dboundary
        TYPE(environ_density) :: vanalytic
        TYPE(environ_function_gaussian) :: test_function
        TYPE(environ_density) :: delta
        TYPE(environ_electrons) :: localelectrons
        TYPE(environ_gradient) :: partial
        !
        CHARACTER(LEN=80) :: sub_name = 'test_energy_derivatives'
        !
        !--------------------------------------------------------------------------------
        !
        cell => bound%scaled%cell
        !
        !--------------------------------------------------------------------------------
        ! Compute the field and the field-aware interface (Need to check)
        !
        CALL update_test_boundary(bound, bound%electrons)
        !
        !--------------------------------------------------------------------------------
        ! Compute functional derivative of the energy wrt interface
        !
        CALL de_dboundary%init(cell)
        !
        localpressure = -0.35
        !
        CALL bound%devolume_dboundary(localpressure, de_dboundary)
        !
        localsurface_tension = 0.0
        !
        ! CALL bound%desurface_dboundary(localsurface_tension, de_dboundary)
        !
        IF (ideriv .EQ. 1) THEN
            !
            !----------------------------------------------------------------------------
            ! Compute functional derivative wrt electronic density
            !
            ! IF (bound%mode .EQ. 'fa-ionic') THEN
            !     CALL compute_dion_field_drho(bound%ions%number, bound%local_spheres, &
            !                                  bound%dion_field_drho, bound%core%fft)
            ! END IF
            !
            CALL vanalytic%init(cell)
            !
            IF (bound%field_aware) THEN
                !
                ! CALL field_aware_de_drho(bound, de_dboundary, vanalytic)
                !
                CALL io%error('field-aware6', 'Option not yet implimented ', 1)
                !
            ELSE
                vanalytic%of_r = bound%dscaled%of_r * de_dboundary%of_r
            END IF
            !
            !----------------------------------------------------------------------------
            ! Loop over gridpoints with rhoelec + or - a delta function
            !
            test_function%f_type = 1
            test_function%dim = 0
            test_function%axis = 3
            test_function%spread = 0.4D0
            test_function%width = 0.D0
            test_function%volume = 1.D0
            !
            epsilon = 0.000008
            !
            CALL delta%init(cell)
            !
            CALL localelectrons%density%init(cell)
            !
            !----------------------------------------------------------------------------
            ! We are only going to check delta functions along the z-axis
            ! passing throught !the O atom
            !
            ALLOCATE (test_function%pos(3))
            test_function%pos(1) = 0.0
            test_function%pos(2) = 0.0
            !
            DO i = 1, cell%dfft%nr3
                de_fd = 0.D0
                test_function%pos(3) = DBLE(i - 1) * cell%at(3, 3) / DBLE(cell%dfft%nr3)
                !
                IF (io%lnode) &
                    WRITE (io%unit, '(a,f14.7)') ' z = ', test_function%pos(3)
                !
                ! IF (test_function%pos(3) .LE. 0.365) CYCLE
                !
                ! IF (test_function%pos(3) .LE. 4.5 .OR. &
                !     test_function%pos(3) .GE. 7.0) CYCLE
                !
                CALL test_function%density(delta, .TRUE.)
                !
                localelectrons%density%of_r = bound%electrons%density%of_r + &
                                              epsilon * delta%of_r
                !
                CALL update_test_boundary(bound, localelectrons)
                !
                CALL bound%evolume(localpressure, etmp)
                de_fd = de_fd + etmp
                !
                CALL bound%esurface(localsurface_tension, etmp)
                de_fd = de_fd + etmp
                !
                ! IF (io%lnode) &
                !     WRITE (io%debug_unit, '(2f20.10)') bound%surface, bound%volume
                !
                localelectrons%density%of_r = bound%electrons%density%of_r - &
                                              epsilon * delta%of_r
                !
                CALL update_test_boundary(bound, localelectrons)
                !
                CALL bound%evolume(localpressure, etmp)
                de_fd = de_fd - etmp
                !
                CALL bound%esurface(localsurface_tension, etmp)
                de_fd = de_fd - etmp
                !
                ! IF (io%lnode) &
                !     WRITE (io%debug_unit, '(2f20.10)') bound%surface, bound%volume
                !
                de_fd = de_fd * 0.5D0 / epsilon
                !
                de_analytic = vanalytic%scalar_product(delta)
                !
                IF (io%lnode) &
                    WRITE (io%debug_unit, '(a,i14.7,3f20.10)') ' z = ', i, &
                    de_analytic, de_fd, de_analytic - de_fd
                !
                FLUSH (io%debug_unit)
                !
                ! STOP
                !
            END DO
            !
            DEALLOCATE (test_function%pos)
            !
            CALL delta%destroy()
            !
            CALL vanalytic%destroy()
            !
        ELSE IF (ideriv .EQ. 2) THEN
            PRINT *, 'ideriv = 2'
            !
            !----------------------------------------------------------------------------
            ! Compute partial derivative with respect to ionic positions
            !
            ! CALCULATE ion field partials in advance
            !
            ! IF (bound%mode .EQ. 'fa-ionic') THEN
            !     !
            !     IF (io%lnode) &
            !         WRITE (io%unit, '(1X,a)') 'outside compute_ion_field_partial'
            !     !
            !     CALL compute_ion_field_partial(bound%ions%number, &
            !                                    bound%local_spheres, &
            !                                    bound%ions, bound%electrons, &
            !                                    bound%ion_field, &
            !                                    bound%partial_of_ion_field, &
            !                                    bound%core%fft)
            !     !
            ! END IF
            !
            dx = 2.0D-3 ! dx expected units: BOHR
            !
            CALL partial%init(cell)
            !
            ALLOCATE (tau0(3, bound%ions%number))
            !
            tau0(:, :) = bound%ions%tau(:, :)
            !
            DO i = 1, bound%ions%number
                !
                bound%ions%tau(:, :) = tau0(:, :)
                ! CALCULATE FORCE FIRST, reset the ionic positions
                !
                CALL bound%ions%update(bound%ions%number, bound%ions%tau)
                !
                CALL update_test_boundary(bound, bound%electrons)
                !
                IF (bound%mode .EQ. 'fa-ionic') THEN
                    CALL bound%dboundary_dions(i, partial)
                    !
                    ssforce = -partial%scalar_product_density(de_dboundary)
                END IF
                !
                ! IF (bound%field_aware) &
                !     CALL field_aware_dboundary_dions(i, bound, partial)
                !
                force = -partial%scalar_product_density(de_dboundary)
                !
                DO ipol = 1, 3
                    !
                    de_fd = 0.D0 ! FINITE DIFFERENCE VALUE STORED HERE
                    bound%ions%tau(:, :) = tau0(:, :) ! RESET IONS
                    !
                    bound%ions%tau(ipol, i) = tau0(ipol, i) - dx ! MINUS dx
                    !
                    CALL bound%ions%update(bound%ions%number, bound%ions%tau)
                    !
                    CALL update_test_boundary(bound, bound%electrons)
                    !
                    CALL bound%evolume(localpressure, etmp)
                    !
                    de_fd = de_fd - etmp
                    !
                    CALL bound%esurface(localsurface_tension, etmp)
                    !
                    de_fd = de_fd - etmp
                    !
                    bound%ions%tau(ipol, i) = tau0(ipol, i) + dx ! PLUS dx
                    !
                    CALL bound%ions%update(bound%ions%number, bound%ions%tau)
                    !
                    CALL update_test_boundary(bound, bound%electrons)
                    !
                    CALL bound%evolume(localpressure, etmp)
                    !
                    de_fd = de_fd + etmp
                    !
                    CALL bound%esurface(localsurface_tension, etmp)
                    !
                    de_fd = de_fd + etmp
                    !
                    de_fd = de_fd * (-0.5D0) / dx
                    ! force is negative of the energy derivative
                    !
                    IF (bound%mode .EQ. 'fa-ionic') THEN
                        !
                        IF (io%lnode) &
                            WRITE (io%debug_unit, '(a,i3,a,i3,4f20.10)') ' i = ', i, &
                            ' ipol != ', ipol, force(ipol), ssforce(ipol), de_fd, &
                            force(ipol) - de_fd
                        !
                    ELSE
                        !
                        IF (io%lnode) &
                            WRITE (io%debug_unit, '(a,i3,a,i3,3f20.10)') ' i = ', i, &
                            ' ipol != ', ipol, force(ipol), de_fd, force(ipol) - de_fd
                        !
                    END IF
                    !
                    ! STOP
                    !
                END DO
                !
            END DO
            !
            CALL partial%destroy()
            !
        END IF
        !
        CALL de_dboundary%destroy()
        !
        STOP
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE test_energy_derivatives
    !------------------------------------------------------------------------------------
    !>
    !! (DEBUGGING)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE extract_boundary_data(bound)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_boundary), INTENT(IN) :: bound
        !
        TYPE(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'extract_boundary_data'
        !
        !--------------------------------------------------------------------------------
        !
        IF (io%lnode) WRITE (io%unit, '(1X,a)') 'extract_boundary_data'
        !
        !--------------------------------------------------------------------------------
        ! Compute the field and the field-aware interface
        !
        CALL bound%electrons%density%write_cube_no_ions(label='e')
        !
        CALL bound%ions%density%write_cube_no_ions(label='i')
        !
        STOP
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE extract_boundary_data
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE environ_debugging
!----------------------------------------------------------------------------------------
