!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
MODULE environ_debugging
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP
    !
    USE physical_types, ONLY: environ_boundary
    USE representation_types, ONLY: environ_density, environ_functions, environ_gradient
    USE cell_types, ONLY: environ_cell
    !
    USE utils_boundary, ONLY: copy_environ_boundary, destroy_environ_boundary
    USE utils_density, ONLY: init_environ_density, destroy_environ_density
    !
    USE utils_gradient, ONLY: init_environ_gradient, update_gradient_modulus, &
                              destroy_environ_gradient
    !
    USE tools_math, ONLY: scalar_product_environ_density, integrate_environ_density
    USE tools_functions, ONLY: density_of_functions, gradient_of_functions
    !
    USE tools_generate_boundary, ONLY: solvent_aware_boundary, solvent_aware_de_dboundary
    !
    USE embedding_volume
    USE embedding_surface
    !
    USE environ_output, ONLY: ionode, environ_unit
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
        TYPE(environ_functions) :: test_function
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
        CALL copy_environ_boundary(boundary, localbound)
        !
        CALL init_environ_density(cell, de_dboundary)
        !
        CALL solvent_aware_boundary(localbound)
        !
        localpressure = 100.D0
        !
        CALL calc_devolume_dboundary(localpressure, localbound, de_dboundary)
        !
        localsurface_tension = 100.D0
        !
        CALL calc_desurface_dboundary(localsurface_tension, localbound, de_dboundary)
        !
        CALL solvent_aware_de_dboundary(localbound, de_dboundary)
        !
        test_function%type_ = 1
        test_function%dim = 0
        test_function%axis = 3
        test_function%spread = 0.3D0
        test_function%width = 0.D0
        test_function%volume = 1.D0
        !
        epsilon = 0.000001
        !
        CALL init_environ_density(cell, delta)
        !
        CALL init_environ_gradient(cell, graddelta)
        !
        ALLOCATE (test_function%pos(3))
        test_function%pos(1) = 11.79D0 / cell%alat
        test_function%pos(2) = 12.05D0 / cell%alat
        !
        DO i = 1, cell%dfft%nr3
            !
            test_function%pos(3) = DBLE(i - 1) * cell%at(3, 3) / DBLE(cell%dfft%nr3)
            !
            CALL density_of_functions(test_function, delta, .TRUE.)
            !
            CALL gradient_of_functions(test_function, graddelta, .TRUE.)
            !
            de_fd = 0.D0
            !
            CALL copy_environ_boundary(boundary, localbound)
            !
            localbound%scaled%of_r = localbound%scaled%of_r + epsilon * delta%of_r
            localbound%volume = integrate_environ_density(localbound%scaled)
            !
            IF (localbound%deriv .GE. 1) THEN
                !
                localbound%gradient%of_r = localbound%gradient%of_r + &
                                           epsilon * graddelta%of_r
                !
                CALL update_gradient_modulus(localbound%gradient)
                !
                localbound%surface = &
                    integrate_environ_density(localbound%gradient%modulus)
                !
            END IF
            !
            CALL solvent_aware_boundary(localbound)
            !
            CALL calc_evolume(localpressure, localbound, eplus)
            !
            de_fd = de_fd + eplus
            !
            CALL calc_esurface(localsurface_tension, localbound, eplus)
            !
            de_fd = de_fd + eplus
            !
            CALL copy_environ_boundary(boundary, localbound)
            !
            localbound%scaled%of_r = localbound%scaled%of_r - epsilon * delta%of_r
            localbound%volume = integrate_environ_density(localbound%scaled)
            !
            IF (localbound%deriv .GE. 1) THEN
                !
                localbound%gradient%of_r = localbound%gradient%of_r - &
                                           epsilon * graddelta%of_r
                !
                CALL update_gradient_modulus(localbound%gradient)
                !
                localbound%surface = &
                    integrate_environ_density(localbound%gradient%modulus)
                !
            END IF
            !
            CALL solvent_aware_boundary(localbound)
            !
            CALL calc_evolume(localpressure, localbound, eminus)
            !
            de_fd = de_fd - eminus
            !
            CALL calc_esurface(localsurface_tension, localbound, eminus)
            !
            de_fd = de_fd - eminus
            !
            de_fd = 0.5D0 * de_fd / epsilon
            !
            de_analytic = scalar_product_environ_density(de_dboundary, delta)
            !
            IF (ionode) &
                WRITE (environ_unit, '(1X,a,f20.10,3f20.10)') &
                ' z = ', test_function%pos(3) * cell%alat, de_analytic, de_fd, &
                de_analytic - de_fd
            !
            FLUSH (environ_unit)
            !
        END DO
        !
        CALL destroy_environ_density(delta)
        !
        CALL destroy_environ_gradient(graddelta)
        !
        CALL destroy_environ_density(de_dboundary)
        !
        CALL destroy_environ_boundary(.TRUE., localbound)
        !
        STOP
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE test_de_dboundary
    !------------------------------------------------------------------------------------
    ! !>
    ! !! (DEBUGGING)
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE test_ion_field_derivatives(ideriv, bound)
    !     !--------------------------------------------------------------------------------
    !     !
    !     USE utils_ions, ONLY: update_environ_ions
    !     USE tools_generate_boundary, ONLY: compute_ion_field, compute_ion_field_partial
    !     !
    !     IMPLICIT NONE
    !     !
    !     INTEGER, INTENT(IN) :: ideriv ! .EQ. 1/2 for electronic/ionic
    !     TYPE(environ_boundary), INTENT(INOUT) :: bound
    !     !
    !     TYPE(environ_cell), POINTER :: cell
    !     TYPE(environ_density) :: rho
    !     TYPE(environ_gradient) :: field
    !     !
    !     INTEGER :: i, ipol
    !     REAL(DP) :: dx, x0, epsilon
    !     REAL(DP), ALLOCATABLE :: fd_partial_of_ion_field(:)
    !     REAL(DP), ALLOCATABLE :: analytic_partial_of_ion_field(:, :, :)
    !     !
    !     INTEGER :: j
    !     REAL(DP) :: tmp
    !     REAL(DP), ALLOCATABLE :: fd_dion_field_drho(:)
    !     TYPE(environ_density), ALLOCATABLE :: analytic_dion_field_drho(:)
    !     TYPE(environ_functions) :: test_function
    !     TYPE(environ_density) :: delta
    !     TYPE(environ_electrons) :: localelectrons
    !     !
    !     cell => bound%scaled%cell
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! Recompute total charge density and field
    !     !
    !     CALL init_environ_density(cell, rho)
    !     rho%of_r = bound%electrons%density%of_r + bound%ions%density%of_r
    !     !
    !     CALL init_environ_gradient(cell, field)
    !     CALL gradv_h_of_rho_r(rho%of_r, field%of_r)
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! Print out individual and global fluxes
    !     !
    !     DO i = 1, bound%ions%number
    !         !
    !         WRITE (environ_unit, '(a,i3,a,f14.7)') &
    !             'flux through soft-sphere number ', i, &
    !             ' equal to ', bound%ion_field(i)
    !         !
    !     END DO
    !     !
    !     rho%of_r = rho%of_r * bound%scaled%of_r
    !     !
    !     WRITE (environ_unit, '(a,f14.7,a,f14.7)') &
    !         'total_charge = ', integrate_environ_density(rho), &
    !         ' total flux throught !soft-spheres = ', SUM(bound%ion_field(:))
    !     !
    !     CALL scalar_product_environ_gradient(field, bound%gradient, rho)

    !     WRITE (environ_unit, '(a,f14.7)') 'actual total flux = ', &
    !         integrate_environ_density!(rho)
    !     !
    !     IF (ideriv .EQ. 1) THEN
    !         !
    !         !----------------------------------------------------------------------------
    !         ! Test functional derivative wrt electronic density of each flux
    !         !
    !         ALLOCATE (analytic_dion_field_drho(bound%ions%number))
    !         !
    !         DO i = 1, bound%ions%number
    !             !
    !             CALL copy_environ_density(bound%dion_field_drho(i), &
    !                                       analytic_dion_field_drho(i))
    !             !
    !         END DO
    !         !
    !         ALLOCATE (fd_dion_field_drho(bound%ions%number))
    !         !
    !         test_function%type_ = 1
    !         test_function%dim = 0
    !         test_function%axis = 3
    !         test_function%spread = 0.3D0
    !         test_function%width = 0.D0
    !         test_function%volume = 1.D0
    !         !
    !         epsilon = 0.000001
    !         !
    !         CALL init_environ_density(cell, delta)
    !         !
    !         CALL init_environ_density(cell, localelectrons%density)
    !         !
    !         !----------------------------------------------------------------------------
    !         ! We are only going to check delta functions along the z axis
    !         ! passing throught the O atom
    !         !
    !         ALLOCATE (test_function%pos(3))
    !         test_function%pos(1) = 11.79D0 / cell%alat
    !         test_function%pos(2) = 12.05D0 / cell%alat
    !         !
    !         DO i = 1, cell%dfft%nr3
    !             !
    !             test_function%pos(3) = DBLE(i - 1) * cell%at(3, 3) / DBLE(cell%dfft%nr3)
    !             CALL density_of_functions(test_function, delta, .TRUE.)
    !             !
    !             fd_dion_field_drho = 0.D0
    !             !
    !             localelectrons%density%of_r = bound%electrons%density%of_r - &
    !                                           epsilon * delta%of_r
    !             !
    !             ! CALL compute_ion_field(bound%ions%number, bound%local_spheres, &
    !             !                        bound%ions, localelectrons, bound%ion_field)
    !             !
    !             fd_dion_field_drho = bound%ion_field
    !             !
    !             localelectrons%density%of_r = bound%electrons%density%of_r + &
    !                                           epsilon * delta%of_r
    !             !
    !             ! CALL compute_ion_field(bound%ions%number, bound%local_spheres, &
    !             !                        bound%ions, localelectrons, bound%ion_field)
    !             !
    !             fd_dion_field_drho = bound%ion_field - fd_dion_field_drho
    !             fd_dion_field_drho = fd_dion_field_drho * 0.5D0 / epsilon
    !             !
    !             IF (ionode) &
    !                 WRITE (environ_unit, *) ' z = ', test_function%pos(3) * cell%alat
    !             !
    !             DO j = 1, bound%ions%number
    !                 tmp = scalar_product_environ_density(analytic_dion_field_drho(j), &
    !                                                      delta)
    !                 !
    !                 IF (ionode) &
    !                     WRITE (environ_unit, '(a,i3,3f20.10)') 'ion = ', j, &
    !                     tmp, fd_dion_field_drho(j), tmp - fd_dion_field_drho(j)
    !                 !
    !             END DO
    !             !
    !         END DO
    !         !
    !         DEALLOCATE (test_function%pos)
    !         CALL destroy_environ_density(delta)
    !         !
    !         DO i = 1, bound%ions%number
    !             CALL destroy_environ_density(analytic_dion_field_drho(i))
    !         END DO
    !         !
    !         DEALLOCATE (analytic_dion_field_drho)
    !         DEALLOCATE (fd_dion_field_drho)
    !         !
    !     ELSE IF (ideriv .EQ. 2) THEN
    !         !
    !         !----------------------------------------------------------------------------
    !         ! Test derivative wrt atomic positions with finite differences
    !         !
    !         CALL compute_ion_field_partial(bound%ions%number, bound%local_spheres, &
    !                                        bound%ions, bound%electrons, &
    !                                        bound%ion_field, bound%partial_of_ion_field, &
    !                                        bound%core%fft)
    !         !
    !         ALLOCATE (analytic_partial_of_ion_field(3, bound%ions%number, &
    !                                                 bound%ions%number))
    !         !
    !         analytic_partial_of_ion_field = bound%partial_of_ion_field
    !         !
    !         dx = 2.0D-3 ! dx expected units: BOHR
    !         !
    !         dx = dx / cell%alat ! convert to alat
    !         !
    !         ALLOCATE (fd_partial_of_ion_field(bound%ions%number))
    !         !
    !         DO i = 1, bound%ions%number
    !             !
    !             DO ipol = 1, 3
    !                 x0 = bound%ions%tau(ipol, i)
    !                 bound%ions%tau(ipol, i) = x0 - dx
    !                 !
    !                 CALL update_environ_ions(bound%ions%number, bound%ions%tau, &
    !                                          bound%ions)
    !                 !
    !                 rho%of_r = bound%electrons%density%of_r + bound%ions%density%of_r
    !                 CALL gradv_h_of_rho_r(rho%of_r, field%of_r)
    !                 !
    !                 ! CALL compute_ion_field(bound%ions%number, bound%local_spheres, &
    !                 !                        bound%ions, bound%electrons, &
    !                 !                        bound%ion_field)
    !                 !
    !                 fd_partial_of_ion_field = bound%ion_field
    !                 !
    !                 bound%ions%tau(ipol, i) = x0 + dx
    !                 !
    !                 CALL update_environ_ions(bound%ions%number, bound%ions%tau, &
    !                                          bound%ions)
    !                 !
    !                 rho%of_r = bound%electrons%density%of_r + bound%ions%density%of_r
    !                 CALL gradv_h_of_rho_r(rho%of_r, field%of_r)
    !                 !
    !                 ! CALL compute_ion_field(bound%ions%number, bound%local_spheres, &
    !                 !                        bound%ions, bound%electrons, &
    !                 !                        bound%ion_field)
    !                 !
    !                 fd_partial_of_ion_field = bound%ion_field - fd_partial_of_ion_field
    !                 fd_partial_of_ion_field = fd_partial_of_ion_field / 2.D0 / dx / cell%alat
    !                 !
    !                 WRITE (environ_unit, *) ' i  = ', i, ' ipol = ', ipol
    !                 !
    !                 WRITE (environ_unit, '(a,10f20.10)') 'analytic     = ', &
    !                     analytic_partial_of_ion_field(ipol, :, i)
    !                 !
    !                 WRITE (environ_unit, '(a,10f20.10)') 'finite-diff  = ', &
    !                     fd_partial_of_ion_field(:)
    !                 !
    !                 WRITE (environ_unit, *) ' '
    !                 !
    !                 bound%ions%tau(ipol, i) = x0
    !                 !
    !             END DO
    !             !
    !         END DO
    !         !
    !         DEALLOCATE (analytic_partial_of_ion_field)
    !         DEALLOCATE (fd_partial_of_ion_field)
    !         !
    !     END IF
    !     !
    !     STOP
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE test_ion_field_derivatives
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !! (DEBUGGING)
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE update_test_boundary(bound, electrons)
    !     !--------------------------------------------------------------------------------
    !     !
    !     USE tools_generate_boundary
    !     !
    !     IMPLICIT NONE
    !     !
    !     TYPE(environ_boundary), INTENT(INOUT) :: bound
    !     TYPE(environ_electrons), INTENT(IN) :: electrons
    !     INTEGER :: debugcubes = 0
    !     CHARACTER(len=100) :: label = 'None'
    !     CHARACTER(len=80) :: sub_name = 'update_test_boundary'
    !     !
    !     SELECT CASE (bound%mode)
    !     CASE ('fa-electronic')
    !         !
    !         ! CALL compute_normal_field(bound%ions, electrons, bound%normal_field)
    !         ! CALL field_aware_density(electrons, bound)
    !         CALL boundary_of_density(bound%density, bound)
    !         !
    !         SELECT CASE (debugcubes)
    !         CASE (2)
    !             label = "standard"
    !             CALL write_cube(bound%scaled, label=label)
    !         CASE (1)
    !             label = "standard"
    !             CALL write_cube(bound%density, label=label)
    !         CASE DEFAULT
    !         END SELECT
    !         !
    !     CASE ('fa-ionic')
    !         !
    !         ! CALL compute_ion_field(bound%ions%number, bound%local_spheres, &
    !         !                        bound%ions, electrons, bound%ion_field)
    !         !
    !         CALL set_soft_spheres(bound, .TRUE.)
    !         CALL boundary_of_functions(bound%ions%number, bound%soft_spheres, bound)
    !         !
    !     CASE ('electronic')
    !         CALL boundary_of_density(electrons%density, bound)
    !     CASE DEFAULT
    !         CALL env_errore(sub_name, 'Unrecognized boundary mode', 1)
    !     END SELECT
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE update_test_boundary
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !! (DEBUGGING)
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE print_ionic(ions)
    !     !--------------------------------------------------------------------------------
    !     !
    !     TYPE(environ_ions), INTENT(IN) :: ions
    !     INTEGER :: i
    !     !
    !     DO i = 1, ions%number
    !         PRINT *, "ATOM", i, ions%tau(1, i), ions%tau(2, i), ions%tau(3, i)
    !     END DO
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE print_ionic
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !! (DEBUGGING)
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE test_energy_derivatives(ideriv, bound)
    !     !--------------------------------------------------------------------------------
    !     !
    !     USE tools_generate_boundary
    !     USE utils_ions, ONLY: update_environ_ions
    !     USE embedding_surface, ONLY: calc_esurface, calc_desurface_dboundary
    !     USE embedding_volume, ONLY: calc_evolume, calc_devolume_dboundary
    !     !
    !     IMPLICIT NONE
    !     !
    !     INTEGER, INTENT(IN) :: ideriv
    !     TYPE(environ_boundary), INTENT(INOUT) :: bound
    !     !
    !     TYPE(environ_cell), POINTER :: cell
    !     !
    !     INTEGER :: i, ipol
    !     REAL(DP) :: epsilon
    !     REAL(DP) :: localpressure, localsurface_tension
    !     REAL(DP) :: de_fd, de_analytic, etmp
    !     REAL(DP) :: dx, x0, force(3), ssforce(3)
    !     REAL(DP) :: flux
    !     REAL(DP), ALLOCATABLE :: tau0(:, :)
    !     TYPE(environ_density) :: de_dboundary
    !     TYPE(environ_density) :: vanalytic
    !     TYPE(environ_functions) :: test_function
    !     TYPE(environ_density) :: delta
    !     TYPE(environ_electrons) :: localelectrons
    !     TYPE(environ_gradient) :: partial
    !     !
    !     CHARACTER(LEN=80) :: sub_name = 'test_energy_derivatives'
    !     !
    !     cell => bound%scaled%cell
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! Compute the field and the field-aware interface (Need to check)
    !     !
    !     CALL update_test_boundary(bound, bound%electrons)
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! Compute functional derivative of the energy wrt interface
    !     !
    !     CALL init_environ_density(cell, de_dboundary)
    !     !
    !     localpressure = -0.35
    !     CALL calc_devolume_dboundary(localpressure, bound, de_dboundary)
    !     !
    !     localsurface_tension = 0.0
    !     ! CALL calc_desurface_dboundary(localsurface_tension, bound, de_dboundary)
    !     !
    !     IF (ideriv .EQ. 1) THEN
    !         !
    !         !----------------------------------------------------------------------------
    !         ! Compute functional derivative wrt electronic density
    !         !
    !         ! IF (bound%mode .EQ. 'fa-ionic') THEN
    !         !     CALL compute_dion_field_drho(bound%ions%number, bound%local_spheres, &
    !         !                                  bound%dion_field_drho, bound%core%fft)
    !         ! END IF
    !         !
    !         CALL init_environ_density(cell, vanalytic)
    !         !
    !         IF (bound%field_aware) THEN
    !             ! CALL field_aware_de_drho(bound, de_dboundary, vanalytic)
    !             CALL env_errore('field-aware6', 'Option not yet implimented ', 1)
    !         ELSE
    !             vanalytic%of_r = bound%dscaled%of_r * de_dboundary%of_r
    !         END IF
    !         !
    !         !----------------------------------------------------------------------------
    !         ! Loop over gridpoints with rhoelec + or - a delta function
    !         !
    !         test_function%type_ = 1
    !         test_function%dim = 0
    !         test_function%axis = 3
    !         test_function%spread = 0.4D0
    !         test_function%width = 0.D0
    !         test_function%volume = 1.D0
    !         !
    !         epsilon = 0.000008
    !         !
    !         CALL init_environ_density(cell, delta)
    !         !
    !         CALL init_environ_density(cell, localelectrons%density)
    !         !
    !         !----------------------------------------------------------------------------
    !         ! We are only going to check delta functions along the z-axis
    !         ! passing throught !the O atom
    !         !
    !         ALLOCATE (test_function%pos(3))
    !         test_function%pos(1) = 0.0 / cell%alat
    !         test_function%pos(2) = 0.0 / cell%alat
    !         !
    !         DO i = 1, cell%dfft%nr3
    !             de_fd = 0.D0
    !             test_function%pos(3) = DBLE(i - 1) * cell%at(3, 3) / DBLE(cell%dfft%nr3)
    !             !
    !             IF (ionode) &
    !                 WRITE (program_unit, '(a,f14.7)') ' z = ', test_function%pos(3)
    !             !
    !             ! IF (test_function%pos(3) .LE. 0.365) CYCLE
    !             !
    !             ! IF (test_function%pos(3) * cell%alat .LE. 4.5 .OR. &
    !             !     test_function%pos(3) * cell%alat .GE. 7.0) CYCLE
    !             !
    !             CALL density_of_functions(test_function, delta, .TRUE.)
    !             !
    !             localelectrons%density%of_r = bound%electrons%density%of_r + &
    !                                           epsilon * delta%of_r
    !             !
    !             CALL update_test_boundary(bound, localelectrons)
    !             !
    !             CALL calc_evolume(localpressure, bound, etmp)
    !             de_fd = de_fd + etmp
    !             !
    !             CALL calc_esurface(localsurface_tension, bound, etmp)
    !             de_fd = de_fd + etmp
    !             !
    !             ! IF (ionode) &
    !             !     WRITE (environ_unit, '(2f20.10)') bound%surface, bound%volume
    !             !
    !             localelectrons%density%of_r = bound%electrons%density%of_r - &
    !                                           epsilon * delta%of_r
    !             !
    !             CALL update_test_boundary(bound, localelectrons)
    !             !
    !             CALL calc_evolume(localpressure, bound, etmp)
    !             de_fd = de_fd - etmp
    !             !
    !             CALL calc_esurface(localsurface_tension, bound, etmp)
    !             de_fd = de_fd - etmp
    !             !
    !             ! IF (ionode) &
    !             !     WRITE (environ_unit, '(2f20.10)') bound%surface, bound%volume
    !             !
    !             de_fd = de_fd * 0.5D0 / epsilon
    !             !
    !             de_analytic = scalar_product_environ_density(vanalytic, delta)
    !             !
    !             IF (ionode) &
    !                 WRITE (environ_unit, '(a,i14.7,3f20.10)') ' z = ', i, &
    !                 de_analytic, de_fd, de_analytic - de_fd
    !             !
    !             FLUSH (environ_unit)
    !             !
    !             ! STOP
    !             !
    !         END DO
    !         !
    !         DEALLOCATE (test_function%pos)
    !         CALL destroy_environ_density(delta)
    !         !
    !         CALL destroy_environ_density(vanalytic)
    !         !
    !     ELSE IF (ideriv .EQ. 2) THEN
    !         PRINT *, 'ideriv = 2'
    !         !
    !         !----------------------------------------------------------------------------
    !         ! Compute partial derivative with respect to ionic positions
    !         !
    !         ! CALCULATE ion field partials in advance
    !         !
    !         ! IF (bound%mode .EQ. 'fa-ionic') THEN
    !         !     !
    !         !     IF (ionode) &
    !         !         WRITE (program_unit, '(1X,a)') 'outside compute_ion_field_partial'
    !         !     !
    !         !     CALL compute_ion_field_partial(bound%ions%number, &
    !         !                                    bound%local_spheres, &
    !         !                                    bound%ions, bound%electrons, &
    !         !                                    bound%ion_field, &
    !         !                                    bound%partial_of_ion_field, &
    !         !                                    bound%core%fft)
    !         !     !
    !         ! END IF
    !         !
    !         dx = 2.0D-3 ! dx expected units: BOHR
    !         !
    !         dx = dx / cell%alat ! convert to alat
    !         !
    !         CALL init_environ_gradient(cell, partial)
    !         ALLOCATE (tau0(3, bound%ions%number))
    !         !
    !         tau0(:, :) = bound%ions%tau(:, :)
    !         !
    !         DO i = 1, bound%ions%number
    !             !
    !             bound%ions%tau(:, :) = tau0(:, :)
    !             ! CALCULATE FORCE FIRST, reset the ionic positions
    !             !
    !             CALL update_environ_ions(bound%ions%number, bound%ions%tau, bound%ions)
    !             CALL update_test_boundary(bound, bound%electrons)
    !             !
    !             IF (bound%mode .EQ. 'fa-ionic') THEN
    !                 CALL calc_dboundary_dions(i, bound, partial)
    !                 !
    !                 ssforce = -scalar_product_environ_gradient_density(partial, &
    !                                                                    de_dboundary)
    !                 !
    !             END IF
    !             !
    !             ! IF (bound%field_aware) &
    !             !     CALL field_aware_dboundary_dions(i, bound, partial)
    !             !
    !             force = -scalar_product_environ_gradient_density(partial, de_dboundary)
    !             !
    !             DO ipol = 1, 3
    !                 !
    !                 de_fd = 0.D0 ! FINITE DIFFERENCE VALUE STORED HERE
    !                 bound%ions%tau(:, :) = tau0(:, :) ! RESET IONS
    !                 !
    !                 bound%ions%tau(ipol, i) = tau0(ipol, i) - dx ! MINUS dx
    !                 !
    !                 CALL update_environ_ions(bound%ions%number, bound%ions%tau, &
    !                                          bound%ions)
    !                 !
    !                 CALL update_test_boundary(bound, bound%electrons)
    !                 !
    !                 CALL calc_evolume(localpressure, bound, etmp)
    !                 de_fd = de_fd - etmp
    !                 CALL calc_esurface(localsurface_tension, bound, etmp)
    !                 de_fd = de_fd - etmp
    !                 !
    !                 ! PLUS dx
    !                 bound%ions%tau(ipol, i) = tau0(ipol, i) + dx
    !                 !
    !                 CALL update_environ_ions(bound%ions%number, bound%ions%tau, &
    !                                          bound%ions)
    !                 !
    !                 CALL update_test_boundary(bound, bound%electrons)
    !                 !
    !                 CALL calc_evolume(localpressure, bound, etmp)
    !                 de_fd = de_fd + etmp
    !                 CALL calc_esurface(localsurface_tension, bound, etmp)
    !                 de_fd = de_fd + etmp
    !                 !
    !                 de_fd = de_fd * (-0.5D0) / dx / cell%alat
    !                 ! force is negative of the energy derivative
    !                 !
    !                 IF (bound%mode .EQ. 'fa-ionic') THEN
    !                     !
    !                     IF (ionode) &
    !                         WRITE (environ_unit, '(a,i3,a,i3,4f20.10)') ' i = ', i, &
    !                         ' ipol != ', ipol, force(ipol), ssforce(ipol), de_fd, &
    !                         force(ipol) - de_fd
    !                     !
    !                 ELSE
    !                     !
    !                     IF (ionode) &
    !                         WRITE (environ_unit, '(a,i3,a,i3,3f20.10)') ' i = ', i, &
    !                         ' ipol != ', ipol, force(ipol), de_fd, force(ipol) - de_fd
    !                     !
    !                 END IF
    !                 !
    !                 ! STOP
    !                 !
    !             END DO
    !             !
    !         END DO
    !         !
    !         CALL destroy_environ_gradient(partial)
    !         !
    !     END IF
    !     !
    !     CALL destroy_environ_density(de_dboundary)
    !     !
    !     STOP
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE test_energy_derivatives
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !! (DEBUGGING)
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE extract_boundary_data(bound)
    !     !------------------------------------------------------------------------------------
    !     !
    !     IMPLICIT NONE
    !     !
    !     TYPE(environ_boundary), INTENT(IN) :: bound
    !     !
    !     TYPE(environ_cell), POINTER :: cell
    !     !
    !     INTEGER :: i, ipol
    !     REAL(DP) :: epsilon
    !     REAL(DP) :: localpressure, localsurface_tension
    !     REAL(DP) :: de_fd, de_analytic, etmp
    !     REAL(DP) :: dx, x0, force(3)
    !     TYPE(environ_density) :: de_dboundary
    !     TYPE(environ_density) :: vanalytic
    !     TYPE(environ_functions) :: test_function
    !     TYPE(environ_density) :: delta
    !     TYPE(environ_electrons) :: localelectrons
    !     TYPE(environ_gradient) :: partial
    !     !
    !     CHARACTER(len=100) :: strg = 'e'
    !     CHARACTER(len=100) :: strl = 'i'
    !     CHARACTER(LEN=80) :: sub_name = 'extract_boundary_data'
    !     !
    !     cell => bound%scaled%cell
    !     !
    !     IF (ionode) WRITE (program_unit, '(1X,a)') 'extract_boundary_data'
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! Compute the field and the field-aware interface
    !     !
    !     CALL write_cube(bound%electrons%density, label=strg)
    !     CALL write_cube(bound%ions%density, label=strl)
    !     !
    !     STOP
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE extract_boundary_data
    ! !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE environ_debugging
!----------------------------------------------------------------------------------------
