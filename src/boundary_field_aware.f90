!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE boundary_field_aware
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE compute_ion_field(nsoft_spheres, soft_spheres, ions, electrons, ion_field)
    !     !--------------------------------------------------------------------------------
    !     !
    !     IMPLICIT NONE
    !     !
    !     INTEGER, INTENT(IN) :: nsoft_spheres
    !     TYPE(environ_functions), INTENT(IN) :: soft_spheres(nsoft_spheres)
    !     TYPE(environ_ions), INTENT(IN) :: ions
    !     TYPE(environ_electrons), INTENT(IN) :: electrons
    !     REAL(DP), INTENT(OUT) :: ion_field(nsoft_spheres)
    !     !
    !     TYPE(environ_cell), POINTER :: cell
    !     !
    !     INTEGER :: i, j
    !     CHARACTER(LEN=80) :: sub_name = 'compute_ion_field'
    !     !
    !     TYPE(environ_density), ALLOCATABLE :: local(:)
    !     !
    !     TYPE(environ_density) :: aux, prod
    !     TYPE(environ_gradient) :: gradaux, field
    !     !
    !     cell => ions%density%cell
    !     !
    !     ALLOCATE (local(nsoft_spheres))
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! Compute field-independent soft spheres and gradients
    !     !
    !     DO i = 1, nsoft_spheres
    !         CALL init_environ_density(cell, local(i))
    !         CALL density_of_functions(soft_spheres(i), local(i), .FALSE.)
    !     END DO
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! Compute field
    !     !
    !     CALL init_environ_density(cell, aux)
    !     aux%of_r = electrons%density%of_r + ions%density%of_r
    !     !
    !     CALL init_environ_gradient(cell, field)
    !     CALL gradv_h_of_rho_r(aux%of_r, field%of_r)
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! Compute field flux
    !     !
    !     ion_field = 0.D0
    !     !
    !     CALL init_environ_density(cell, prod)
    !     CALL init_environ_gradient(cell, gradaux)
    !     !
    !     DO i = 1, nsoft_spheres
    !         !
    !         !----------------------------------------------------------------------------
    !         ! Compute product of other soft-spheres
    !         !
    !         prod%of_r = 1.D0
    !         !
    !         DO j = 1, nsoft_spheres
    !             IF (j .EQ. i) CYCLE
    !             prod%of_r = prod%of_r * local(j)%of_r
    !         END DO
    !         !
    !         !----------------------------------------------------------------------------
    !         ! Compute field flux through soft-sphere interface
    !         !
    !         CALL gradient_of_functions(soft_spheres(i), gradaux, .TRUE.)
    !         !
    !         CALL scalar_product_environ_gradient(field, gradaux, aux) ! here aux is the !normal field
    !         !
    !         aux%of_r = -aux%of_r * prod%of_r
    !         !
    !         ion_field(i) = integrate_environ_density(aux)
    !         !
    !     END DO
    !     !
    !     CALL destroy_environ_gradient(gradaux)
    !     CALL destroy_environ_density(prod)
    !     !
    !     CALL destroy_environ_gradient(field)
    !     CALL destroy_environ_density(aux)
    !     !
    !     DO i = 1, nsoft_spheres
    !         CALL destroy_environ_density(local(i))
    !     END DO
    !     !
    !     DEALLOCATE (local)
    !     !
    !     RETURN
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE compute_ion_field
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE compute_ion_field_partial(nsoft_spheres, soft_spheres, ions, &
    !                                      electrons, ion_field, partial_of_ion_field, fft)
    !     !--------------------------------------------------------------------------------
    !     !
    !     USE core_fft, ONLY: hessv_h_of_rho_r
    !     !
    !     IMPLICIT NONE
    !     !
    !     INTEGER, INTENT(IN) :: nsoft_spheres
    !     TYPE(environ_functions), INTENT(IN) :: soft_spheres(nsoft_spheres)
    !     TYPE(environ_ions), INTENT(IN) :: ions
    !     TYPE(environ_electrons), INTENT(IN) :: electrons
    !     TYPE(fft_core), INTENT(IN) :: fft
    !     REAL(DP), INTENT(OUT) :: ion_field(nsoft_spheres)
    !     REAL(DP), INTENT(OUT) :: partial_of_ion_field(3, nsoft_spheres, nsoft_spheres)
    !     !
    !     INTEGER, POINTER :: nnr, ir_end, deriv
    !     TYPE(environ_cell), POINTER :: cell
    !     !
    !     INTEGER :: i, j, k, ipol
    !     CHARACTER(LEN=80) :: sub_name = 'compute_ion_field'
    !     !
    !     TYPE(environ_density), ALLOCATABLE :: local(:)
    !     TYPE(environ_gradient), ALLOCATABLE :: gradlocal(:)
    !     TYPE(environ_hessian) :: hesslocal
    !     !
    !     TYPE(environ_density) :: aux, prod
    !     TYPE(environ_gradient) :: gradaux, field
    !     TYPE(environ_hessian) :: hessaux
    !     !
    !     cell => ions%density%cell
    !     nnr => cell%nnr
    !     ir_end => cell%ir_end
    !     !
    !     ALLOCATE (local(nsoft_spheres))
    !     ALLOCATE (gradlocal(nsoft_spheres))
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! Compute field-independent soft spheres and gradients
    !     !
    !     DO i = 1, nsoft_spheres
    !         CALL init_environ_density(cell, local(i))
    !         CALL init_environ_gradient(cell, gradlocal(i))
    !         CALL density_of_functions(soft_spheres(i), local(i), .FALSE.)
    !         CALL gradient_of_functions(soft_spheres(i), gradlocal(i), .FALSE.)
    !     END DO
    !     !
    !     CALL init_environ_hessian(cell, hesslocal)
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! Compute field
    !     !
    !     CALL init_environ_density(cell, aux)
    !     aux%of_r = electrons%density%of_r + ions%density%of_r
    !     !
    !     CALL init_environ_gradient(cell, field)
    !     CALL gradv_h_of_rho_r(aux%of_r, field%of_r)
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! Compute field flux
    !     !
    !     ion_field = 0.D0
    !     partial_of_ion_field = 0.D0
    !     !
    !     CALL init_environ_density(cell, prod)
    !     CALL init_environ_gradient(cell, gradaux)
    !     CALL init_environ_hessian(cell, hessaux)
    !     !
    !     DO i = 1, nsoft_spheres
    !         !
    !         !----------------------------------------------------------------------------
    !         ! Compute product of other soft-spheres
    !         !
    !         prod%of_r = 1.D0
    !         DO j = 1, nsoft_spheres
    !             IF (j .EQ. i) CYCLE
    !             prod%of_r = prod%of_r * local(j)%of_r
    !         END DO
    !         !
    !         !----------------------------------------------------------------------------
    !         ! Compute field flux through soft-sphere interface
    !         !
    !         CALL scalar_product_environ_gradient(field, gradlocal(i), aux)
    !         ! here aux is the normal field
    !         !
    !         aux%of_r = -aux%of_r * prod%of_r
    !         !
    !         ion_field(i) = integrate_environ_density(aux)
    !         !
    !         !----------------------------------------------------------------------------
    !         ! Compute partial derivatives of field flux wrt ionic positions
    !         !
    !         DO j = 1, nsoft_spheres
    !             !
    !             !------------------------------------------------------------------------
    !             ! Compute hessian of poisson potential of individual nuclei #TODO
    !             !
    !             ! DEBUG try changing indices.. was j
    !             CALL density_of_functions(ions%smeared_ions(j), aux, .TRUE.)
    !             ! THIS STEP SHOULD BE MOVED OUT OF THIS LOOP
    !             !
    !             CALL hessv_h_of_rho_r(aux%of_r, hesslocal%of_r, fft)
    !             ! THIS STEP SHOULD BE MOVED OUT OF THIS LOOP
    !             !
    !             CALL scalar_product_environ_hessian(hesslocal, gradlocal(i), gradaux)
    !             !
    !             partial_of_ion_field(:, i, j) = &
    !                 partial_of_ion_field(:, i, j) - &
    !                 scalar_product_environ_gradient_density(gradaux, prod)
    !             !
    !             IF (j .EQ. i) THEN
    !                 !
    !                 !--------------------------------------------------------------------
    !                 ! Hessian of soft-sphere times the field
    !                 !
    !                 CALL hessian_of_functions(soft_spheres(i), hessaux, .TRUE.)
    !                 !
    !                 CALL scalar_product_environ_hessian(hessaux, field, gradaux)
    !                 !
    !                 partial_of_ion_field(:, i, j) = &
    !                     partial_of_ion_field(:, i, j) + &
    !                     scalar_product_environ_gradient_density(gradaux, prod)
    !                 !
    !             ELSE
    !                 !
    !                 !--------------------------------------------------------------------
    !                 ! ion field times gradient of different soft-sphere
    !                 !
    !                 ! DEBUG try changing indices...was i
    !                 CALL scalar_product_environ_gradient(gradlocal(i), field, aux)
    !                 ! here aux !is the normal field
    !                 !
    !                 DO k = 1, nsoft_spheres
    !                     IF (k .EQ. j .OR. k .EQ. i) CYCLE
    !                     aux%of_r = aux%of_r * local(k)%of_r
    !                 END DO
    !                 !
    !                 ! DEBUG try changing indices... was j
    !                 partial_of_ion_field(:, i, j) = &
    !                     partial_of_ion_field(:, i, j) + &
    !                     scalar_product_environ_gradient_density(gradlocal(j), aux)
    !                 !
    !             END IF
    !             !
    !         END DO
    !         !
    !     END DO
    !     !
    !     CALL destroy_environ_gradient(field)
    !     CALL destroy_environ_density(prod)
    !     CALL destroy_environ_density(aux)
    !     CALL destroy_environ_gradient(gradaux)
    !     CALL destroy_environ_hessian(hessaux)
    !     !
    !     CALL destroy_environ_hessian(hesslocal)
    !     !
    !     DO i = 1, nsoft_spheres
    !         CALL destroy_environ_gradient(gradlocal(i))
    !         CALL destroy_environ_density(local(i))
    !     END DO
    !     !
    !     DEALLOCATE (gradlocal)
    !     DEALLOCATE (local)
    !     !
    !     RETURN
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE compute_ion_field_partial
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE compute_dion_field_drho(nsoft_spheres, soft_spheres, dion_field_drho, fft)
    !     !--------------------------------------------------------------------------------
    !     !
    !     USE core_fft, ONLY: field_of_gradrho
    !     !
    !     IMPLICIT NONE
    !     !
    !     INTEGER, INTENT(IN) :: nsoft_spheres
    !     TYPE(environ_functions), INTENT(IN) :: soft_spheres(nsoft_spheres)
    !     TYPE(fft_core), INTENT(IN) :: fft
    !     TYPE(environ_density), INTENT(INOUT) :: dion_field_drho(nsoft_spheres)
    !     !
    !     TYPE(environ_cell), POINTER :: cell
    !     !
    !     INTEGER :: i, j, ipol
    !     CHARACTER(LEN=80) :: sub_name = 'compute_dion_field_drho'
    !     !
    !     TYPE(environ_density), ALLOCATABLE :: local(:)
    !     !
    !     TYPE(environ_density) :: prod
    !     TYPE(environ_gradient) :: gradaux
    !     !
    !     IF (nsoft_spheres .LT. 1) CALL errore(sub_name, 'Missing soft-spheres', 1)
    !     cell => dion_field_drho(1)%cell
    !     !
    !     ALLOCATE (local(nsoft_spheres))
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! Compute field-independent soft spheres and gradients
    !     !
    !     DO i = 1, nsoft_spheres
    !         CALL init_environ_density(cell, local(i))
    !         CALL density_of_functions(soft_spheres(i), local(i), .FALSE.)
    !     END DO
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! Compute field flux
    !     !
    !     CALL init_environ_density(cell, prod)
    !     CALL init_environ_gradient(cell, gradaux)
    !     !
    !     DO i = 1, nsoft_spheres
    !         !
    !         !----------------------------------------------------------------------------
    !         ! Compute product of other soft-spheres
    !         !
    !         prod%of_r = 1.D0
    !         !
    !         DO j = 1, nsoft_spheres
    !             IF (j .EQ. i) CYCLE
    !             prod%of_r = prod%of_r * local(j)%of_r
    !         END DO
    !         !
    !         !----------------------------------------------------------------------------
    !         ! Compute functional derivative of field flux wrt electronic density
    !         !
    !         CALL gradient_of_functions(soft_spheres(i), gradaux, .TRUE.)
    !         !
    !         DO ipol = 1, 3
    !             gradaux%of_r(ipol, :) = gradaux%of_r(ipol, :) * prod%of_r(:)
    !         END DO
    !         !
    !         CALL field_of_gradrho(gradaux%of_r, dion_field_drho(i)%of_r, fft)
    !     END DO
    !     !
    !     CALL destroy_environ_gradient(gradaux)
    !     !
    !     DO i = 1, nsoft_spheres
    !         CALL destroy_environ_density(local(i))
    !     END DO
    !     !
    !     DEALLOCATE (local)
    !     !
    !     RETURN
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE compute_dion_field_drho
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE compute_normal_field(ions, electrons, normal_field)
    !     !--------------------------------------------------------------------------------
    !     !
    !     IMPLICIT NONE
    !     !
    !     TYPE(environ_ions), INTENT(IN) :: ions
    !     TYPE(environ_electrons), INTENT(IN) :: electrons
    !     TYPE(environ_density), INTENT(INOUT) :: normal_field
    !     !
    !     TYPE(environ_cell), POINTER :: cell
    !     TYPE(environ_density) :: rho
    !     TYPE(environ_gradient) :: field, gradrho
    !     INTEGER :: i
    !     !
    !     cell => normal_field%cell
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! Compute field
    !     !
    !     CALL init_environ_density(cell, rho)
    !     rho%of_r = electrons%density%of_r + ions%density%of_r
    !     !
    !     CALL init_environ_gradient(cell, field)
    !     CALL gradv_h_of_rho_r(rho%of_r, field%of_r)
    !     !
    !     CALL init_environ_gradient(cell, gradrho)
    !     CALL external_gradient(electrons%density%of_r, gradrho%of_r)
    !     !
    !     ! temporary factor for testing
    !     ! CALL update_gradient_modulus(gradrho)
    !     ! need to divide by the modulus factor ( but don't touch if zero )
    !     ! loop elementwise to be safe
    !     ! DO i = 1, cell%nnr
    !     !     !
    !     !     IF (gradrho%modulus%of_r(i) < 1E-10) THEN
    !     !         gradrho%modulus%of_r(i) = 1.D0
    !     !     END IF
    !     !     !
    !     ! END DO
    !     ! !
    !     CALL scalar_product_environ_gradient(field, gradrho, normal_field)
    !     !
    !     ! temporary modification of the normal field to remove density dependency
    !     ! normal_field%of_r = normal_field%of_r / gradrho%modulus%of_r
    !     !
    !     CALL destroy_environ_gradient(gradrho)
    !     CALL destroy_environ_gradient(field)
    !     CALL destroy_environ_density(rho)
    !     !
    !     RETURN
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE compute_normal_field
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !!
    ! !------------------------------------------------------------------------------------
    ! FUNCTION scaling_of_field(field_factor, charge_asymmetry, field_max, field_min, &
    !                           ion_field)
    !     !--------------------------------------------------------------------------------
    !     !
    !     IMPLICIT NONE
    !     !
    !     REAL(DP) :: scaling_of_field
    !     REAL(DP) :: field_factor, charge_asymmetry, field_max, field_min, ion_field
    !     !
    !     REAL(DP) :: field, fact, arg
    !     !
    !     field = ABS(ion_field)
    !     !
    !     fact = (charge_asymmetry - SIGN(1.D0, ion_field))**2 * field_factor
    !     !
    !     IF (field .LE. field_min) THEN
    !         scaling_of_field = 0.D0
    !     ELSE IF (field .LE. field_max) THEN
    !         arg = tpi * (field - field_min) / (field_max - field_min)
    !         scaling_of_field = (arg - SIN(arg)) / tpi
    !     ELSE
    !         scaling_of_field = 1.D0
    !     END IF
    !     !
    !     scaling_of_field = 1.D0 - scaling_of_field * fact
    !     !
    !     RETURN
    !     !
    !     !--------------------------------------------------------------------------------
    ! END FUNCTION scaling_of_field
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !!
    ! !------------------------------------------------------------------------------------
    ! FUNCTION dscaling_of_field(field_factor, charge_asymmetry, field_max, &
    !                            field_min, ion_field)
    !     !--------------------------------------------------------------------------------
    !     !
    !     IMPLICIT NONE
    !     !
    !     REAL(DP) :: dscaling_of_field
    !     REAL(DP) :: field_factor, charge_asymmetry, field_max, field_min, ion_field
    !     !
    !     REAL(DP) :: field, fact, arg
    !     !
    !     field = ABS(ion_field)
    !     !
    !     fact = (charge_asymmetry - SIGN(1.D0, ion_field))**2 * field_factor
    !     !
    !     IF (field .LE. field_min .OR. field .GT. field_max) THEN
    !         dscaling_of_field = 0.D0
    !     ELSE
    !         arg = tpi * (field - field_min) / (field_max - field_min)
    !         dscaling_of_field = (1.D0 - COS(arg)) / (field_max - field_min)
    !     END IF
    !     !
    !     dscaling_of_field = -dscaling_of_field * fact * SIGN(1.D0, ion_field)
    !     !
    !     RETURN
    !     !
    !     !--------------------------------------------------------------------------------
    ! END FUNCTION dscaling_of_field
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !! @brief Outputs a density object representing part of the field-aware scaling
    ! !!
    ! !! Scales the density in the field-aware self-consistent model
    ! !! the function is defined in e^-f(E_norm), which multiplies the electronic
    ! !! density to represent the new effective density. NOTE that the exponent is
    ! !! applied here! The f function is the same switching function introduced
    ! !! by Andreussi et al. for the original SCCS
    ! !!
    ! !! @param[in]      field_factor      the muliplicity factor
    ! !! @param[in]      field_max         maximum electric field cutoff
    ! !! @param[in]      field_min         minimum electric field cutoff
    ! !! @param[in]      normal_field      the field to define the scaling object
    ! !! @param[inout]   scaled           scaling
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE scaling_of_density(field_factor, charge_asymmetry, field_max, &
    !                               field_min, normal_field, scaled)
    !     !--------------------------------------------------------------------------------
    !     !
    !     IMPLICIT NONE
    !     !
    !     REAL(DP), INTENT(IN) :: field_factor, field_max, field_min, charge_asymmetry
    !     REAL(DP) :: arg, fact
    !     REAL(DP) :: field
    !     !
    !     TYPE(environ_density), INTENT(IN) :: normal_field
    !     TYPE(environ_density), INTENT(INOUT) :: scaled
    !     !
    !     INTEGER, POINTER :: nnr
    !     INTEGER :: i
    !     !
    !     nnr => normal_field%cell%nnr
    !     !
    !     DO i = 1, nnr
    !         field = normal_field%of_r(i)
    !         fact = (charge_asymmetry - SIGN(1.D0, field))**2 * field_factor
    !         !
    !         IF (ABS(field) .LE. field_min) THEN
    !             arg = 0.D0
    !         ELSE IF (ABS(field) .LE. field_max) THEN
    !             arg = tpi * (ABS(field) - (field_min)) / ((field_max) - (field_min))
    !             arg = (arg - SIN(arg)) / tpi
    !         ELSE
    !             arg = 1.D0
    !         END IF
    !         !
    !         scaled%of_r(i) = EXP(arg * fact * (-1.D0))
    !         !
    !     END DO
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE scaling_of_density
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !! @brief Outputs a density object representing part of the field-aware
    ! !! scaling derivative
    ! !!
    ! !! Applies the derivative of the scaling of density function
    ! !! Note that scaling of density is an exponential function, exp(-fx)
    ! !! and the derivative function here only returns the derivative of f,
    ! !! which should then be multiplied with the result of scaling of density
    ! !! for the correct result
    ! !!
    ! !! @param[in]      field_factor      the muliplicity factor
    ! !! @param[in]      field_max         maximum electric field cutoff
    ! !! @param[in]      field_min         minimum electric field cutoff
    ! !! @param[in]      normal_field      the field to define the scaling object
    ! !! @param[inout]   dscaled           scaling derivative (wrt normal field)
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE dscaling_of_density(field_factor, charge_asymmetry, field_max, &
    !                                field_min, normal_field, dscaled)
    !     !--------------------------------------------------------------------------------
    !     !
    !     IMPLICIT NONE
    !     !
    !     REAL(DP), INTENT(IN) :: field_factor, field_max, field_min, charge_asymmetry
    !     REAL(DP) :: arg, fact
    !     REAL(DP) :: field
    !     !
    !     TYPE(environ_density), INTENT(IN) :: normal_field
    !     TYPE(environ_density), INTENT(INOUT) :: dscaled
    !     !
    !     INTEGER, POINTER :: nnr
    !     INTEGER :: i
    !     !
    !     nnr => normal_field%cell%nnr
    !     !
    !     DO i = 1, nnr
    !         field = normal_field%of_r(i)
    !         fact = (charge_asymmetry - SIGN(1.D0, field))**2 * field_factor
    !         !
    !         IF (ABS(field) .LE. field_min) THEN
    !             arg = 0.D0
    !         ELSE IF (ABS(field) .LE. field_max) THEN
    !             arg = tpi * (ABS(field) - (field_min)) / ((field_max) - (field_min))
    !             arg = (1.D0 - COS(arg)) / ((field_max) - (field_min))
    !             arg = arg * SIGN(1.D0, field)
    !         ELSE
    !             arg = 0.D0
    !         END IF
    !         !
    !         dscaled%of_r(i) = arg * fact * (-1.D0)
    !         !
    !     END DO
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE dscaling_of_density
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE field_aware_density(electrons, boundary)
    !     !--------------------------------------------------------------------------------
    !     !
    !     IMPLICIT NONE
    !     !
    !     TYPE(environ_boundary), INTENT(INOUT) :: boundary
    !     TYPE(environ_density) :: scaling
    !     !
    !     TYPE(environ_cell), POINTER :: cell
    !     TYPE(environ_electrons), INTENT(IN) :: electrons
    !     REAL(DP) :: field_factor, charge_asymmetry, field_max, field_min
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! call scaling of density to get the scaling factor and add this to
    !     ! the density
    !     !
    !     cell => boundary%electrons%density%cell
    !     field_factor = boundary%field_factor
    !     charge_asymmetry = boundary%charge_asymmetry
    !     field_max = boundary%field_max
    !     field_min = boundary%field_min
    !     !
    !     CALL init_environ_density(cell, scaling)
    !     !
    !     CALL scaling_of_density(field_factor, charge_asymmetry, field_max, field_min, &
    !                             boundary%normal_field, scaling)
    !     !
    !     boundary%density%of_r = scaling%of_r * electrons%density%of_r
    !     !
    !     CALL destroy_environ_density(scaling)
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE field_aware_density
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !! @brief d normal_field / d rho
    ! !!
    ! !! Intermediate function that calculates the functional derivative of
    ! !! Normal electric field with respect to the electric density. To
    ! !! evaluate one needs an arbitrary function. Implemented here, the
    ! !! function is only evaluated at call, and discarded afterwards.
    ! !! Note that this value is only useful for debugging.
    ! !!
    ! !! @param[in]      boundary       the boundary object contains the normal field
    ! !! @param[in]      funct          to reduce this functional, one needs to define
    ! !!                                some arbitrary function to integrate over.
    ! !! @param[inout]   dfield_drho    functional derivative
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE delectronic_field_drho(boundary, funct, dfield_drho)
    !     !--------------------------------------------------------------------------------
    !     !
    !     USE core_fft, ONLY: field_of_gradrho
    !     !
    !     IMPLICIT NONE
    !     !
    !     TYPE(environ_boundary), INTENT(IN) :: boundary
    !     TYPE(environ_functions), INTENT(IN) :: funct
    !     TYPE(environ_density), INTENT(OUT) :: dfield_drho
    !     !
    !     TYPE(environ_cell), POINTER :: cell
    !     TYPE(fft_core), POINTER :: fft
    !     TYPE(environ_gradient) :: gradaux
    !     TYPE(environ_density) :: gaussian
    !     INTEGER :: ipol
    !     !
    !     cell => boundary%electrons%density%cell
    !     fft => boundary%core%fft
    !     !
    !     CALL init_environ_density(cell, dfield_drho)
    !     CALL init_environ_density(cell, gaussian)
    !     CALL init_environ_gradient(cell, gradaux)
    !     !
    !     CALL external_gradient(boundary%electrons%density%of_r, gradaux%of_r)
    !     CALL density_of_functions(funct, gaussian, .TRUE.)
    !     !
    !     DO ipol = 1, 3
    !         gradaux%of_r(ipol, :) = gradaux%of_r(ipol, :) * gaussian%of_r(:)
    !     END DO
    !     !
    !     CALL field_of_gradrho(gradaux%of_r, dfield_drho%of_r, fft)
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! use gaussian container to store auxiliary density to be added to
    !     ! dfield_drho
    !     !
    !     gaussian%of_r = gaussian%of_r * 2.D0 * tpi
    !     gaussian%of_r = gaussian%of_r * (boundary%electrons%density%of_r + &
    !                                      boundary%ions%density%of_r)
    !     dfield_drho%of_r = dfield_drho%of_r + gaussian%of_r
    !     !
    !     CALL destroy_environ_density(gaussian)
    !     CALL destroy_environ_gradient(gradaux)
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE delectronic_field_drho
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE field_aware_de_drho(boundary, de_dboundary, de_drho)
    !     !--------------------------------------------------------------------------------
    !     !
    !     USE tools_functions, ONLY: derivative_of_functions
    !     USE core_fft, ONLY: field_of_gradrho
    !     !
    !     IMPLICIT NONE
    !     !
    !     TYPE(environ_boundary), INTENT(IN) :: boundary
    !     TYPE(environ_density), INTENT(IN) :: de_dboundary
    !     TYPE(environ_density), INTENT(INOUT) :: de_drho
    !     !
    !     TYPE(environ_cell), POINTER :: cell
    !     TYPE(fft_core), POINTER :: fft
    !     INTEGER, POINTER :: nsoft_spheres, nnr
    !     !
    !     INTEGER :: i, j, ipol
    !     REAL(DP) :: df
    !     CHARACTER(LEN=80) :: sub_name = 'field_aware_de_drho'
    !     CHARACTER(LEN=100) :: l_daux = 'daux'
    !     CHARACTER(LEN=100) :: l_aux = 'aux'
    !     !
    !     TYPE(environ_density), ALLOCATABLE :: local(:)
    !     !
    !     TYPE(environ_density) :: aux, daux, rhotot
    !     TYPE(environ_gradient) :: gradaux
    !     !
    !     cell => de_drho%cell
    !     fft => boundary%core%fft
    !     nnr => cell%nnr
    !     !
    !     IF (boundary%mode .EQ. 'fa-ionic') THEN
    !         nsoft_spheres => boundary%ions%number
    !         !
    !         IF (nsoft_spheres .LE. 0) &
    !             CALL errore(sub_name, 'Inconsistent number of soft-spheres', 1)
    !         !
    !         ALLOCATE (local(nsoft_spheres))
    !         !
    !         !----------------------------------------------------------------------------
    !         ! Compute field-independent soft spheres and gradients
    !         !
    !         DO i = 1, nsoft_spheres
    !             CALL init_environ_density(cell, local(i))
    !             CALL density_of_functions(boundary%soft_spheres(i), local(i), .FALSE.)
    !         END DO
    !         !
    !         CALL init_environ_density(cell, aux)
    !         !
    !         DO i = 1, nsoft_spheres
    !             CALL derivative_of_functions(boundary%soft_spheres(i), aux, .TRUE.)
    !             !
    !             DO j = 1, nsoft_spheres
    !                 IF (j .EQ. i) CYCLE
    !                 aux%of_r = aux%of_r * local(j)%of_r
    !             END DO
    !             !
    !             df = dscaling_of_field(boundary%field_factor, &
    !                                    boundary%charge_asymmetry, &
    !                                    boundary%field_max, &
    !                                    boundary%field_min, &
    !                                    boundary%ion_field(i))
    !             !
    !             df = df * boundary%ions%iontype(boundary%ions%ityp(i))%solvationrad * &
    !                  boundary%alpha * scalar_product_environ_density(aux, de_dboundary)
    !             !
    !             de_drho%of_r = de_drho%of_r + boundary%dion_field_drho(i)%of_r * df
    !         END DO
    !         !
    !         CALL destroy_environ_density(aux)
    !         !
    !         DO i = 1, nsoft_spheres
    !             CALL destroy_environ_density(local(i))
    !         END DO
    !         !
    !         DEALLOCATE (local)
    !     ELSE IF (boundary%mode .EQ. 'fa-electronic') THEN
    !         !
    !         ! Try loading in components that need to be reused, use aux for the field aware function, which
    !         ! is the same as the ionic function (scaling_of_field), only taken as an exponent. daux is the
    !         ! result from dscaling_of_field, both of these components are used twice.
    !         !
    !         ! Note that this implementation is not speed optimal due to the fact that many of the components
    !         ! are calculated on the fly and stored temporarily in order to save global storage.
    !         !
    !         CALL init_environ_density(cell, aux)
    !         CALL init_environ_density(cell, daux)
    !         !
    !         CALL scaling_of_density(boundary%field_factor, boundary%charge_asymmetry, &
    !                                 boundary%field_max, boundary%field_min, &
    !                                 boundary%normal_field, aux)
    !         !
    !         ! COMMENTED TO ISOLATE ANOTHER TERM
    !         CALL dscaling_of_density(boundary%field_factor, boundary%charge_asymmetry, &
    !                                  boundary%field_max, boundary%field_min, &
    !                                  boundary%normal_field, daux)
    !         !
    !         ! CALL write_cube(boundary%normal_field)
    !         ! CALL write_cube(daux, label=l_daux)
    !         ! CALL write_cube(aux, label=l_aux)
    !         !
    !         !----------------------------------------------------------------------------
    !         ! compute the green's function part first
    !         !
    !         CALL init_environ_density(cell, rhotot)
    !         CALL init_environ_gradient(cell, gradaux)
    !         !
    !         ! assume boundary_of_density has already been called, this contains
    !         ! the sccs !calls that create dscaled ( to be used here... )
    !         !
    !         CALL external_gradient(boundary%electrons%density%of_r, gradaux%of_r)
    !         !
    !         ! COMMENTED THIS TERM TO ISOLATE ANOTHER TERM
    !         !
    !         DO ipol = 1, 3
    !             gradaux%of_r(ipol, :) = gradaux%of_r(ipol, :) * boundary%electrons%density%of_r(:)
    !             gradaux%of_r(ipol, :) = gradaux%of_r(ipol, :) * boundary%dscaled%of_r(:)
    !             gradaux%of_r(ipol, :) = gradaux%of_r(ipol, :) * aux%of_r(:)
    !             gradaux%of_r(ipol, :) = gradaux%of_r(ipol, :) * daux%of_r(:)
    !             gradaux%of_r(ipol, :) = gradaux%of_r(ipol, :) * de_dboundary%of_r(:)
    !         END DO
    !         !
    !         CALL field_of_gradrho(gradaux%of_r, de_drho%of_r, fft)
    !         !
    !         ! add remaining terms, using daux to update new term temporarily
    !         ! COMMENTED TO ISOLATE ANOTHER TERM
    !         !
    !         ! rhotot%of_r = boundary%electrons%density%of_r + &
    !         !               boundary%ions%density%of_r
    !         !
    !         ! daux%of_r = daux%of_r * boundary%electrons%density%of_r
    !         ! daux%of_r = daux%of_r * rhotot%of_r
    !         ! daux%of_r = daux%of_r * (-2.D0) * tpi
    !         ! daux%of_r = daux%of_r + 1.D0
    !         ! daux%of_r = daux%of_r * aux%of_r
    !         ! FOLLOWING TERM IS TEMPORARY
    !         ! daux%of_r = aux%of_r
    !         ! THIS STUFF IS NOT TEMPORARY
    !         ! daux%of_r = daux%of_r * boundary%dscaled%of_r
    !         ! daux%of_r = daux%of_r * de_dboundary%of_r
    !         ! de_drho%of_r = daux%of_r - de_drho%of_r
    !         ! FOLLOWING TERM IS TEMPORARY
    !         ! de_drho%of_r = daux%of_r
    !         ! deallocate
    !         !
    !         CALL destroy_environ_density(aux)
    !         CALL destroy_environ_density(daux)
    !         CALL destroy_environ_density(rhotot)
    !         CALL destroy_environ_gradient(gradaux)
    !         !
    !     END IF
    !     !
    !     RETURN
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE field_aware_de_drho
    ! !------------------------------------------------------------------------------------
    ! !>
    ! !!
    ! !------------------------------------------------------------------------------------
    ! SUBROUTINE field_aware_dboundary_dions(index, boundary, partial)
    !     !--------------------------------------------------------------------------------
    !     !
    !     USE tools_functions, ONLY: derivative_of_functions
    !     USE core_fft, ONLY: hessv_h_of_rho_r
    !     !
    !     IMPLICIT NONE
    !     !
    !     INTEGER, INTENT(IN) :: index
    !     TYPE(environ_boundary), INTENT(IN) :: boundary
    !     TYPE(environ_gradient), INTENT(INOUT) :: partial
    !     !
    !     TYPE(environ_cell), POINTER :: cell
    !     TYPE(fft_core), POINTER :: fft
    !     INTEGER, POINTER :: nsoft_spheres, nnr
    !     !
    !     INTEGER :: i, j, ipol
    !     REAL(DP) :: df
    !     CHARACTER(LEN=80) :: sub_name = 'field_aware_dboundary_dions'
    !     CHARACTER(len=100) :: strh1 = 'h1'
    !     CHARACTER(len=100) :: strh2 = 'h2'
    !     CHARACTER(len=100) :: strdh1 = 'dh1'
    !     CHARACTER(len=100) :: strdh2 = 'dh2'
    !     !
    !     TYPE(environ_density), ALLOCATABLE :: local(:)
    !     !
    !     TYPE(environ_density) :: aux, daux
    !     TYPE(environ_gradient) :: gradaux, gradlocal
    !     TYPE(environ_hessian) :: hesslocal
    !     !
    !     cell => boundary%scaled%cell
    !     fft => boundary%core%fft
    !     nnr => cell%nnr
    !     !
    !     IF (ionode) &
    !         WRITE (program_unit, '(1X,a)') 'in function field_aware_dboundary_dions'
    !     !
    !     IF (boundary%mode .EQ. 'fa-ionic') THEN
    !         nsoft_spheres => boundary%ions%number
    !         !
    !         IF (nsoft_spheres .LE. 0) &
    !             CALL errore(sub_name, 'Inconsistent number of soft-spheres', 1)
    !         !
    !         ALLOCATE (local(nsoft_spheres))
    !         !
    !         !----------------------------------------------------------------------------
    !         ! Compute field-independent soft spheres and gradients
    !         !
    !         DO i = 1, nsoft_spheres
    !             CALL init_environ_density(cell, local(i))
    !             CALL density_of_functions(boundary%soft_spheres(i), local(i), .FALSE.)
    !         END DO
    !         !
    !         ! CALL write_cube(local(1), label=strh1)
    !         ! CALL write_cube(local(2), label=strh2)
    !         !
    !         CALL init_environ_density(cell, aux)
    !         CALL init_environ_gradient(cell, gradaux)
    !         !
    !         DO i = 1, nsoft_spheres
    !             CALL derivative_of_functions(boundary%soft_spheres(i), aux, .TRUE.)
    !             ! IF (i .EQ. 1) CALL write_cube(aux, label=strdh1)
    !             ! IF (i .EQ. 2) CALL write_cube(aux, label=strdh2)
    !             !
    !             DO j = 1, nsoft_spheres
    !                 IF (j .EQ. i) CYCLE
    !                 aux%of_r = aux%of_r * local(j)%of_r
    !             END DO
    !             !
    !             df = dscaling_of_field(boundary%field_factor, &
    !                                    boundary%charge_asymmetry, &
    !                                    boundary%field_max, &
    !                                    boundary%field_min, &
    !                                    boundary%ion_field(i))
    !             !
    !             PRINT *, "df", i, boundary%ion_field(i), df
    !             !
    !             df = df * boundary%ions%iontype(boundary%ions%ityp(i))%solvationrad * &
    !                  boundary%alpha
    !             !
    !             aux%of_r = aux%of_r * df
    !             !
    !             DO ipol = 1, 3
    !                 !
    !                 gradaux%of_r(ipol, :) = &
    !                     gradaux%of_r(ipol, :) + &
    !                     aux%of_r(:) * boundary%partial_of_ion_field(ipol, i, index)
    !                 !
    !             END DO
    !             !
    !         END DO
    !         !
    !         ! STOP
    !         !
    !         partial%of_r = partial%of_r + gradaux%of_r
    !         !
    !         CALL destroy_environ_gradient(gradaux)
    !         !
    !         CALL destroy_environ_density(aux)
    !         !
    !         DO i = 1, nsoft_spheres
    !             CALL destroy_environ_density(local(i))
    !         END DO
    !         !
    !         DEALLOCATE (local)
    !     ELSE IF (boundary%mode .EQ. 'fa-electronic') THEN
    !         !
    !         ! Some of these terms may be able to be calculated outside of this function, which is repeated per
    !         ! atom. These terms do not depend on atom, but do take up space so are currently calculated on the
    !         ! fly. It may be preferable to change this (since it's only 1 environ density object...)
    !         !
    !         CALL init_environ_density(cell, aux)
    !         CALL init_environ_density(cell, daux)
    !         !
    !         CALL init_environ_gradient(cell, gradaux)
    !         CALL init_environ_gradient(cell, gradlocal)
    !         !
    !         CALL init_environ_hessian(cell, hesslocal)
    !         !
    !         CALL scaling_of_density(boundary%field_factor, boundary%charge_asymmetry, &
    !                                 boundary%field_max, boundary%field_min, &
    !                                 boundary%normal_field, aux)
    !         !
    !         ! CALL write_cube(aux, label=stra)
    !         !
    !         CALL dscaling_of_density(boundary%field_factor, boundary%charge_asymmetry, &
    !                                  boundary%field_max, boundary%field_min, &
    !                                  boundary%normal_field, daux)

    !         ! CALL write_cube(boundary%electrons%density, label=strd)
    !         !
    !         daux%of_r = daux%of_r * boundary%electrons%density%of_r
    !         ! CALL write_cube(daux, label=strd)
    !         daux%of_r = daux%of_r * aux%of_r
    !         daux%of_r = daux%of_r * boundary%dscaled%of_r
    !         ! CALL write_cube(boundary%dscaled, label=stra)
    !         !
    !         ! STOP
    !         !
    !         !----------------------------------------------------------------------------
    !         ! now calculate gradaux, which will build the partial normal field,
    !         ! and then the partial interface
    !         !
    !         CALL density_of_functions(boundary%ions%smeared_ions(index), aux, .TRUE.)
    !         CALL hessv_h_of_rho_r(aux%of_r, hesslocal%of_r, fft)
    !         CALL external_gradient(boundary%electrons%density%of_r, gradlocal%of_r)
    !         CALL scalar_product_environ_hessian(hesslocal, gradlocal, gradaux)
    !         !
    !         DO ipol = 1, 3
    !             gradaux%of_r(ipol, :) = gradaux%of_r(ipol, :) * daux%of_r(:)
    !         END DO
    !         !
    !         partial%of_r = gradaux%of_r
    !         !
    !         CALL destroy_environ_density(aux)
    !         CALL destroy_environ_density(daux)
    !         !
    !         CALL destroy_environ_gradient(gradaux)
    !         CALL destroy_environ_gradient(gradlocal)
    !         !
    !         CALL destroy_environ_hessian(hesslocal)
    !         !
    !     END IF
    !     !
    !     RETURN
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE field_aware_dboundary_dions
    ! !------------------------------------------------------------------------------------
    ! !
    ! !------------------------------------------------------------------------------------
END MODULE boundary_field_aware
!----------------------------------------------------------------------------------------
