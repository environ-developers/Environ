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
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Ismaila Dabo       (DMSE, Penn State)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------------------
!>
!! Module containing the main routines to handle environ_dielectric
!! derived data types.
!!
!! Environ_dielectric is the type to store the details of the dielectric
!! embedding. It contains the specifics of externally-defined dielectric
!! regions and it links the boundary details. Starting from these quantities,
!! It builds the dielectric function in space (stored in %epsilon component)
!! and the factors derived from it that are required by the generalized
!! Poisson solver (gradient of the logarithm, sqrt factor needed by
!! preconditioned conjugate gradient, etc.).
!!
!----------------------------------------------------------------------------------------
MODULE tools_dielectric
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: e2, fpi
    !
    USE physical_types, ONLY: environ_dielectric
    USE representation_types, ONLY: environ_density, environ_gradient
    USE cell_types, ONLY: environ_cell
    !
    USE utils_density, ONLY: init_environ_density, destroy_environ_density
    !
    USE utils_gradient, ONLY: init_environ_gradient, update_gradient_modulus, &
                              destroy_environ_gradient
    !
    USE tools_math, ONLY: scalar_product_environ_gradient, integrate_environ_density
    USE core_fft, ONLY: gradient_fft
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: dielectric_of_potential, calc_dedielectric_dboundary, &
              calc_dvdielectric_dboundary
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE dielectric_of_potential(charges, potential, dielectric)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: charges, potential
        !
        TYPE(environ_dielectric), INTENT(INOUT) :: dielectric
        !
        TYPE(environ_cell), POINTER :: cell
        !
        TYPE(environ_gradient) :: gradient
        !
        CHARACTER(LEN=80) :: sub_name = 'dielectric_of_potential'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(potential%cell, charges%cell)) &
            CALL env_errore(sub_name, &
                            'Mismatch in domains of potential and charges', 1)
        !
        IF (.NOT. ASSOCIATED(potential%cell, dielectric%density%cell)) &
            CALL env_errore(sub_name, &
                            'Mismatch in domains of potential and dielectric', 1)
        !
        cell => charges%cell
        !
        CALL init_environ_gradient(cell, gradient)
        !
        CALL gradient_fft(dielectric%boundary%core%fft, potential, gradient)
        !
        CALL scalar_product_environ_gradient(dielectric%gradlog, gradient, &
                                             dielectric%density)
        !
        dielectric%density%of_r = dielectric%density%of_r / fpi / e2 + charges%of_r * &
                                  (1.D0 - dielectric%epsilon%of_r) / &
                                  dielectric%epsilon%of_r
        !
        CALL destroy_environ_gradient(gradient)
        !
        dielectric%charge = integrate_environ_density(dielectric%density)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE dielectric_of_potential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dedielectric_dboundary(dielectric, velectrostatic, de_dboundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_dielectric), INTENT(IN) :: dielectric
        TYPE(environ_density), INTENT(IN) :: velectrostatic
        !
        TYPE(environ_density), INTENT(INOUT) :: de_dboundary
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_gradient) :: gradient
        !
        !--------------------------------------------------------------------------------
        !
        cell => de_dboundary%cell
        !
        CALL init_environ_gradient(cell, gradient)
        !
        CALL gradient_fft(dielectric%boundary%core%fft, velectrostatic, gradient)
        !
        CALL update_gradient_modulus(gradient)
        !
        de_dboundary%of_r = de_dboundary%of_r - &
                            gradient%modulus%of_r**2 * dielectric%depsilon%of_r * &
                            0.5D0 / fpi / e2
        !
        CALL destroy_environ_gradient(gradient)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dedielectric_dboundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dvdielectric_dboundary(dielectric, velectrostatic, &
                                           dvelectrostatic, dv_dboundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_dielectric), INTENT(IN) :: dielectric
        TYPE(environ_density), INTENT(IN) :: velectrostatic
        TYPE(environ_density), INTENT(IN) :: dvelectrostatic
        !
        TYPE(environ_density), INTENT(INOUT) :: dv_dboundary
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density) :: aux
        TYPE(environ_gradient) :: gradient, dgradient
        !
        !--------------------------------------------------------------------------------
        !
        cell => dv_dboundary%cell
        !
        CALL init_environ_gradient(cell, gradient)
        !
        CALL gradient_fft(dielectric%boundary%core%fft, velectrostatic, gradient)
        !
        CALL init_environ_gradient(cell, dgradient)
        !
        CALL gradient_fft(dielectric%boundary%core%fft, dvelectrostatic, dgradient)
        !
        CALL init_environ_density(cell, aux)
        !
        CALL scalar_product_environ_gradient(gradient, dgradient, aux)
        !
        CALL destroy_environ_gradient(gradient)
        !
        CALL destroy_environ_gradient(dgradient)
        !
        dv_dboundary%of_r = dv_dboundary%of_r - &
                            aux%of_r * dielectric%depsilon%of_r / (fpi * e2)
        !
        CALL destroy_environ_density(aux)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dvdielectric_dboundary
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE tools_dielectric
!----------------------------------------------------------------------------------------
