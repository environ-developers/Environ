! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
! Copyright (C) 2006-2010 Quantum ESPRESSO group
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
!! Module containing the main drivers to compute Environ contributions
!! to Kohn-Sham potential, total energy and inter-atomic forces
!!
!----------------------------------------------------------------------------------------
MODULE environ_main
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP, e2, fpi
    !
    USE physical_types, ONLY: environ_charges, environ_electrons
    USE representation_types, ONLY: environ_density, environ_gradient
    !
    USE environ_base, ONLY: vzero, solvent, lelectrostatic, velectrostatic, &
                            vreference, dvtot, lsoftcavity, vsoftcavity, lsurface, &
                            env_surface_tension, lvolume, env_pressure, lconfine, &
                            env_confine, vconfine, lstatic, static, lexternals, &
                            lelectrolyte, electrolyte, lsoftsolvent, lsoftelectrolyte, &
                            system_cell, environment_cell, system_charges, &
                            environment_charges, mapping, system_electrons, &
                            environment_electrons, niter, lrigidcavity, lrigidsolvent, &
                            lrigidelectrolyte, loptical, optical, &
                            system_response_charges, environment_response_charges
    !
    USE electrostatic_base, ONLY: reference, outer
    !
    USE utils_density, ONLY: init_environ_density, destroy_environ_density
    USE utils_gradient, ONLY: init_environ_gradient, destroy_environ_gradient
    ! USE utils_charges, ONLY: update_environ_charges ! #TODO keep for now until external tests are fully debugged
    USE utils_electrons
    !
    USE tools_charges, ONLY: charges_of_potential
    USE tools_dielectric, ONLY: calc_dedielectric_dboundary, calc_dvdielectric_dboundary
    USE tools_electrolyte, ONLY: calc_deelectrolyte_dboundary, calc_eelectrolyte
    !
    USE tools_math, ONLY: scalar_product_environ_density, &
                          scalar_product_environ_gradient_density
    !
    USE tools_mapping, ONLY: map_large_to_small
    !
    USE tools_generate_boundary, ONLY: solvent_aware_de_dboundary, calc_dboundary_dions
    !
    USE embedding_electrostatic, ONLY: calc_velectrostatic, calc_eelectrostatic, &
                                       calc_felectrostatic
    !
    USE embedding_confine, ONLY: calc_vconfine, calc_deconfine_dboundary
    USE embedding_surface, ONLY: calc_desurface_dboundary, calc_esurface
    USE embedding_volume, ONLY: calc_devolume_dboundary, calc_evolume
    !
    USE environ_output, ONLY: print_environ_density, print_environ_charges
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the Environ contribution to the local potential. All
    !! the Environ modules need to be called here. The potentials are
    !! all computed on the dense real-space grid and added to vtot.
    !
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_venviron(update, nnr, vtot, local_verbose)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: update
        INTEGER, INTENT(IN) :: nnr
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        !
        REAL(DP), INTENT(OUT) :: vtot(nnr)
        !
        TYPE(environ_density) :: aux
        TYPE(environ_density) :: de_dboundary
        !
        !--------------------------------------------------------------------------------
        ! If not updating the potentials, add old potentials and exit
        !
        vtot = vzero%of_r
        !
        IF (.NOT. update) THEN
            vtot = vtot + dvtot%of_r
            !
            IF (PRESENT(local_verbose)) THEN
                !
                CALL print_environ_density(dvtot, local_verbose)
                !
                IF (lelectrostatic) THEN
                    !
                    CALL print_environ_density(vreference, local_verbose)
                    !
                    CALL print_environ_density(velectrostatic, local_verbose)
                    !
                    CALL print_environ_charges(system_charges, local_verbose, &
                                               local_depth=0)
                END IF
                !
                IF (lconfine) &
                    CALL print_environ_density(vconfine, local_verbose)
                !
                IF (lsoftcavity) &
                    CALL print_environ_density(vsoftcavity, local_verbose)
                !
            END IF
            !
            RETURN
            !
        END IF
        !
        dvtot%of_r = 0.D0
        !
        CALL init_environ_density(system_cell, aux)
        !
        !--------------------------------------------------------------------------------
        ! If any form of electrostatic embedding is present, calculate its contribution
        !
        IF (lelectrostatic) THEN
            !
            !----------------------------------------------------------------------------
            ! Electrostatics is also computed inside the calling program,
            ! need to remove the reference #TODO to-be-decided
            !
            CALL calc_velectrostatic(reference, system_charges, vreference)
            !
            CALL print_environ_density(vreference)
            !
            CALL calc_velectrostatic(outer, environment_charges, velectrostatic)
            !
            ! IF (lexternals) CALL update_environ_charges(environment_charges, lexternals)
            ! #TODO keep for now until external tests are fully debugged
            !
            CALL print_environ_density(velectrostatic)
            !
            CALL map_large_to_small(mapping, velectrostatic, aux)
            !
            dvtot%of_r = aux%of_r - vreference%of_r
            !
            CALL charges_of_potential(velectrostatic, environment_charges)
            !
            CALL print_environ_charges(environment_charges)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! #TODO add brief description of confine
        !
        IF (lconfine) THEN
            !
            CALL calc_vconfine(env_confine, solvent, vconfine)
            !
            CALL print_environ_density(vconfine)
            !
            CALL map_large_to_small(mapping, vconfine, aux)
            !
            dvtot%of_r = dvtot%of_r + aux%of_r
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Compute the total potential depending on the boundary
        !
        IF (lsoftcavity) THEN
            vsoftcavity%of_r = 0.D0
            !
            CALL init_environ_density(environment_cell, de_dboundary)
            !
            IF (lsoftsolvent) THEN
                de_dboundary%of_r = 0.D0
                !
                ! if surface tension greater than zero, calculate cavity contribution
                IF (lsurface) &
                    CALL calc_desurface_dboundary(env_surface_tension, solvent, &
                                                  de_dboundary)
                !
                ! if external pressure different from zero, calculate PV contribution
                IF (lvolume) &
                    CALL calc_devolume_dboundary(env_pressure, solvent, de_dboundary)
                !
                ! if confinement potential not zero, calculate confine contribution
                IF (lconfine) &
                    CALL calc_deconfine_dboundary(env_confine, &
                                                  environment_charges%electrons%density, &
                                                  de_dboundary)
                !
                ! if dielectric embedding, calcultes dielectric contribution
                IF (lstatic) &
                    CALL calc_dedielectric_dboundary(static, velectrostatic, &
                                                     de_dboundary)
                !
                ! if solvent-aware interface correct the potential
                IF (solvent%solvent_aware) &
                    CALL solvent_aware_de_dboundary(solvent, de_dboundary)
                !
                IF (solvent%field_aware) THEN
                    !
                    ! CALL field_aware_de_drho(solvent, de_dboundary, vsoftcavity) ! #TODO field-aware
                    ! if field-aware interface use a more cumbersome formula
                    !
                    CALL errore('field-aware1', 'Option not yet implimented ', 1)
                    !
                ELSE
                    !
                    vsoftcavity%of_r = de_dboundary%of_r * solvent%dscaled%of_r
                    ! multiply by derivative of the boundary w.r.t electronic density
                    !
                END IF
                !
            END IF
            !
            IF (lsoftelectrolyte) THEN
                de_dboundary%of_r = 0.D0
                !
                CALL calc_deelectrolyte_dboundary(electrolyte, de_dboundary)
                ! if electrolyte is present add its non-electrostatic contribution
                !
                ! if solvent-aware interface correct the potential
                IF (electrolyte%boundary%solvent_aware) &
                    CALL solvent_aware_de_dboundary(electrolyte%boundary, de_dboundary)
                !
                IF (electrolyte%boundary%field_aware) THEN
                    !
                    !--------------------------------------------------------------------
                    ! CALL field_aware_de_drho(electrolyte%boundary, de_dboundary, vsoftcavity) ! #TODO field-aware
                    ! if field-aware, correct the derivative of the interface function
                    !
                    CALL errore('field-aware2', 'Option not yet implimented ', 1)
                    !
                ELSE
                    !
                    ! multiply for the derivative of the boundary w.r.t electronic density
                    vsoftcavity%of_r = &
                        vsoftcavity%of_r + &
                        de_dboundary%of_r * electrolyte%boundary%dscaled%of_r
                    !
                END IF
                !
            END IF
            !
            CALL print_environ_density(vsoftcavity)
            !
            CALL map_large_to_small(mapping, vsoftcavity, aux)
            !
            dvtot%of_r = dvtot%of_r + aux%of_r
            !
            CALL destroy_environ_density(de_dboundary)
            !
        END IF
        !
        CALL destroy_environ_density(aux)
        !
        CALL print_environ_density(dvtot, local_verbose)
        !
        vtot = vtot + dvtot%of_r
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_venviron
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the Environ contribution to the energy. We must remove
    !! int v_environ * rhoelec that is automatically included in the
    !! energy computed as the sum of Kohn-Sham eigenvalues.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_eenviron(deenviron, eelectrostatic, esurface, evolume, econfine, &
                             eelectrolyte)
        !--------------------------------------------------------------------------------
        !
        ! USE embedding_confine, ONLY: calc_econfine #TODO to-be-decided
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(OUT) :: deenviron, eelectrostatic, esurface, &
                                 evolume, econfine, eelectrolyte
        !
        REAL(DP) :: ereference
        !
        !--------------------------------------------------------------------------------
        !
        eelectrostatic = 0.D0
        esurface = 0.D0
        evolume = 0.D0
        econfine = 0.D0
        eelectrolyte = 0.D0
        !
        niter = niter + 1
        !
        deenviron = -scalar_product_environ_density(system_electrons%density, dvtot)
        ! calculates the energy corrections
        !
        !--------------------------------------------------------------------------------
        ! If electrostatic is on, compute electrostatic energy
        !
        IF (lelectrostatic) THEN
            !
            CALL calc_eelectrostatic(reference%core, system_charges, vreference, &
                                     ereference)
            !
            CALL calc_eelectrostatic(outer%core, environment_charges, velectrostatic, &
                                     eelectrostatic)
            !
            eelectrostatic = eelectrostatic - ereference
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        IF (lsurface) CALL calc_esurface(env_surface_tension, solvent, esurface)
        ! if surface tension not zero, compute cavitation energy
        !
        IF (lvolume) CALL calc_evolume(env_pressure, solvent, evolume)
        ! if pressure not zero, compute PV energy
        !
        ! if confinement potential not zero compute confine energy
        IF (lconfine) &
            econfine = scalar_product_environ_density(environment_electrons%density, &
                                                      vconfine)
        !
        IF (lelectrolyte) CALL calc_eelectrolyte(electrolyte, eelectrolyte)
        ! if electrolyte is present, calculate its non-electrostatic contribution
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_eenviron
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the Environ contribution to the forces. Due to
    !! Hellman-Feynman only a few of the Environ modules have an
    !! effect on the atomic forces.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_fenviron(nat, force_environ)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        !
        REAL(DP), INTENT(INOUT) :: force_environ(3, nat)
        !
        INTEGER :: i
        TYPE(environ_density) :: de_dboundary
        TYPE(environ_gradient) :: partial
        !
        !--------------------------------------------------------------------------------
        !
        force_environ = 0.D0
        !
        ! compute the electrostatic embedding contribution to the interatomic forces
        IF (lelectrostatic) &
            CALL calc_felectrostatic(outer, nat, environment_charges, force_environ)
        !
        !--------------------------------------------------------------------------------
        ! Compute the total forces depending on the boundary
        !
        IF (lrigidcavity) THEN
            !
            CALL init_environ_density(environment_cell, de_dboundary)
            !
            CALL init_environ_gradient(environment_cell, partial)
            !
            IF (lrigidsolvent) THEN
                de_dboundary%of_r = 0.D0
                !
                ! if surface tension greater than zero, calculate cavity contribution
                IF (lsurface) &
                    CALL calc_desurface_dboundary(env_surface_tension, &
                                                  solvent, de_dboundary)
                !
                ! if external pressure not zero, calculate PV contribution
                IF (lvolume) &
                    CALL calc_devolume_dboundary(env_pressure, solvent, de_dboundary)
                !
                ! if confinement potential not zero, calculate confine contribution
                IF (lconfine) &
                    CALL calc_deconfine_dboundary(env_confine, &
                                                  environment_charges%electrons%density, &
                                                  de_dboundary)
                !
                ! if dielectric embedding, calculate dielectric contribution
                IF (lstatic) &
                    CALL calc_dedielectric_dboundary(static, velectrostatic, &
                                                     de_dboundary)
                !
                ! if solvent-aware, correct the potential
                IF (solvent%solvent_aware) &
                    CALL solvent_aware_de_dboundary(solvent, de_dboundary)
                !
                ! if field-aware, compute partial derivatives of field fluxes w.r.t ionic positions ! #TODO field-aware
                ! IF (solvent%mode .EQ. 'fa-ionic') &
                !     CALL compute_ion_field_partial(solvent%ions%number, &
                !                                    solvent%soft_spheres, &
                !                                    solvent%ions, solvent%electrons, &
                !                                    solvent%ion_field, &
                !                                    solvent%partial_of_ion_field, &
                !                                    solvent%core%fft)
                !
                !------------------------------------------------------------------------
                ! Multiply by derivative of the boundary w.r.t ionic positions
                !
                DO i = 1, nat
                    !
                    CALL calc_dboundary_dions(i, solvent, partial)
                    !
                    ! if field-aware, correct the derivative of the interface function
                    ! IF (solvent%field_aware) &
                    !     CALL field_aware_dboundary_dions(i, solvent, partial) ! #TODO field-aware
                    !
                    force_environ(:, i) = &
                        force_environ(:, i) - &
                        scalar_product_environ_gradient_density(partial, de_dboundary)
                    !
                END DO
                !
            END IF
            !
            IF (lrigidelectrolyte) THEN
                de_dboundary%of_r = 0.D0
                !
                ! if electrolyte is present, add its non-electrostatic contribution
                CALL calc_deelectrolyte_dboundary(electrolyte, de_dboundary)
                !
                ! if solvent-aware, correct the potential
                IF (electrolyte%boundary%solvent_aware) &
                    CALL solvent_aware_de_dboundary(electrolyte%boundary, de_dboundary)
                !
                ! if field-aware, compute partial derivatives of field fluxes w.r.t ionic positions ! #TODO field-aware
                ! IF (electrolyte%boundary%mode .EQ. 'fa-ionic') &
                !     CALL compute_ion_field_partial(electrolyte%boundary%ions%number, &
                !                                    electrolyte%boundary%soft_spheres, &
                !                                    electrolyte%boundary%ions, &
                !                                    electrolyte%boundary%electrons, &
                !                                    electrolyte%boundary%ion_field, &
                !                                    electrolyte%boundary%partial_of_ion_field, &
                !                                    electrolyte%boundary%core%fft)
                !
                !------------------------------------------------------------------------
                ! Multiply by derivative of the boundary w.r.t ionic positions
                !
                DO i = 1, nat
                    !
                    CALL calc_dboundary_dions(i, electrolyte%boundary, partial)
                    !
                    ! if field-aware interface, correct the derivative of the interface function ! #TODO field-aware
                    ! IF (electrolyte%boundary%field_aware) &
                    !     CALL field_aware_dboundary_dions(i, electrolyte%boundary, &
                    !                                      partial)
                    !
                    force_environ(:, i) = &
                        force_environ(:, i) - &
                        scalar_product_environ_gradient_density(partial, de_dboundary)
                    !
                END DO
                !
            END IF
            !
            CALL destroy_environ_gradient(partial)
            !
            CALL destroy_environ_density(de_dboundary)
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_fenviron
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the Environ contribution to the response potential in TD calculations
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dvenviron(nnr, dvtot)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(OUT) :: dvtot(nnr)
        !
        TYPE(environ_density) :: aux
        TYPE(environ_density) :: dvreference
        TYPE(environ_density) :: dvelectrostatic
        TYPE(environ_density) :: dvsoftcavity
        TYPE(environ_density) :: dv_dboundary
        !
        CALL init_environ_density(system_cell, aux)
        !
        !--------------------------------------------------------------------------------
        ! If any form of electrostatic embedding is present, calculate its contribution
        !
        IF (lelectrostatic) THEN
            !
            !----------------------------------------------------------------------------
            ! Electrostatics is also computed inside the calling program,
            ! need to remove the reference
            !
            CALL init_environ_density(system_cell, dvreference)
            !
            CALL calc_velectrostatic(reference, system_response_charges, dvreference)
            !
            CALL init_environ_density(environment_cell, dvelectrostatic)
            !
            CALL calc_velectrostatic(outer, environment_response_charges, dvelectrostatic)
            !
            CALL map_large_to_small(mapping, dvelectrostatic, aux)
            !
            dvtot(:) = dvtot(:) + aux%of_r(:) - dvreference%of_r(:)
            !
            CALL destroy_environ_density(dvreference)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Compute the response potential depending on the boundary
        !
        IF (lsoftcavity) THEN
            !
            CALL init_environ_density(environment_cell, dvsoftcavity)
            !
            CALL init_environ_density(environment_cell, dv_dboundary)
            !
            IF (lsoftsolvent) THEN
                dv_dboundary%of_r = 0.D0
                !
                ! if dielectric embedding, calcultes dielectric contribution
                IF (loptical) &
                     CALL calc_dvdielectric_dboundary(optical, velectrostatic, &
                     dvelectrostatic, dv_dboundary)
                !
                dvsoftcavity%of_r = dv_dboundary%of_r * solvent%dscaled%of_r
                !
            END IF
            !
            CALL map_large_to_small(mapping, dvsoftcavity, aux)
            !
            dvtot(:) = dvtot(:) + aux%of_r(:)
            !
            CALL destroy_environ_density(dv_dboundary)
            !
            CALL destroy_environ_density(dvsoftcavity)
            !
        END IF
        !
        IF (lelectrostatic) CALL destroy_environ_density(dvelectrostatic)
        !
        CALL destroy_environ_density(aux)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dvenviron
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE environ_main
!----------------------------------------------------------------------------------------
