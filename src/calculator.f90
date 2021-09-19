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
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Ismaila Dabo       (DMSE, Penn State)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_calculator
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
    USE class_environ
    USE class_gradient
    USE class_mapping
    !
    USE env_write_cube
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
    TYPE, PUBLIC :: environ_calculator
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, NOPASS :: force => calc_fenviron
        PROCEDURE, NOPASS :: energy => calc_eenviron
        PROCEDURE, NOPASS :: denergy => calc_deenviron
        PROCEDURE, NOPASS :: potential => calc_venviron
        PROCEDURE, NOPASS :: dpotential => calc_dvenviron
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_calculator
    !------------------------------------------------------------------------------------
    !
    TYPE(environ_calculator), PUBLIC :: calc
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  GENERAL METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the Environ contribution to the local potential. All
    !! the Environ modules need to be called here. The potentials are
    !! all computed on the dense real-space grid and added to vtot.
    !
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_venviron(env, update, local_verbose)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: update
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        !
        CLASS(environ_obj), TARGET, INTENT(INOUT) :: env
        !
        TYPE(environ_density) :: aux
        TYPE(environ_density) :: de_dboundary
        !
        TYPE(environ_mapping), POINTER :: mapping
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_venviron'
        !
        !--------------------------------------------------------------------------------
        !
        mapping => env%setup%mapping
        !
        !--------------------------------------------------------------------------------
        ! If not updating the potentials, add old potentials and exit
        !
        IF (.NOT. update) THEN
            !
            IF (PRESENT(local_verbose)) THEN
                !
                CALL write_cube(env%dvtot, env%system_ions, local_verbose)
                !
                IF (env%setup%lelectrostatic) THEN
                    !
                    CALL write_cube(env%vreference, env%system_ions, local_verbose)
                    !
                    CALL write_cube(env%velectrostatic, env%system_ions, local_verbose)
                    !
                    CALL env%system_charges%printout(local_verbose)
                    !
                END IF
                !
                IF (env%setup%lconfine) &
                    CALL write_cube(env%vconfine, env%system_ions, local_verbose)
                !
                IF (env%setup%lsoftcavity) &
                    CALL write_cube(env%vsoftcavity, env%system_ions, local_verbose)
                !
            END IF
            !
            RETURN
            !
        END IF
        !
        env%dvtot%of_r = 0.D0
        !
        CALL aux%init(env%setup%system_cell)
        !
        !--------------------------------------------------------------------------------
        ! If any form of electrostatic embedding is present, calculate its contribution
        !
        IF (env%setup%lelectrostatic) THEN
            !
            !----------------------------------------------------------------------------
            ! Electrostatics is also computed inside the calling program,
            ! need to remove the reference #TODO to-be-decided
            !
            CALL env%setup%reference%calc_v(env%system_charges, env%vreference)
            !
            CALL write_cube(env%vreference, env%system_ions)
            !
            CALL env%setup%outer%calc_v(env%environment_charges, env%velectrostatic)
            !
            ! IF (this%setup%lexternals) CALL this%environment_charges%update()
            ! #TODO keep for now until external tests are fully debugged
            !
            CALL write_cube(env%velectrostatic, env%system_ions)
            !
            CALL mapping%to_small(env%velectrostatic, aux)
            !
            env%dvtot%of_r = aux%of_r - env%vreference%of_r
            !
            CALL env%environment_charges%of_potential(env%velectrostatic)
            !
            CALL env%environment_charges%printout()
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! #TODO add brief description of confine
        !
        IF (env%setup%lconfine) THEN
            !
            CALL env%solvent%vconfine(env%setup%confine, env%vconfine)
            !
            CALL write_cube(env%vconfine, env%system_ions)
            !
            CALL mapping%to_small(env%vconfine, aux)
            !
            env%dvtot%of_r = env%dvtot%of_r + aux%of_r
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Compute the total potential depending on the boundary
        !
        IF (env%setup%lsoftcavity) THEN
            env%vsoftcavity%of_r = 0.D0
            !
            CALL de_dboundary%init(env%setup%environment_cell)
            !
            IF (env%setup%lsoftsolvent) THEN
                de_dboundary%of_r = 0.D0
                !
                ! if surface tension greater than zero, calculate cavity contribution
                IF (env%setup%lsurface) &
                    CALL env%solvent%desurface_dboundary(env%setup%surface_tension, &
                                                         de_dboundary)
                !
                ! if external pressure different from zero, calculate PV contribution
                IF (env%setup%lvolume) &
                    CALL env%solvent%devolume_dboundary(env%setup%pressure, &
                                                        de_dboundary)
                !
                ! if confinement potential not zero, calculate confine contribution
                IF (env%setup%lconfine) &
                    CALL env%solvent%deconfine_dboundary(env%setup%confine, &
                                                         env%environment_charges%electrons%density, &
                                                         de_dboundary)
                !
                ! if dielectric embedding, calculate dielectric contribution
                IF (env%setup%lstatic) &
                    CALL env%static%de_dboundary(env%velectrostatic, de_dboundary)
                !
                ! if solvent-aware interface correct the potential
                IF (env%solvent%solvent_aware) &
                    CALL env%solvent%sa_de_dboundary(de_dboundary)
                !
                IF (env%solvent%field_aware) THEN
                    CALL io%error(sub_name, 'field-aware not yet implimented ', 1)
                ELSE
                    !
                    env%vsoftcavity%of_r = de_dboundary%of_r * env%solvent%dscaled%of_r
                    ! multiply by derivative of the boundary w.r.t electronic density
                    !
                END IF
                !
            END IF
            !
            IF (env%setup%lsoftelectrolyte) THEN
                de_dboundary%of_r = 0.D0
                !
                CALL env%electrolyte%de_dboundary(de_dboundary)
                ! if electrolyte is present add its non-electrostatic contribution
                !
                ! if solvent-aware interface correct the potential
                IF (env%electrolyte%boundary%solvent_aware) &
                    CALL env%electrolyte%boundary%sa_de_dboundary(de_dboundary)
                !
                IF (env%electrolyte%boundary%field_aware) THEN
                    CALL io%error(sub_name, 'field-aware not yet implimented ', 1)
                ELSE
                    !
                    ! multiply for the derivative of the boundary w.r.t electronic density
                    env%vsoftcavity%of_r = &
                        env%vsoftcavity%of_r + &
                        de_dboundary%of_r * env%electrolyte%boundary%dscaled%of_r
                    !
                END IF
                !
            END IF
            !
            CALL write_cube(env%vsoftcavity, env%system_ions)
            !
            CALL mapping%to_small(env%vsoftcavity, aux)
            !
            env%dvtot%of_r = env%dvtot%of_r + aux%of_r
            !
            CALL de_dboundary%destroy()
            !
        END IF
        !
        CALL aux%destroy()
        !
        CALL write_cube(env%dvtot, env%system_ions, local_verbose)
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
    SUBROUTINE calc_eenviron(env, total_energy)
        !--------------------------------------------------------------------------------
        !
        ! USE embedding_confine, ONLY: calc_econfine #TODO to-be-decided
        !
        IMPLICIT NONE
        !
        CLASS(environ_obj), INTENT(INOUT) :: env
        REAL(DP), INTENT(INOUT) :: total_energy
        !
        REAL(DP) :: ereference
        !
        !--------------------------------------------------------------------------------
        !
        env%eelectrostatic = 0.D0
        env%esurface = 0.D0
        env%evolume = 0.D0
        env%econfine = 0.D0
        env%eelectrolyte = 0.D0
        !
        env%setup%niter = env%setup%niter + 1
        !
        !--------------------------------------------------------------------------------
        ! If electrostatic is on, compute electrostatic energy
        !
        IF (env%setup%lelectrostatic) THEN
            !
            CALL env%setup%reference%calc_e(env%system_charges, env%vreference, &
                                            ereference)
            !
            CALL env%setup%outer%calc_e(env%environment_charges, &
                                        env%velectrostatic, env%eelectrostatic)
            !
            env%eelectrostatic = env%eelectrostatic - ereference
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        ! if surface tension not zero, compute cavitation energy
        IF (env%setup%lsurface) &
            CALL env%solvent%esurface(env%setup%surface_tension, env%esurface)
        !
        IF (env%setup%lvolume) CALL env%solvent%evolume(env%setup%pressure, env%evolume)
        ! if pressure not zero, compute PV energy
        !
        ! if confinement potential not zero compute confine energy
        IF (env%setup%lconfine) &
            env%econfine = env%environment_electrons%density%scalar_product(env%vconfine)
        !
        ! if electrolyte is present, calculate its non-electrostatic contribution
        IF (env%setup%lelectrolyte) CALL env%electrolyte%energy(env%eelectrolyte)
        !
        total_energy = total_energy + env%eelectrostatic + env%esurface + &
                       env%evolume + env%econfine + env%eelectrolyte
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
    SUBROUTINE calc_fenviron(env, nat, force_environ)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        !
        CLASS(environ_obj), TARGET, INTENT(INOUT) :: env
        REAL(DP), INTENT(INOUT) :: force_environ(3, nat)
        !
        INTEGER :: i
        TYPE(environ_density) :: de_dboundary
        TYPE(environ_gradient) :: partial
        !
        TYPE(environ_cell), POINTER :: environment_cell
        !
        REAL(DP), DIMENSION(3, nat) :: freference, felectrostatic
        !
        !--------------------------------------------------------------------------------
        !
        environment_cell => env%setup%environment_cell
        !
        !--------------------------------------------------------------------------------
        ! Compute the electrostatic embedding contribution to the interatomic forces.
        ! If using doublecell, take the difference of contributions from the full charge
        ! density computed in the environment cell (outer solver) and those from the
        ! charge density of ions/electrons computed in the system cell (reference solver)
        !
        force_environ = 0.D0
        freference = 0.D0
        !
        IF (env%setup%lelectrostatic) THEN
            !
            CALL env%setup%outer%calc_f(nat, env%environment_charges, felectrostatic, &
                                        env%setup%ldoublecell)
            !
            IF (env%setup%ldoublecell) THEN
                !
                CALL env%setup%reference%calc_f(nat, env%system_charges, freference, &
                                                env%setup%ldoublecell)
                !
            END IF
            !
            force_environ = felectrostatic - freference
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Compute the total forces depending on the boundary
        !
        IF (env%setup%lrigidcavity) THEN
            !
            CALL de_dboundary%init(environment_cell)
            !
            CALL partial%init(environment_cell)
            !
            IF (env%setup%lrigidsolvent) THEN
                !
                de_dboundary%of_r = 0.D0
                !
                ! if surface tension greater than zero, calculate cavity contribution
                IF (env%setup%lsurface) &
                    CALL env%solvent%desurface_dboundary(env%setup%surface_tension, &
                                                         de_dboundary)
                !
                ! if external pressure not zero, calculate PV contribution
                IF (env%setup%lvolume) &
                    CALL env%solvent%devolume_dboundary(env%setup%pressure, de_dboundary)
                !
                ! if confinement potential not zero, calculate confine contribution
                IF (env%setup%lconfine) &
                    CALL env%solvent%deconfine_dboundary(env%setup%confine, &
                                                         env%environment_charges%electrons%density, &
                                                         de_dboundary)
                !
                ! if dielectric embedding, calculate dielectric contribution
                IF (env%setup%lstatic) &
                    CALL env%static%de_dboundary(env%velectrostatic, de_dboundary)
                !
                ! if solvent-aware, correct the potential
                IF (env%solvent%solvent_aware) &
                    CALL env%solvent%sa_de_dboundary(de_dboundary)
                !
                !------------------------------------------------------------------------
                ! Multiply by derivative of the boundary w.r.t ionic positions
                !
                DO i = 1, nat
                    !
                    CALL env%solvent%dboundary_dions(i, partial)
                    !
                    force_environ(:, i) = force_environ(:, i) - &
                                          partial%scalar_product_density(de_dboundary)
                    !
                END DO
                !
            END IF
            !
            IF (env%setup%lrigidelectrolyte) THEN
                !
                de_dboundary%of_r = 0.D0
                !
                ! if electrolyte is present, add its non-electrostatic contribution
                CALL env%electrolyte%de_dboundary(de_dboundary)
                !
                ! if solvent-aware, correct the potential
                IF (env%electrolyte%boundary%solvent_aware) &
                    CALL env%electrolyte%boundary%sa_de_dboundary(de_dboundary)
                !
                !------------------------------------------------------------------------
                ! Multiply by derivative of the boundary w.r.t ionic positions
                !
                DO i = 1, nat
                    !
                    CALL env%electrolyte%boundary%dboundary_dions(i, partial)
                    !
                    force_environ(:, i) = force_environ(:, i) - &
                                          partial%scalar_product_density(de_dboundary)
                    !
                END DO
                !
            END IF
            !
            CALL partial%destroy()
            !
            CALL de_dboundary%destroy()
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_fenviron
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the Environ contribution to the response potential in TD calculations
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dvenviron(env, nnr, dv)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        !
        CLASS(environ_obj), TARGET, INTENT(INOUT) :: env
        REAL(DP), INTENT(INOUT) :: dv(nnr)
        !
        TYPE(environ_density) :: aux
        TYPE(environ_density) :: dvreference
        TYPE(environ_density) :: dvelectrostatic
        TYPE(environ_density) :: dvsoftcavity
        TYPE(environ_density) :: dv_dboundary
        !
        TYPE(environ_cell), POINTER :: system_cell, environment_cell
        TYPE(environ_mapping), POINTER :: mapping
        !
        !--------------------------------------------------------------------------------
        !
        system_cell => env%setup%system_cell
        environment_cell => env%setup%environment_cell
        mapping => env%setup%mapping
        !
        !--------------------------------------------------------------------------------
        !
        CALL aux%init(system_cell)
        !
        !--------------------------------------------------------------------------------
        ! If any form of electrostatic embedding is present, calculate its contribution
        !
        IF (env%setup%lelectrostatic) THEN
            !
            !----------------------------------------------------------------------------
            ! Electrostatics is also computed inside the calling program,
            ! need to remove the reference
            !
            CALL dvreference%init(system_cell)
            !
            CALL env%setup%reference%calc_v(env%system_response_charges, dvreference)
            !
            CALL dvelectrostatic%init(environment_cell)
            !
            CALL env%setup%outer%calc_v(env%environment_response_charges, &
                                        dvelectrostatic)
            !
            CALL mapping%to_small(dvelectrostatic, aux)
            !
            dv(:) = dv(:) + aux%of_r(:) - dvreference%of_r(:)
            !
            CALL dvreference%destroy()
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Compute the response potential depending on the boundary
        !
        IF (env%setup%lsoftcavity) THEN
            !
            CALL dvsoftcavity%init(environment_cell)
            !
            CALL dv_dboundary%init(environment_cell)
            !
            IF (env%setup%lsoftsolvent) THEN
                dv_dboundary%of_r = 0.D0
                !
                ! if dielectric embedding, calcultes dielectric contribution
                IF (env%setup%loptical) &
                    CALL env%optical%dv_dboundary(env%velectrostatic, &
                                                  dvelectrostatic, dv_dboundary)
                !
                dvsoftcavity%of_r = dv_dboundary%of_r * env%solvent%dscaled%of_r
            END IF
            !
            CALL mapping%to_small(dvsoftcavity, aux)
            !
            dv(:) = dv(:) + aux%of_r(:)
            !
            CALL dv_dboundary%destroy()
            !
            CALL dvsoftcavity%destroy()
            !
        END IF
        !
        IF (env%setup%lelectrostatic) CALL dvelectrostatic%destroy()
        !
        CALL aux%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dvenviron
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the energy corrections in PW calculations
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_deenviron(env, total_energy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_obj), INTENT(INOUT) :: env
        REAL(DP), INTENT(INOUT) :: total_energy
        !
        !--------------------------------------------------------------------------------
        !
        env%deenviron = -env%system_electrons%density%scalar_product(env%dvtot)
        total_energy = total_energy + env%deenviron
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_deenviron
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_calculator
!----------------------------------------------------------------------------------------
