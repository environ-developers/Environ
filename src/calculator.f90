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
    USE class_environ
    USE class_setup
    !
    USE class_cell
    USE class_density
    USE class_gradient
    !
    USE class_boundary_electronic
    USE class_boundary_ionic
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
        TYPE(environ_main), POINTER :: main => NULL()
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: force => calc_fenviron
        PROCEDURE :: energy => calc_eenviron
        PROCEDURE :: denergy => calc_deenviron
        PROCEDURE :: potential => calc_venviron
        PROCEDURE :: dpotential => calc_dvenviron
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_calculator
    !------------------------------------------------------------------------------------
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
    SUBROUTINE calc_venviron(this, update, local_verbose)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: update
        INTEGER, OPTIONAL, INTENT(IN) :: local_verbose
        !
        CLASS(environ_calculator), TARGET, INTENT(INOUT) :: this
        !
        TYPE(environ_density) :: aux
        TYPE(environ_density) :: de_dboundary
        !
        TYPE(environ_main), POINTER :: main
        TYPE(environ_setup), POINTER :: setup
        !
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_venviron'
        !
        !--------------------------------------------------------------------------------
        !
        main => this%main
        setup => main%setup
        !
        !--------------------------------------------------------------------------------
        ! If not updating, add old potentials and exit
        !
        IF (.NOT. update) THEN
            !
            IF (PRESENT(local_verbose)) THEN
                !
                CALL write_cube(main%dvtot, main%system_ions, local_verbose)
                !
                IF (setup%lelectrostatic) THEN
                    !
                    CALL write_cube(main%vreference, main%system_ions, local_verbose)
                    !
                    CALL write_cube(main%velectrostatic, main%system_ions, local_verbose)
                    !
                    CALL main%system_charges%printout(local_verbose)
                    !
                END IF
                !
                IF (setup%lconfine) &
                    CALL write_cube(main%vconfine, main%system_ions, local_verbose)
                !
                IF (setup%lsoftcavity) &
                    CALL write_cube(main%vsoftcavity, main%system_ions, local_verbose)
                !
            END IF
            !
            RETURN
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! If updating, calculate new potentials
        !
        main%dvtot%of_r = 0.D0
        !
        CALL aux%init(setup%system_cell)
        !
        !--------------------------------------------------------------------------------
        ! If any form of electrostatic embedding is present, calculate its contribution
        !
        IF (setup%lelectrostatic) THEN
            !
            !----------------------------------------------------------------------------
            ! Electrostatics is also computed inside the calling program,
            ! need to remove the reference #TODO to-be-decided
            !
            CALL setup%reference%calc_v(main%system_charges, main%vreference)
            !
            CALL write_cube(main%vreference, main%system_ions)
            !
            CALL setup%outer%calc_v(main%environment_charges, main%velectrostatic)
            !
            ! IF (this%setup%lexternals) CALL this%environment_charges%update()
            ! #TODO keep for now until external tests are fully debugged
            !
            CALL write_cube(main%velectrostatic, main%system_ions)
            !
            CALL setup%mapping%to_small(main%velectrostatic, aux)
            !
            main%dvtot%of_r = aux%of_r - main%vreference%of_r
            !
            CALL main%environment_charges%of_potential(main%velectrostatic)
            !
            CALL main%environment_charges%printout()
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! #TODO add brief description of confine
        !
        IF (setup%lconfine) THEN
            !
            CALL main%solvent%vconfine(setup%confine, main%vconfine)
            !
            CALL write_cube(main%vconfine, main%system_ions)
            !
            CALL setup%mapping%to_small(main%vconfine, aux)
            !
            main%dvtot%of_r = main%dvtot%of_r + aux%of_r
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Compute the total potential depending on the boundary
        !
        IF (setup%lsoftcavity) THEN
            main%vsoftcavity%of_r = 0.D0
            !
            CALL de_dboundary%init(setup%environment_cell)
            !
            IF (setup%lsoftsolvent) THEN
                de_dboundary%of_r = 0.D0
                !
                ! if surface tension greater than zero, calculate cavity contribution
                IF (setup%lsurface) &
                    CALL main%solvent%desurface_dboundary(setup%surface_tension, &
                                                          de_dboundary)
                !
                ! if external pressure different from zero, calculate PV contribution
                IF (setup%lvolume) &
                    CALL main%solvent%devolume_dboundary(setup%pressure, &
                                                         de_dboundary)
                !
                ! if confinement potential not zero, calculate confine contribution
                IF (setup%lconfine) &
                    CALL main%solvent%deconfine_dboundary(setup%confine, &
                                                          main%environment_charges%electrons%density, &
                                                          de_dboundary)
                !
                ! if dielectric embedding, calculate dielectric contribution
                IF (setup%lstatic) &
                    CALL main%static%de_dboundary(main%velectrostatic, de_dboundary)
                !
                ! if solvent-aware interface correct the potential
                IF (main%solvent%solvent_aware) &
                    CALL main%solvent%sa_de_dboundary(de_dboundary)
                !
                SELECT TYPE (solvent => main%solvent)
                    !
                TYPE IS (environ_boundary_electronic)
                    main%vsoftcavity%of_r = de_dboundary%of_r * solvent%dscaled%of_r
                    !
                END SELECT
                !
            END IF
            !
            IF (setup%lsoftelectrolyte) THEN
                de_dboundary%of_r = 0.D0
                !
                CALL main%electrolyte%de_dboundary(de_dboundary)
                ! if electrolyte is present add its non-electrostatic contribution
                !
                ! if solvent-aware interface correct the potential
                IF (main%electrolyte%boundary%solvent_aware) &
                    CALL main%electrolyte%boundary%sa_de_dboundary(de_dboundary)
                !
                SELECT TYPE (boundary => main%electrolyte%boundary)
                    !
                TYPE IS (environ_boundary_electronic)
                    !
                    main%vsoftcavity%of_r = main%vsoftcavity%of_r + &
                                            de_dboundary%of_r * boundary%dscaled%of_r
                    !
                END SELECT
                !
            END IF
            !
            CALL write_cube(main%vsoftcavity, main%system_ions)
            !
            CALL setup%mapping%to_small(main%vsoftcavity, aux)
            !
            main%dvtot%of_r = main%dvtot%of_r + aux%of_r
            !
            CALL de_dboundary%destroy()
            !
        END IF
        !
        CALL aux%destroy()
        !
        CALL write_cube(main%dvtot, main%system_ions, local_verbose)
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
    SUBROUTINE calc_eenviron(this, total_energy)
        !--------------------------------------------------------------------------------
        !
        ! USE embedding_confine, ONLY: calc_econfine #TODO to-be-decided
        !
        IMPLICIT NONE
        !
        CLASS(environ_calculator), TARGET, INTENT(INOUT) :: this
        REAL(DP), INTENT(INOUT) :: total_energy
        !
        REAL(DP) :: ereference
        !
        TYPE(environ_main), POINTER :: main
        TYPE(environ_setup), POINTER :: setup
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_eenviron'
        !
        !--------------------------------------------------------------------------------
        !
        main => this%main
        setup => main%setup
        !
        !--------------------------------------------------------------------------------
        !
        main%eelectrostatic = 0.D0
        main%esurface = 0.D0
        main%evolume = 0.D0
        main%econfine = 0.D0
        main%eelectrolyte = 0.D0
        !
        setup%niter = setup%niter + 1
        !
        !--------------------------------------------------------------------------------
        ! If electrostatic is on, compute electrostatic energy
        !
        IF (setup%lelectrostatic) THEN
            !
            CALL setup%reference%calc_e(main%system_charges, main%vreference, ereference)
            !
            CALL setup%outer%calc_e(main%environment_charges, main%velectrostatic, &
                                    main%eelectrostatic)
            !
            main%eelectrostatic = main%eelectrostatic - ereference
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        ! if surface tension not zero, compute cavitation energy
        IF (setup%lsurface) &
            CALL main%solvent%esurface(setup%surface_tension, main%esurface)
        !
        IF (setup%lvolume) CALL main%solvent%evolume(setup%pressure, main%evolume)
        ! if pressure not zero, compute PV energy
        !
        ! if confinement potential not zero compute confine energy
        IF (setup%lconfine) &
            main%econfine = main%environment_electrons%density%scalar_product(main%vconfine)
        !
        ! if electrolyte is present, calculate its non-electrostatic contribution
        IF (setup%lelectrolyte) CALL main%electrolyte%energy(main%eelectrolyte)
        !
        total_energy = total_energy + main%eelectrostatic + main%esurface + &
                       main%evolume + main%econfine + main%eelectrolyte
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
    SUBROUTINE calc_fenviron(this, nat, force_environ)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        !
        CLASS(environ_calculator), TARGET, INTENT(INOUT) :: this
        REAL(DP), INTENT(INOUT) :: force_environ(3, nat)
        !
        INTEGER :: i
        TYPE(environ_density) :: de_dboundary
        TYPE(environ_gradient) :: partial
        !
        REAL(DP), DIMENSION(3, nat) :: freference, felectrostatic
        !
        TYPE(environ_main), POINTER :: main
        TYPE(environ_setup), POINTER :: setup
        TYPE(environ_cell), POINTER :: environment_cell
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_fenviron'
        !
        !--------------------------------------------------------------------------------
        !
        main => this%main
        setup => main%setup
        environment_cell => setup%environment_cell
        !
        !--------------------------------------------------------------------------------
        ! Compute the electrostatic embedding contribution to the interatomic forces.
        ! If using doublecell, take the difference of contributions from the full charge
        ! density computed in the environment cell (outer solver) and those from the
        ! charge density of ions/electrons computed in the system cell (reference solver)
        !
        force_environ = 0.D0
        !
        IF (setup%lelectrostatic) THEN
            !
            IF (setup%ldoublecell) THEN
                !
                CALL setup%reference%calc_f(nat, main%system_charges, freference, &
                                            setup%ldoublecell)
                !
            ELSE
                freference = 0.D0
            END IF
            !
            CALL setup%outer%calc_f(nat, main%environment_charges, felectrostatic, &
                                    setup%ldoublecell)
            !
            force_environ = felectrostatic - freference
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Compute the total forces depending on the boundary
        !
        IF (setup%lrigidcavity) THEN
            !
            CALL de_dboundary%init(environment_cell)
            !
            CALL partial%init(environment_cell)
            !
            IF (setup%lrigidsolvent) THEN
                de_dboundary%of_r = 0.D0
                !
                ! if surface tension greater than zero, calculate cavity contribution
                IF (setup%lsurface) &
                    CALL main%solvent%desurface_dboundary(setup%surface_tension, &
                                                          de_dboundary)
                !
                ! if external pressure not zero, calculate PV contribution
                IF (setup%lvolume) &
                    CALL main%solvent%devolume_dboundary(setup%pressure, de_dboundary)
                !
                ! if confinement potential not zero, calculate confine contribution
                IF (setup%lconfine) &
                    CALL main%solvent%deconfine_dboundary(setup%confine, &
                                                          main%environment_charges%electrons%density, &
                                                          de_dboundary)
                !
                ! if dielectric embedding, calculate dielectric contribution
                IF (setup%lstatic) &
                    CALL main%static%de_dboundary(main%velectrostatic, de_dboundary)
                !
                ! if solvent-aware, correct the potential
                IF (main%solvent%solvent_aware) &
                    CALL main%solvent%sa_de_dboundary(de_dboundary)
                !
                IF (main%solvent%field_aware) CALL main%solvent%ion_field_partial()
                !
                !------------------------------------------------------------------------
                ! Multiply by derivative of the boundary w.r.t ionic positions
                !
                DO i = 1, nat
                    !
                    CALL main%solvent%dboundary_dions(i, partial)
                    !
                    IF (main%solvent%field_aware) &
                        CALL main%solvent%fa_dboundary_dions(i, partial)
                    !
                    force_environ(:, i) = force_environ(:, i) - &
                                          partial%scalar_product_density(de_dboundary)
                    !
                END DO
                !
            END IF
            !
            IF (setup%lrigidelectrolyte) THEN
                de_dboundary%of_r = 0.D0
                !
                ! if electrolyte is present, add its non-electrostatic contribution
                CALL main%electrolyte%de_dboundary(de_dboundary)
                !
                ! if solvent-aware, correct the potential
                IF (main%electrolyte%boundary%solvent_aware) &
                    CALL main%electrolyte%boundary%sa_de_dboundary(de_dboundary)
                !
                IF (main%electrolyte%boundary%field_aware) &
                    CALL main%electrolyte%boundary%ion_field_partial()
                !
                !------------------------------------------------------------------------
                ! Multiply by derivative of the boundary w.r.t ionic positions
                !
                DO i = 1, nat
                    !
                    CALL main%electrolyte%boundary%dboundary_dions(i, partial)
                    !
                    IF (main%electrolyte%boundary%field_aware) &
                        CALL main%electrolyte%boundary%fa_dboundary_dions(i, partial)
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
    SUBROUTINE calc_dvenviron(this, nnr, dv)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        !
        CLASS(environ_calculator), TARGET, INTENT(INOUT) :: this
        REAL(DP), INTENT(INOUT) :: dv(nnr)
        !
        TYPE(environ_density) :: aux
        TYPE(environ_density) :: dvreference
        TYPE(environ_density) :: dvelectrostatic
        TYPE(environ_density) :: dvsoftcavity
        TYPE(environ_density) :: dv_dboundary
        !
        TYPE(environ_main), POINTER :: main
        TYPE(environ_setup), POINTER :: setup
        TYPE(environ_cell), POINTER :: system_cell, environment_cell
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_dvenviron'
        !
        !--------------------------------------------------------------------------------
        !
        main => this%main
        setup => main%setup
        system_cell => setup%system_cell
        environment_cell => setup%environment_cell
        !
        !--------------------------------------------------------------------------------
        !
        IF (setup%optical_permittivity == 1.D0) RETURN
        !
        !--------------------------------------------------------------------------------
        !
        CALL aux%init(system_cell)
        !
        !--------------------------------------------------------------------------------
        ! If any form of electrostatic embedding is present, calculate its contribution
        !
        IF (setup%lelectrostatic) THEN
            !
            !----------------------------------------------------------------------------
            ! Electrostatics is also computed inside the calling program,
            ! need to remove the reference
            !
            CALL dvreference%init(system_cell)
            !
            CALL setup%reference%calc_v(main%system_response_charges, dvreference)
            !
            CALL dvelectrostatic%init(environment_cell)
            !
            CALL setup%outer%calc_v(main%environment_response_charges, dvelectrostatic)
            !
            CALL setup%mapping%to_small(dvelectrostatic, aux)
            !
            dv = dv + aux%of_r - dvreference%of_r
            !
            CALL dvreference%destroy()
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Compute the response potential depending on the boundary
        !
        IF (setup%lsoftcavity) THEN
            !
            CALL dvsoftcavity%init(environment_cell)
            !
            CALL dv_dboundary%init(environment_cell)
            !
            IF (setup%lsoftsolvent) THEN
                dv_dboundary%of_r = 0.D0
                !
                ! if dielectric embedding, calcultes dielectric contribution
                IF (setup%loptical) &
                    CALL main%optical%dv_dboundary(main%velectrostatic, &
                                                   dvelectrostatic, dv_dboundary)
                !
                SELECT TYPE (solvent => main%solvent)
                    !
                TYPE IS (environ_boundary_electronic)
                    dvsoftcavity%of_r = dv_dboundary%of_r * solvent%dscaled%of_r
                    !
                END SELECT
                !
            END IF
            !
            CALL setup%mapping%to_small(dvsoftcavity, aux)
            !
            dv = dv + aux%of_r
            !
            CALL dv_dboundary%destroy()
            !
            CALL dvsoftcavity%destroy()
            !
        END IF
        !
        IF (setup%lelectrostatic) CALL dvelectrostatic%destroy()
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
    SUBROUTINE calc_deenviron(this, total_energy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_calculator), TARGET, INTENT(INOUT) :: this
        REAL(DP), INTENT(INOUT) :: total_energy
        !
        TYPE(environ_main), POINTER :: main
        !
        !--------------------------------------------------------------------------------
        !
        main => this%main
        !
        !--------------------------------------------------------------------------------
        !
        main%deenviron = -main%system_electrons%density%scalar_product(main%dvtot)
        total_energy = total_energy + main%deenviron
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_deenviron
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_calculator
!----------------------------------------------------------------------------------------
