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
!
!----------------------------------------------------------------------------------------
!>
!! Module containing the main routines to handle environ_ions
!! derived data types.
!!
!! Environ_ions contains all the details of the ions of the QM system,
!! including the atomic type, mass and charge, the solvation radii,
!! the sizes of the ionic cores and the spread to be used when treating
!! the ions as gaussians. The module also contains all the routines
!! to handle environ_ions and to generate smooth ionic density from the
!! ions specification.
!!
!----------------------------------------------------------------------------------------
MODULE utils_ions
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: e2, pi, tpi, fpi, bohr_radius_angs
    USE environ_types
    USE environ_base, ONLY: potential_shift
    USE environ_output
    USE utils_functions
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=2) :: elements(92)
    !
    REAL(DP), DIMENSION(92) :: pauling_radii, bondi_radii, UFF_diameters, &
                               MUFF_diameters, weights
    !
    DATA elements/ &
        'h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne', 'na', 'mg', 'al', 'si', &
        'p', 's', 'cl', 'ar', 'k', 'ca', 'sc', 'ti', 'v', 'cr', 'mn', 'fe', 'co', &
        'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr', 'y', 'zr', &
        'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 'te', 'i', &
        'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', &
        'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w', 're', 'os', 'ir', 'pt', 'au', &
        'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn', 'fr', 'ra', 'ac', 'th', 'pa', 'u'/
    !
    DATA pauling_radii/1.20_DP, 0.00_DP, & ! H, -
        0.00_DP, 0.00_DP, 0.00_DP, 1.50_DP, 1.50_DP, 1.40_DP, 1.35_DP, 0.00_DP, & ! -, -, -, C, N, O, F, -
        0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.90_DP, 1.85_DP, 1.80_DP, 0.00_DP, & ! -, -, -, -, P, S, Cl, -
        0.00_DP, 0.00_DP, 10*0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.95_DP, 0.00_DP, & ! -,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,Br,-
        0.00_DP, 0.00_DP, 10*0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 2.15_DP, 0.00_DP, & ! -,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,I,-
        38*0.00_DP/ ! ...
    !
    DATA bondi_radii/1.20_DP, 1.40_DP, & ! H, He
        1.82_DP, 1.45_DP, 1.80_DP, 1.70_DP, 1.55_DP, 1.52_DP, 1.47_DP, 1.54_DP, & ! Li, Be, B, C, N, O, F, Ne
        2.27_DP, 1.73_DP, 2.30_DP, 2.10_DP, 1.80_DP, 1.80_DP, 1.75_DP, 1.88_DP, & ! Na, Mg, Al, Si, P, S, Cl, Ar
        2.75_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.63_DP, 0.00_DP, & ! K, -, -, -, -, -, -, Ni, -
        0.00_DP, 1.40_DP, 1.39_DP, 1.87_DP, 2.19_DP, 1.85_DP, 1.90_DP, 1.85_DP, 2.02_DP, & ! -, Cu, Zn, Ga, Ge, As, Se, Be, Kr
        0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, & ! -, -, -, -, -, -, -, -, -
        1.63_DP, 1.72_DP, 1.58_DP, 1.93_DP, 2.17_DP, 0.00_DP, 2.06_DP, 1.98_DP, 2.16_DP, & ! Pd, Ag, Cd, In, Sn, -, Te, I, Xe
        0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, & ! -, -, -, -, -, -, -, -
        0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, & ! -, -, -, -, -, -, -, -
        0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.75_DP, & ! -, -, -, -, -, -, -, Pt
        1.66_DP, 1.55_DP, 1.96_DP, 2.02_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, & ! Au, Hg, Tl, Pb, -, -, -, -
        0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.86_DP/ ! -,-,-,-,-,U
    !
    DATA UFF_diameters/2.886_DP, 2.362_DP, & ! H, He
        2.451_DP, 2.745_DP, 4.083_DP, 3.851_DP, 3.660_DP, 3.500_DP, 3.364_DP, 3.243_DP, & ! Li, Be, B, C, N, O, F, Ne
        2.983_DP, 3.021_DP, 4.499_DP, 4.295_DP, 4.147_DP, 4.035_DP, 3.947_DP, 3.868_DP, & ! Na, Mg, Al, Si, P, S, Cl, Ar
        3.812_DP, 3.399_DP, 3.295_DP, 3.175_DP, 3.144_DP, 3.023_DP, 2.961_DP, 2.912_DP, 2.872_DP, & ! K, Ca, Sc, Ti, V, Cr, Mn, Ni, Fe
        2.834_DP, 3.495_DP, 2.763_DP, 4.383_DP, 4.280_DP, 4.230_DP, 4.205_DP, 4.189_DP, 4.141_DP, & ! Co, Cu, Zn, Ga, Ge, As, Se, Br, Kr
        4.114_DP, 3.641_DP, 3.345_DP, 3.124_DP, 3.165_DP, 3.052_DP, 2.998_DP, 2.963_DP, 2.929_DP, & ! Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh
        2.899_DP, 3.148_DP, 2.848_DP, 4.463_DP, 4.392_DP, 4.420_DP, 4.470_DP, 4.500_DP, 4.404_DP, & ! Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe
        4.517_DP, 3.703_DP, 3.522_DP, 3.556_DP, 3.606_DP, 3.575_DP, 3.547_DP, 3.520_DP, & ! Cs, Ba, La, Ce, Pr, Nd, Pm, Sm
        3.493_DP, 3.368_DP, 3.451_DP, 3.428_DP, 3.409_DP, 3.391_DP, 3.374_DP, 3.355_DP, & ! Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb
        3.640_DP, 3.141_DP, 3.170_DP, 3.069_DP, 2.954_DP, 3.120_DP, 2.840_DP, 2.754_DP, & ! Lu, Hf, Ta, W, Re, Os, Ir, Pt
        3.293_DP, 2.705_DP, 4.337_DP, 4.297_DP, 4.379_DP, 4.709_DP, 4.750_DP, 4.765_DP, & ! Au, Hg, Tl, Pb, Bi, Po, At, Rn
        4.900_DP, 3.677_DP, 3.478_DP, 3.396_DP, 3.424_DP, 3.395_DP/ ! Fr, Ra, Ac, Th, Pa, U
    !
    DATA MUFF_diameters/2.886_DP, 2.362_DP, & ! H, He
        2.451_DP, 2.745_DP, 4.083_DP, 3.851_DP, 3.100_DP, 3.500_DP, 3.364_DP, 3.243_DP, & ! Li, Be, B, C, N, O, F, Ne
        2.983_DP, 3.021_DP, 4.499_DP, 4.295_DP, 4.147_DP, 4.035_DP, 3.947_DP, 3.868_DP, & ! Na, Mg, Al, Si, P, S, Cl, Ar
        3.812_DP, 3.399_DP, 3.295_DP, 3.175_DP, 3.144_DP, 3.023_DP, 2.961_DP, 2.912_DP, 2.872_DP, & ! K, Ca, Sc, Ti, V, Cr, Mn, Ni, Fe
        2.834_DP, 3.495_DP, 2.763_DP, 4.383_DP, 4.280_DP, 4.230_DP, 4.205_DP, 4.189_DP, 4.141_DP, & ! Co, Cu, Zn, Ga, Ge, As, Se, Br, Kr
        4.114_DP, 3.641_DP, 3.345_DP, 3.124_DP, 3.165_DP, 3.052_DP, 2.998_DP, 2.963_DP, 2.929_DP, & ! Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh
        2.899_DP, 3.148_DP, 2.848_DP, 4.463_DP, 4.392_DP, 4.420_DP, 4.470_DP, 4.500_DP, 4.404_DP, & ! Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe
        4.517_DP, 3.703_DP, 3.522_DP, 3.556_DP, 3.606_DP, 3.575_DP, 3.547_DP, 3.520_DP, & ! Cs, Ba, La, Ce, Pr, Nd, Pm, Sm
        3.493_DP, 3.368_DP, 3.451_DP, 3.428_DP, 3.409_DP, 3.391_DP, 3.374_DP, 3.355_DP, & ! Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb
        3.640_DP, 3.141_DP, 3.170_DP, 3.069_DP, 2.954_DP, 3.120_DP, 2.840_DP, 2.754_DP, & ! Lu, Hf, Ta, W, Re, Os, Ir, Pt
        3.293_DP, 2.705_DP, 4.337_DP, 4.297_DP, 4.379_DP, 4.709_DP, 4.750_DP, 4.765_DP, & ! Au, Hg, Tl, Pb, Bi, Po, At, Rn
        4.900_DP, 3.677_DP, 3.478_DP, 3.396_DP, 3.424_DP, 3.395_DP/ ! Fr, Ra, Ac, Th, Pa, U
    !
    DATA weights/1.00794_DP, 4.002602_DP, & ! H, He
        6.941_DP, 9.012182_DP, 10.811_DP, 12.0107_DP, 14.0067_DP, 15.9994_DP, 18.9984032_DP, 20.1797_DP, & ! Li, Be, B, C, N, O, F, Ne
        22.98977_DP, 24.3050_DP, 26.981538_DP, 28.0855_DP, 30.973761_DP, 32.065_DP, 35.453_DP, 39.948_DP, & ! Na, Mg, Al, Si, P, S, Cl, Ar
        39.0983_DP, 40.078_DP, 44.955910_DP, 47.867_DP, 50.9415_DP, 51.9961_DP, 54.938049_DP, 58.6934_DP, 55.845_DP, & ! K, Ca, Sc, Ti, V, Cr, Mn, Ni, Fe
        58.9332_DP, 63.546_DP, 65.409_DP, 69.723_DP, 72.64_DP, 74.9216_DP, 78.96_DP, 79.904_DP, 83.798_DP, & ! Co, Cu, Zn, Ga, Ge, As, Se, Br, Kr
        85.4678_DP, 87.62_DP, 88.90585_DP, 91.224_DP, 92.90638_DP, 95.94_DP, 98.0_DP, 101.07_DP, 102.90550_DP, & ! Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh
        106.42_DP, 107.8682_DP, 112.411_DP, 114.818_DP, 118.710_DP, 121.760_DP, 127.60_DP, 126.90447_DP, 131.293_DP, & ! Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe
        132.90545_DP, 137.327_DP, 138.9055_DP, 140.116_DP, 140.90766_DP, 144.24_DP, 145.0_DP, 150.36_DP, & ! Cs, Ba, La, Ce, Pr, Nd, Pm, Sm
        151.964_DP, 157.25_DP, 158.925354_DP, 162.500_DP, 164.930328_DP, 167.259_DP, 168.934218_DP, 173.045_DP, & ! Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb
        174.9668_DP, 178.49_DP, 180.9479_DP, 183.84_DP, 186.207_DP, 190.23_DP, 192.217_DP, 195.078_DP, & ! Lu, Hf, Ta, W, Re, Os, Ir, Pt
        196.96655_DP, 200.59_DP, 204.3833_DP, 207.2_DP, 208.98038_DP, 209.0_DP, 210.0_DP, 222.0_DP, & ! Au, Hg, Tl, Pb, Bi, Po, At, Rn
        223.0_DP, 226.0_DP, 227.03_DP, 232.04_DP, 231.04_DP, 238.02891_DP/ ! Fr, Ra, Ac, Th, Pa, U
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: create_environ_ions, init_environ_ions_first, init_environ_ions_second, &
              update_environ_ions, destroy_environ_ions
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_ions(ions)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_ions), INTENT(INOUT) :: ions
        !
        CHARACTER(LEN=80) :: label = ' '
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_ions'
        !
        !--------------------------------------------------------------------------------
        !
        ions%update = .FALSE.
        !
        IF (ALLOCATED(ions%ityp)) &
            CALL errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        IF (ALLOCATED(ions%iontype)) &
            CALL errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        NULLIFY (ions%tau)
        !
        ions%use_smeared_ions = .FALSE.
        !
        IF (ALLOCATED(ions%smeared_ions)) &
            CALL errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        label = 'smeared_ions'
        !
        CALL create_environ_density(ions%density, label)
        !
        ions%use_core_electrons = .FALSE.
        !
        IF (ALLOCATED(ions%core_electrons)) &
            CALL errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        label = 'core_electrons'
        !
        CALL create_environ_density(ions%core, label)
        !
        IF (ALLOCATED(ions%vloc)) &
            CALL errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_ions
    !------------------------------------------------------------------------------------
    !>
    !! First step of ions initialization, cannot initialize everything
    !! because some information is missing when this routine is called
    !! NEED TO REVISE THE POSITION OF THE CALL INSIDE QE TO MERGE FIRST AND SECOND
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_ions_first(nat, ntyp, lsoftcavity, lcoredensity, &
                                       lsmearedions, radius_mode, atom_label, &
                                       atomicspread, corespread, solvationrad, ions)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat, ntyp
        LOGICAL, INTENT(IN) :: lsoftcavity, lcoredensity, lsmearedions
        CHARACTER(LEN=80), INTENT(IN) :: radius_mode
        CHARACTER(LEN=3), INTENT(IN) :: atom_label(ntyp)
        REAL(DP), DIMENSION(ntyp), INTENT(IN) :: atomicspread, corespread, solvationrad
        !
        TYPE(environ_ions), INTENT(INOUT) :: ions
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_ions_first'
        !
        !--------------------------------------------------------------------------------
        !
        ions%center = 0.D0
        ions%dipole = 0.D0
        ions%quadrupole_pc = 0.D0
        ions%quadrupole_gauss = 0.D0
        ions%quadrupole_correction = 0.D0
        ions%selfenergy_correction = 0.D0
        !
        ions%number = nat
        ions%ntyp = ntyp
        !
        !--------------------------------------------------------------------------------
        ! Allocate the basic vectors, cannot initialize them here
        !
        ions%initialized = .FALSE.
        ALLOCATE (ions%tau(3, nat))
        ions%tau = 0.D0
        ALLOCATE (ions%ityp(nat))
        ions%ityp = 0
        ions%alat = 0.D0
        !
        !--------------------------------------------------------------------------------
        ! Set ions types, note that also valence charges cannot be initialized here
        !
        ALLOCATE (ions%iontype(ntyp))
        !
        ALLOCATE (ions%vloc(ntyp))
        !
        DO i = 1, ntyp
            !
            CALL set_iontype_defaults(i, atom_label(i), radius_mode, ions%iontype(i))
            ! given the label we could set some of the properties with defaults
            !
            !----------------------------------------------------------------------------
            ! Check if values were provided in input and overwrite them
            !
            IF (atomicspread(i) > 0) ions%iontype(i)%atomicspread = atomicspread(i)
            !
            IF (corespread(i) > 0) ions%iontype(i)%corespread = corespread(i)
            !
            IF (solvationrad(i) > 0) ions%iontype(i)%solvationrad = solvationrad(i)
            !
            !----------------------------------------------------------------------------
            ! If need cavity defined exclusively on ions, check radius is not zero
            !
            IF (.NOT. lsoftcavity .AND. (ions%iontype(i)%solvationrad == 0.D0)) &
                CALL errore(sub_name, &
                            'Missing solvation radius for one of the atom types', 1)
            !
            !----------------------------------------------------------------------------
            ! If need smeared ions, check spread is not zero
            !
            IF (lsmearedions .AND. (ions%iontype(i)%atomicspread == 0.D0)) &
                CALL errore(sub_name, &
                            'Missing atomic spread for one of the atom types', 1)
            !
            CALL create_environ_density(ions%vloc(i))
            !
        END DO
        !
        ions%use_smeared_ions = lsmearedions
        ions%use_core_electrons = lcoredensity
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_ions_first
    !------------------------------------------------------------------------------------
    !>
    !! Second step of initialization, passing the information on types, atomic valence
    !! charges and whether we need to compute the smeared density or not
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_ions_second(nat, ntyp, nnr, ityp, zv, cell, ions, vloc)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat, ntyp, nnr
        INTEGER, INTENT(IN) :: ityp(nat)
        REAL(DP), INTENT(IN) :: zv(ntyp)
        TYPE(environ_cell), INTENT(IN) :: cell
        REAL(DP), OPTIONAL, INTENT(IN) :: vloc(nnr, ntyp)
        !
        TYPE(environ_ions), INTENT(INOUT) :: ions
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_ions_second'
        !
        !--------------------------------------------------------------------------------
        ! Check on dimensions, can skip if merged with first step
        !
        IF (ions%number /= nat) &
            CALL errore(sub_name, 'Mismatch in number of atoms', 1)
        !
        IF (ions%ntyp /= ntyp) &
            CALL errore(sub_name, 'Mismatch in number of atom types', 1)
        !
        ions%alat = cell%alat ! needed because the ionic positions are scaled by alat
        ions%ityp = ityp
        !
        DO i = 1, ions%ntyp
            ions%iontype(i)%zv = -zv(i)
            !
            !----------------------------------------------------------------------------
            ! Store pseudopotential in R-space in ions%vloc%of_r
            !
            IF (PRESENT(vloc)) THEN
                !
                !------------------------------------------------------------------------
                ! Input/reference potentials in system cell
                !
                IF (.NOT. ASSOCIATED(ions%vloc(i)%cell)) &
                    CALL init_environ_density(cell, ions%vloc(i))
                !
                ions%vloc(i)%of_r = vloc(:, i)
            END IF
            !
        END DO
        !
        ions%charge = 0.D0
        !
        DO i = 1, ions%number
            ions%charge = ions%charge + ions%iontype(ions%ityp(i))%zv
        END DO
        !
        IF (ions%use_smeared_ions) THEN
            !
            !----------------------------------------------------------------------------
            ! #TODO The following test on allocation is only there because of when this
            ! initialization is called. If merged with the first step, remove the test.
            !
            IF (.NOT. ALLOCATED(ions%density%of_r)) THEN
                CALL init_environ_density(cell, ions%density)
            ELSE
                ions%density%of_r = 0.D0
            END IF
            !
            !----------------------------------------------------------------------------
            ! Build smeared ions from iontype data
            !
            IF (.NOT. ALLOCATED(ions%smeared_ions)) THEN
                ALLOCATE (ions%smeared_ions(ions%number))
                !
                DO i = 1, ions%number
                    ions%smeared_ions(i) = &
                        environ_functions(1, 1, 0, 0.0_DP, &
                                          ions%iontype(ions%ityp(i))%atomicspread, &
                                          ions%iontype(ions%ityp(i))%zv, &
                                          ions%tau(:, i))
                END DO
                !
            END IF
            !
        END IF
        !
        IF (ions%use_core_electrons) THEN
            !
            !----------------------------------------------------------------------------
            ! #TODO The following test on allocation is only there because of when this
            ! initialization is called. If merged with the first step, remove the test.
            !
            IF (.NOT. ALLOCATED(ions%core%of_r)) THEN
                CALL init_environ_density(cell, ions%core)
            ELSE
                ions%core%of_r = 0.D0
            END IF
            !
            !----------------------------------------------------------------------------
            ! Build core electrons from iontype data
            !
            IF (.NOT. ALLOCATED(ions%core_electrons)) THEN
                ALLOCATE (ions%core_electrons(ions%number))
                !
                DO i = 1, ions%number
                    !
                    ions%core_electrons(i) = &
                        environ_functions(1, 1, 0, 0.0_DP, &
                                          ions%iontype(ions%ityp(i))%corespread, &
                                          -ions%iontype(ions%ityp(i))%zv, &
                                          ions%tau(:, i))
                    !
                END DO
                !
            END IF
            !
        END IF
        !
        ions%initialized = .TRUE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_ions_second
    !------------------------------------------------------------------------------------
    !>
    !! Update ionic positions and compute derived quantities
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_ions(nat, tau, ions)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        REAL(DP), INTENT(IN) :: tau(3, nat)
        !
        TYPE(environ_ions), INTENT(INOUT) :: ions
        !
        INTEGER :: i
        INTEGER :: dim, axis
        REAL(DP) :: charge, spread
        REAL(DP) :: pos(3)
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_ions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ions%number /= nat) &
            CALL errore(sub_name, 'Mismatch in number of atoms', 1)
        !
        ions%tau = tau ! update positions
        !
        !--------------------------------------------------------------------------------
        ! Center of ionic charge used by three sub-modules
        !
        ions%center = 0.D0
        !
        DO i = 1, ions%number
            !
            ions%center(:) = ions%center(:) + &
                             ions%tau(:, i) * ions%iontype(ions%ityp(i))%zv
            !
        END DO
        !
        IF (ABS(ions%charge) < 1.D-8) &
            CALL errore(sub_name, 'Ionic charge equal to zero', 1)
        !
        ions%center = ions%center / ions%charge
        !
        !--------------------------------------------------------------------------------
        ! If needed, generate a fictitious ion density using gaussians
        !
        IF (ions%use_smeared_ions) &
            CALL density_of_functions(ions%number, ions%smeared_ions, ions%density, &
                                      .TRUE.)
        !
        !--------------------------------------------------------------------------------
        ! Compute quadrupole moment of point-like (and gaussian) nuclei
        !
        ions%dipole = 0.D0 ! this is due to the choice of ionic center
        ions%quadrupole_pc = 0.D0
        ions%quadrupole_correction = 0.D0
        ions%selfenergy_correction = 0.D0
        !
        DO i = 1, ions%number
            !
            ions%quadrupole_pc(:) = ions%quadrupole_pc(:) + &
                                    ions%iontype(ions%ityp(i))%zv * &
                                    ((ions%tau(:, i) - ions%center(:)) * ions%alat)**2
            !
            IF (ions%use_smeared_ions) THEN
                !
                ions%quadrupole_correction = &
                    ions%quadrupole_correction + ions%iontype(ions%ityp(i))%zv * &
                    ions%iontype(ions%ityp(i))%atomicspread**2 * 0.5D0
                !
                ions%selfenergy_correction = &
                    ions%selfenergy_correction + ions%iontype(ions%ityp(i))%zv**2 / &
                    ions%iontype(ions%ityp(i))%atomicspread * SQRT(2.D0 / pi)
                !
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Calculate potential shift due to Gaussian nuclei
        !
        IF (ions%use_smeared_ions) &
            potential_shift = ions%quadrupole_correction * &
                              tpi * e2 / ions%density%cell%omega
        !
        IF (ions%use_smeared_ions) &
            ions%quadrupole_gauss(:) = ions%quadrupole_pc(:) + &
                                       ions%quadrupole_correction
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_ions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_ions(lflag, ions)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(environ_ions), INTENT(INOUT) :: ions
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_ions'
        !
        !--------------------------------------------------------------------------------
        ! ityp, tau and iontype should have been allocated
        ! raise an error if they are not
        !
        IF (ions%initialized) THEN
            !
            IF (ions%use_smeared_ions) THEN
                !
                CALL destroy_environ_density(ions%density)
                !
                CALL destroy_environ_functions(ions%number, ions%smeared_ions)
                !
            END IF
            !
            ions%charge = 0.D0
            !
            DO i = 1, ions%ntyp
                CALL destroy_environ_density(ions%vloc(i))
            END DO
            !
            ions%initialized = .FALSE.
            !
        END IF
        !
        IF (lflag) THEN
            !
            IF (.NOT. ALLOCATED(ions%ityp)) &
                CALL errore(sub_name, 'Trying to destroy a non allocated object', 1)
            !
            DEALLOCATE (ions%ityp)
            !
            IF (.NOT. ALLOCATED(ions%iontype)) &
                CALL errore(sub_name, 'Trying to destroy a non allocated object', 1)
            !
            DEALLOCATE (ions%iontype)
            !
            IF (.NOT. ASSOCIATED(ions%tau)) &
                CALL errore(sub_name, 'Trying to destroy a non associated object', 1)
            !
            DEALLOCATE (ions%tau)
            DEALLOCATE (ions%vloc)
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_ions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_iontype_defaults(index, label, radius_mode, iontype)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: index
        CHARACTER(LEN=3), INTENT(IN) :: label
        CHARACTER(LEN=80), INTENT(IN) :: radius_mode
        !
        TYPE(environ_iontype), INTENT(INOUT) :: iontype
        !
        CHARACTER(LEN=80) :: sub_name = 'set_iontype_defaults'
        !
        !--------------------------------------------------------------------------------
        !
        iontype%index = index
        iontype%label = label
        iontype%atmnum = get_atmnum(label)
        !
        IF (iontype%atmnum == 0) &
            CALL errore(sub_name, &
                        'Can not assign the atom type associated with input label', 1)
        !
        ! iontype%weight = weights(iontype%atmnum) ! #TODO future work
        !
        iontype%atomicspread = 0.5D0
        iontype%corespread = 0.5D0
        !
        SELECT CASE (radius_mode)
        CASE ('pauling')
            iontype%solvationrad = pauling_radii(iontype%atmnum)
        CASE ('bondi')
            iontype%solvationrad = bondi_radii(iontype%atmnum)
        CASE ('uff')
            iontype%solvationrad = UFF_diameters(iontype%atmnum) * 0.5_DP
        CASE ('muff')
            iontype%solvationrad = MUFF_diameters(iontype%atmnum) * 0.5_DP
        CASE DEFAULT
            CALL errore(sub_name, 'Unknown radius_mode', 1)
        END SELECT
        !
        iontype%solvationrad = iontype%solvationrad / bohr_radius_angs
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_iontype_defaults
    !------------------------------------------------------------------------------------
    !>
    !! original version by O. Andreussi (MIT)
    !! modified by Edan Bainglass (UNT)
    !!
    !------------------------------------------------------------------------------------
    FUNCTION get_atmnum(label)
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=3), INTENT(IN) :: label
        !
        INTEGER :: i, get_atmnum
        CHARACTER(LEN=2) :: clean_label, lower
        !
        !--------------------------------------------------------------------------------
        !
        clean_label = TRIM(ADJUSTL(label))
        !
        CALL lowcase(clean_label, lower)
        !
        get_atmnum = 0
        !
        DO i = 1, SIZE(elements)
            IF (lower == elements(i)) THEN
                get_atmnum = i
                EXIT
            END IF
        END DO
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_atmnum
    !------------------------------------------------------------------------------------
    !>
    !! Uses ASCII table values to shift uppercase (65-90) to lowercase (97-122)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE lowcase(string, lower_str)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=2), INTENT(IN) :: string
        !
        CHARACTER(LEN=2), INTENT(OUT) :: lower_str
        !
        CHARACTER(LEN=1) :: c
        INTEGER :: i, ci
        !
        !--------------------------------------------------------------------------------
        !
        lower_str = string
        !
        DO i = 1, LEN(string)
            c = string(i:i)
            ci = ICHAR(c)
            !
            IF (ci >= 65 .AND. ci <= 90) THEN
                lower_str(i:i) = CHAR(ci + 32)
            END IF
            !
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE lowcase
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_ions
!----------------------------------------------------------------------------------------
