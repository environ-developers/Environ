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
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_iontype
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_char_ops, ONLY: env_lowercase

    USE environ_param, ONLY: DP, BOHR_RADIUS_ANGS
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: print_environ_iontypes, get_element
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE, PUBLIC :: environ_iontype
        !--------------------------------------------------------------------------------
        !
        INTEGER :: index
        INTEGER :: atmnum
        CHARACTER(LEN=3) :: label
        REAL(DP) :: zv
        REAL(DP) :: atomicspread
        REAL(DP) :: corespread
        REAL(DP) :: solvationrad
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: init => init_environ_iontype
        PROCEDURE :: set_defaults => set_iontype_defaults
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_iontype
    !------------------------------------------------------------------------------------
    !
    INTERFACE get_element
        MODULE PROCEDURE get_element_by_number, get_element_by_weight
    END INTERFACE get_element
    !
    !------------------------------------------------------------------------------------
    !
    CHARACTER(LEN=3) :: elements(92)
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
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   ADMIN METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_iontype(this, index, atom_label, zv, radius_mode, &
                                    atomicspread, corespread, solvationrad, &
                                    lsoftcavity, lsmearedions)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: index
        CHARACTER(LEN=*), INTENT(IN) :: atom_label
        CHARACTER(LEN=*), INTENT(IN) :: radius_mode
        REAL(DP), INTENT(IN) :: zv, atomicspread, corespread, solvationrad
        LOGICAL, INTENT(IN) :: lsoftcavity, lsmearedions
        !
        CLASS(environ_iontype), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_iontype'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%set_defaults(index, atom_label, radius_mode)
        !
        this%zv = -zv
        !
        IF (atomicspread > 0) this%atomicspread = atomicspread
        !
        IF (this%label == 'H') THEN
            this%corespread = 1.D-10
        ELSE IF (corespread > 0) THEN
            this%corespread = corespread
        END IF
        !
        IF (solvationrad > 0) this%solvationrad = solvationrad
        !
        !--------------------------------------------------------------------------------
        ! If cavity is defined exclusively on ions, check that radius is not zero
        !
        IF (.NOT. lsoftcavity .AND. (this%solvationrad == 0.D0)) &
            CALL io%error(sub_name, &
                          'Missing solvation radius for one of the atom types', 1)
        !
        !--------------------------------------------------------------------------------
        ! If using smeared ions, check that spread is not zero
        !
        IF (lsmearedions .AND. (this%atomicspread == 0.D0)) &
            CALL io%error(sub_name, &
                          'Missing atomic spread for one of the atom types', 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_iontype
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               PRIVATE HELPER METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_iontype_defaults(this, index, label, radius_mode)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: index
        CHARACTER(LEN=*), INTENT(IN) :: label
        CHARACTER(LEN=*), INTENT(IN) :: radius_mode
        !
        CLASS(environ_iontype), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'set_iontype_defaults'
        !
        !--------------------------------------------------------------------------------
        !
        this%index = index
        this%label = label
        this%atmnum = get_atmnum(label)
        !
        IF (this%atmnum == 0) &
            CALL io%error(sub_name, &
                          'Cannot assign the atom type associated with input label', 1)
        !
        ! this%weight = weights(this%atmnum) ! #TODO future work
        !
        this%atomicspread = 0.5D0
        this%corespread = 0.5D0
        !
        SELECT CASE (radius_mode)
            !
        CASE ('pauling')
            this%solvationrad = pauling_radii(this%atmnum)
            !
        CASE ('bondi')
            this%solvationrad = bondi_radii(this%atmnum)
            !
        CASE ('uff')
            this%solvationrad = UFF_diameters(this%atmnum) * 0.5_DP
            !
        CASE ('muff')
            this%solvationrad = MUFF_diameters(this%atmnum) * 0.5_DP
            !
        CASE DEFAULT
            CALL io%error(sub_name, "Unknown radius_mode", 1)
            !
        END SELECT
        !
        this%solvationrad = this%solvationrad / BOHR_RADIUS_ANGS
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
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: label
        !
        INTEGER :: i, get_atmnum
        CHARACTER(LEN=2) :: lowcase_label
        !
        !--------------------------------------------------------------------------------
        !
        lowcase_label = env_lowercase(TRIM(ADJUSTL(label)))
        !
        get_atmnum = 0
        !
        DO i = 1, SIZE(elements)
            !
            IF (lowcase_label == elements(i)) THEN
                get_atmnum = i
                !
                EXIT
                !
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_atmnum
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               PUBLIC HELPER METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    CHARACTER(LEN=3) FUNCTION get_element_by_number(number) RESULT(label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: number
        !
        CHARACTER(LEN=80) :: sub_name = 'get_element_by_number'
        !
        !--------------------------------------------------------------------------------
        !
        IF (number > SIZE(elements)) THEN
            WRITE (io%unit, '(/, 5X, "Atomic number = ", I10)') number
            !
            CALL io%error(sub_name, "Atomic number out of bounds", 1)
            !
        END IF
        !
        label = elements(number)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_element_by_number
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    CHARACTER(LEN=3) FUNCTION get_element_by_weight(weight) RESULT(label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: weight
        !
        INTEGER :: i, index
        !
        CHARACTER(LEN=80) :: sub_name = 'get_element_by_weight'
        !
        !--------------------------------------------------------------------------------
        !
        index = -1
        !
        DO i = 1, SIZE(weights)
            !
            IF (ABS(weights(i) - weight) < 1.D-1) THEN
                index = i
                !
                EXIT
                !
            END IF
            !
        END DO
        !
        IF (index < 0) THEN
            WRITE (io%unit, '(/, 5X, "Atomic weight = ", F9.5)') weight
            !
            CALL io%error(sub_name, "Wrong atomic weight", 1)
            !
        END IF
        !
        label = elements(index)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_element_by_weight
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Prints the details of the iontypes
    !!
    !! If called by a parent object, prints details in block format
    !!
    !! @param verbose       : (INTEGER) adds verbosity to global verbose
    !! @param debug_verbose : (INTEGER) replaces global verbose for debugging
    !! @param unit          : (INTEGER) output target (default = io%debug_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_iontypes(this, ntyp, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ntyp
        TYPE(environ_iontype), INTENT(IN) :: this(ntyp)
        INTEGER, OPTIONAL, INTENT(IN) :: verbose, debug_verbose, unit
        !
        INTEGER :: i
        INTEGER :: base_verbose, local_verbose, local_unit
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_iontypes'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        IF (PRESENT(debug_verbose)) THEN
            base_verbose = debug_verbose
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = verbose
            ELSE
                local_verbose = debug_verbose
            END IF
            !
        ELSE IF (io%verbosity > 0) THEN
            base_verbose = io%verbosity
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = base_verbose + verbose
            ELSE
                local_verbose = base_verbose
            END IF
            !
        ELSE
            RETURN
        END IF
        !
        IF (PRESENT(unit)) THEN
            local_unit = unit
        ELSE
            local_unit = io%debug_unit
        END IF
        !
        IF (local_verbose >= 1) THEN
            !
            IF (local_verbose >= base_verbose) THEN ! header
                WRITE (local_unit, 1000)
            ELSE
                !
                CALL io%block_divider(local_verbose, base_verbose, local_unit)
                !
                WRITE (local_unit, 1001)
            END IF
            !
            IF (local_verbose < 3) THEN
                WRITE (local_unit, 1002) ! table headers
                !
                DO i = 1, ntyp
                    !
                    WRITE (local_unit, 1003) &
                        this(i)%index, this(i)%label, &
                        this(i)%atmnum, this(i)%zv
                    !
                END DO
                !
            ELSE
                WRITE (local_unit, 1004) ! table headers
                !
                DO i = 1, ntyp
                    !
                    WRITE (local_unit, 1005) &
                        this(i)%index, this(i)%label, &
                        this(i)%atmnum, this(i)%zv, &
                        this(i)%atomicspread, this(i)%corespread, &
                        this(i)%solvationrad
                    !
                END DO
                !
            END IF
            !
            IF (local_verbose < base_verbose) &
                CALL io%block_divider(local_verbose, base_verbose, local_unit)
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), " IONTYPES ", 66('%'))
1001    FORMAT(/, " IONTYPES", /, " ========")
        !
1002    FORMAT(/, "               atomic", /, &
                "   i | label | number | charge", /, 1X, 29('-'))
        !
1003    FORMAT(1X, I3, " | ", A5, " | ", I6, " | ", F6.2)
!
1004    FORMAT(/, "               atomic            atomic    core    solvation", /, &
                "   i | label | number | charge | spread | spread |    radius", /, &
                1X, 59('-'))
!
1005    FORMAT(1X, I3, " | ", A5, " | ", I6, 3(" | ", F6.2), " |", F10.2)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_iontypes
    !------------------------------------------------------------------------------------
    !
    !
    !------------------------------------------------------------------------------------
END MODULE class_iontype
!----------------------------------------------------------------------------------------
