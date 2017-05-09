!
! Copyright (C) Oliviero Andreussi
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE ions_utils

  USE environ_types
  USE functions

  PRIVATE

  PUBLIC :: create_environ_ions, init_environ_ions_first, init_environ_ions_second, &
       & update_environ_ions, destroy_environ_ions

CONTAINS

  SUBROUTINE create_environ_ions(ions)

    IMPLICIT NONE

    TYPE( environ_ions ), INTENT(INOUT) :: ions
    CHARACTER (LEN=80) :: sub_name = 'create_environ_ions'

    ions%update = .FALSE.

    IF ( ALLOCATED( ions%ityp ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)
    IF ( ALLOCATED( ions%iontype ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)

    NULLIFY( ions%tau )

    ions%use_smeared_ions = .FALSE.
    IF ( ALLOCATED( ions%smeared_ions ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)
    CALL create_environ_density( ions%density )

    ions%use_core_electrons = .FALSE.
    IF ( ALLOCATED( ions%core_electrons ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)
    CALL create_environ_density( ions%core )

    RETURN

  END SUBROUTINE create_environ_ions

  SUBROUTINE init_environ_ions_first( nat, ntyp, lsoftcavity, lcoredensity, &
       &  lsmearedions, radius_mode, atom_label, atomicspread, corespread, solvationrad, ions )

    ! First step of ions initialization, cannot initialize everything
    ! because some infos are missing when this routine is called
    ! NEED TO REVISE THE POSITION OF THE CALL INSIDE QE TO MERGE FIRST AND SECOND

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nat, ntyp
    LOGICAL, INTENT(IN) :: lsoftcavity
    LOGICAL, INTENT(IN) :: lcoredensity
    LOGICAL, INTENT(IN) :: lsmearedions
    CHARACTER( LEN=80 ), INTENT(IN) :: radius_mode
    CHARACTER( LEN=3 ), DIMENSION(ntyp), INTENT(IN) :: atom_label
    REAL( DP ), DIMENSION(ntyp), INTENT(IN) :: atomicspread
    REAL( DP ), DIMENSION(ntyp), INTENT(IN) :: corespread
    REAL( DP ), DIMENSION(ntyp), INTENT(IN) :: solvationrad
    TYPE( environ_ions ), INTENT(INOUT) :: ions

    CHARACTER( LEN = 80 ) :: sub_name = 'init_environ_ions_first'

    INTEGER :: i

    ions%center = 0.D0

    ions%number = nat
    ions%ntyp = ntyp

    ! Allocate the basic vectors, cannot initialize them here

    ALLOCATE( ions%tau( 3, nat ) )
    ALLOCATE( ions%ityp( nat ) )

    ! Set ions types, note that also valence charges cannot be initialized here

    ALLOCATE( ions%iontype( ntyp ) )

    DO i = 1, ntyp

       ! Given the label we could set some of the properties with defaults

       CALL set_iontype_defaults( i, atom_label(i), radius_mode, ions%iontype(i) )

       ! Check if values were provided in input and overwrite them

       IF ( atomicspread(i) .GT. 0 ) ions%iontype(i)%atomicspread = atomicspread(i)
       IF ( corespread(i) .GT. 0 ) ions%iontype(i)%corespread = corespread(i)
       IF ( solvationrad(i) .GT. 0 ) ions%iontype(i)%solvationrad = solvationrad(i)

       ! If need cavity defined exclusively on ions, check radius is not zero

       IF ( .NOT. lsoftcavity .AND. ( ions%iontype(i)%solvationrad .EQ. 0.D0 ) ) &
            & CALL errore(sub_name,'Missing solvation radius for one of the atom types',1)

       ! If need smeared ions, check spread is not zero

       IF ( lsmearedions .AND. ( ions%iontype(i)%atomicspread .EQ. 0.D0 ) ) &
            & CALL errore(sub_name,'Missing atomic spread for one of the atom types',1)

    END DO

    ions%use_smeared_ions = lsmearedions

    ions%use_core_electrons = lcoredensity

    RETURN

  END SUBROUTINE init_environ_ions_first

  SUBROUTINE init_environ_ions_second( nat, ntyp, ityp, zv, cell, ions )

    ! Second step of initialization, passing the information on types,
    ! atomic valence charges and whether we need to compute the smeared density or not

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nat, ntyp
    INTEGER, DIMENSION(nat), INTENT(IN) :: ityp
    REAL(DP), DIMENSION(ntyp), INTENT(IN) :: zv
    TYPE( environ_cell), INTENT(IN) :: cell
    TYPE( environ_ions ), INTENT(INOUT) :: ions

    CHARACTER( LEN = 80 ) :: sub_name = 'init_environ_ions_second'

    INTEGER :: i

    ! Check on dimensions, can skip if merged with first step

    IF ( ions%number .NE. nat ) CALL errore(sub_name,'Mismatch in number of atoms',1)
    IF ( ions%ntyp .NE. ntyp ) CALL errore(sub_name,'Mismatch in number of atom types',1)

    ions%ityp = ityp

    DO i = 1, ions%ntyp

       ions%iontype(i)%zv = zv(i)

    ENDDO

    ions%charge = 0.D0

    DO i = 1, ions%number

       ions%charge = ions%charge + ions%iontype(ions%ityp(i))%zv

    ENDDO

    IF ( ions%use_smeared_ions ) THEN
       ! THE FOLLOWING TEST ON ALLOCATION IS ONLY THERE BECAUSE OF WHEN THIS INITIALIZATION
       ! IS CALLED, IF MERGED WITH THE FIRST STEP REMOVE THE TEST
       IF ( .NOT. ALLOCATED( ions%density%of_r ) ) THEN
          CALL init_environ_density( cell, ions%density )
       ELSE
          ions%density%of_r = 0.D0
       ENDIF

       ! Build smeared ions from iontype data

       IF ( .NOT. ALLOCATED( ions%smeared_ions ) ) THEN
          ALLOCATE( ions%smeared_ions( ions%number ) )
          DO i = 1, ions%number
             ions%smeared_ions(i) = environ_functions(1,1,0,0.0_DP,&
                  & ions%iontype(ions%ityp(i))%atomicspread,&
                  & ions%iontype(ions%ityp(i))%zv,ions%tau(:,i))
          ENDDO
       ENDIF

    ENDIF

    IF ( ions%use_core_electrons ) THEN
       ! THE FOLLOWING TEST ON ALLOCATION IS ONLY THERE BECAUSE OF WHEN THIS INITIALIZATION
       ! IS CALLED, IF MERGED WITH THE FIRST STEP REMOVE THE TEST
       IF ( .NOT. ALLOCATED( ions%core%of_r ) ) THEN
          CALL init_environ_density( cell, ions%core )
       ELSE
          ions%core%of_r = 0.D0
       ENDIF

       ! Build core electrons from iontype data

       IF ( .NOT. ALLOCATED( ions%core_electrons ) ) THEN
          ALLOCATE( ions%core_electrons( ions%number ) )
          DO i = 1, ions%number
             ions%core_electrons(i) = environ_functions(1,1,0,0.0_DP,&
                  & ions%iontype(ions%ityp(i))%corespread,&
                  & ions%iontype(ions%ityp(i))%zv,ions%tau(:,i))
          ENDDO
       ENDIF

    ENDIF

    RETURN

  END SUBROUTINE init_environ_ions_second

  SUBROUTINE update_environ_ions( nat, tau, ions )

    ! Update ionic positions and compute derived quantities

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nat
    REAL(DP), DIMENSION(3,nat), INTENT(IN) :: tau
    TYPE( environ_ions ), INTENT(INOUT) :: ions

    INTEGER :: i
    INTEGER :: dim, axis
    REAL(DP) :: charge, spread
    REAL(DP), DIMENSION(3) :: pos
    CHARACTER(LEN=80) :: sub_name = 'update_environ_ions'

    ! Check on dimensions

    IF ( ions%number .NE. nat ) CALL errore(sub_name,'Mismatch in number of atoms',1)

    ! Update positions

    ions%tau = tau

    ! Center of ionic charge used by three sub-modules

    ions%center = 0.D0
    DO i = 1, ions%number

       ions%center(:) = ions%center(:) + ions%tau(:,i)*ions%iontype(ions%ityp(i))%zv

    ENDDO
    ions%center = ions%center / ions%charge

    ! If needed, generate a fictitious ion density using gaussians

    IF ( ions%use_smeared_ions ) CALL density_of_functions(ions%number,ions%smeared_ions,ions%density)

    RETURN

  END SUBROUTINE update_environ_ions

  SUBROUTINE destroy_environ_ions(lflag, ions)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_ions ), INTENT(INOUT) :: ions
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_ions'

    ! ityp, tau and iontype should have been allocated
    ! raise an error if they are not

    IF ( lflag ) THEN

       ! These components were allocated first, only destroy if lflag = .TRUE.

       IF (.NOT.ALLOCATED(ions%ityp)) &
            & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
       DEALLOCATE( ions%ityp )
       IF (.NOT.ALLOCATED(ions%iontype)) &
            & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
       DEALLOCATE( ions%iontype )

       IF (.NOT.ASSOCIATED(ions%tau)) &
            & CALL errore(sub_name,'Trying to destroy a non associated object',1)
       DEALLOCATE( ions%tau )

    ENDIF

    IF ( ions%use_smeared_ions ) THEN
       CALL destroy_environ_density( ions%density )
       CALL destroy_environ_functions( ions%number, ions%smeared_ions )
    ENDIF

    RETURN

  END SUBROUTINE destroy_environ_ions

  SUBROUTINE set_iontype_defaults(index,label,radius_mode,iontype)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: index
    CHARACTER( LEN=3 ), INTENT(IN) :: label
    CHARACTER( LEN=80 ), INTENT(IN) :: radius_mode
    TYPE(environ_iontype), INTENT(INOUT) :: iontype

    CHARACTER( LEN=80 ) :: sub_name = 'set_iontype_defaults'
    REAL( DP ), DIMENSION(92) :: pauling_radii
    DATA pauling_radii/ 1.20_DP, 0.00_DP, & ! H, -
         & 0.00_DP, 0.00_DP, 0.00_DP, 1.50_DP, 1.50_DP, 1.40_DP, 1.35_DP, 0.00_DP, & ! -, -, -, C, N, O, F, -
         & 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.90_DP, 1.85_DP, 1.80_DP, 0.00_DP, & ! -, -, -, -, P, S, Cl, -
         & 0.00_DP, 0.00_DP, 10*0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.95_DP, 0.00_DP, & ! -,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,Br,-
         & 0.00_DP, 0.00_DP, 10*0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 2.15_DP, 0.00_DP, & ! -,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,I,-
         & 38*0.00_DP / ! ...
    REAL( DP ), DIMENSION(92) :: bondi_radii
    DATA bondi_radii/ 1.20_DP, 1.40_DP, & ! H, He
         & 1.82_DP, 1.45_DP, 1.80_DP, 1.70_DP, 1.55_DP, 1.52_DP, 1.47_DP, 1.54_DP, & ! Li, Be, B, C, N, O, F, Ne
         & 2.27_DP, 1.73_DP, 2.30_DP, 2.10_DP, 1.80_DP, 1.80_DP, 1.75_DP, 1.88_DP, & ! Na, Mg, Al, Si, P, S, Cl, Ar
         & 2.75_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.63_DP, 0.00_DP, & ! K, -, -, -, -, -, -, Ni, -
         & 0.00_DP, 1.40_DP, 1.39_DP, 1.87_DP, 2.19_DP, 1.85_DP, 1.90_DP, 1.85_DP, 2.02_DP, & ! -, Cu, Zn, Ga, Ge, As, Se, Be, Kr
         & 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, & ! -, -, -, -, -, -, -, -, -
         & 1.63_DP, 1.72_DP, 1.58_DP, 1.93_DP, 2.17_DP, 0.00_DP, 2.06_DP, 1.98_DP, 2.16_DP, & ! Pd, Ag, Cd, In, Sn, -, Te, I, Xe
         & 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, & ! -, -, -, -, -, -, -, -
         & 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, & ! -, -, -, -, -, -, -, -
         & 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.75_DP, & ! -, -, -, -, -, -, -, Pt
         & 1.66_DP, 1.55_DP, 1.96_DP, 2.02_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, & ! Au, Hg, Tl, Pb, -, -, -, -
         & 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.86_DP / ! -,-,-,-,-,U
    REAL( DP ), DIMENSION(92) :: UFF_diameters
    DATA UFF_diameters/  2.886_DP, 2.362_DP, & ! H, He
         & 2.451_DP, 2.745_DP, 4.083_DP, 3.851_DP, 3.660_DP, 3.500_DP, 3.364_DP, 3.243_DP, & ! Li, Be, B, C, N, O, F, Ne
         & 2.983_DP, 3.021_DP, 4.499_DP, 4.295_DP, 4.147_DP, 4.035_DP, 3.947_DP, 3.868_DP, & ! Na, Mg, Al, Si, P, S, Cl, Ar
         & 3.812_DP, 3.399_DP, 3.295_DP, 3.175_DP, 3.144_DP, 3.023_DP, 2.961_DP, 2.912_DP, 2.872_DP, & ! K, Ca, Sc, Ti, V, Cr, Mn, Ni, Fe
         & 2.834_DP, 3.495_DP, 2.763_DP, 4.383_DP, 4.280_DP, 4.230_DP, 4.205_DP, 4.189_DP, 4.141_DP, & ! Co, Cu, Zn, Ga, Ge, As, Se, Br, Kr
         & 4.114_DP, 3.641_DP, 3.345_DP, 3.124_DP, 3.165_DP, 3.052_DP, 2.998_DP, 2.963_DP, 2.929_DP, & ! Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh
         & 2.899_DP, 3.148_DP, 2.848_DP, 4.463_DP, 4.392_DP, 4.420_DP, 4.470_DP, 4.500_DP, 4.404_DP, & ! Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe
         & 4.517_DP, 3.703_DP, 3.522_DP, 3.556_DP, 3.606_DP, 3.575_DP, 3.547_DP, 3.520_DP, & ! Cs, Ba, La, Ce, Pr, Nd, Pm, Sm
         & 3.493_DP, 3.368_DP, 3.451_DP, 3.428_DP, 3.409_DP, 3.391_DP, 3.374_DP, 3.355_DP, & ! Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb
         & 3.640_DP, 3.141_DP, 3.170_DP, 3.069_DP, 2.954_DP, 3.120_DP, 2.840_DP, 2.754_DP, & ! Lu, Hf, Ta, W, Re, Os, Ir, Pt
         & 3.293_DP, 2.705_DP, 4.337_DP, 4.297_DP, 4.379_DP, 4.709_DP, 4.750_DP, 4.765_DP, & ! Au, Hg, Tl, Pb, Bi, Po, At, Rn
         & 4.900_DP, 3.677_DP, 3.478_DP, 3.396_DP, 3.424_DP, 3.395_DP / ! Fr, Ra, Ac, Th, Pa, U

    iontype%index = index
    iontype%label = label
    iontype%zv = 0.D0 ! this cannot be initialized here at this time
    iontype%atmnum = get_atmnum(label)

    IF ( iontype%atmnum .EQ. 0 ) &
         & CALL errore(sub_name,'Can not assign the atom type associated with input label',1)

    iontype%atomicspread = 0.5D0
    iontype%corespread = 0.5D0

    SELECT CASE ( radius_mode )

    CASE ( 'pauling' )

       iontype%solvationrad = pauling_radii(iontype%atmnum)

    CASE ( 'bondi' )

       iontype%solvationrad = bondi_radii(iontype%atmnum)

    CASE ( 'uff' )

       iontype%solvationrad = UFF_diameters(iontype%atmnum) * 0.5_DP

    CASE DEFAULT

       CALL errore(sub_name,'Unknown radius_mode',1)

    END SELECT

    RETURN

  END SUBROUTINE set_iontype_defaults

!--------------------------------------------------------------------
  FUNCTION get_atmnum(label)
!wgt  FUNCTION get_atmwgt(label)
!
! original version by O. Andreussi (MIT)
!
!--------------------------------------------------------------------

    INTEGER :: get_atmnum
!wgt    REAL*8 :: get_atmwgt

    CHARACTER*(*), INTENT(IN) :: label
    CHARACTER*3 :: tmplab

    INTEGER :: num
    REAL*8 :: weigth

    tmplab=TRIM(ADJUSTL(label))
    CALL lowcase(tmplab)
    IF (tmplab(1:1).EQ.'a')  THEN
      IF (tmplab(2:2).EQ.'c') THEN
        num=89
        weigth=227.03
      ELSE IF (tmplab(2:2).EQ.'l') THEN
        num=13
        weigth=26.981538
      ELSE IF (tmplab(2:2).EQ.'r') THEN
        num=18
        weigth=39.948
      ELSE IF (tmplab(2:2).EQ.'g') THEN
        num=47
        weigth=107.8682
      ELSE IF (tmplab(2:2).EQ.'s') THEN
        num=33
        weigth=74.9216
      ELSE IF (tmplab(2:2).EQ.'u') THEN
        num=79
        weigth=196.96655
      ELSE IF (tmplab(2:2).EQ.'t') THEN
        num=85
        weigth=210
      ELSE
        num=13
        weigth=26.981538
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'b')  THEN
      IF (tmplab(2:2).EQ.'e') THEN
        num=4
        weigth=9.012182
      ELSE IF (tmplab(2:2).EQ.'r') THEN
        num=35
        weigth=79.904
      ELSE IF (tmplab(2:2).EQ.'a') THEN
        num=56
        weigth=137.327
      ELSE IF (tmplab(2:2).EQ.'i') THEN
        num=83
        weigth=208.98038
      ELSE
        num=5
        weigth=10.811
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'c')  THEN
      IF (tmplab(2:2).EQ.'a') THEN
        num=20
        weigth=40.078
      ELSE IF (tmplab(2:2).EQ.'r') THEN
        num=24
        weigth=51.9961
      ELSE IF (tmplab(2:2).EQ.'o') THEN
        num=27
        weigth=58.9332
      ELSE IF (tmplab(2:2).EQ.'l') THEN
        num=17
        weigth=35.453
      ELSE IF (tmplab(2:2).EQ.'s') THEN
        num=55
        weigth=132.90545
      ELSE IF (tmplab(2:2).EQ.'d') THEN
        num=48
        weigth=112.411
      ELSE IF (tmplab(2:2).EQ.'u') THEN
        num=29
        weigth=63.546
      ELSE IF (tmplab(2:2).EQ.'e') THEN
        num=58
        weigth=140.116
      ELSE
        num=6
        weigth=12.0107
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'f')  THEN
      IF (tmplab(2:2).EQ.'e') THEN
        num=26
        weigth=55.845
      ELSE IF (tmplab(2:2).EQ.'r') THEN
        num=87
        weigth=223
      ELSE
        num=9
        weigth=18.9984032
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'g')  THEN
      IF (tmplab(2:2).EQ.'a') THEN
        num=31
        weigth=69.723
      ELSE IF (tmplab(2:2).EQ.'e') THEN
        num=32
        weigth=72.64
      ELSE
        num=31
        weigth=69.723
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'h')  THEN
      IF (tmplab(2:2).EQ.'e') THEN
        num=2
        weigth=4.002602
      ELSE IF (tmplab(2:2).EQ.'g') THEN
        num=80
        weigth=200.59
      ELSE IF (tmplab(2:2).EQ.'f') THEN
        num=72
        weigth=178.49
      ELSE
        num=1
        weigth=1.00794
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'i')  THEN
      IF (tmplab(2:2).EQ.'n') THEN
        num=49
        weigth=114.818
      ELSE IF (tmplab(2:2).EQ.'r') THEN
        num=77
        weigth=192.217
      ELSE
        num=53
        weigth=126.90447
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'k')  THEN
      IF (tmplab(2:2).EQ.'r') THEN
        num=36
        weigth=83.798
      ELSE
        num=19
        weigth=39.0983
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'l')  THEN
      IF (tmplab(2:2).EQ.'i') THEN
        num=3
        weigth=6.941
      ELSE IF (tmplab(2:2).EQ.'a') THEN
        num=57
        weigth=138.9055
      ELSE
        num=3
        weigth=6.941
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'m')  THEN
      IF (tmplab(2:2).EQ.'g') THEN
        num=12
        weigth=24.3050
      ELSE IF (tmplab(2:2).EQ.'n') THEN
        num=25
        weigth=54.938049
      ELSE IF (tmplab(2:2).EQ.'o') THEN
        num=42
        weigth=95.94
      ELSE
        num=25
        weigth=54.938049
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'n')  THEN
      IF (tmplab(2:2).EQ.'a') THEN
        num=11
        weigth=22.98977
      ELSE IF (tmplab(2:2).EQ.'i') THEN
        num=28
        weigth=58.6934
      ELSE IF (tmplab(2:2).EQ.'b') THEN
        num=41
        weigth=92.90638
      ELSE IF (tmplab(2:2).EQ.'e') THEN
        num=10
        weigth=20.1797
      ELSE IF (tmplab(2:2).EQ.'d') THEN
        num=60
        weigth=144.24
      ELSE
        num=7
        weigth=14.0067
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'o')  THEN
      IF (tmplab(2:2).EQ.'s') THEN
        num=76
        weigth=190.23
      ELSE
        num=8
        weigth=15.9994
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'p')  THEN
      IF (tmplab(2:2).EQ.'a') THEN
        num=91
        weigth=231.04
      ELSE IF (tmplab(2:2).EQ.'d') THEN
        num=46
        weigth=106.42
      ELSE IF (tmplab(2:2).EQ.'t') THEN
        num=78
        weigth=195.078
      ELSE IF (tmplab(2:2).EQ.'b') THEN
        num=82
        weigth=207.2
      ELSE IF (tmplab(2:2).EQ.'o') THEN
        num=84
        weigth=209
      ELSE
        num=15
        weigth=30.973761
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'r')  THEN
      IF (tmplab(2:2).EQ.'b') THEN
        num=37
        weigth=85.4678
      ELSE IF (tmplab(2:2).EQ.'u') THEN
        num=44
        weigth=101.07
      ELSE IF (tmplab(2:2).EQ.'h') THEN
        num=45
        weigth=102.90550
      ELSE IF (tmplab(2:2).EQ.'e') THEN
        num=75
        weigth=186.207
      ELSE IF (tmplab(2:2).EQ.'n') THEN
        num=86
        weigth=222
      ELSE
        num=37
        weigth=85.4678
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'s')  THEN
      IF (tmplab(2:2).EQ.'i') THEN
        num=14
        weigth=28.0855
      ELSE IF (tmplab(2:2).EQ.'e') THEN
        num=34
        weigth=78.96
      ELSE IF (tmplab(2:2).EQ.'c') THEN
        num=21
        weigth=44.955910
      ELSE IF (tmplab(2:2).EQ.'r') THEN
        num=38
        weigth=87.62
      ELSE IF (tmplab(2:2).EQ.'n') THEN
        num=50
        weigth=118.710
      ELSE IF (tmplab(2:2).EQ.'b') THEN
        num=51
        weigth=121.760
      ELSE
        num=16
        weigth=32.065
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'t')  THEN
      IF (tmplab(2:2).EQ.'a') THEN
        num=73
        weigth=180.9479
      ELSE IF (tmplab(2:2).EQ.'l') THEN
        num=81
        weigth=204.3833
      ELSE IF (tmplab(2:2).EQ.'c') THEN
        num=43
        weigth=98
      ELSE IF (tmplab(2:2).EQ.'h') THEN
        num=90
        weigth=232.04
      ELSE IF (tmplab(2:2).EQ.'i') THEN
        num=22
        weigth=47.867
      ELSE IF (tmplab(2:2).EQ.'e') THEN
        num=52
        weigth=127.60
      ELSE
        num=22
        weigth=47.867
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'u')  THEN
      num=92
      weigth=238.02891
    ELSE IF (tmplab(1:1).EQ.'v')  THEN
      num=23
      weigth=50.9415
    ELSE IF (tmplab(1:1).EQ.'w')  THEN
      num=74
      weigth=183.84
    ELSE IF (tmplab(1:1).EQ.'x')  THEN
      num=54
      weigth=131.293
    ELSE IF (tmplab(1:1).EQ.'y')  THEN
      num=39
      weigth=88.90585
    ELSE IF (tmplab(1:1).EQ.'z')  THEN
      IF (tmplab(2:2).EQ.'n') THEN
        num=30
        weigth=65.409
      ELSE IF (tmplab(2:2).EQ.'r') THEN
        num=40
        weigth=91.224
      ELSE
        num=30
        weigth=65.409
      ENDIF
    ELSE
      num=0
      weigth=0
    ENDIF

    get_atmnum=num
!    get_atmwgt=weigth

!--------------------------------------------------------------------
  END FUNCTION get_atmnum
!wgt  END FUNCTION get_atmwgt
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE lowcase(string)
!--------------------------------------------------------------------

    CHARACTER*(*), INTENT(inout) :: string

    INTEGER :: i, length
    CHARACTER(1) :: letter

    length=LEN(string)

    DO i = 1,min(255,length)
      letter = string(i:i)
      IF(letter.eq.'A') THEN
        letter = 'a'
      ELSE IF(letter.eq.'B') THEN
        letter = 'b'
      ELSE IF(letter.eq.'C') THEN
        letter = 'c'
      ELSE IF(letter.eq.'D') THEN
        letter = 'd'
      ELSE IF(letter.eq.'E') THEN
        letter = 'e'
      ELSE IF(letter.eq.'F') THEN
        letter = 'f'
      ELSE IF(letter.eq.'G') THEN
        letter = 'g'
      ELSE IF(letter.eq.'H') THEN
        letter = 'h'
      ELSE IF(letter.eq.'I') THEN
        letter = 'i'
      ELSE IF(letter.eq.'J') THEN
        letter = 'j'
      ELSE IF(letter.eq.'K') THEN
        letter = 'k'
      ELSE IF(letter.eq.'L') THEN
        letter = 'l'
      ELSE IF(letter.eq.'M') THEN
        letter = 'm'
      ELSE IF(letter.eq.'N') THEN
        letter = 'n'
      ELSE IF(letter.eq.'O') THEN
        letter = 'o'
      ELSE IF(letter.eq.'P') THEN
        letter = 'p'
      ELSE IF(letter.eq.'Q') THEN
        letter = 'q'
      ELSE IF(letter.eq.'R') THEN
        letter = 'r'
      ELSE IF(letter.eq.'S') THEN
        letter = 's'
      ELSE IF(letter.eq.'T') THEN
        letter = 't'
      ELSE IF(letter.eq.'U') THEN
        letter = 'u'
      ELSE IF(letter.eq.'V') THEN
        letter = 'v'
      ELSE IF(letter.eq.'W') THEN
        letter = 'w'
      ELSE IF(letter.eq.'X') THEN
        letter = 'x'
      ELSE IF(letter.eq.'Y') THEN
        letter = 'y'
      ELSE IF(letter.eq.'Z') THEN
        letter = 'z'
      END IF
      string(i:i) = letter
    END DO

    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE lowcase
!--------------------------------------------------------------------

END MODULE ions_utils
