
! Copyright (C) 2017 Environ group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by Q. Campbell
!
!----------------------------------------------------------------------------
MODULE correction_ms_gcs
!----------------------------------------------------------------------------
  !
  ! ...
  !
  USE environ_types
  USE electrostatic_types
  USE environ_output,    ONLY : environ_unit
  USE environ_base,      ONLY : e2, semiconductor
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: calc_vms_gcs, calc_gradvms_gcs
  !
CONTAINS
!---------------------------------------------------------------------------
  SUBROUTINE calc_vms_gcs( oned_analytic, electrolyte, semiconductor_in, charges, potential )
!---------------------------------------------------------------------------
    !
    ! ... Given the total explicit charge, the value of the field at the boundary
    !     is obtained by Gauss's law
    ! (1) ez = - tpi * e2 * charge * axis_length / omega / env_static_permittivity
    !
    ! ... By integrating the Mott Schottky model, we can find that the schottky barrier
    !     is equal to
    ! (2) vms = fact*(ez)**2 + kbt
    !
    !     where
    ! (3) fact = sc_permittivity/tpi / e2 /2.D0 /carrier_density
    !
    !  You can then find the potential as a parabolic function of distance
    ! (4) v_analytic = (distance)**2/fact/4 - ez*(distance)
    !
    !  This parabolic function is only accurate until the maximum is reached
    !  We can quantify when this will happen as the depletion length
    !
    ! (5) depletion_length = 2.D0 *fact*ez
    !
    !  After the depletion length, the potential will remain flat
    !  For more details and derivation see the appendix in Ch11 of Schmickler
    !  Interfacial Electrochemistry book
    !
    IMPLICIT NONE
!    SAVE
    !
    TYPE( oned_analytic_core ), TARGET, INTENT(IN) :: oned_analytic
    TYPE( environ_semiconductor ), TARGET, INTENT(INOUT) :: semiconductor_in
    TYPE( environ_electrolyte ), TARGET, INTENT(IN) :: electrolyte
    TYPE( environ_density ), TARGET, INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential
    !
    INTEGER, POINTER :: nnr
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER, POINTER :: env_periodicity
    INTEGER, POINTER :: slab_axis
    REAL( DP ), POINTER :: alat, omega, axis_length
    REAL( DP ), DIMENSION(:), POINTER :: origin
    REAL( DP ), DIMENSION(:,:), POINTER :: axis
    !
    REAL( DP ), POINTER :: cion, zion, permittivity_ms, permittivity_gcs
    REAL( DP ), POINTER :: xstern_ms, xstern_gcs, electrode_charge, carrier_density
    REAL( DP ) :: kbt, invkbt
    !
    TYPE( environ_density ), TARGET :: local
    REAL( DP ), DIMENSION(:), POINTER :: v
    !
    INTEGER :: i, icount
    REAL( DP ) :: ez, ez_ms,ez_gcs, fact, vms, vstern
    REAL( DP ) :: arg, const, depletion_length
    REAL( DP ) :: dv, vbound, v_cut, v_edge
    REAL( DP ) :: asinh, coth, acoth
    REAL( DP ) :: f1, f2, max_axis
    REAL( DP ) :: area, vtmp, distance
    REAL(DP) :: dipole(0:3), quadrupole(3)
    REAL(DP) :: tot_charge, tot_dipole(3), tot_quadrupole(3)
    CHARACTER( LEN = 80 ) :: sub_name = 'calc_vms_gcs'
    !
    CALL start_clock ('calc_vms_gcs')
    !
    ! ... Aliases and sanity checks
    !
    IF ( .NOT. ASSOCIATED( potential%cell, charges%cell ) ) &
         & CALL errore(sub_name,'Missmatch in domains of potential and charges',1)
    IF ( potential % cell % nnr .NE. oned_analytic % n ) &
         & CALL errore(sub_name,'Missmatch in domains of potential and solver',1)
    cell => potential % cell
    nnr => cell % nnr
    !


    alat => oned_analytic % alat
    omega => oned_analytic % omega
    env_periodicity => oned_analytic % d
    slab_axis => oned_analytic % axis
    axis_length => oned_analytic % size
    origin => oned_analytic % origin
    axis => oned_analytic % x

    !
    ! ... Get parameters of semiconductor to compute analytic correction
    !

    IF ( electrolyte % ntyp .NE. 2 ) &
         & CALL errore(sub_name,'Unexpected number of counterionic species, different from two',1)
    cion => electrolyte % ioncctype(1) % cbulk
    zion => electrolyte % ioncctype(1) % z
    permittivity_gcs => electrolyte%permittivity
    xstern_gcs => electrolyte%boundary%simple%width


    permittivity_ms => semiconductor_in%permittivity
    carrier_density => semiconductor_in%carrier_density
    electrode_charge => semiconductor_in%electrode_charge
    xstern_ms => semiconductor_in%simple%width

    WRITE( environ_unit, * )"carrier density: ",carrier_density

    !
    WRITE( environ_unit, * )"MS xstern: ",xstern_ms
    WRITE( environ_unit, * )"carrier density: ",carrier_density
    ! ... Set Boltzmann factors
    !
    kbt = semiconductor_in % temperature * k_boltzmann_ry
    invkbt = 1.D0 / kbt
    !
    IF ( env_periodicity .NE. 2 ) &
         & CALL errore(sub_name,'Option not yet implemented: 1D Poisson-Boltzmann solver only for 2D systems',1)
    !
    CALL init_environ_density( cell, local )
    v => local % of_r
   !
    ! ... Compute multipoles of the system wrt the chosen origin
    !
    CALL compute_dipole( nnr, charges%of_r, origin, dipole, quadrupole )
    !
    tot_charge = dipole(0)

    tot_dipole = dipole(1:3)
    tot_quadrupole = quadrupole
    area = omega / axis_length
    semiconductor_in%surf_area_per_sq_cm = area* 2.8002D-17     ! Value of 1 square bohr in cm^2
    semiconductor%surf_area_per_sq_cm = area* 2.8002D-17
    !
    ! ... First apply parabolic correction
    !
    fact = e2 * tpi / omega
    const = - pi / 3.D0 * tot_charge / axis_length * e2 - fact * tot_quadrupole(slab_axis)
    v(:) = - tot_charge * axis(1,:)**2 + 2.D0 * tot_dipole(slab_axis) * axis(1,:)
    v(:) = fact * v(:) + const
    dv = - fact * 4.D0 * tot_dipole(slab_axis) * xstern_gcs
    !
    ! ... Compute the physical properties of the interface
    !
    ! Starting with the GCS props

    zion = ABS(zion)
    IF (ABS(tot_charge) < 1.D-6) THEN
      ez_gcs = 0.D0
    ELSE
      ez_gcs =  tpi * e2 * electrode_charge / area ! / permittivity
    END IF

    ez = - tpi * e2 * tot_charge / area ! / permittivity
    WRITE (environ_unit, *)"ez: ",ez
    fact = - e2 * SQRT( 8.D0 * fpi * cion * kbt / e2 ) !/ permittivity )
    arg = ez_gcs/fact
    asinh = LOG(arg + SQRT( arg**2 + 1 ))
    vstern = 2.D0 * kbt / zion * asinh
    arg = vstern * 0.25D0 * invkbt * zion
    coth = ( EXP( 2.D0 * arg ) + 1.D0 ) / ( EXP( 2.D0 * arg ) - 1.D0 )
    const = coth * EXP( zion * fact * invkbt * 0.5D0 * xstern_gcs )


    !
    ! Gathering the potential at the boundary condition
    vbound = 0.D0
    icount = 0
    DO i = 1, nnr
       !
       IF ( ABS(axis(1,i)) .GE. xstern_gcs ) THEN
          !
          icount = icount + 1
          vbound = vbound + potential % of_r(i) + v(i) - ez * (ABS(axis(1,i)) - xstern_gcs)
          !
       ENDIF
       !
    ENDDO
    CALL mp_sum(icount,cell%comm)
    CALL mp_sum(vbound,cell%comm)
    vbound = vbound / DBLE(icount)
    WRITE (environ_unit, *)"vbound: ",vbound
    !
    ! ... Compute some constants needed for the gcs calculation
    !
    f1 = - fact * zion * invkbt * 0.5D0
    f2 = 4.D0 * kbt / zion
    !
    ! ... Compute the analytic potential and charge
    ! ... adding in gcs effect first, only on the positive side of xstern
    !
    WRITE (environ_unit, *)"vstern: ",vstern
    v = v - vbound + vstern
    DO i = 1, nnr
       !
       IF ( axis(1,i) .GE. xstern_gcs ) THEN
          !
          ! ... Gouy-Chapmann-Stern analytic solution on the outside
          !
          arg = const * EXP( ABS(axis(1,i)) * f1 )
          IF ( ABS(arg) .GT. 1.D0 ) THEN
             acoth = 0.5D0 * LOG( (arg + 1.D0) / (arg - 1.D0) )
          ELSE
             acoth = 0.D0
          END IF
          vtmp =  f2 * acoth


          ! Having to add extra handling for electrode charge

          IF ( ISNAN(vtmp) ) THEN
            vtmp = 0.D0
          END IF
          !
          ! ... Remove source potential (linear) and add analytic one
          !
          WRITE( environ_unit, *)"v_gcs corr: ",vtmp
          v(i) =  v(i) + vtmp - vstern - ez * (ABS(axis(1,i))-xstern_gcs) !+ ez_gcs * xstern_gcs ! vtmp - potential % of_r(i)
          !
          WRITE( environ_unit, *)"v_i: ",v(i)
       ENDIF
       !
    ENDDO

    !
    !   Now adding in ms contribution on negative side of  axis
    !

    ! Now moving on to the ms props
    WRITE( environ_unit, *)"charge: ",tot_charge
    IF (ABS(tot_charge) < 1.D-6) THEN
      ez_ms = 0.D0
    ELSE
      ez_ms= tpi * e2 * (-electrode_charge-tot_charge) / area ! / permittivity !in units of Ry/bohr
    END IF
    WRITE( environ_unit, * )"Mott Schottky electric field: ",ez_ms
    fact = 1.D0/tpi / e2 /2.D0 /carrier_density !*permittivity
    WRITE(  environ_unit, *)"MS Prefactor: ",fact
    arg = fact* (ez_ms**2.D0)
    vms =  arg ! +kbt
    !Finds the total length of the depletion region
    depletion_length = ABS(2.D0 *fact*ez_ms)
    WRITE ( environ_unit, * )"depletion length: ",depletion_length
    WRITE ( environ_unit, * )"vms: ",vms


    !
    ! Finding the v cutoff potential
    v_cut= 0.D0
    icount = 0
    DO i = 1, nnr

       !
       !WRITE (environ_unit, *)"axis : ",ABS(axis(1,i) - xstern_ms)

       IF (( axis(1,i) .LT. 0.D0 ) .AND. ( (ABS(axis(1,i)) - xstern_ms) .LE. 0.1D0 )) THEN
          !
          icount = icount + 1
          v_cut = v_cut + potential % of_r(i) + v(i)
          !
       ENDIF
       !
    ENDDO

    CALL mp_sum(icount,cell%comm)
    CALL mp_sum(v_cut,cell%comm)
    WRITE (environ_unit, *)"v_cut: ",v_cut
    WRITE (environ_unit, *)"icount: ",icount
    v_cut = v_cut / DBLE(icount)
    WRITE (environ_unit, *)"v_cut: ",v_cut

    semiconductor_in%bulk_sc_fermi = v_cut+ vms+ semiconductor_in%flatband_fermi
    semiconductor%bulk_sc_fermi = v_cut+ vms+ semiconductor%flatband_fermi
    WRITE ( environ_unit, * )"bulk semiconductor fermi level: ",semiconductor_in%bulk_sc_fermi

    DO i = 1, nnr

       IF ( -axis(1,i) .GE. xstern_ms ) THEN
          distance = ABS(axis(1,i)) - xstern_ms
          !Only applies parabolic equation if still within the depletion width
          !
          IF ( distance <= depletion_length) THEN
             !
             ! ... Mott Schottky analytic solution on the outside
             !
             IF (ez_ms < 0) THEN
                vtmp = -(distance)**2.D0 / fact/4.D0 + ez_ms*(distance)
             ELSE
                vtmp = (distance)**2.D0 / fact/4.D0 - ez_ms*(distance)
             END IF
          ELSE
             vtmp = 0.D0
          END IF
          ! WRITE (environ_unit, *)"This is the axis value: ",axis(1,i)
          ! WRITE (environ_unit, *) "Distance: ", distance
          ! WRITE (environ_unit, *) "ms correction: ", vtmp
          !
          ! ... Remove source potential (linear) and add analytic one
          !
          v(i) =  v(i) + vtmp  -ez*distance!-vms ! vtmp - potential % of_r(i)
          !WRITE( environ_unit, *)"This is the vi: ",ez*distance
          !
       ENDIF
       !
    ENDDO


    ! Adjust the potential to always be 0 on left side
    ! First determine max axis point
    max_axis = 0.0
    DO i = 1, nnr
       !
       IF ( axis(1,i) > max_axis) THEN
          !
          max_axis = axis(1,i)
          !
       ENDIF
       !
    ENDDO



    v_edge= 0.D0
    icount = 0
    DO i = 1, nnr
       !
       IF ( axis(1,i) .EQ. max_axis ) THEN
          !
          icount = icount + 1
          v_edge = v_edge + potential % of_r(i) + v(i)
          !
       ENDIF
       !
    ENDDO
    v_edge = v_edge / DBLE(icount)
    v = v - v_edge
    !
    potential % of_r = potential % of_r + v

    semiconductor = semiconductor_in
    !
    CALL destroy_environ_density(local)
    !
    CALL stop_clock ('calc_vms_gcs')
    !




    RETURN



    !
!---------------------------------------------------------------------------
END SUBROUTINE calc_vms_gcs
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  SUBROUTINE calc_gradvms_gcs( oned_analytic, electrolyte, semiconductor_in, charges, gradv )
!---------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( oned_analytic_core ), TARGET, INTENT(IN) :: oned_analytic
    TYPE( environ_electrolyte ), TARGET, INTENT(IN) :: electrolyte
    TYPE( environ_semiconductor ), TARGET, INTENT(IN) :: semiconductor_in
    TYPE( environ_density ), TARGET, INTENT(IN) :: charges
    TYPE( environ_gradient ), INTENT(INOUT) :: gradv
    !
    INTEGER, POINTER :: nnr
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER, POINTER :: env_periodicity
    INTEGER, POINTER :: slab_axis
    REAL( DP ), POINTER :: alat, omega, axis_length
    REAL( DP ), DIMENSION(:), POINTER :: origin
    REAL( DP ), DIMENSION(:,:), POINTER :: axis
    !
    REAL( DP ), POINTER ::  permittivity_ms, permittivity_gcs, xstern_ms, carrier_density
    REAL( DP ), POINTER :: cion, zion, xstern_gcs, electrode_charge
    REAL( DP ) :: kbt, invkbt
    !
    TYPE( environ_gradient ), TARGET :: glocal
    REAL( DP ), DIMENSION(:,:), POINTER :: gvstern
    !
    INTEGER :: i
    REAL( DP ) :: ez, ez_ms, ez_gcs, fact, vms, vtmp
    REAL( DP ) :: vstern, const
    REAL( DP ) :: lin_k, lin_e, lin_c
    REAL( DP ) :: arg, asinh, coth, acoth
    REAL( DP ) :: f1, f2, depletion_length
    REAL( DP ) :: area, dvtmp_dx, distance
    REAL(DP) :: dipole(0:3), quadrupole(3)
    REAL(DP) :: tot_charge, tot_dipole(3), tot_quadrupole(3)
    CHARACTER( LEN = 80 ) :: sub_name = 'calc_gradvms_gcs'
    !
    CALL start_clock ('calc_gvms_gcs')
    !
    ! ... Aliases and sanity checks
    !
    IF ( .NOT. ASSOCIATED( gradv%cell, charges%cell ) ) &
         & CALL errore(sub_name,'Missmatch in domains of potential and charges',1)
    IF ( gradv % cell % nnr .NE. oned_analytic % n ) &
         & CALL errore(sub_name,'Missmatch in domains of potential and solver',1)
    cell => gradv % cell
    nnr => cell % nnr
    !
    alat => oned_analytic % alat
    omega => oned_analytic % omega
    env_periodicity => oned_analytic % d
    slab_axis => oned_analytic % axis
    axis_length => oned_analytic % size
    origin => oned_analytic % origin
    axis => oned_analytic % x
    !
    ! ... Get parameters of semiconductor to compute analytic correction
    !
    permittivity_ms => semiconductor_in%permittivity
    electrode_charge => semiconductor_in%electrode_charge
    carrier_density => semiconductor_in%carrier_density
    xstern_ms => semiconductor_in%simple%width

    cion => electrolyte % ioncctype(1) % cbulk
    zion => electrolyte % ioncctype(1) % z
    permittivity_gcs => electrolyte%permittivity
    xstern_gcs => electrolyte%boundary%simple%width
    !
    ! ... Set Boltzmann factors
    !
    kbt = semiconductor_in % temperature * k_boltzmann_ry
    invkbt = 1.D0 / kbt
    !
    IF ( env_periodicity .NE. 2 ) &
         & CALL errore(sub_name,'Option not yet implemented: 1D Poisson-Boltzmann solver only for 2D systems',1)
    !

    CALL init_environ_gradient( cell, glocal )
    gvstern => glocal % of_r
    !
    ! ... Compute multipoles of the system wrt the chosen origin
    !
    CALL compute_dipole( nnr, charges%of_r, origin, dipole, quadrupole )
    !
    tot_charge = dipole(0)
    tot_dipole = dipole(1:3)
    area = omega / axis_length
    !
    ! ... First compute the gradient of parabolic correction
    !
    fact = e2 * fpi / omega
    gvstern(slab_axis,:) = tot_dipole(slab_axis) - tot_charge * axis(1,:)
    gvstern = gvstern * fact

    !
    ! ... Compute the physical properties of the interface
    !
    ez = - tpi * e2 * tot_charge / area !/ permittivity !! the input charge density includes explicit and
                                                        !! polarization charges, so tot_charge already accounts
                                                        !! for the dielectric screening. permittivity needs not
                                                        !! to be included
    ez_gcs = - tpi * e2 * electrode_charge / area ! / permittivity
    fact = - e2 * SQRT( 8.D0 * fpi * cion * kbt / e2 )!/ permittivity )
    arg = ez/fact
    asinh = LOG(arg + SQRT( arg**2 + 1 ))
    vstern = 2.D0 * kbt / zion * asinh
    arg = vstern * 0.25D0 * invkbt * zion
    coth = ( EXP( 2.D0 * arg ) + 1.D0 ) / ( EXP( 2.D0 * arg ) - 1.D0 )
    const = coth * EXP( zion * fact * invkbt * 0.5D0 * xstern_gcs )
    !
    ! ... Compute linearized quantities
    !
    ! note that this is simplified due to same c_ion
    lin_k = SQRT(electrolyte%k2)
    lin_e = SQRT(permittivity_gcs)
    lin_c = -1.D0 * ez * lin_e / lin_k * EXP(lin_k * xstern_gcs / lin_e)
    !
    ! ... Compute the analytic gradient of potential
    !     Note that the only contribution different from the parabolic
    !     correction is in the region of the diffuse layer
    !
    ! Will first do this for the guoy chapman system
    IF ( electrolyte % linearized ) THEN
       !
       ! ... Compute some constants needed for the calculation
       !
       f1 = -1.D0 * lin_k / lin_e
       !
       DO i = 1, nnr
          !
          IF ( axis(1,i) .GE. xstern_gcs ) THEN
             !
             ! ... Linearized Gouy-Chapmann-Stern analytic solution on the outside
             !
             arg = f1 * ABS(axis(1,i))
             dvtmp_dx = lin_c * f1 * EXP( arg )
             !
             ! ... Remove source potential and add analytic one
             !
             gvstern(slab_axis,i) = -gradv % of_r(slab_axis,i) + &
                         & ( dvtmp_dx - ez ) * ABS(axis(1,i))/axis(1,i)
             !
          ENDIF
          !
       ENDDO
       !
    ELSE
       !
       ! ... Compute some constants needed for the calculation
       !
       f1 = - fact * zion * invkbt * 0.5D0
       f2 = 4.D0 * kbt / zion
       !
       DO i = 1, nnr
          !
          IF ( axis(1,i) .GE. xstern_gcs ) THEN
             !
             ! ... Gouy-Chapmann-Stern analytic solution on the outside
             !
             arg = const * EXP( ABS(axis(1,i)) * f1 )
             dvtmp_dx = f1 * f2 * arg / ( 1.D0 - arg ** 2 )
             !
             ! ... Remove source potential (linear) and add analytic one
             !
             gvstern(slab_axis,i) = -gradv % of_r(slab_axis,i) + &
                         & ( dvtmp_dx - ez ) * ABS(axis(1,i))/axis(1,i)
             !
          ENDIF
          !
       ENDDO
       !
    ENDIF
    !
    ! Now applying the mott schottky side

    ez_ms= tpi * e2 * (electrode_charge-tot_charge) / area ! / permittivity !in units of Ry/bohr
    fact = 1.D0/tpi / e2 /2.D0 /carrier_density !*permittivity
    arg = fact* (ez**2.D0)
    vms =  arg ! +kbt
    !Finds the total length of the depletion region
    depletion_length = 2.D0 *fact*ez_ms

    DO i = 1, nnr

       IF ( -axis(1,i) .GE. xstern_ms ) THEN
          distance = ABS(axis(1,i)) - xstern_ms
          !Only applies parabolic equation if still within the depletion width
          !
          IF ( distance <= depletion_length) THEN
             !
             ! ... Mott Schottky analytic solution on the outside
             !
             vtmp = 2.D0*(distance) / fact/4.D0 + ez_ms
          ELSE
             vtmp = 0.D0
          END IF
          !WRITE (environ_unit, *)"This is the axis value: ",axis(1,i)
          !WRITE (environ_unit, *) "Distance: ", distance
          !
          ! ... Remove source potential (linear) and add analytic one
          !
          gvstern(slab_axis,i) =  gvstern(slab_axis,i) + vtmp -ez!-vms ! vtmp - potential % of_r(i)
          !WRITE( environ_unit, *)"This is the vi: ",ez*distance
          !
       ENDIF
       !
    ENDDO


    gradv % of_r = gradv % of_r + gvstern

    !
    CALL destroy_environ_gradient(glocal)
    !
    CALL stop_clock ('calc_gvms_gcs')
    !
    RETURN
!---------------------------------------------------------------------------
END SUBROUTINE calc_gradvms_gcs
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
END MODULE correction_ms_gcs
!---------------------------------------------------------------------------
