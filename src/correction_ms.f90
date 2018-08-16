!
! Copyright (C) 2017 Environ group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by Q. Campbell 
!
!----------------------------------------------------------------------------
MODULE correction_ms
!----------------------------------------------------------------------------
  !
  ! ...
  !
  USE environ_types
  USE electrostatic_types
  USE environ_output,    ONLY : environ_unit
  USE environ_base,      ONLY : e2
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: calc_vms, calc_gradvms
  !
CONTAINS
!---------------------------------------------------------------------------
  SUBROUTINE calc_vms( oned_analytic, semiconductor, charges, potential )
!---------------------------------------------------------------------------
    !
    ! ... Given the total explicit charge, the value of the field at the boundary
    !     is obtained by Gauss's law
    ! (1) ez = - tpi * e2 * charge * axis_length / omega / env_static_permittivity
    !
    ! ... By integrating the Guy-Chapman model one has a relation that links the derivative
    !     of the potential (the field) to the value of the potential
    ! (2) dv/dz = - ez = fact * SINH( v(z) * zion / 2 / kbt )
    !     where
    ! (3) fact = - e2 * SQRT( 32 * pi * cd * kbt / e2 / env_static_permittivity )
    !
    ! ... By combining (1) and (2) one can derive the analytic charge from the knowledge of the potential
    !     at the boundary,
    ! (4) charge_ext = fact * env_static_permittivity * omega / axis_lenght / tpi / e2 * SINH( vskin * zion / 2 / kbt )
    !
    ! ... or one can compute the value of the potential at the interface corresponding to a certain
    !     explicit charge density,
    ! (5) vskin_analytic = 2.D0 * kbt / zion * ASINH( ez / fact )
    !
    ! ... Eventually, by integrating Eq. (2) the analytic form of the potential is found
    ! (6) vanalytic(z) = 4.D0 * kbt / zion * ACOTH( const * EXP( - z * fact * zion / kbt / 2.D0  ) )
    !     were the constant is determined by the condition that vanalytic(xskin) = vskin
    !
    IMPLICIT NONE
    !
    TYPE( oned_analytic_core ), TARGET, INTENT(IN) :: oned_analytic
    TYPE( environ_semiconductor ), TARGET, INTENT(IN) :: semiconductor
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
    REAL( DP ), POINTER :: permittivity, xstern, carrier_density
    REAL( DP ) :: kbt, invkbt
    !
    TYPE( environ_density ), TARGET :: local
    REAL( DP ), DIMENSION(:), POINTER :: v
    !
    INTEGER :: i, icount
    REAL( DP ) :: ez, fact, vms
    REAL( DP ) :: arg, const
    REAL( DP ) :: area, vtmp, vbound
    REAL(DP) :: dipole(0:3), quadrupole(3)
    REAL(DP) :: tot_charge, tot_dipole(3), tot_quadrupole(3)
    CHARACTER( LEN = 80 ) :: sub_name = 'calc_vms'
    !
    CALL start_clock ('calc_vms')
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
    permittivity => semiconductor%permittivity
    carrier_density => semiconductor%carrier_density
    xstern => semiconductor%simple%width
    
    !   convert carrier density to units of (bohr)^-3
    carrier_density = carrier_density *6.448D24
    !
    WRITE( environ_unit, * )"xstern: ",xstern
    WRITE( environ_unit, * )"carrier density: ",carrier_density
    ! ... Set Boltzmann factors
    !
    kbt = semiconductor % temperature * k_boltzmann_ry
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
    CALL compute_dipole( nnr, 1, charges%of_r, origin, dipole, quadrupole )
    !
    tot_charge = dipole(0)
    tot_dipole = dipole(1:3)
    tot_quadrupole = quadrupole
    area = omega / axis_length
    !
    ! ... First apply parabolic correction
    !
    fact = e2 * tpi / omega
    const = - pi / 3.D0 * tot_charge / axis_length * e2 - fact * tot_quadrupole(slab_axis)
    v(:) = - tot_charge * axis(1,:)**2 + 2.D0 * tot_dipole(slab_axis) * axis(1,:)
    v(:) = fact * v(:) + const
    !dv = - fact * 4.D0 * tot_dipole(slab_axis) * xstern
    !
    ! ... Compute the physical properties of the interface
    !
    !zion = ABS(zion)
    WRITE( environ_unit, *)"charge: ",tot_charge
    ez = - tpi * e2 * tot_charge / area  / permittivity !in units of Ry/bohr 
    WRITE( environ_unit, * )"electric field: ",ez
    fact = permittivity /tpi / e2 /2.D0 /carrier_density
    WRITE(  environ_unit, *)"Prefactor: ",fact
    arg = fact* (e2**2.D0)
    vms =  arg ! +kbt
    WRITE ( environ_unit, * )"vms: ",vms
    !
    vbound = 0.D0
    icount = 0
    DO i = 1, nnr
       !
       IF ( ABS(axis(1,i)) .GE. xstern ) THEN
          !
          icount = icount + 1
          vbound = vbound + potential % of_r(i) + v(i) - ez * (ABS(axis(1,i)) - xstern)
          !
       ENDIF
       !
    ENDDO
    CALL mp_sum(icount,cell%comm)
    CALL mp_sum(vbound,cell%comm)
    vbound = vbound / DBLE(icount)
    !
    ! ... Compute some constants needed for the calculation
    !
    !f1 =  fact * zion * invkbt * 0.5D0
    !
    ! ... Compute the analytic potential and charge
    !
    v = v - vbound + vms
    DO i = 1, nnr
       !IF MOD(i,15) .EQ. 0 THEN
       WRITE(environ_unit,*)axis(1,i)
       !ENDIF
       IF ( ABS(axis(1,i)) .GE. xstern ) THEN
          !
          ! ... Gouy-Chapmann-Stern analytic solution on the outside
          !
          vtmp = (ABS(axis(1,i))-xstern)**2.D0 / fact - ez*(ABS(axis(1,i))-xstern)
          !
          ! ... Remove source potential (linear) and add analytic one
          !
          v(i) =  v(i) + vtmp - vms - ez * ABS(axis(1,i)) + ez * xstern ! vtmp - potential % of_r(i)
          !
       ENDIF
       !
    ENDDO
    !
    potential % of_r = potential % of_r + v
    !
    CALL destroy_environ_density(local)
    !
    CALL stop_clock ('calc_vms')
    !




    RETURN



    !
!---------------------------------------------------------------------------
  END SUBROUTINE calc_vms
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  SUBROUTINE calc_gradvms( oned_analytic, semiconductor, charges, gradv )
!---------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( oned_analytic_core ), TARGET, INTENT(IN) :: oned_analytic
    TYPE( environ_semiconductor ), TARGET, INTENT(IN) :: semiconductor
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
    REAL( DP ), POINTER ::  permittivity, xstern, carrier_density
    REAL( DP ) :: kbt, invkbt
    !
    TYPE( environ_gradient ), TARGET :: glocal
    REAL( DP ), DIMENSION(:,:), POINTER :: gvms
    !
    INTEGER :: i
    REAL( DP ) :: ez, fact, vms
    REAL( DP ) :: arg
    REAL( DP ) :: area, dvtmp_dx
    REAL(DP) :: dipole(0:3), quadrupole(3)
    REAL(DP) :: tot_charge, tot_dipole(3), tot_quadrupole(3)
    CHARACTER( LEN = 80 ) :: sub_name = 'calc_gradvms'
    !
    CALL start_clock ('calc_gvms')
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
    permittivity => semiconductor%permittivity
    carrier_density => semiconductor%carrier_density
    xstern => semiconductor%simple%width
    !
    ! ... Set Boltzmann factors
    !
    kbt = semiconductor % temperature * k_boltzmann_ry
    invkbt = 1.D0 / kbt
    !
    IF ( env_periodicity .NE. 2 ) &
         & CALL errore(sub_name,'Option not yet implemented: 1D Poisson-Boltzmann solver only for 2D systems',1)
    !
    CALL init_environ_gradient( cell, glocal )
    gvms => glocal % of_r
    !
    ! ... Compute multipoles of the system wrt the chosen origin
    !
    CALL compute_dipole( nnr, 1, charges%of_r, origin, dipole, quadrupole )
    !
    tot_charge = dipole(0)
    tot_dipole = dipole(1:3)
    !
    ! ... First compute the gradient of parabolic correction
    !
    fact = e2 * fpi / omega
    gvms(slab_axis,:) = tot_dipole(slab_axis) - tot_charge * axis(1,:)
    gvms = gvms * fact
    !
    ! ... Compute the physical properties of the interface
    !
!    zion = ABS(zion)
!    ez = - tpi * e2 * tot_charge / area / permittivity
!    fact = - e2 * SQRT( 8.D0 * fpi * cion * kbt / e2 / permittivity )
!    arg = ez/fact
!    asinh = LOG(arg + SQRT( arg**2 + 1 ))
!    vstern = 2.D0 * kbt / zion * asinh
!    arg = vstern * 0.25D0 * invkbt * zion
!    coth = ( EXP( 2.D0 * arg ) + 1.D0 ) / ( EXP( 2.D0 * arg ) - 1.D0 )
!    const = coth * EXP( zion * fact * invkbt * 0.5D0 * xstern )
!    !
!    ! ... Compute some constants needed for the calculation
!    !
!    f1 = - fact * zion * invkbt * 0.5D0
!    f2 = 4.D0 * kbt / zion
!    !
!    ! ... Compute the analytic gradient of potential
!    !     Note that the only contribution different from the parabolic
!    !     correction is in the region of the diffuse layer
!    !
!    DO i = 1, nnr
!       !
!       IF ( ABS(axis(1,i)) .GE. xstern ) THEN
!          !
!          ! ... Gouy-Chapmann-Stern analytic solution on the outside
!          !
!          arg = const * EXP( ABS(axis(1,i)) * f1 )
!          dvtmp_dx = f1 * f2 * arg / ( 1.D0 - arg ** 2 )
!          !
!          ! ... Remove source potential (linear) and add analytic one
!          !
!          gvstern(slab_axis,i) =  gvstern(slab_axis,i) + dvtmp_dx - ez
!          !
!       ENDIF
!       !
!    ENDDO
    !
    ! Seems like 
    gradv % of_r = gradv % of_r + gvms
    !
    CALL destroy_environ_gradient(glocal)
    !
    CALL stop_clock ('calc_gvms')
    !
    RETURN
!---------------------------------------------------------------------------
  END SUBROUTINE calc_gradvms
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
END MODULE correction_ms
!---------------------------------------------------------------------------

