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
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------
MODULE correction_gcs
!----------------------------------------------------------------------------
  !
  ! ...
  !
  USE modules_constants, ONLY : e2, k_boltzmann_ry, pi, tpi, fpi
  USE environ_types
  USE core_types
  USE environ_output
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: calc_vgcs, calc_gradvgcs
  !
CONTAINS
    !  Subroutine: calc_vgcs
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
    !> Function that calculates the Gouy-Chapman correction and adds it to the
    !! potential.
    !! 
    !! Procedure is as follows:
    !! -# execute parabolic correction
    !! -# shift to fit open boundary conditions
    !! -# add analytic solution outside stern boundary
    !!
    !! The following outlines the Gouy-Chapman correction for the non-linear case.
    !!
    !! Given the total explicit charge, the value of the field at the boundary
    !! is obtained by Gauss's law
    !! \f[
    !!    E_z = \frac{-2\pi q}{A \epsilon_0}
    !! \f]
    !!
    !! By integrating the Gouy-Chapman model one has a relation that links the derivative
    !! of the potential (the field) to the value of the potential
    !! \f[
    !!    \frac{dv}{dz} = -E_z = f\sinh( v(x) z_d / 2k_BT )
    !! \f]
    !! where \f$z_d\f$ is the electrolyte charge and 
    !! \f[
    !!    f = \sqrt{\frac{ 32\pi c_d k_BT}{\epsilon_0}}
    !! \f]
    !!
    !! By combining the equations for the electric field, one can derive the 
    !! analytic charge from the knowledge of the potential at the boundary,
    !! \f[
    !!    Q_{ext} = \frac{f\epsilon_0 A}{2\pi}\sinh\left(\frac{v_s z_d}
    !!    {2k_B T}\right)
    !! \f]
    !! where \f$v_s\f$ is vstern
    !!
    !! or one can compute the value of the potential at the interface corresponding 
    !! to a certain explicit charge density,
    !! \f[
    !!    v_a = \frac{2 k_BT}{z_d} \sinh^{-1}(E_z / f)
    !! \f]
    !! where \f$v_a\f$ is the analytic value of vstern.
    !!
    !! Eventually, by integrating the equation for potential, the analytic form 
    !! of the potential is found
    !! \f[
    !!    v(x) = \frac{4k_BT}{z_d}\coth(c \exp( -x f z_d / 2k_BT  ) )
    !! \f]
    !! where c is determined by the condition that v(xstern) = vstern
    !!
    !! For the linearized case, the analytic solution is calculated as follows:
    !! 
    !! The electric field can be calculated as in the non-linear case
    !! 
    !! The equation for the potential is
    !! \f[
    !!    v(x) = c \exp(-k\left|x\right|/\sqrt{\epsilon_0})
    !! \f]
    !!
    !! where c is determined by relating the derivative of the potential with the
    !! value of the electric field at xstern, and k is a quantity that depends on the
    !! electrolyte.
!---------------------------------------------------------------------------
  SUBROUTINE calc_vgcs( oned_analytic, electrolyte, charges, potential )
!---------------------------------------------------------------------------
    IMPLICIT NONE
    !
    TYPE( oned_analytic_core ), TARGET, INTENT(IN) :: oned_analytic
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
    REAL( DP ), POINTER :: cion, permittivity, xstern
    REAL( DP ) :: kbt, invkbt
    !
    TYPE( environ_density ), TARGET :: local
    REAL( DP ), DIMENSION(:), POINTER :: v
    !
    INTEGER :: i, icount
    REAL( DP ) :: ez, fact, vstern, const
    REAL( DP ) :: dv, vbound, zion
    REAL( DP ) :: arg, asinh, coth, acoth
    REAL( DP ) :: f1, f2
    REAL( DP ) :: area, vtmp
    REAL( DP ) :: lin_k, lin_e, lin_c
    REAL(DP) :: dipole(0:3), quadrupole(3)
    REAL(DP) :: tot_charge, tot_dipole(3), tot_quadrupole(3)
    CHARACTER( LEN = 80 ) :: sub_name = 'calc_vgcs'
    !
    CALL start_clock ('calc_vgcs')
    !
    ! ... Aliases and sanity checks
    !
    IF ( .NOT. ASSOCIATED( potential%cell, charges%cell ) ) &
         & CALL errore(sub_name,'Missmatch in domains of potential and charges',1)
    IF ( .NOT. ASSOCIATED( potential%cell, oned_analytic%cell ) ) &
         & CALL errore(sub_name,'Missmatch in domains of potential and solver',1)
    cell => potential % cell
    nnr => cell % nnr
    alat => cell % alat
    omega => cell % omega
    !
    env_periodicity => oned_analytic % d
    slab_axis => oned_analytic % axis
    axis_length => oned_analytic % size
    origin => oned_analytic % origin
    axis => oned_analytic % x
    !
    ! ... Get parameters of electrolyte to compute analytic correction
    !
    IF ( electrolyte % ntyp .NE. 2 ) &
         & CALL errore(sub_name,'Unexpected number of counterionic species, different from two',1)
    cion => electrolyte % ioncctype(1) % cbulk
    zion = ABS(electrolyte % ioncctype(1) % z)
    permittivity => electrolyte%permittivity
    xstern => electrolyte%boundary%simple%width
    !
    ! ... Set Boltzmann factors
    !
    kbt = electrolyte % temperature * k_boltzmann_ry
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
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X QE-6.2.X QE-6.3.X
!    CALL compute_dipole( nnr, 1, charges%of_r, origin, dipole, quadrupole )
! Compatible with QE-6.4.X, and QE-GIT
    CALL compute_dipole( nnr, charges%of_r, origin, dipole, quadrupole )
! END BACKWARD COMPATIBILITY
    !
    tot_charge = dipole(0)
    tot_dipole = dipole(1:3)
    tot_quadrupole = quadrupole
    area = omega / axis_length
    !
    ! ... First apply parabolic correction
    fact = e2 * tpi / omega
    const = - pi / 3.D0 * tot_charge / axis_length * e2 - fact * tot_quadrupole(slab_axis)
    v(:) = - tot_charge * axis(1,:)**2 + 2.D0 * tot_dipole(slab_axis) * axis(1,:)
    v = fact * v + const
    !
    ! ... Compute the physical properties of the interface
    !
    ez = - tpi * e2 * tot_charge / area !/ permittivity !! the input charge density includes explicit and
                                                        !! polarization charges, so tot_charge already accounts
                                                        !! for the dielectric screening. permittivity needs not
                                                        !! to be included
    fact = - e2 * SQRT( 8.D0 * fpi * cion * kbt / e2 / permittivity )
    arg = ez/fact
    asinh = LOG(arg + SQRT( arg**2 + 1 ))
    vstern = 2.D0 * kbt / zion * asinh
    arg = vstern * 0.25D0 * invkbt * zion
    coth = ( EXP( 2.D0 * arg ) + 1.D0 ) / ( EXP( 2.D0 * arg ) - 1.D0 )
    const = coth * EXP( zion * fact * invkbt * 0.5D0 * xstern )
    !
    ! ... Compute linearized quantities
    !
    ! note that this is simplified due to same c_ion
    lin_k = SQRT(electrolyte%k2)
    lin_e = SQRT(permittivity)
    lin_c = -1.D0 * ez * lin_e / lin_k * EXP(lin_k * xstern / lin_e)
    IF ( electrolyte % linearized ) THEN
       vstern = lin_c * EXP(-1.D0 * lin_k * xstern / lin_e)
    ENDIF
    !
    ! ... Compute value of the reference potential at the boundary with electrolyte
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
    v = v - vbound + vstern
    !
    ! ... Compute some constants needed for the calculation
    !
    f1 = - fact * zion * invkbt * 0.5D0
    f2 = 4.D0 * kbt / zion
    !
    ! ... Compute the analytic potential
    !
    IF ( electrolyte % linearized ) THEN
       !
       DO i = 1, nnr
          !
          IF ( ABS(axis(1,i)) .GE. xstern ) THEN
             !
             ! ... Linearized Gouy-Chapmann-Stern analytic solution on the outside
             !
             vtmp = lin_c * EXP( -1.D0 * lin_k * ABS(axis(1,i)) / lin_e)
             !
             ! ... Remove source potential and add analytic one
             !
             !v(i) = v(i) + vtmp - vstern - ez * ABS(axis(1,i)) + ez * xstern
             v(i) = vtmp - potential % of_r(i)
             !
          ENDIF
          !
       ENDDO
       !
    ELSE
       !
       DO i = 1, nnr
          !
          IF ( ABS(axis(1,i)) .GE. xstern ) THEN
             !
             ! ... Gouy-Chapmann-Stern analytic solution on the outside
             !
             arg = const * EXP( ABS(axis(1,i)) * f1 )
             IF ( ABS(arg) .GT. 1.D0 ) THEN
                acoth = 0.5D0 * LOG( (arg + 1.D0) / (arg - 1.D0) )
             ELSE
                acoth = 0.D0
             ENDIF
             vtmp =  f2 * acoth
             !
             ! ... Remove source potential and add analytic one
             !
!             v(i) =  v(i) + vtmp - vstern - ez * ABS(axis(1,i)) + ez * xstern
             v(i) = vtmp - potential % of_r(i)
             !
          ENDIF
          !
       ENDDO
    ENDIF
    !
    potential % of_r = potential % of_r + v
    !
    CALL destroy_environ_density(local)
    !
    CALL stop_clock ('calc_vgcs')
    !
    RETURN
    !
!---------------------------------------------------------------------------
  END SUBROUTINE calc_vgcs
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  SUBROUTINE calc_gradvgcs( oned_analytic, electrolyte, charges, gradv )
!---------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( oned_analytic_core ), TARGET, INTENT(IN) :: oned_analytic
    TYPE( environ_electrolyte ), TARGET, INTENT(IN) :: electrolyte
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
    REAL( DP ), POINTER :: cion, permittivity, xstern
    REAL( DP ) :: zion, kbt, invkbt
    !
    TYPE( environ_gradient ), TARGET :: glocal
    REAL( DP ), DIMENSION(:,:), POINTER :: gvstern
    !
    INTEGER :: i
    REAL( DP ) :: ez, fact, vstern, const
    REAL( DP ) :: arg, asinh, coth, acoth
    REAL( DP ) :: lin_k, lin_e, lin_c
    REAL( DP ) :: f1, f2
    REAL( DP ) :: area, dvtmp_dx
    REAL(DP) :: dipole(0:3), quadrupole(3)
    REAL(DP) :: tot_charge, tot_dipole(3), tot_quadrupole(3)
    CHARACTER( LEN = 80 ) :: sub_name = 'calc_gradvgcs'
    !
    CALL start_clock ('calc_gvst')
    !
    ! ... Aliases and sanity checks
    !
    IF ( .NOT. ASSOCIATED( gradv%cell, charges%cell ) ) &
         & CALL errore(sub_name,'Missmatch in domains of potential and charges',1)
    IF ( .NOT. ASSOCIATED ( gradv%cell, oned_analytic%cell ) ) &
         & CALL errore(sub_name,'Missmatch in domains of potential and solver',1)
    cell => gradv % cell
    nnr => cell % nnr
    alat => cell % alat
    omega => cell % omega
    !
    env_periodicity => oned_analytic % d
    slab_axis => oned_analytic % axis
    axis_length => oned_analytic % size
    origin => oned_analytic % origin
    axis => oned_analytic % x
    !
    ! ... Get parameters of electrolyte to compute analytic correction
    !
    IF ( electrolyte % ntyp .NE. 2 ) &
         & CALL errore(sub_name,'Unexpected number of counterionic species, different from two',1)
    cion => electrolyte % ioncctype(1) % cbulk
    zion = ABS( electrolyte % ioncctype(1) % z )
    permittivity => electrolyte%permittivity
    xstern => electrolyte%boundary%simple%width
    !
    ! ... Set Boltzmann factors
    !
    kbt = electrolyte % temperature * k_boltzmann_ry
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
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X QE-6.2.X QE-6.3.X
!    CALL compute_dipole( nnr, 1, charges%of_r, origin, dipole, quadrupole )
! Compatible with QE-6.4.X, and QE-GIT
    CALL compute_dipole( nnr, charges%of_r, origin, dipole, quadrupole )
! END BACKWARD COMPATIBILITY
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
    fact = - e2 * SQRT( 8.D0 * fpi * cion * kbt / e2 / permittivity )
    arg = ez/fact
    asinh = LOG(arg + SQRT( arg**2 + 1 ))
    vstern = 2.D0 * kbt / zion * asinh
    arg = vstern * 0.25D0 * invkbt * zion
    coth = ( EXP( 2.D0 * arg ) + 1.D0 ) / ( EXP( 2.D0 * arg ) - 1.D0 )
    const = coth * EXP( zion * fact * invkbt * 0.5D0 * xstern )
    !
    ! ... Compute linearized quantities
    !
    ! note that this is simplified due to same c_ion
    lin_k = SQRT(electrolyte%k2)
    lin_e = SQRT(permittivity)
    lin_c = -1.D0 * ez * lin_e / lin_k * EXP(lin_k * xstern / lin_e)
    !
    ! ... Compute the analytic gradient of potential
    !     Note that the only contribution different from the parabolic
    !     correction is in the region of the diffuse layer
    !
    IF ( electrolyte % linearized ) THEN
       !
       ! ... Compute some constants needed for the calculation
       !
       f1 = -1.D0 * lin_k / lin_e
       !
       DO i = 1, nnr
          !
          IF ( ABS(axis(1,i)) .GE. xstern ) THEN
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
          IF ( ABS(axis(1,i)) .GE. xstern ) THEN
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
    gradv % of_r = gradv % of_r + gvstern
    !
    CALL destroy_environ_gradient(glocal)
    !
    CALL stop_clock ('calc_gvst')
    !
    RETURN
!---------------------------------------------------------------------------
  END SUBROUTINE calc_gradvgcs
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
END MODULE correction_gcs
!---------------------------------------------------------------------------
