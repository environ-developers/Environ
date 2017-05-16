!--------------------------------------------------------------------
MODULE electrolyte_utils
!--------------------------------------------------------------------

  USE environ_types

  USE dielectric
  USE boundary
  
!  USE constants,      ONLY: k_boltzmann_ry, pi, tpi, fpi
!  USE io_global,      ONLY: stdout
!  USE mp,             ONLY: mp_sum
!  USE mp_bands,       ONLY: intra_bgrp_comm
!  USE environ_cell,   ONLY: domega, alat, omega, ntot, at
!  USE environ_ions,   ONLY: avg_pos, rhoions
!  USE environ_base,   ONLY: verbose, e2, ir_end, environ_unit,           &
!                            env_ioncc_level, nrep, cion, cionmax, zion,  &
!                            solvent_temperature, rhoioncc, rhopolcc,     &
!                            env_static_permittivity, rhopol, slab_axis,  &
!                            rhomin, rhopb
!  USE generate_function, ONLY: planar_average
!  USE periodic,       ONLY: calc_vperiodic
!  USE environ_debug,  ONLY: write_cube
!  USE generate_f_of_rho, ONLY: epsilonfunct

  IMPLICIT NONE

!  INTEGER :: naxis, shift, cell_max, cell_min, n1d, nr1, nr2, nr3
!  REAL( DP ) :: area, darea, axis_length, dx, xstern, deltastern
!
!  REAL( DP ) :: kbt, invkbt, fsw
!
!  REAL( DP ), ALLOCATABLE :: axis(:), switch(:), step(:)
!  REAL( DP ), ALLOCATABLE :: axis1d(:), switch1d(:), step1d(:)
!
!  REAL( DP ) :: charge_fix, dipole_fix, dv_fix, vfix_avg, vfix1d_avg
!  REAL( DP ), ALLOCATABLE :: rhofix(:), vfix(:), vzero(:), gradvzero(:,:), rhozero(:), laplvzero(:)
!  REAL( DP ), ALLOCATABLE :: rhofix1d(:), vfix1d(:)
!  REAL( DP ), ALLOCATABLE :: czero1d(:,:), rhozero1d(:), vzero1d(:)
!
!  INTEGER :: nnrlocal, nspinlocal
!  REAL( DP ), ALLOCATABLE :: rhoinlocal(:), rhooutlocal(:), vtotlocal(:,:)
!
!  INTEGER :: ifdtype, nfdpoint, ncfd, ncfd2
!  INTEGER, ALLOCATABLE :: icfd(:), icfd2(:)
!  REAL( DP ), ALLOCATABLE :: eps1d(:), deps1d(:)
!
!  INTEGER :: nsp
!  REAL( DP ), ALLOCATABLE :: cb(:), z(:), mu(:)

  PRIVATE

  PUBLIC :: create_environ_electrolyte, init_environ_electrolyte_first, &
       & init_environ_electrolyte_second, destroy_environ_electrolyte

CONTAINS

  SUBROUTINE create_environ_electrolyte( electrolyte )

    IMPLICIT NONE

    TYPE( environ_electrolyte ), INTENT(INOUT) :: electrolyte

    CHARACTER( LEN=80 ) :: sub_name = 'create_environ_electrolyte'
    CHARACTER( LEN=80 ) :: label = 'electrolyte'

    CALL create_environ_boundary( electrolyte%boundary )
    CALL create_environ_density( electrolyte%density, label )

    IF ( ALLOCATED( electrolyte%ioncctype ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)

    RETURN

  END SUBROUTINE create_environ_electrolyte

  SUBROUTINE init_environ_electrolyte_first( ntyp, mode, stype, rhomax, rhomin, &
             & tbeta, distance, spread, alpha, softness, electrons, ions, &
             & temperature, cbulk, cmax, radius, z, electrolyte )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ntyp, stype
    CHARACTER( LEN=80 ), INTENT(IN) :: mode
    REAL( DP ), INTENT(IN) :: rhomax, rhomin, tbeta, distance, spread, alpha, softness, temperature
    REAL( DP ), DIMENSION(ntyp), INTENT(IN) :: cbulk, cmax, radius, z
    TYPE( environ_electrons ), INTENT(IN) :: electrons
    TYPE( environ_ions ), INTENT(IN) :: ions
    TYPE( environ_electrolyte ), INTENT(INOUT) :: electrolyte

    INTEGER :: ityp
    REAL( DP ) :: neutral
    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_electrolyte_first'

    electrolyte%ntyp = ntyp

    electrolyte%temperature = temperature

    CALL init_environ_boundary_first( mode, stype, rhomax, rhomin, tbeta, .FALSE., &
         & 0.D0, alpha, softness, electrons, ions, electrolyte%boundary )

    ALLOCATE( electrolyte%ioncctype(ntyp) )

    neutral = 0.D0
    DO ityp = 1, ntyp
       ! If the radius is provided in input, compute cmax from it
       electrolyte%ioncctype(ityp)%cmax = cmax(ityp) * bohr_radius_si**3 / amu_si
       IF ( cmax(ityp) .EQ. 0.D0 .AND. radius(ityp) .GT. 0.D0 ) &
            & electrolyte%ioncctype(ityp)%cmax  = 0.64D0 * 3.D0 / fpi / radius(ityp)**3
       ! Double check that the bulk and max concentrations in input are compatible
       IF ( cbulk(ityp) .GT. 0.D0 .AND. cmax(ityp) .LT. cbulk(ityp) ) &
            & call errore (sub_name,'cmax should be at least greater than cbulk',1)
       electrolyte%ioncctype(ityp)%cbulk = cbulk(ityp) * bohr_radius_si**3 / amu_si
       electrolyte%ioncctype(ityp)%radius = radius(ityp)
       electrolyte%ioncctype(ityp)%z = z(ityp)
       neutral = neutral + cbulk(ityp)*z(ityp)
    END DO

    IF ( neutral .GT. 1.D-8 ) CALL errore(sub_name,'Bulk electrolyte is not neutral',1)

    RETURN

  END SUBROUTINE init_environ_electrolyte_first

  SUBROUTINE init_environ_electrolyte_second( cell, electrolyte )

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_electrolyte ), INTENT(INOUT) :: electrolyte

    CALL init_environ_boundary_second( cell, electrolyte%boundary )

    CALL init_environ_density( cell, electrolyte%density )

    RETURN

  END SUBROUTINE init_environ_electrolyte_second

  SUBROUTINE destroy_environ_electrolyte( lflag, electrolyte )

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_electrolyte ), INTENT(INOUT) :: electrolyte
    CHARACTER( LEN=80 ) :: sub_name = 'destroy_environ_electrolyte'

    IF ( lflag ) THEN

       ! These components were allocated first, destroy only if lflag = .TRUE.

       IF ( .NOT. ALLOCATED( electrolyte%ioncctype ) ) &
            & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
       DEALLOCATE( electrolyte%ioncctype )

    ENDIF

    CALL destroy_environ_boundary( lflag, electrolyte%boundary )
    CALL destroy_environ_density( electrolyte%density )

    RETURN

  END SUBROUTINE destroy_environ_electrolyte

!!--------------------------------------------------------------------
!  SUBROUTINE ioncc_initbase( ifdtype_, nfdpoint_, nnr )
!!--------------------------------------------------------------------
!    !
!    USE fd_gradient,       ONLY : init_fd_gradient, init_fd_laplacian
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: ifdtype_, nfdpoint_, nnr
!    !
!    INTEGER :: i
!    !
!    IF (ALLOCATED(axis)) DEALLOCATE(axis)
!    ALLOCATE(axis(nnr))
!    axis = 0.D0
!    !
!    IF (ALLOCATED(switch)) DEALLOCATE(switch)
!    ALLOCATE(switch(nnr))
!    switch = 0.D0
!    !
!    IF (ALLOCATED(step)) DEALLOCATE(step)
!    ALLOCATE(step(nnr))
!    step = 0.D0
!    !
!    kbt = k_boltzmann_ry * solvent_temperature
!    invkbt = 1.D0/kbt
!    IF ( verbose .GE. 2 ) THEN
!       WRITE(environ_unit,*)'kbt =',kbt
!       WRITE(environ_unit,*)'invkbt = ',invkbt
!       WRITE(environ_unit,*)'cionmax = ',cionmax
!       WRITE(environ_unit,*)'kbt * cionmax = ',kbt*cionmax
!    ENDIF
!    !
!    fsw = LOG( rhomin / rhopb )
!    !
!    ! ... Set the finite difference derivative coefficients
!    !
!    ifdtype = ifdtype_
!    nfdpoint = nfdpoint_
!    IF (ALLOCATED(icfd)) DEALLOCATE(icfd)
!    ALLOCATE(icfd(-nfdpoint:nfdpoint))
!    CALL init_fd_gradient( ifdtype, nfdpoint, ncfd, icfd )
!    !
!    ! ... Set the central finite difference second derivative coefficients
!    !
!    IF (ALLOCATED(icfd2)) DEALLOCATE(icfd2)
!    ALLOCATE(icfd2(-nfdpoint:nfdpoint))
!    CALL init_fd_laplacian( 1, nfdpoint, ncfd2, icfd2 )
!    !
!    nsp = 2
!    ALLOCATE(z(nsp)) ! ionic charge
!    ALLOCATE(cb(nsp)) ! bulk concentration
!    z(1)  =  zion ! charge
!    z(2)  = -zion
!    cb(1) = cion ! bulk concentration
!    cb(2) = cion
!    !
!    ALLOCATE(mu(nsp))
!    CALL calc_muioncc(nsp,cb,mu)
!    !
!    WRITE(environ_unit,*)nsp,cb,mu
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE ioncc_initbase
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE ioncc_initcell( nnr, n1, n2, n3, at )
!!--------------------------------------------------------------------
!    !
!    USE generate_function, ONLY: generate_axis
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: nnr, n1, n2, n3
!    REAL(DP), INTENT(IN) :: at(3,3)
!    !
!    INTEGER :: ia, irep, imin, imax
!    REAL(DP) :: xmin, xmax, arg
!    !
!    nr1 = n1
!    nr2 = n2
!    nr3 = n3
!    !
!    avg_pos = 0.D0
!    !
!    ! ... Set the properties of the axis and generate it
!    !
!    axis_length = at( slab_axis, slab_axis ) * alat
!    SELECT CASE (slab_axis)
!    CASE (1)
!      naxis = n1
!    CASE (2)
!      naxis = n2
!    CASE (3)
!      naxis = n3
!    END SELECT
!    CALL generate_axis( nnr, slab_axis, avg_pos, axis )
!    IF ( verbose .GE. 2 ) CALL write_cube( nnr, axis, 'axis.cube' )
!    !
!    ! ... Set the main indexes to convert 3D quantities in the cell
!    !     into 1D quantities in extended cell (unit cell plus 2*nrep replicas)
!    !
!    n1d = (2*nrep+1)*naxis
!    cell_min = nrep*naxis + 1
!    cell_max = (nrep+1)*naxis
!    dx = axis_length/DBLE(naxis)
!    area = omega/axis_length
!    darea = dx*area
!    shift = naxis/2 - NINT(avg_pos(slab_axis)*alat/dx)
!    IF ( verbose .GE. 1 ) THEN
!      WRITE(environ_unit,8301)naxis,naxis/2,nrep,n1d
!      WRITE(environ_unit,8302)cell_min,cell_max,dx,area
!      WRITE(environ_unit,8303)darea,avg_pos(slab_axis)*alat,shift
!    ENDIF
!    !
!    ! ... Properties of the ioncc density (only temporary)
!    !
!    xstern = 3.D0 ! only temporary, PUT THIS IN INPUT KEYWORDS
!    deltastern = 1.D0 ! only temporary, PUT THIS IN INPUT KEYWORDS
!    !
!    ! ... Generate the 1D axis: first compute it in the original cell
!    !
!    IF (ALLOCATED(axis1d)) DEALLOCATE(axis1d)
!    ALLOCATE(axis1d(n1d))
!    axis1d = 0.D0
!    CALL planar_average( nnr, naxis, slab_axis, shift, .false., axis, axis1d(cell_min:cell_max) )
!    !
!    ! ... then generate the extended replicas
!    !
!    DO irep = 1, 2*nrep+1
!      imin = (irep-1)*naxis + 1
!      imax = imin + naxis - 1
!      axis1d(imin:imax) = axis1d(cell_min:cell_max) + ( irep - 2*nrep ) * axis_length
!    ENDDO
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, axis1d, 'axis1d.dat' )
!    !
!    ! ... Generate the step function, to allow finite step in potential
!    !
!    CALL generate_step(nnr,axis,step)
!    IF (ALLOCATED(step1d)) DEALLOCATE(step1d)
!    ALLOCATE(step1d(n1d))
!    CALL generate_step(n1d,axis1d,step1d)
!    !
!    ! ... Generate the 1D switching function, for Level 2 or higher ioncc density
!    !
!    IF (ALLOCATED(switch1d)) DEALLOCATE(switch1d)
!    ALLOCATE(switch1d(n1d))
!    CALL generate_switch(n1d,xstern,deltastern,axis1d,switch1d)
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, switch1d, 'switch1d.dat' )
!    !
!    RETURN
!    !
!8301 FORMAT(1X,'naxis = ',i4,' naxis/2 = ',i4,' nrep = ',i3,' n1d = ',i4)
!8302 FORMAT(1X,'cell_min = ',i4,' cell_max = ',i4,' dx = ',F14.6,' area = ',F14.6)
!8303 FORMAT(1X,' darea = ',F16.10,' avg_pos = ',F14.6,' shift = ',i5)
!!--------------------------------------------------------------------
!  END SUBROUTINE ioncc_initcell
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE ioncc_clean( )
!!--------------------------------------------------------------------
!    !
!    ! ... Clean up of local variables
!    !
!    IMPLICIT NONE
!    !
!    IF (ALLOCATED(axis)) DEALLOCATE(axis)
!    IF (ALLOCATED(switch)) DEALLOCATE(switch)
!    IF (ALLOCATED(step)) DEALLOCATE(step)
!    IF (ALLOCATED(axis1d)) DEALLOCATE(axis1d)
!    IF (ALLOCATED(switch1d)) DEALLOCATE(switch1d)
!    IF (ALLOCATED(step1d)) DEALLOCATE(step1d)
!    !
!    IF (ALLOCATED(icfd)) DEALLOCATE(icfd)
!    IF (ALLOCATED(icfd2)) DEALLOCATE(icfd2)
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE ioncc_clean
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE calc_vioncc(  nnr, nspin, rhoelec, vioncc )
!!--------------------------------------------------------------------
!    !
!    USE generate_function, ONLY: generate_gaussian
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: nnr, nspin
!    REAL( DP ), INTENT(IN) :: rhoelec(nnr)
!    REAL( DP ), INTENT(OUT) :: vioncc(nnr)
!    !
!    INTEGER :: ia, ir, irep, iter
!    REAL( DP ), ALLOCATABLE :: rhotot(:)
!    REAL( DP ), ALLOCATABLE :: vzero1d(:), vioncc1d(:), switch1dtmp(:)
!    REAL( DP ), ALLOCATABLE :: rhostern1d(:), vstern1d(:), cstern1d(:,:)
!    REAL( DP ), ALLOCATABLE :: rhoioncc1d(:), rhotot1d(:), cioncc1d(:,:)
!    !
!    REAL( DP ) :: charge, xtmp, spread, pos(3)
!    INTEGER :: dim
!    !
!    !
!    CALL start_clock( 'get_ioncc' )
!    !
!    ! ... Initialize the potential
!    !
!    vioncc = 0.D0
!    IF ( env_ioncc_level .LT. 1 ) RETURN
!    ALLOCATE( vioncc1d(n1d) )
!    vioncc1d = 0.D0
!    !
!    ! ... Set the total density, only exclude the ioncc densities
!    !
!    ALLOCATE( rhotot(nnr) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! TESTING ONLY
!    !
!    ! ... Test desity, we build a uniform charge density thick 2.D0 a.u. around avg_pos
!    !     The total chage of the analytic system is 1 e
!    !
!    rhotot = 0.D0
!    charge = -1.d0
!    spread = 0.5D0
!    dim = 2
!    pos = avg_pos
!    CALL generate_gaussian( nnr, dim, slab_axis, charge, spread, pos, rhotot )
!    charge = SUM(rhotot)*domega
!    CALL mp_sum( charge, intra_bgrp_comm )
!    WRITE(environ_unit,*)'input charge = ',charge
!!!!OLD-POSSIBLY-DIFFICULT-TO-CONVERGE-WITH-FFT    DO ir = 1, ir_end
!!!!OLD-POSSIBLY-DIFFICULT-TO-CONVERGE-WITH-FFT      IF ( ABS(axis(ir)) .LT. 1.D0 ) rhotot(ir) = -1.D0
!!!!OLD-POSSIBLY-DIFFICULT-TO-CONVERGE-WITH-FFT    ENDDO
!!!!OLD-POSSIBLY-DIFFICULT-TO-CONVERGE-WITH-FFT    charge = SUM(rhotot)*domega
!!!!OLD-POSSIBLY-DIFFICULT-TO-CONVERGE-WITH-FFT    CALL mp_sum( charge, intra_bgrp_comm )
!!!!OLD-POSSIBLY-DIFFICULT-TO-CONVERGE-WITH-FFT    rhotot = rhotot / ABS(charge)
!!!! TESTING ONLY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    rhotot = rhoelec + rhoions
!!!!    IF ( env_static_permittivity .GT. 1.D0 ) rhotot = rhotot + rhopol
!    IF ( verbose .GE. 3 ) CALL write_cube( nnr, rhotot, 'rhotot.cube' )
!    !
!    ! ... Compute and store some basic quantities of the source charge density
!    !
!    CALL setfix( nnr, nspin, n1d, rhotot )
!    !
!    ! ... Generate the 3D switching function, for Level 3 or higher ioncc density
!    !
!    IF ( env_ioncc_level .GE. 3 ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! TESTING ONLY !!!!!!!!!!!
!      ALLOCATE(switch1dtmp(n1d))
!      xtmp = xstern - 0.5D0
!      CALL generate_switch(n1d, xtmp, deltastern, axis1d, switch1dtmp )
!      IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, switch1dtmp, 'switch1dtmp.dat' )
!      switch = 0.D0
!      CALL planar_average( nnr, naxis, slab_axis, shift, .true., switch, switch1dtmp(cell_min:cell_max) )
!!!! TESTING ONLY !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!      DO ir = 1, nnr
!!!!        switch( ir ) = ( 1.D0 - switch( ir ) ) *    &
!!!!          ( epsilonfunct( rhoelec( ir ), rhomin, rhopb, fsw, 2.D0, 1 ) - 1.D0 )
!!!!      END DO
!    ENDIF
!    IF ( verbose .GE. 3 ) CALL write_cube( nnr, switch, 'switch.cube' )
!    !
!    ! ... Compute Level 1 ioncc (analytic 1D diffuse layer with sharp boundary)
!    !
!    ALLOCATE( rhostern1d(n1d) )
!    rhostern1d = 0.D0
!    CALL calc_rhostern( n1d, charge_fix, dipole_fix, xstern, axis1d, step1d, vfix1d, vioncc1d, rhostern1d )
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, rhostern1d, 'rhostern1d.dat' )
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, vioncc1d, 'vstern1d.dat' )
!    ALLOCATE( cstern1d(nsp,n1d) )
!    cstern1d = 0.D0
!    CALL calc_cstern( n1d, nsp, charge_fix, dipole_fix, xstern, axis1d, step1d, vfix1d, vioncc1d, cstern1d )
!    IF ( ALLOCATED(switch1dtmp) ) DEALLOCATE(switch1dtmp)
!    ALLOCATE(switch1dtmp(n1d))
!    switch1dtmp = 1.D0
!    CALL sum_cioncc( n1d, nsp, switch1dtmp, z, cstern1d, rhostern1d )
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, rhostern1d, 'rhostern1d_new.dat' )
!    !
!    ALLOCATE( rhoioncc1d(n1d) )
!    rhoioncc1d = 0.D0
!    ALLOCATE( cioncc1d(nsp,n1d) )
!    cioncc1d = 0.D0
!    IF ( env_ioncc_level .EQ. 1 ) THEN
!      CALL planar_average( nnr, naxis, slab_axis, shift, .true., vioncc, vioncc1d(cell_min:cell_max) )
!      IF ( verbose .GE. 2 ) CALL write_cube( nnr, vioncc, 'vioncc.cube')
!      CALL planar_average( nnr, naxis, slab_axis, shift, .true., rhoioncc, rhostern1d(cell_min:cell_max) )
!      IF ( verbose .GE. 2 ) CALL write_cube( nnr, rhoioncc, 'rhoioncc.cube')
!    !
!    ! ... Compute Level 2 ioncc (self-consistent 1D diffuse layer with switching function)
!    !
!    ELSE IF ( env_ioncc_level .EQ. 2 ) THEN
!      CALL calc_rhoioncc1d_lbfgs( n1d, nsp, vioncc1d, cioncc1d )
!!      CALL calc_vioncc1d_lbfgs( n1d, nsp, vioncc1d, cioncc1d )
!!      CALL calc_fullioncc1d_lbfgs( n1d, nsp, vioncc1d, cioncc1d )
!      CALL sum_cioncc( n1d, nsp, switch1d, z, cioncc1d, rhoioncc1d )
!      IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, rhoioncc1d, 'rhoioncc1d.dat' )
!      IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, vioncc1d, 'vioncc1d.dat' )
!      CALL planar_average( nnr, naxis, slab_axis, shift, .true., vioncc, vioncc1d(cell_min:cell_max) )
!      CALL planar_average( nnr, naxis, slab_axis, shift, .true., rhoioncc, rhoioncc1d(cell_min:cell_max) )
!    !
!    ! ... Compute Level 3 ioncc (self-consistent internal 3D diffuse layer with electrons-based
!    !      switching function plus self-consistent external 1D diffuse layer)
!    !
!    ELSE IF ( env_ioncc_level .EQ. 3 ) THEN
!      CALL calc_vioncc1d_lbfgs( n1d, nsp, vioncc1d, cioncc1d )
!      IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, rhoioncc1d, 'rhoioncc1d.dat' )
!      IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, vioncc1d, 'vioncc1d.dat' )
!      CALL planar_average( nnr, naxis, slab_axis, shift, .true., vioncc, vioncc1d(cell_min:cell_max) )
!      vioncc = vioncc + vfix
!      CALL calc_vioncc_lbfgs( nnr, nsp, rhoioncc, vioncc )
!    ELSE
!      WRITE(stdout,*)'ERROR: specified ioncc level is not implemented'
!    END IF
!    !
!    CALL stop_clock( 'get_ioncc' )
!    !
!    CALL write_cube( nnr, vioncc(:), 'dump.cube' )
!    !
!    stop
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE calc_vioncc
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE calc_rhostern( n, charge, dipole, x0, axis, step, v0, v, rho )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: n
!    !
!    REAL( DP ), INTENT(IN) :: charge, dipole
!    REAL( DP ), INTENT(IN) :: x0
!    REAL( DP ), INTENT(IN) :: axis(n)
!    REAL( DP ), INTENT(IN) :: step(n)
!    REAL( DP ), INTENT(IN) :: v0(n)
!    REAL( DP ), INTENT(OUT) :: v(n)
!    REAL( DP ), INTENT(OUT) :: rho(n)
!    !
!    INTEGER :: i, i1, i2
!    REAL( DP ) :: ez, fact, vstern, const
!    REAL( DP ) :: v1, v2, dv, vbound
!    REAL( DP ) :: arg, asinh, coth, acoth
!    REAL( DP ) :: f1, f2, f3, f4
!    !
!    ! ... Given the total explicit charge, the value of the field at the boundary
!    !     is obtained by Gauss's law
!    ! (1) ez = - tpi * e2 * charge * axis_length / omega / env_static_permittivity
!    !
!    ! ... By integrating the Guy-Chapman model one has a relation that links the derivative
!    !     of the potential (the field) to the value of the potential
!    ! (2) dv/dz = - ez = fact * SINH( v(z) * zion / 2 / kbt )
!    !     where
!    ! (3) fact = - e2 * SQRT( 32 * pi * cd * kbt / e2 / env_static_permittivity )
!    !
!    ! ... By combining (1) and (2) one can derive the analytic charge from the knowledge of the potential
!    !     at the boundary,
!    ! (4) charge_ext = fact * env_static_permittivity * omega / axis_lenght / tpi / e2 * SINH( vskin * zion / 2 / kbt )
!    !
!    ! ... or one can compute the value of the potential at the interface corresponding to a certain
!    !     explicit charge density,
!    ! (5) vskin_analytic = 2.D0 * kbt / zion * ASINH( ez / fact )
!    !
!    ! ... Eventually, by integrating Eq. (2) the analytic form of the potential is found
!    ! (6) vanalytic(z) = 4.D0 * kbt / zion * ACOTH( const * EXP( - z * fact * zion / kbt / 2.D0  ) )
!    !     were the constant is determined by the condition that vanalytic(xskin) = vskin
!    !
!    ! ... Compute the physical properties of the interface
!    !
!    ez = - tpi * e2 * charge / area / env_static_permittivity
!    fact = - e2 * SQRT( 8.D0 * fpi * cion * kbt / e2 / env_static_permittivity )
!    arg = ez/fact
!    asinh = LOG(arg + SQRT( arg**2 + 1 ))
!    vstern = 2.D0 * kbt / zion * asinh
!    arg = vstern * 0.25D0 * invkbt * zion
!    coth = ( EXP( 2.D0 * arg ) + 1.D0 ) / ( EXP( 2.D0 * arg ) - 1.D0 )
!    const = coth * EXP( zion * fact * invkbt * 0.5D0 * x0 )
!    !
!!    IF ( verbose .GE. 2 ) WRITE(environ_unit,8001)ez,fact,vstern,const
!    !
!    ! ... Compute the potential of the source charge at the boundary, to set the constant shift
!    !     Linear interpolate between gridpoints around x0 and average over both sides of slab
!    !
!    i1 = cell_min + naxis/2 + INT( x0 / dx )
!    i2 = cell_min + naxis/2 - 1 - INT( x0 / dx )
!    IF ( i1 .LE. 0 .OR. i1 .GT. n-1 .OR. i2 .LE. 1 .OR. i2 .GT. n ) THEN
!      WRITE(stdout,8004)
!      STOP
!    ENDIF
!    v1 = v0(i1)   + ( x0 - axis(i1) )   * ( v0(i1+1) - v0(i1) ) / ( axis(i1+1)-axis(i1) )
!    v2 = v0(i2-1) + ( - x0 - axis(i2-1) ) * ( v0(i2) - v0(i2-1) ) / ( axis(i2)-axis(i2-1) )
!    dv = 2.D0 * fpi * dipole / area
!    vbound = ( v1 + v2 ) * 0.5D0
!    !
!!    IF ( verbose .GE. 3 ) WRITE(environ_unit,8003)i1,axis(i1),axis(i1+1),v0(i1),v0(i1+1),v1
!!    IF ( verbose .GE. 3 ) WRITE(environ_unit,8003)i2,axis(i2-1),axis(i2),v0(i2-1),v0(i2),v2
!!    IF ( verbose .GE. 3 ) WRITE(environ_unit,8005)vbound,dv,v1-dv/2.D0,v2+dv/2.D0
!    !
!    ! ... Compute some constants needed for the calculation
!    !
!    f1 = - fact * zion * invkbt * 0.5D0
!    f2 = 4.D0 * kbt / zion
!    f3 = -2.D0 * zion * cion
!    f4 = zion * invkbt
!    !
!!    IF ( verbose .GE. 3 ) WRITE(environ_unit,8002)f1,f2,f3,f4
!    !
!    ! ... Compute the analytic potential and charge
!    !
!    v = 0.D0
!    rho = 0.D0
!    DO i = 1, n
!      IF ( ABS(axis(i)) .LT. x0 ) THEN
!        !
!        ! ... Constant shift of the potential in the inside
!        !
!        v(i) = vstern - vbound
!        !
!      ELSE
!        !
!        arg = const * EXP( ABS(axis(i)) * f1)
!        IF ( ABS(arg) .GT. 1.D0 ) THEN
!          acoth = 0.5D0 * LOG( (arg + 1.D0) / (arg - 1.D0) )
!        ELSE
!          acoth = 0.D0
!        END IF
!        v(i) =  f2 * acoth
!        rho(i) = f3 * SINH( f4 * v(i) )
!        !
!        ! ... Need to remove source potential from the outside
!        !
!        v(i) = v(i) + dv*step(i)  - v0(i)
!        !
!      ENDIF
!    ENDDO
!    !
!    RETURN
!    !
!8001 FORMAT(1X,'ez = ',E14.6,' fact =',E14.6,' v0 = ',E14.6,' const = ',E14.6)
!8002 FORMAT(1X,'f1 = ',E14.6,' f2 = ',E14.6,' f3 = ',E14.6,' f4 = ',E14.6)
!8003 FORMAT(1X,'i = ',i4,' x1 = ',F14.6,' x2 = ',F14.6,' v1 = ',F14.6,' v2 = ',F14.6,' v = ',F14.6)
!8004 FORMAT(1X,'ERROR: stern layer is outside of 1d cell, need to use more nrep')
!8005 FORMAT(1X,'vbound = ',F14.6,' dv = ',F14.6,' v1_shifted = ',F14.6,' v2_shifted = ',F14.6)
!!--------------------------------------------------------------------
!  END SUBROUTINE calc_rhostern
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE calc_cstern( n, nsp, charge, dipole, x0, axis, step, v0, v, c )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: n, nsp
!    !
!    REAL( DP ), INTENT(IN) :: charge, dipole
!    REAL( DP ), INTENT(IN) :: x0
!    REAL( DP ), INTENT(IN) :: axis(n)
!    REAL( DP ), INTENT(IN) :: step(n)
!    REAL( DP ), INTENT(IN) :: v0(n)
!    REAL( DP ), INTENT(OUT) :: v(n)
!    REAL( DP ), INTENT(OUT) :: c(nsp,n)
!    !
!    INTEGER :: i, i1, i2
!    REAL( DP ) :: ez, fact, vstern, const
!    REAL( DP ) :: v1, v2, dv, vbound
!    REAL( DP ) :: arg, asinh, coth, acoth
!    REAL( DP ) :: f1, f2, f3, f4
!    !
!    ! ... Given the total explicit charge, the value of the field at the boundary
!    !     is obtained by Gauss's law
!    ! (1) ez = - tpi * e2 * charge * axis_length / omega / env_static_permittivity
!    !
!    ! ... By integrating the Guy-Chapman model one has a relation that links the derivative
!    !     of the potential (the field) to the value of the potential
!    ! (2) dv/dz = - ez = fact * SINH( v(z) * zion / 2 / kbt )
!    !     where
!    ! (3) fact = - e2 * SQRT( 32 * pi * cd * kbt / e2 / env_static_permittivity )
!    !
!    ! ... By combining (1) and (2) one can derive the analytic charge from the knowledge of the potential
!    !     at the boundary,
!    ! (4) charge_ext = fact * env_static_permittivity * omega / axis_lenght / tpi / e2 * SINH( vskin * zion / 2 / kbt )
!    !
!    ! ... or one can compute the value of the potential at the interface corresponding to a certain
!    !     explicit charge density,
!    ! (5) vskin_analytic = 2.D0 * kbt / zion * ASINH( ez / fact )
!    !
!    ! ... Eventually, by integrating Eq. (2) the analytic form of the potential is found
!    ! (6) vanalytic(z) = 4.D0 * kbt / zion * ACOTH( const * EXP( - z * fact * zion / kbt / 2.D0  ) )
!    !     were the constant is determined by the condition that vanalytic(xskin) = vskin
!    !
!    ! ... Compute the physical properties of the interface
!    !
!    ez = - tpi * e2 * charge / area / env_static_permittivity
!    fact = - e2 * SQRT( 8.D0 * fpi * cion * kbt / e2 / env_static_permittivity )
!    arg = ez/fact
!    asinh = LOG(arg + SQRT( arg**2 + 1 ))
!    vstern = 2.D0 * kbt / zion * asinh
!    arg = vstern * 0.25D0 * invkbt * zion
!    coth = ( EXP( 2.D0 * arg ) + 1.D0 ) / ( EXP( 2.D0 * arg ) - 1.D0 )
!    const = coth * EXP( zion * fact * invkbt * 0.5D0 * x0 )
!    !
!!    IF ( verbose .GE. 2 ) WRITE(environ_unit,8001)ez,fact,vstern,const
!    !
!    ! ... Compute the potential of the source charge at the boundary, to set the constant shift
!    !     Linear interpolate between gridpoints around x0 and average over both sides of slab
!    !
!    i1 = cell_min + naxis/2 + INT( x0 / dx )
!    i2 = cell_min + naxis/2 - 1 - INT( x0 / dx )
!    IF ( i1 .LE. 0 .OR. i1 .GT. n-1 .OR. i2 .LE. 1 .OR. i2 .GT. n ) THEN
!      WRITE(stdout,8004)
!      STOP
!    ENDIF
!    v1 = v0(i1)   + ( x0 - axis(i1) )   * ( v0(i1+1) - v0(i1) ) / ( axis(i1+1)-axis(i1) )
!    v2 = v0(i2-1) + ( - x0 - axis(i2-1) ) * ( v0(i2) - v0(i2-1) ) / ( axis(i2)-axis(i2-1) )
!    dv = 2.D0 * fpi * dipole / area
!    vbound = ( v1 + v2 ) * 0.5D0
!    !
!!    IF ( verbose .GE. 3 ) WRITE(environ_unit,8003)i1,axis(i1),axis(i1+1),v0(i1),v0(i1+1),v1
!!    IF ( verbose .GE. 3 ) WRITE(environ_unit,8003)i2,axis(i2-1),axis(i2),v0(i2-1),v0(i2),v2
!!    IF ( verbose .GE. 3 ) WRITE(environ_unit,8005)vbound,dv,v1-dv/2.D0,v2+dv/2.D0
!    !
!    ! ... Compute some constants needed for the calculation
!    !
!    f1 = - fact * zion * invkbt * 0.5D0
!    f2 = 4.D0 * kbt / zion
!    f3 = cion / cionmax
!    f4 = zion * invkbt
!    !
!!    IF ( verbose .GE. 3 ) WRITE(environ_unit,8002)f1,f2,f3,f4
!    !
!    ! ... Compute the analytic potential and charge
!    !
!    v = 0.D0
!    c = 0.D0
!    DO i = 1, n
!      IF ( ABS(axis(i)) .LT. x0 ) THEN
!        !
!        ! ... Constant shift of the potential in the inside
!        !
!        v(i) = vstern - vbound
!        !
!      ELSE
!        !
!        arg = const * EXP( ABS(axis(i)) * f1)
!        IF ( ABS(arg) .GT. 1.D0 ) THEN
!          acoth = 0.5D0 * LOG( (arg + 1.D0) / (arg - 1.D0) )
!        ELSE
!          acoth = 0.D0
!        END IF
!        v(i) =  f2 * acoth
!        c(1,i) = f3 * ( EXP( - f4 * v(i) ) ) ! - 1.D0 )
!        c(2,i) = f3 * ( EXP( f4 * v(i) ) ) ! - 1.D0 )
!        !
!        ! ... Need to remove source potential from the outside
!        !
!        v(i) = v(i) + dv*step(i)  - v0(i)
!        !
!      ENDIF
!    ENDDO
!    !
!    RETURN
!    !
!8001 FORMAT(1X,'ez = ',E14.6,' fact =',E14.6,' v0 = ',E14.6,' const = ',E14.6)
!8002 FORMAT(1X,'f1 = ',E14.6,' f2 = ',E14.6,' f3 = ',E14.6,' f4 = ',E14.6)
!8003 FORMAT(1X,'i = ',i4,' x1 = ',F14.6,' x2 = ',F14.6,' v1 = ',F14.6,' v2 = ',F14.6,' v = ',F14.6)
!8004 FORMAT(1X,'ERROR: stern layer is outside of 1d cell, need to use more nrep')
!8005 FORMAT(1X,'vbound = ',F14.6,' dv = ',F14.6,' v1_shifted = ',F14.6,' v2_shifted = ',F14.6)
!!--------------------------------------------------------------------
!  END SUBROUTINE calc_cstern
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE setfix( nnr, nspin, n1d, rho )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT( IN ) :: nnr, nspin, n1d
!    REAL( DP ), DIMENSION( nnr ), INTENT( IN ) :: rho
!    !
!    ! ... Local variables
!    !
!    REAL( DP ) :: charge, ehart
!    REAL( DP ), DIMENSION( nnr ) :: vperiodic
!    REAL( DP ), DIMENSION( nnr, nspin ) :: vtot
!    !
!    ! ... Store full source charge
!    !
!    IF ( ALLOCATED(rhofix) ) DEALLOCATE(rhofix)
!    ALLOCATE( rhofix( nnr ) )
!    rhofix = rho
!    !
!    ! ... Compute the non-periodic potential of the source charge
!    !
!    IF ( ALLOCATED( vfix ) ) DEALLOCATE( vfix )
!    ALLOCATE( vfix( nnr ) )
!    vtot = 0.D0
!    CALL v_h_of_rho_r( rhofix, ehart, charge, vtot )
!    vperiodic = 0.D0
!    CALL calc_vperiodic( nnr, nspin, .FALSE., rhofix, vperiodic )
!    vfix( : ) = vtot( :, 1 ) + vperiodic(:)
!    CALL write_cube( nnr, vfix, 'vfix.cube' )
!    !
!    ! ... Convert the total source charge into 1D ...
!    !
!    IF ( ALLOCATED(rhofix1d) ) DEALLOCATE( rhofix1d )
!    ALLOCATE( rhofix1d(n1d) )
!    rhofix1d = 0.D0
!    CALL planar_average( nnr, naxis, slab_axis, shift, .false., rhofix, rhofix1d(cell_min:cell_max) )
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, rhofix1d, 'rhofix1d.dat' )
!    !
!    ! ... and compute its non-periodic potential
!    !
!    IF ( ALLOCATED( vfix1d ) ) DEALLOCATE( vfix1d )
!    ALLOCATE( vfix1d(n1d) )
!    vfix1d = 0.D0
!    CALL v1d_of_rho1d( n1d, dx, rhofix1d, vfix1d )
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, vfix1d, 'vfix1d.dat' )
!    !
!    ! ... Compute and store average total potentials
!    !
!    vfix_avg = SUM(vfix)/DBLE(ntot)
!    CALL mp_sum( vfix_avg, intra_bgrp_comm )
!    vfix1d_avg = SUM(vfix1d(cell_min:cell_max))/DBLE(naxis)
!    IF ( verbose .GE. 1 ) WRITE(environ_unit,*)'vfix_avg = ',vfix_avg,' vfix1d_avg = ',vfix1d_avg
!    !
!    ! ... Compute and store total charge and dipole and internal potential drop
!    !
!    charge_fix = SUM(rhofix1d(:))*darea
!    dipole_fix = SUM(rhofix1d(:)*axis1d(:))*darea
!    dv_fix = 2.D0 * fpi * dipole_fix / area
!    IF ( verbose .GE. 1 ) WRITE( environ_unit, * )'charge_fix = ',charge_fix,' dipole_fix = ',dipole_fix
!    IF ( verbose .GE. 1 ) WRITE( environ_unit, * )'dv_fix = ',dv_fix
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE setfix
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE cleanfix( )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    ! ... Clean local quantities
!    !
!    IF ( ALLOCATED( rhofix   ) ) DEALLOCATE( rhofix   )
!    IF ( ALLOCATED( rhofix1d ) ) DEALLOCATE( rhofix1d )
!    IF ( ALLOCATED( vfix     ) ) DEALLOCATE( vfix     )
!    IF ( ALLOCATED( vfix1d   ) ) DEALLOCATE( vfix1d   )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE cleanfix
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE vioncc1d_fgeval(n,v,f,g)
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: n
!    REAL( DP ), DIMENSION(n), INTENT(IN) :: v
!    REAL( DP ), INTENT(OUT) :: f
!    REAL( DP ), DIMENSION(n), INTENT(OUT) :: g
!    !
!    INTEGER :: i, j
!    REAL( DP ) :: df, d2f
!    REAL( DP ), DIMENSION(n) :: residue, dresidue, rtmp
!    !
!    ! ...The error is computed from the norm of the residue
!    !    which is taken to be the difference of d/dx(eps(x)*d/dx\phi(x))
!    !    from -4*pi*e2*(rhofix+rhoioncc)
!!    CALL write_1d( n, axis1d, v, 'v_input.dat' )
!    !
!    CALL generate_rhoioncc( n, 1, switch1d, v, residue )
!    !
!    residue = fpi * e2 * ( residue + rhofix1d )
!!    CALL write_1d( n, axis1d, residue, 'residue_rho.dat' )
!    !
!    ! ...The gradient of the error is 2 * residue \cdot (d/d\phi residue)
!    !
!    CALL generate_drhoioncc( n, 1, switch1d, v, dresidue )
!    !
!    dresidue = fpi * e2 * dresidue
!!    CALL write_1d( n, axis1d, dresidue, 'dresidue.dat' )
!    !
!    ! ...Compute the first part of the residue vector and the total error
!    !
!    rtmp = 0.D0
!    DO i = nfdpoint+1, n-nfdpoint
!       !
!!       df = 0.D0
!       d2f = 0.D0
!       DO j = -nfdpoint, nfdpoint
!!          df = df + icfd(j) * v(i+j)
!          d2f = d2f + icfd2(j) * v(i+j)
!       ENDDO
!!       df = df / DBLE(ncfd) / dx
!       d2f = d2f / DBLE(ncfd2) / dx**2
!       !
!       rtmp(i) = d2f
!       residue(i) = residue(i) + d2f !*eps1d(i) + deps1d(i) * df
!       !
!    ENDDO
!!    CALL write_1d( n, axis1d, rtmp, 'residue_deriv.dat' )
!!    CALL write_1d( n, axis1d, residue, 'residue_full.dat' )
!    !
!    f = SUM(residue**2)
!    !
!    ! ... Compute the gradient (the derivative of the first part is constant)
!    !
!    g = 0.D0
!    DO i = nfdpoint+1, n-nfdpoint
!       DO j = -nfdpoint, nfdpoint
!!          df = icfd(j) / DBLE(ncfd) / dx
!          d2f = icfd2(j) / DBLE(ncfd2) / dx**2
!          g(i) = g(i) + residue(i+j) * d2f ! (d2f * eps(i+j) + df * deps(i+j))
!       ENDDO
!       g(i) = g(i) + residue(i) * dresidue(i)
!    ENDDO
!    g = g * 2.D0
!!    CALL write_1d( n, axis1d, g, 'gradient.dat' )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE vioncc1d_fgeval
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE vioncc1d_int_fgeval(n,v,f,g)
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: n
!    REAL( DP ), DIMENSION(n), INTENT(INOUT) :: v
!    REAL( DP ), INTENT(OUT) :: f
!    REAL( DP ), DIMENSION(n), INTENT(OUT) :: g
!    !
!    INTEGER :: i, j, isp
!    REAL( DP ) :: df, d2f
!    REAL( DP ), DIMENSION(n) :: rhotmp, residue, dresidue, l, s, vtmp
!    REAL( DP ), DIMENSION(nsp,n) :: c
!    !
!    ! ...The error is computed from the variational expression of the
!    !    free energy, written in terms of the potential
!    !    F = \int ( -0.5*|dv/dx|**2 + 4*pi*e2*(rhofix+rhoioncc)*v ) dx
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, v, 'v_input.dat' )
!    !
!    ! ... Compute the first derivative of the potential
!    !
!    g = 0.D0
!    DO i = nfdpoint+1, n-nfdpoint
!       DO j = -nfdpoint, nfdpoint
!          g(i) = g(i) + icfd(j) * v(i+j)
!       ENDDO
!       g(i) = g(i) / DBLE(ncfd) / dx
!    ENDDO
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, g, 'gradv.dat' )
!    !
!    vtmp = v ! + vzero1d
!    CALL generate_cioncc( n, nsp, 1, switch1d, vtmp, c )
!    !
!    CALL sum_cioncc( n, nsp, switch1d, z, c, rhotmp )
!    !rhotmp = rhotmp - rhozero1d
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, rhotmp, 'rhotest.dat' )
!    !
!    CALL calc_sioncc( n, nsp, switch1d, c, s )
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, s, 'entropy.dat' )
!    !
!    residue = 0.D0
!    CALL sum_cioncc( n, nsp, switch1d, mu, c, residue )
!    !
!    residue(:) = residue(:) + 0.5D0 / fpi / e2 * g(:)**2 - ( rhofix1d(:) + rhotmp(:) ) * v(:) + s(:)
!!    residue(:) = residue(:) + 0.5D0 / fpi / e2 * g(:)**2 - rhotmp(:) * vtmp(:) + s(:)
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, residue, 'residue.dat' )
!    !
!    f = SUM(residue)*dx
!    !
!    ! ...The gradient of the free energy: gradient of gradient (laplacian)
!    !    plus density
!    !
!    l = 0.D0
!    DO i = nfdpoint+1, n-nfdpoint
!       DO j = -nfdpoint, nfdpoint
!          l(i) = l(i) + icfd(j) * g(i+j)
!       ENDDO
!       l(i) = l(i) / DBLE(ncfd) / dx
!    ENDDO
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, l, 'laplv.dat' )
!    !
!    g(:) = l(:) / fpi / e2 + ( rhotmp(:) + rhofix1d(:) )
!!    g(:) = l(:) / fpi / e2 + rhotmp(:)
!    g = - g * dx
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, g, 'gradf.dat' )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE vioncc1d_int_fgeval
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE rhoioncc1d_fgeval(n,nsp,x,f,g)
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: n, nsp
!    REAL( DP ), DIMENSION(nsp,n), INTENT(IN) :: x
!    REAL( DP ), INTENT(OUT) :: f
!    REAL( DP ), DIMENSION(nsp,n), INTENT(OUT) :: g
!    !
!    INTEGER :: isp
!    REAL( DP ), DIMENSION(nsp, n) :: residue, c, ctmp
!    REAL( DP ), DIMENSION(n) :: rhotmp, vtmp
!    !
!    ! ... The error is computed from the variational expression of the
!    !    free energy, written in terms of the countercharge concentration
!    !    F = \int ( -0.5*|dv/dx|**2 + 4*pi*e2*(rhofix+rhoioncc)*v ) dx
!    !
!    DO isp = 1, nsp
!       c(isp,:) = x(isp,:)**2
!!       c(isp,:) = ( ATAN(x(isp,:)) / pi + 0.5D0 )
!    ENDDO
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, c(1,:), 'c1_input.dat' )
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, c(2,:), 'c2_input.dat' )
!    !
!    ! ... Compute the charge density from the concentrations
!    !
!    CALL sum_cioncc( n, nsp, switch1d, z, c, rhotmp )
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, rhotmp, 'rho_input.dat' )
!    !
!    ! ... Compute the potential
!    !
!    rhotmp = rhotmp + rhofix1d
!    !
!    CALL v1d_of_rho1d( n, dx, rhotmp, vtmp )
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, vtmp, 'v_input.dat' )
!    !
!    CALL generate_cioncc( n, nsp, 1, switch1d, vtmp, ctmp )
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, ctmp(1,:), 'ctmp1.dat' )
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, ctmp(2,:), 'ctmp2.dat' )
!    !
!    residue = c - ctmp
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, residue(1,:), 'residue1.dat' )
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, residue(2,:), 'residue2.dat' )
!    !
!    f = 0.D0
!    DO isp = 1, nsp
!       f = f + SUM( residue(isp, :)**2 )
!    ENDDO
!    !
!    ! ... Compute the gradient
!    !
!    g = 0.D0
!    CALL generate_dcioncc( n, nsp, 1, switch1d, vtmp, g )
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, g(1,:), 'dc1.dat' )
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, g(2,:), 'dc2.dat' )
!    !
!    rhotmp = 0.D0
!    !
!    DO isp = 1, nsp
!       rhotmp = rhotmp + g(isp,:) * residue(isp,:)
!    ENDDO
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, rhotmp, 'rhotmp.dat' )
!    !
!    CALL v1d_of_rho1d( n, dx, rhotmp, vtmp )
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, vtmp, 'vtmp.dat' )
!    !
!    DO isp = 1, nsp
!       !
!       g(isp,:) = residue(isp, :) - vtmp(:) * switch1d(:) * cionmax * z(isp)
!       g(isp,:) = g(isp,:) * 2.D0 * x(isp,:)
!!       g(isp,:) = g(isp,:) * 1.D0 / ( 1.D0 + x(isp,:)**2 ) / pi
!       !
!    ENDDO
!    !
!    g = g * 2.D0
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, g(1,:), 'gradf.dat' )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE rhoioncc1d_fgeval
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE rhoioncc1d_int_fgeval(n,nsp,x,f,g)
!!--------------------------------------------------------------------
!    !
!    USE numerical_recipes, ONLY: zbrent, zbrac
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: n, nsp
!    REAL( DP ), DIMENSION(nsp,n), INTENT(IN) :: x
!    REAL( DP ), INTENT(OUT) :: f
!    REAL( DP ), DIMENSION(nsp,n), INTENT(OUT) :: g
!    !
!    INTEGER :: i, j, isp
!    REAL( DP ) :: df, d2f, slope, vavg
!    REAL( DP ), DIMENSION(nsp,n) :: c
!    REAL( DP ), DIMENSION(n) :: rhotmp, residue, dresidue, l, s, vtmp, gtmp, verr
!    !
!    REAL( DP ) :: dlambda, lambda, x1, x2, f1, f2, flambda, tol
!    LOGICAL :: success
!    !
!    ! ... The error is computed from the variational expression of the
!    !    free energy, written in terms of the countercharge concentration
!    !    F = \int ( -0.5*|dv/dx|**2 + 4*pi*e2*(rhofix+rhoioncc)*v ) dx
!    !
!    DO isp = 1, nsp
!       c(isp,:) = x(isp,:)**2
!    ENDDO
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, c(1,:), 'c1_input.dat' )
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, c(2,:), 'c2_input.dat' )
!    !
!    ! ... Compute the charge density from the concentrations
!    !
!    rhotmp = 0.D0
!    CALL sum_cioncc( n, nsp, switch1d, z, c, rhotmp )
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, rhotmp, 'rho_input.dat' )
!    !
!    ! ... Compute the potential
!    !
!    rhotmp = rhotmp !+ rhofix1d
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, rhotmp, 'rhotmp.dat' )
!    !
!!    slope = SUM(rhotmp) * dx * fpi * e2 / 2.D0
!    IF ( verbose .GE. 2 ) WRITE(environ_unit,*)'charge = ',SUM(rhotmp)*darea
!    !
!    vtmp = 0.D0
!    CALL v1d_of_rho1d( n, dx, rhotmp, vtmp )
!    !
!    vtmp = vtmp + vfix1d
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, vtmp, 'vtmp.dat' )
!    !
!    ! ... Compute the gradient of the potential
!    !
!!    gtmp = - slope * step1d(:) * 2.D0
!    gtmp = 0.D0
!    DO i = nfdpoint + 1, n - nfdpoint ! 2*nfdpoint+1, n-2*nfdpoint
!!       gtmp(i) = 0.D0
!       DO j = -nfdpoint, nfdpoint
!          gtmp(i) = gtmp(i) + icfd(j) * vtmp(i+j)
!       ENDDO
!       gtmp(i) = gtmp(i) / DBLE(ncfd) / dx
!    ENDDO
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, gtmp, 'gtmp.dat' )
!    !
!    ! ... Compute the entropy
!    !
!    CALL calc_sioncc( n, nsp, switch1d, c, s )
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, s, 'entropy.dat' )
!    !
!    ! ... Constraint the countercharge, compute the multiplier
!    !
!    vzero1d = vtmp
!    dlambda = 1.D-2
!    x1 = dlambda
!    x2 = -dlambda
!    CALL zbrac(ioncc_charge, x1, x2, success)
!!    f1 = ioncc_charge(x1)
!!    f2 = ioncc_charge(x2)
!!    IF ( verbose .EQ. 4 ) WRITE(environ_unit,*)'x1 = ',x1,' f1 = ',f1
!!    IF ( verbose .EQ. 4 ) WRITE(environ_unit,*)'x2 = ',x2,' f2 = ',f2
!    IF ( .NOT. success ) THEN
!       WRITE(environ_unit,*)'ERROR: could not find a zero'
!       STOP
!    ENDIF
!    tol=1.D-8
!    CALL zbrent(ioncc_charge, x1, x2, tol, lambda)
!    flambda = ioncc_charge(lambda)
!    IF ( verbose .EQ. 2 ) WRITE(environ_unit,*)'lambda = ',lambda,' flambda = ',flambda
!    !
!    ! ... Compute the chemical potential terms
!    !
!    residue = 0.D0
!    CALL sum_cioncc( n, nsp, switch1d, mu, c, residue )
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, residue, 'sumctimesmu.dat')
!    !
!    residue(:) = residue(:) + 0.5D0 / fpi / e2 * gtmp(:)**2 - ( rhotmp(:) + rhofix1d(:) ) * vtmp(:) + s(:) &
!         & + lambda * ( rhotmp(:) + rhofix1d(:) )
!    !
!    f =  SUM(residue)*dx
!    !
!!!!TMP!!!    l = 0.D0
!!!!TMP!!!    DO i = nfdpoint+1, n-nfdpoint
!!!!TMP!!!       DO j = -nfdpoint, nfdpoint
!!!!TMP!!!          l(i) = l(i) + icfd(j) * gtmp(i+j)
!!!!TMP!!!       ENDDO
!!!!TMP!!!       l(i) = l(i) / DBLE(ncfd) / dx
!!!!TMP!!!    ENDDO
!!!!TMP!!!    gtmp(:) = l(:) / fpi / e2 + rhotmp(:)
!!!!TMP!!!    !
!!!!TMP!!!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, gtmp, 'delta.dat' )
!!!!TMP!!!    !
!!!!TMP!!!    CALL v1d_of_rho1d( n, dx, gtmp, verr )
!!!!TMP!!!    !
!!!!TMP!!!    verr(:) = verr(:) * z(1) * cionmax * switch1d(:)
!!!!TMP!!!    !
!!!!TMP!!!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, verr, 'verr.dat' )
!    !
!    ! ...The gradient of the free energy: ds/dc + mu + z * v
!    !
!    g = 0.D0
!    CALL calc_dsioncc( n, nsp, switch1d, c, g )
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, g(1,:), 'ds1.dat' )
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, g(2,:), 'ds2.dat' )
!    !
!    DO isp = 1, nsp
!       g(isp,:) = g(isp,:) + ( mu(isp) - z(isp) * ( vtmp(:) - lambda ) ) * cionmax * switch1d(:)
!       g(isp,:) = g(isp,:) * 2.D0 * x(isp,:)
!    ENDDO
!    !
!    g = g * dx
!    !
!    gtmp(:) = mu(1) * cionmax * switch1d(:)
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, gtmp(:), 'mu.dat' )
!    !
!    gtmp(:) = - z(1) * vtmp(:) * cionmax * switch1d(:)
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, gtmp(:), 'zvtmp.dat' )
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, g(1,:), 'gradf1.dat' )
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, g(2,:), 'gradf2.dat' )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE rhoioncc1d_int_fgeval
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE fullioncc1d_int_fgeval(n,nsp,x,f,g)
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: n, nsp
!    REAL( DP ), DIMENSION(nsp+1,n), INTENT(IN) :: x
!    REAL( DP ), INTENT(OUT) :: f
!    REAL( DP ), DIMENSION(nsp+1,n), INTENT(OUT) :: g
!    !
!    INTEGER :: i, j, isp
!    REAL( DP ) :: df, d2f, slope, vavg
!    REAL( DP ), DIMENSION(nsp,n) :: c
!    REAL( DP ), DIMENSION(n) :: rhotmp, residue, dresidue, l, s, vtmp, gtmp, verr
!    !
!    ! ... The error is computed from the variational expression of the
!    !    free energy, written in terms of the countercharge concentration
!    !    F = \int ( -0.5*|dv/dx|**2 + 4*pi*e2*(rhofix+rhoioncc)*v ) dx
!    !
!    DO isp = 1, nsp
!       c(isp,:) = x(isp,:)**2 ! + cb(isp)/cionmax
!    ENDDO
!    !
!    vtmp(:) = x(nsp+1,:)
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, c(1,:), 'c1_input.dat' )
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, c(2,:), 'c2_input.dat' )
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, vtmp(:), 'v_input.dat' )
!    !
!    ! ... Compute the charge density from the concentrations
!    !
!    CALL sum_cioncc( n, nsp, switch1d, z, c, rhotmp )
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, rhotmp, 'rho_input.dat' )
!    !
!    ! ... Compute the total charge
!    !
!    rhotmp = rhofix1d + rhotmp
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, rhotmp, 'rhotmp.dat' )
!    !
!    ! ... Compute the gradient of the potential
!    !
!    gtmp = 0.D0
!    DO i = nfdpoint+1, n-nfdpoint
!       DO j = -nfdpoint, nfdpoint
!          gtmp(i) = gtmp(i) + icfd(j) * vtmp(i+j)
!       ENDDO
!       gtmp(i) = gtmp(i) / DBLE(ncfd) / dx
!    ENDDO
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, gtmp, 'gtmp.dat' )
!    !
!    ! ... Compute the entropy
!    !
!    CALL calc_sioncc( n, nsp, switch1d, c, s )
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, s, 'entropy.dat' )
!    !
!    ! ... Compute the chemical potential terms
!    !
!    residue = 0.D0
!    CALL sum_cioncc( n, nsp, switch1d, mu, c, residue )
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, residue, 'sumctimesmu.dat')
!    !
!    residue(:) = residue(:) + 0.5D0 / fpi / e2 * gtmp(:)**2 - rhotmp(:) * vtmp(:) + s(:)
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, residue, 'residue.dat' )
!    !
!    f = SUM(residue)*dx
!    !
!    ! ...The gradient of the free energy wrt the potential: lapl(v)/fpi + rhotmp
!    !
!    g = 0.D0
!    !
!    l = 0.D0
!    DO i = nfdpoint+1, n-nfdpoint
!       DO j = -nfdpoint, nfdpoint
!          l(i) = l(i) + icfd(j) * gtmp(i+j)
!       ENDDO
!       l(i) = l(i) / DBLE(ncfd) / dx
!    ENDDO
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, l(:), 'lapl.dat' )
!    !
!    g(nsp+1,:) = - l(:) / fpi / e2 - rhotmp(:)
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, g(nsp+1,:), 'gv.dat' )
!    !
!    ! ...The gradient of the free energy wrt the ioncc concentration: ds/dc + mu + z * v
!    !
!    CALL calc_dsioncc( n, nsp, switch1d, c, g(1:nsp,:) )
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, g(1,:), 'ds1.dat' )
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, g(2,:), 'ds2.dat' )
!    !
!    DO isp = 1, nsp
!       g(isp,:) = g(isp,:) + ( mu(isp) - z(isp) * vtmp(:) ) * cionmax * switch1d(:)
!       g(isp,:) = g(isp,:) * 2.D0 * x(isp,:)
!    ENDDO
!    !
!    g = g * dx
!    !
!    gtmp(:) = mu(1) * cionmax * switch1d(:)
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, gtmp(:), 'mu.dat' )
!    !
!    gtmp(:) = - z(1) * vtmp(:) * cionmax * switch1d(:)
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, gtmp(:), 'zvtmp.dat' )
!    !
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, g(1,:), 'gradf1.dat' )
!    IF ( verbose .EQ. 4 ) CALL write_1d( n, axis1d, g(2,:), 'gradf2.dat' )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE fullioncc1d_int_fgeval
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE calc_fullioncc1d_lbfgs(  n1d, nsp, v, c )
!!--------------------------------------------------------------------
!    IMPLICIT NONE
!    !
!    REAL( DP ), PARAMETER :: tol = 1.D-8, xtol = 1.D-14
!    !
!    INTEGER, INTENT(IN) :: n1d, nsp
!    !
!    REAL( DP ), INTENT(INOUT) :: v(n1d)
!    REAL( DP ), INTENT(INOUT) :: c(nsp, n1d)
!    !
!    ! ... Local variables
!    !
!    LOGICAL :: diagc0
!    INTEGER :: n, m, length, iflag, iprint(2), nfdpoint, i, j, isp
!    REAL( DP ) :: f, gfd, delta, xtmp, ftmp
!    REAL( DP ), DIMENSION(n1d) :: rho
!    REAL( DP ), DIMENSION( : , : ), ALLOCATABLE :: g, gtmp
!    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: diag, work
!    REAL( DP ), DIMENSION( :, : ), ALLOCATABLE :: x
!    !
!    ! ... Initialize L-BFGS
!    !
!    diagc0 = .FALSE.
!    m = 20
!    n = ( nsp + 1 ) * n1d
!    length = n*(2*m+1)+2*m
!    iprint(1) = -1
!    iprint(2) = 0
!    IF ( verbose .EQ. 1 ) THEN
!      iprint(1) = 0
!    ENDIF
!    ALLOCATE( x( nsp+1, n1d ) )
!    ALLOCATE( g( nsp+1, n1d ) )
!    ALLOCATE( diag( n ) )
!    ALLOCATE( work( length ) )
!    !
!    ! ... Initial function and gradient evaluations
!    !
!    DO isp = 1, nsp
!       x(isp,:) = SQRT(ABS(c(isp,:)))
!    ENDDO
!    x(nsp+1,:) = v(:)
!    CALL write_1d( n1d, axis1d, x(1,:), 'c1_start.dat' )
!    CALL write_1d( n1d, axis1d, x(2,:), 'c2_start.dat' )
!    CALL write_1d( n1d, axis1d, x(3,:), 'v_start.dat' )
!    verbose = 4
!    CALL fullioncc1d_int_fgeval(n1d,nsp,x,f,g)
!    verbose = 2
!    !
!    ! ... Check gradient via finite differences
!    !
!    delta = 0.000001D0
!    ALLOCATE( gtmp( nsp+1, n1d ) )
!    DO isp = 1, nsp+1
!       WRITE(environ_unit,*)
!       DO i = 1, n1d
!          xtmp = x(isp,i)
!          x(isp,i) = xtmp + delta
!          CALL fullioncc1d_int_fgeval(n1d,nsp,x,f,gtmp)
!          gfd = f
!          x(isp,i) = xtmp - delta
!          CALL fullioncc1d_int_fgeval(n1d,nsp,x,f,gtmp)
!          gfd = ( gfd - f ) / 2.d0 / delta
!          x(isp,i) = xtmp
!          WRITE(environ_unit,'(4f20.10)')(DBLE(i)-naxis*(DBLE(nrep)+0.5D0))*dx,g(isp,i)-gfd,g(isp,i),gfd
!       ENDDO
!    ENDDO
!    CALL flush(environ_unit)
!    DEALLOCATE( gtmp )
!    !
!    i = 0
!    iflag = 0
!    DO
!      CALL lbfgs( n, m, x, f, g, diagc0, diag, iprint, tol, xtol, work, iflag)
!      IF ( iflag .LT. 0 ) THEN ! linesearch failed for some reason
!         WRITE(stdout,*)&
!           &'ERROR: lbfgs failed to converge in rhoioncc1d with iflag = ',iflag
!         EXIT
!         STOP
!      ELSE IF ( iflag .EQ. 0 ) THEN ! Calculation converged, exit the loop
!         EXIT
!      ELSE IF ( iflag .EQ. 1 ) THEN ! lbfgs is asking to evaluate f and g
!         i = i + 1
!         IF ( i .EQ. 1 ) EXIT
!         CALL fullioncc1d_int_fgeval(n1d,nsp,x,f,g)
!         write(environ_unit,*)'fullioncc_fgeval = ',i,f
!         CALL flush(environ_unit)
!      ELSE IF ( iflag .GT. 1 ) THEN ! not supposed to happen with diagc0=.FALSE.
!         WRITE(stdout,*)&
!           &'ERROR: unexpected flag from lbfgs optimization, iflag = ',iflag
!         STOP
!      ENDIF
!    ENDDO
!    !
!    verbose = 4
!    CALL fullioncc1d_int_fgeval(n1d,nsp,x,f,g)
!    verbose = 2
!    !
!    ! ... Check gradient via finite differences
!    !
!    delta = 0.000001D0
!    ALLOCATE( gtmp( nsp+1, n1d ) )
!    DO isp = 1, nsp + 1
!       WRITE(environ_unit,*)
!       DO j = 1, n1d
!          xtmp = x(isp,j)
!          x(isp,j) = xtmp + delta
!          CALL fullioncc1d_int_fgeval(n1d,nsp,x,f,gtmp)
!          gfd = f
!          x(isp,j) = xtmp - delta
!          CALL fullioncc1d_int_fgeval(n1d,nsp,x,f,gtmp)
!          gfd = ( gfd - f ) / 2.d0 / delta
!          x(isp,j) = xtmp
!          WRITE(environ_unit,'(4f20.10)')(DBLE(j)-naxis*(DBLE(nrep)+0.5D0))*dx,g(isp,j)-gfd,g(isp,j),gfd
!       ENDDO
!    ENDDO
!    CALL flush(environ_unit)
!    DEALLOCATE( gtmp )
!    !
!    c(1:nsp,:) = x(1:nsp,:)
!    v(:) = x(nsp+1,:)
!    !
!    rho = 0.D0
!    CALL sum_cioncc( n1d, nsp, switch1d, z, c, rho )
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, rho, 'rhoioncc1d.dat' )
!    !
!    CALL v1d_of_rho1d( n1d, dx, rho, v )
!    v = v + vfix1d
!    !
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, v, 'v1d.dat' )
!    v = v - vfix1d - vfix_avg + vfix1d_avg
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, v, 'vext1d.dat' )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE calc_fullioncc1d_lbfgs
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE calc_rhoioncc1d_lbfgs(  n1d, nsp, v, c )
!!--------------------------------------------------------------------
!    IMPLICIT NONE
!    !
!    REAL( DP ), PARAMETER :: tol = 1.D-6, xtol = 1.D-14
!    !
!    INTEGER, INTENT(IN) :: n1d, nsp
!    !
!    REAL( DP ), INTENT(OUT) :: v(n1d)
!    REAL( DP ), INTENT(INOUT) :: c(nsp, n1d)
!    !
!    ! ... Local variables
!    !
!    LOGICAL :: diagc0
!    INTEGER :: n, m, length, iflag, iprint(2), nfdpoint, i, j, isp
!    REAL( DP ) :: charge_tot, f, gfd, delta, ctmp, ftmp, fv, fvtmp, fx
!    REAL( DP ) :: charge(nsp), fact(nsp), tot, integral
!    REAL( DP ), DIMENSION(n1d) :: rho
!    REAL( DP ), DIMENSION( : , : ), ALLOCATABLE :: g, gtmp, x, gx, c2
!    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: diag, work, gv
!    !
!    ! ... Build the constant approximation of rhoioncc, rescale the guess
!    !
!    WRITE(environ_unit,*)'cionmax ',cionmax,'cion ',cion,'cion/cionmax ',cion/cionmax
!    tot = 0.D0
!    fact = 0.D0
!    DO isp = 1, nsp
!       IF ( z(isp)*charge_fix .LT. 0.D0 ) THEN
!          tot = tot + 1.D0
!          fact(isp) = - charge_fix / z(isp) / cionmax
!       ENDIF
!    ENDDO
!    IF ( tot .EQ. 0.D0 ) THEN
!       WRITE(environ_unit,*)'ERROR: missing compensating ionic charge in ioncc'
!       STOP
!    ENDIF
!    WRITE(environ_unit,*)' tot = ',tot
!    integral = SUM(switch1d(:)**2*(1.D0-switch1d(:))) * darea
!    WRITE(environ_unit,*)'integral = ',integral
!    fact = fact / tot / integral
!    WRITE(environ_unit,*)'fact =', fact
!    !
!    IF ( ALLOCATED(czero1d) ) DEALLOCATE(czero1d)
!    ALLOCATE(czero1d(nsp,n1d))
!    charge = 0.D0
!    DO isp = 1, nsp
!       czero1d(isp,:) = ( cb(isp) + fact(isp) * ( 1.D0 - switch1d(:) ) ) * switch1d(:)
!       charge(isp) = SUM(czero1d(isp,:)*switch1d(:)) * darea * cionmax * z(isp)
!       WRITE(environ_unit,*)'charge(',isp,') = ',charge(isp)
!    ENDDO
!    charge_tot = SUM(charge)
!    WRITE(environ_unit,*)'charge = ',charge_tot,' charge_fix = ',charge_fix
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, czero1d(1,:), 'c1_zero.dat' )
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, czero1d(2,:), 'c2_zero.dat' )
!    !
!    IF ( ALLOCATED(rhozero1d) ) DEALLOCATE(rhozero1d)
!    ALLOCATE(rhozero1d(n1d))
!    rhozero1d = 0.D0
!    CALL sum_cioncc( n1d, nsp, switch1d, z, czero1d, rhozero1d )
!    IF ( verbose .GE. 3 ) CALL write_1d( n1d, axis1d, rhozero1d, 'rhozero1d.dat' )
!    rho = rhozero1d + rhofix1d
!    IF ( ALLOCATED(vzero1d) ) DEALLOCATE(vzero1d)
!    ALLOCATE(vzero1d(n1d))
!    CALL v1d_of_rho1d( n1d, dx, rho, vzero1d )
!    IF ( verbose .GE. 3 ) CALL write_1d( n1d, axis1d, vzero1d, 'vzero1d.dat' )
!    !
!    ! ... Initialize L-BFGS
!    !
!    diagc0 = .FALSE.
!    m = 20
!    n = nsp * n1d
!    length = n*(2*m+1)+2*m
!    iprint(1) = -1
!    iprint(2) = 0
!    IF ( verbose .EQ. 1 ) THEN
!      iprint(1) = 0
!    ENDIF
!    ALLOCATE( g( nsp, n1d ) )
!    ALLOCATE( diag( n ) )
!    ALLOCATE( work( length ) )
!    !
!    ! ... Initial function and gradient evaluations
!    !
!    DO isp = 1, nsp
!       c(isp,:) = SQRT(ABS(czero1d(isp,:)))
!    ENDDO
!    verbose = 4
!    CALL rhoioncc1d_int_fgeval(n1d,nsp,c,f,g)
!    verbose = 2
!    !
!    ! ... Check gradient via finite differences
!    !
!    delta = 0.0000001D0
!    ALLOCATE( gtmp( nsp, n1d ) )
!    DO isp = 1, nsp
!       WRITE(environ_unit,*)
!       DO i = 1, n1d
!          ctmp = c(isp,i)
!          c(isp,i) = ctmp + delta
!          CALL rhoioncc1d_int_fgeval(n1d,nsp,c,f,gtmp)
!          gfd = f
!          c(isp,i) = ctmp - delta
!          CALL rhoioncc1d_int_fgeval(n1d,nsp,c,f,gtmp)
!          gfd = ( gfd - f ) / 2.d0 / delta
!          c(isp,i) = ctmp
!          WRITE(environ_unit,'(4f20.10)')(DBLE(i)-naxis*(DBLE(nrep)+0.5D0))*dx,g(isp,i)-gfd,g(isp,i),gfd
!       ENDDO
!    ENDDO
!    CALL flush(environ_unit)
!    DEALLOCATE( gtmp )
!    !
!!    CALL rhoioncc1d_fgeval(n1d,nsp,c,f,g)
!    !
!    i = 0
!    iflag = 0
!    ALLOCATE( gv( n1d ) )
!    ALLOCATE( x( nsp+1, n1d ) )
!    ALLOCATE( gx( nsp+1, n1d ) )
!    ALLOCATE( gtmp( nsp, n1d ) )
!    ALLOCATE( c2( nsp, n1d ) )
!    DO
!      CALL lbfgs( n, m, c, f, g, diagc0, diag, iprint, tol, xtol, work, iflag)
!      IF ( iflag .LT. 0 ) THEN ! linesearch failed for some reason
!         WRITE(stdout,*)&
!           &'ERROR: lbfgs failed to converge in rhoioncc1d with iflag = ',iflag
!         EXIT
!         STOP
!      ELSE IF ( iflag .EQ. 0 ) THEN ! Calculation converged, exit the loop
!         EXIT
!      ELSE IF ( iflag .EQ. 1 ) THEN ! lbfgs is asking to evaluate f and g
!         i = i + 1
!         CALL rhoioncc1d_fgeval(n1d,nsp,c,ftmp,gtmp)
!         CALL rhoioncc1d_int_fgeval(n1d,nsp,c,f,g)
!         !
!         DO isp = 1 , nsp
!            c2(isp,:) = c(isp,:)**2
!         ENDDO
!         rho = 0.D0
!         CALL sum_cioncc( n1d, nsp, switch1d, z, c2, rho )
!         v = 0.D0
!         CALL v1d_of_rho1d( n1d, dx, rho, v )
!         v = v + vfix1d
!         !
!         CALL vioncc1d_fgeval( n1d, v, fv, gv)
!         CALL vioncc1d_int_fgeval( n1d, v, fvtmp, gv )
!         !
!         DO isp = 1 , nsp
!            x(isp,:) = c(isp,:)
!         ENDDO
!         x(nsp+1,:) = v(:)
!         CALL fullioncc1d_int_fgeval(n1d,nsp,x,fx,gx)
!         !
!         write(environ_unit,'(x,a,i8,5g20.10)')'rhoioncc_fgeval = ',i,f,ftmp,fv,fvtmp,fx
!         CALL flush(environ_unit)
!      ELSE IF ( iflag .GT. 1 ) THEN ! not supposed to happen with diagc0=.FALSE.
!         WRITE(stdout,*)&
!           &'ERROR: unexpected flag from lbfgs optimization, iflag = ',iflag
!         STOP
!      ENDIF
!    ENDDO
!    DEALLOCATE( c2 )
!    DEALLOCATE( gv )
!    DEALLOCATE( gtmp )
!    DEALLOCATE( gx, x )
!    !
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, c(1,:), 'c1_final.dat' )
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, c(2,:), 'c2_final.dat' )
!    !
!    verbose = 4
!    CALL rhoioncc1d_int_fgeval(n1d,nsp,c,f,g)
!    verbose = 2
!    !
!    ! ... Check gradient via finite differences
!    !
!    delta = 0.00000001D0
!    ALLOCATE( gtmp( nsp, n1d ) )
!    DO isp = 1, nsp
!       WRITE(environ_unit,*)
!       DO j = 1, n1d
!          ctmp = c(isp,j)
!          c(isp,j) = ctmp + delta
!          CALL rhoioncc1d_int_fgeval(n1d,nsp,c,f,gtmp)
!          gfd = f
!          c(isp,j) = ctmp - delta
!          CALL rhoioncc1d_int_fgeval(n1d,nsp,c,f,gtmp)
!          gfd = ( gfd - f ) / 2.d0 / delta
!          c(isp,j) = ctmp
!          WRITE(environ_unit,'(4f20.10)')(DBLE(j)-naxis*(DBLE(nrep)+0.5D0))*dx,g(isp,j)-gfd,g(isp,j),gfd
!       ENDDO
!    ENDDO
!    CALL flush(environ_unit)
!    DEALLOCATE( gtmp )
!    !
!    ALLOCATE( c2( nsp, n ) )
!    DO isp = 1, nsp
!       c2(isp,:) = c(isp,:)**2
!    ENDDO
!    rho = 0.D0
!    CALL sum_cioncc( n1d, nsp, switch1d, z, c2, rho )
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, rho, 'rhoioncc1d.dat' )
!    DEALLOCATE( c2 )
!    !
!    CALL v1d_of_rho1d( n1d, dx, rho, v )
!    v = v + vfix1d
!    !
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, v, 'v1d.dat' )
!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, v, 'v1d_first.dat' )
!!    v = v - vfix1d - vfix_avg + vfix1d_avg
!!    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, v, 'vext1d.dat' )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE calc_rhoioncc1d_lbfgs
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE calc_vioncc1d_lbfgs(  n, nsp, v, c )
!!--------------------------------------------------------------------
!    IMPLICIT NONE
!    !
!    REAL( DP ), PARAMETER :: tol = 1.D-8, xtol = 1.D-14
!    !
!    INTEGER, INTENT(IN) :: n, nsp
!    !
!    REAL( DP ), INTENT(OUT) :: v(n)
!    REAL( DP ), INTENT(INOUT) :: c(nsp,n)
!    !
!    ! ... Local variables
!    !
!    LOGICAL :: diagc0
!    INTEGER :: m, length, iflag, iprint(2), nfdpoint, i, j, isp
!    REAL( DP ) :: charge, f, gfd, delta, vtmp, ftmp, fc, fctmp, fx
!    REAL( DP ), DIMENSION( :, : ), ALLOCATABLE :: gc, gx, x, c2
!    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: g, gtmp, diag, work
!    !
!    REAL( DP ), DIMENSION( n ) :: rho
!    !
!    ! ... Build the constant approximation of rhoioncc, rescale the guess
!    !
!!    rho = 0.D0
!!    ALLOCATE( c2( nsp, n ) )
!!    DO isp = 1, nsp
!!       c2( isp, : ) = c(isp,:)**2
!!    ENDDO
!!    CALL sum_cioncc( n, nsp, switch1d, z, c2, rho )
!!    !
!!    IF ( ALLOCATED(rhozero1d) ) DEALLOCATE(rhozero1d)
!!    ALLOCATE(rhozero1d(n))
!!    rhozero1d(:) = rho(:)
!!    charge = SUM(rhozero1d(:))*darea
!!    rhozero1d = rhozero1d*ABS(charge_fix/charge)
!!    IF ( verbose .GE. 3 ) CALL write_1d( n, axis1d, rhozero1d, 'rhozero1d.dat' )
!!    rho = rhozero1d + rhofix1d
!!    IF ( ALLOCATED(vzero1d) ) DEALLOCATE(vzero1d)
!!    ALLOCATE(vzero1d(n))
!!    CALL v1d_of_rho1d( n, dx, rho, vzero1d )
!!    IF ( verbose .GE. 3 ) CALL write_1d( n, axis1d, vzero1d, 'vzero1d.dat' )
!    !
!    ! ... Initialize L-BFGS
!    !
!    diagc0 = .FALSE.
!    m = 20
!    length = n*(2*m+1)+2*m
!    iprint(1) = -1
!    iprint(2) = 0
!    IF ( verbose .EQ. 1 ) THEN
!      iprint(1) = 0
!!    ELSE IF ( verbose .GT. 1 ) THEN
!!      iprint(1) = 1
!    ENDIF
!    ALLOCATE( g( n ) )
!    ALLOCATE( diag( n ) )
!    ALLOCATE( work( length ) )
!    !
!    ! ... Initial function and gradient evaluations
!    !
!    v = 0.D0 !v ! 0.D0 !v - vfix1d
!    CALL write_1d( n, axis1d, v, 'v_start.dat' )
!!DEBUG!    verbose = 4
!    CALL vioncc1d_fgeval(n,v,f,g)
!!DEBUG!    verbose = 2
!!DEBUG!    !
!!DEBUG!    ! ... Check gradient via finite differences
!!DEBUG!    !
!!DEBUG!    delta = 0.000001D0
!!DEBUG!    ALLOCATE( gtmp( n ) )
!!DEBUG!    DO i = 1, n
!!DEBUG!      vtmp = v(i)
!!DEBUG!      v(i) = vtmp + delta
!!DEBUG!      CALL vioncc1d_int_fgeval(n,v,f,gtmp)
!!DEBUG!      gfd = f
!!DEBUG!      v(i) = vtmp - delta
!!DEBUG!      CALL vioncc1d_int_fgeval(n,v,f,gtmp)
!!DEBUG!      gfd = ( gfd - f ) / 2.d0 / delta
!!DEBUG!      v(i) = vtmp
!!DEBUG!      WRITE(environ_unit,'(4f20.10)')(DBLE(i)-naxis*(DBLE(nrep)+0.5D0))*dx,g(i)-gfd,g(i),gfd
!!DEBUG!    ENDDO
!!DEBUG!    CALL flush(environ_unit)
!!DEBUG!    DEALLOCATE( gtmp )
!    !
!    i = 0
!    iflag = 0
!    ALLOCATE( gtmp( n ) )  ! DEBUG !
!!DEBUG!    ALLOCATE( gc( nsp, n ) )  ! DEBUG !
!!DEBUG!    ALLOCATE( x( nsp+1, n ) )  ! DEBUG !
!!DEBUG!    ALLOCATE( gx( nsp+1, n ) )  ! DEBUG !
!    DO
!      CALL lbfgs( n, m, v, f, g, diagc0, diag, iprint, tol, xtol, work, iflag)
!      IF ( iflag .LT. 0 ) THEN ! linesearch failed for some reason
!         WRITE(stdout,*)&
!           &'ERROR: lbfgs failed to converge in rhoioncc1d with iflag = ',iflag
!         EXIT!STOP ! DEBUG !
!      ELSE IF ( iflag .EQ. 0 ) THEN ! Calculation converged, exit the loop
!         EXIT
!      ELSE IF ( iflag .EQ. 1 ) THEN ! lbfgs is asking to evaluate f and g
!         i = i + 1
!         CALL vioncc1d_fgeval(n,v,f,g) ! DEBUG !
!         CALL vioncc1d_int_fgeval(n,v,ftmp,gtmp)
!         !
!!DEBUG!         CALL generate_cioncc( n, nsp, 1, switch1d, v, c ) ! DEBUG !
!!DEBUG!         DO isp = 1 , nsp ! DEBUG !
!!DEBUG!            c(isp,:) = SQRT(ABS(c(isp,:))) ! DEBUG !
!!DEBUG!            x(isp,:) = c(isp,:) ! DEBUG !
!!DEBUG!         ENDDO ! DEBUG !
!!DEBUG!         CALL rhoioncc1d_fgeval(n, nsp, c, fc, gc) ! DEBUG !
!!DEBUG!         CALL rhoioncc1d_int_fgeval(n, nsp, c, fctmp, gc)   ! DEBUG !
!!DEBUG!         x(nsp+1,:) = v(:) ! DEBUG !
!!DEBUG!         CALL fullioncc1d_int_fgeval(n,nsp,x,fx,gx) ! DEBUG !
!         !
!!DEBUG!         WRITE(environ_unit,'(x,a,i8,5g20.10)')'vioncc_fgeval = ',i,fc,fctmp,ftmp,f,fx
!         WRITE(environ_unit,'(1X,a,i8,2g20.10)')'vioncc_fgeval = ',i,f,ftmp
!         CALL flush(environ_unit)
!      ELSE IF ( iflag .GT. 1 ) THEN ! not supposed to happen with diagc0=.FALSE.
!         WRITE(stdout,*)&
!           &'ERROR: unexpected flag from lbfgs optimization, iflag = ',iflag
!         RETURN!STOP
!      ENDIF
!    ENDDO
!    DEALLOCATE( gtmp ) ! DEBUG !
!!DEBUG!    DEALLOCATE( gc ) ! DEBUG !
!!DEBUG!    DEALLOCATE( gx ) ! DEBUG !
!!DEBUG!    DEALLOCATE( x ) ! DEBUG !
!!DEBUG!    !
!!DEBUG!    verbose = 4
!!DEBUG!    CALL vioncc1d_int_fgeval(n,v,f,g)
!!DEBUG!    verbose = 2
!!DEBUG!    !
!!DEBUG!    ! ... Check gradient via finite differences
!!DEBUG!    !
!!DEBUG!    delta = 0.000001D0
!!DEBUG!    ALLOCATE( gtmp( n ) )
!!DEBUG!    DO i = 1, n
!!DEBUG!      vtmp = v(i)
!!DEBUG!      v(i) = vtmp + delta
!!DEBUG!      CALL vioncc1d_int_fgeval(n,v,f,gtmp)
!!DEBUG!      gfd = f
!!DEBUG!      v(i) = vtmp - delta
!!DEBUG!      CALL vioncc1d_int_fgeval(n,v,f,gtmp)
!!DEBUG!      gfd = ( gfd - f ) / 2.d0 / delta
!!DEBUG!      v(i) = vtmp
!!DEBUG!      WRITE(environ_unit,'(4f20.10)')(DBLE(i)-naxis*(DBLE(nrep)+0.5D0))*dx,g(i)-gfd,g(i),gfd
!!DEBUG!    ENDDO
!!DEBUG!    CALL flush(environ_unit)
!!DEBUG!    DEALLOCATE( gtmp )
!!DEBUG!    !
!!    v = v + vzero1d
!!    !
!    CALL generate_cioncc( n, nsp, 1, switch1d, v, c )
!    !
!    IF ( verbose .GE. 2 ) CALL write_1d( n, axis1d, v, 'v1d.dat' )
!    IF ( verbose .GE. 2 ) CALL write_1d( n, axis1d, v, 'v1d_final.dat' )
!    v = v - vfix1d - vfix_avg + vfix1d_avg
!    IF ( verbose .GE. 2 ) CALL write_1d( n, axis1d, v, 'vext1d.dat' )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE calc_vioncc1d_lbfgs
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE setlocal( nnr, nspin, rhoin )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT( IN ) :: nnr, nspin
!    REAL( DP ), DIMENSION( nnr ), INTENT( IN ) :: rhoin
!    !
!    ! ... Store local quantities
!    !
!    nnrlocal = nnr
!    nspinlocal = nspin
!    IF ( ALLOCATED(rhoinlocal) ) DEALLOCATE(rhoinlocal)
!    ALLOCATE( rhoinlocal( nnr ) )
!    rhoinlocal = rhoin
!    IF ( ALLOCATED(rhooutlocal) ) DEALLOCATE(rhooutlocal)
!    ALLOCATE( rhooutlocal( nnr ) )
!    rhooutlocal = 0.D0
!    IF ( ALLOCATED(vtotlocal) ) DEALLOCATE(vtotlocal)
!    ALLOCATE( vtotlocal( nnr, nspin ) )
!    vtotlocal = 0.D0
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE setlocal
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE cleanlocal( )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    ! ... Clean local quantities
!    !
!    IF ( ALLOCATED( rhoinlocal  ) ) DEALLOCATE( rhoinlocal  )
!    IF ( ALLOCATED( rhooutlocal ) ) DEALLOCATE( rhooutlocal )
!    IF ( ALLOCATED( vtotlocal   ) ) DEALLOCATE( vtotlocal   )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE cleanlocal
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE ioncc_fgeval(n,v,f,g)
!!--------------------------------------------------------------------
!    !
!    USE fd_gradient,   ONLY : calc_fd_gradient, calc_fd_laplacian
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: n
!    REAL( DP ), DIMENSION(n), INTENT(IN) :: v
!    REAL( DP ), INTENT(OUT) :: f
!    REAL( DP ), DIMENSION(n), INTENT(INOUT) :: g
!    !
!    INTEGER :: i, j
!    REAL( DP ) :: df, d2f
!    REAL( DP ), DIMENSION(n) :: residue, dresidue, laplv, vtmp
!    !
!!    CALL write_cube( n, v, 'v_in_fgeval.cube' )
!    !
!    ! ...The error is computed from the norm of the residue
!    !    which is taken to be the difference of d/dx(eps(x)*d/dx\phi(x))
!    !    from -4*pi*e2*(rhofix+rhoioncc)
!    vtmp = v !+ vfix
!    !
!    CALL generate_rhoioncc( n, 1, switch, vtmp, residue )
!    !
!    residue = fpi * e2 * residue
!    residue = residue + fpi * e2 * rhofix
!!    CALL write_cube( n, residue, 'residue_1.cube' )
!    !
!    ! ...The gradient of the error is 2 * residue \cdot (d/d\phi residue)
!    !
!    CALL generate_drhoioncc( n, 1, switch, vtmp, dresidue )
!    !
!    dresidue = fpi * e2 * dresidue
!!    CALL write_cube( n, dresidue, 'dresidue.cube' )
!    !
!    ! ...Compute the first part of the residue vector and the total error
!    !
!    ! need lapl(v) * eps + grad(eps) \cdot grad(v)
!    CALL external_laplacian( v, laplv )
!!    CALL calc_fd_laplacian( nfdpoint, icfd2, ncfd2, n, v, laplv)
!!    CALL write_cube( n, laplv, 'residue_2.cube' )
!    !
!    residue = residue + laplv !*eps + gradepsgradv
!!    CALL write_cube( n, residue, 'residue_tot.cube' )
!    !
!    f = SUM(residue**2)
!    CALL mp_sum( f, intra_bgrp_comm )
!    !
!    ! ... Compute the gradient of the residue wrt the components of v
!    !
!    CALL external_laplacian( residue, g )
!!    CALL calc_fd_laplacian( nfdpoint, icfd2, ncfd2, n, residue, g )
!!    CALL write_cube( n, g, 'g.cube' )
!    g(:) = g(:) + residue(:) * dresidue(:)
!    g = g * 2.D0
!    !
!    CALL flush(environ_unit)
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE ioncc_fgeval
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE ioncc_heval(n,v,diag)
!!--------------------------------------------------------------------
!    !
!    USE fd_gradient,   ONLY : calc_fd_gradient, calc_fd_laplacian
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: n
!    REAL( DP ), DIMENSION(n), INTENT(IN) :: v
!    REAL( DP ), DIMENSION(n), INTENT(OUT) :: diag
!    !
!    INTEGER :: i, j
!    REAL( DP ) :: df, d2f, fact
!    REAL( DP ), DIMENSION(n) :: residue, dresidue, d2residue, laplv
!    !
!!    CALL write_cube( n, v, 'v_in_fgeval.cube' )
!    !
!    ! ...The error is computed from the norm of the residue
!    !    which is taken to be the difference of d/dx(eps(x)*d/dx\phi(x))
!    !    from -4*pi*e2*(rhofix+rhoioncc)
!    !
!    CALL generate_rhoioncc( n, 1, switch, v, residue )
!    !
!    residue = fpi * e2 * ( residue + rhofix )
!!    CALL write_cube( n, residue, 'residue_1.cube' )
!    !
!    ! ...Compute the first part of the residue vector and the total error
!    !
!    ! need lapl(v) * eps + grad(eps) \cdot grad(v)
!!    CALL calc_fd_gradient()
!    CALL calc_fd_laplacian( nfdpoint, icfd2, ncfd2, n, v, laplv)
!!    CALL write_cube( n, laplv, 'residue_2.cube' )
!    !
!    residue = residue + laplv !*eps + gradepsgradv
!!    CALL write_cube( n, residue, 'residue_tot.cube' )
!    !
!    ! ...The hessian of the error is approximately 2 * residue \cdot (d/d\phi residue)
!    !
!    CALL generate_drhoioncc( n, 1, switch, v, dresidue )
!    !
!    dresidue = fpi * e2 * dresidue
!!    CALL write_cube( n, dresidue, 'dresidue.cube' )
!    !
!    ! ... second derivative
!    !
!    CALL generate_d2rhoioncc( n, 1, switch, v, d2residue )
!    !
!    d2residue = fpi * e2 * d2residue
!!    CALL write_cube( n, d2residue, 'd2residue.cube' )
!    !
!    fact = (nr1/at(1,1)/alat)**2 + (nr2/at(2,2)/alat)**2 + (nr3/at(3,3)/alat)**2
!    diag(:) = (icfd2(0)/DBLE(ncfd2)*fact + dresidue(:))**2 + &
!         & residue(:)*d2residue(:)
!    !
!    diag = 0.5D0 / diag(i)
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE ioncc_heval
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE ioncc_int_fgeval(n,nsp,v,f,g)
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: n, nsp
!    REAL( DP ), DIMENSION(n), INTENT(IN) :: v
!    REAL( DP ), INTENT(OUT) :: f
!    REAL( DP ), DIMENSION(n), INTENT(OUT) :: g
!    !
!    REAL( DP ), DIMENSION(n) :: residue, laplv, rho, s, vtmp
!    REAL( DP ), DIMENSION(3, n) :: gradv
!    REAL( DP ), DIMENSION(nsp, n) :: c
!    !
!    INTEGER :: i
!    !
!    CALL external_gradient( v, gradv )
!    !
!    vtmp = v !+ vzero
!    CALL generate_cioncc( n, nsp, 1, switch, vtmp, c )
!    !
!    CALL sum_cioncc( n, nsp, switch, z, c, rho )
!    !
!    CALL calc_sioncc( n, nsp, switch, c, s )
!    !
!    CALL sum_cioncc( n, nsp, switch, mu, c, residue )
!    !
!    DO i = 1, n
!       residue(i) = residue(i) + 0.5D0 / fpi / e2 * SUM(gradv(:,i)**2) - ( rhofix(i) + rho(i) ) * v(i) + s(i)
!      ! residue(i) = residue(i) + 1.D0 / fpi / e2 * SUM(gradv(:,i)*gradvzero(:,i)) - (rho(i)-rhozero(i)) * vzero(i)
!    ENDDO
!    !
!    f = SUM(residue)*domega
!    CALL mp_sum( f, intra_bgrp_comm )
!    !
!    CALL external_laplacian( v, laplv )
!    !
!    IF ( verbose .EQ. 4 ) THEN
!       WRITE(environ_unit,*)'gradv**2'
!       WRITE(environ_unit,*)
!       DO i = 1, n, nr1*nr2
!          WRITE(environ_unit,*)i/nr1/nr2*dx,0.5D0 / fpi / e2 * SUM(gradv(:,i)**2)
!       ENDDO
!       WRITE(environ_unit,*)
!       !
!       WRITE(environ_unit,*)'rho'
!       WRITE(environ_unit,*)
!       DO i = 1, n, nr1*nr2
!          WRITE(environ_unit,*)i/nr1/nr2*dx,rho(i)
!       ENDDO
!       WRITE(environ_unit,*)
!       !
!       WRITE(environ_unit,*)'residue'
!       WRITE(environ_unit,*)
!       DO i = 1, n, nr1*nr2
!          WRITE(environ_unit,*)i/nr1/nr2*dx,residue(i)
!       ENDDO
!       WRITE(environ_unit,*)
!       !
!       WRITE(environ_unit,*)'entropy'
!       WRITE(environ_unit,*)
!       DO i = 1, n, nr1*nr2
!          WRITE(environ_unit,*)i/nr1/nr2*dx,s(i)
!       ENDDO
!       WRITE(environ_unit,*)
!       !
!       WRITE(environ_unit,*)'laplv'
!       WRITE(environ_unit,*)
!       DO i = 1, n, nr1*nr2
!          WRITE(environ_unit,*)i/nr1/nr2*dx,laplv(i)
!       ENDDO
!       WRITE(environ_unit,*)
!       FLUSH(environ_unit)
!    ENDIF
!    !
!    g(:) = laplv(:) / fpi / e2 + ( rho(:) + rhofix(:) )
!   ! g(:) = g(:) + laplvzero(:) / fpi / e2
!    g = - g * domega
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE ioncc_int_fgeval
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE ioncc_mix_fgeval(n,nsp,alpha,v,f,g)
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: n, nsp
!    REAL( DP ), INTENT(IN) :: alpha
!    REAL( DP ), DIMENSION(n), INTENT(IN) :: v
!    REAL( DP ), INTENT(OUT) :: f
!    REAL( DP ), DIMENSION(n), INTENT(OUT) :: g
!    !
!    REAL( DP ), DIMENSION(n) :: residue, dresidue, laplv, rho, s
!    REAL( DP ), DIMENSION(3, n) :: gradv
!    REAL( DP ), DIMENSION(nsp, n) :: c
!    !
!    INTEGER :: i
!    !
!    CALL external_gradient( v, gradv )
!    !
!    CALL generate_cioncc( n, nsp, 1, switch, v, c )
!    !
!    CALL sum_cioncc( n, nsp, switch, z, c, rho )
!    !
!    CALL calc_sioncc( n, nsp, switch, c, s )
!    !
!    CALL sum_cioncc( n, nsp, switch, mu, c, residue )
!    !
!    DO i = 1, n
!       residue(i) = residue(i) + 0.5D0 / fpi / e2 * SUM(gradv(:,i)**2) - ( rhofix(i) + rho(i) ) * v(i) + s(i)
!    ENDDO
!    !
!    f = SUM(residue)*domega
!    !
!    CALL external_laplacian( v, laplv )
!    !
!    residue = laplv(:) / fpi / e2 + ( rho(:) + rhofix(:) )
!    !
!    g = 0.D0
!    !
!    IF ( alpha .GT. 0.D0 ) THEN
!       !
!       f = f + SUM(residue**2) * alpha
!       !
!       CALL generate_drhoioncc( n, 1, switch, v, dresidue )
!       !
!       CALL external_laplacian( residue, g )
!       !
!       g(:) = 2.D0 * alpha * ( g(:) + residue * dresidue )
!       !
!    ENDIF
!    !
!    g = g - residue * domega
!    !
!    CALL mp_sum( f, intra_bgrp_comm )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE ioncc_mix_fgeval
!!--------------------------------------------------------------------
!!!!OLD!!!!--------------------------------------------------------------------
!!!!OLD!!!  SUBROUTINE ioncc_int_fgeval(n,v,f,g)
!!!!OLD!!!!--------------------------------------------------------------------
!!!!OLD!!!    !
!!!!OLD!!!    USE fd_gradient,   ONLY : calc_fd_gradient, calc_fd_laplacian
!!!!OLD!!!    !
!!!!OLD!!!    IMPLICIT NONE
!!!!OLD!!!    !
!!!!OLD!!!    INTEGER, INTENT(IN) :: n
!!!!OLD!!!    REAL( DP ), DIMENSION(n), INTENT(IN) :: v
!!!!OLD!!!    REAL( DP ), INTENT(OUT) :: f
!!!!OLD!!!    REAL( DP ), DIMENSION(n), INTENT(OUT) :: g
!!!!OLD!!!    !
!!!!OLD!!!    INTEGER :: i, j
!!!!OLD!!!    REAL( DP ) :: df, d2f, charge, vavg
!!!!OLD!!!    REAL( DP ), DIMENSION(n) :: rhotmp, residue, dresidue, vtmp
!!!!OLD!!!    REAL( DP ), DIMENSION(3,n) :: gradv, gradtmp
!!!!OLD!!!    !
!!!!OLD!!!!    CALL write_cube( n, v, 'v_in_fgeval.cube' )
!!!!OLD!!!    !
!!!!OLD!!!    ! ...The error is computed from the total free energy (only the
!!!!OLD!!!    !    terms that depend on the electrostatic potential)
!!!!OLD!!!    !
!!!!OLD!!!    vtmp = v !+ vfix
!!!!OLD!!!    !
!!!!OLD!!!    vavg = SUM(vtmp)/n
!!!!OLD!!!    CALL mp_sum( vavg, intra_bgrp_comm )
!!!!OLD!!!    WRITE(environ_unit,*)'vavg = ',vavg
!!!!OLD!!!    !
!!!!OLD!!!    CALL calc_fd_gradient( nfdpoint, icfd, ncfd, n, vtmp, gradv)
!!!!OLD!!!    !
!!!!OLD!!!    CALL generate_rhoioncc( n, 1, switch, vtmp, rhotmp )
!!!!OLD!!!    !
!!!!OLD!!!    charge = SUM(rhofix+rhotmp)*domega
!!!!OLD!!!    CALL mp_sum( charge, intra_bgrp_comm )
!!!!OLD!!!    WRITE(environ_unit,*)'charge = ',charge
!!!!OLD!!!    DO i = 1, n
!!!!OLD!!!       residue(i) = -0.5D0*SUM(gradv(:,i)**2)/fpi/e2 + &
!!!!OLD!!!            & (rhofix(i)+rhotmp(i))*(vtmp(i)-vavg)
!!!!OLD!!!!            & rhofix(i)*(vtmp(i)-vavg)
!!!!OLD!!!    ENDDO
!!!!OLD!!!    !
!!!!OLD!!!    f = SUM(residue)*domega
!!!!OLD!!!    CALL mp_sum( f, intra_bgrp_comm )
!!!!OLD!!!    WRITE(environ_unit,*)'fgeval: f = ',f
!!!!OLD!!!    !
!!!!OLD!!!    ! ... Compute the gradient wrt the components of v
!!!!OLD!!!    !
!!!!OLD!!!    ! CALL calc_fd_laplacian( nfdpoint, icfd2, ncfd2, n, vtmp, g )
!!!!OLD!!!    g = 0.D0
!!!!OLD!!!    DO i = 1, 3
!!!!OLD!!!       CALL calc_fd_gradient( nfdpoint, icfd, ncfd, n, gradv(i,:), gradtmp )
!!!!OLD!!!       g(:) = g(:) + gradtmp(i,:)
!!!!OLD!!!    ENDDO
!!!!OLD!!!    !
!!!!OLD!!!    CALL generate_drhoioncc( n, 1, switch, vtmp, dresidue )
!!!!OLD!!!    !
!!!!OLD!!!    g = g/fpi/e2 + rhofix - charge / omega + rhotmp + dresidue*(vtmp-vavg)
!!!!OLD!!!    g = g*domega
!!!!OLD!!!    !
!!!!OLD!!!    CALL flush(environ_unit)
!!!!OLD!!!    !
!!!!OLD!!!    RETURN
!!!!OLD!!!    !
!!!!OLD!!!!--------------------------------------------------------------------
!!!!OLD!!!  END SUBROUTINE ioncc_int_fgeval
!!!!OLD!!!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE calc_vioncc_lbfgs( n, nsp, rho, v )
!!--------------------------------------------------------------------
!    !
!    USE scatter_mod,    ONLY: scatter_grid, gather_grid
!    USE fft_base,       ONLY: dfftp
!    USE random_numbers, ONLY: randy
!    USE mp,             ONLY: mp_bcast
!    USE mp_bands,       ONLY: me_bgrp, root_bgrp
!    !
!    IMPLICIT NONE
!    !
!    REAL( DP ), PARAMETER :: tol = 1.D-10, xtol = 1.D-10
!    !
!    INTEGER, INTENT(IN) :: n, nsp
!    !
!    REAL( DP ), INTENT(INOUT) :: rho(n)
!    REAL( DP ), INTENT(OUT) :: v(n)
!    !
!    ! ... Local variables
!    !
!    LOGICAL :: diagc0
!    INTEGER :: m, length, iflag, iprint(2)
!    REAL( DP ) :: f
!    REAL( DP ), DIMENSION( n ) :: g, diag
!    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: gtot, diagtot, vtot, work
!    !
!    INTEGER :: i, ntest, itest, istep
!    REAL( DP ) :: vtmp, delta, gfd, hfd, ftmp, alpha
!    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: gtmp
!    !
!    ! ... Initialize L-BFGS
!    !
!    diagc0 = .FALSE.
!    m = 10
!    length = ntot*(2*m+1)+2*m
!    iprint(1) = -1
!    iprint(2) = 0
!    IF ( verbose .EQ. 1 ) THEN
!      iprint(1) = 0
!    ENDIF
!    ALLOCATE( work( length ) )
!    ALLOCATE( vtot( ntot ) )
!    ALLOCATE( gtot( ntot ) )
!    ALLOCATE( diagtot( ntot ) )
!    !
!    ! ... Initial function and gradient evaluations
!    !
!    IF ( ALLOCATED(vzero) ) DEALLOCATE(vzero)
!    ALLOCATE(vzero(n))
!    vzero = v
!    !
!    IF ( ALLOCATED(gradvzero) ) DEALLOCATE(gradvzero)
!    ALLOCATE(gradvzero(3,n))
!    CALL external_gradient( vzero, gradvzero )
!    !
!    IF ( ALLOCATED(rhozero) ) DEALLOCATE(rhozero)
!    ALLOCATE(rhozero(n))
!    CALL generate_rhoioncc( n, 1, switch, vzero, rhozero )
!    !
!    IF ( ALLOCATED(laplvzero) ) DEALLOCATE(laplvzero)
!    ALLOCATE(laplvzero(n))
!    CALL external_laplacian( vzero, laplvzero )
!    !
!    v = 0.D0! v + vfix
!    verbose = 4
!    alpha = 1.D-6
!    CALL ioncc_mix_fgeval(n,nsp,alpha,v,f,g)
!    verbose = 2
!#if defined (__MPI)
!    CALL gather_grid ( dfftp, v, vtot )
!    CALL gather_grid ( dfftp, g, gtot )
!    CALL gather_grid ( dfftp, diag, diagtot )
!#else
!    vtot = v
!    gtot = g
!    diagtot = diag
!#endif
!    WRITE(environ_unit,*)'vzero'
!    WRITE(environ_unit,*)SUM(v)/n
!    DO i = 1, n, nr1*nr2
!       WRITE(environ_unit,*)i/nr1/nr2*dx,vzero(i)
!    ENDDO
!    WRITE(environ_unit,*)
!    WRITE(environ_unit,*)'rhozero'
!    WRITE(environ_unit,*)
!    DO i = 1, n, nr1*nr2
!       WRITE(environ_unit,*)i/nr1/nr2*dx,rhozero(i)
!    ENDDO
!    WRITE(environ_unit,*)
!    WRITE(environ_unit,*)'v1'
!    WRITE(environ_unit,*)SUM(v)/n
!    DO i = 1, n, nr1*nr2
!       WRITE(environ_unit,*)i/nr1/nr2*dx,v(i)
!    ENDDO
!    WRITE(environ_unit,*)
!!DEBUG!SINGLE-PROC~    !
!!DEBUG!SINGLE-PROC~    ! ... Check gradient via finite differences
!!DEBUG!SINGLE-PROC~    !
!!DEBUG!SINGLE-PROC~    delta = 0.000001D0
!!DEBUG!SINGLE-PROC~    ALLOCATE( gtmp( n ) )
!!DEBUG!SINGLE-PROC~    DO i = 1, n, nr1*nr2
!!DEBUG!SINGLE-PROC~       vtmp = v(i)
!!DEBUG!SINGLE-PROC~       v(i) = vtmp + delta
!!DEBUG!SINGLE-PROC~       CALL ioncc_int_fgeval(n,nsp,v,f,gtmp)
!!DEBUG!SINGLE-PROC~       gfd = f
!!DEBUG!SINGLE-PROC~       v(i) = vtmp - delta
!!DEBUG!SINGLE-PROC~       CALL ioncc_int_fgeval(n,nsp,v,f,gtmp)
!!DEBUG!SINGLE-PROC~       gfd = ( gfd - f ) / 2.d0 / delta
!!DEBUG!SINGLE-PROC~       v(i) = vtmp
!!DEBUG!SINGLE-PROC~       WRITE(environ_unit,'(4f20.10)')(i-1)/nr1/nr2*dx,g(i)-gfd,g(i),gfd
!!DEBUG!SINGLE-PROC~       FLUSH(environ_unit)
!!DEBUG!SINGLE-PROC~    ENDDO
!!DEBUG!SINGLE-PROC~    DEALLOCATE( gtmp )
!    !
!!    IF ( diagc0 ) CALL ioncc_int_heval(n,v,diag)
!    !
!    ! ... Check diagonal inverse hessian via finite differences
!    !
!!!!!    delta = 0.00001D0
!!!!!    ALLOCATE( gtmp( n ) )
!!!!!    ntest = 1000
!!!!!    DO i = 1, ntest
!!!!!      itest = INT( n * randy() ) + 1
!!!!!      vtmp = v(itest)
!!!!!      v(itest) = vtmp + delta
!!!!!      CALL ioncc_fgeval(n,v,f,gtmp)
!!!!!      hfd = gtmp(itest)
!!!!!      v(itest) = vtmp - delta
!!!!!      CALL ioncc_fgeval(n,v,f,gtmp)
!!!!!      hfd = ( hfd - gtmp(itest) ) / 2.d0 / delta
!!!!!      hfd = 1.D0 / hfd
!!!!!      v(itest) = vtmp
!!!!!      WRITE(environ_unit,'(i8,3f20.10)')itest,diag(itest)-hfd,diag(itest),hfd
!!!!!      CALL flush(environ_unit)
!!!!!    ENDDO
!!!!!    DEALLOCATE( gtmp )
!    !
!    iflag = 0
!    istep = 0
!    ALLOCATE(gtmp(n))
!    DO
!      IF ( me_bgrp .EQ. root_bgrp ) CALL lbfgs( ntot, m, vtot, f, gtot, diagc0, diagtot, iprint, tol, xtol, work, iflag)
!#if defined (__MPI)
!      CALL mp_bcast( iflag, root_bgrp, intra_bgrp_comm )
!      CALL scatter_grid ( dfftp, vtot, v )
!#else
!      v = vtot
!#endif
!      IF ( iflag .LT. 0 ) THEN ! linesearch failed for some reason
!         WRITE(stdout,*)&
!           &'ERROR: lbfgs failed to converge in rhoioncc1d with iflag = ',iflag
!         EXIT!RETURN!STOP
!      ELSE IF ( iflag .EQ. 0 ) THEN ! Calculation converged, exit the loop
!         EXIT
!      ELSE IF ( iflag .EQ. 1 ) THEN ! lbfgs is asking to evaluate f and g
!         istep = istep + 1
!         CALL ioncc_fgeval(n,v,ftmp,gtmp)
!         CALL ioncc_mix_fgeval(n,nsp,alpha,v,f,g)
!         WRITE(environ_unit,'(1X,a,i10,2g30.20)')'fgeval = ',istep,f,ftmp
!         FLUSH(environ_unit)
!      ELSE IF ( iflag .EQ. 2 ) THEN ! lbfgs is asking to estimate h^-1
!!         CALL ioncc_heval(n,v,diag)
!      ELSE IF ( iflag .GT. 2 ) THEN ! not supposed to happen
!         WRITE(stdout,*)&
!           &'ERROR: unexpected flag from lbfgs optimization, iflag = ',iflag
!         RETURN!STOP
!      ENDIF
!#if defined (__MPI)
!      CALL gather_grid ( dfftp, g, gtot )
!      IF ( iflag .EQ. 2 ) CALL gather_grid ( dfftp, diag, diagtot )
!#else
!      gtot = g
!      IF ( iflag .EQ. 2 ) diagtot = diag
!#endif
!    ENDDO
!    DEALLOCATE(gtmp)
!    CALL write_cube( n, v, 'vtot.cube' )
!    !
!    ! | | | ... DEBUG ... | | |
!    ! V V V               V V V
!    WRITE(environ_unit,*)'vloc'
!    WRITE(environ_unit,*)SUM(v)/n
!    DO i = 1, n, nr1*nr2
!       WRITE(environ_unit,*)(i-1)/nr1/nr2*dx,v(i)
!    ENDDO
!    WRITE(environ_unit,*)
!    FLUSH(environ_unit)
!    !
!    verbose = 4
!    CALL ioncc_mix_fgeval(n,nsp,alpha,v,f,g)
!    verbose = 2
!    !
!    delta = 0.000001D0
!    ALLOCATE( gtmp( n ) )
!    DO i = 1, n, nr1*nr2
!       vtmp = v(i)
!       v(i) = vtmp + delta
!       CALL ioncc_mix_fgeval(n,nsp,alpha,v,f,gtmp)
!       gfd = f
!       v(i) = vtmp - delta
!       CALL ioncc_mix_fgeval(n,nsp,alpha,v,f,gtmp)
!       gfd = ( gfd - f ) / 2.d0 / delta
!       v(i) = vtmp
!       WRITE(environ_unit,'(4f20.10)')(i-1)/nr1/nr2*dx,g(i)-gfd,g(i),gfd
!       FLUSH(environ_unit)
!    ENDDO
!    DEALLOCATE( gtmp )
!    ! ^ ^ ^                   ^ ^ ^
!    ! | | | ... END DEBUG ... | | |
!    !
!    CALL generate_rhoioncc( n, 1, switch, v, rho )
!    CALL write_cube( n, rho, 'rhoioncc.cube' )
!    !
!    v = v - vfix
!    CALL write_cube( n, v, 'vioncc.cube' )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE calc_vioncc_lbfgs
!!--------------------------------------------------------------------
!!!!!!!--------------------------------------------------------------------
!!!!!!  FUNCTION ioncc_feval(x)
!!!!!!!--------------------------------------------------------------------
!!!!!!    USE kinds, ONLY: DP
!!!!!!    !
!!!!!!    IMPLICIT NONE
!!!!!!    !
!!!!!!    REAL( DP ), INTENT(IN) :: x
!!!!!!    REAL( DP ) :: ioncc_feval
!!!!!!    !
!!!!!!    ! ... Local variables
!!!!!!    !
!!!!!!    REAL( DP ) :: charge, ehart
!!!!!!    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: vioncc1d, rhoioncc1d
!!!!!!    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: vzero, rhozero, vperiodic
!!!!!!    !
!!!!!!    ! ... Allocate
!!!!!!    !
!!!!!!    ALLOCATE( vioncc1d( n1d ) )
!!!!!!    vioncc1d = 0.D0
!!!!!!    ALLOCATE( rhoioncc1d( n1d ) )
!!!!!!    rhoioncc1d = 0.D0
!!!!!!    ALLOCATE( vzero( nnrlocal ) )
!!!!!!    ALLOCATE( rhozero( nnrlocal ) )
!!!!!!    ALLOCATE( vperiodic( nnrlocal ) )
!!!!!!    !
!!!!!!    ! ... Given x, compute the 1d analytic solution
!!!!!!    !
!!!!!!    CALL calc_rhostern( n1d, charge_fix, dipole_fix, x, axis1d, step1d, vfix1d, vioncc1d, rhoioncc1d )
!!!!!!    !
!!!!!!    ! ... Compute the new switching function
!!!!!!    !
!!!!!!    CALL generate_switch( n1d, x, deltastern, axis1d, switch1d )
!!!!!!    !
!!!!!!    ! ... Compute the 1d numerical solution
!!!!!!    !
!!!!!!    rhoioncc1d = rhoioncc1d * switch1d
!!!!!!    CALL calc_rhoioncc1d_lbfgs( n1d, vioncc1d, rhoioncc1d )
!!!!!!    !
!!!!!!    ! ... Bring the solution to 3D
!!!!!!    !
!!!!!!    CALL planar_average( nnrlocal, naxis, slab_axis, shift, .true., rhozero, rhoioncc1d(cell_min:cell_max) )
!!!!!!    CALL planar_average( nnrlocal, naxis, slab_axis, shift, .true., vzero, vioncc1d(cell_min:cell_max) )
!!!!!!    CALL write_1d( n1d, axis1d, rhoioncc1d, 'rhoioncc1d.dat' )
!!!!!!    CALL write_1d( n1d, axis1d, vioncc1d, 'vioncc1d.dat' )
!!!!!!    !
!!!!!!    ! ... Compute the error on 1D density with 3D switch
!!!!!!    !
!!!!!!    vtotlocal( :, 1 ) = vfix(:) + vzero(:)
!!!!!!    CALL generate_rhoioncc( nnrlocal, 1, switch, vtotlocal, rhooutlocal )
!!!!!!    rhooutlocal = rhooutlocal - rhozero
!!!!!!    charge = SUM(rhooutlocal**2)*domega
!!!!!!    CALL mp_sum( charge, intra_bgrp_comm )
!!!!!!    WRITE(environ_unit,*)'error of 1D charge = ',charge
!!!!!!    ioncc_feval = 1000.D0*charge**2
!!!!!!    !
!!!!!!    ! ... Compute the 3D correction
!!!!!!    !
!!!!!!    vtotlocal = 0.D0
!!!!!!    CALL v_h_of_rho_r( rhoinlocal, ehart, charge, vtotlocal )
!!!!!!    vperiodic = 0.D0
!!!!!!    CALL calc_vperiodic( nnrlocal, nspinlocal, .FALSE., rhoinlocal, vperiodic )
!!!!!!    vtotlocal( :, 1 ) = vtotlocal( :, 1 ) + vperiodic(:) + vfix(:) + vzero(:)
!!!!!!    !
!!!!!!    ! ... Get the new 3D countercharge density
!!!!!!    !
!!!!!!    CALL generate_rhoioncc( nnrlocal, 1, switch, vtotlocal(:,1), rhooutlocal )
!!!!!!    !
!!!!!!    ! ... Integrate the corrective charge: should sum up to zero
!!!!!!    !
!!!!!!    rhooutlocal = (rhooutlocal - rhozero)
!!!!!!    charge = SUM(rhooutlocal)*domega
!!!!!!    CALL mp_sum( charge, intra_bgrp_comm )
!!!!!!    WRITE(environ_unit,*)'error of 3D charge = ',charge
!!!!!!    !
!!!!!!    ioncc_feval = ioncc_feval + charge**2
!!!!!!    !
!!!!!!    DEALLOCATE( vperiodic )
!!!!!!    DEALLOCATE( rhozero )
!!!!!!    DEALLOCATE( vzero )
!!!!!!    DEALLOCATE( rhoioncc1d )
!!!!!!    DEALLOCATE( vioncc1d )
!!!!!!    !
!!!!!!    RETURN
!!!!!!    !
!!!!!!!--------------------------------------------------------------------
!!!!!! END FUNCTION ioncc_feval
!!!!!!!--------------------------------------------------------------------
!!!!!!!!--------------------------------------------------------------------
!!!!!!!  SUBROUTINE calc_rhoioncc( nnr, nspin, rhoioncc, vioncc )
!!!!!!!!--------------------------------------------------------------------
!!!!!!!    !
!!!!!!!    IMPLICIT NONE
!!!!!!!    !
!!!!!!!    INTEGER, PARAMETER :: maxiter = 100000
!!!!!!!    REAL( DP ), PARAMETER :: mix = 0.1
!!!!!!!    REAL( DP ), PARAMETER :: tol = 1.D-6
!!!!!!!    !
!!!!!!!    INTEGER, INTENT(IN) :: nnr, nspin
!!!!!!!    !
!!!!!!!    REAL( DP ), INTENT(INOUT) :: rhoioncc(nnr)
!!!!!!!    REAL( DP ), INTENT(OUT) :: vioncc(nnr)
!!!!!!!    !
!!!!!!!    INTEGER :: iter
!!!!!!!    REAL( DP ) :: charge, deltarho
!!!!!!!    REAL( DP ) :: ax, bx, cx, fa, fb, fc, xout
!!!!!!!    REAL( DP ), ALLOCATABLE, DIMENSION(:) :: residue
!!!!!!!    !
!!!!!!!    REAL( DP ) :: golden
!!!!!!!    !
!!!!!!!    CALL setlocal( nnr, nspin, rhoioncc )
!!!!!!!    !
!!!!!!!    ax = xstern - 1.0D0
!!!!!!!    bx = xstern + deltastern + 0.5D0
!!!!!!!    !
!!!!!!!    fa = ioncc_feval(ax)
!!!!!!!    WRITE(environ_unit,*)'ax = ',ax,' fa = ',fa
!!!!!!!    fb = ioncc_feval(bx)
!!!!!!!    WRITE(environ_unit,*)'bx = ',bx,' fb = ',fb
!!!!!!!    CALL mnbrak( ax, bx, cx, fa, fb, fc, ioncc_feval )
!!!!!!!    WRITE(environ_unit,*)'ax = ',ax,' fa = ',fa
!!!!!!!    WRITE(environ_unit,*)'bx = ',bx,' fb = ',fb
!!!!!!!    WRITE(environ_unit,*)'cx = ',cx,' fc = ',fc
!!!!!!!    !
!!!!!!!    ! ... 3D Refinement, iterating on the difference of the 1D solution and the 3D solution
!!!!!!!    !     assume there is no difference (wishfull tihnking)
!!!!!!!    !
!!!!!!!    rhoioncc = 0.D0
!!!!!!!    !
!!!!!!!    ALLOCATE(residue(nnr))
!!!!!!!!    DO iter = 1, 1!maxiter
!!!!!!!      !
!!!!!!!      rhoinlocal = rhoioncc
!!!!!!!      !
!!!!!!!      charge = golden(ax,bx,cx,ioncc_feval,tol,xout)
!!!!!!!      WRITE(environ_unit,*)'xout = ',xout,' charge = ',charge
!!!!!!!!!!      !
!!!!!!!!!!      residue = rhooutlocal - rhoinlocal
!!!!!!!!!!      deltarho = SUM(residue(:)*residue(:))
!!!!!!!!!!      CALL mp_sum( deltarho, intra_bgrp_comm )
!!!!!!!!!!      deltarho = SQRT(deltarho)/DBLE(ntot)
!!!!!!!!!!      IF ( verbose .GE. 1 ) WRITE(environ_unit,8101)iter,deltarho,charge,charge_fix
!!!!!!!!!!      !
!!!!!!!!!!!      CALL mix_rhoioncc_gspace( nnr, 1, mix, rhoout, rhoin, rhoioncc )
!!!!!!!!!!      rhoioncc = mix * rhooutlocal + ( 1.D0 - mix ) * rhoinlocal
!!!!!!!!!!      !
!!!!!!!!!!      IF (deltarho.LT.tol) THEN
!!!!!!!!!!        IF ( verbose .GE. 1 ) WRITE(environ_unit,8102)
!!!!!!!!!!        EXIT
!!!!!!!!!!      ELSE IF ( iter .EQ. maxiter ) THEN
!!!!!!!!!!        WRITE(stdout,8103)
!!!!!!!!!!      ENDIF
!!!!!!!!!!      !
!!!!!!!!!!    ENDDO
!!!!!!!!!!    WRITE(stdout,8104)deltarho,iter
!!!!!!!!!!    DEALLOCATE(residue)
!!!!!!!!!!    !
!!!!!!!!!!    vioncc(:) = vtotlocal(:,1) - vfix(:)
!!!!!!!!!!    !
!!!!!!!!!!    CALL cleanlocal()
!!!!!!!    !
!!!!!!!    RETURN
!!!!!!!    !
!!!!!!!8101 FORMAT(1X,'iter = ',I5,' deltarho = ',E14.6,' charge_out = ',E14.6,' charge_in = ',E14.6)
!!!!!!!8102 FORMAT(1X,'Level 3 ioncc charges are converged!')
!!!!!!!8103 FORMAT(1X,'Warning: Level 3 ioncc charges are not converged')
!!!!!!!8104 FORMAT(1X,'        L3 ioncc accuracy =',1PE8.1,', # of iterations = ',i5)
!!!!!!!8105 FORMAT(1X,'iter = ',I5,' charge_tmp = ',E14.6,' charge_ext = ',E14.6)
!!!!!!!!--------------------------------------------------------------------
!!!!!!!  END SUBROUTINE calc_rhoioncc
!!!!!!!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE generate_rhoioncc( n, type, switch, v, rho )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    REAL, PARAMETER :: max_arg = 3.D1 !1.D2
!    !
!    INTEGER, INTENT(IN) :: n, type
!    REAL( DP ), INTENT(IN) :: switch(n)
!    REAL( DP ), INTENT(IN) :: v(n)
!    REAL( DP ), INTENT(INOUT) :: rho(n)
!    !
!    INTEGER :: i
!    REAL( DP ) :: arg, fact, denom
!    !
!    fact = -2.D0 * zion * cion
!    DO i = 1, n
!      IF (ABS(switch(i)).LT.1.D-30) THEN
!        rho(i) = 0.D0
!      ELSE
!        arg = zion * invkbt * v(i)
!        IF ( ABS(arg) .GT. max_arg ) arg = arg/ABS(arg)*max_arg ! try to avoid exploding charge
!        IF ( type .EQ. 1 ) THEN
!          IF ( cionmax .EQ. 0.D0 ) THEN
!            rho(i) = fact * sinh( arg ) * switch(i)
!          ELSE
!            IF ( ABS(arg) .LT. max_arg ) THEN
!              denom = 1.D0 + cion / cionmax * 2.D0 * ( cosh(arg) - 1.D0 )
!              rho(i) = fact * sinh( arg ) / denom * switch(i)
!            ELSE
!              rho(i) = - cionmax * zion * switch(i) * arg/ABS(arg)
!            ENDIF
!          ENDIF
!        ELSE IF ( type .EQ. 2 ) THEN
!          rho(i) = fact * ( arg ) * switch(i)
!        END IF
!      ENDIF
!    ENDDO
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE generate_rhoioncc
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE generate_drhoioncc( n, type, switch, v, drho )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    REAL, PARAMETER :: max_arg = 2.D1
!    !
!    INTEGER, INTENT(IN) :: n, type
!    REAL( DP ), INTENT(IN) :: switch(n)
!    REAL( DP ), INTENT(IN) :: v(n)
!    REAL( DP ), INTENT(INOUT) :: drho(n)
!    !
!    INTEGER :: i
!    REAL( DP ) :: arg, fact, denom, numer, alpha
!    !
!    drho = 0.D0
!    fact = - 2.D0 * zion * cion * zion * invkbt
!    IF ( cionmax .NE. 0.D0 ) alpha = 2.D0 * cion / cionmax
!    DO i = 1, n
!      IF (ABS(switch(i)).LT.1.D-30) THEN
!        drho(i) = 0.D0
!      ELSE
!        arg = zion * invkbt * v(i)
!        IF ( ABS(arg) .GT. max_arg ) arg = arg/ABS(arg)*max_arg ! try to avoid exploding charge
!        IF ( type .EQ. 1 ) THEN
!          IF ( cionmax .EQ. 0.D0 ) THEN
!            drho(i) = fact * cosh( arg ) * switch(i)
!          ELSE
!            IF ( ABS(arg) .LT. max_arg ) THEN
!              denom = 1.D0 + alpha * ( cosh(arg) - 1.D0 )
!              numer = alpha + ( 1.D0 - alpha ) * cosh(arg)
!              drho(i) = fact * switch(i) * numer / denom ** 2
!            ELSE
!              drho(i) = 0.D0
!            ENDIF
!          ENDIF
!        ELSE IF ( type .EQ. 2 ) THEN
!          drho(i) = fact * switch(i)
!        END IF
!      ENDIF
!    ENDDO
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE generate_drhoioncc
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE generate_d2rhoioncc( n, type, switch, v, d2rho )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    REAL, PARAMETER :: max_arg = 2.D1
!    !
!    INTEGER, INTENT(IN) :: n, type
!    REAL( DP ), INTENT(IN) :: switch(n)
!    REAL( DP ), INTENT(IN) :: v(n)
!    REAL( DP ), INTENT(INOUT) :: d2rho(n)
!    !
!    INTEGER :: i
!    REAL( DP ) :: arg, fact, denom, numer, alpha
!    !
!    d2rho = 0.D0
!    fact = - 2.D0 * zion * cion * zion * invkbt * zion * invkbt
!    IF ( cionmax .NE. 0.D0 ) alpha = 2.D0 * cion / cionmax
!    DO i = 1, n
!      IF (ABS(switch(i)).LT.1.D-30) THEN
!        d2rho(i) = 0.D0
!      ELSE
!        arg = zion * invkbt * v(i)
!        IF ( ABS(arg) .GT. max_arg ) arg = arg/ABS(arg)*max_arg ! try to avoid exploding charge
!        IF ( type .EQ. 1 ) THEN
!          IF ( cionmax .EQ. 0.D0 ) THEN
!            d2rho(i) = fact * sinh( arg ) * switch(i)
!          ELSE
!            IF ( ABS(arg) .LT. max_arg ) THEN
!              denom = 1.D0 + alpha * ( cosh(arg) - 1.D0 )
!              numer = 1.D0 - alpha * ( alpha + 2.D0 ) + alpha * ( alpha - 1.D0 ) * cosh(arg)
!              d2rho(i) = fact * switch(i) * sinh(arg) * numer / ( denom )**3
!            ELSE
!              d2rho(i) = 0.D0
!            ENDIF
!          ENDIF
!        ELSE IF ( type .EQ. 2 ) THEN
!          d2rho(i) = 0.D0
!        END IF
!      ENDIF
!    ENDDO
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE generate_d2rhoioncc
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE generate_cioncc( n, nsp, type, switch, v, c )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    REAL, PARAMETER :: max_arg = 3.D1 !1.D2
!    !
!    INTEGER, INTENT(IN) :: n, nsp, type
!    REAL( DP ), INTENT(IN) :: switch(n), v(n)
!    REAL( DP ), INTENT(INOUT) :: c(nsp,n)
!    !
!    INTEGER :: i, isp
!    REAL( DP ) :: arg, fact, denom
!    !
!    ! For the moment, assume nsp = 2, isp = 1 is +zion, isp = 2 is -zion
!    ! this parameters will be among the input keywords
!    ! c in output is the ioncc concentration divided by cionmax(r)=cionmax*switch(r)
!    !
!    c = 0.D0
!    DO isp = 1, nsp
!       fact = cb(isp) / cionmax
!       DO i = 1, n
!          IF (ABS(switch(i)).LT.1.D-30) THEN
!             c(isp,i) = 0.D0
!          ELSE
!             arg = - z(isp) * invkbt * v(i)
!             IF ( arg .LT. -max_arg ) THEN
!                c(isp,i) = 0.D0 ! - fact
!                GOTO 20
!             ELSE IF ( arg .GT. max_arg ) THEN
!                arg = max_arg           ! try to avoid exploding charge
!             ENDIF
!             IF ( type .EQ. 1 ) THEN
!                IF ( cionmax .EQ. 0.D0 ) THEN
!                   c(isp,i) = fact *  exp(arg) ! ( exp( arg ) - 1.D0 )
!                ELSE
!                   IF ( arg .LT. max_arg ) THEN
!                      denom = 1.D0 + cb(isp) / cionmax * 2.D0 * ( cosh(arg) - 1.D0 )
!                      c(isp,i) = fact * exp( arg ) / denom ! ( exp( arg ) / denom - 1.D0 )
!                   ELSE
!                      c(isp,i) = 1.D0 ! - fact
!                   ENDIF
!                ENDIF
!             ELSE IF ( type .EQ. 2 ) THEN
!                c(isp,i) = fact * ( arg + 1.D0 ) ! arg
!             END IF
!          ENDIF
!20        CONTINUE
!       ENDDO
!    ENDDO
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE generate_cioncc
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE generate_dcioncc( n, nsp, type, switch, v, dc )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    REAL, PARAMETER :: max_arg = 3.D1 !1.D2
!    !
!    INTEGER, INTENT(IN) :: n, nsp, type
!    REAL( DP ), INTENT(IN) :: switch(n), v(n)
!    REAL( DP ), INTENT(INOUT) :: dc(nsp,n)
!    !
!    INTEGER :: i, isp
!    REAL( DP ) :: arg, fact, denom, numer, alpha
!    !
!    ! For the moment, assume nsp = 2, isp = 1 is +zion, isp = 2 is -zion
!    ! this parameters will be among the input keywords
!    ! c in output is the ioncc concentration divided by cionmax(r)=cionmax*switch(r)
!    !
!    dc = 0.D0
!    DO isp = 1, nsp
!       fact = - cb(isp) / cionmax * z(isp) * invkbt
!       IF ( cionmax .NE. 0.D0 ) alpha = 2.D0 * cb(isp) / cionmax
!       DO i = 1, n
!          IF (ABS(switch(i)).LT.1.D-30) THEN
!             dc(isp,i) = 0.D0
!          ELSE
!             arg = - z(isp) * invkbt * v(i)
!             IF ( arg .LT. -max_arg ) THEN
!                dc(isp,i) = 0.D0
!                GOTO 20
!             ELSE IF ( arg .GT. max_arg ) THEN
!                arg = max_arg           ! try to avoid exploding charge
!             ENDIF
!             IF ( type .EQ. 1 ) THEN
!                IF ( cionmax .EQ. 0.D0 ) THEN
!                   dc(isp,i) = fact *  ( exp( arg ) )
!                ELSE
!                   IF ( arg .LT. max_arg ) THEN
!                      denom = 1.D0 + alpha * ( cosh(arg) - 1.D0 )
!                      numer = alpha + ( 1.D0 - alpha ) * exp( arg )
!                      dc(isp,i) = fact * numer / denom ** 2
!                   ELSE
!                      dc(isp,i) = 0.D0
!                   ENDIF
!                ENDIF
!             ELSE IF ( type .EQ. 2 ) THEN
!                dc(isp,i) = fact
!             END IF
!          ENDIF
!20        CONTINUE
!       ENDDO
!    ENDDO
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE generate_dcioncc
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE calc_sioncc(  n, nsp, switch, c, s )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: n, nsp
!    REAL( DP ), INTENT(IN) :: switch(n)
!    REAL( DP ), INTENT(IN) :: c(nsp,n)
!    REAL( DP ), INTENT(OUT) :: s(n)
!    !
!    INTEGER :: i, isp
!    REAL( DP ) :: ctmp
!    !
!    s = 0.D0
!    DO i = 1, n
!       IF ( switch(i) .LT. 1.D-30 ) CYCLE
!       DO isp = 1, nsp
!          IF ( c(isp,i) .LE. 0.D0 ) CYCLE
!          s(i) = s(i) + c(isp,i)*LOG(c(isp,i))
!!          IF ( i .EQ. 1 ) write(environ_unit,*)'s(1) ',s(i),c(isp,1),LOG(c(isp,1))
!       ENDDO
!       ctmp = 1.D0 - SUM( c(:,i) )
!       IF ( ctmp .LE. 0.D0 ) CYCLE
!       s(i) = s(i) + ctmp * LOG( ctmp )
!!       IF ( i .EQ. 1 ) write(environ_unit,*)'s(1) ',s(i),ctmp,LOG(ctmp)
!    ENDDO
!    s = - kbt * cionmax * switch * s
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE calc_sioncc
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE calc_dsioncc(  n, nsp, switch, c, ds )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: n, nsp
!    REAL( DP ), INTENT(IN) :: switch(n)
!    REAL( DP ), INTENT(IN) :: c(nsp,n)
!    REAL( DP ), INTENT(OUT) :: ds(nsp,n)
!    !
!    INTEGER :: i, isp
!    REAL( DP ) :: ctmp
!    !
!    ds = 0.D0
!    DO i = 1, n
!       IF ( switch(i) .LT. 1.D-30 ) CYCLE
!       DO isp = 1, nsp
!          IF ( c(isp,i) .LE. 0.D0 ) CYCLE
!          ds(isp,i) = LOG( c(isp,i) ) + 1.D0
!       ENDDO
!       ctmp = 1.D0 - SUM( c(:,i) )
!       IF ( ctmp .LE. 0.D0 ) CYCLE
!       DO isp = 1, nsp
!          ds(isp,i) = ds(isp,i) - ( LOG( ctmp ) + 1.D0 )
!       ENDDO
!    ENDDO
!    DO isp = 1, nsp
!       ds(isp,:) = - kbt * cionmax * switch(:) * ds(isp,:)
!    ENDDO
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE calc_dsioncc
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE sum_cioncc(  n, nsp, switch, weight, c, csum )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: n, nsp
!    REAL( DP ), INTENT(IN) :: switch(n)
!    REAL( DP ), INTENT(IN) :: weight(nsp)
!    REAL( DP ), INTENT(IN) :: c(nsp,n)
!    REAL( DP ), INTENT(OUT) :: csum(n)
!    !
!    INTEGER :: isp
!    !
!    csum = 0.D0
!    DO isp = 1, nsp
!       csum(:) = csum(:) + weight(isp) * c(isp,:) * cionmax * switch(:)
!    ENDDO
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE sum_cioncc
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE calc_muioncc(  nsp, cb, mu )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: nsp
!    REAL( DP ), INTENT(IN) :: cb(nsp)
!    REAL( DP ), INTENT(OUT) :: mu(nsp)
!    !
!    INTEGER :: isp
!    REAL( DP ) :: ctot
!    !
!    ctot = SUM(cb(:))
!    IF ( cionmax .LT. ctot ) THEN
!       WRITE( stdout, * )&
!         &'ERROR: cionmax is lower that total bulk concentration =',ctot
!       STOP
!    ENDIF
!    DO isp = 1, nsp
!      mu(isp) = - kbt * LOG( (cionmax - ctot)/cb(isp) )
!    ENDDO
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE calc_muioncc
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE calc_eioncc(  nnr, rho, eioncc )
!!--------------------------------------------------------------------
!    !
!    USE environ_base,  ONLY : vioncc
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: nnr
!    REAL( DP ), INTENT(IN) :: rho(nnr)
!    REAL( DP ), INTENT(OUT) :: eioncc
!    !
!    REAL( DP ), ALLOCATABLE :: rhotot(:)
!    !
!    ALLOCATE(rhotot(nnr))
!    rhotot = rhoions + rho
!    !
!    eioncc = 0.5D0 * SUM( vioncc(:) * rhotot( : ) ) * domega
!    !
!    DEALLOCATE(rhotot)
!    !
!    CALL mp_sum( eioncc, intra_bgrp_comm )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE calc_eioncc
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!      SUBROUTINE write_1d( n1d, axis1d, f1d, filename )
!!--------------------------------------------------------------------
!      !
!      USE kinds,          ONLY : DP
!      USE io_global,      ONLY : ionode
!      !
!      IMPLICIT NONE
!      !
!      INTEGER, INTENT(IN)      :: n1d
!      REAL( DP ), INTENT(IN)   :: axis1d( n1d )
!      REAL( DP ), INTENT(IN)   :: f1d( n1d )
!      !
!      CHARACTER(*), INTENT(IN) :: filename
!      !
!      INTEGER                  :: ia
!      !
!      IF( ionode ) THEN
!        !
!        OPEN( 300, file = TRIM( filename ), status = 'unknown' )
!        !
!        DO ia = 1, n1d
!          WRITE( 300, '(f20.10,x,e20.10)') axis1d(ia), f1d(ia)
!        ENDDO
!        !
!        CLOSE( 300 )
!        !
!      END IF
!      !
!      RETURN
!      !
!!--------------------------------------------------------------------
!      END SUBROUTINE write_1d
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!      SUBROUTINE generate_step( n, axis, step )
!!--------------------------------------------------------------------
!      !
!      IMPLICIT NONE
!      !
!      INTEGER, INTENT(IN) :: n
!      !
!      REAL( DP ), INTENT(IN) :: axis(n)
!      REAL( DP ), INTENT(OUT) :: step(n)
!      !
!      INTEGER :: i
!      !
!      DO i = 1, n
!        IF ( axis(i) .GT. 0.D0 ) THEN
!          step(i) = 0.5D0
!        ELSE
!          step(i) = -0.5D0
!        ENDIF
!      ENDDO
!      !
!      RETURN
!      !
!!--------------------------------------------------------------------
!      END SUBROUTINE generate_step
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!      SUBROUTINE generate_switch( n, x0, delta, axis, switch )
!!--------------------------------------------------------------------
!      !
!      IMPLICIT NONE
!      !
!      INTEGER, INTENT(IN) :: n
!      !
!      REAL( DP ), INTENT(IN) :: x0, delta
!      REAL( DP ), INTENT(IN) :: axis(n)
!      REAL( DP ), INTENT(OUT) :: switch(n)
!      !
!      INTEGER :: i
!      !
!      REAL( DP ) :: xmin, xmax, arg
!      !
!      xmin = x0
!      xmax = x0 + delta
!      DO i = 1, n
!        IF ( ABS(axis(i)) .LT. xmin ) THEN
!          switch(i) = 0.D0
!        ELSE IF ( ABS(axis(i)) .LT. xmax ) THEN
!          arg = tpi * ( ABS(axis(i)) - xmin)/(xmax - xmin)
!          switch(i) = ( arg - SIN( arg ) ) / tpi
!        ELSE
!          switch(i) = 1.D0
!        ENDIF
!      ENDDO
!      !
!      RETURN
!      !
!!--------------------------------------------------------------------
!      END SUBROUTINE generate_switch
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!      SUBROUTINE v1d_of_rho1d( n, dx, rho, v )
!!--------------------------------------------------------------------
!      !
!      IMPLICIT NONE
!      !
!      INTEGER, INTENT(IN) :: n
!      !
!      REAL( DP ), INTENT(IN) :: dx
!      REAL( DP ), INTENT(IN) :: rho(n)
!      REAL( DP ), INTENT(OUT) :: v(n)
!      !
!      INTEGER :: i
!      REAL( DP ) :: e(n)
!      REAL( DP ) :: vlast, slope, vshift
!      !
!      CALL integral_function_1d(n,dx,rho,e)
!      CALL integral_function_1d(n,dx,e,v)
!      v = - fpi * e2 * v
!      !
!!      vlast = v(nfdpoint+1)
!!      slope = SUM(rho) * dx * fpi * e2 / 2.D0
!!      DO i = 1, nfdpoint
!!         v(i) = vlast + slope * dx * (i - nfdpoint - 1)
!!      ENDDO
!!      vlast = v(n-nfdpoint)
!!      DO i = n, n-nfdpoint+1
!!         v(i) = vlast + slope * dx * (( n - nfdpoint) - i )
!!      ENDDO
!!      vshift = ( v(1) + v(n) ) / 2.D0
!!      v = v - vshift
!      !
!      RETURN
!      !
!!--------------------------------------------------------------------
!      END SUBROUTINE v1d_of_rho1d
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!      SUBROUTINE integral_function_1d( n, dx, f, fint )
!!--------------------------------------------------------------------
!      !
!      USE kinds,          ONLY : DP
!      !
!      IMPLICIT NONE
!      !
!      INTEGER, INTENT(IN)      :: n
!      REAL( DP ), INTENT(IN)   :: dx
!      REAL( DP ), INTENT(IN)   :: f( n )
!      REAL( DP ), INTENT(OUT)  :: fint( n )
!      !
!      INTEGER                  :: i, j
!      !
!      REAL( DP ) :: c(3)
!      !
!      REAL( DP ) :: ftmp(n)
!      !
!      c(1) = 7.D0/6.D0
!      c(2) = 23.D0/24.D0 - 7.D0/6.D0
!      c(3) = 1.D0 - 23.D0/24.D0
!      !
!      fint = 0.d0
!      !
!      DO i = 2, n
!         fint(i) = fint(i-1)
!         DO j = 1, MIN(i-1,3)
!            fint(i) = fint(i) + c(j) * f(i-j)
!         ENDDO
!      ENDDO
!      !
!      ftmp = 0.D0
!      !
!      DO i = n-1, 1, -1
!         ftmp(i) = ftmp(i+1)
!         DO j = 1, MIN(n-i,3)
!            ftmp(i) = ftmp(i) + c(j) * f(i+j)
!         ENDDO
!      ENDDO
!      !
!      fint = ( fint - ftmp ) / 2.D0 * dx
!      !
!      RETURN
!      !
!!--------------------------------------------------------------------
!      END SUBROUTINE integral_function_1d
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!      SUBROUTINE integral_function_1d_original( n, dx, f, fint )
!!--------------------------------------------------------------------
!      !
!      USE kinds,          ONLY : DP
!      !
!      IMPLICIT NONE
!      !
!      INTEGER, INTENT(IN)      :: n
!      REAL( DP ), INTENT(IN)   :: dx
!      REAL( DP ), INTENT(IN)   :: f( n )
!      REAL( DP ), INTENT(OUT)  :: fint( n )
!      !
!      INTEGER                  :: i
!      !
!      REAL( DP ) :: c1, c2, c3, c4, shift
!      !
!      REAL( DP ) :: ftmp(n)
!      !
!      c1 = 3.D0/8.D0
!      c2 = 7.D0/6.D0 - 3.D0/8.D0
!      c3 = 23.D0/24.D0 - 7.D0/6.D0
!      c4 = 1.D0 - 23.D0/24.D0
!      !
!      fint = 0.d0
!      !
!      fint(2) = 0.5*(f(1)+f(2))
!      !
!      fint(3) = (f(1)+f(3)+4.D0*f(2))/3.D0
!      !
!      fint(4) = (f(1)+f(4)+3.D0*(f(2)+f(3)))*3.D0/8.D0
!      !
!      fint(5) = (14.D0*(f(1)+f(5))+64.D0*(f(2)+f(4))+24*f(3))/45.D0
!      !
!      fint(6) = 3.D0/8.D0*(f(1)+f(6))+7.D0/6.D0*(f(2)+f(5))+23.D0/24.D0*(f(3)+f(4))
!      DO i = 7, n
!        fint(i) = fint(i-1) + f(i)*c1 + f(i-1)*c2 + f(i-2)*c3 + f(i-3)*c4
!      ENDDO
!      !
!      ftmp = 0.D0
!      !
!      ftmp(n-1) = 0.5*(f(n)+f(n-1))
!      !
!      ftmp(n-2) = (f(n)+f(n-2)+4.D0*f(n-1))/3.D0
!      !
!      ftmp(n-3) = (f(n)+f(n-3)+3.D0*(f(n-1)+f(n-2)))*3.D0/8.D0
!      !
!      ftmp(n-4) = (14.D0*(f(n)+f(n-4))+64.D0*(f(n-1)+f(n-3))+24*f(n-2))/45.D0
!      !
!      ftmp(n-5) = 3.D0/8.D0*(f(n)+f(n-5))+7.D0/6.D0*(f(n-1)+f(n-4))+23.D0/24.D0*(f(n-2)+f(n-3))
!      DO i = n-6, 1, -1
!        ftmp(i) = ftmp(i+1) + f(i)*c1 + f(i+1)*c2 + f(i+2)*c3 + f(i+3)*c4
!      ENDDO
!      !
!      fint = ( fint - ftmp )/2.D0
!      !
!!      shift = (fint(1) + fint(n))/2.D0
!!      fint = ( fint - shift ) * dx
!      fint = fint * dx
!      !
!      RETURN
!      !
!!--------------------------------------------------------------------
!      END SUBROUTINE integral_function_1d_original
!!--------------------------------------------------------------------
!!!!!!!      SUBROUTINE mix_rhoioncc_gspace( nnr, nspin, mixlocal, rhoout, rhoin, rhoioncc )
!!!!!!!
!!!!!!!      USE fft_base,       ONLY : dfftp
!!!!!!!      USE fft_interfaces, ONLY : fwfft, invfft
!!!!!!!      USE gvect,          ONLY : nl, nlm, ngm, gg, gstart
!!!!!!!      USE control_flags,  ONLY : gamma_only
!!!!!!!
!!!!!!!      IMPLICIT NONE
!!!!!!!
!!!!!!!      INTEGER, INTENT(IN) :: nnr, nspin
!!!!!!!      REAL*8, INTENT(IN) :: mixlocal
!!!!!!!
!!!!!!!      REAL*8, DIMENSION(nnr,nspin), INTENT(IN) :: rhoout, rhoin
!!!!!!!      REAL*8, DIMENSION(nnr,nspin), INTENT(OUT) :: rhoioncc
!!!!!!!      !
!!!!!!!      COMPLEX( DP ), ALLOCATABLE :: rhooutg(:,:), rhoing(:,:), rhoionccg(:,:)
!!!!!!!      COMPLEX( DP ), ALLOCATABLE :: aux( : )
!!!!!!!      INTEGER :: is, ig
!!!!!!!      REAL*8 :: fac
!!!!!!!
!!!!!!!      ALLOCATE( rhooutg(ngm,nspin), rhoing(ngm,nspin), rhoionccg(ngm,nspin))
!!!!!!!      ALLOCATE( aux(nnr) )
!!!!!!!
!!!!!!!      DO is = 1, nspin
!!!!!!!        aux(:) = CMPLX(rhoout( : , is ),0.D0,kind=dp)
!!!!!!!        CALL fwfft ('Dense', aux, dfftp)
!!!!!!!        rhooutg(:,is) = aux(nl(:))
!!!!!!!        aux(:) = CMPLX(rhoin( : , is ),0.D0,kind=dp)
!!!!!!!        CALL fwfft ('Dense', aux, dfftp)
!!!!!!!        rhoing(:,is) = aux(nl(:))
!!!!!!!      END DO
!!!!!!!
!!!!!!!      rhoionccg = 0.D0
!!!!!!!      DO ig = gstart, ngm
!!!!!!!        fac = mixlocal / (1.D0+gg(ig))
!!!!!!!        rhoionccg(ig,1) = fac * rhooutg(ig,1) + ( 1.D0 - fac ) * rhoing(ig,1)
!!!!!!!        IF ( nspin == 2 ) &
!!!!!!!          & rhoionccg(ig,2) = fac * rhooutg(ig,2) + ( 1.D0 - fac ) * rhoing(ig,2)
!!!!!!!      ENDDO
!!!!!!!
!!!!!!!      DO is = 1, nspin
!!!!!!!        aux(nl(1:ngm)) = rhoionccg(1:ngm,is)
!!!!!!!        IF ( gamma_only ) &
!!!!!!!          & aux(nlm(1:ngm)) = CMPLX(REAL(rhoionccg(1:ngm,is)),-AIMAG(rhoionccg(1:ngm,is)),KIND=dp)
!!!!!!!        CALL invfft ('Dense', aux, dfftp)
!!!!!!!        rhoioncc(:,is) = DBLE(aux(:))
!!!!!!!      ENDDO
!!!!!!!
!!!!!!!      DEALLOCATE( aux )
!!!!!!!      DEALLOCATE( rhooutg, rhoing, rhoionccg )
!!!!!!!
!!!!!!!      RETURN
!!!!!!!
!!!!!!!      END SUBROUTINE mix_rhoioncc_gspace
!!!!!!!
!!--------------------------------------------------------------------
!  FUNCTION ioncc_charge(lambda)
!!--------------------------------------------------------------------
!    !
!    USE kinds, ONLY: DP
!    !
!    IMPLICIT NONE
!    !
!    REAL( DP ), INTENT(IN) :: lambda
!    REAL( DP ) :: ioncc_charge
!    !
!    ! ... Local variables
!    !
!    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: rholocal, vlocal
!    REAL( DP ), DIMENSION( :, : ), ALLOCATABLE :: clocal
!    !
!    ! ... Allocate
!    !
!    ALLOCATE( clocal( nsp, n1d ) )
!    ALLOCATE( rholocal( n1d ) )
!    ALLOCATE( vlocal( n1d ) )
!    vlocal = vzero1d - lambda
!    !
!    CALL generate_cioncc( n1d, nsp, 1, switch1d, vlocal, clocal )
!    !
!    CALL sum_cioncc( n1d, nsp, switch1d, z, clocal, rholocal )
!    !
!    ioncc_charge = SUM( rholocal ) * darea + charge_fix
!    !
!    DEALLOCATE( vlocal )
!    DEALLOCATE( clocal )
!    DEALLOCATE( rholocal )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
! END FUNCTION ioncc_charge
!!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE electrolyte_utils
!--------------------------------------------------------------------
