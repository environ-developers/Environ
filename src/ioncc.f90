!--------------------------------------------------------------------
MODULE ioncc
!--------------------------------------------------------------------

  USE kinds,          ONLY: DP
  USE constants,      ONLY: k_boltzmann_ry, pi, tpi, fpi, e2
  USE io_global,      ONLY: stdout
  USE mp,             ONLY: mp_sum
  USE mp_bands,       ONLY: intra_bgrp_comm
  USE environ_cell,   ONLY: domega, alat, omega, ntot
  USE environ_ions,   ONLY: avg_pos, rhoions
  USE environ_base,   ONLY: verbose, ir_end, environ_unit,               & 
                            env_ioncc_level, nrep, cion, zion,           &
                            solvent_temperature, rhoioncc, rhopolcc,     & 
                            env_static_permittivity, rhopol, slab_axis,  &
                            rhomin, rhopb
  USE generate_function, ONLY: planar_average
  USE periodic,       ONLY: calc_vperiodic
  USE environ_debug,  ONLY: write_cube
  USE generate_f_of_rho, ONLY: epsilonfunct

  IMPLICIT NONE

  INTEGER :: naxis, shift, cell_max, cell_min, n1d
  REAL( DP ) :: area, darea, axis_length, dx, xstern, deltastern, vavg

  REAL( DP ) :: kbt, invkbt, fsw
  
  REAL( DP ), ALLOCATABLE :: axis(:), switch(:), step(:)
  REAL( DP ), ALLOCATABLE :: axis1d(:), switch1d(:), step1d(:)
  REAL( DP ), ALLOCATABLE :: rhoioncc1dref(:), rhotot1dref(:)

  SAVE

  PRIVATE

  PUBLIC :: ioncc_initbase, ioncc_initcell, ioncc_clean, calc_vioncc, &
            calc_eioncc

CONTAINS

!--------------------------------------------------------------------
  SUBROUTINE ioncc_initbase( nnr )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nnr
    !
    IF (ALLOCATED(axis)) DEALLOCATE(axis)
    ALLOCATE(axis(nnr))
    axis = 0.D0
    !    
    IF (ALLOCATED(switch)) DEALLOCATE(switch)
    ALLOCATE(switch(nnr))
    switch = 0.D0
    !
    IF (ALLOCATED(step)) DEALLOCATE(step)
    ALLOCATE(step(nnr))
    step = 0.D0
    !
    kbt = k_boltzmann_ry * solvent_temperature
    invkbt = 1.D0/kbt
    !
    fsw = LOG( rhomin / rhopb )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE ioncc_initbase
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE ioncc_initcell( nnr, n1, n2, n3, at )
!--------------------------------------------------------------------
    !
    USE generate_function, ONLY: generate_axis
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nnr, n1, n2, n3
    REAL(DP), INTENT(IN) :: at(3,3)
    !
    INTEGER :: ia, irep, imin, imax
    REAL(DP) :: xmin, xmax, arg
    !
    avg_pos = 0.D0
    !
    ! ... Set the properties of the axis and generate it
    !
    axis_length = at( slab_axis, slab_axis ) * alat
    SELECT CASE (slab_axis)
    CASE (1)
      naxis = n1
    CASE (2) 
      naxis = n2
    CASE (3)
      naxis = n3
    END SELECT
    CALL generate_axis( nnr, slab_axis, avg_pos, axis )
    IF ( verbose .GE. 2 ) CALL write_cube( nnr, axis, 'axis.cube' )
    !
    ! ... Set the main indexes to convert 3D quantities in the cell
    !     into 1D quantities in extended cell (unit cell plus 2*nrep replicas)
    !
    n1d = (2*nrep+1)*naxis
    cell_min = nrep*naxis + 1
    cell_max = (nrep+1)*naxis
    dx = axis_length/DBLE(naxis)
    area = omega/axis_length
    darea = dx*area
    shift = naxis/2 - NINT(avg_pos(slab_axis)*alat/dx) 
    IF ( verbose .GE. 1 ) THEN 
      WRITE(environ_unit,8301)naxis,naxis/2,nrep,n1d
      WRITE(environ_unit,8302)cell_min,cell_max,dx,area
      WRITE(environ_unit,8303)darea,avg_pos(slab_axis)*alat,shift
    ENDIF
    !
    ! ... Properties of the ioncc density (only temporary)
    !
    xstern = 3.D0 ! only temporary, PUT THIS IN INPUT KEYWORDS
    deltastern = 4.D0 ! only temporary, PUT THIS IN INPUT KEYWORDS
    !
    ! ... Generate the 1D axis: first compute it in the original cell
    !
    IF (ALLOCATED(axis1d)) DEALLOCATE(axis1d)
    ALLOCATE(axis1d(n1d))
    axis1d = 0.D0
    CALL planar_average( nnr, naxis, slab_axis, shift, .false., axis, axis1d(cell_min:cell_max) ) 
    !
    ! ... then generate the extended replicas
    !
    DO irep = 1, 2*nrep+1
      imin = (irep-1)*naxis + 1
      imax = imin + naxis - 1
      axis1d(imin:imax) = axis1d(cell_min:cell_max) + ( irep - 2*nrep ) * axis_length
    ENDDO  
    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, axis1d, 'axis1d.dat' )
    !
    ! ... Generate the step function, to allow finite step in potential
    !
    CALL generate_step(nnr,axis,step)
    IF (ALLOCATED(step1d)) DEALLOCATE(step1d)
    ALLOCATE(step1d(n1d))
    CALL generate_step(n1d,axis1d,step1d)
    !
    ! ... Generate the 1D switching function, for Level 2 or higher ioncc density
    !
    IF (ALLOCATED(switch1d)) DEALLOCATE(switch1d)
    ALLOCATE(switch1d(n1d))
    xmax = xstern + deltastern
    CALL generate_switch( n1d, 2, xstern, xmax, 1.D0, axis1d, switch1d )
    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, switch1d, 'switch1d.dat' )    
    !
    RETURN
    !
8301 FORMAT(1X,'naxis = ',i4,' naxis/2 = ',i4,' nrep = ',i3,' n1d = ',i4)
8302 FORMAT(1X,'cell_min = ',i4,' cell_max = ',i4,' dx = ',F14.6,' area = ',F14.6)
8303 FORMAT(1X,' darea = ',F16.10,' avg_pos = ',F14.6,' shift = ',i5)
!--------------------------------------------------------------------
  END SUBROUTINE ioncc_initcell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE ioncc_clean( )
!--------------------------------------------------------------------
    !
    ! ... Clean up of local variables
    !
    IMPLICIT NONE
    !
    IF (ALLOCATED(axis)) DEALLOCATE(axis)
    IF (ALLOCATED(switch)) DEALLOCATE(switch)
    IF (ALLOCATED(step)) DEALLOCATE(step)
    IF (ALLOCATED(axis1d)) DEALLOCATE(axis1d)
    IF (ALLOCATED(switch1d)) DEALLOCATE(switch1d)
    IF (ALLOCATED(step1d)) DEALLOCATE(step1d)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE ioncc_clean
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_vioncc(  nnr, nspin, rhoelec, vioncc )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nnr, nspin
    REAL( DP ), INTENT(IN) :: rhoelec(nnr)
    REAL( DP ), INTENT(OUT) :: vioncc(nnr)
    !
    INTEGER :: ia, ir, irep, iter
    REAL( DP ) :: tot_charge, tot_dipole
    REAL( DP ), ALLOCATABLE :: rhotot(:)
    REAL( DP ), ALLOCATABLE :: vzero1d(:), vioncc1d(:), switch1dtmp(:)
    REAL( DP ), ALLOCATABLE :: rhostern1d(:), vstern1d(:)
    REAL( DP ), ALLOCATABLE :: rhoioncc1d(:), rhotot1d(:)
    !
    REAL( DP ) :: xmax, xmin, tmp, arg
    !
    CALL start_clock( 'get_ioncc' ) 
    !
    ! ... Initialize the potential
    !
    vioncc = 0.D0
    IF ( env_ioncc_level .LT. 1 ) RETURN
    ALLOCATE( vioncc1d(n1d) ) 
    vioncc1d = 0.D0
    !
    ! ... Set the total density, only exclude the ioncc densities
    !
    ALLOCATE( rhotot(nnr) )
!!!    rhotot = rhoelec + rhoions 
!!!    IF ( env_static_permittivity .GT. 1.D0 ) rhotot = rhotot + rhopol
    !
    ! ... Test desity, we build a uniform charge density thick 2.D0 a.u. around avg_pos
    !     The total chage of the analytic system is 1 e 
    !
    rhotot = 0.D0
    DO ir = 1, ir_end
      IF ( ABS(axis(ir)) .LT. 1.D0 ) rhotot(ir) = -1.D0
    ENDDO
    tot_charge = SUM(rhotot)*domega 
    CALL mp_sum( tot_charge, intra_bgrp_comm ) 
    rhotot = rhotot / ABS(tot_charge) 
    IF ( verbose .GE. 3 ) CALL write_cube( nnr, rhotot, 'rhotot.cube' )
    !
    ! ... Generate the 3D switching function, for Level 3 or higher ioncc density
    !
    IF ( env_ioncc_level .GE. 3 ) THEN
      switch = 0.D0
!!!      DO ir = 1, nnr
!!!        switch( ir ) = ( 1.D0 - switch( ir ) ) *    &
!!!          ( epsilonfunct( rhoelec( ir ), rhomin, rhopb, fsw, 2.D0, 1 ) - 1.D0 )
!!!      END DO
      ! ... test 3D switching function, converted from an homogeneous one in 1D, but shifted wrt level 2 switching function
      xmin = xstern - 1.0 
      xmax = xmin + deltastern
      CALL generate_switch( nnr, 2, xmin, xmax, 1.D0, axis, switch )
      IF ( verbose .GE. 1 ) CALL write_cube( nnr, switch, 'switch.cube' ) ! DEBUG, reset it to 3 or 2
    ENDIF
    !
    ! ... Convert the total source charge into 1D ... 
    !
    ALLOCATE( rhotot1d(n1d) )
    rhotot1d = 0.D0
    CALL planar_average( nnr, naxis, slab_axis, shift, .false., rhotot, rhotot1d(cell_min:cell_max) ) 
    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, rhotot1d, 'rhotot1d.dat' )
    !
    ! ... and compute charge, dipole and non-periodic potential
    !
    tot_charge = SUM(rhotot1d(:))*darea
    tot_dipole = SUM(rhotot1d(:)*axis1d(:))*darea
    ALLOCATE( vzero1d(n1d) ) 
    vzero1d = 0.D0
    CALL v1d_of_rho1d( n1d, dx, rhotot1d, vzero1d )
    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, vzero1d, 'vzero1d.dat' )
    !
    ! ... Compute Level 1 ioncc (analytic 1D diffuse layer with sharp boundary)
    !
    ALLOCATE( rhostern1d(n1d) )
    rhostern1d = 0.D0
    CALL calc_rhostern( n1d, tot_charge, tot_dipole, xstern, axis1d, step1d, vzero1d, vioncc1d, rhostern1d )
    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, rhostern1d, 'rhostern1d.dat' )   
    IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, vioncc1d, 'vstern1d.dat' ) 
    !
    ALLOCATE( rhoioncc1d(n1d) ) 
    rhoioncc1d = 0.D0
    IF ( env_ioncc_level .EQ. 1 ) THEN
      CALL planar_average( nnr, naxis, slab_axis, shift, .true., vioncc, vioncc1d(cell_min:cell_max) )
      CALL planar_average( nnr, naxis, slab_axis, shift, .true., rhoioncc, rhostern1d(cell_min:cell_max) )
    !
    ! ... Compute Level 2 ioncc (self-consistent 1D diffuse layer with switching function)
    !
    ELSE IF ( env_ioncc_level .EQ. 2 ) THEN
      rhoioncc1d = rhostern1d * switch1d ! remove and start from previous guess, unless at the fisrt iteration
      CALL calc_rhoioncc1d( n1d, 0.d0, axis1d, step1d, switch1d, rhotot1d, vioncc1d, rhoioncc1d )
      IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, rhoioncc1d, 'rhoioncc1d.dat' )   
      IF ( verbose .GE. 2 ) CALL write_1d( n1d, axis1d, vioncc1d, 'vioncc1d.dat' ) 
      CALL planar_average( nnr, naxis, slab_axis, shift, .true., vioncc, vioncc1d(cell_min:cell_max) )
      CALL planar_average( nnr, naxis, slab_axis, shift, .true., rhoioncc, rhoioncc1d(cell_min:cell_max) )
    !
    ! ... Compute Level 3 ioncc (self-consistent internal 3D diffuse layer with electrons-based
    !      switching function plus self-consistent external 1D diffuse layer)
    !
    ELSE IF ( env_ioncc_level .EQ. 3 ) THEN  
      rhoioncc1d = rhostern1d
      CALL calc_rhoioncc( nnr, nspin, n1d,  xstern-2, xstern+1, 1.D0, xstern, &
           & axis1d, step1d, rhoioncc1d, switch, rhotot, rhoioncc, vioncc ) 
    ELSE
      WRITE(stdout,*)'ERROR: specified ioncc level is not implemented'
    END IF 
    IF ( verbose .GE. 2 ) CALL write_cube( nnr, vioncc, 'vioncc.cube')
    IF ( verbose .GE. 2 ) CALL write_cube( nnr, rhoioncc, 'rhoioncc.cube')
    !
    CALL stop_clock( 'get_ioncc' ) 
    !
    CALL write_cube( nnr, vioncc(:), 'dump.cube' )
    !
    stop
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_vioncc
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_rhostern( n, charge, dipole, x0, axis, step, v0, v, rho )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n
    !
    REAL( DP ), INTENT(IN) :: charge, dipole
    REAL( DP ), INTENT(IN) :: x0
    REAL( DP ), INTENT(IN) :: axis(n)
    REAL( DP ), INTENT(IN) :: step(n)
    REAL( DP ), INTENT(IN) :: v0(n)
    REAL( DP ), INTENT(OUT) :: v(n)
    REAL( DP ), INTENT(OUT) :: rho(n)
    !
    INTEGER :: i, i1, i2
    REAL( DP ) :: ez, fact, vstern, const
    REAL( DP ) :: v1, v2, dv, vbound
    REAL( DP ) :: arg, asinh, coth, acoth 
    REAL( DP ) :: f1, f2, f3, f4
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
    ! ... Compute the physical properties of the interface
    !
    ez = - tpi * e2 * charge / area / env_static_permittivity 
    fact = - e2 * SQRT( 8.D0 * fpi * cion * kbt / e2 / env_static_permittivity )
    arg = ez/fact 
    asinh = LOG(arg + SQRT( arg**2 + 1 ))
    vstern = 2.D0 * kbt / zion * asinh
    arg = vstern * 0.25D0 * invkbt * zion
    coth = ( EXP( 2.D0 * arg ) + 1.D0 ) / ( EXP( 2.D0 * arg ) - 1.D0 )
    const = coth * EXP( zion * fact * invkbt * 0.5D0 * x0 )
    !
    IF ( verbose .GE. 1 ) WRITE(environ_unit,8001)ez,fact,vstern,const
    !
    ! ... Compute the potential of the source charge at the boundary, to set the constant shift
    !     Linear interpolate between gridpoints around x0 and average over both sides of slab
    !
    i1 = cell_min + naxis/2 + INT( x0 / dx )
    i2 = cell_min + naxis/2 - 1 - INT( x0 / dx ) 
    IF ( i1 .LE. 0 .OR. i1 .GT. n-1 .OR. i2 .LE. 1 .OR. i2 .GT. n ) THEN
      WRITE(stdout,8004)
      STOP
    ENDIF
    v1 = v0(i1)   + ( x0 - axis(i1) )   * ( v0(i1+1) - v0(i1) ) / ( axis(i1+1)-axis(i1) )
    v2 = v0(i2-1) + ( - x0 - axis(i2-1) ) * ( v0(i2) - v0(i2-1) ) / ( axis(i2)-axis(i2-1) )
    dv = 2.D0 * fpi * dipole / area
    vbound = ( v1 + v2 ) * 0.5D0
    !
    IF ( verbose .GE. 2 ) WRITE(environ_unit,8003)i1,axis(i1),axis(i1+1),v0(i1),v0(i1+1),v1
    IF ( verbose .GE. 2 ) WRITE(environ_unit,8003)i2,axis(i2-1),axis(i2),v0(i2-1),v0(i2),v2
    IF ( verbose .GE. 2 ) WRITE(environ_unit,8005)vbound,dv,v1-dv/2.D0,v2+dv/2.D0
    !
    ! ... Compute some constants needed for the calculation
    !
    f1 = - fact * zion * invkbt * 0.5D0 
    f2 = 4.D0 * kbt / zion
    f3 = -2.D0 * zion * cion 
    f4 = zion * invkbt
    !
    IF ( verbose .GE. 2 ) WRITE(environ_unit,8002)f1,f2,f3,f4
    !
    ! ... Compute the analytic potential and charge
    !
    v = 0.D0
    rho = 0.D0
    DO i = 1, n
      IF ( ABS(axis(i)) .LT. x0 ) THEN
        !
        ! ... Constant shift of the potential in the inside
        !
        v(i) = vstern - vbound
        !
      ELSE 
        !
        arg = const * EXP( ABS(axis(i)) * f1)
        IF ( ABS(arg) .GT. 1.D0 ) THEN
          acoth = 0.5D0 * LOG( (arg + 1.D0) / (arg - 1.D0) )
        ELSE 
          acoth = 0.D0
        END IF 
        v(i) =  f2 * acoth
        rho(i) = f3 * SINH( f4 * v(i) )
        !
        ! ... Need to remove source potential from the outside
        !
        v(i) = v(i) + dv*step(i)  - v0(i)
        !
      ENDIF
    ENDDO
    !
    RETURN
    !
8001 FORMAT(1X,'ez = ',E14.6,' fact =',E14.6,' v0 = ',E14.6,' const = ',E14.6)
8002 FORMAT(1X,'f1 = ',E14.6,' f2 = ',E14.6,' f3 = ',E14.6,' f4 = ',E14.6)
8003 FORMAT(1X,'i = ',i4,' x1 = ',F14.6,' x2 = ',F14.6,' v1 = ',F14.6,' v2 = ',F14.6,' v = ',F14.6)
8004 FORMAT(1X,'ERROR: stern layer is outside of 1d cell, need to use more nrep')
8005 FORMAT(1X,'vbound = ',F14.6,' dv = ',F14.6,' v1_shifted = ',F14.6,' v2_shifted = ',F14.6) 
!--------------------------------------------------------------------
  END SUBROUTINE calc_rhostern
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_rhoioncc1d(  n, v0, axis, step, switch, rhofix, v, rho )
!--------------------------------------------------------------------
    IMPLICIT NONE
    !
    INTEGER, PARAMETER :: maxiter = 1000000
    REAL( DP ), PARAMETER :: mix = 0.1
    REAL( DP ), PARAMETER :: tol = 1.D-10
    !
    INTEGER, INTENT(IN) :: n
    !
    REAL( DP ), INTENT(IN) :: v0, switch(n), step(n), axis(n)
    REAL( DP ), INTENT(IN) :: rhofix(n)
    !
    REAL( DP ), INTENT(OUT) :: v(n)
    REAL( DP ), INTENT(INOUT) :: rho(n)
    !
    INTEGER :: iter, i
    REAL( DP ) :: charge_in, charge_out, dipole_in, v1d0, dv, deltarho, mixlocal
    !
    REAL( DP ) :: rhozero(n), rhotot(n)
    REAL( DP ) :: vfix(n)
    !
    REAL( DP ) :: rhoin(n), rhoout(n), residue(n)
    !
    ! ... Compute total charge and potential of the source density
    !
    charge_in = SUM(rhofix)*darea
    dipole_in = SUM(rhofix(:)*axis(:))*darea
    CALL v1d_of_rho1d( n, dx, rhofix, vfix )
    v1d0 = SUM(vfix(cell_min:cell_max))/DBLE(naxis)
    IF ( verbose .GE. 2 ) write(environ_unit,*)' v1d0 ',v1d0,' shift ',v0-v1d0
    IF ( verbose .GE. 3 ) CALL write_1d( n, axis, vfix, 'vfix1d.dat' )
    dv = 2.D0 * fpi * dipole_in / area
    !
    ! ... Build the constant approximation of rhoioncc, rescale the guess
    !
    rhozero(:) = rho(:)
    charge_out = SUM(rhozero(:))*darea
    rhozero = rhozero*ABS(charge_in/charge_out)
    IF ( verbose .GE. 3 ) CALL write_1d( n, axis, rhozero, 'rhozero1d.dat' )
    !
    ! ... Self-consistent iterations
    !
    DO iter = 1, maxiter
      rhoin = rho
      rhotot = rhofix + rhozero + rhoin
      CALL v1d_of_rho1d( n, dx, rhotot, v ) 
      v(:) = v(:) + dv * step(:)
      IF ( verbose .GE. 4 ) CALL write_1d( n, axis, v, 'vtmp1d.dat' )
      CALL generate_rhoioncc( n, 1, switch, v, rhoout )
      IF ( verbose .GE. 4 ) CALL write_1d( n, axis, rhoout, 'rhoout1d.dat' )
      charge_out = SUM(rhoout)*darea 
      mixlocal = mix
      IF ( ABS(charge_out-charge_in) .GT. 1.D-3 ) &
        mixlocal = mix*ABS(charge_in/charge_out)*1.d-3
      rhoout = rhoout - rhozero
      residue = rhoout - rhoin
      rho = mixlocal * rhoout + (1.D0 - mixlocal) * rhoin
      deltarho = SUM(residue(:)*residue(:))
      deltarho = SQRT(deltarho) / DBLE(n1d)
      IF ( verbose .GE. 3 ) WRITE(environ_unit,8201)iter,deltarho,charge_out
      IF (deltarho.LT.tol) THEN
        IF ( verbose .GE. 1 ) WRITE(environ_unit,8202)
        EXIT
      ELSE IF ( iter .EQ. maxiter ) THEN
        WRITE(stdout,8203)
      ENDIF
    ENDDO  
    IF ( verbose .GE. 1 ) WRITE(environ_unit,8204)deltarho,iter
    !
    rho = rho + rhozero + rhofix
    CALL v1d_of_rho1d( n, dx, rho, v ) 
    IF ( verbose .GE. 3 ) CALL write_1d( n, axis, v, 'v1d.dat' )
    rho = rho - rhofix
    v = v - vfix - v0 + v1d0 
    IF ( verbose .GE. 3 ) CALL write_1d( n, axis, v, 'vext1d.dat' )
    !
    RETURN
    !
8201 FORMAT(1X,'iter = ',i10,' deltarho = ',E14.6,' charge_out = ',E14.6)
8202 FORMAT(1X,'Level 2 ioncc charges converged!')
8203 FORMAT(1X,'Warning: Level 2 ioncc charges are not converged')
8204 FORMAT(1X,'        L2 ioncc accuracy =',1PE8.1,', # of iterations = ',i6)
!--------------------------------------------------------------------
  END SUBROUTINE calc_rhoioncc1d
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_rhoioncc3d( nnr, nspin, mix, tol, switch, vfix, rhozero, rhoioncc, vioncc, failed )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, PARAMETER :: maxiter = 100000
    !
    INTEGER, INTENT(IN) :: nnr, nspin
    !
    REAL( DP ), INTENT(IN) :: mix, tol
    REAL( DP ), INTENT(IN) :: switch(nnr), vfix(nnr), rhozero(nnr)
    REAL( DP ), INTENT(INOUT) :: rhoioncc(nnr)
    REAL( DP ), INTENT(OUT) :: vioncc(nnr)
    !
    LOGICAL, INTENT(OUT) :: failed
    !
    INTEGER :: iter
    !
    REAL( DP ) :: charge, ehart, charge_iter, charge_out, mixlocal, deltarho, deltaold
    REAL( DP ), ALLOCATABLE, DIMENSION(:) :: rhoin, rhoout, residue
    REAL( DP ), ALLOCATABLE, DIMENSION(:) :: vperiodic 
    REAL( DP ), ALLOCATABLE, DIMENSION(:,:) :: rhotot, vtot
    !
    failed = .FALSE.
    !
    ALLOCATE(rhotot(nnr,nspin))
    ALLOCATE(vtot(nnr,nspin))
    ALLOCATE(vperiodic(nnr))
    rhotot = 0.D0
    vtot = 0.D0
    !
    ! ... 3D Refinement, iterating on the difference of the 1D solution and the 3D solution
    !     assume there is no difference (wishfull tihnking)
    !
    rhoioncc = 0.D0 
    !
    ALLOCATE(rhoin(nnr))
    ALLOCATE(rhoout(nnr))
    ALLOCATE(residue(nnr))
    deltarho = 0.D0
    DO iter = 1, maxiter
      !
      rhoin = rhoioncc
      !
      ! ... Compute the potential coming from the iterative charge
      !
      rhotot( :, 1 ) = rhoin( : )
      IF ( nspin .EQ. 2 ) rhotot( :, 2 ) = 0.D0
      vtot = 0.D0
      CALL v_h_of_rho_r( rhotot, ehart, charge, vtot )
      vperiodic = 0.D0
      CALL calc_vperiodic( nnr, nspin, .FALSE., rhoin, vperiodic )
      vtot( :, 1 ) = vtot( :, 1 ) + vperiodic( : ) + vfix( : ) 
      IF ( verbose .GE. 4 ) CALL write_cube( nnr, vtot(:,1), 'vtottmp.cube') 
      !
      CALL generate_rhoioncc( nnr, 1, switch, vtot(:,1), rhoout )
      !
      rhoout = (rhoout - rhozero)
      IF ( verbose .GE. 4 ) CALL write_cube( nnr, rhoout, 'rhoout.cube') ! DEBUG ONLY, CHANGE TO VERBOSE .GE. 4
      !
      residue = rhoout - rhoin
      deltaold = deltarho
      deltarho = SUM(residue(:)*residue(:))
      CALL mp_sum( deltarho, intra_bgrp_comm )
      deltarho = SQRT(deltarho)/DBLE(ntot)
      charge_out = SUM(rhoout)*domega
      CALL mp_sum( charge_out, intra_bgrp_comm )
      ! 
      ! ... If the predicted charge is very large, use a tiny mixing
      !
      mixlocal = mix * MIN(1.D0,ABS(1.D0/charge_out))
      IF ( deltarho .GT. 1.D10 ) mixlocal = MIN(mixlocal*10.D0,0.1)
      !
      CALL mix_rhoioncc_gspace( nnr, 1, mixlocal, 1.D0, rhoout, rhoin, rhoioncc ) 
!      rhoioncc = mixlocal * rhoout + ( 1.D0 - mixlocal ) * rhoin ! seems to be working only if passing through G-space...
      IF ( verbose .GE. 1 ) WRITE(environ_unit,8101)iter,deltarho,deltarho/deltaold,charge_out,mixlocal
      IF ( verbose .GE. 4 ) CALL write_cube( nnr, rhoioncc, 'rhoioncc.cube') ! DEBUG ONLY, CHANGE TO VERBOSE .GE. 4
      !
      IF (deltarho.LT.tol) THEN
        IF ( verbose .GE. 1 ) WRITE(environ_unit,8102)
        EXIT
      ELSE IF ( iter .GT. 1 .AND. (deltarho/deltaold .GT. 1.D1) ) THEN
        WRITE(stdout,8106)
        failed = .TRUE.
        RETURN
      ELSE IF ( iter .EQ. maxiter ) THEN
        WRITE(stdout,8103)
      ENDIF
      !
    ENDDO
    WRITE(stdout,8104)deltarho,iter 
    DEALLOCATE(rhoin,rhoout,residue)
    !
    IF ( verbose .GE. 1 ) CALL write_cube( nnr, rhoioncc(:), 'rhoioncc.cube') ! DEBUG ONLY CHANGE TO VERBOSE .GE. 3
    vioncc(:) = vtot( :, 1 )
    !
    DEALLOCATE(rhotot,vtot,vperiodic)
    !
    RETURN
    !
8101 FORMAT(1X,'iter = ',I5,' deltarho = ',E14.6,' drho/dold = ',E14.6,' charge_out = ',E14.6,' mixlocal = ',E14.6)
8102 FORMAT(1X,'Level 3 ioncc charges are converged!')
8103 FORMAT(1X,'Warning: Level 3 ioncc charges are not converged')
8104 FORMAT(1X,'        L3 ioncc accuracy =',1PE8.1,', # of iterations = ',i5)
8105 FORMAT(1X,'iter = ',I5,' charge_tmp = ',E14.6,' charge_ext = ',E14.6)
8106 FORMAT(1X,'Warning: Level 3 ioncc charges are exploding, trying with a different mix')
!--------------------------------------------------------------------
  END SUBROUTINE calc_rhoioncc3d
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_rhoioncc( nnr, nspin, n1d,  x1dmin, x1dmax, spread1d, x1d, &
                          & axis1d, step1d, rhoioncc1d, switch, rhofix, rhoioncc, vioncc )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, PARAMETER :: maxiter = 100000
    REAL( DP ), PARAMETER :: mix = 0.1
    REAL( DP ), PARAMETER :: tol = 1.D-10
    !
    INTEGER, INTENT(IN) :: nnr, nspin, n1d
    !
    REAL( DP ), INTENT(IN) :: x1dmin, x1dmax, spread1d
    REAL( DP ), INTENT(INOUT) :: x1d
    !
    REAL( DP ), INTENT(IN) :: axis1d(n1d), step1d(n1d)
    REAL( DP ), INTENT(INOUT) :: rhoioncc1d(n1d)
    REAL( DP ), INTENT(IN) :: switch(nnr), rhofix(nnr)
    REAL( DP ), INTENT(INOUT) :: rhoioncc(nnr)
    REAL( DP ), INTENT(OUT) :: vioncc(nnr)
    !
    INTEGER :: imix
    ! 
    LOGICAL :: failed
    LOGICAL, SAVE :: first = .TRUE.
    !
    REAL( DP ) :: charge_in, charge_out, charge_zero, mixlocal, charge, ehart, deltarho, fact
    REAL( DP ), ALLOCATABLE, DIMENSION(:) :: rhozero
    REAL( DP ), ALLOCATABLE, DIMENSION(:) :: vperiodic, vfix, vzero
    REAL( DP ), ALLOCATABLE, DIMENSION(:,:) :: rhotot, vtot
    !
    REAL( DP ) :: rhotot1d(n1d), vioncc1d(n1d)
    !
    REAL( DP ) :: x1dtest, xmax, tollocal, minx1derror
    !
    ALLOCATE(rhotot(nnr,nspin))
    ALLOCATE(rhozero(nnr))
    ALLOCATE(vtot(nnr,nspin))
    ALLOCATE(vperiodic(nnr))
    ALLOCATE(vzero(nnr))
    ALLOCATE(vfix(nnr))
    vtot = 0.D0
    vzero = 0.D0
    rhozero = 0.D0
    vfix = 0.D0
    !
    ! ... Compute the non-periodic potential of the source charge
    !
    charge_in = SUM(rhofix)*domega
    CALL mp_sum( charge_in, intra_bgrp_comm )
    rhotot(:,1) = rhofix(:)
    IF ( nspin .EQ. 2 ) rhotot(:,2) = 0.D0
    CALL v_h_of_rho_r( rhotot, ehart, charge, vtot )
    vperiodic = 0.D0
    CALL calc_vperiodic( nnr, nspin, .FALSE., rhotot(:,1), vperiodic )
    vfix( : ) = vtot( :, 1 ) + vperiodic(:)  
    vavg = SUM(vfix)/DBLE(ntot)
    CALL mp_sum( vavg, intra_bgrp_comm )
    write(environ_unit,*)'vavg ',vavg
    IF ( verbose .GE. 3 ) CALL write_cube( nnr, rhofix(:), 'rhofix.cube')
    IF ( verbose .GE. 1 ) CALL write_cube( nnr, vfix(:), 'vfix.cube') ! DEBUG, reset to 3
    !
    rhotot1d = 0.D0
    CALL planar_average( nnr, naxis, slab_axis, shift, .false., rhotot(:,1), rhotot1d(cell_min:cell_max) )
    !
    ! ... If first step, need to find the optimal 1D solution
    !
    IF ( first ) THEN
      first = .FALSE.
      !
      ! ... the initial guess is in the golden ratio between x1dmax and x1dmin
      !
      x1dtest = x1dmin + .61803399 * ( x1dmax - x1dmin ) 
      xmax = x1dtest + deltastern
      WRITE( environ_unit, * )'x1dmax = ',x1dmax,' x1dmin = ',x1dmin,' x1dtest = ',x1dtest
      CALL generate_switch( n1d, 2, x1dtest, xmax, 1.D0, axis1d, switch1d )
      IF ( verbose .GE. 1 ) CALL write_1d( n1d, axis1d, switch1d, 'switch1d_first.dat' )
      CALL calc_rhoioncc1d( n1d, vavg, axis1d, step1d, switch1d, rhotot1d, vioncc1d, rhoioncc1d )
      CALL planar_average( nnr, naxis, slab_axis, shift, .true., rhozero, rhoioncc1d(cell_min:cell_max) )
      CALL planar_average( nnr, naxis, slab_axis, shift, .true., vzero, vioncc1d(cell_min:cell_max) )
      IF ( verbose .GE. 1 ) CALL write_cube( nnr, rhozero(:), 'rhozero.cube') ! DEBUG, reset it to 3
      IF ( verbose .GE. 1 ) CALL write_cube( nnr, vzero(:), 'vzero.cube') ! DEBUG, reset it to 3
      charge_zero = SUM(rhozero)*domega
      CALL mp_sum( charge_zero, intra_bgrp_comm )
      !
      ! ... try to converge level 3 with lower tol
      !
      vfix = vfix + vzero
      mixlocal = mix * 5.D-2
      tollocal = 1.D-6
      DO imix = 1, 10
        IF ( verbose .GE. 1 ) WRITE( environ_unit, 7001 ) mixlocal 
        CALL calc_rhoioncc3d( nnr, nspin, mixlocal, tollocal, switch, vfix, rhozero, rhoioncc, vioncc, failed )
        IF ( failed ) THEN
          mixlocal = mixlocal * 0.5D0
        ELSE
          GOTO 201
        ENDIF
      ENDDO
      WRITE( stdout, * )"ERROR: calc_rhoioncc failed to converged even with the smallest mix" 
      STOP
201   CONTINUE 
      vfix = vfix - vzero
      !
      ! ... use the converged density to get the optimal x1d
      !
      ALLOCATE(rhoioncc1dref(n1d))
      ALLOCATE(rhotot1dref(n1d))
      rhoioncc1dref = 0.D0
      CALL planar_average( nnr, naxis, slab_axis, shift, .false., rhoioncc, rhoioncc1dref(cell_min:cell_max) )
      rhoioncc1dref = rhoioncc1dref + rhoioncc1d
      rhotot1dref = rhotot1d
      minx1derror = golden(x1dmin,x1dtest,x1dmax,tol,x1d)
      DEALLOCATE(rhoioncc1dref)
      DEALLOCATE(rhotot1dref)
      !
    END IF 
    !
    !
    !
    xmax = x1d + deltastern
    CALL generate_switch(n1d,2,x1d,xmax,spread1d,axis1d,switch1d)
    IF ( verbose .GE. 1 ) CALL write_1d( n1d, axis1d, switch1d, 'switch1d.dat' ) ! DEBUG, reset it to 3
    !
    ! ... Compute the 1D solution 
    ! 
    CALL calc_rhoioncc1d( n1d, vavg, axis1d, step1d, switch1d, rhotot1d, vioncc1d, rhoioncc1d )
    CALL planar_average( nnr, naxis, slab_axis, shift, .true., rhozero, rhoioncc1d(cell_min:cell_max) )
    CALL planar_average( nnr, naxis, slab_axis, shift, .true., vzero, vioncc1d(cell_min:cell_max) )
    charge_zero = SUM(rhozero)*domega
    CALL mp_sum( charge_zero, intra_bgrp_comm )
    IF ( verbose .GE. 1 ) CALL write_cube( nnr, rhozero(:), 'rhozero.cube') ! DEBUG, reset it to 3
    IF ( verbose .GE. 1 ) CALL write_cube( nnr, vzero(:), 'vzero.cube') ! DEBUG, reset it to 3
    !
    ! ... Compute the 3D refinement, decreasing mixing if it explodes
    !
    vfix = vfix + vzero
    mixlocal = mix * 5.D-2
    DO imix = 1, 10
      IF ( verbose .GE. 1 ) WRITE( environ_unit, 7001 ) mixlocal 
      CALL calc_rhoioncc3d( nnr, nspin, mixlocal, tol, switch, vfix, rhozero, rhoioncc, vioncc, failed )
      IF ( failed ) THEN
        mixlocal = mixlocal * 0.5D0
      ELSE
        GOTO 200
      ENDIF
    ENDDO
    WRITE( stdout, * )"ERROR: calc_rhoioncc failed to converged even with the smallest mix" 
    STOP
    !
    ! ... Converged 3D refinement
    !
200 CONTINUE 
    charge_out = SUM(rhoioncc)*domega
    CALL mp_sum( charge_out, intra_bgrp_comm )
    !
    vioncc = vioncc + vzero + vperiodic
    !
    RETURN
    !
7001 FORMAT(1X,'Trying to compute level 3 ioncc density with mixing =',E14.6)
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_rhoioncc
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_rhoioncc( n, type, switch, v, rho ) 
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL( DP ), PARAMETER :: max_arg = 2.D1
    !
    INTEGER, INTENT(IN) :: n, type
    REAL( DP ), INTENT(IN) :: switch(n)
    REAL( DP ), INTENT(IN) :: v(n)
    REAL( DP ), INTENT(INOUT) :: rho(n)
    !
    INTEGER :: i
    REAL( DP ) :: arg, fact
    !
    fact = -2.D0 * zion * cion 
    DO i = 1, n
      IF (ABS(switch(i)).LT.1.D-10) THEN
        rho(i) = 0.D0
      ELSE 
        arg = zion * invkbt * v(i)
        IF ( arg .GT. max_arg ) arg = max_arg
        IF ( type .EQ. 1 ) THEN
          rho(i) = fact * sinh( arg ) * switch(i)
        ELSE IF ( type .EQ. 2 ) THEN
          rho(i) = fact * ( arg ) * switch(i)
        ELSE IF ( type .EQ. 3 ) THEN
          rho(i) = fact * ( arg  + arg**3 / 6.D0 ) * switch(i)
        ELSE IF ( type .EQ. 4 ) THEN
          rho(i) = fact * ( arg  + arg**3 / 6.D0 + arg**5 / 120D0 ) * switch(i)
        END IF 
      ENDIF
    ENDDO
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_rhoioncc
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_eioncc(  nnr, rho, eioncc )
!--------------------------------------------------------------------
    !
    USE environ_base,  ONLY : vioncc
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nnr
    REAL( DP ), INTENT(IN) :: rho(nnr)
    REAL( DP ), INTENT(OUT) :: eioncc
    !
    REAL( DP ), ALLOCATABLE :: rhotot(:)
    !
    ALLOCATE(rhotot(nnr))
    rhotot = rhoions + rho
    !
    eioncc = 0.5D0 * SUM( vioncc(:) * rhotot( : ) ) * domega
    !
    DEALLOCATE(rhotot)
    !
    CALL mp_sum( eioncc, intra_bgrp_comm )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_eioncc
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE write_1d( n1d, axis1d, f1d, filename )
!--------------------------------------------------------------------
      !
      USE kinds,          ONLY : DP
      USE io_global,      ONLY : ionode
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)      :: n1d
      REAL( DP ), INTENT(IN)   :: axis1d( n1d )
      REAL( DP ), INTENT(IN)   :: f1d( n1d )
      !
      CHARACTER(*), INTENT(IN) :: filename
      !
      INTEGER                  :: ia
      !
      IF( ionode ) THEN
        !
        OPEN( 300, file = TRIM( filename ), status = 'unknown' )
        !
        DO ia = 1, n1d
          WRITE( 300, '(f20.10,x,e20.10)') axis1d(ia), f1d(ia)
        ENDDO
        !
        CLOSE( 300 )
        !
      END IF
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE write_1d
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE generate_step( n, axis, step )
!--------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: n
      !
      REAL( DP ), INTENT(IN) :: axis(n)
      REAL( DP ), INTENT(OUT) :: step(n)
      !
      INTEGER :: i
      !
      DO i = 1, n
        IF ( axis(i) .GT. 0.D0 ) THEN
          step(i) = 0.5D0
        ELSE
          step(i) = -0.5D0
        ENDIF
      ENDDO
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE generate_step
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE generate_switch( n, ityp, xmin, xmax, smoothness, axis, switch )
!--------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: n, ityp
      !
      REAL( DP ), INTENT(IN) :: xmin, xmax, smoothness
      REAL( DP ), INTENT(IN) :: axis(n)
      REAL( DP ), INTENT(OUT) :: switch(n)
      !
      INTEGER :: i
      !
      REAL( DP ) :: arg, fact
      !
      switch = 0.D0
      SELECT CASE( ityp )
      CASE ( 1 ) 
        fact = tpi / (xmax - xmin)
        DO i = 1, n
          IF ( ABS(axis(i)) .LE. xmin ) THEN
            switch(i) = 0.D0
          ELSE IF ( ABS(axis(i)) .LT. xmax ) THEN
            arg = fact * ( ABS(axis(i)) - xmin )
            switch(i) = ( arg - SIN( arg ) ) / tpi
          ELSE
            switch(i) = 1.D0
          ENDIF   
        ENDDO
      CASE ( 2 )
        DO i = 1, n
          IF ( ABS(axis(i)) .GT. xmin ) THEN
            arg = -smoothness/(ABS(axis(i))-xmin)**2
            switch(i) = EXP(arg)
          ENDIF
        ENDDO
      CASE ( 3 )
        DO i = 1, n
          IF ( ABS(axis(i)) .LE. xmin ) THEN
            switch(i) = 0.D0
          ELSE IF ( ABS(axis(i)) .LT. xmax ) THEN
            arg = (ABS(axis(i))-xmax)*smoothness
            switch(i) = EXP(arg)
          ELSE 
            switch(i) = 1.D0
          ENDIF
        ENDDO
      END SELECT
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE generate_switch
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE v1d_of_rho1d( n, dx, rho, v )
!--------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: n
      !
      REAL( DP ), INTENT(IN) :: dx
      REAL( DP ), INTENT(IN) :: rho(n)
      REAL( DP ), INTENT(OUT) :: v(n)
      !
      REAL( DP ) :: e(n)
      !
      CALL integral_function_1d(n,dx,rho,e)
      CALL integral_function_1d(n,dx,e,v)
      v = - fpi * e2 * v
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE v1d_of_rho1d
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE integral_function_1d( n, dx, f, fint )
!--------------------------------------------------------------------
      !
      USE kinds,          ONLY : DP
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)      :: n
      REAL( DP ), INTENT(IN)   :: dx
      REAL( DP ), INTENT(IN)   :: f( n )
      REAL( DP ), INTENT(OUT)  :: fint( n )
      !
      INTEGER                  :: i
      !
      REAL( DP ) :: c1, c2, c3, c4, shift
      !
      REAL( DP ) :: ftmp(n)
      !
      c1 = 3.D0/8.D0
      c2 = 7.D0/6.D0 - 3.D0/8.D0
      c3 = 23.D0/24.D0 - 7.D0/6.D0
      c4 = 1.D0 - 23.D0/24.D0
      !
      fint = 0.d0
      !
      fint(2) = 0.5*(f(1)+f(2))
      !
      fint(3) = (f(1)+f(3)+4.D0*f(2))/3.D0
      !
      fint(4) = (f(1)+f(4)+3.D0*(f(2)+f(3)))*3.D0/8.D0
      !
      fint(5) = (14.D0*(f(1)+f(5))+64.D0*(f(2)+f(4))+24*f(3))/45.D0
      !
      fint(6) = 3.D0/8.D0*(f(1)+f(6))+7.D0/6.D0*(f(2)+f(5))+23.D0/24.D0*(f(3)+f(4))
      DO i = 7, n
        fint(i) = fint(i-1) + f(i)*c1 + f(i-1)*c2 + f(i-2)*c3 + f(i-3)*c4
      ENDDO
      !
      ftmp = 0.D0
      !
      ftmp(n-1) = 0.5*(f(n)+f(n-1))
      !
      ftmp(n-2) = (f(n)+f(n-2)+4.D0*f(n-1))/3.D0
      !
      ftmp(n-3) = (f(n)+f(n-3)+3.D0*(f(n-1)+f(n-2)))*3.D0/8.D0
      !
      ftmp(n-4) = (14.D0*(f(n)+f(n-4))+64.D0*(f(n-1)+f(n-3))+24*f(n-2))/45.D0
      !
      ftmp(n-5) = 3.D0/8.D0*(f(n)+f(n-5))+7.D0/6.D0*(f(n-1)+f(n-4))+23.D0/24.D0*(f(n-2)+f(n-3))
      DO i = n-6, 1, -1
        ftmp(i) = ftmp(i+1) + f(i)*c1 + f(i+1)*c2 + f(i+2)*c3 + f(i+3)*c4
      ENDDO
      !
      fint = ( fint - ftmp )/2.D0
      !
!      shift = (fint(1) + fint(n))/2.D0
!      fint = ( fint - shift ) * dx
      fint = fint * dx
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE integral_function_1d
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE mix_rhoioncc_gspace( nnr, nspin, alpha, beta, rhoout, rhoin, rhoioncc )
!--------------------------------------------------------------------

      USE fft_base,       ONLY : dfftp
      USE fft_interfaces, ONLY : fwfft, invfft
      USE gvect,          ONLY : nl, nlm, ngm, g, gg, gstart
      USE control_flags,  ONLY : gamma_only
 
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nnr, nspin
      REAL(DP), INTENT(IN) :: alpha, beta
      
      REAL(DP), DIMENSION(nnr,nspin), INTENT(IN) :: rhoout, rhoin
      REAL(DP), DIMENSION(nnr,nspin), INTENT(OUT) :: rhoioncc
      !
      COMPLEX( DP ), ALLOCATABLE :: rhooutg(:,:), rhoing(:,:), rhoionccg(:,:)
      COMPLEX( DP ), ALLOCATABLE :: aux( : )
      INTEGER :: is, ig
      REAL(DP) :: fac
 
      ALLOCATE( rhooutg(ngm,nspin), rhoing(ngm,nspin), rhoionccg(ngm,nspin))
      ALLOCATE( aux(nnr) )

      DO is = 1, nspin
        aux(:) = CMPLX(rhoout( : , is ),0.D0,kind=dp) 
        CALL fwfft ('Dense', aux, dfftp)
        rhooutg(:,is) = aux(nl(:))
        aux(:) = CMPLX(rhoin( : , is ),0.D0,kind=dp) 
        CALL fwfft ('Dense', aux, dfftp)
        rhoing(:,is) = aux(nl(:))
      END DO

      rhoionccg = 0.D0
      DO ig = gstart, ngm
!        fac = alpha * (beta*gg(ig))
        fac = alpha !/ (1.D0+beta*gg(ig))
!        fac = alpha / (1.D0+beta*g(3,ig)**2)
        rhoionccg(ig,1) = fac * rhooutg(ig,1) + ( 1.D0 - fac ) * rhoing(ig,1)
        IF ( nspin == 2 ) &
          & rhoionccg(ig,2) = fac * rhooutg(ig,2) + ( 1.D0 - fac ) * rhoing(ig,2)   
      ENDDO
      
      DO is = 1, nspin
        aux(nl(1:ngm)) = rhoionccg(1:ngm,is)
        IF ( gamma_only ) &
          & aux(nlm(1:ngm)) = CMPLX(REAL(rhoionccg(1:ngm,is)),-AIMAG(rhoionccg(1:ngm,is)),KIND=dp)
        CALL invfft ('Dense', aux, dfftp)
        rhoioncc(:,is) = DBLE(aux(:))   
      ENDDO

      DEALLOCATE( aux )
      DEALLOCATE( rhooutg, rhoing, rhoionccg )

      RETURN
!--------------------------------------------------------------------
      END SUBROUTINE mix_rhoioncc_gspace
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE x1derror( x1d , f )
!--------------------------------------------------------------------
      IMPLICIT NONE

      REAL(DP) :: x1d
      REAL(DP) :: f
      REAL(DP) :: xmax
      REAL(DP), DIMENSION(n1d) :: switchloc, rhoionccloc, vionccloc

      xmax = x1d + deltastern
      CALL generate_switch( n1d, 2, x1d, xmax, 1.D0, axis1d, switchloc )
      CALL calc_rhoioncc1d( n1d, vavg, axis1d, step1d, switchloc, rhotot1dref, vionccloc, rhoionccloc )
      
      rhoionccloc = rhoionccloc - rhoioncc1dref
      f = SQRT(SUM(rhoionccloc(:)**2))

      RETURN
!--------------------------------------------------------------------
      END SUBROUTINE x1derror
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION golden(ax,bx,cx,tol,xmin)
!--------------------------------------------------------------------
      IMPLICIT NONE
      REAL(DP) :: golden,ax,bx,cx,tol,xmin,R,C
      PARAMETER (R=.61803399,C=1.-R)
!Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
!between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine performs
!a golden section search for the minimum, isolating it to a fractional precision of about
!tol. The abscissa of the minimum is returned as xmin, and the minimum function value
!is returned as golden, the returned function value.
!Parameters: The golden ratios.
      REAL(DP) f1,f2,x0,x1,x2,x3
        write( environ_unit, *)'golden search'
      x0=ax
      x3=cx
      if(abs(cx-bx).gt.abs(bx-ax))then ! Make x0 to x1 the smaller segment,
        x1=bx
        x2=bx+C*(cx-bx) ! and fill in the new point to be tried.
      else
        x2=bx
        x1=bx-C*(bx-ax)
      endif
      CALL x1derror(x1,f1) ! The initial function evaluations. Note that we never need to
      CALL x1derror(x2,f2) ! evaluate the function at the original endpoints.
1     if(abs(x3-x0).gt.tol*(abs(x1)+abs(x2))) then ! Do-while loop: we keep returning here.
        write( environ_unit, *)x0,x1,x2,x3
        write( environ_unit, *)f1,f2
        if(f2.lt.f1)then ! One possible outcome,
          x0=x1 ! its housekeeping,
          x1=x2
          x2=R*x1+C*x3
          f1=f2
          CALL x1derror(x2,f2) ! and a new function evaluation.
        else ! The other outcome,
          x3=x2
          x2=x1
          x1=R*x2+C*x0
          f2=f1
          CALL x1derror(x1,f1) ! and its new function evaluation.
        endif
        goto 1 ! Back to see if we are done.
      endif
      if(f1.lt.f2)then ! We are done. Output the best of the two current values.
        golden=f1
        xmin=x1
      else
        golden=f2
        xmin=x2
      endif
      stop
      return 
!--------------------------------------------------------------------
      END FUNCTION golden
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE ioncc
!--------------------------------------------------------------------
