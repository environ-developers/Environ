!!
!! Copyright (C) 2007-2008 Quantum-ESPRESSO group
!! This file is distributed under the terms of the
!! GNU General Public License. See the file `License'
!! in the root directory of the present distribution,
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! original version by O. Andreussi and N. Marzari
!!
!!--------------------------------------------------------------------
!MODULE extcharge
!!--------------------------------------------------------------------
!
!  USE kinds,          ONLY: DP
!  USE io_global,      ONLY: stdout
!  USE mp,             ONLY: mp_sum
!  USE mp_bands,       ONLY: intra_bgrp_comm
!  USE environ_cell,   ONLY: domega
!  USE environ_ions,   ONLY: avg_pos, rhoions
!  USE environ_base,   ONLY: verbose, environ_unit,                   &
!                            system_pos, system_width,                &
!                            env_external_charges, extcharge_origin,  &
!                            extcharge_dim, extcharge_axis,           &
!                            extcharge_pos, extcharge_spread,         &
!                            extcharge_charge, rhoexternal
!  USE environ_debug,  ONLY: write_cube
!  !
!  IMPLICIT NONE
!  !
!  LOGICAL :: first = .TRUE.
!  !
!  REAL(DP) :: eextself = 0.D0
!  !
!  SAVE
!  !
!  PRIVATE
!  !
!  PUBLIC :: generate_extcharge, calc_vextcharge, calc_eextcharge, calc_fextcharge
!  !
!CONTAINS
!  !
!!--------------------------------------------------------------------
!  SUBROUTINE generate_extcharge( nnr, alat )
!!--------------------------------------------------------------------
!    !
!    USE generate_function, ONLY : generate_gaussian
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: nnr
!    REAL(DP) :: alat
!    !
!    LOGICAL :: shift
!    INTEGER :: iextcharge
!    INTEGER :: axis, dim
!    REAL(DP) :: charge, spread, norm, pos(3), pos0(3)
!    !
!    rhoexternal = 0.D0
!    !
!    IF ( env_external_charges .LE. 0 ) RETURN
!    !
!    first = .TRUE.
!    !
!    ! Use the origin of the cell as default origin, removed other options
!    !
!    SELECT CASE ( extcharge_origin )
!    CASE ( 0 )
!       shift = .false.
!       pos0 = 0.D0
!    CASE ( 1 )
!       shift = .false.
!       pos0 = system_pos
!    CASE ( 2 )
!       shift = .true.
!       pos0 = system_pos
!    END SELECT
!    !
!    ! Generate gaussian densities (planar,lines or points, depending on dim)
!    !
!    DO iextcharge = 1, env_external_charges
!      !
!      norm = SQRT(SUM(extcharge_pos(:,iextcharge)**2))
!      IF ( .NOT. shift ) THEN
!        pos(:) = extcharge_pos(:,iextcharge)/alat + pos0(:)
!      ELSE IF ( norm .NE. 0.D0 ) THEN
!        pos(:) = extcharge_pos(:,iextcharge)/alat*(1.D0+system_width/norm) + pos0(:)
!      ELSE
!        WRITE(stdout,*)'ERROR: ill-defined position of external charge'
!        STOP
!      ENDIF
!      charge = extcharge_charge(iextcharge)
!      spread = extcharge_spread(iextcharge)
!      dim = extcharge_dim(iextcharge)
!      axis = extcharge_axis(iextcharge)
!      !
!      ! Extcharge densities only for points, lines or planes (DIM = 0, 1 or 2)
!      !
!      IF ( dim .GE. 3 .OR. dim .LT. 0 ) THEN
!        WRITE(stdout,*)'Warning: wrong extcharge_dim, doing nothing'
!        CYCLE
!      END IF
!      !
!      CALL generate_gaussian( nnr, dim, axis, charge, spread, pos, rhoexternal )
!      !
!    ENDDO
!    !
!    IF ( verbose .GE. 2 ) CALL write_cube( nnr, rhoexternal, 'rhoexternal.cube' )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE generate_extcharge
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE calc_vextcharge(  nnr, nspin, vextcharge )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: nnr, nspin
!    REAL( DP ), INTENT(OUT) :: vextcharge(nnr)
!    !
!    REAL( DP ) :: ehart, charge
!    REAL( DP ), ALLOCATABLE :: vaux(:,:), rhoaux(:,:)
!
!    CALL start_clock( 'get_extcharge' )
!
!    IF ( .NOT. first ) RETURN
!    first = .FALSE.
!    vextcharge = 0.D0
!
!    ALLOCATE( rhoaux( nnr, nspin ) )
!    ALLOCATE( vaux( nnr, nspin ) )
!    rhoaux( :, 1 ) = rhoexternal(:)
!    IF ( nspin .EQ. 2 ) rhoaux( :, 2 ) = 0.D0
!    vaux = 0.D0
!    CALL v_h_of_rho_r( rhoaux, ehart, charge, vaux )
!    vextcharge(:) = vaux( :, 1 )
!    DEALLOCATE( rhoaux )
!    IF ( verbose .GE. 2 ) CALL write_cube( nnr, vextcharge, 'vextcharge.cube' )
!    DEALLOCATE( vaux )
!
!    eextself = 0.5 * SUM( vextcharge( : ) * rhoexternal( : ) ) * domega
!    CALL mp_sum( eextself, intra_bgrp_comm )
!    IF ( verbose .GE. 1 ) WRITE(environ_unit,*)&
!      & 'External charge self energy',eextself
!
!    CALL stop_clock( 'get_extcharge' )
!
!    RETURN
!
!!--------------------------------------------------------------------
!  END SUBROUTINE calc_vextcharge
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE calc_eextcharge(  nnr, rho, eextcharge )
!!--------------------------------------------------------------------
!    !
!    USE environ_base,  ONLY : vextcharge
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: nnr
!    REAL( DP ), INTENT(IN) :: rho(nnr)
!    REAL( DP ), INTENT(OUT) :: eextcharge
!    !
!    REAL( DP ), ALLOCATABLE :: rhotot(:)
!    !
!    ALLOCATE(rhotot(nnr))
!    rhotot = rhoions + rho
!    !
!    eextcharge = SUM( vextcharge(:) * rhotot( : ) ) * domega
!    !
!    DEALLOCATE(rhotot)
!    !
!    CALL mp_sum( eextcharge, intra_bgrp_comm )
!    !
!    eextcharge = eextcharge + eextself
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE calc_eextcharge
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE calc_fextcharge(  nnr, nspin, nat, f )
!!--------------------------------------------------------------------
!    !
!    USE generate_function, ONLY : generate_gradgaussian
!    USE environ_cell,      ONLY : alat
!    USE environ_base,      ONLY : env_static_permittivity, rhopol
!    USE environ_ions,      ONLY : ntyp, ityp, zv, tau, avg_pos, zvtot, rhoions
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN)       :: nnr, nspin, nat
!    REAL(DP), INTENT(INOUT)   :: f( 3, nat )
!    !
!    ! Since the external charges can only be placed independently of
!    ! atomic positions there is no explicit contributions to the
!    ! inter-atomic forces
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE calc_fextcharge
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!END MODULE extcharge
!!--------------------------------------------------------------------
