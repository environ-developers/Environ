!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by O. Andreussi and N. Marzari
!
!--------------------------------------------------------------------
MODULE externals_utils
!--------------------------------------------------------------------

  USE environ_types
  USE functions
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
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: create_environ_externals, init_environ_externals_first, &
       & init_environ_externals_second, update_environ_externals, destroy_environ_externals
  !
  CONTAINS
  SUBROUTINE create_environ_externals(externals)

    IMPLICIT NONE

    TYPE( environ_externals ), INTENT(INOUT) :: externals
    CHARACTER (LEN=80) :: sub_name = 'create_environ_externals'

    externals%update = .FALSE.
    externals%number = 0
    IF ( ALLOCATED( externals%functions ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)
    CALL create_environ_density( externals%density, "externals" )
    externals%charge = 0.D0

    RETURN

  END SUBROUTINE create_environ_externals

  SUBROUTINE init_environ_externals_first( nexternals, dims, axis, pos, &
       & spreads, charge, externals )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nexternals
    INTEGER, DIMENSION(nexternals), INTENT(IN) :: dims, axis
    REAL( DP ), DIMENSION(3,nexternals), INTENT(IN) :: pos
    REAL( DP ), DIMENSION(nexternals), INTENT(IN) :: spreads, charge
    TYPE( environ_externals ), INTENT(INOUT) :: externals

    INTEGER :: i

    externals%number = nexternals
    ALLOCATE(externals%functions(externals%number))
    DO i = 1, externals%number
       ALLOCATE(externals%functions(i)%pos(3))
       externals%functions(i)%dim    = dims(i)
       externals%functions(i)%axis   = axis(i)
       externals%functions(i)%pos(:) = pos(:,i)
       externals%functions(i)%spread = spreads(i)
       externals%functions(i)%width  = spreads(i)
       externals%functions(i)%volume = charge(i)
    ENDDO

    RETURN

  END SUBROUTINE init_environ_externals_first

  SUBROUTINE init_environ_externals_second( cell, externals )

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_externals ), INTENT(INOUT) :: externals

    CALL init_environ_density( cell, externals%density )

    RETURN

  END SUBROUTINE init_environ_externals_second

  SUBROUTINE update_environ_externals( externals )

    IMPLICIT NONE

    TYPE( environ_externals ), INTENT(INOUT) :: externals

    CALL density_of_functions( externals%number, externals%functions, externals%density )

    RETURN

  END SUBROUTINE update_environ_externals

  SUBROUTINE destroy_environ_externals( lflag, externals )

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_externals ), INTENT(INOUT) :: externals

    IF ( lflag ) CALL destroy_environ_functions( externals%number, externals%functions )

    CALL destroy_environ_density( externals%density )

    RETURN

  END SUBROUTINE destroy_environ_externals
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
!--------------------------------------------------------------------
END MODULE externals_utils
!--------------------------------------------------------------------
