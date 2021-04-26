!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE cell_types
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP, tpi
    !
    USE fft_types, ONLY: fft_type_descriptor, fft_type_init, fft_type_deallocate
    !
    USE stick_base, ONLY: sticks_map, sticks_map_deallocate
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !------------------------------------------------------------------------------------
    !>
    !! A simulation cell
    !!
    !! Notes:
    !!
    !! 1. ir_end can be different from nnr due to FFT-grid optimization yielding
    !!    additional, unphysical grid points
    !!
    !! 2. cell corners are utilized in minimum_image()
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_cell
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: update = .FALSE.
        LOGICAL :: cubic = .FALSE.
        REAL(DP) :: omega, domega ! volume quantities
        REAL(DP) :: origin(3)
        REAL(DP) :: at(3, 3) ! real-space lattice vectors
        REAL(DP) :: bg(3, 3) ! reciprocal lattice vectors
        REAL(DP) :: corners(3, 8)
        !
        !--------------------------------------------------------------------------------
        ! Units needed to scale real and reciprocal space cells
        ! #TODO redesign Environ units
        !
        REAL(DP) :: alat ! lattice parameter
        REAL(DP) :: tpiba, tpiba2 ! tpi = 2*pi
        !
        !--------------------------------------------------------------------------------
        ! Properties of the grid
        !
        TYPE(fft_type_descriptor) :: dfft
        INTEGER :: ntot ! total number of grid points
        INTEGER :: nnr ! number of grid points allocated in every processor
        INTEGER :: ir_end ! actual number grid points accessed by each processor
        INTEGER :: j0, k0 ! starting indexes of processor-specific boxes of grid points
        REAL(DP) :: in1, in2, in3 ! inverse number of grid points
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_cell
    !------------------------------------------------------------------------------------
    !>
    !! The mapping between system and environment simulation cells
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_mapping
        !--------------------------------------------------------------------------------
        !
        INTEGER :: nrep(3)
        ! number of system cells in environment cell along a_i = 2 * nrep_i + 1
        !
        TYPE(environ_cell), POINTER :: small ! system cell
        TYPE(environ_cell), POINTER :: large ! environment cell
        !
        INTEGER, ALLOCATABLE :: map(:)
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_mapping
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE cell_types
!----------------------------------------------------------------------------------------
