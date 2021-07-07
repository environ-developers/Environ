# Environ Utilities

    This library implements various basic tasks such as timing, tracing,
    optimized memory accesses and an abstraction layer for the MPI subroutines

    Library borrowed from QE 6.7 (UtilXlib) and simplified for Environ purposes

### Modifications

---

- Applied formatting standards
- Discarded unused files
- Renamed files for consistency
- Renamed modules/routines using env_ prefix to force uniqueness
- Packaged all sorting routines in sorting.f90
- Packaged all string operations in char_ops.f90
- Packaged input/output and error routines in io.f90
- Moved kinds into kinds.f90

- Merged GPU/non-GPU functionality into single routines in base_mp
  - mp routines reflect this merge
  - NOTE: GPU is not yet tested in Environ and is experimental in QE
  - NOTE: base_mp routines MUST be placed outside of a module to avoid compiler type/rank errors
