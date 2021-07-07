# Environ FFTs

    Implements real space grid parallelization of FFTs.

    Library borrowed from QE 6.7 (FFTXlib) and simplified for Environ purposes.

### Modifications

---

- Applied formatting standards
- Discarded unused files
- Renamed files for consistency
- Renamed modules/routines using env\_ prefix to force uniqueness
- Now only handling potentials and charge densities ('Rho')
- fft_interfaces is now fft_main

  - Primary driver, containing interfaces and fw/inv fft routines

- Split fft_parallel into GPU and non-GPU modules

  - Merged 2d functionality into respective modules

- Merged fft_error into utils/errors.f90
