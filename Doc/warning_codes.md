# Environ Warning Codes

1001. Ignored card

      An unexpected input file card has been ignored. Check input file

</br>

1002. Strange lattice parameter

      Host lattice parameter unusually small

</br>

1003. Missing FFT-grid information

      Environ is missing some information required to establish both the FFT-grid and its parallelization scheme

      - If only a G-vector cutoff value (`gcutm`) is provided, Environ generates its FFT-grid and parallelization stick map
      - Alternatively, a user may provide `env_ecut` in the Environ input (`gcutm = env_ecut / pi**2`)
        - Note that `env_ecut` is in Rydberg and is equivalent, for example, to Quantum ESPRESSO's `ecutrho`
      - If either are provided together with an FFT-grid (`nr`), Environ sets its internal FFT-grid to `nr` and uses `gcutm` for the stick map
      - If only `nr` is provided, Environ will derive a conservative overestimate of `gcutm` from `nr` and the lattice (`at`)
        - Note that this overestimate may introduce errors. Caution is advised

</br>

1004. Left-handed lattice vectors

</br>

1005. Wrong integral of ERFC function

      Error in numerical integral w.r.t analytic solution exceeds threshold

      Reasons:

      - the grid is too coarse with respect to the spread of the ERFC (or vice-versa)
        - consider increasing grid spacing and/or ERFC function spread
      - the ERFC function spills over the cell boundary
        - consider increasing cell size
