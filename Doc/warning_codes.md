# Environ Warning Codes

1001. Ignored card

      An unexpected input file card has been ignored. Check input file

1002. No free unit

      IO module was unable to find an empty unit for an IO operation

1003. Strange lattice parameter

      Host lattice parameter unusually small

1004. Left-handed lattice vectors

1005. Wrong integral of ERFC function

      Error in numerical integral w.r.t analytic solution exceeds threshold

      Reasons:

      - the grid is too coarse with respect to the spread of the ERFC (or vice-versa)
        - consider increasing grid spacing and/or ERFC function spread
      - the ERFC function spills over the cell boundary
        - consider increasing cell size
