# Environ Examples

    Example QE+Environ runs

example01: pw.x

    Calculation of total solvation free energy of an isolated molecule
    including electrostatic/cavitation/PV contributions. The solvation
    cavity is defined according to the revised-SCCS.

example02: pw.x

    Same as example01, but the solvation cavity is constructed from
    atomic-centered interlocking spheres ("soft-spheres").

example03: pw.x

    Use of different periodic boundary correction schemes for charged
    isolated systems in vacuum and in solution: Martyna-Tuckerman or
    Point-Counter-Charge (parabolic).

example04: pw.x

    Use of parabolic correction for periodic boundary conditions in
    neutral and charged 2D systems (slab).

example05: pw.x

    Use of external classical charge distribution (planar) to model
    Helmholtz layer above a charged 2D system.

example06: pw.x

    Same as example05, but the electrolyte counter charge layer is
    modeled according to the linearized Poisson-Boltzmann model.

example07: pw.x

    Use of the solvent-aware interface to prevent the dielectric to
    penetrate regions that should remain solvent-free.

example08: pw.x turbo_lanczos.x

    Calculation of the optical spectrum of a molecule in vacuum and
    in solution (only dielectric)

example09: pw.x turbo_davidson.x

    Calculation of the optical spectrum of a molecule in vacuum and
    in solution (only dielectric)

example10: cp.x (only for QE-5.3.0 and later versions)

    Calculation of the total solvation free energy of an isolated
    molecule including electrostatic/cavitation/PV contributions

example11: neb.x

    Calculation of the minimum energy path (MEP) of a simple 
    activated reaction i.e. the collinear proton transfer reaction
