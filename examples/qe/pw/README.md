# Environ PW Examples

sccs: pw.x

    Calculation of total solvation free energy of an isolated molecule
    including electrostatic/cavitation/PV contributions. The solvation
    cavity is defined according to the revised-SCCS.

sscs: pw.x

    Same as example01, but the solvation cavity is constructed from
    atomic-centered interlocking spheres ("soft-spheres").

pbc: pw.x

    Use of different periodic boundary correction schemes for charged
    isolated systems in vacuum and in solution: Martyna-Tuckerman or
    Point-Counter-Charge (parabolic).

slab: pw.x

    Use of parabolic correction for periodic boundary conditions in
    neutral and charged 2D systems (slab).

helmholtz: pw.x

    Use of external classical charge distribution (planar) to model
    Helmholtz layer above a charged 2D system.

helmholtz_linpb: pw.x

    Same as example05, but the electrolyte counter charge layer is
    modeled according to the linearized Poisson-Boltzmann model.

mott_schottky: pw.x

    Use of a Mott-Schottky counter charge layer to model a charged
    slab embedded in a semiconducting region with different defect
    densities.

solvent_aware: pw.x

    Use of the solvent-aware interface to prevent the dielectric to
    penetrate regions that should remain solvent-free.

field_aware: pw.x

    Use of the field-aware interface to automatically adjust the
    interface for charged solutes.
