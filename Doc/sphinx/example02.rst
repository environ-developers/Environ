.. Environ documentation example02 file, created by
   Matthew Truscott on Fri Mar 29 2019.

Example: Solvation Energy (soft-sphere)
=======================================

In the previous example, one detail that was not explained was the choice of interface model.
As stated previously, Environ is designed to be as modular as is reasonable, which means the user is able
to define settings that enable different models, solvers, and corrections to different subsystems. Previously
the interface model, that defines the separation of the quantum mechanical system (set in the PW input file),
with the solvent environment, was set to SCCS. The SCCS model sets a smooth boundary according to the electronic
density at any particular point. This function is in fact a piece-wise function that is chosen to be easily
differentiable [1]_.

This example repeats the previous setup, instead with a soft-spheres implementation. Since is a boundary
parameter, one can see on execution of the bash script that the environ input files here differ, now possessing
parameters related to the &BOUNDARY keyword.

Simply by setting the solvent_mode parameter one can switch from SCCS to soft-sphere model, environ_type presets
can take care of the required (tuned) parameters associated with the model. Note that the solver (and auxiliary)
parameters are not required

.. [1] O. Andreussi, I. Dabo, N. Marzari, J. Chem. Phys. 136, 064102 (2012) 
