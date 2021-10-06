.. Environ documentation example02 file, created by
   Matthew Truscott on Fri Mar 29 2019.

Example: Solvation Energy (SSCS)
================================

In the previous example, one detail that was not explained was the choice of interface model.
As stated previously, Environ is designed to be as modular as is reasonable, which means the user is able
to define settings that enable different models, solvers, and corrections to different subsystems. Previously
the interface model, that defines the separation of the quantum mechanical system (set in the PW input file),
with the solvent environment, was set to SCCS. The SCCS model sets a smooth boundary according to the electronic
density at any particular point. This function is in fact a piece-wise function that is chosen to be easily
differentiable [1]_.

An alternative model that Environ has implemented is the soft-sphere (SSCS) model [2]_. This model sets
a smooth boundary according to the ionic positions. Interlocking spheres are taken around each ion, whose
radii are by default consistent with the UFF force fields [3]_ (note that in the paper, there is an exception 
for nitrogen, which is based off Bondi's radius for nitrogen [4]_, an option which is not a default in Environ)
in a nature similar to PCM, only these are smooth functions (spherical error functions).

This example repeats the previous setup, instead with a soft-spheres implementation. Since this is a boundary
parameter, one can see on execution of the bash script that the environ input files here differ, now possessing
parameters related to the &BOUNDARY keyword.

.. literalinclude:: example02.in
   :language: fortran
   :lines: 9, 11

Simply by setting the solvent_mode parameter one can switch from SCCS to soft-sphere model, environ_type presets
can take care of the required (tuned) parameters associated with the model. Note that the solver (and auxiliary)
parameters are not required for the SSCS.

Calculating the solvation energy is exactly the same procedure as before, and results in a value of
-5.786kcal/mol, within chemical accuracy of the experimental results, as seen in the first example.

.. [1] O. Andreussi, I. Dabo, N. Marzari, J. Chem. Phys. 136, 064102 (2012) 
.. [2] G. Fisicaro et al., J. Chem. Comput. 2017, 13, 8, 3829-3845
.. [3] A. K. Rappe et al., J. Am. Chem. Soc. (1992), 114 (25) 114 10024-35
.. [4] A. Bondi, J. Phys. Chem. (1964), 68 (3) 68 441-51
