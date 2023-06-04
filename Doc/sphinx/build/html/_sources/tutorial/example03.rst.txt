.. Environ documentation example03 file.
   Created by Matthew Truscott on Mon Apr 8 2019.

Example: Isolated Systems
=========================

This example shows how to use pw.x to simulate isolated charge systems in vacuum or immersed in a continuum
dielectric constant. Note that this can be thought of as an alternative to Environ's own periodic correction
as seen in Example 1. Since this setting resides in the pw input file, both input files should be checked for
consistence. In the first example we removed all periodicity from the system with a correction that results
in a fairly good approximation. 

.. literalinclude:: example01.in
   :language: fortran
   :lines: 10-12

In this example, we demonstrate a system without periodicity, thus removing the parameters stated before and
keeping to the defaults. We also demonstrate an alternative to the pbc correction, which is the
martyna-tuckerman [1]_ correction. This correction is seen as a more accurate method, at the cost of a larger
required cellsize. For reliable results, a cell length of twice the molecule length is recommended, in any
direction. This can be added into the pw input file, under the SYSTEM keyword.

.. code-block:: fortran

   &SYSTEM
      assume_isolated = 'martyna-tuckerman'

The tolerance is different from previous examples. This value is chosen to facilitate the self-consistent
process and will directly affect the speed and accuracy of the convergence. One should set this according to
the type of simulation being run, these values in the examples serving as guidelines for that decision.

Clearly for larger systems, the parabolic correction is a more feasible choice for imposing isolation.
This is a real-space quadratic correction of the electrostatic potential [2]_, which has been shown to
provide energy accuracy improvements of almost 1 order of magnitude in comparison with the alternative PCC
and GCC schemes for large cell sizes. This correction scheme implemented in Environ is limited to 0D and 2D
systems (which account for the majority of systems that Environ typically works with).

.. [1] G. J. Martyna and M. E. Tuckerman, J. Chem. Phys. 110, 2810 (1999)
.. [2] \I. Dabo et al. Phys. Rev. B 77, 115139 (2008)
