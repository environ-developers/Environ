.. Environ documentation diffuse layer models file, created by
   Matthew Truscott on Mon Apr 8 2019. Contains general description
   and comparison of diffuse layer models.

Diffuse Layer Models
====================

This page presents a more in-depth description of the diffuse layer models implemented in Environ, along with
all revelant input parameters.

Explicit Charge
---------------

Environ has the ability to represent external charges of various types. As with most of the examples on this
page, we will assume a setup consisting of some noble metal slab (i.e. a two-dimensional system) that can
possess some surface charge. To compensate we then consider the addition of countercharge, either explicitly or
via an electrolyte solution. Explicit charges come in a variety of implemented forms depending on the specified
dimensionality. Suppose we want a point charge electron at some arbitrary position. One could enter the 
following into the Environ input file,

.. code-block:: fortran

   EXTERNAL_CHARGES (bohr)
   -1 0. 0. 0. 1.0 0 0

Remember to update the env_external_charges parameter appropriately (here we could have something like

.. code-block:: fortran

   &ENVIRON
      env_external_charges = 1

This places a point-charge (since the dimensionality is set to 0), at the origin. The charge is given a negative
value whose magnitude equals that of an electron. The spread is set to 1.0. The dimensionality is given as 0, 
and the axis setting is irrelevent but set to 0 anyway. The definition of the spread is
based off the following 1D definition for the gaussian, 

.. math::

   ae^{\frac{(x-b)^2}{c^2}}

Note the lack of the factor of two that some definitions use. In three dimensions, the equation is effectively
the same, only radial around the 3-dimensional point specified in the input. The gaussian is normalized to the
specified charge. This convention has a lot of transferability, and by simply changing the dimensionality of the
input, can be extended to line charges, and planar charges. For these objects, the axis is relevant. For a line 
charge, the chosen axis corresponds to a the axis parallel to the line charge, and for a plane, it corresponds 
to the axis normal to the plane. Note that this definition clearly does not account for all possible definitions
of line and plane charge distributions, but this can be achieved instead by changing the atomic positions of 
the slab.

For an example of plane countercharges, see ex05_


