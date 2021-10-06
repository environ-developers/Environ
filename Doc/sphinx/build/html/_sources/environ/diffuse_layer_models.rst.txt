.. Environ documentation diffuse layer models file.
   Created by Matthew Truscott on Mon Apr 8 2019. 
   Contains general description and comparison of diffuse layer models.
   Updated by Edan Bainglass on Mon Oct 4 2021.

Diffuse Layer Models
====================

This page presents a more in-depth description of the diffuse layer models implemented in Environ, along with
all revelant input parameters.

.. _expq:

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

For an example of plane countercharges, see :ref:`ex05`. 

For models with charge distribution, the free energy functional can be written as

.. math::

   F^{\text{PC}}[\rho(\mathbf{r}), \phi(\mathbf{r}] = \int\left[-\frac{\epsilon(\mathbf{r})}{8\pi}\lvert\nabla\phi(\mathbf{r})\rvert^2 + \rho(\mathbf{r})\phi(\mathbf{r}) + \rho^{\text{ions}}(\mathbf{r})\phi(\mathbf{r})\right]d\mathbf{r}

where :math:`\rho(\mathbf{r})` is the total charge density of the solute, :math:`\phi(\mathbf{r})` is the
electrostatic potential, and :math:`\rho^{\text{ions}}(\mathbf{r})` is the external charge density that
mimics the counterion accumulation.

Poisson-Boltzmann model
-----------------------

The Poisson-Boltzmann model accounts for the chemical potential and the entropy of the ions in solution.
Adding these consequently explicit concentration-dependent terms results in a different free energy functional
than before,

.. math::

   F^{\text{PC}}[\rho(\mathbf{r}), \phi(\mathbf{r}, \{c_i(\mathbf{r}\}] = \int\left[-\frac{\epsilon(\mathbf{r})}{8\pi}\lvert\nabla\phi(\mathbf{r})\rvert^2 + \rho(\mathbf{r})\phi(\mathbf{r}) + \rho^{\text{ions}}(\mathbf{r})\phi(\mathbf{r})\right.
   
   \left.-\sum^{\text{p}}_{i=1}\mu_i(c_i(\mathbf{r})-c_i^0)-T(s[\{c_i(\mathbf{r})\}]-s[\{c_i^0\}])\right]d\mathbf{r}.

Here, :math:`\mu_i` is the chemical potential of the ith electrolyte species and T is the temperature, that
is included as part of the entropy term.

We assume that these electrolytes have point-charge, and there is ideal mixing, leading us to the following
entropy expression,

.. math::

   s[\{c_i\}] = -k_B\sum^{\text{p}}_{i=1}c_i(\mathbf{r})\ln\frac{c_i(\mathbf{r})}{\gamma(\mathbf{r})}

where :math:`\gamma(\mathbf{r})` is the exclusion function and sets the boundary between the electrolyte
and solute regions.

The functional above is
minimized with respect to the concentrations in order to find the equilibrium for electrolyte concentrations.
This leads to the well-known Poisson-Boltzmann equation, which can be numerically solved using Newton's
iterative algorithm and a preconditioned conjugate gradient algorithm, which works with linear equations. 

There are times when it is beneficial to approximate the Poisson-Boltzmann equation, which is non-linear and
thus non-trivial to solve. For :math:`z_i\phi(\mathbf{r}) \ll k_BT`, that is, when the electrostatic potential
is low, or in the high-temperature limit, one can approximate the Poisson-Boltzmann by its linear form, which
is derived by Taylor expanding the appropriate terms in the full Poisson-Boltzmann equation. The result is
a less computationally expensive approach to representing the electrolyte, while still maintaining the same
parameter set. 

In the case where we deal with linear slabs with two-dimensional periodicity, a different approximation can
be used. The idea here is that the one-dimensional Poisson-Boltzmann equation can be analytically solved and
thus this result can be applied as a PBC (periodic boundary condition) correction. This approach is also
known as the Gouy-Chapman Stern model for an electrolyte, and compares well to the numerical approach for
our examples. 

Modified Poisson-Boltzmann model
--------------------------------

The Poisson-Boltzmann can be improved by dropping the assumpion of point-like ions. This assumption leads to
an overestimate of the electrolyte countercharge accumulation at electrode surfaces, and thus, by accounting
for the steric repulsion between ions, which opposes the electrostatic attraction between the counterions and
the electrode surface, one can improve on the full Poisson-Boltzmann model. The size-modified
Poisson-Boltzmann (MPB) can be derived from the same free energy functional as before, only with a modified
entropy expression

.. math::

   s[\{c_i\}] = -k_B\sum^{\text{p}}_{i=1}c_i(\mathbf{r})\ln\frac{c_i(\mathbf{r})}{c_{\text{max}}\gamma(\mathbf{r})}

   -k_B\left(c_{\text{max}}\gamma(\mathbf{r} - \sum^{\text{p}}_{i=1}c_i(\mathbf{r})\right)\ln\left(1 - \sum^{\text{p}}_{i=1}\frac{c_i(\mathbf{r})}{c_{\text{max}}\gamma(\mathbf{r})}\right).

The idea is to essentially impose a space dependent maximum ionic concentration, and the result is a better
representation of the electrolyte, verified by a comparison to experimental differential capacitance.

Additional Interactions
-----------------------

Environ
-------

Environ has implemented all of the above models in a modular way that allows one to mix and match models and
correction methods where reasonable. 
   




