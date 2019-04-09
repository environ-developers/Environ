.. Environ documentation example03 file, created by
   Matthew Truscott on Mon Apr 8 2019.

Example: Electrolyte Solutions
==============================

This example, like the previous two, models a Pt (111) slab with a CO molecule adsorbed on the surface. The
slab is charged and an electrolyte solution is added. This solution can be thought of as a solvent (whose
characteristics we have now implemented a number of times), and electrolyte ions, that are defined based off
a set of parameters. The parameters describing the ions are placed in the ENVIRON keyword

.. code-block:: fortran

   &ENVIRON
      env_electrolyte_ntyp = 2
      zion(1) = 1
      zion(2) = -1
      cion(1) = 5.0
      cion(2) = 5.0
      cionmax = 15.0

The env_electrolyte_ntyp parameter sets the number of ionic species that are to be defined, and for each species
the charge (zion) and the concentration (cion) are to be specified. Finally the maximum concentration can be
set (respresenting the maximum concentration of ions at any point). 

A number of electrolyte models can be implemented in Environ. In this example, the linearized Poisson-Boltzmann
model is implemented. By default, the Poisson-Boltzmann model is set. By specifying the parameter

.. code-block:: fortran

   &ENVIRON
      electrolyte_linearized = .true.

the model is changed to the linear equivalent. Refer to the publications [1]_, [2]_ for more details on these 
models. A pbc correction (as explained in example 1) is strongly recommended here. As in previous examples, the
solvent_mode parameter is set to 'full', due to a hole close to the Pt ions. A solvent-accessible but ion-free
region in the proximity of the surface (Stern layer) is modeled by introducing an ion-exclusion function as
described here [3]_. 

The interface function respresenting the boundary between the quantum mechanical region and the electrolyte
embedding environment is set in the BOUNDARY keyword

.. code-block:: fortran

   &BOUNDARY
      electrolyte_mode = 'system'
      electrolyte_distance = 10.0
      electrolyte_spread = 0.5

Here, the electrolyte_mode parameter is set to system, corresponding to a simple analytical function. 

.. [1] F. Nattino et al., J. Chem. Phys. 150, 041722 (2019)
.. [2] G. Fisicaro et al., J. Chem. Phys. 144, 014103 (2016)
.. [3] I. Dabo et al., arXiv 0901.0096 (2008)
