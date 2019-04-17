.. Environ documentation example03 file, created by
   Matthew Truscott on Mon Apr 8 2019.

Example: Electrolyte Solutions
==============================

This example, like the previous two, models a 2-D metallic structure, this time, an Ag (100) slab, in
water solution. The
slab is charged and an electrolyte solution is added. This solution can be thought of as a solvent (whose
characteristics we have now implemented a number of times), and electrolyte ions, that are defined based off
a set of parameters. The parameters describing the ions are placed in the ENVIRON keyword

.. code-block:: fortran

   &ENVIRON
      env_electrolyte_ntyp = 2
      zion(1) = 1
      zion(2) = -1
      cion(1) = 1.0
      cion(2) = 1.0
      cionmax = 10.0

The env_electrolyte_ntyp parameter sets the number of ionic species that are to be defined, and for each species
the charge (zion) and the concentration (cion) are to be specified. Finally the maximum concentration can be
set (respresenting the maximum concentration of ions at any point). 

A number of electrolyte models can be implemented in Environ. This example iterates through the models based
the numerical method, which utilizes the Poisson-Boltzmann equation or some variant. It also includes some
analytical analogues that, rather than solve the Poisson-Boltzmann equation, add a correction to the potential.
By default, the Poisson-Boltzmann model is set. By specifying the parameter

.. code-block:: fortran

   &ENVIRON
      electrolyte_linearized = .true.

the model is changed to the linear equivalent. Refer to the publications [1]_, [2]_ for more details on these 
models. Some pbc correction (as explained in example 1) is necessary here, since calculations involving the
electrolyte require open boundary conditions. In the example, this is set to
parabolic, referring to the numerical approach, whereas a change to the gcs correction,

.. code-block:: fortran

   &ELECTROSTATIC
      pbc_correction = 'gcs'

will switch to the analytic approach to solving the relevant Poisson-Boltzmann equation.
Note that this does not include the size-modified equation. The page on diffuse-layer models will describe
each of the models and how they differ in more detail. In particular, the required parameters for each models
do vary slightly. For the linearized Poisson-Boltzmann model and the Poisson-Boltzmann model, the charge
and concentrations need to be specified. The additional max concentration (cionmax) parameter should be
included when choosing the modified Poisson-Boltzmann model. Notice that the model need not be specified
explicitly and is instead inferred by the parameters supplied to the input. The physical models all require
at least the concentrations and the charges of the ions. By setting the maximum concentration to be non-zero,
Environ will assume the user wishes to use the modified Poisson-Boltzmann model. 

.. note::

   The gcs correction (Gouy-Chapman Stern) is a 1-D analytical solution to the Poisson-Boltzmann equation
   and thus is only valid for 2 dimensional systems. The correction has been shown to produce close results
   to the numerical equivalents and should therefore be considered if the user is looking to save computational
   time. The parabolic correction (as seen in previous examples) is a general correction that works in any
   number of dimensions, however in order for the gcs correction to be valid, there should be no parabolic
   correction.

As in previous examples, the
solvent_mode parameter is set to 'full', due to a hole close to the Ag ions. A solvent-accessible but ion-free
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
.. [3] \I. Dabo et al., arXiv 0901.0096 (2008)
