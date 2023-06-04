.. Environ documentation api objects file.
   Created by Edan Bainglass on Sun Jul 3 2022.

Objects
=======

Environ exposes the following class instances to provide host programs
with greater flexibility of calculation. For more information on each
class, please reference the provided links.

----

``setup``
#########

The setup class (`environ_setup`_) contains calculation flags,
simulation cells, and the numerical solvers and cores.

``main``
########

The main class (`environ_main`_) contains the physical quantities
of Environ, as well as calculated potentials and energies.

``calc``
########

The calculator class (`environ_calculator`_) contains the primary
calculators to compute energies, potentials, and forces.

``clean``
#########

The destructor class (`environ_destructor`_) contains methods to
destroy any objects created in the process of a calculation.

.. _environ_setup: https://github.com/environ-developers/Environ/blob/master/src/setup.f90
.. _environ_main: https://github.com/environ-developers/Environ/blob/master/src/main.f90
.. _environ_calculator: https://github.com/environ-developers/Environ/blob/master/src/calculator.f90
.. _environ_destructor: https://github.com/environ-developers/Environ/blob/master/src/destructor.f90
