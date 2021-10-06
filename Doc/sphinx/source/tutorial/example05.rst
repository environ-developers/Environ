.. Environ documentation example03 file.
   Created by Matthew Truscott on Mon Apr 8 2019.

.. _ex05:

Example: Charge Distributions
=============================

This example follows the previous, now adding a fixed, planar and smooth distribution of charge (specifically a 
helmholtz plane). The environ input file for this kind of system is a little more complicated. 

First is the use of the keyword solvent_mode = 'full'. For this calculation, this choice is mandatory due to
the Pt valence density, which results in a hole close to the ionic position, due to the missing core electrons.
This option adds an additional charge density (gaussians) at the location of the ions, that represent the core
electrons and the nuclei. This results in a better description of the interface between vacuum and solvent.
This choice is useful for particular types of atoms where this effect can occur, for example halogens and
transition metals.

.. code-block:: fortran

   &BOUNDARY
      solvent_mode = 'full'

To represent the helmholtz plane, it is necessary to specify the number of external charges in the ENVIRON
keyword and then specify the the charges themselves in a block labelled EXTERNAL_CHARGES (at the end of the
input file). Notice how this structure is similar to how the pw input files handle atomic positions.

.. code-block:: fortran

   EXTERNAL_CHARGES (bohr)
   -0.5 0. 0. 25.697 1.0 2 3
   -0.5 0. 0. -10.303 1.0 2 3

Each line specifies a charge, and possesses 7 values separated by spaces. In order these describe the charge
(in units of electron charge), x position, y position, z position, spread, dimensionality, and axis. In this
example we wish to define charge planes parallel to the Pt slab. Since we are defining two-dimensional objects
that lie perpendicular to the z-axis, the x and y positions are irrelevent and are thus set to 0.0. The charge
functions are represented by gaussians, and thus the parameters are simply position and spread. Finally the
charges can represent basic shapes depending on the dimensionality defined. For a slab, we specify a
dimensionality of 2. For the direction of the slab, we specify 3, which represents the z-axis. This convention
is consistent with the boundary correction options (explained in the previous example) that define the system
to be periodic in two-dimensions as opposed to the three by default.

We set this system up so that there exists two countercharge planes on either side of the Pt slab. Since the
system is periodic, imposing symmetry in this manner allows the simulation to perform better and is generally
advised. Since there are now two planes of half charge, the resultant electric field is half what it would have
otherwise been. For detailed information on these planar charge setups, see [1]_, [2]_, or :ref:`expq`. 

.. [1] C. L. Fu and K. M. Ho, Phys. Rev. Lett. 63, 1617 (1989)
.. [2] F. Nattino et al., J. Chem. Phys. 041722 (2019)

