.. Environ documentation example03 file.
   Created by Matthew Truscott on Mon Apr 8 2019.

Example: Two-Dimensional Systems
================================

This example shows how to use pw.x to model two-dimensional periodic systems in contact with a continuum solvent.
This is particurly useful for diffuse layer systems and other surface investigations. For this calculation, a
Pt (111) slab is chosen with a CO molecule adsorbed on the surface. 
By running the bash file and inspecting the environ input files, one will notice that
most of the input remains unchanged. The SCCS (self-consistent charge solvation) model is chosen for the solvent
by default and the solvent type is set manually, for example for vacuum,

.. code-block:: fortran
   
   &ENVIRON
      environ_type = 'input'
      env_static_permittivity = 1
      env_surface_tension = 0.D0
      env_pressure = 0.D0

The change is in the pbc_dim parameter. This needs to be set accordingly and the axis should be subsequently
chosen. This axis parameter determines the direction of the periodicity (that is the plane normal to the axis
will be parallel to the surface or slab). This parameter expects an integer value, 1 for x, 2 for y, and 3 for
z.

.. code-block:: fortran

   &ELECTROSTATIC
      pbc_dim = 2
      pbc_axis = 3

The tolerance is also different from previous examples. This parameter is typically set in order to optimize
the ability for the self-consistent process to converge, and should be picked appropriately depending on the
system. These examples serve as good guidelines for each type of simulation.

Along with the electrostatic solvation energy, which calculated similarly to Example 1, where one takes the
difference between the total energy in vacuum with total energy in electrolyte solvent, this example
calculates the Fermi energy both in vacuum and in solution. Unlike other energy terms, the printed total
Fermi energy is decomposed into terms, one from QE and a correction (delta energy) from Environ. Hence, the
user should sum these terms during post-processing in order to get the Fermi energy of the system (this
is done for the user in this example by the bash script).
