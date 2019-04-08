.. Environ documentation example03 file, created by
   Matthew Truscott on Mon Apr 8 2019.

Example: Two-Dimensional Systems
================================

This example shows how to use pw.x to model two-dimensional periodic systems in contact with a continuum solvent.
This is particurly useful for diffuse layer systems and other surface investigations. For this calculation, a
Pt (111) slab is chosen. By running the bash file and inspecting the environ input files, one will notice that
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
