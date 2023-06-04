.. Environ documentation example01 file.
   Created by Matthew Truscott on Tue Mar 26 2019.

Example: Solvation Energy (SCCS)
================================

This first example will demonstrate using the pw.x executable with Environ to calculate the solvation energy
and other related solvation quantities for a water molecule in a water solvent. A more in-depth explanation 
of the demonstration can be found in the README. The overall idea is as follows: solvation energy is the
difference in energy between a solute in vacuum and in solvent. By default, the pw.x program calculates
the energy of a system in vacuum. The Environ module allows us to account for additional effects or
correction factors. It also allows us to add some non-vacuum environment, in this case, a water solvent.
Hence by running pw.x twice with different Environ parameters, one can calculate the energy of the system
in vacuum and in water, thus obtaining the solvation energy. 

In general, execution of pw.x requires a correctly formatted input file with all the descriptive properties of
the desired system. With the Environ addon, an additional input file is required that contains simulation
parameters for the environment. Unlike the pw input file that can be arbitrarily named and then specified on
execution of the program, the Environ file must be named ‘environ.in’. To enable the environ addon on a pw
calculation, the ``--environ`` modifier should be added. Hence the command, ``./pw.x --environ < example.in`` 
would 
feed in a pw input file named example.in into pw.x and run a serial calculation via whichever FORTRAN compiler 
Quantum ESPRESSO has been configured with (as with any program, call ``pw.x`` from its actual position, or add
it to the PATH). Since the ``--environ`` modifier is added, ``pw.x`` now expects 
an environ.in file. Failure to do so will result in an error. 

The bash script run_example.sh generates all the necessary files described above and executes multiple pw 
calculations successively. Due to the constraint of the environ.in name, a change in environment should require
files to be renamed on the fly (something that is taken care of here in the bash script). Alternatively one
can set up multiple runs to execute in separate directories thus avoiding this issue entirely.

The general format of the environ.in file can be seen starting line 154 of the bash script. Alternatively one
can run the bash script and view the generated environ files in the results folder (these have the environ
prefix and the .in suffix). Refer to the input documentation (which can be viewed on the browser from the local
html file, Doc/INPUT_Environ.html) for in-depth information on any of the variables as well as the expected
structure. The general structure splits parameters into a number of keywords, &ENVIRON, &BOUNDARY, and
&ELECTROSTATIC. This file has fortran-like formatting and as such, '!' can prepend comments.

.. literalinclude:: example01.in
   :language: fortran
   :lines: 1-3

Environ is designed to be modular and as a result, many of the parameters are interchangeable. Each
parameter belongs to a corresponding keyword and therefore should be placed accordingly. The parameters
used in the program are also helpfully described in the bash file. By default, the verbosity parameter (verbose)
is set to zero, which limits output to the standard pw output location (which is to the terminal, and therefore
is often piped into an output file). 

Some helpful environment specific quantities can also be printed out. These
are typically used on the development side for debugging. The option ``verbose=1`` (or any higher integer value)
will output these in a file named environ.debug. The option ``verbose=2`` will additionally output gaussian cube
files that contain the density objects of important quantities. Specifying this option will impact performance
due to the heavy I/O operations, and for the majority of simulations, a verbosity of zero ought to be
sufficient, which is the default.

.. literalinclude:: example01.in
   :language: fortran
   :lines: 2,4

The environ threshold specifies when to start environ contributions to the SCF iterative step. The default
value is suitable primarily for small systems. Note that since these input files follow fortran convention,
double precision is notated as above. 

.. literalinclude:: example01.in
   :language: fortran
   :lines: 2,5

The environ type simplifies the requirement of some of the physical parameters (such as pressure, surface
tension, and static permittivity), which arise depending on the interface model applied. One can manually add
these, but these have been tuned to optimal values, and therefore it is generally preferable to use the preset
values by setting this parameter accordingly (vacuum, water, water-anions, water-cations). 

.. note::

   This option is equivalent to

   .. code-block:: fortran

      &ENVIRON
         environ_type = 'input'
         env_surface_tension = 0
         env_pressure = 0
         env_static_permittivity = 1

.. literalinclude:: example01.in
   :language: fortran
   :lines: 10-12
         
By design, pw assumes a simulation cell with 3D periodicity. This can be overcome by setting a simulation cell
with enough space and enabling the pbc_correction parameter (which in turn requires the env_electrostatic
parameter to be set as true). This correction is an approximation but should be sufficient for most simulation
runs. For more accurate results, one may want to refer to martyna-tuckerman (see Q-E input documentation under
‘assume-isolated’), which does require a larger simulation cell (and therefore more physical memory). For
simulations that require one to retain periodicity in 1 or 2 dimensions (for example, a diffuse-layer
simulation), the pbc_dim parameter can be set appropriately. More details on this quadratic correction method
can be found in the publication [1]_. 

.. literalinclude:: example01.in
   :language: fortran
   :lines: 10,13-14

The electrostatic calculation has its own accuracy that, for the most part, should be set manually in the input
file (since the value itself is reliant on the solver picked). 

On running a pw calculation with Environ while reading the output, a few differences are apparent. Firstly
before the SCF loop, Environ Module will be printed, followed by a brief description of the module, and some of
the settings (as defined in the input file). On each step of the SCF iteration, there are additional energy
corrections printed that are a result of the Environ module (the electrostatic embedding term and the correction
to one-el term). These corrections should be non-zero once the threshold specified in Environ’s input file is
met. Finally one can inspect the computation cost of the environ processes at the end of the file. 

.. literalinclude:: example01.in
   :language: fortran
   :lines: 8-9

In this example, the boundary is left empty. By default, Environ will use the SCCS model to define the
interface. This model uses the electrostatic density to define the separation between solute and solvent. 

.. literalinclude:: example01.in
   :language: fortran
   :lines: 10,15-16

A number of different solvers can be employed to calculate for electrostatic potential. Typically these default
to sensible values according to the setup, however, one can manually set these along with the auxiliary
parameter (in the event of say, polarization charge).

For small investigations such as this example, where only 3 simulations are run, it is sufficient to process
the result directly from the command line. When scaling up to a higher number of simulations, bash or
python scripting can be exploited to automate the process of calculating the solvation energy. Suppose the
user is in the directory containing the newly created PW output files (in the results folder), the output files
of Quantum ESPRESSO are designed for convenient extraction of commonly used results via grep. Calling
``grep ! *.out`` will return the total energy values. Since this is a PW relax calculation, the energy is
calculated multiple times, while updating the positions of the ions depending on the calculated forces in the
system. One should expect to get a result like:

::
     
   h2o_vacuum_direct.out:!    total energy              =     -34.26804813 Ry
   h2o_vacuum_direct.out:!    total energy              =     -34.26815510 Ry
   h2o_vacuum_direct.out:!    total energy              =     -34.26825783 Ry
   h2o_vacuum_direct.out:!    total energy              =     -34.26826238 Ry
   h2o_water_cg.out:!    total energy              =     -34.28845991 Ry
   h2o_water_cg.out:!    total energy              =     -34.28851745 Ry
   h2o_water_cg.out:!    total energy              =     -34.28858148 Ry
   h2o_water_iterative.out:!    total energy              =     -34.28845478 Ry
   h2o_water_iterative.out:!    total energy              =     -34.28851656 Ry
   h2o_water_iterative.out:!    total energy              =     -34.28857886 Ry

Note that the values may not match exactly due to the compiler used. It is useful to bear in mind when
analysing energies in Rydberg, differences of E-2 are typical for solvation energies, whereas differences of 
less than 3E-3 are less than chemical accuracy. We see that the final energy of the system depends on the 
solver used, but the difference is 3E-5, well below chemical accuracy and therefore not of much concern for
regular applications. Hence the solvation energy can be calculated (as the difference between the solvation
and vacuum energies) as -6.374kcal/mol (or 2.03E-2Ry) via the conjugative gradient solver, and -6.373kcal/mol
(or 2.03E-2Ry) via the iterative solver. This compares closely to the experimental solvation energy of
water (-6.3kcal/mol [2]_).

If the user is just interested in the final energy in a relax calculation, a command like ``grep Final *.out``
will achieve this. This example conveniently automates the calculation of these values via bash (along with
other helpful energy values), and prints them out into ``results.txt`` for reference. 

.. [1] I. Dabo et al., Phys. Rev. B 77, 115139 (2008)
.. [2] A. V. Marenich et al., Minnesota Solvation Database - version 2012; University of Minnesota: Minneapolis
