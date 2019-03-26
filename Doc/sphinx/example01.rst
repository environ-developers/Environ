.. Environ documentation example01 file, created by
   Matthew Truscott on Tue Mar 26 2019.

Example 1
=========

This first example will demonstrate using the pw.x executable with Environ to calculate the solvation energy
and other related solvation quantities for a water molecule in a water solvent. A more in-depth explanation 
of the demonstration can be found in the README.

In general, execution of pw.x requires a correctly formatted input file with all the descriptive properties of
the desired system. With the Environ addon, an additional input file is required that contain simulation
parameters for the environment. Unlike the pw input file that can be arbitrarily named and then specified on
execution of the program, the Environ file must be named ‘environ.in’. To enable the environ addon on a pw
calculation, the --environ modifier should be added. Hence the command, ./pw.x --environ < example.in would 
feed in a pw input file named example.in into pw.x and run a serial calculation via whichever FORTRAN compiler 
Quantum ESPRESSO has been configured with. Since the --environ modifier added, pw now expects an environ.in 
file. Failure to do so will result in an error. 

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
are typically used on the development side for debugging. The option verbose=1 (or any higher integer value)
will output these in a file named environ.debug. The option verbose=2 will additionally output gaussian cube
files that contain the density objects of important quantities. Specifying this option will impact performance
due to the heavy I/O operations, and for the majority of simulations, a verbosity of zero ought to be
sufficient.

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

By design, pw assumes a simulation cell with 3D periodicity. This can be overcome by setting a simulation cell
with enough space and enabling the pbc_correction parameter (which in turn requires the env_electrostatic
parameter to be set as true). This correction is an approximation but should be sufficient for most simulation
runs. For more accurate results, one may want to refer to martyna-tuckerman (see Q-E input documentation under
‘assume-isolated’), which does require a larger simulation cell (and therefore more physical memory). For
simulations that require one to retain periodicity in 1 or 2 dimensions (for example, a diffuse-layer
simulation), the pbc_dim parameter can be set appropriately. 

The electrostatic calculation has its own accuracy that, for the most part, should be set manually in the input
file (since the value itself is reliant on the solver picked). 

On running a pw calculation with Environ while reading the output, a few differences are apparent. Firstly
before the SCF loop, Environ Module will be printed, followed by a brief description of the module, and some of
the settings (as defined in the input file). On each step of the SCF iteration, there are additional energy
corrections printed that are a result of the Environ module (the electrostatic embedding term and the correction
to one-el term). These corrections should be non-zero once the threshold specified in Environ’s input file is
met. Finally one can inspect the computation cost of the environ processes at the end of the file. 

In this example, the boundary is left empty. By default, Environ will use the SCCS model to define the
interface. This model uses the electrostatic density to define the separation between solute and solvent. 

A number of different solvers can be employed to calculate for electrostatic potential. Typically these default
to sensible values according to the setup, however, one can manually set these along with the auxiliary
parameter (in the event of say, polarization charge). 

