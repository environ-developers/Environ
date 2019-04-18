.. Environ FAQ file, created by Matthew
   Truscott on Thu Apr 18 2019.

FAQ
===

This section attempts to address Frequently Asked Questions regarding Environ. If your query does not
relate to any of these sections, check out the Q&A section on the website, which has instructions on joining
the mailing group where you can chat with one of the developers. 

Converge Issues
---------------

There are physical and/or numerical reasons for why a calculation may explode. The latter is easier to solve,
by tuning some of the parameters in the input file. The former is a little more tricky to handle. 

Roughly speaking, in a calculation with Environ, one needs to include some extra terms to the Kohn-Sham
potential of the system in vacuum. These terms depend on the electronic density itself, in particular due to
the fact that the interface between the system and the environment is defined in terms of the electronic
density. Moreover, to compute one of these terms (e.g. the electrostatic interaction with the dielectric
continuum), Environ uses an iterative approach, so at every SCF step, there is an additional iterative procedure
performed by Environ, let's call this the polarization iteration, in contrast to the SCF iteration.

There are two Environ parameters which are crucial to help convergence of the total SCF calculation, 
``tol``, which controls the convergence of the polarization iteration, and ``environ_thr``, which determines
when the polarization iteration starts. 

When you decrease ``tol``, the extra term coming from the interaction with the dielectric will be computed more
accurately, but it will require more time to be computed. On the other hand, if the parameter is not set to
be sufficiently accurate, the SCF procedure may not convergence as quickly, thus slowing the simulation down
noticeably. Furthermore, a low convergence speed in the SCF may leads to oscillation around the SCF threshold,
leading to convergence failure. Due to the fact that in general, the polarization iteration is less expensive
than the SCF iteration, it is recommended that the user opts for a decreased ``tol`` when encountering these
issues. Hence if in doubt, start with a generic value for the tolerance, be it the default or some value taken
from the examples (depending on what you're trying to simulate), and decrease if necessary. 

The existence of ``environ_thr`` allows for better convergence in general, since it is advisable not to add
polarization correction terms in the first few SCF steps, due to the fact the accuracy of these energy terms
will still be poor (although if one is starting from a good set of input parameters, it is probably fine to
enable the polarization correction from the start of the iteration, in which case the 
``environ_restart = .true.`` will enable this). On the other hand, waiting too long to
start adding Environ corrections can lead to poor convergence due to the fact that before the polarization
iteration corrections kick in, the system is converging for vacuum, not a solution. Hence one can analyze the
output file and see when the correction terms begin, in order to judge whether ``environ_thr`` is good or not.
If the SCF accuracy explodes or converges poorly after the polarization iterations begin, it is recommended
that one increases ``environ_thr``. One rule of thumb to go by would be to make sure that Environ corrections
skip only around 3-5 SCF steps.

A further source of errors may arise due to the use of pseudopotentials for handling core electrons. This
problem is more common in transition metals of halogens, and it comes from the neglect of these core
electrons in the charge density used to build the interface with the environment. In some cases there appears
to be a hole in the charge density on the nuclei, which is then filled with the environment (in fact, a
couple of the examples in the tutorial address this fact and account for it in the input file). The solution
is to use ``solvent_mode = 'full'``. This option adds additional gaussians centered on each atom, effectively
describing the core electrons and thus not allowing the environment to permeate inside the atoms, while
not compromising the definition of the interface. 

Another suggestion is to simplify the physics in the problem. In particular, since in most applications the 
most important effects of the environment are the electrostatic ones, you may want to switch off all the other
non-electrostatic terms. This is possible by setting ``environ_type = 'input'`` and manually setting the
parameters that you want. By only setting ``env_static_permittivity = 78.3``, one only retains the dielectric
continuum in the calculation, since all other contributions are off by default.

When simulating a 2D system, convergence of the polarization potential is made more difficult by the artificial
finite electric field coming from periodic boundary conditions. To avoid this, it would be better to use a PBC
correction, in particular for slabs in Environ, a parabolic correction is implemented, and alternatively, one
can utilize the Quantum ESPRESSO isolated options. See example 3 for a demonstration on isolated systems. 
Increasing the cellsize perpendicular to the slab will decrease periodic effects but increase convergence times
due to the larger cell.

One reason of poor convergence is simply a complex shape of the interface. New releases implement models
designed to deal with some of these complex shapes. 

Harris-Foulkes Estimate
-----------------------

The Harris-Foulkes estimate does not take into account the extra terms coming from Environ, as it would cost
more than necessary for a quantity that is only used to have an estimate. As a matter of fact, also the SCF
accuracy does not take into full account of the environment, but again it would need a secondary correction
that would cost more than its utility.

Environ Total-Energy Sum
------------------------

The energy output in an Environ calculation may be confusing in that the total energy is not the sum of the
reported terms and it is not immediately clear what solvation energy means, but before going into the details
of why, one point should be addressed. ``dGsol``, and in particular ``dGel``, are not quantities that you
get from a single calculation in solution, but rather are obtained as the difference in energy from a 
calculation fully optimized (nuclei and electrons) in solution minus a calculation fully optimized (nuclei and
electrons) in vacuum. This is explained in Example 1 of the Tutorial. In general, for these values, one needs
two geometry optimization (relax) calculations for these energy values that are typically reported in our 
publications. Hence solvation energy requires two calculations. 

The different terms do not add to the same energy due to the fact that there is a spurious extra term which is
subtracted from the one-electron contribution, but is not reported in the output. The reason is that total 
energy is computed in pw making use of the secular expression (see for example W. E. Pickett, Computer Physics
Reports, 9(3), 115-198, 1989, equations 5.21-5.22 etc.), thus from the sum of occupied states eigenvalues
(so called band structure energy, "eband" in pw.x), in which some contributions are included correctly (the
kinetic and the electron-nuclei interaction), some are double counted (hartree), some need to be removed and
added back with the correct expression (xc term), and some are missing (nuclei-nuclei interaction via ewald sum).
What PW reports in the output as the one-electron contribution is in fact ebands minus the double counted
hartree and minus the wrong xc term (these spurious terms are summed up in pw.x into a term named "deband").

Similarly, when performing the solvated calculation, there are some terms which are double counted (the
electronic part of the solvation energy), some that are missing (the ionic part of the solvation energy), some
that have the wrong expression (pressure and cavitation) and some that need to be removed (for example a term
coming from the rho-dependence of the dielectric constant). All these spurious terms are collected in a term
named "deenviron", which has the same meaning as "deband" seen above. To get to the point, what is wrong in the
output is the one-electron contribution, which is reported including the spurious deenviron term. As the Environ
module was designed not to affect the standard printout of QE, this output was not modified, and consequently,
the terms do not sum to the correct total energy. The correction deenviron and deband terms were not reported,
due to their insignifance in anything alone. However, this may well change in future Environ releases. 
