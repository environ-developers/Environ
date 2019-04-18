.. Environ documentation installation instructions, created
   by Matthew Truscott on Tue Mar 26 2019.
   Contains installation instructions.

Installation instructions
=========================

Environ requires a Quantum ESPRESSO installation. At the very simplest, one can install the pw.x program
by cloning quantum-espresso `here`_ and checking out the tag qe-6.4 for the most up to date version. 
Alternatively one can download the archive at the `release`_ page. Then,

1. configure QE following the standard procedure, running::

      ./configure

2. compile QE without the Environ module (i.e. standard pw installation)::
   
      make pw

If there are problems with the preliminary steps, look up for solution on the `PW-forum`_ or refer to
the Quantum ESPRESSO `documentation`_ and website. To install Environ, the following steps are
necessary. All commands should be executed in the root directory of Quantum ESPRESSO, in which the
Environ files should be extracted. The Environ files are assumed to exist (see step 1 below) in
a folder named Environ that sits in the Quantum ESPRESSO root directory

1. If not done already, clone (or extract) the Environ module to the root of the Quantum ESPRESSO installation
   ::
   
      git clone https://gitlab.com/olivieroandreussi/Environ.git

2. run the QE script addonpatch.sh with the -patch option::

      ./install/addsonpatch.sh Environ Environ/src Modules -patch

3. run the Environ installation script with the -patch option::

      ./Environ/patches/environpatch.sh -patch

4. run the QE script to regenerate modules' dependencies::

      ./install/makedeps.sh

5. re-compile, e.g.::

      make pw

To clean up the compilation of Environ, follow these steps:

1. run the QE script addsonpatch.sh with the -revert option::

      ./install/addsonpatch.sh Environ Environ/src Modules -revert

2. run the Environ installation script with the -revert option::

      ./Environ/patches/environpatch.sh -revert

3. run the QE script to regenerate modules' dependencies::

      ./install/makedeps.sh

4. remove objects, modules, and exectuables::

      make clean

Installation Issues
-------------------

Linker Problem
~~~~~~~~~~~~~~

In Environ v1.0, there is an issue with the installation procedure for codes different from pw.x. The problem
seems to depend on the compiler, but it is present in the most common Intel Fortran compiler. The solution
for this problem requires some work, which is described here, and will be fixed in future releases.

There is some circular dependencies between Environ modules and PW modules, but since the linkers only look 
for the dependencies written at the time of the Environ installation, they will not be able to function
as intended.

To fix this, one can edit the Makefile manually, to add these PW dependencies, so that instead of::

   PWOBJS = some-path/PW/src/libpw.a
   QEMODS = some-path/Modules/libquemod.a some-more-libraries.a

we have::

   PWOBJS = some-path/PW/src/libpw.a
   QEMODS = some-path/Modules/libquemod.a some-more-libraries.a some-path/PW/src/libpw.a

This may be necessary for the following files::

   PHonon/FD/Makefile
   PHonon/Gamma/Makefile
   PHonon/PH/Makefile
   PWCOND/src/Makefile
   TDDFPT/src/Makefile
   XSpectra/src/Makefile
   GWW/bse/Makefile
   GWW/head/Makefile
   GWW/pw4gww/Makefile

There is an exception, for the CPV Makefile, instead of this::

   QEMODS=../../Modules/libqemod.a

do this::

   QEMODS=../../Modules/libqemod.a libcp.a

.. _here: https://gitlab.com/QEF/q-e
.. _release: https://github.com/QEF/q-e/releases
.. _PW-forum: https://www.quantum-espresso.org/forum
.. _documentation: https://www.quantum-espresso.org/Doc/user_guide/
