.. Environ documentation installation instructions.
   Created by Matthew Truscott on Tue Mar 26 2019.
   Contains installation instructions.
   Updated by Edan Bainglass on Mon Oct 5 2021.


Environ :guilabel:`v1.0` up to :guilabel:`v2.0`
===============================================

.. note:: Run the following commands from the QE root directory.

Installing
----------

1. run the script :command:`addonpatch.sh` with the -patch option ::

      ./install/addsonpatch.sh Environ Environ/src Modules -patch

2. patch the necessary QE files with ::

      ./Environ/patches/environpatch.sh -patch

3. update QE dependencies ::

      ./install/makedeps.sh

4. re-compile the chosen QE package (pw, cp, tddfpt, or xspectra) ::

      make <package>


Uninstalling
------------

1. run the QE script addsonpatch.sh with the -revert option ::

      ./install/addsonpatch.sh Environ Environ/src Modules -revert

2. run the Environ installation script with the -revert option ::

      ./Environ/patches/environpatch.sh -revert

3. run the QE script to regenerate modules' dependencies ::

      ./install/makedeps.sh

4. remove objects, modules, and exectuables ::

      make clean

Known Issues
------------

In Environ :guilabel:`v1.0`, there is an issue with the installation procedure for codes different from pw.x. The problem seems to depend on the compiler, but it is present in the most common Intel Fortran compiler. The solution for this problem requires some work, which is described here, and will be fixed in future releases.

There is some circular dependencies between Environ modules and PW modules, but since the linkers only look for the dependencies written at the time of the Environ installation, they will not be able to function as intended.

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
