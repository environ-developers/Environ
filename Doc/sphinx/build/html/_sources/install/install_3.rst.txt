.. Environ documentation installation instructions.
   Created by Edan Bainglass on Sun Jun 26 2022.


Environ :guilabel:`v3.0` and up
===============================

Installing
----------

From Environ root ::

      make -jN compile

.. note::
      `N` is the number of cores to use in parallel for compilation. If issues arise, use `N` = 1

From QE root, or QE/build if configured with `cmake` ::

      make -jN <package>

.. note::
      Environ supports the following QE packages: PW, CP, TDDFPT, XSpectra, and NEB

Uninstalling
------------

From Environ root ::

      make clean

To also remove the configuration files ::

      make veryclean
