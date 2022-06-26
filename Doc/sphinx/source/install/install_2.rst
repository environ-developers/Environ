.. Environ documentation installation instructions.
   Created by Edan Bainglass on Mon Oct 5 2021.
   Contains installation instructions.


Environ :guilabel:`v2.0` - :guilabel:`v3.0`
===========================================

.. note:: Run the following commands from the Environ root directory.

Installing
----------

1. install Environ ::

      make install

   1. select the QE package you wish to work with

   2. select the number of cores for parallel compilation (default = 1)

.. note:: If issues arise, try compiling in serial.


Uninstalling
------------

To uninstall Environ, run ::

      make uninstall

1. enter :command:`y` to proceed

2. when prompted for, enter :command:`y` to uninstall QE (optional)
