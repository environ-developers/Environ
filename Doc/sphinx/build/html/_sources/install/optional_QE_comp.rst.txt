.. Environ documentation installation instructions.
   Created by Edan Bainglass on Mon Oct 5 2021.
   Contains installation instructions.


Pre-compilation
===============

.. note:: Optional for Environ :guilabel:`v2.0` and up.

The installation process of Environ :guilabel:`v2.0` (or later) will automatically pre-compile QE and will abort if any issues arise. However, if using an older version of Environ, QE pre-compilation is required at this stage. You may do so from the QE root directory with ::

      make <package>

where *package* is any of the packages supported by Environ: pw, cp, tddfpt, or xspectra.
If issues arise, solutions may be found in the `QE docs`_ or on the `QE forum`_.

.. note:: If QE is manually pre-compiled, Environ :guilabel:`v2.0` (or later) re-compilation of QE will proceed quickly.

.. _QE docs: https://www.quantum-espresso.org/Doc/user_guide/
.. _QE forum: https://www.quantum-espresso.org/forum
