.. Environ documentation installation instructions, created
   by Matthew Truscott on Tue Mar 26 2019.
   Contains installation instructions.

Installation Instructions
=========================

Environ requires a Quantum ESPRESSO installation. At the very simplest, one can install the pw.x program
by cloning quantum-espresso `here`_ and checking out the tag qe-6.4 for the most up to date version. 
Alternatively one can download the archive at the `release`_ page. Then,

1. configure QE following the standard procedure, running::

   ./configure

2. compile QE without the Environ module (i.e. standard pw installation)::
   
   make pw

.. _here: https://gitlab.com/QEF/q-e.git 
.. _release: https://github.com/QEF/q-e/releases
