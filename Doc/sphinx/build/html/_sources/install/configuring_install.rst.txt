.. Environ documentation installation instructions.
   Created by Edan Bainglass on Mon Oct 5 2021.
   Contains installation instructions.


Configuration
=============

QE
--

Quantum-ESPRESSO uses a :command:`configure` script to detect compilers and libraries necessary for compilation. To configure the environment, run the following command from the QE root directory ::

      ./configure

For additional flags for fine-tuning the configuration process, check out the `QE docs`_.

Environ :guilabel:`v2.0` and up
-------------------------------

The :command:`configure` script was adopted in the recent Environ :guilabel:`v2.0` release. To configure Environ, run the same command from the Environ root directory ::

      ./configure

making sure to use the same compiler flags and libraries.

.. _QE docs: https://www.quantum-espresso.org/Doc/user_guide/
