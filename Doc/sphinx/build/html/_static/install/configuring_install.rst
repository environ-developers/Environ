.. Environ documentation installation instructions.
   Created by Edan Bainglass on Mon Oct 5 2021.
   Updated by Edan Bainglass on Sun Jun 26 2022.
   Contains installation instructions.


Configuration
=============

QE
--

make
####

QE can use its :command:`configure` script to detect compilers and libraries
necessary for compilation. To configure the environment, run the following
command from the QE root directory ::

      ./configure

For additional flags for fine-tuning the configuration process,
check out the `QE docs`_.

cmake
#####

.. note::
      When using `cmake`, Environ must be pre-compiled. See :doc:`install_3` compilation instructions

Alternatively, one can use QE's `cmake` installation process. If chosen,
create a new build directory in QE root, and inside it, run ::

      cmake -DCMAKE_Fortran_COMPILER=<...> -DCMAKE_C_COMPILER=<...> -DENVIRON_ROOT=<absolute_path_to_Environ> ..

For additional information on the cmake installation process,
check out the `QE cmake guide`_

Environ :guilabel:`v2.0` and up
-------------------------------

The :command:`configure` script was adopted in the recent Environ
:guilabel:`v2.0` release. To configure Environ, run the same command from
the Environ root directory ::

      ./configure

making sure to use the same compiler flags and libraries.

.. note::
      If using cmake to compile QE, modify the ESPRESSO_ROOT variable in the following files to reflect the build directory: Environ/tests/environment, Environ/examples/qe/environment


.. _QE docs: https://www.quantum-espresso.org/Doc/user_guide/
.. _QE cmake guide: https://gitlab.com/QEF/q-e/-/wikis/Developers/CMake-build-system
