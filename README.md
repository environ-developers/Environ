# Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)

### This file is part of Environ version 3.0

    Environ 3.0 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    Environ 3.0 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more detail, either the file
    `License' in the root directory of the present distribution, or
    online at <http://www.gnu.org/licenses/>.

#

# Authors

| Name               | Institution                                                     |
| ------------------ | --------------------------------------------------------------- |
| Oliviero Andreussi | Department of Physics, University of North Texas                |
| Francesco Nattino  | THEOS and NCCR-MARVEL, Ecole Polytechnique Federale de Lausanne |
| Iurii Timrov       | THEOS and NCCR-MARVEL, Ecole Polytechnique Federale de Lausanne |
| Ismaila Dabo       | Department of Materials Science and Engineering, Penn State     |
| Nicola Marzari     | THEOS and NCCR-MARVEL, Ecole Polytechnique Federale de Lausanne |
| Quinn Campbell     | Sandia National Labs                                            |
| Edan Bainglass     | Department of Physics, University of North Texas                |
| Matthew Truscott   | Department of Physics, University of North Texas                |
| Gabriel Medrano    | Department of Physics, University of North Texas                |

#

# Compilation instructions for Environ module

    This distribution of the Environ module is intended for Quantum ESPRESSO.
    The steps below will cover the following:

    - configuration of QE and Environ
    - QE compilation
    - Environ compilation
    - QE patching, dependency updating, and recompilation

    Note that Environ now uses its own FFTXlib and UtilXlib libraries modeled
    respectively after the FFTXlib and UtilXlib libraries of QE 6.7

    If there are any problems with the QE compilation step, look for solutions
    on the PW-forum, or refer to the Quantum ESPRESSO documentation and/or website:

      http://www.quantum-espresso.org/

# QE + Environ 3.0 Installation

    The following instructions refer to Environ 3.0 installed on QE >= 7.1
    For previous versions, please refer to the instructions on the website:

      https://environ.readthedocs.io/en/latest/install/install.html

#

From Environ root:

1. configure Environ with

   - `./configure --with-qe=<absolute_path_to_QE>`
   - the script will attempt to find the necessary libraries/compilers
     - if it fails, the user may need to manually modify `make.inc`

2. compile Environ with

   - `make -jN compile`
   - `N` = number of cores for compilation (default = 1)
     - if issues arise, retry with `N` = 1

#

From QE root:

# using make

1. configure QE with

   - `./configure --with-environ=<absolute_path_to_Environ>`
   - the script will attempt to find the necessary libraries/compilers
     - if it fails, the user may need to manually modify `make.inc`

2. compile desired QE package with

   - `make -jN <package>`
   - `N` = number of cores for compilation (default = 1)
     - if issues arise, retry with `N = 1`
   - Environ supports the following QE packages: PW, CP, TDDFPT, XSpectra, NEB

# using cmake

1. create and navigate to a new directory (`build/` is a common choice)
2. generate build environment with

   - `cmake -DCMAKE_Fortran_COMPILER=<...> -DCMAKE_C_COMPILER=<...> -DENVIRON_ROOT=<absolute_path_to_Environ> ..`
   - more information can be found at https://gitlab.com/QEF/q-e/-/wikis/Developers/CMake-build-system

3. compile desired QE package as described above

#

If using cmake, modify the ESPRESSO_ROOT variable in the following files to reflect the build dir:

- Environ/tests/environment
- Environ/examples/qe/environment

#

# Uninstall Environ 3.0

From Environ root:

- `make clean`

- (OPTIONAL) `make veryclean` to also remove configuration files

---
