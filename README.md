# Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)

### This file is part of Environ version 2.0

    Environ 2.0 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    Environ 2.0 is distributed in the hope that it will be useful,
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

    Note that Environ now uses its own FFTs and utils libraries modeled
    respectively after the FFTXlib and UtilXlib libraries of QE 6.7

    If there are any problems with the QE compilation step, look for solutions
    on the PW-forum, or refer to the Quantum ESPRESSO documentation and/or website:

# http://www.quantum-espresso.org/

# QE + Environ 2.0 Installation

    The following instructions refer to Environ 2.0 installed on QE >= 6.3
    For previous versions of QE (namely QE-6.1.X and QE-6.2.X), please refer
    to the instructions on the website:

# http://www.quantum-environ.org/installation.html

From QE root:

1. configure QE installation with

   - `./configure [--prefix=]`
   - no prefix should be enough in most cases

From Environ root: (placed in the QE root directory)

1. configure Environ with

   - `./configure`
   - in most cases, no prefix should be enough

2. install QE + Environ 2.0 with

   - `make install-QE+Environ`

     a) select QE package (pw, tddfpt, xspectra, or all)

     - see QE documentation for package detail

     b) select number of cores (for parallel compilation)

     - if issues arise in parallel, run in serial (core=1)"

   * NOTE: installation pre-compiles QE, compiles Environ, applies
           Environ patches to QE, and re-compiles QE. This scheme is
           not only safer in general, but also required for QE < 6.7
           due to the use of iotk (see QE documentation for details)

   * compilation will write logs to Environ/install at each step
   * compilation will terminate if errors are detected at any step

3. run executables with -environ flag, e.g.

   - pw.x -environ < pw.in > pw.out

# Uninstall QE + Environ 2.0

From QE root:

1. uninstall QE + Environ 2.0 with

   - `make uninstall-QE+Environ`
   - (optional) uninstall QE

---
