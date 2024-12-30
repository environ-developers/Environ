.. Environ documentation installation instructions.
   Created by Edan Bainglass on Mon Oct 5 2021.
   Updated by Edan Bainglass on Sun Jun 26 2022.
   Contains installation instructions.


Obtaining source files
======================

QE
--

Clone the repository and checkout the relevant QE version (see note below) ::

      git clone https://gitlab.com/QEF/q-e
      git checkout <tag>

or download the relevant version's archive from the `QE releases`_ page
and unpack the archive using ::

      tar zxvf q-e-qe-X.Y.Z.tar.gz

or, if your :command:`tar` doesn't recognize the :guilabel:`z` flag ::

      gunzip -c q-e-qe-X.Y.Z.tar.gz | tar xvf -

+-------------------+-----------------------------------------+
| Environ Version   | Supported QE Version                    |
+===================+=========================================+
| :guilabel:`v3.1+` | :guilabel:`v7.1+`                       |
+-------------------+-----------------------------------------+
| :guilabel:`v3.0`  | :guilabel:`v7.1` - :guilabel:`v7.2`     |
+-------------------+-----------------------------------------+
| :guilabel:`<v3.0` | :guilabel:`v6.3` - :guilabel:`v7.0`     |
+-------------------+-----------------------------------------+

.. note::
    Quantum Espresso versions :guilabel:`>v7.3` require Environ version :guilabel:`>v3.1`

Environ
-------

Clone the Environ repository ::

      git clone https://github.com/environ-developers/Environ.git

or download an archive from the `Environ releases`_ page
and follow the above directions to unpack the archive

.. note::
      Environ versions prior to :guilabel:`v3.0` require Environ root to be placed inside QE root
      Environ versions :guilabel:`v3.0+` MUST NOT be placed inside the QE root directory

.. _QE releases: https://gitlab.com/QEF/q-e/-/releases
.. _Environ releases: https://github.com/environ-developers/Environ/releases
