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

      tar zxvf qe-X.Y.Z.tar.gz

or, if your :command:`tar` doesn't recognize the :guilabel:`z` flag ::

      gunzip -c qe-X.Y.Z.tar.gz | tar xvf -

.. note::
      Environ :guilabel:`v3.0` supports QE :guilabel:`qe-7.1` and up. Prior Environ versions support QE :guilabel:`qe-6.3` - :guilabel:`qe-7.0`

Environ
-------

Clone the Environ repository ::

      git clone https://github.com/environ-developers/Environ.git

or download an archive from the `Environ releases`_ page
and follow the above directions to unpack the archive

.. note::
      Environ versions prior to :guilabel:`v3.0` require Environ root to be placed inside QE root

.. _QE releases: https://github.com/QEF/q-e/releases
.. _Environ releases: https://github.com/environ-developers/Environ/releases
