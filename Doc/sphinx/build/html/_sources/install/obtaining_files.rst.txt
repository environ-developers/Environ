.. Environ documentation installation instructions.
   Created by Edan Bainglass on Mon Oct 5 2021.
   Contains installation instructions.


Obtaining source files
======================

QE
--

1. clone the repository and checkout :guilabel:`qe-6.3` (or later) ::

      git clone https://gitlab.com/QEF/q-e
      git checkout qe-6.3

or

1. download a qe-6.3 archive (or later) from the `QE releases`_ page
2. unpack the archive ::

      tar zxvf qe-X.Y.Z.tar.gz

   or, if your :command:`tar` doesn't recognize the :guilabel:`z` flag ::

      gunzip -c qe-X.Y.Z.tar.gz | tar xvf -

Environ
-------

1. clone the Environ repository into the QE root directory ::

      git clone https://github.com/environ-developers/Environ.git

or

1. download an archive from the `Environ releases`_ page
2. unpack the archive inside the QE root directory ::

      tar zxvf qe-X.Y.Z.tar.gz

   or, if your :command:`tar` doesn't recognize the :guilabel:`z` flag ::

      gunzip -c qe-X.Y.Z.tar.gz | tar xvf -

.. _QE releases: https://github.com/QEF/q-e/releases
.. _Environ releases: https://github.com/environ-developers/Environ/releases
