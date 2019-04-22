# DOXYGEN README
#### Notes For Users
Environ uses doxygen to auto-generate documentation. In order to locally generate documentation, one should do the following: 

1. Download doxygen from the provided [link](http://www.stack.nl/~dimitri/doxygen/download.html). Note that doxygen has a number of prerequisites. 
2. run 'doxygen doxy-environ' from the Environ/Docs (this) folder.
3. Two folders should be generated, html and latex.
4. To generate the pdf manual, run make in the latex directory.

#### Notes For Developers
Doxygen uses a particular commenting scheme for fortran by default. Checkout the documentation which should be provided by the doxygen installation, or alternatively can be accessed from this [link](http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html). 