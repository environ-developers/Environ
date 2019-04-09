.. Environ documentation interface models file, created by
   Matthew Truscott on Mon Apr 8 2019. Contains general
   description and comparison of interface models.

Interface Models
================

This page presents a more in-depth description of the interface models implemented in Environ, along
with all revelant input parameters. 

Continuum Models
----------------

As opposed to the more commonly used atomistic description of solvents in chemistry simulations, Environ
implements a set of continuum models, with the goal of providing a method of describing systems in
environments with the primary goal of maximizing the focus on the quantum mechanical system. The result is
a relatively succinct description of the environment unsuitable for capturing certain effects within this
regime, yet the advantages of reduced computational cost make this approach invaluable for many applications.



Self-Consistent Continuum Solvation (SCCS)
------------------------------------------

The SCCS model is described comprehensively in the 2012 publication [1]_, we present here a summary of the
theory and methodology behind the model.

The idea of this interface model 

.. [1] O. Andreussi, I. Dabo, and N. Marzari, J. Chem. Phys. 136, 064102 (2012)
