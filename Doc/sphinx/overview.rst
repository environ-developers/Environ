
Overview
--------

As opposed to the more commonly used atomistic description of solvents in chemistry simulations, Environ
implements a set of continuum models, with the goal of providing a method of describing systems in
environments with the primary goal of maximizing the focus on the quantum mechanical system. The result is
a relatively succinct description of the environment unsuitable for capturing certain effects within this
regime, yet the advantages of reduced computational cost make this approach invaluable for many applications.

These models belong to a family of hierarchical approaches that separate the treatment of the quantum
mechanical system from the environment. Varying choices for representing this environment result in varying
groups of methods, each with their own ideal applications. For example Quantum Mechanics has been coupled 
with Molecular Mechanics (QM/MM) thus retaining atomistic details of the embedding system while still 
treating it with a less computationally intensive scheme. The drawback of such an approach is perhaps not 
going far enough in terms of cost reduction, since a statistical sampling is still required for the 
environment configurations. Additionally the accuracy of such an approach is strictly dependent on the quality
of the classical parameterization used to describe the components within the environment. 

Continuum models go one step further, now removing atomistic details in the environment under the assumption
that a statistical average of the environment results in a homogenous continuum. This can be justified for 
many solvent systems (with the exception of perhaps solvent atoms close to the QM system, in certain scenarios).
One popular continuum model within the chemistry community is the Polarizable Continuum Model (PCM), which 
has been progressively extended and optimized for simulation stability and speed over the years. The 
breakthrough in the field of condensed matter was the Fattebert and Gygi (FG) model, which kickstarted
substantial interest and many publications in the following years. The continuum models implemented in Environ
are in some ways derived from the ideas presented in the PCM and the FG model. These have been successfully
applied to neutral liquids and electrolyte solutions, with periodic crystal structures and small isolated
molecules. 

More recent developments have been made to improve these methods, especially towards better transferability,
and there has been work in proposed multi-scale approaches that combine multiple different embedding schemes,
with the same goals of effectively distributing computational time towards more important areas. 
