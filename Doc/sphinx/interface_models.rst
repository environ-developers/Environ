.. Environ documentation interface models file, created by
   Matthew Truscott on Mon Apr 8 2019. Contains general
   description and comparison of interface models.

Interface Models
================

This page presents a more in-depth description of the interface models implemented in Environ, along
with all revelant input parameters. Much of this text is structured to mirror a tutorial review on
Continuum Embeddings, [1]_. 

In general the interface is a 3D continuous, differentiable scalar field with a range from 0 to 1. 
A value of 1 (or within some tolerance) refers to the system, and a value of 0 (or within some tolerance)
refers to the environment. The boundary is a function of the ensemble of the positions of the environment
molecule. The strategy taken by continuum embedding models is to define this boundary from the system
information only. This avoids the need for statistical sampling and/or other computationally expensive
approaches. 

The tutorial review summarizes the reasons for choosing a smooth boundary over a sharp one. Environ
implements smooth boundaries since the QM simulation of the embedded system works in 3D (and reciprocal)
space. The choice to retain 3D when working with the interface bypasses the possible issues encountered
in discretization, where singularities and discontinuities are possible. We argue that from a computational
standpoint, it is more reasonable to work with a set of vector and scalar fields that share the same
numerical domain. 


Self-Consistent Continuum Solvation (SCCS)
------------------------------------------

The SCCS model is described comprehensively in the 2012 publication [2]_, we present here a summary of the
theory and methodology behind the model.

This model chooses to base the definition of the interface on the electronic density, which is a scalar
field that varies from a higher magnitude close to the ions, to a lower magnitude as we move into the
environment. This approach is based on the original model proposed by Fattebert and Gygi [3]_. 

.. math::

   s(\mathbf{r}) = \frac{1}{2}\left(1 - \frac{1 - (\rho^{\text{el}}(\mathbf{r})/\rho_0)^{2\beta}}{1 + (\rho^{\text{el}}(\mathbf{r})/\rho_0)^{2\beta}}\right)

which has two parameters, :math:`\rho_0` and :math:`\beta`. Instead it uses a piece-wise definition,

.. math::

   s(\mathbf{r}) = t(\ln(\rho^{\text{el}}(\mathbf{r})))

where :math:`t(x)` is a smooth function, and the two parameters are :math:`\rho_{\text{min}}` and
:math:`\rho_{\text{max}}` that bound the above function, the function returns 1 with an input above 
the max density bound and 0 with an input below the min density bound.

Soft-Sphere Continuum Solvation (SSCS)
--------------------------------------

The SSCS model is described comprehensively in the 2017 publication [4]_, again, we present here a
summary of the theory and methodology behind the model. 

This model choose to base the definition of the interface on the ionic positions, and takes a similar
approach to the PCM model [5]_. In the same vein as the SCCS model, rather than sticking with a 2D
definition of the boundary, the model defines its smooth 3D interface function on interlocking smooth
spheres centered on the ionic positions, and scaled depending on the atom, with an additional global
scaling for parameterization. 

Non-local interfaces
--------------------

The SCCS model can be thought of as a local interface, the value of the interface function is solely
dependent on the electronic density at that position. Likewise, the interface function of the soft-sphere
is agnostic of the 'bigger picture'. There are a number of limitations of this assumption. Firstly, 
it is easy to imagine that more complex system with artifacts such as cavities devoid of molecule 
where a solvent should
not access, will be misrepresented by a local model, which will by design fill any cavities without
regard for the physical implications. Secondly, in the SCCS, the net system charge will influence the size
of the QM cavity. Physically this is desirable since for charged systems, a solvent such as water would
interact more closely with the system, but only for positive charges, where the SCCS will shrink
the QM cavity. For negative charges, the model scales the system in the wrong direction and therefore
produces poor results without a reparameterization. 

There have been a number of recent techniques proposed in the literature that have been implemented
(or still under development) in Environ. First is the solvent-aware model, which can be toggled in
conjunction with either implemented interface model.


.. [1] O. Andreussi and G. Fisicaro, Wiley DOI: 10.1002/qua.25725 (2018)
.. [2] O. Andreussi, I. Dabo, and N. Marzari, J. Chem. Phys. 136, 064102 (2012)
.. [3] 
