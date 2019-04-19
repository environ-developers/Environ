.. Environ pressure file, created by
   Matthew Truscott on Fri 19 2019

Pressure
========

There is a pressure field caused by the embedding environment. We represent this pressure field in a continuum
fashion, following the concept of electronic enthalpy introduced by Cococcioni et al [1]_. The approach exploits
the notion of quantum volume, which can be straightforwardly defined in terms of a continuum interface function
as

.. math::
   
   V \equiv V[s] = \int s(\mathbf{r})d^3\mathbf{r}.

The contribution to the enthalpy of the system is expressed as

.. math::

   G^{\text{PV}}[s] = P^{\text{ext}}V[s],

with :math:`P^{\text{ext}}` the external pressure of the environment. This approach is based on an electronic
interface, but equivalent arguments may be adopted for ionic or mixed interfaces. The key ingredient for the 
derivation of the enthalpy contribution to the Kohn-Sham potential or the inter-atomic forces is the functional
derivative of the additional energetic term with respect to the interface function, which for the equation
defining the continuum interface function above, is simply

.. math::

   \frac{\delta F^{\text{PV}}[s]}{\delta s}(\mathbf{r}) = 1.

.. [1] Cococcioni M, Mauri F, Ceder G, Marzari N, Phys. Rev. Lett. 2005 4;94(14):145501
