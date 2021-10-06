.. Environ tension file.
   Created by Matthew Truscott on Fri Apr 19 2019.

Tension
=======

Cococcioni et al. [1]_ proposed an environment effect related to the quantum surface of the embedded system. The
macroscopic physical picture is the one of surface tension, which controls the size of the boundary between the
system and environment. An embedding interaction based on the quantum surface was proposed by Scherlis et al. [2]_
to estimate the free energy penalty of creating a void into the continuum liquid solution, that is, the
cavitation free energy,

.. math::

   F^{\text{cav}}[s] = \gamma S[s],

where the surface is defined in terms of the interface function,

.. math::

   S[s] = \int\lvert\nabla s(\mathbf{r})\rvert d^3\mathbf{r}.

Various publications have generalized and revised this cavitation free energy function [3]_, [4]_, the 
implementation of the surface tension being a result of this work. The functional derivative of the energy 
contribution with respect to the interface function is

.. math::

   \frac{\delta F^{\text{cav}}[s]}{\delta s}(\mathbf{r}) = -\nabla\cdot\left(\frac{\nabla s(\mathbf{r})}{\lvert\nabla s(\mathbf{r}\rvert}\right),

which is used in the calculation of the contribution to the Kohn-Sham potential or the interatomic forces.

.. [1] Cococcioni M, Mauri F, Ceder G, Marzari N, Phys. Rev. Lett. 2005 4;94(14):145501
.. [2] D. A. Scherlis, J. L. Fattebert, F. Gygi, M. Cococcioni, N. Marzari, J. Chem. Phys. 2006, 124(7), 74103
.. [3] Andreussi O, Dabo I, Marzari N, J. Chem. Phys. 2012 2;136(6):064102
.. [4] G. Fisicaro, L. Genovese, O. Andreussi, S. Mandal, N. N. Nair, N. Marzari, et al., J. Chem. Theory Comput. 2017, 13(8), 3829
