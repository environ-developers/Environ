
Dielectric Medium
=================

Electrostatic solvation effects can be described in terms of a reaction field, which provides a stabilization
effect to the system. In classical electrostatics, this effect is represented in terms of a dielectric medium
with some boundary and a dielectric permittivity, which has different representations depending on the
application. Here we consider the permittivity as a function in space, thus making the definition of a 
sharp boundary surface redundant, which itself is consistent with our definition of the interface function. 

Fattebert and Gygi presented an electrostatic picture in terms of 3D fields, with the total electrostatic
free energy of the embedded system as

.. math::

   F[\rho^{\text{el}}, \{\mathbf{R}_i\}] = \int\rho^{\text{sys}}(\mathbf{r})\phi(\mathbf{r})d^3\mathbf{r}-\int\frac{1}{8\pi}\epsilon(\mathbf{r})\lvert\nabla\phi(\mathbf{r})\rvert^2d^3\mathbf{r}

where 
