
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

where the electrostatic density term can be decomposed into electronic and ionic terms,

.. math::

   \rho^{\text{sys}}(\mathbf{r}) = \left(\rho^{\text{el}}(\mathbf{r}) + \sum_i z_i\delta(\lvert\mathbf{r}-\mathbf{R}_i\rvert)\right),

:math:`z_i` are the ionic charges, :math:`\phi(\mathbf{r})` is the electrostatic potential, and
:math:`\epsilon(\mathbf{r})` is the dielectric permittivity, which reflects the system/continuum separation
and is thus a function of the continuum interface, for example,

.. math::

   \epsilon(\mathbf{r}) \equiv \epsilon(s(\mathbf{r})) = 1 + (\epsilon_0 - 1)(1 - s(\mathbf{r})),

where :math:`\epsilon_0` is the bulk dielectric permittivity of the environment. 

.. note::

   This dielectric permittivity is the constant set in the Environ input file. Hence Environ takes this
   value, and calculates an interface function that scales from 0 (environment) to 1 (system). As one can
   see from the equation, the result is a smoothly scaling dielectric permittivity function that transitions
   at the computed interface, from 1 (system) representing vacuum, to the permittivity specified by the
   environment, say 80 for water.

The functional derivative of the total electrostatic free energy functional defined above with respect to the
electrostatic potential is simply the generalized Poisson equation (GPE), 

.. math::

   \nabla\cdot\epsilon(\mathbf{r})\nabla\phi(\mathbf{r}) = -4\pi\rho^{\text{sys}}(\mathbf{r}),

which is similar to the more popularly known Poisson equation in vacuum, 

.. math::

   \nabla^2\phi^0(\mathbf{r}) = -4\pi\rho^{\text{sys}}(\mathbf{r}),

following the convention that the dielectric permittivity of vacuum is simply 1. The polarization potential
is simply the difference between the vacuum potential and the generalized potential. The generalized Poisson
equation is non-trivial to solve, and an alternative apporach is to define an auxiliary polarization density,
spread at the continuum interface, and expressed as

.. math::

   \rho^{\text{pol}}(\mathbf{r}) = \frac{1}{4\pi}\nabla\ln\epsilon(\mathbf{r})\cdot\nabla\phi(\mathbf{r}) - \frac{\epsilon(\mathbf{r}) - 1}{\epsilon(\mathbf{r})}\rho^{\text{sys}}(\mathbf{r}).

This depends on the gradient of the logarithm of the dielectric function, and thus Andreussi et al. proposed
an alternative formulation of the dielectric in terms of the continuum interface,

.. math::

   \epsilon(\mathbf{r}) = \exp(\log(\epsilon_0[1 - s(\mathbf{r}))).

When optimizing the embedded system degrees of freedom, the functional derivatives of the energy with
respect to the electronic density or atomic positions need to be computed. It is possible to break down these
derivatives into more transferable forms,

.. math::

   V^{\text{interface}}_{\text{KS}}(\mathbf{r}) = \int\frac{\delta s(\mathbf{r}^{\prime})}{\delta\rho^{\text{el}}(\mathbf{r})}\frac{\delta F[s]}{\delta s(\mathbf{r}^{\prime})}d^3(\mathbf{r})^{\prime},

and

.. math::

   f^{\text{interface}}_i = -\int\frac{\delta s(\mathbf{r})}{\delta\mathbf{R}_a}\frac{\delta F[s]}{\delta s(\mathbf{r})}d^3(\mathbf{r}).
