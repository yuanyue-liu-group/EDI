Introduction
===============

Background
------------

In recent years, there has been an increasing demand for the development of software tools that can perform electron defect interaction calculations.
Such tools are essential in materials science and engineering as they help in understanding the properties and behavior of materials at the atomic and molecular levels.
The interaction between electrons and defects can have a significant impact on the properties of materials, including their mechanical, electrical, and optical properties.
To develop advanced materials and devices, it is crucial to understand the processes involving electrons and atomic defects in the solid state.
We develop a software tool that can perform electron defect interaction calculations and further study its potential applications.


Softare
------------

EDI(electron-defect interaction) uses established first-principles methods such as density functional theory (DFT) as starting points for computing electron dynamics.
The current distribution of EDI focuses on electron-defect (e-d) interactions and related transport properties, including electrical conductivity, scattering rate, and mobility.
It also includes routines for computing spin-related information.
The transport module enables accurate calculations of charge transport in a wide range of functional materials.
EDI is computationally efficient with good scalabilities.
It is suitable for use on both high-performance supercomputers and smaller computer clusters.


Methodology
-------------

The e-d interaction matrix element :math:`M_{ij}=<\phi_i|V|\phi_j>` could be calculated with the main engine of edi.
For neutral defects, we employ a supercell method. The perturbation potential is calculated from DFT.  
The scattering matrix element is calculated by integrating the pertubration potential in supercell and the wavefunctions,
which consists of local component integrated in real space and non-local component integrated in reciprocal space.

.. math::
  M_{ij}= M_{ij}^L+ M_{ij}^{NL}

The local part is integrated in real space:

.. math::
  M_{ij}^L= \int d^3r \phi_i(r)^* \phi_j(r) \Delta H(r) 

And the non-local part comes from the pseudo-potential, in separable Kleinman-Bylander form, the non-local pseudo-potential is:

.. math::
  V^{NL}= \Sigma_{mn} | \beta_m> D_{mn} < \beta_n |

The corresponding matrix element is:

.. math::
  M_{ij}^{NL}= \Sigma_{mn} <\phi_i | \beta_m^P> D_{mn} < \beta_n^P | \phi_j>- \Sigma_{mn} <\phi_i | \beta_m^D> D_{mn} < \beta_n^D | \phi_j>

Here the superscript `P` and `D` denotes pristine and defect structures respectively.
Here the overlap integral between wavefunction and projector :math:`<\phi | \beta>` is calculated in reciprocal space.
The details of the methods could be found `here <https://pubs.acs.org/doi/10.1021/acsnano.4c01033>`_.



The scattering rate is calculated from matrix element based on Fermi's golden rule.
The equation of the momentum relaxation time approximation (MRTA) is as follows.

.. math::
  \gamma_i= \frac{2\pi n_d}{ \hbar } \Sigma_j|M_{ij}|^2 (1-cos(\theta_{ij})) \delta(E_i-E_j)

The self-energy relaxation time approximation (SERTA) is similar, without angle term:

.. math::
  \gamma_i= \frac{2\pi n_d}{ \hbar } \Sigma_j|M_{ij}|^2  \delta(E_i-E_j)

The mobility is calculated using Boltzmann transport theory.

.. math::
  \mu= -\frac{e}{ 2 n_c } \Sigma_i \gamma_i^{-1} v_i^2 f'(E_i) 

In the sampling of k point to calculate mobility, 2 methods are implemented:

- Uniform grid: suitable for all systems, with gaussian smearing representing the delta function.

- Triangular integral algorithm: suitable for 2D system.



Capabilities
-------------

Currently, the following properties are calculated by EDI:

- EDI matrix
- scattering rate/relaxation time
- carrier mobility/conductivity 


Calculation or following situations currently under development:


- charged defects
- surfaces
- grain boundaries

Reference
----------

The details in the above sections could be found in `ACS Nano 2024, 18, 11, 8511–8516 <https://pubs.acs.org/doi/10.1021/acsnano.4c01033>`_.



