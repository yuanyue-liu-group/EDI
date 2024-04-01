Introduction
===============

Background
------------

In recent years, there has been an increasing demand for the development of software tools that can perform electron defect interaction calculations.
Such tools are essential in materials science and engineering as they help in understanding the properties and behavior of materials at the atomic and molecular levels.
We develop a software tool that can perform electron defect interaction calculations and further study its potential applications.

Electron defect interaction is a fundamental concept in materials science and engineering that describes the interaction between electrons and defects in a material.
Defects are any irregularities or imperfections in the structure of a material, including but not limited to cases such as vacancies, interstitials, dislocations, and grain boundaries.
The interaction between electrons and defects can have a significant impact on the properties of materials, including their mechanical, electrical, and optical properties.
To develop advanced materials and devices, it is crucial to understand the processes involving electrons and atomic defects in the solid state.
This requires computational tools that can predict the physical properties of materials by taking into account their atomic and electronic structure.

Softare
------------

EDI(electron-defect interaction) is a validated code that provides a unified platform for computing electron interactions, transport, and ultrafast dynamics in materials.
It uses established first-principles methods such as density functional theory (DFT) as starting points for computing electron dynamics.
The current distribution of EDI focuses on electron-defect (e-d) interactions and related transport properties, including electrical conductivity, scattering rate, and mobility.
It also includes routines for computing spin-related information.
The transport module enables accurate calculations of charge transport in a wide range of functional materials.
EDI is computationally efficient with good scalabilities.

It is suitable for use on both high-performance supercomputers and smaller computer clusters.
Target users include experts in first-principles calculations and materials theory, as well as experimental researchers investigating charge transport, advanced functional materials, and semiconductor or solid-state devices.
The code is efficient with MPI parallelization, and scales linearly with the supercell size, thus only limited by the DFT code.


Methodology
-------------

For neutral defects, we employ a supercell method. The perturbation potential is calculated from DFT.  
The scattering matrix element is calculated by integrating the pertubration potential in supercell and the wavefunctions,
which consists of local component integrated in real space and non-local component integrated in reciprocal space.

For charged defects, they are treated as a point charge and its perturbation potential as screened Coulomb potential. 
The dielectric function is obtained from first-principles calculations.
Several types of model are supported in the code, including dielectric constant, scalar dielectric function, and dielectric matrix.
The details of the methods could be found `here <https://pubs.acs.org/doi/10.1021/acsnano.4c01033>`_.

The e-d interaction matrix element :math:`M_{ij}=<\phi_i|V|\phi_j>` could be calculated with the main engine of edi, which facilitates calculation of scattering rate based on Fermi's golden rule under the momentum relaxation time approximation (MRTA).

.. math::
  \gamma_i= \frac{2\pi n_d}{ \hbar } \Sigma_j|M_{ij}|^2 (1-cos(\theta_{ij})) \delta(E_i-E_j)

The mobility is calculated using Boltzmann transport theory.

.. math::
  \mu= -\frac{e}{ 2 n_c } \Sigma_i \gamma_i^{-1} v_i^2 f'(E_i) 

In the sampling of k point to calculate mobility, 2 methods are implemented:

- Uniform grid: suitable for all systems, with gaussian smearing representing the delta function.

- Triangular integral algorithm: suitable for 2D system.



Capabilities
-------------

Currently, the following functions are supported by EDI:

- Calculate matrix element of electrons scattered by point defect

   * Atomic neutral defect: vacancy, substitution, interstitial, etc

   * Point-charge charged defect: various screening model

- Calculate transport property such as carrier mobility 

- Calculate spin relaxation time

Functions currently under development:

- Calculate atomic charged defect

- Calculate line defect scattering process

- Calculate plane defect scattering process

Performance 
-------------

EDI uses optimal algorithm to calculate different part of scattering matrix element, giving the optimal performance of accuracy and effeciency.
The scalability of EDI is very good for system size. 
The calcualtion cost scales linearly with the volum of super cell, making it capable of calculating large systems easily.
Calculation of matrix element is parallelized over k point pairs. 
EDI could easily run on large HPCs and utilize the full capacity.


Reference
----------

The details in the above sections could be found in `ACS Nano 2024, 18, 11, 8511–8516 <https://pubs.acs.org/doi/10.1021/acsnano.4c01033>`_.



