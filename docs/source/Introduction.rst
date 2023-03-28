Introduction
=====

.. _installation:

Background
------------

In recent years, there has been an increasing demand for the development of software tools that can perform electron defect interaction calculations.
Such tools are essential in materials science and engineering as they help in understanding the properties and behavior of materials at the atomic and molecular levels.
We develop a software tool that can perform electron defect interaction calculations and its potential applications.

Electron defect interaction is a fundamental concept in materials science and engineering that describes the interaction between electrons and defects in a material.
Defects are any irregularities or imperfections in the structure of a material, such as vacancies, interstitials, and dislocations.
The interaction between electrons and defects can have a significant impact on the properties of materials, including their mechanical, electrical, and optical properties.
To develop advanced materials and devices, it is crucial to understand the processes involving electrons and atomic defects in the solid state.
This requires computational tools that can predict the physical properties of materials by taking into account their atomic and electronic structure.

EDI(electron-defect interaction) is a validated code that provides a unified platform for computing electron interactions, transport, and ultrafast dynamics in materials.
It uses established first-principles methods such as density functional theory (DFT) as starting points for computing electron dynamics.
The current distribution of EDI focuses on electron-defect (e-d) interactions and related transport properties, including electrical conductivity, scattering rate, and mobility.
It also includes routines for computing spin-related infromation.
The transport module enables accurate calculations of charge transport in a wide range of functional materials.
The code is efficient with MPI parallelization, and scales linearly with the supercell size, thus only limited by the DFT code.

It is suitable for use on both high-performance supercomputers and smaller computer clusters.
Target users include experts in first-principles calculations and materials theory, as well as experimental researchers investigating charge transport, advanced functional materials, and semiconductor or solid-state devices.
The code is efficient with MPI parallelization, and scales linearly with the supercell size, thus only limited by the DFT code.


Methodology
----------


The mobility is calculated using Boltzmann transport theory under the momentum relaxation time approximation (MRTA).
In the sampling of k point to calculate mobility, triangular integral algorithm is used.
The triangular integral is a commonly used mathematical tool to calculate physical properties of crystalline materials.
It involves integrating over a triangular region in the Brillouin zone, which is a mathematical construct that represents the periodicity of the crystal lattice.
It is used to calculate quantities such as the density of states and the electronic band structure, which are important for understanding the behavior of electrons in solids.
The integral takes into account the periodicity of the crystal lattice, and is used to calculate the dispersion relation of electrons, which determines their momentum and energy.
The triangular integral is a fundamental tool in the study of electronic properties of materials, and is widely used in both theoretical and experimental research.
For neutral defects, we employ a supercell method. The perturbation potential is calculated from DFT.  
The charged defect are treated as a point charge and its perturbation potential as screened Coulomb potential. 


Capabilities
----
EDI is capable of the following functions:

- Calculate electron-defect scattering matrix element 

- Find reduced k-point sampling using triangular integral algorithm

- Calculate carrier mobility limited by both neutral and charged defect as a function of carrier concentration

Performance 
----

EDI uses optimal algorithm to calculate different part of scattering matrix element, giving the optimal performance of accuracy and effeciency.
The scalability of EDI is very good for system size. 
The calcualtion cost scales linearly with the volum of super cell, making it capable of calculating large systems easily.
Calculation of matrix element is parallelizede over k point pairs. 
EDI could easily run on large HPCs and utilize the full capacity.


Reference
----

The details in the above sections could be found in reference to be released soon.



