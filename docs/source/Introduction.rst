Introduction
=====

.. _installation:

Background
----------

In recent years, there has been an increasing demand for the development of software tools that can perform electron defect interaction calculations. Such tools are essential in materials science and engineering as they help in understanding the properties and behavior of materials at the atomic and molecular levels. The purpose of this research paper is to discuss a software tool that can perform electron defect interaction calculations and its potential applications. Electron defect interaction is a fundamental concept in materials science and engineering that describes the interaction between electrons and defects in a material. Defects are any irregularities or imperfections in the structure of a material, such as vacancies, interstitials, and dislocations. The interaction between electrons and defects can have a significant impact on the properties of materials, including their mechanical, electrical, and optical properties.
To develop advanced materials and devices, it is crucial to understand the processes involving electrons and atomic defects in the solid state. This requires computational tools that can predict the physical properties of materials by taking into account their atomic and electronic structure. EDI(electron-defect interaction) is a validated code that provides a unified platform for computing electron interactions, transport, and ultrafast dynamics in materials. It uses established first-principles methods such as density functional theory (DFT) as starting points for computing electron dynamics. The current distribution of EDI focuses on electron-defect (e-d) interactions and related transport properties, including electrical conductivity, scattering rate, and mobility. It also includes routines for computing spin-related infromation. The transport module enables accurate calculations of charge transport in a wide range of functional materials. The code is efficient with MPI parallelization, and scales linearly with the supercell size, thus only limited by the DFT code. 
It is suitable for use on both high-performance supercomputers and smaller computer clusters. Target users include experts in first-principles calculations and materials theory, as well as experimental researchers investigating charge transport, advanced functional materials, and semiconductor or solid-state devices. The paper provides a detailed description of the theory and numerical methods implemented in the code, its capabilities and workflows, technical aspects, example calculations, parallelization strategy, and future development plans.

The code is efficient with MPI parallelization, and scales linearly with the supercell size, thus only limited by the DFT code.


Methodology
----------


Boltzmann transport
^^^^
The mobility is calculated using Boltzmann transport theory under the momentum relaxation time approximation (MRTA). The carrier mobility is:, where the relaxation time can be obtained as:  
,			(1)
where nd is the concentration of defect (in our case, chalcogen vacancy concentration, nV), ΩBZis the Brillouin zone (BZ) volume, and  is the angle between the velocities of the initial state ik and final state jk+q for the scattering process. M is the matrix element:
, 					(2)
where ΔH is the perturbation potential induced by the defect, and φik and φik+q are the wavefunctions of the initial and final states. The integral of delta function is evaluated with triangular integral. 

Triangular integral
^^^^
The triangular integral is a commonly used mathematical tool to calculate physical properties of crystalline materials. It involves integrating over a triangular region in the Brillouin zone, which is a mathematical construct that represents the periodicity of the crystal lattice. The triangular integral is used to calculate quantities such as the density of states and the electronic band structure, which are important for understanding the behavior of electrons in solids. The integral takes into account the periodicity of the crystal lattice, and is used to calculate the dispersion relation of electrons, which determines their momentum and energy. The triangular integral is a fundamental tool in the study of electronic properties of materials, and is widely used in both theoretical and experimental research.

Electron interaction with neutral defect
^^^^
For neutral defects, we employ a supercell method. The perturbation potential is calculated from DFT as: 
, 					(5)
where Hd and Hp are the KS potentials of the supercells with and without the defect. The KS potential is a sum of the Hartree potential, exchange-correlation potential, and pseudopotential (PP). The first two potentials are local in real space, while the pseudopotential has a non-local (NL) component. With the norm-conserving pseudopotential, the non-local component has a separable Kleinmann-Bylander format:
,				(6)
where s is the index of atoms, l1 and l2 are indices of the projector function β for the respective atoms, σ1 and σ2 are the spin indices for the non-colinear wavefunction representation For the local matrix, we calculate it in real space as:
, 		(9)
where u is the periodic part of the Bloch wavefunction, which is calculated from the primitive cell and extended to the supercell, and Ωsup is the volume of the supercell. 
The non-local matrix can be calculated as:
.	           (10)
Here sd and sp are the atom indices of the supercell for the defect-containing structure and the pristine structure respectively. The overlap integral of wavefunction and projector  is calculated with the implementation: , where Ωuc is unit cell volume, and G is the reciprocal lattice vector of the unit cell. The hat in  and  represent the Fourier transform of u(r) and β(r) respectively30. Note that the phase factor eikr in φ is absorbed into , which can be obtained from Quantum Espresso with simple adaptions, so that we can calculate the overlap integral easily, which facilitates the implementation of Eq. (10). 

Electrons interaction with charged defect
^^^^
As mentioned before, following the widely-used approximation in literature 16–18,22,23, we treat the charged defect as a point charge and its perturbation potential as screened Coulomb potential. The matrix element can be computed as:
.		    (11)
Here  is the Fourier transform of the screened potential generated by a point charge located at r0 :
, 			(12)
where  is the fourier transform of density of point charge located at r0,  is the Coulomb kernel truncated in the z-direction (perpendicular to the basal plane) to avoid fictitious interactions between periodic images31:
,   		 (13)
where the subscript xy denotes the component on x-y plane, Lz is the periodic length of the unit cell along the vacuum direction, and Z is the number of elementary charges carried by the defect. 
ε in Eq. (12) is the dielectric function, which can be obtained from the density response function χ as
.				(14)
χ can be computed from first principles as follows: 
, (15)
where Nk is the total number of k points in the summation. Note that the Fermi distribution f in Eq. (15) depends on the Fermi level, thereby incorporating the effect of free carrier screening. However, in practice, it is computationally very expensive to converge χ against the k grid, especially for moderate to low nc. To reduce the computational cost, we follow the approximation by separating the total response function into the “intrinsic” part and “free carrier” part:
,					(16)
The intrinsic part is contributed by inter-band transition and can be calculated relatively easily using Eq. (15). The free carrier part is contributed by the intra-band transition. Assuming that the carriers are electron gas distributed on a 2D (zero-thickness) plane located at z=z0, the density response function will have the following form (see SI for derivation):
. 		(17)
For MX2, the band edge states are mainly contributed by the M atoms. Hence we put the electron gas plane at the same position as the M plane. 



Capabilities
----

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

The reference could be found at [arxiv](./another-page.html).



