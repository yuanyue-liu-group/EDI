---
layout: default
---


  {:refdef: style="text-align: center;"}
 ![fig1](./figs/TOC.png){:height="394px" width="429px"}
  {: refdef}

It is crucial to understand the processes involving electrons and atomic defects in the solid state, in order to develop advanced materials and devices. 
This requires computational tools that can predict the physical properties of materials by taking into account their atomic and electronic structure. 
**EDI** (**e**lectron-**d**efect **i**nteraction) is a validated code that provides a unified platform for computing electron interactions, transport, and ultrafast dynamics in materials. 
It uses established first-principles methods such as density functional theory (DFT) as starting points for computing electron dynamics. 
The current distribution of EDI focuses on electron-defect (e-d) interactions and related transport properties, including electrical conductivity, scattering rate, and mobility. 
It also includes routines for computing spin-related infromation. 
The transport module enables accurate calculations of charge transport in a wide range of functional materials. 
The code is efficient with MPI parallelization, and scales linearly with the supercell size, thus only limited by the DFT code.

# Capabilities

EDI has the following functions:

- Calculate electron-defect scattering matrix element 

- K-point sampling for 2D system based on triangular integral algorithm

- Calculate transport property such as carrier mobility 




  {:refdef: style="text-align: center;"}
   ![fig1](./figs/fig1.png){:height="394px" width="429px"}
  {: refdef}

# Performance 

EDI uses optimized algorithm to calculate different part of scattering matrix element, giving the optimal performance of accuracy and effeciency.
The scalability of EDI is very good for system size. 
The calcualtion cost scales linearly with the volum of super cell, making it capable of calculating large systems easily.
Calculation of matrix element is parallelizede over k point pairs. 
EDI could easily run on large HPCs and utilize the full capacity.


# Reference

 Zhongcan Xiao, Rongjing Guo, Chenmu Zhang, and Yuanyue Liu. Point Defect Limited Carrier Mobility in 2D Transition Metal Dichalcogenides.  ACS Nano.  [DOI: 10.1021/acsnano.4c01033](https://pubs.acs.org/doi/10.1021/acsnano.4c01033)

