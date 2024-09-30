Input
=====


The fortran namelist format is used for the input of EDI.
The default input file name is `calcmdefect.dat`
Required data fils should be specified in the input file.
Different calculation options are set in input file as well.
Below is the detailed discussion of the input variables and their meanings.


Input file
------------

A brief table showing the meaning of all variables are listed below in the table:

====================      ======================================
variable                     meaning                            
====================      ======================================
prefix                       qe prefix
outdir                       qe outdir
lvacalign                    vacuum Energy aligment
vac_idx                     vacuum  location for  Energy aligment
lconstalign                   constant Energy alignment
Eshift_d                    energy shift for defect
Eshift_p                    energy shift for pristine
wt_filename                  weight file for mobility calculation
degauss                     gaussian smearing in delta function
noncolin                     non-colinear calculation
lspinorb                     spin-orbit calculation
calcmlocal                   calculate local part M
calcmnonlocal                   calculate non-local part M
V_d_filename                 defect system local potential 
Bxc_1_d_filename             defect system magnetic field along x
Bxc_2_d_filename             defect system magnetic field along y
Bxc_3_d_filename             defect system magnetic field along z
V_p_filename                 pristine system local potential
Bxc_1_p_filename             pristine system magnetic field along x
Bxc_2_p_filename             pristine system magnetic field along y
Bxc_3_p_filename             pristine system magnetic field along z
====================      ======================================


An example input file is shown below:

.. code-block:: console

    &calcmcontrol
    prefix='mos2',
    outdir='dout/'
    lvacalign=.true.
    vac_idx=0
    lconstalign=.false.
    Eshift_d=0.0
    Eshift_p=0.0 
    wt_filename='wt.dat'
    noncolin =.true.
    lspinorb =.true.
    calcmlocal = .true.
    calcmnonlocal = .true.
    V_d_filename='./V_d.dat'
    Bxc_1_d_filename='./Bxc_1_d.dat'
    Bxc_2_d_filename='./Bxc_2_d.dat'
    Bxc_3_d_filename='./Bxc_3_d.dat'
    V_p_filename='./V_p.dat'
    Bxc_1_p_filename='./Bxc_1_p.dat'
    Bxc_2_p_filename='./Bxc_2_p.dat'
    Bxc_3_p_filename='./Bxc_3_p.dat'
    /

Input parameters
----------------

A detailed description of the input parameters is as follows:

QE parameters 
^^^^^^^^^^^^^^^^^^

These parameters should be the same as used in QE:
 ``prefix``, ``outdir``, ``noncolin``,  ``lspinorb``  

Energy alignment
^^^^^^^^^^^^^^^^
The energies calculated from different systems may not be able to directly compare. 
In order to obtain the correct perturbation potential, we need to choose proper energy alignment methods.
EDI provides 2 types of energy alignment algorithms:

* vacuum alignment
* constant alignment

Vacuum alignment could be used for materials confined along at least one direction, where a 2D plane in vacuum with location set in the input file will be used to calculate an averaged energy as a reference point.
Currently, only plane perpendicular to z direction is supported.
To use vacuum alignment, set ``lvacalign`` to ``.true.``.
``vac_idx`` also needs to be set.
This parameter sets the location of the vacuum plane, in the form of the FFT grid number index from the DFT calculation.


Constant alignment could be used for all materials, applying a constant shift to the two structures, the value should be the core level energies of proper element, or any other chioces see fit by the users. 
To use constatn alignment, set ``lconstalign`` to ``.true.``.
``Eshift_d`` and ``Eshift_p`` needs to be set for the constatn level corrections in this option.
They represent the reference energy of defect and pristine structures respectively.


K point sampling
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The initial and final wavefunctions for the scattering process are needed for the calcualtion of matrix elment.
The band number and k points are the index for the wavefunctions.  
Two methods are provided for the sampling of k points:
The index of the wavefunction pairs are given in the weight file, which is set by parameter ``wt_filename``.
The weight file can be obtained with the provided scripts.

* If uniform grid is used: 
  
   A complete list of :math:`C_N^2` kpoint pairs with the gaussian smearing is needed in the weight file.

* If triangular integral method for 2D system is used:

   The wavefunctions pairs are determined using triangular algorithm from the energy conservation term in the Fermi's golden rule.


Neutral defect perturbtation potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The neutral defect perturbation potential is separated into local and non-local parts. 
To calculate matrix element from it, set ``calcmlocal`` and ``calcmnonlocal`` to ``.true.``. 
Additionally, the following parameters should be set to determine the files for the potentials.

.. code::

  V_d_filename          
  Bxc_1_d_filename      
  Bxc_2_d_filename      
  Bxc_3_d_filename      
  V_p_filename          
  Bxc_1_p_filename      
  Bxc_2_p_filename      
  Bxc_3_p_filename      

.. note::
  The Bxc files are needed for spin-orbital coupling (SOC) calculations, and could be obtained with QE postprocessing tool.





Transport module input parameters
-------------------------------------------


The transport calculation is performed by python script ``mu.py``. 
In the script file, one can modify the parameters for respective materials.
Here is a list of the important parameters. 

====================      ======================================
variable                     meaning                            
====================      ======================================
   Ngrid                    k point mesh grid size
   withangle                Option of angle term in RTA model
   Efermi                   Fermi level
   Eb                       Band edge level
   Cd                       Relative defect concentration
   eh                       Type of carrrier
====================      ======================================

More details of the meaning are given in the following.
``Ngrid`` needs to be the same as the k point mesh grid size used in the wannier interpolation.
If the option ``withangle`` is true, then the MRTA model is used, otherwise the SERTA model is used;
``Efermi`` is the absolute location relative to the band edge, which is set by ``Eb``.
The relative defect concentration is set by ``Cd``.
The option ``eh`` determines the type of carrrier, ``e`` for electron, ``h`` for hole.
