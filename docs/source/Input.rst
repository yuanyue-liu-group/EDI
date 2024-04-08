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
lvacalign                    vacuum aligment
vac_idx                     vacuum alignment location
lcorealign                   core alignment
core_v_d                    core alignment value for defect
core_v_p                    core alignment value for pristine
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
calcmcharge                 calculate charged defect
mcharge_dolfa               use LFA approximation in charged calculation
qeh_eps_filename            dielectric function file from QEH
doqeh                       use QEH dielectric function 
dogw                        use BGW dielectric function
k0screen_read               Lindhard model carrier screening
gw_epsmat_filename          BGW dielectric function file for grid q
gw_eps0mat_filename          BGW dielectric function file for small q
====================      ======================================


An example input file is shown below:

.. code-block:: console

    &calcmcontrol
    prefix='mos2',
    outdir='dout/'
    lvacalign=.true.
    vac_idx=0
    lcorealign=.false.
    core_v_d=0.0
    core_v_p=0.0 
    wt_filename='wt.dat'
    klist_filename='scfklist.dat'
    ev_filename='v.dat'
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
    calcmcharge=.true.
    mcharge_dolfa=.true.
    qeh_eps_filename='./eps.dat'
    eps_type='gw'
    dogw=.true.
    chidat='./chi.dat'
    !eps_type='qeh'
    !eps_type='tf'
    !doqeh=.true.
    !k0screen_read=0.27
    gw_epsmat_filename='./epsmat.h5'
    gw_eps0mat_filename='./eps0mat.h5'
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
* core alignment

Vacuum alignment could be used for materials confined along at least one direction, where a 2D plane in vacuum with location set in the input file will be used to calculate an averaged energy as a reference point.
Currently, only plane perpendicular to z direction is supported.
To use vacuum alignment, set ``lvacalign`` to ``.true.``.
``vac_idx`` also needs to be set.
This parameter sets the location of the vacuum plane, in the form of the FFT grid number index from the DFT calculation.


Core alignment could be used for all materials, the value should be the core level energies of proper element. 
To use core alignment, set ``lcorealign`` to ``.true.``.
``core_v_d`` and ``core_v_p`` needs to be set for the core level corrections in this option.
The represent the core level energy of defect and pristine structures respectively.


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




Charged defect perturbtation potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If defect is charged, the perturbation potential is represented with a different model from neutral ones.
Currently, supported model is Coulomb potential of a point charge, screened by the material. 
Various screening model is supported by EDI.

To perform this calculation, set the parameter ``calcmcharge`` to ``.true.``.

Local Fielad Approximation (LFA) is supported for the charged defect systems.
To turn on, set the parameter ``mcharge_dolfa`` to ``.true.``.

Currently, the supported screening models include:

* Thomas-Fermi model with dielectric constant

    Set ``k0screen_read`` to use dielectric constant

* Quantum Electrostatic Hetereostructure model (scalar dielectric function)

    Set ``doqeh`` to use QEH dielectric function.

    Set ``qeh_eps_filename`` for the dielectric function obtainedfrom QEH model

* Lindhard model (matrix dielectric function)

    Set ``dogw`` to use BGW dielectric function

    Set ``gw_epsmat_filename`` for the full dielectric matrix for q grid obtained from BGW

    Set ``gw_eps0mat_filename`` for the full dielectric matrix for small q obtained from BGW
