Tutorial
========

In the example, we provide the instructions and files needed to calculate carrier mobility limited by S-vacancy of MoS\ :sub:`2`\ .
After downloading the package, the files could be found in ``example`` folder.

Setup
------------

This example consists of the input files for QE to generate the wavefunction and the perturbation potential modelling a S-vacancy of MoS\ :sub:`2`\  in a 6-by-6 supercell.
The obtained files could then be used to calculate the e-d interaction matrix elements with EDI.
Further instructions are then provided to obtain the carrier mobility limited defect. 

When running the example, first follow the instructions in the Installation section and compile ``pw.x``, ``pp.x``, and ``edi.x``.
Assume they are located in ``$BIN`` folder,
then copy the files in the example to a working directory to perform the calculation.


 
DFT calculation
----------------

1. Scf calculation to obtain the perturbation potential:

   Use the input files ``pristine.scf.in`` and ``defect.scf.in`` to perform scf calculations for the pristine and defect structures.

.. code-block:: console

    $ mpirun -np $N $BIN/pw.x <prinstine.scf.in  >prinstine.scf.out

    $ mpirun -np $N $BIN/pw.x <defect.scf.in     >defect.scf.out


2. Postprocessing to generate the files of the perturbation potential:

   Use the input files ``pristine.pp.in`` and ``defect.pp.in`` to generate potential data files for the pristine and defect structures.

.. code-block:: console

    $ mpirun -np $N $BIN/pp.x <prinstine.pp.in  >prinstine.pp.out

    $ mpirun -np $N $BIN/pp.x <defect.pp.in     >defect.pp.out

EDI preprocessing
-----------------

3. Calculate the band structure for full k grid

   Use the input files ``unitcell.scf.in``, ``unitcell.w90.in`` to generate a coarse grid band structure.

.. code-block:: console

    $ mpirun -np $N $BIN/pw.x <unitcell.scf.in  >unitcell.scf.out

    $ mpirun -np $N $BIN/pw.x <unitcell.w90.in  >unitcell.w90.out

..

   Use the input files  ``w90_mos2.win``, ``pw2w90.in``  to generate a fine grid band structure.  This example is using ``wannier90.x``, assuming under the same directory ``$BIN``.

.. code-block:: console

    $ mpirun -np $N $BIN/wannier90.x -pp w90_mos2

    $ mpirun -np $N $BIN/pw2wannier90.x <pw2w90.in >pw2w90.out 

    $ mpirun -np $N $BIN/wannier90.x     w90_mos2

    $ mpirun -np $N $BIN/postw90.x     w90_mos2


..

   The generated velocity file ``w90_mos2_geninterp.dat`` is needed for next step.

4. Generate the needed k points for e-d interaction matrix elements with triangular integral method

   Use the script ``wmat_mos2.py`` to obtain a list of k points in file ``kpt.dat``, and a list of k point pairs to calculate the e-d interaction matrix element together with the weight ``wt.dat``.


.. code-block:: console

    $ python wmat_mos2.py

..

    Prepare the input files for nscf calculation with the data from ``kpt.dat`` file.
    This resulted file is provided as ``unitcell.nscf.in``.


EDI calculation
----------------

5. Generate the wavefunctions for the needed k points

   Use the input files ``unitcell.scf.in`` and ``unitcell.nscf.in`` to generate wavefunctions with proper k points.

.. code-block:: console

    $ mpirun -np $N $BIN/pw.x         <unitcell.scf.in  >unitcell.scf.out

    $ mpirun -np $N $BIN/pw.x -nk $NK <unitcell.nscf.in  >unitcell.nscf.out


6. Calculate e-d interaction matrix element

   Use the input files ``calcmdefect.dat`` and prepared weight file ``wt.dat`` to perform matrix element calculation with ``edi.x``.


.. code-block:: console

    $ mpirun -np $N $BIN/edi.x -ni $N  >output

..

   This step will generate an output files containing the caluclation setup information and matrix element, and a postprocessing file ``pp.dat``.

   Following is an example of the output file:

.. code-block:: console

     Program EDI v.1.1 starts on  1Apr2024 at 17: 7:33 

     This program is part of the open-source Quantum ESPRESSO suite

     Parallel version (MPI), running on    16 processors

     MPI processes distributed on     1 nodes
     path-images division:  nimage    =      16
     242322 MiB available memory on the printing compute node when the environment starts


     Reading xml data from directory:

     dout/mx2.save/

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation= PBE
                           (   1   4   3   4   0   0   0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want


     G-vector sticks info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Sum         397     397    151                51529    51529   12137

     Using Slab Decomposition

     ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D
      The code is running with the 2D cutoff
      Please refer to:
      Sohier, T., Calandra, M., & Mauri, F. (2017), 
      Density functional perturbation theory for gated two-dimensional heterostructures:
      Theoretical developments and application to flexural phonons in graphene.
      Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448
     ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D
     WT files:scfwt.dat  MD5 sum:1912e3e22fde8d136b8032730343bf21
     Potential files:V_d.dat  MD5 sum:d83e0171d8167d2cee7ac12b11a4ba15
     Potential files:V_p.dat  MD5 sum:f6e9538660ee2e3a792a297b12eb2aad
     V_d_shift,V_p_shift   1.5199210366590288        1.5232846408985841     
    
                                ----------------------------
    
                                 Start M calculation k loop
    
                                ----------------------------
                                    Neutral defect
    
                           The matrix elements are in the following format:
     Mif, band and k point index of |phi_i>,  band and k point index of |phi_j>,  value of <phi_i|M|phi_f>
    
     Mif          14           1  ->              14           1                     (0.246589541,0.00000000)  
     Mif          14           1  ->              14           2                   (0.123558328,-0.213141322)  
     Mif          14           1  ->              14           3              (-4.693769291E-02,-0.240326583)  
     Mif          14           1  ->              14           4               (0.237022981,-4.940619692E-02)  
     Mif          14           1  ->              14           5                (0.223071113,8.340778202E-02)  
     Mif          14           1  ->              14           6               (8.643350005E-02,-0.216425717)  
     Mif          14           1  ->              14           7               (0.234236583,-6.557287276E-02)  
     Mif          14           1  ->              14           8              (-0.227120012,-9.330487996E-02)  
     Mif          14           1  ->              14           9                   (-0.126966208,0.209361851)  
     Mif          14           1  ->              14          10                (0.231537133,6.921063364E-02)  
     Mif          14           1  ->              14          11                    (0.104719497,0.212930799)  
     Mif          14           1  ->              14          12                  (-0.162798613,-0.165025607)  
     Mif          14           1  ->              14          13               (-4.728108644E-03,0.225269243)  
     Mif          14           1  ->              14          14              (-0.234879598,-3.820278868E-02)  
     Mif          14           1  ->              14          15              (-0.224445209,-9.087809175E-02)  
     Mif          14           1  ->              14          16                   (0.187783465,-0.157618567)  

..

   The first part is header information. Followed by the run time information, including parallelization, DFT system setup, functional, FFT grid size, and other optional DFT parameters.
   The hash value of the EDI required files are also printed. 
   The main components of the output file consists of matrix element calculation .



Mobility calculation
--------------------

7. Calcualte the carrier mobility

   Previous calculation gives ``pp.dat`` file, use this file and the postprocessing script ``mu.py`` to calculate the carrier mobility.
   Current supported model is MRTA. Other models such as iterative BTE methods are under development. 


.. code-block:: console

    $ python mu.py 
..


   This step will generate an output file containing the mobility, as well as the scattering rate which is ready to be plotted.
