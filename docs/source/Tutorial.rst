Tutorial
========

In the example, we provide the instructions and files needed to calculate carrier mobility limited by S-vacancy of MoS\ :sub:`2`\ .
After downloading the package, the files could be found in ``example`` folder.

A general step-by-step workflow of using EDI is given below: 

.. .. graphviz::
   digraph {
      "From" -> "To";
   }



-     Create the supercell structures

-     Calculate the perturbation potential from supercell structures

-     Generate needed k points

-     Calculate wavefunctions from unitcell structures

-     Calculate EDI matrix element

-     Postprocessing to get transport properties

When running the example, first follow the instructions in the Installation section and compile ``pw.x``, ``pp.x``, and ``edi.x``.
Assume they are located in ``$BIN`` folder,
then copy the files in the example to a working directory to perform the calculation.



Create the supercell structures
------------------------------------

The first step of studying a defect in a material is setting up the supercell models for both the defect and pristine structure. 
In this example, we provide the input files modelling a S-vacancy of MoS\ :sub:`2`\  in a 6-by-6 supercell.


.. note::
   For higher accuracy, the defect should be located in the center region of the supercell. 
   And in the relaxation of the defect system, it is advisable to fix some atoms in the region far away from the defect.

Calculate the perturbation potential from supercell structures
------------------------------------------------------------------------

 
After we obtain the relaxed structures, we can calculate the potential files.
The structures provided in the example are relaxed.
Below are command lines one could use for the calculations.

The first step is scf calculation.
   Use the input files ``pristine.scf.in`` and ``defect.scf.in`` to perform scf calculations for the pristine and relaxed defect structures.

.. code-block:: console

    $ mpirun -np $N $BIN/pw.x <prinstine.scf.in  >prinstine.scf.out

    $ mpirun -np $N $BIN/pw.x <defect.scf.in     >defect.scf.out


The following step is postprocessing.
   Use the input files ``pristine.pp.in`` and ``defect.pp.in`` to generate potential data files for the pristine and defect structures, which are ``V_p.dat`` and ``V_d.dat`` respectively.

.. code-block:: console

    $ mpirun -np $N $BIN/pp.x <prinstine.pp.in  >prinstine.pp.out

    $ mpirun -np $N $BIN/pp.x <defect.pp.in     >defect.pp.out

.. note::
   For SOC calculations, one also need to generate the exchange-correlation magnetization files.
   The needed input file could be adapted from the above files by changing the ``plot_num`` to 13.
   Note there are 3 directions so 3 data files would be generated.

Generate the needed k points  
----------------------------------



To generate the needed k points, one needs to first calculate the band structure for full k grid.
   Use the input files ``unitcell.scf.in``, ``unitcell.w90.in`` to generate a coarse grid band structure.

.. code-block:: console

    $ mpirun -np $N $BIN/pw.x <unitcell.scf.in  >unitcell.scf.out

    $ mpirun -np $N $BIN/pw.x <unitcell.w90.in  >unitcell.w90.out

..

Then use the input files  ``w90_mos2.win``, ``pw2w90.in``  to generate a fine grid band structure.  This example is using ``wannier90.x``, assuming under the same directory ``$BIN``.

.. code-block:: console

    $ mpirun -np $N $BIN/wannier90.x -pp w90_mos2

    $ mpirun -np $N $BIN/pw2wannier90.x <pw2w90.in >pw2w90.out 

    $ mpirun -np $N $BIN/wannier90.x     w90_mos2

    $ mpirun -np $N $BIN/postw90.x     w90_mos2


..

   The generated velocity file ``w90_mos2_geninterp.dat`` is needed for next step.

With the data, we can continue to generate the needed k points for e-d interaction matrix elements with triangular integral method.
   Use the script ``wmat.py`` to obtain a list of k points in file ``kpt.dat``. This file contains the k points of all the initial and final wavefunctions combined. Another file containing the list of k point pairs between initial and final k points, that are needed to calculate the e-d interaction matrix element together with the weight is generated named ``wt.dat``.


.. code-block:: console

    $ python wmat.py

..

    Prepare the input files for nscf calculation with the data from ``kpt.dat`` file, by inserting the k point lists from ``kpt.dat`` file..
    This resulted file is provided as ``unitcell.nscf.in``.

.. note::
   The wmat.py script is provided in the package. Here are a few important paramters in the file. 
   ``velocity_fn``: the interpolated band structure file name;
   ``triangular_wt``: True to use triangular integral for 2D systems, false to use Gaussian smearing.
   ``nkf``: the fine grid number;
   ``nbnd_valence``: valence band number in scf code;
   ``nbnd_valence_w90``: valence band number in Wannier results.
   The parameters should be set correctly to generate the 



Calculate wavefunctions from unitcell structures
------------------------------------------------

To generate the wavefunctions for the needed k points, one may use the input files ``unitcell.scf.in`` and ``unitcell.nscf.in`` to generate wavefunctions with proper k points.

.. code-block:: console

    $ mpirun -np $N $BIN/pw.x         <unitcell.scf.in  >unitcell.scf.out

    $ mpirun -np $N $BIN/pw.x -nk $NK <unitcell.nscf.in  >unitcell.nscf.out


Calculate EDI matrix element
--------------------------------

After all the above data are prepared, we may calculate e-d interaction matrix element.
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




Postprocessing to get transport properties
------------------------------------------------

Finally, we  can calcualte the carrier mobility and conductivity.
Previous calculation gives ``pp.dat`` file, use this file and the postprocessing script ``mu.py`` to calculate the carrier transport property.
Current supported model is MRTA and SERTA. Other models such as iterative BTE methods are under development. 


.. code-block:: console

    $ python mu.py 
..


   This step will generate an output file containing the mobility, conductivity, as well as the scattering rate, velocity, matrix element etc, which is ready to be plotted.

   An example of the data file showing the final result is shown below:

.. code-block:: console

     mu: 38.039856, sigma: 2.50257e6, nc: 1.829376e+11
     ni ki: kxyz(crystal), gamma (s^-1), E(eV) ,  vx vy(cm/s),vx vy(au), f, df(eV^-1), dos, sum m, sum m*angleterm,   Nkf
    14 274 3.05555552000000e-01 3.12500000000000e-01  7.36583698498118e+13 -2.30667785000000e-01 -3.76223686598431e+06 -5.62481368991860e+06 -8.84246469000000e-01 -1.32201183000000e+00  8.89812147616810e-03   3.44490035561073e-01  1.21347127229965e+00     9.21033388684503e+03   8.35614505146145e+03   350  5.73365571094954e-01   
    14 17 2.77777791000000e-01 3.40277791000000e-01  7.43338596702315e+13 -2.07558706000000e-01 -5.70217550230501e+06 -3.12028566625926e+06 -1.34019434000000e+00 -7.33367324000000e-01  3.62714556241616e-03   1.41171460058016e-01  1.26777099679015e+00     9.92550586468186e+03   8.52836643272525e+03   407  6.71131911250348e-01   
    14 385 2.91666657000000e-01 3.47222209000000e-01  7.48096831276885e+13 -2.53816575000000e-01 -5.32297846410817e+06 -1.22312930764192e+06 -1.25107086000000e+00 -2.87474662000000e-01  2.16952573261703e-02   8.29084888114122e-01  1.16154924399571e+00     8.01721768830764e+03   7.40936120862637e+03   299  5.47680126140007e-01   
    14 479 3.12500000000000e-01 3.26388896000000e-01  9.20002036066427e+13 -2.82599092000000e-01 -2.63754005336711e+06 -2.71704457723223e+06 -6.19906604000000e-01 -6.38592720000000e-01  6.38995500741894e-02   2.33657803025412e+00  1.08348767520257e+00     6.80040648561869e+03   6.55335519504726e+03   215  5.58218390563172e-01   
    14 2 2.70833343000000e-01 3.54166657000000e-01  7.11520247389928e+13 -2.01504454000000e-01 -6.19011619258427e+06 -1.57203385516821e+06 -1.45487607000000e+00 -3.69478434000000e-01  2.86542298264522e-03   1.11609856788115e-01  1.28053214666691e+00     1.01619076182392e+04   8.47963673539922e+03   419  6.79864609646184e-01   

..


   The first line of the file shows the mobility, conductivity, and carrier concentration at the current Fermi level.
   The user can tune the Fermi level and calculate different conductivities.
   The second line indicates what data are stored in the following data: from left to right, they represent the index and coordinates of the initial k points, the scattering rate calculated from MRTA, band structure information including energy and velocity, Fermi velocity, summation of the matrix element over final k points, and total number of final k points.

   Another example of the data file for each individual initial k point is given beblow:


.. code-block:: console

    # ni, ki: kxyz(crystal), gamma (s^-1), E(eV) ,  vx vy(cm/s),vx vy(au), f, df(eV^-1), dos, sum m, sum m*angleterm,   Nkf
    # 14 217 : 2.77777791000000e-01 3.47222209000000e-01  7.55370034990726e+13 -2.17326492000000e-01 -5.75554761388175e+06 -2.12220390635853e+06 -1.35273850000000e+00 -4.98786062000000e-01  5.30323243605497e-03   2.06058912569693e-01  1.24619904485609e+00     9.53615703575223e+03   8.25602260716569e+03   383 
    # nf kf: kfx, kfy(crystal), Wt, M1 M2, |M|, arg(M),vf1,vf2, E(eV), angleterm(1-cos(theta)), thetai, thetaf
    14 132 6.25000000000000e-01 6.52777791000000e-01     4.37359465000000e-03    2.66175210440000e-02 -2.35831102736000e-01  2.37328467412162e-01 -1.45840508330738e+00     -1.23103797000000e+00 -1.03775382000000e+00     -2.09001586000000e-01   5.96553526372161e-02 3.49484899460205e+00  3.84200342668222e+00 
    14 46 3.33333343000000e-01 3.88888896000000e-01     4.28812673000000e-05    -4.23372957200000e+00 8.69354682240000e+00  9.66965472196854e+00 2.02398737565657e+00     -1.61909744000000e-01 1.62098110000000e+00     -1.94189295000000e-01   1.25098994182196e+00 3.49484899460205e+00  1.67034986998397e+00 
    14 285 3.81944448000000e-01 3.12500000000000e-01     2.88865599000000e-03    -6.79999668160000e-02 8.10520752576000e-04  6.80047971165818e-02 3.12967378936845e+00     1.58664262000000e+00 1.04926139000000e-01     -2.30619013000000e-01   1.95903465875036e+00 3.49484899460205e+00  6.60347705282734e-02 

..

   The first two lines shows the information of the initial k point, which is the same as the overall data.
   The third line indicates the conten stored in the following: the index and coordinates of the final k points, the weight used in calculation of scattering rate, and the matrix element of scattering between the initial and final states.

.. note::
   For details on the input parameters of the ``mu.py`` script, see the Input section.
