Tutorial
========

In the example, we provide the instructions and files needed to calcualte carrier mobility limited by S-vacancy of MoS\ :sub:`2`\ .

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

EDI calculation
----------------

3. Generate the needed k points for e-d interaction matrix elements with triangular integral method

   Use the script ``wmat_mos2.py`` to obtain a list of k points in file ``kpt.dat``, and a list of k point pairs to calculate the e-d interaction matrix element together with the weight ``wt.dat``.

   Prepare the input files for nscf calculation with the data from ``kpt.dat`` file.
   This resulted file is provided as ``unitcell.nscf.in``.

.. code-block:: console

    $ python wmat_mos2.py




4. Generate the wavefunctions for the needed k points

   Use the input files ``unitcell.scf.in`` and ``unitcell.nscf.in`` to generate wavefunctions with proper k points.

.. code-block:: console

    $ mpirun -np $N $BIN/pw.x         <unitcell.scf.in  >unitcell.scf.out

    $ mpirun -np $N $BIN/pw.x -nk $NK <unitcell.nscf.in  >unitcell.nscf.out

5. Calculate e-d interaction matrix element

   Use the input files ``calcmdefect.dat`` and prepared weight file ``wt.dat`` to perform matrix element calculation with ``edi.x``.


.. code-block:: console

    $ mpirun -np $N $BIN/edi.x -ni $N  >output


Mobility calculation
--------------------

6. Use MRTA model to calcualted the carrier mobility

   Previous calculation gives ``pp.dat`` file, use this file and the postprocessing script ``mu.py`` to calculate the carrier mobility.


.. code-block:: console

    $ python mu.py 

This step will generate an output file containing the mobility, as well as the scattering rate which is ready to be plot.
