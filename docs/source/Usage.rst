Usage
=====

.. _installation:

Running
------------

EDI requires the wavefunction and perturbation potential obtained from first principle code as input data for calculations.
Currently, QE is supported DFT engine.
Here is a step-by-step instruction for generating the required data.

Preparation
^^^^^

.. code-block:: console

    $ mpirun -np $N edi.x -ni $N  >output


calculated
 
To run edi, use the following command:

.. code-block:: console

    $ mpirun -np $N edi.x -ni $N  >output

The input file is predefined text file. 
The output file contains calculation procedure information.

Parallelization level
------------

EDI is parallelized over k point pairs.
EDI uses optimized algorithm to calculate different part of scattering matrix element, giving the optimal performance of accuracy and effeciency.

Optimization and scaling
----------------
The scalability of EDI is linear with system size. 
Namelhy the calcualtion cost scales linearly with the volum of super cell, making it capable of calculating large systems easily.
Calculation of matrix element is parallelizede over k point pairs.  
EDI could easily run on large HPCs and utilize the full capacity.
Here the number of core is `$N`, it could be any value divisible by the total number of kpoint pairs for optimal performance.


Postprocessing
------------

After a successful run, the output files will be generated. 
The default file name is pp.dat.
It could be used to calculate the carrier mobility.

Calculating mobility
------------

To calcualted the transport mobility, the postprocessing script mu.py could be used.


.. code-block:: console

    $ python mu.py 


Out put data containing the transport properties, including scattering rate and mobility, will be written automatically.
