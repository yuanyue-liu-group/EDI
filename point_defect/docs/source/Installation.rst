Installation
===============

Download
--------

EDI is a freeware licensed under GPL v3.0.
Its stable version can be downloaded as zip file from the `main website <https://zhongcanxiao.github.io/EDI>`_.
Alternatively, the most up to date version could be accessed from `github page <https://github.com/zhongcanxiao/EDI>`_. 
However, stable version is encouraged to use since develop version may not be fully tested.

Prerequisites
-------------

The following packages are required to install EDI:

* A GCC Fortran compiler

* MPI 

* HDF5

* Quantum Espresso

Compile Quantum Espresso
---------------------------

Currently, Quantum Espresso package is required to run EDI.
Additional required library hdf5 should be compiled with QE.
Supported version of Quantum Espresso is  `6.8 <https://gitlab.com/QEF/q-e/-/archive/qe-6.8/q-e-qe-6.8.tar.gz>`_. 
To compile edi, first download the respective version. 


To compiled QE with hdf5 support:

.. code-block:: console

   $ tar xvf q-e-qe-6.8.tar.gz 
   $ cd q-e-qe-6.8
   $ configure --with-hdf5  
   $ make --with-hdf5 pw

.. note::
    When encountering problem compiling hdf5, try to use ``h5fc -show`` to find hdf5 library.

Compilation
------------

After QE is compiled, copy the EDI into QE root directory.

.. code-block:: console

   $ cp -r $src/EDI . 
 
where ``$src`` is the directory containing the EDI main folder.
Change into the main directory of EDI and compile.
EDI makes use of HDF5 to access the dielectric function data.
Note the gfortran compiler with Fortran 2008 support is required.

.. code-block:: console

   $ cd EDI
   $ make -j $(nproc)

where you should substitute ``nproc`` with the number of cores available for parallel compilation. 

Uninstall
----------------

To clean the compiled files, run:


.. code-block:: console

   $ make clean

.. _installation:

Compiling the documentation
---------------------------

The documentation could be downloaded from main website, as well as compiled locally.
To do this you need to have the following available on your machine:

* sphinx

* pdflatex

Then type:

.. code-block:: console

   $ cd docs
   $ make latexpdf



Installation instructions for specific systems
--------------------------------------------------------------------

Ubuntu
^^^^^^

Under QE root folder::

   $ configure --with-hdf5  --with-hdf5-include=/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include
   $ make --with-hdf5 pw

Note that paths to the HDF5 library may need to be updated.
Tested on Ubuntu 20.04.

Lonestar 6
^^^^^^^^^^

Under QE root folder::

   $ module load hdf5 fftw3 gcc mkl 
   $ ./configure --with-hdf5=$TACC_HDF5_DIR
   $ make --with-hdf5 pw

Anvil
^^^^^

Under QE root folder::

   $ module load hdf5 fftw gcc intel-mkl 
   $ ./configure --with-hdf5 --with-hdf5-libs="-lhdf5_fortran -lhdf5" --with-hdf5-include="$HDF5_HOME/include"
   $ make --with-hdf5 pw

.. note::
    When running under Anvil, if the wait time is extensively long at gw_bcast routine, it's likely the memory is out.



