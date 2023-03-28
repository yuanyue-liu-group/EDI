Installation
=====

.. _installation:

Prerequisites 
------------

To compile edi, first install Quantum Espresso. It needs to be compiled with hdf5 support.


.. code-block:: console

   $ make --with-hdf5 pw

Compilation
------------

Then change into the main directory of EDI and compile.

.. code-block:: console

   $ cd EDI
   $ make edi.x

Uninstall
----------------

To clean the compiled files, run:


.. code-block:: console

   $ make clean


