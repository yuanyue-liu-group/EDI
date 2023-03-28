Usage
=====

.. _installation:

Installation
------------

To compile edi, first install Quantum Espresso:

.. code-block:: console

   (.venv) $ make --with-hdf5 pw
   (.venv) $ cd EDI
   (.venv) $ make edi.x

Creating recipes
----------------

To retrieve a list of random ingredients,
you can use the ``lumache.get_random_ingredients()`` function:

.. autofunction:: lumache.get_random_ingredients

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
will raise an exception.

.. autoexception:: lumache.InvalidKindError

For example:

>>> import lumache
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']

