Python-MIC is a Python MD demo with energy (and forces) computed in C++.

Dependencies
------------

Recent versions of the following: Python, Cython, Numpy, a C++ compiler.


How to build the extension (in-place)
-------------------------------------

.. code-block:: bash

    ./setup.py build_ext -i


How to run the example NVE MD simulation
----------------------------------------

.. code-block:: bash

    ./verlet.py
