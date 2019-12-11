A Python wrapper of ``minpack``'s ``hybrd1`` function.

******
Notice
******

This product includes software developed by the University of Chicago, as Operator of Argonne National Laboratory.

More information can be found in "NOTICE".

************
Installation
************

Pip
===

This package can be installed with ``pip`` via the setup function provided by ``numpy``.

Simply execute

.. code-block:: bash

    $ pip install .

within this directory.

**NOTE**: While this worked nicely and out-of-the-box on Linux, on Windows some additional environment variables have to be set. However, I never got it working properly, so I relied on the manual compilation explained below.

Linux
=====

1) Install ``python`` (and ``pip``)
2) Install ``numpy``:

.. code-block:: bash

    $ pip install numpy

3) Run Makefile:

.. code-block:: bash

    $ make linux

Windows
=======

1) Download and install `Anaconda`_

.. _Anaconda: https://www.anaconda.com/distribution/

2) Download and install `MinGW`_ to ``C:\mingw`` with the following options:

   - Architecture: ``x86_64``
   - Threads: ``posix``
   - Exception: ``seh``

.. _MinGW: https://sourceforge.net/projects/mingw-w64/

3) Add the following to the ``Path`` environment variable:
    - ``C:\mingw\mingw64\bin``
    - ``C:\Users\USERNAME\Anaconda3``
    - ``C:\Users\USERNAME\Anaconda3\Scripts``

4) Make sure that you have the correct C++ Redistributable installed:

   +------------+----------+
   | Visual C++ | CPython  |
   +============+==========+
   | 14.X       | 3.5+     |
   +------------+----------+
   | 10.0       | 3.3, 3.4 |
   +------------+----------+
   | 9.0        | 3.2-     |
   +------------+----------+

See `Microsoft Visual Cpp Build Tools`_

.. _Microsoft Visual Cpp Build Tools: https://visualstudio.microsoft.com/downloads/

5) Run Makefile:

.. code-block:: bash

    $ make windows

*******
Testing
*******

The following script can be used to test the ``optimize`` module:

.. code-block:: python

    #!/usr/bin/python3

    import numpy as np
    import optimize as opt

    def f(x, a):
        if isinstance(a, list):
            # scipy.optimize.root passes a list
            a = a[0]

        fvec = np.array([a*np.sin(x[0]), a*np.cos(x[1])])
        return fvec

    print(opt.root(f, np.array([.5, .5]), args=[1.]))

    # For comparison
    import scipy.optimize as sscipy_opt

    print(scipy_opt.root(f, np.array([.5, .5]), args=[1.]))

**********************
Additional information
**********************

See also:

- `Configuring your compilers`_
- `Testing your compilers`_
- `Example of setup file`_

.. _Configuring your compilers: https://python-at-risoe.pages.windenergy.dtu.dk/compiling-on-windows/configuration.html#configuration
.. _Testing your compilers: https://python-at-risoe.pages.windenergy.dtu.dk/compiling-on-windows/testing.html#testing
.. _Example of setup file: https://github.com/scipy/scipy/blob/81c096001974f0b5efe29ec83b54f725cc681540/scipy/fftpack/setup.py
