========
Building
========

I recommend you build the code in a separate directory, e.g., in the project
directory run

.. code-block:: bash

    $ mkdir build
    $ cd build

When inside your build directory run

.. code-block:: bash

    cmake -B . -S Path_to_project_directory
    make gsl-asw

This will build the gsl-asw executable in the bin directory.

In order to build tests you can add the CMake flag

::

    -DBUILD_TESTING
